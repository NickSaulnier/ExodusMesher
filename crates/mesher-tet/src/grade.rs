//! Grading a finished tetrahedral mesh so that neighbouring cells' sizes
//! stay within a ratio (Milestone 18): where a cell meets one more than
//! `limit` times smaller across a face (sizes measured as in
//! `quality::fv_quality`, each cell's longest edge), the larger cell's
//! longest edge is split at its midpoint, every tet around the edge split
//! in two (longest-edge bisection), until no face is over the limit or the
//! budget is spent.
//!
//! Small cells are where the geometry needs them (two feature edges passing
//! close, a sliver input triangle); the kernel's refinement follows the size
//! field away from them, which knows nothing of them, so a cell of the
//! field's size can sit beside one ten times smaller. Bisection grades the
//! larger cells down towards them instead, a wave that stops where the
//! ratio is within the limit.
//!
//! A split boundary edge splits its two boundary faces, each child keeping
//! its parent's tag, so side sets follow; the midpoint is on the edge to
//! rounding, so the boundary is the same surface. A split is skipped if a
//! new tet is not positive exactly and in floating point. Deterministic:
//! cells are taken largest first, then by index, and new vertices and tets
//! are numbered in that order.

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{norm, sub};
use mesher_geom::{Point3, Sign};

use crate::fxhash::FxHashMap;
use crate::quality::float_volume;

/// Faces of a tet, each with the vertex left out.
const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// The six edges of a tet.
const EDGES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

/// How grading works.
#[derive(Debug, Clone, PartialEq)]
pub struct GradeParams {
    /// Largest ratio between face neighbours' sizes left alone.
    pub limit: f64,
    /// Most tets grading may add, as a share of the mesh's (at least
    /// `min_budget`).
    pub budget: f64,
    pub min_budget: usize,
}

impl Default for GradeParams {
    fn default() -> GradeParams {
        GradeParams {
            limit: 4.0,
            budget: 0.5,
            min_budget: 10_000,
        }
    }
}

/// What grading did.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct GradeStats {
    /// Edges split (interior and boundary).
    pub splits: usize,
    pub boundary_splits: usize,
    /// Cells over the limit whose split was skipped (a new tet would not be
    /// positive).
    pub skipped: usize,
    /// Rounds of splitting.
    pub rounds: usize,
    /// Whether the budget ran out with faces still over the limit.
    pub out_of_budget: bool,
}

/// A graded mesh: the input's vertices first (none removed), then the new
/// ones; the tets positively oriented; the tagged boundary faces, children
/// of a split face in place of it.
#[derive(Debug, Clone)]
pub struct Graded {
    pub points: Vec<Point3>,
    pub tets: Vec<[u32; 4]>,
    pub faces: Vec<([u32; 3], u32)>,
    pub stats: GradeStats,
}

fn length(points: &[Point3], a: u32, b: u32) -> f64 {
    norm(sub(points[a as usize], points[b as usize]))
}

fn sorted3(mut f: [u32; 3]) -> [u32; 3] {
    f.sort_unstable();
    f
}

struct Mesh {
    points: Vec<Point3>,
    tets: Vec<[u32; 4]>,
    alive: Vec<bool>,
    longest: Vec<f64>,
    star: Vec<Vec<u32>>,
    /// Boundary edges (each vertex pair once, smaller first) and how many
    /// boundary faces hold each.
    boundary_edges: FxHashMap<[u32; 2], u32>,
    /// Tagged faces: oriented vertices, tag, alive; by sorted vertices.
    faces: Vec<([u32; 3], u32, bool)>,
    face_index: FxHashMap<[u32; 3], usize>,
}

impl Mesh {
    fn longest_of(&self, t: [u32; 4]) -> f64 {
        EDGES
            .iter()
            .map(|&[i, j]| length(&self.points, t[i], t[j]))
            .fold(0.0, f64::max)
    }

    fn add(&mut self, t: [u32; 4]) -> u32 {
        let id = self.tets.len() as u32;
        self.longest.push(self.longest_of(t));
        self.tets.push(t);
        self.alive.push(true);
        for v in t {
            self.star[v as usize].push(id);
        }
        id
    }

    fn remove(&mut self, id: u32) {
        self.alive[id as usize] = false;
        for v in self.tets[id as usize] {
            self.star[v as usize].retain(|&s| s != id);
        }
    }

    /// The tet across face `f` of tet `t`.
    fn across(&self, t: u32, f: [u32; 3]) -> Option<u32> {
        self.star[f[0] as usize].iter().copied().find(|&s| {
            s != t && {
                let tet = self.tets[s as usize];
                tet.contains(&f[1]) && tet.contains(&f[2])
            }
        })
    }

    fn is_boundary_edge(&self, a: u32, b: u32) -> bool {
        self.boundary_edges.contains_key(&[a.min(b), a.max(b)])
    }

    fn add_boundary_edge(&mut self, a: u32, b: u32) {
        *self.boundary_edges.entry([a.min(b), a.max(b)]).or_insert(0) += 1;
    }

    fn drop_boundary_edge(&mut self, a: u32, b: u32) {
        let key = [a.min(b), a.max(b)];
        if let Some(n) = self.boundary_edges.get_mut(&key) {
            *n -= 1;
            if *n == 0 {
                self.boundary_edges.remove(&key);
            }
        }
    }

    /// Split edge `ab` at its midpoint: the new tets, or `None` (and the
    /// mesh unchanged) if one would not be positive.
    fn split(&mut self, a: u32, b: u32) -> Option<Vec<u32>> {
        let ring: Vec<u32> = self.star[a as usize]
            .iter()
            .copied()
            .filter(|&t| self.tets[t as usize].contains(&b))
            .collect();
        if ring.is_empty() {
            return None;
        }
        let (pa, pb) = (self.points[a as usize], self.points[b as usize]);
        let mid = [0, 1, 2].map(|k| 0.5 * (pa[k] + pb[k]));
        let m = self.points.len() as u32;
        let at = |v: u32| {
            if v == m { mid } else { self.points[v as usize] }
        };
        let mut new = Vec::with_capacity(2 * ring.len());
        for &t in &ring {
            let tet = self.tets[t as usize];
            for child in [
                tet.map(|x| if x == a { m } else { x }),
                tet.map(|x| if x == b { m } else { x }),
            ] {
                let p = child.map(at);
                if orient3d(p[0], p[1], p[2], p[3]) != Sign::Positive || float_volume(p) <= 0.0 {
                    return None;
                }
                new.push(child);
            }
        }
        // Tagged boundary faces on the edge split in two, keeping the tag.
        let mut children = Vec::new();
        if self.is_boundary_edge(a, b) {
            for &t in &ring {
                let tet = self.tets[t as usize];
                for c in tet {
                    if c == a || c == b {
                        continue;
                    }
                    if let Some(&i) = self.face_index.get(&sorted3([a, b, c])) {
                        if self.faces[i].2 {
                            children.push(i);
                        }
                    }
                }
            }
            children.sort_unstable();
            children.dedup();
        }
        self.points.push(mid);
        self.star.push(Vec::new());
        for &t in &ring {
            self.remove(t);
        }
        let ids = new.into_iter().map(|t| self.add(t)).collect();
        for i in children {
            let (f, tag, _) = self.faces[i];
            self.faces[i].2 = false;
            self.face_index.remove(&sorted3(f));
            for k in 0..3 {
                let (u, v) = (f[k], f[(k + 1) % 3]);
                self.drop_boundary_edge(u, v);
            }
            for child in [
                f.map(|x| if x == a { m } else { x }),
                f.map(|x| if x == b { m } else { x }),
            ] {
                self.face_index.insert(sorted3(child), self.faces.len());
                self.faces.push((child, tag, true));
                for k in 0..3 {
                    self.add_boundary_edge(child[k], child[(k + 1) % 3]);
                }
            }
        }
        Some(ids)
    }

    /// Whether a face of tet `t` meets a cell more than `limit` times
    /// smaller or larger; the larger cells found are pushed to `large`.
    fn check(&self, t: u32, limit: f64, large: &mut Vec<u32>) {
        let tet = self.tets[t as usize];
        let lt = self.longest[t as usize];
        for f in FACES {
            let face = f.map(|k| tet[k]);
            let Some(u) = self.across(t, face) else {
                continue;
            };
            let lu = self.longest[u as usize];
            if lt > limit * lu && lu > 0.0 {
                large.push(t);
            } else if lu > limit * lt && lt > 0.0 {
                large.push(u);
            }
        }
    }
}

/// Grade mesh `tets` on `points` (positively oriented), with its boundary
/// faces tagged by `faces` (outward or not; a face not listed has no tag and
/// is split untagged).
#[must_use]
pub fn grade(
    points: &[Point3],
    tets: &[[u32; 4]],
    faces: &[([u32; 3], u32)],
    params: &GradeParams,
) -> Graded {
    let _prof = crate::prof::scope("grade::grade");
    let n = points.len();
    let mut mesh = Mesh {
        points: points.to_vec(),
        tets: Vec::with_capacity(tets.len()),
        alive: Vec::with_capacity(tets.len()),
        longest: Vec::with_capacity(tets.len()),
        star: vec![Vec::new(); n],
        boundary_edges: FxHashMap::default(),
        faces: Vec::new(),
        face_index: FxHashMap::default(),
    };
    let mut stats = GradeStats::default();
    if tets.iter().flatten().any(|&v| v as usize >= n) {
        return Graded {
            points: points.to_vec(),
            tets: tets.to_vec(),
            faces: faces.to_vec(),
            stats,
        };
    }
    for &t in tets {
        mesh.add(t);
    }
    // The boundary: faces with one tet.
    let mut count: FxHashMap<[u32; 3], u32> = FxHashMap::default();
    for t in tets {
        for f in FACES {
            *count.entry(sorted3(f.map(|k| t[k]))).or_insert(0) += 1;
        }
    }
    let tags: FxHashMap<[u32; 3], ([u32; 3], u32)> =
        faces.iter().map(|&(f, g)| (sorted3(f), (f, g))).collect();
    for t in tets {
        for f in FACES {
            let face = f.map(|k| t[k]);
            let key = sorted3(face);
            if count.get(&key) == Some(&1) {
                for k in 0..3 {
                    mesh.add_boundary_edge(face[k], face[(k + 1) % 3]);
                }
                // The caller's orientation for listed faces.
                let (oriented, tag) = tags.get(&key).copied().unwrap_or((face, u32::MAX));
                mesh.face_index.insert(key, mesh.faces.len());
                mesh.faces.push((oriented, tag, true));
            }
        }
    }
    let budget = ((params.budget * tets.len() as f64) as usize).max(params.min_budget);
    let start = tets.len();

    let mut work: Vec<u32> = (0..tets.len() as u32).collect();
    while !work.is_empty() {
        let mut large = Vec::new();
        for &t in &work {
            if mesh.alive[t as usize] {
                mesh.check(t, params.limit, &mut large);
            }
        }
        if large.is_empty() {
            break;
        }
        large.sort_unstable_by(|&x, &y| {
            mesh.longest[y as usize]
                .total_cmp(&mesh.longest[x as usize])
                .then(x.cmp(&y))
        });
        large.dedup();
        stats.rounds += 1;
        let mut next = Vec::new();
        for t in large {
            if !mesh.alive[t as usize] {
                continue;
            }
            if mesh.tets.len() - start >= budget {
                stats.out_of_budget = true;
                break;
            }
            let tet = mesh.tets[t as usize];
            let [i, j] = EDGES
                .iter()
                .copied()
                .max_by(|&[i, j], &[k, l]| {
                    length(&mesh.points, tet[i], tet[j])
                        .total_cmp(&length(&mesh.points, tet[k], tet[l]))
                        .then_with(|| {
                            let key = |p: u32, q: u32| (p.min(q), p.max(q));
                            key(tet[k], tet[l]).cmp(&key(tet[i], tet[j]))
                        })
                })
                .unwrap_or([0, 1]);
            let (a, b) = (tet[i], tet[j]);
            let boundary = mesh.is_boundary_edge(a, b);
            match mesh.split(a, b) {
                Some(ids) => {
                    stats.splits += 1;
                    stats.boundary_splits += usize::from(boundary);
                    next.extend(ids);
                }
                None => stats.skipped += 1,
            }
        }
        if stats.out_of_budget {
            break;
        }
        next.retain(|&t| mesh.alive[t as usize]);
        work = next;
    }

    let tets = mesh
        .tets
        .iter()
        .zip(&mesh.alive)
        .filter(|(_, a)| **a)
        .map(|(t, _)| *t)
        .collect();
    let faces = mesh
        .faces
        .iter()
        .filter(|f| f.2 && f.1 != u32::MAX)
        .map(|f| (f.0, f.1))
        .collect();
    Graded {
        points: mesh.points,
        tets,
        faces,
        stats,
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

    use super::*;
    use crate::quality::fv_quality;

    fn volume(p: &[Point3], t: &[[u32; 4]]) -> f64 {
        t.iter()
            .map(|t| float_volume(t.map(|v| p[v as usize])))
            .sum()
    }

    /// Boundary faces (sorted) of a mesh.
    fn boundary(t: &[[u32; 4]]) -> Vec<[u32; 3]> {
        let mut count: std::collections::BTreeMap<[u32; 3], u32> = Default::default();
        for t in t {
            for f in FACES {
                *count.entry(sorted3(f.map(|k| t[k]))).or_insert(0) += 1;
            }
        }
        count
            .into_iter()
            .filter(|(_, n)| *n == 1)
            .map(|(f, _)| f)
            .collect()
    }

    /// The Delaunay tetrahedralization of a unit cube's corners and a
    /// cluster of size `size` at `at`: cells of the cube's size beside
    /// cells of the cluster's.
    fn jump_at(at: Point3, size: f64) -> (Vec<Point3>, Vec<[u32; 4]>) {
        let mut p: Vec<Point3> = (0..8)
            .map(|i| [f64::from(i & 1), f64::from((i >> 1) & 1), f64::from(i >> 2)])
            .collect();
        for i in 0..8 {
            let c = [f64::from(i & 1), f64::from((i >> 1) & 1), f64::from(i >> 2)];
            p.push([0, 1, 2].map(|k| at[k] + size * (c[k] - 0.5 + 0.1 * f64::from(i))));
        }
        let t = crate::delaunay::delaunay(&p).unwrap().tets;
        (p, t)
    }

    /// The cluster at the cube's centre.
    fn jump() -> (Vec<Point3>, Vec<[u32; 4]>) {
        jump_at([0.5; 3], 0.01)
    }

    /// A cluster at the centre is graded from inside, no boundary edge
    /// split; one with a vertex on a face has its boundary edges split too,
    /// and the faces' tags follow.
    #[test]
    fn grading_brings_neighbours_within_the_limit() {
        for (at, size, on_boundary) in [([0.5; 3], 0.01, false), ([0.5, 0.45, 0.0005], 0.001, true)]
        {
            grades(at, size, on_boundary);
        }
    }

    fn grades(at: Point3, size: f64, on_boundary: bool) {
        let (p, t) = jump_at(at, size);
        let before = fv_quality(&p, &t).max_size_ratio;
        assert!(before > 4.0, "{before}");
        let faces: Vec<([u32; 3], u32)> = boundary(&t)
            .into_iter()
            .enumerate()
            .map(|(i, f)| (f, i as u32 % 3))
            .collect();
        let area = |p: &[Point3], f: &[([u32; 3], u32)], g: u32| -> f64 {
            f.iter()
                .filter(|x| x.1 == g)
                .map(|(f, _)| {
                    let [a, b, c] = f.map(|v| p[v as usize]);
                    0.5 * norm(mesher_geom::vec3::cross(sub(b, a), sub(c, a)))
                })
                .sum()
        };
        let g = grade(&p, &t, &faces, &GradeParams::default());
        let after = fv_quality(&g.points, &g.tets);
        assert!(after.max_size_ratio <= 4.0, "{}", after.max_size_ratio);
        assert!(
            g.stats.splits > 0 && !g.stats.out_of_budget,
            "{:?}",
            g.stats
        );
        assert_eq!(g.stats.boundary_splits > 0, on_boundary, "{:?}", g.stats);
        for t in &g.tets {
            let q = t.map(|v| g.points[v as usize]);
            assert_eq!(orient3d(q[0], q[1], q[2], q[3]), Sign::Positive);
        }
        assert!((volume(&g.points, &g.tets) - volume(&p, &t)).abs() < 1e-12);
        // Every boundary face is tagged, and each tag keeps its area.
        assert_eq!(boundary(&g.tets).len(), g.faces.len());
        for tag in 0..3 {
            assert!((area(&g.points, &g.faces, tag) - area(&p, &faces, tag)).abs() < 1e-12);
        }
        // Deterministic.
        let again = grade(&p, &t, &faces, &GradeParams::default());
        assert_eq!(again.tets, g.tets);
        assert_eq!(again.faces, g.faces);
    }

    #[test]
    fn a_graded_mesh_is_unchanged() {
        let (p, t) = jump();
        let g = grade(&p, &t, &[], &GradeParams::default());
        let h = grade(&g.points, &g.tets, &[], &GradeParams::default());
        assert_eq!(h.stats.splits, 0);
        assert_eq!(h.tets, g.tets);
    }

    #[test]
    fn the_budget_stops_grading() {
        let (p, t) = jump();
        let g = grade(
            &p,
            &t,
            &[],
            &GradeParams {
                budget: 0.0,
                min_budget: 1,
                ..GradeParams::default()
            },
        );
        assert!(g.stats.out_of_budget);
        assert_eq!(g.stats.splits, 1);
    }
}
