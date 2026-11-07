//! Isosurface stuffing (Labelle and Shewchuk, "Isosurface Stuffing: Fast
//! Tetrahedral Meshes with Good Dihedral Angles", SIGGRAPH 2007), uniform
//! lattice version.
//!
//! 1. Fill the bounding box with a body-centred cubic (BCC) lattice: black
//!    vertices at cube corners, red vertices at cube centres, and four tets
//!    around every red edge (between centres of face-adjacent cubes).
//! 2. Evaluate the level function at every lattice vertex (negative inside).
//! 3. Cut points: on every lattice edge whose ends have strictly opposite
//!    signs, find the zero crossing by bisection.
//! 4. Warp: a vertex closer to a cut point on one of its edges than α × the
//!    edge's length (α = 0.24999 for long edges, 0.41189 for short) moves onto
//!    the nearest such cut point and becomes a zero vertex; cut points on its
//!    other edges are discarded.
//! 5. Stencils: a tet with no outside vertex is kept whole if it has an
//!    inside vertex (all-zero tets are dropped: they are slivers lying in the
//!    surface). Otherwise its inside part (vertices with negative or zero
//!    level, plus its cut points) is a small polytope. Its boundary is
//!    triangulated face by face (a quadrilateral, which only appears on a face
//!    with two inside vertices and one outside vertex, is split by the
//!    diagonal from its lower-numbered inside vertex), and the polytope is
//!    filled by coning from its lowest-numbered inside vertex. Both rules
//!    depend only on the shared face, so neighbouring tets conform.
//!
//! The level function is any `Fn(Point3) -> f64`; the mesher uses
//! `0.5 − winding number`.

use std::collections::{BTreeMap, HashMap};

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{add, distance_squared, scale, sub};
use mesher_geom::{Aabb, Point3, Sign};
use rayon::prelude::*;

/// Warping thresholds from the paper.
pub const ALPHA_LONG: f64 = 0.24999;
pub const ALPHA_SHORT: f64 = 0.41189;
/// The dihedral-angle bounds the paper proves for these thresholds.
pub const DIHEDRAL_MIN_DEGREES: f64 = 10.7;
pub const DIHEDRAL_MAX_DEGREES: f64 = 164.8;

#[derive(Debug, Clone, PartialEq)]
pub struct StuffingParams {
    /// Lattice cube edge length (absolute).
    pub cell_size: f64,
    pub alpha_long: f64,
    pub alpha_short: f64,
    /// Refuse lattices with more vertices than this.
    pub max_lattice_vertices: u64,
    /// Bisection steps per cut point (each halves the error).
    pub bisection_steps: u32,
}

impl StuffingParams {
    #[must_use]
    pub fn with_cell_size(cell_size: f64) -> StuffingParams {
        StuffingParams {
            cell_size,
            alpha_long: ALPHA_LONG,
            alpha_short: ALPHA_SHORT,
            max_lattice_vertices: 40_000_000,
            bisection_steps: 52,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StuffingError {
    /// Cell size not positive and finite.
    BadCellSize,
    /// Nothing to mesh (empty bounding box).
    EmptyBounds,
    /// The lattice would exceed `max_lattice_vertices`.
    Budget { lattice_vertices: u64, limit: u64 },
}

impl std::fmt::Display for StuffingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StuffingError::BadCellSize => write!(f, "cell size must be positive and finite"),
            StuffingError::EmptyBounds => write!(f, "empty bounding box"),
            StuffingError::Budget {
                lattice_vertices,
                limit,
            } => write!(
                f,
                "lattice of {lattice_vertices} vertices exceeds the limit of {limit}"
            ),
        }
    }
}

/// The stuffed mesh: positively oriented tets.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Stuffed {
    pub points: Vec<Point3>,
    pub tets: Vec<[u32; 4]>,
    pub lattice_vertices: u64,
    pub cut_points: u64,
    pub warped: u64,
    /// Output tets of zero volume, dropped (the paper's guarantee says none).
    pub dropped: u64,
}

struct Lattice {
    origin: Point3,
    h: f64,
    /// Cubes per axis.
    n: [usize; 3],
}

impl Lattice {
    fn black_count(&self) -> usize {
        (self.n[0] + 1) * (self.n[1] + 1) * (self.n[2] + 1)
    }

    fn red_count(&self) -> usize {
        self.n[0] * self.n[1] * self.n[2]
    }

    fn black(&self, i: [usize; 3]) -> u32 {
        (i[0] + (self.n[0] + 1) * (i[1] + (self.n[1] + 1) * i[2])) as u32
    }

    fn red(&self, i: [usize; 3]) -> u32 {
        (self.black_count() + i[0] + self.n[0] * (i[1] + self.n[1] * i[2])) as u32
    }

    fn is_red(&self, id: u32) -> bool {
        id as usize >= self.black_count()
    }

    fn point(&self, id: u32) -> Point3 {
        let id = id as usize;
        let (coords, half) = if id < self.black_count() {
            let (nx, ny) = (self.n[0] + 1, self.n[1] + 1);
            ([id % nx, (id / nx) % ny, id / (nx * ny)], 0.0)
        } else {
            let r = id - self.black_count();
            let (nx, ny) = (self.n[0], self.n[1]);
            ([r % nx, (r / nx) % ny, r / (nx * ny)], 0.5)
        };
        [0, 1, 2].map(|k| self.origin[k] + self.h * (coords[k] as f64 + half))
    }

    /// Every BCC tet: the two centres of face-adjacent cubes plus one edge of
    /// their shared face.
    fn for_each_tet(&self, mut f: impl FnMut([u32; 4])) {
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            for k in 0..self.n[2] {
                for j in 0..self.n[1] {
                    for i in 0..self.n[0] {
                        let cube = [i, j, k];
                        if cube[a] + 1 >= self.n[a] {
                            continue;
                        }
                        let mut next = cube;
                        next[a] += 1;
                        let (r1, r2) = (self.red(cube), self.red(next));
                        let corner = |db: usize, dc: usize| {
                            let mut v = cube;
                            v[a] += 1;
                            v[b] += db;
                            v[c] += dc;
                            self.black(v)
                        };
                        let q = [corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)];
                        for e in 0..4 {
                            f([r1, r2, q[e], q[(e + 1) % 4]]);
                        }
                    }
                }
            }
        }
    }
}

/// A vertex of the output before numbering: a lattice vertex or the cut
/// point on a lattice edge (by index into the sorted cut-edge list).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Node {
    Lattice(u32),
    Cut(u32),
}

const TET_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// Stuff the region where `level < 0` inside `bounds`.
pub fn stuff<F>(level: F, bounds: &Aabb, params: &StuffingParams) -> Result<Stuffed, StuffingError>
where
    F: Fn(Point3) -> f64 + Sync,
{
    let h = params.cell_size;
    if !(h > 0.0 && h.is_finite()) {
        return Err(StuffingError::BadCellSize);
    }
    if bounds.is_empty() {
        return Err(StuffingError::EmptyBounds);
    }
    // One and a half cells of margin keep the outer lattice layer outside.
    let origin = bounds.min.map(|c| c - 1.5 * h);
    let n: [usize; 3] = [0, 1, 2].map(|k| {
        let span = bounds.max[k] - bounds.min[k] + 3.0 * h;
        ((span / h).ceil() as usize).max(2)
    });
    let lattice_vertices = n
        .iter()
        .map(|&c| c as u64 + 1)
        .product::<u64>()
        .saturating_add(n.iter().map(|&c| c as u64).product::<u64>());
    if lattice_vertices > params.max_lattice_vertices {
        return Err(StuffingError::Budget {
            lattice_vertices,
            limit: params.max_lattice_vertices,
        });
    }
    let lat = Lattice { origin, h, n };
    let total = lat.black_count() + lat.red_count();

    // Signs at every lattice vertex.
    let mut sign: Vec<i8> = (0..total as u32)
        .into_par_iter()
        .map(|id| {
            let l = level(lat.point(id));
            if l < 0.0 {
                -1
            } else if l > 0.0 {
                1
            } else {
                0
            }
        })
        .collect();

    // Cut edges: lattice edges whose ends have strictly opposite signs.
    let mut cut_edges: Vec<(u32, u32)> = Vec::new();
    lat.for_each_tet(|t| {
        for (i, j) in TET_EDGES {
            let (u, v) = (t[i], t[j]);
            if sign[u as usize] * sign[v as usize] < 0 {
                cut_edges.push((u.min(v), u.max(v)));
            }
        }
    });
    cut_edges.sort_unstable();
    cut_edges.dedup();
    let cut_point: Vec<Point3> = cut_edges
        .par_iter()
        .map(|&(u, v)| {
            let (inside, outside) = if sign[u as usize] < 0 { (u, v) } else { (v, u) };
            let (p, q) = (lat.point(inside), lat.point(outside));
            let (mut lo, mut hi) = (0.0f64, 1.0f64);
            for _ in 0..params.bisection_steps {
                let mid = 0.5 * (lo + hi);
                if level(add(p, scale(sub(q, p), mid))) < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            add(p, scale(sub(q, p), 0.5 * (lo + hi)))
        })
        .collect();
    let cut_index = |u: u32, v: u32| -> Option<u32> {
        cut_edges
            .binary_search(&(u.min(v), u.max(v)))
            .ok()
            .map(|i| i as u32)
    };

    // Warp vertices that a cut point violates.
    let mut active = vec![true; cut_edges.len()];
    let mut incident: Vec<(u32, u32)> = cut_edges
        .iter()
        .enumerate()
        .flat_map(|(e, &(u, v))| [(u, e as u32), (v, e as u32)])
        .collect();
    incident.sort_unstable();
    let mut warped: BTreeMap<u32, Point3> = BTreeMap::new();
    let mut i = 0;
    while i < incident.len() {
        let v = incident[i].0;
        let j = (i..incident.len())
            .find(|&j| incident[j].0 != v)
            .unwrap_or(incident.len());
        let pv = lat.point(v);
        let mut best: Option<(f64, u32)> = None;
        for &(_, e) in &incident[i..j] {
            if !active[e as usize] {
                continue;
            }
            let (a, b) = cut_edges[e as usize];
            let long = lat.is_red(a) == lat.is_red(b);
            let (length, alpha) = if long {
                (h, params.alpha_long)
            } else {
                (h * 3f64.sqrt() / 2.0, params.alpha_short)
            };
            let d2 = distance_squared(cut_point[e as usize], pv);
            let limit = alpha * length;
            if d2 < limit * limit && best.is_none_or(|(bd, be)| (d2, e) < (bd, be)) {
                best = Some((d2, e));
            }
        }
        if let Some((_, e)) = best {
            warped.insert(v, cut_point[e as usize]);
            sign[v as usize] = 0;
            for &(_, e2) in &incident[i..j] {
                active[e2 as usize] = false;
            }
        }
        i = j;
    }
    let position = |node: Node| match node {
        Node::Lattice(v) => warped.get(&v).copied().unwrap_or_else(|| lat.point(v)),
        Node::Cut(e) => cut_point[e as usize],
    };

    // Stencils.
    let mut out = Stuffed {
        lattice_vertices,
        cut_points: active.iter().filter(|&&a| a).count() as u64,
        warped: warped.len() as u64,
        ..Stuffed::default()
    };
    let mut numbering: HashMap<Node, u32> = HashMap::new();
    let mut emit = |nodes: [Node; 4], out: &mut Stuffed| {
        let p = nodes.map(position);
        let mut t = nodes;
        match orient3d(p[0], p[1], p[2], p[3]) {
            Sign::Positive => {}
            Sign::Negative => t.swap(0, 1),
            Sign::Zero => {
                out.dropped += 1;
                return;
            }
        }
        let ids = t.map(|node| {
            *numbering.entry(node).or_insert_with(|| {
                out.points.push(position(node));
                (out.points.len() - 1) as u32
            })
        });
        out.tets.push(ids);
    };
    let s = |v: u32| sign[v as usize];
    lat.for_each_tet(|t| {
        let neg = t.iter().filter(|&&v| s(v) < 0).count();
        let pos = t.iter().filter(|&&v| s(v) > 0).count();
        if pos == 0 {
            // Whole tets need an inside vertex: a tet whose four vertices were
            // all warped onto the surface is a sliver lying in it, and its
            // faces close the boundary just as well.
            if neg > 0 {
                emit(t.map(Node::Lattice), &mut out);
            }
            return;
        }
        if neg == 0 {
            return;
        }
        let cut = |u: u32, v: u32| -> Option<Node> {
            (s(u) * s(v) < 0)
                .then(|| cut_index(u, v))
                .flatten()
                .map(Node::Cut)
        };
        let apex = t
            .iter()
            .copied()
            .filter(|&v| s(v) < 0)
            .min()
            .unwrap_or(t[0]);
        let mut triangles: Vec<[Node; 3]> = Vec::new();

        // Inside part of each lattice face.
        for omit in 0..4 {
            let face: Vec<u32> = (0..4).filter(|&k| k != omit).map(|k| t[k]).collect();
            let mut poly: Vec<Node> = Vec::with_capacity(4);
            for k in 0..3 {
                let (v, w) = (face[k], face[(k + 1) % 3]);
                if s(v) <= 0 {
                    poly.push(Node::Lattice(v));
                }
                if let Some(c) = cut(v, w) {
                    poly.push(c);
                }
            }
            match poly.len() {
                3 => triangles.push([poly[0], poly[1], poly[2]]),
                4 => {
                    // Two inside vertices: split from the lower-numbered one.
                    let lo = poly
                        .iter()
                        .filter_map(|n| match n {
                            Node::Lattice(v) => Some(*v),
                            Node::Cut(_) => None,
                        })
                        .min()
                        .unwrap_or(apex);
                    let k = poly
                        .iter()
                        .position(|&n| n == Node::Lattice(lo))
                        .unwrap_or(0);
                    let (a, b, c, d) = (
                        poly[k],
                        poly[(k + 1) % 4],
                        poly[(k + 2) % 4],
                        poly[(k + 3) % 4],
                    );
                    triangles.push([a, b, c]);
                    triangles.push([a, c, d]);
                }
                _ => {}
            }
        }

        // The cap: the surface's section through the tet.
        let zeros: Vec<Node> = t
            .iter()
            .filter(|&&v| s(v) == 0)
            .map(|&v| Node::Lattice(v))
            .collect();
        let mut cuts: Vec<(u32, u32, Node)> = Vec::new();
        for (i, j) in TET_EDGES {
            if let Some(c) = cut(t[i], t[j]) {
                let (inside, outside) = if s(t[i]) < 0 {
                    (t[i], t[j])
                } else {
                    (t[j], t[i])
                };
                cuts.push((inside, outside, c));
            }
        }
        if cuts.len() == 4 {
            // Two inside, two outside: order the cut points around the quad.
            let ins: Vec<u32> = t.iter().copied().filter(|&v| s(v) < 0).collect();
            let outs: Vec<u32> = t.iter().copied().filter(|&v| s(v) > 0).collect();
            let at = |a: u32, p: u32| {
                cuts.iter()
                    .find(|c| c.0 == a && c.1 == p)
                    .map_or(Node::Cut(0), |c| c.2)
            };
            let (a, b, p, q) = (ins[0], ins[1], outs[0], outs[1]);
            let quad = [at(a, p), at(a, q), at(b, q), at(b, p)];
            let d02 = distance_squared(position(quad[0]), position(quad[2]));
            let d13 = distance_squared(position(quad[1]), position(quad[3]));
            if d02 <= d13 {
                triangles.push([quad[0], quad[1], quad[2]]);
                triangles.push([quad[0], quad[2], quad[3]]);
            } else {
                triangles.push([quad[0], quad[1], quad[3]]);
                triangles.push([quad[1], quad[2], quad[3]]);
            }
        } else {
            let cap: Vec<Node> = zeros
                .iter()
                .copied()
                .chain(cuts.iter().map(|c| c.2))
                .collect();
            if cap.len() == 3 {
                triangles.push([cap[0], cap[1], cap[2]]);
            }
        }

        // Cone from the apex over every boundary triangle not touching it.
        for tri in triangles {
            if tri.contains(&Node::Lattice(apex)) {
                continue;
            }
            emit([Node::Lattice(apex), tri[0], tri[1], tri[2]], &mut out);
        }
    });
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sphere(c: Point3, r: f64) -> impl Fn(Point3) -> f64 + Sync {
        move |p: Point3| distance_squared(p, c).sqrt() - r
    }

    #[test]
    fn lattice_counts() {
        let lat = Lattice {
            origin: [0.0; 3],
            h: 1.0,
            n: [3, 3, 3],
        };
        let mut count = 0;
        lat.for_each_tet(|_| count += 1);
        // 3 axes × (2 × 3 × 3 cube pairs) × 4 edges.
        assert_eq!(count, 3 * 18 * 4);
        assert_eq!(lat.point(lat.red([1, 1, 1])), [1.5, 1.5, 1.5]);
        assert_eq!(lat.point(lat.black([3, 2, 1])), [3.0, 2.0, 1.0]);
    }

    #[test]
    fn unit_sphere_volume_and_validity() {
        let bounds = Aabb::of_points(&[[-1.0; 3], [1.0; 3]]);
        let out = stuff(
            sphere([0.0; 3], 1.0),
            &bounds,
            &StuffingParams::with_cell_size(0.1),
        )
        .unwrap();
        assert_eq!(out.dropped, 0);
        let volume: f64 = out
            .tets
            .iter()
            .map(|t| {
                let [a, b, c, d] = t.map(|v| out.points[v as usize]);
                mesher_geom::vec3::dot(sub(b, a), mesher_geom::vec3::cross(sub(c, a), sub(d, a)))
                    / 6.0
            })
            .sum();
        let exact = 4.0 / 3.0 * std::f64::consts::PI;
        assert!((volume - exact).abs() / exact < 0.01, "volume {volume}");
        let (lo, hi) = crate::dihedral_range(&out.points, &out.tets).unwrap();
        assert!(
            lo >= DIHEDRAL_MIN_DEGREES && hi <= DIHEDRAL_MAX_DEGREES,
            "{lo} {hi}"
        );
    }

    #[test]
    fn errors() {
        let b = Aabb::of_points(&[[0.0; 3], [1.0; 3]]);
        let f = sphere([0.5; 3], 0.4);
        assert_eq!(
            stuff(&f, &b, &StuffingParams::with_cell_size(0.0)),
            Err(StuffingError::BadCellSize)
        );
        assert_eq!(
            stuff(&f, &Aabb::EMPTY, &StuffingParams::with_cell_size(0.1)),
            Err(StuffingError::EmptyBounds)
        );
        assert!(matches!(
            stuff(&f, &b, &StuffingParams::with_cell_size(1e-4)),
            Err(StuffingError::Budget { .. })
        ));
    }
}
