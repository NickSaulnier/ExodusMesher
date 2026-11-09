//! Surface analysis shared by diagnosis, repair and the surface gate: one
//! function computes every topological and geometric fact about a triangle
//! surface, and each stage turns the facts into its own records.

use std::collections::BTreeMap;

use mesher_geom::vec3::{cross, dot, norm, sub};
use mesher_geom::{
    Bvh, EdgeId, EdgeKind, FastWindingNumber, GeomError, HalfEdge, Point3, SurfaceMesh,
    TriangleSoup, faces_intersect,
};

use super::Tolerances;
use super::parts::centroid;
use super::uf::DisjointSets;

/// An orientation component whose faces disagree on winding.
#[derive(Debug, Clone)]
pub(crate) struct Disagreement {
    /// Faces to flip (the minority of the component).
    pub minority: Vec<u32>,
    pub inconsistent_edges: Vec<EdgeId>,
    /// The component cannot be oriented at all (a Möbius-like strip).
    pub non_orientable: bool,
}

/// Connected boundary edges of one shell, directed as the oriented faces
/// traverse them.
#[derive(Debug, Clone)]
pub(crate) struct BoundaryLoop {
    pub shell: usize,
    /// `([from, to], face)`, in edge order.
    pub edges: Vec<([u32; 2], u32)>,
}

#[derive(Debug)]
pub(crate) struct Analysis {
    pub mesh: SurfaceMesh,
    pub nonmanifold_edges: Vec<EdgeId>,
    /// `(vertex, fans)` for vertices joining more than one fan.
    pub nonmanifold_vertices: Vec<(u32, usize)>,
    pub shell_faces: Vec<Vec<u32>>,
    pub disagreements: Vec<Disagreement>,
    /// Faces as majority-oriented (minorities flipped).
    pub oriented: Vec<[u32; 3]>,
    pub loops: Vec<BoundaryLoop>,
    /// Per shell: signed volume with holes closed, area, and volume / area^1.5.
    pub volume: Vec<f64>,
    pub area: Vec<f64>,
    pub ratio: Vec<f64>,
    pub sheet: Vec<bool>,
    /// Face pairs of one shell that intersect beyond shared features.
    pub self_pairs: Vec<(u32, u32)>,
    /// Intersecting face pairs between two shells, keyed by shell pair.
    pub shell_overlaps: BTreeMap<(usize, usize), Vec<(u32, u32)>>,
    /// Per shell: how many other disjoint shells contain it, and whether its
    /// volume sign is wrong for that depth.
    pub depth: Vec<u64>,
    pub inverted: Vec<bool>,
    /// `(face, min angle in degrees)` below the sliver threshold.
    pub slivers: Vec<(u32, f64)>,
}

/// Analyse a surface. `contact` gives, per face, the vertex ids used to decide
/// which contacts between faces are shared features (pass `tris` itself, or
/// provenance ids after vertices have been split).
pub(crate) fn analyze_surface(
    points: &[Point3],
    tris: &[[u32; 3]],
    contact: &[[u32; 3]],
    tol: &Tolerances,
) -> Result<Analysis, GeomError> {
    let soup = TriangleSoup::new(points, tris)?;
    let mesh = SurfaceMesh::new(&soup);
    let nf = tris.len();

    let nonmanifold_edges: Vec<EdgeId> = mesh
        .edges()
        .filter(|&e| mesh.edge_kind(e) == Some(EdgeKind::NonManifold))
        .collect();
    let nonmanifold_vertices: Vec<(u32, usize)> = (0..mesh.vertex_count())
        .map(|v| (v, mesh.vertex_fan_count(v)))
        .filter(|&(_, fans)| fans > 1)
        .collect();

    // Shells: faces connected through edges.
    let mut sets = DisjointSets::new(nf);
    for e in mesh.edges() {
        for w in mesh.edge_half_edges(e).windows(2) {
            sets.union(w[0] / 3, w[1] / 3);
        }
    }
    let mut shell_index: BTreeMap<u32, usize> = BTreeMap::new();
    let mut shell_of = vec![0usize; nf];
    for (f, s) in shell_of.iter_mut().enumerate() {
        let root = sets.find(f as u32);
        let next = shell_index.len();
        *s = *shell_index.entry(root).or_insert(next);
    }
    let shells = shell_index.len();
    let mut shell_faces: Vec<Vec<u32>> = vec![Vec::new(); shells];
    for f in 0..nf {
        shell_faces[shell_of[f]].push(f as u32);
    }

    // Orientation: propagate across edges with exactly two faces; the minority
    // of each component is to be flipped.
    let origin = |h: u32| mesh.endpoints(HalfEdge(h)).map_or(0, |e| e.0);
    let mut flag: Vec<Option<bool>> = vec![None; nf];
    let mut flip = vec![false; nf];
    let mut disagreements = Vec::new();
    for seed in 0..nf {
        if flag[seed].is_some() {
            continue;
        }
        flag[seed] = Some(false);
        let mut members = vec![seed as u32];
        let mut stack = vec![seed as u32];
        let mut non_orientable = false;
        let mut inconsistent_edges = Vec::new();
        while let Some(f) = stack.pop() {
            let ff = flag[f as usize].unwrap_or(false);
            for c in 0..3 {
                let h = f * 3 + c;
                let Some(e) = mesh.edge(HalfEdge(h)) else {
                    continue;
                };
                let [x, y] = match *mesh.edge_half_edges(e) {
                    [x, y] => [x, y],
                    _ => continue,
                };
                let other = if x == h { y } else { x };
                let same_direction = origin(other) == origin(h);
                if same_direction && h < other {
                    inconsistent_edges.push(e);
                }
                let want = ff ^ same_direction;
                let g = (other / 3) as usize;
                match flag[g] {
                    None => {
                        flag[g] = Some(want);
                        members.push(g as u32);
                        stack.push(g as u32);
                    }
                    Some(have) => non_orientable |= have != want,
                }
            }
        }
        let ones = members
            .iter()
            .filter(|&&f| flag[f as usize] == Some(true))
            .count();
        let minority_flag = ones * 2 <= members.len();
        let mut minority: Vec<u32> = members
            .iter()
            .copied()
            .filter(|&f| flag[f as usize] == Some(minority_flag))
            .collect();
        minority.sort_unstable();
        for &f in &minority {
            flip[f as usize] = true;
        }
        if !minority.is_empty() || non_orientable {
            disagreements.push(Disagreement {
                minority,
                inconsistent_edges,
                non_orientable,
            });
        }
    }
    let oriented: Vec<[u32; 3]> = tris
        .iter()
        .zip(&flip)
        .map(|(&t, &f)| if f { [t[0], t[2], t[1]] } else { t })
        .collect();

    // Boundary loops: boundary edges of one shell joined through shared vertices.
    let boundary: Vec<(usize, [u32; 2], u32)> = mesh
        .edges()
        .filter(|&e| mesh.edge_kind(e) == Some(EdgeKind::Boundary))
        .filter_map(|e| {
            let h = *mesh.edge_half_edges(e).first()?;
            let (i, j) = mesh.endpoints(HalfEdge(h))?;
            let f = h / 3;
            let dir = if flip[f as usize] { [j, i] } else { [i, j] };
            Some((shell_of[f as usize], dir, f))
        })
        .collect();
    let mut joined = DisjointSets::new(boundary.len());
    let mut first_at: BTreeMap<(usize, u32), u32> = BTreeMap::new();
    for (k, &(s, [i, j], _)) in boundary.iter().enumerate() {
        for v in [i, j] {
            match first_at.get(&(s, v)) {
                Some(&other) => joined.union(other, k as u32),
                None => {
                    first_at.insert((s, v), k as u32);
                }
            }
        }
    }
    let mut by_loop: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for k in 0..boundary.len() {
        by_loop.entry(joined.find(k as u32)).or_default().push(k);
    }
    let loops: Vec<BoundaryLoop> = by_loop
        .into_values()
        .map(|ks| BoundaryLoop {
            shell: boundary[ks[0]].0,
            edges: ks.iter().map(|&k| (boundary[k].1, boundary[k].2)).collect(),
        })
        .collect();

    // Per shell: signed volume with holes closed by a fan from each loop's
    // centroid, and area.
    let mut volume = vec![0.0f64; shells];
    let mut area = vec![0.0f64; shells];
    let reference: Vec<Point3> = shell_faces
        .iter()
        .map(|fs| points[tris[fs[0] as usize][0] as usize])
        .collect();
    for (f, t) in oriented.iter().enumerate() {
        let s = shell_of[f];
        let [a, b, c] = t.map(|v| points[v as usize]);
        volume[s] += tet(reference[s], a, b, c);
        area[s] += norm(cross(sub(b, a), sub(c, a))) / 2.0;
    }
    for l in &loops {
        let mut verts: Vec<u32> = l.edges.iter().flat_map(|e| e.0).collect();
        verts.sort_unstable();
        verts.dedup();
        let n = verts.len() as f64;
        let c = [0, 1, 2].map(|k| verts.iter().map(|&v| points[v as usize][k]).sum::<f64>() / n);
        for &([i, j], _) in &l.edges {
            volume[l.shell] += tet(
                reference[l.shell],
                c,
                points[j as usize],
                points[i as usize],
            );
        }
    }
    let ratio: Vec<f64> = (0..shells)
        .map(|s| {
            if area[s] > 0.0 {
                volume[s].abs() / libm::pow(area[s], 1.5)
            } else {
                0.0
            }
        })
        .collect();
    let sheet: Vec<bool> = ratio
        .iter()
        .map(|&r| r <= tol.open_sheet_volume_ratio)
        .collect();

    // Intersections beyond shared features (by `contact` ids).
    let bvh = Bvh::new(soup);
    let mut self_pairs = Vec::new();
    let mut shell_overlaps: BTreeMap<(usize, usize), Vec<(u32, u32)>> = BTreeMap::new();
    for (f, g) in bvh.self_overlaps() {
        let pos = |t: [u32; 3]| t.map(|v| points[v as usize]);
        let (tf, tg) = (tris[f as usize], tris[g as usize]);
        if !faces_intersect(contact[f as usize], pos(tf), contact[g as usize], pos(tg)) {
            continue;
        }
        let (sf, sg) = (shell_of[f as usize], shell_of[g as usize]);
        if sf == sg {
            self_pairs.push((f, g));
        } else {
            shell_overlaps
                .entry((sf.min(sg), sf.max(sg)))
                .or_default()
                .push((f, g));
        }
    }

    // Nesting and inversion: outward shells need positive volume, cavities
    // (inside an odd number of other disjoint shells) negative.
    let solid: Vec<usize> = (0..shells).filter(|&s| !sheet[s]).collect();
    let shell_tris: Vec<Vec<[u32; 3]>> = shell_faces
        .iter()
        .map(|fs| fs.iter().map(|&f| oriented[f as usize]).collect())
        .collect();
    let mut depth = vec![0u64; shells];
    let mut inverted = vec![false; shells];
    for &s in &solid {
        if solid.len() > 1 {
            let probe = centroid(points, shell_tris[s][0]);
            for &t in &solid {
                if t == s || shell_overlaps.contains_key(&(s.min(t), s.max(t))) {
                    continue;
                }
                let soup = TriangleSoup::new(points, &shell_tris[t])?;
                if FastWindingNumber::new(soup).winding_number(probe).abs() > 0.5 {
                    depth[s] += 1;
                }
            }
        }
        inverted[s] = (volume[s] > 0.0) != (depth[s] % 2 == 0);
    }

    let slivers: Vec<(u32, f64)> = tris
        .iter()
        .enumerate()
        .map(|(f, t)| (f as u32, min_angle_degrees(t.map(|v| points[v as usize]))))
        .filter(|&(_, a)| a < tol.sliver_angle_degrees)
        .collect();

    Ok(Analysis {
        mesh,
        nonmanifold_edges,
        nonmanifold_vertices,
        shell_faces,
        disagreements,
        oriented,
        loops,
        volume,
        area,
        ratio,
        sheet,
        self_pairs,
        shell_overlaps,
        depth,
        inverted,
        slivers,
    })
}

pub(crate) fn tet(r: Point3, a: Point3, b: Point3, c: Point3) -> f64 {
    dot(sub(a, r), cross(sub(b, r), sub(c, r))) / 6.0
}

fn min_angle_degrees(p: [Point3; 3]) -> f64 {
    (0..3)
        .map(|i| {
            let u = sub(p[(i + 1) % 3], p[i]);
            let v = sub(p[(i + 2) % 3], p[i]);
            libm::atan2(norm(cross(u, v)), dot(u, v))
        })
        .fold(f64::INFINITY, f64::min)
        .to_degrees()
}
