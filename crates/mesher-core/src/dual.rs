//! The vertex-centred median dual of a tet mesh (Milestone 21): one
//! polyhedral cell around each vertex the tets use (one per connected piece
//! of the vertex's star, where the tets around a vertex fall apart, as
//! where two cubes meet at a corner or an edge).
//!
//! Built from the finite-volume literature's description of the median dual
//! (Barth's vertex-centred schemes; polyhedral cells from tetrahedra as
//! Perić and Jasak describe them), not from any implementation. Its points
//! are the tets' centroids, the tet faces' centroids, the boundary edges'
//! midpoints and the boundary vertices, in that order:
//!
//! - The face between the cells of vertices `a` and `b` is the ring around
//!   edge `ab` of the centroids of the tets on it and of the faces between
//!   them, turning about `a → b` by the right-hand rule, so its normal points
//!   from `a`'s cell into `b`'s. For a boundary edge the ring opens at the
//!   surface and closes through the edge's midpoint and the centroids of its
//!   two boundary faces. An edge whose tets form several fans (joined by no
//!   face) has a face per fan.
//! - Each boundary triangle gives each of its vertices' cells a
//!   quadrilateral: the vertex, its two edges' midpoints and the triangle's
//!   centroid, in the triangle's plane and its outward order, in the
//!   triangle's patch. So the dual's boundary is the tet mesh's surface
//!   exactly, with its patches.
//!
//! Cells are numbered by block, vertex and first tet, so each block is a zone;
//! faces are in [`FaceMesh`]'s order. Deterministic. The construction is
//! combinatorial: which way a ring turns follows from each tet's positive
//! orientation, so no floating-point decision is made; [`check`] then
//! tests the result's geometry as a solver computes it.

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{add, cross, dot, norm, scale, sub};
use mesher_geom::{Point3, Sign};

use crate::mesh::{CellShape, FaceMesh, MeshProblem, Patch, Polygons, VolumeMesh, Zone};

/// A tet mesh's dual: its points, its faces (cells as polyhedra), and each
/// cell's vertex in the tet mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Dual {
    pub points: Vec<Point3>,
    pub faces: FaceMesh,
    /// Per cell, the tet mesh's vertex it surrounds.
    pub vertices: Vec<u32>,
}

/// The median dual of `mesh`.
///
/// # Errors
/// The mesh's [`MeshProblem`] when it has no face view, when a vertex is in
/// tets of two blocks ([`MeshProblem::SharedVertex`]), or when the tets
/// around an edge do not form fans, open ones between two boundary faces or
/// closed ones ([`MeshProblem::NonManifoldEdge`]); a mesh the volume gate
/// accepts has none of these.
pub fn dual(mesh: &VolumeMesh) -> Result<Dual, MeshProblem> {
    let tf = mesh.faces()?;
    let tets: Vec<[u32; 4]> = mesh.elements().collect();
    let np = mesh.points.len();

    let mut block_of = vec![u32::MAX; np];
    for (b, block) in mesh.blocks.iter().enumerate() {
        for t in &block.tets {
            for &v in t {
                let Some(slot) = block_of.get_mut(v as usize) else {
                    return Err(MeshProblem::NodeOutOfRange {
                        element: 0,
                        node: v,
                    });
                };
                if *slot == u32::MAX {
                    *slot = b as u32;
                } else if *slot != b as u32 {
                    return Err(MeshProblem::SharedVertex { node: v });
                }
            }
        }
    }

    // Each tet's faces in the face view, by the local vertex opposite.
    let mut tet_face = vec![[u32::MAX; 4]; tets.len()];
    for (k, face) in tf.faces.iter().enumerate() {
        for t in [Some(tf.owner[k]), tf.neighbour.get(k).copied()]
            .into_iter()
            .flatten()
        {
            let tet = tets[t as usize];
            if let Some(i) = (0..4).find(|&i| !face.contains(&tet[i])) {
                tet_face[t as usize][i] = k as u32;
            }
        }
    }
    if let Some(t) = tet_face.iter().position(|f| f.contains(&u32::MAX)) {
        return Err(MeshProblem::BadFaceView {
            face: t as u32,
            problem: "a tet's face missing from the face view",
        });
    }
    let boundary = |k: u32| k as usize >= tf.internal();
    let local =
        |t: u32, v: u32| -> usize { tets[t as usize].iter().position(|&x| x == v).unwrap_or(0) };

    // Cells: the pieces of each vertex's star, its tets joined through the
    // faces they share at the vertex; a piece is a tet corner `4 t + i`.
    let mut star = UnionFind::new(4 * tets.len());
    for k in 0..tf.internal() {
        let (o, n) = (tf.owner[k], tf.neighbour[k]);
        for &u in tf.faces.get(k) {
            star.union(4 * o as usize + local(o, u), 4 * n as usize + local(n, u));
        }
    }
    let mut corners: Vec<(u32, u32, u32, usize)> = Vec::with_capacity(4 * tets.len());
    for (t, tet) in tets.iter().enumerate() {
        for (i, &v) in tet.iter().enumerate() {
            corners.push((block_of[v as usize], v, t as u32, i));
        }
    }
    corners.sort_unstable();
    let mut cell_of_root = vec![u32::MAX; 4 * tets.len()];
    let mut vertices: Vec<u32> = Vec::new();
    let mut zone_counts = vec![0usize; mesh.blocks.len()];
    for &(b, v, t, i) in &corners {
        let root = star.find(4 * t as usize + i);
        if cell_of_root[root] == u32::MAX {
            cell_of_root[root] = vertices.len() as u32;
            vertices.push(v);
            zone_counts[b as usize] += 1;
        }
    }
    let cell = |star: &mut UnionFind, t: u32, v: u32| -> u32 {
        cell_of_root[star.find(4 * t as usize + local(t, v))]
    };

    // Points: tet centroids, face centroids, boundary edge midpoints,
    // boundary vertices.
    let p = |v: u32| mesh.points[v as usize];
    let mut boundary_vertices: Vec<u32> = Vec::new();
    let mut boundary_edges: Vec<(u32, u32)> = Vec::new();
    for k in tf.internal()..tf.faces.len() {
        let face = tf.faces.get(k);
        for (i, &a) in face.iter().enumerate() {
            let b = face[(i + 1) % face.len()];
            boundary_vertices.push(a);
            boundary_edges.push((a.min(b), a.max(b)));
        }
    }
    boundary_vertices.sort_unstable();
    boundary_vertices.dedup();
    boundary_edges.sort_unstable();
    boundary_edges.dedup();
    let tet_point = |t: u32| t;
    let face_point = |k: u32| tets.len() as u32 + k;
    let edges_from = (tets.len() + tf.faces.len()) as u32;
    let vertices_from = edges_from + boundary_edges.len() as u32;
    let edge_point = |a: u32, b: u32| -> u32 {
        let key = (a.min(b), a.max(b));
        edges_from + boundary_edges.binary_search(&key).unwrap_or(0) as u32
    };
    let vertex_point =
        |v: u32| -> u32 { vertices_from + boundary_vertices.binary_search(&v).unwrap_or(0) as u32 };
    let mean = |nodes: &[u32]| -> Point3 {
        let sum = nodes.iter().fold([0.0; 3], |s, &v| add(s, p(v)));
        scale(sum, 1.0 / nodes.len() as f64)
    };
    let mut points = Vec::with_capacity(vertices_from as usize + boundary_vertices.len());
    points.extend(tets.iter().map(|t| mean(t)));
    points.extend(tf.faces.iter().map(mean));
    points.extend(
        boundary_edges
            .iter()
            .map(|&(a, b)| scale(add(p(a), p(b)), 0.5)),
    );
    points.extend(boundary_vertices.iter().map(|&v| p(v)));

    // Each tet's six edges, `a < b`, with the faces a ring around `a → b`
    // enters and leaves the tet by. In a positive tet (a, b, c, d) (an even
    // permutation of its nodes), turning about a → b by the right-hand rule
    // goes from face abc (opposite d) to face abd (opposite c).
    let mut around: Vec<(u32, u32, u32, u32, u32)> = Vec::with_capacity(6 * tets.len());
    for (t, tet) in tets.iter().enumerate() {
        for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            let (i, j) = if tet[i] < tet[j] { (i, j) } else { (j, i) };
            let mut rest = (0..4).filter(|&m| m != i && m != j);
            let (k, l) = (rest.next().unwrap_or(0), rest.next().unwrap_or(0));
            let (enter, leave) = if even([i, j, k, l]) {
                (tet_face[t][l], tet_face[t][k])
            } else {
                (tet_face[t][k], tet_face[t][l])
            };
            around.push((tet[i], tet[j], t as u32, enter, leave));
        }
    }
    around.sort_unstable();

    let mut internal: Vec<(u32, u32, Vec<u32>)> = Vec::new();
    let mut s = 0;
    while s < around.len() {
        let (a, b) = (around[s].0, around[s].1);
        let e = (s..around.len())
            .find(|&e| (around[e].0, around[e].1) != (a, b))
            .unwrap_or(around.len());
        let fan = &around[s..e];
        let bad = || MeshProblem::NonManifoldEdge { a, b };
        let mut used = vec![false; fan.len()];
        // Open fans first (from each boundary face they start at), then the
        // closed ones, each from its first tet.
        let starts = (0..fan.len())
            .filter(|&i| boundary(fan[i].3))
            .chain(0..fan.len());
        for first in starts {
            if used[first] {
                continue;
            }
            let open = boundary(fan[first].3);
            let mut ring = Vec::with_capacity(2 * fan.len() + 2);
            if open {
                ring.push(edge_point(a, b));
                ring.push(face_point(fan[first].3));
            }
            let mut at = first;
            loop {
                if used[at] {
                    return Err(bad());
                }
                used[at] = true;
                let (_, _, t, _, leave) = fan[at];
                ring.push(tet_point(t));
                ring.push(face_point(leave));
                if boundary(leave) {
                    if open {
                        break;
                    }
                    return Err(bad());
                }
                match fan.iter().position(|x| x.3 == leave) {
                    Some(next) if next == first && !open => break,
                    Some(next) => at = next,
                    None => return Err(bad()),
                }
            }
            let t = fan[first].2;
            internal.push((cell(&mut star, t, a), cell(&mut star, t, b), ring));
        }
        s = e;
    }
    internal.sort_unstable_by_key(|&(o, n, _)| (o, n));

    let quads: usize = 3 * (tf.faces.len() - tf.internal());
    let mut faces = FaceMesh {
        faces: Polygons::with_capacity(internal.len() + quads, 0),
        owner: Vec::with_capacity(internal.len() + quads),
        neighbour: Vec::with_capacity(internal.len()),
        patches: Vec::new(),
        zones: Vec::new(),
        cells: vertices.len(),
        shape: CellShape::Polyhedron,
    };
    for (o, n, ring) in &internal {
        faces.faces.push(ring);
        faces.owner.push(*o);
        faces.neighbour.push(*n);
    }
    for patch in &tf.patches {
        let mut quads: Vec<(u32, [u32; 4])> = Vec::with_capacity(3 * patch.count);
        for k in patch.start..patch.start + patch.count {
            let tri = tf.faces.get(k);
            let owner = tf.owner[k];
            for i in 0..tri.len() {
                let (v, next, prev) = (tri[i], tri[(i + 1) % 3], tri[(i + 2) % 3]);
                quads.push((
                    cell(&mut star, owner, v),
                    [
                        vertex_point(v),
                        edge_point(v, next),
                        face_point(k as u32),
                        edge_point(prev, v),
                    ],
                ));
            }
        }
        quads.sort_by_key(|&(c, _)| c);
        faces.patches.push(Patch {
            name: patch.name.clone(),
            start: faces.faces.len(),
            count: quads.len(),
        });
        for (c, q) in quads {
            faces.faces.push(&q);
            faces.owner.push(c);
        }
    }
    let mut start = 0;
    for (block, count) in mesh.blocks.iter().zip(zone_counts) {
        faces.zones.push(Zone {
            name: block.name.clone(),
            start,
            count,
        });
        start += count;
    }
    faces.check(points.len())?;
    Ok(Dual {
        points,
        faces,
        vertices,
    })
}

/// Disjoint sets with path halving, the smaller root kept.
struct UnionFind(Vec<usize>);

impl UnionFind {
    fn new(n: usize) -> UnionFind {
        UnionFind((0..n).collect())
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        let (lo, hi) = (a.min(b), a.max(b));
        self.0[hi] = lo;
    }
}

/// Whether a permutation of 0..4 is even.
fn even(p: [usize; 4]) -> bool {
    let mut inversions = 0;
    for i in 0..4 {
        for j in i + 1..4 {
            if p[i] > p[j] {
                inversions += 1;
            }
        }
    }
    inversions % 2 == 0
}

/// A polyhedral mesh's geometry, computed as a cell-centred finite-volume
/// solver computes it: each face split into triangles about its nodes'
/// average; its area vector the triangles' sum, its centre their
/// area-weighted centroid; each cell's centre and volume from the pyramids
/// its faces make with a first estimate of its centre (the average of its
/// faces' centres).
#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    pub face_centres: Vec<Point3>,
    pub face_areas: Vec<Point3>,
    pub cell_centres: Vec<Point3>,
    pub cell_volumes: Vec<f64>,
}

/// The geometry of the cells of `faces` on `points` (a layout that passes
/// [`FaceMesh::check`]).
#[must_use]
pub fn geometry(points: &[Point3], faces: &FaceMesh) -> Geometry {
    let n = faces.faces.len();
    let mut face_centres = Vec::with_capacity(n);
    let mut face_areas = Vec::with_capacity(n);
    for face in faces.faces.iter() {
        let (centre, area) = polygon(points, face);
        face_centres.push(centre);
        face_areas.push(area);
    }
    let cells = faces.cells;
    let mut estimate = vec![[0.0; 3]; cells];
    let mut count = vec![0u32; cells];
    for (k, &centre) in face_centres.iter().enumerate() {
        for c in sides(faces, k).map(|(c, _)| c) {
            estimate[c] = add(estimate[c], centre);
            count[c] += 1;
        }
    }
    for (e, &m) in estimate.iter_mut().zip(&count) {
        *e = scale(*e, 1.0 / f64::from(m.max(1)));
    }
    let mut volume = vec![0.0; cells];
    let mut moment = vec![[0.0; 3]; cells];
    for k in 0..n {
        for (c, out) in sides(faces, k) {
            let v = out * dot(face_areas[k], sub(face_centres[k], estimate[c])) / 3.0;
            let centroid = add(scale(face_centres[k], 0.75), scale(estimate[c], 0.25));
            volume[c] += v;
            moment[c] = add(moment[c], scale(centroid, v));
        }
    }
    let cell_centres = (0..cells)
        .map(|c| {
            if volume[c].abs() > 0.0 {
                scale(moment[c], 1.0 / volume[c])
            } else {
                estimate[c]
            }
        })
        .collect();
    Geometry {
        face_centres,
        face_areas,
        cell_centres,
        cell_volumes: volume,
    }
}

/// A face's centre and area vector, as [`Geometry`] computes them.
fn polygon(points: &[Point3], face: &[u32]) -> (Point3, Point3) {
    let a = average(points, face);
    let mut area = [0.0; 3];
    let mut weighted = [0.0; 3];
    let mut total = 0.0;
    for (i, &u) in face.iter().enumerate() {
        let w = face[(i + 1) % face.len()];
        let (pu, pw) = (points[u as usize], points[w as usize]);
        let t = scale(cross(sub(pu, a), sub(pw, a)), 0.5);
        let m = norm(t);
        area = add(area, t);
        weighted = add(weighted, scale(add(add(a, pu), pw), m / 3.0));
        total += m;
    }
    let centre = if total > 0.0 {
        scale(weighted, 1.0 / total)
    } else {
        a
    };
    (centre, area)
}

/// The cells on face `k`: its owner (the face out of it, sign 1) and its
/// neighbour (sign -1), if any.
fn sides(faces: &FaceMesh, k: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
    [
        Some((faces.owner[k] as usize, 1.0)),
        faces.neighbour.get(k).map(|&n| (n as usize, -1.0)),
    ]
    .into_iter()
    .flatten()
}

fn average(points: &[Point3], nodes: &[u32]) -> Point3 {
    let sum = nodes
        .iter()
        .fold([0.0; 3], |s, &v| add(s, points[v as usize]));
    scale(sum, 1.0 / nodes.len().max(1) as f64)
}

/// At most this many cells are listed per finding of [`check`].
pub const LISTED_CELLS: usize = 20;

/// What [`check`] finds in a polyhedral mesh.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CellCheck {
    pub cells: usize,
    /// Cells whose faces do not close: some edge of theirs is not on
    /// exactly two of their faces, once each way (the first
    /// [`LISTED_CELLS`], and how many).
    pub open: (Vec<u32>, usize),
    /// Cells with a face triangle (about the face's average point) that is
    /// not strictly on the outer side of the cell's centre, by exact
    /// orientation on the floating-point points: the decomposition a solver
    /// uses has a flat or inverted tet. A measure of warped faces, not of
    /// validity: the median dual's rings are not planar, and some of their
    /// triangles lean past the centre even in the duals of clean meshes
    /// that `checkMesh` passes.
    pub inverted_triangles: (Vec<u32>, usize),
    /// Cells with a face whose pyramid on the cell's centre (the face's area
    /// vector against the line from the centre to the face's centre) is not
    /// positive: `checkMesh`'s face pyramids check.
    pub inverted_pyramids: (Vec<u32>, usize),
    /// The cells' volumes summed, and the smallest.
    pub volume: f64,
    pub min_volume: f64,
}

impl CellCheck {
    /// Whether every cell is closed and every face pyramid positive (what
    /// `checkMesh` requires; [`CellCheck::inverted_triangles`] is not
    /// counted).
    #[must_use]
    pub fn valid(&self) -> bool {
        self.open.1 == 0 && self.inverted_pyramids.1 == 0
    }
}

fn listed(mut cells: Vec<u32>) -> (Vec<u32>, usize) {
    cells.sort_unstable();
    cells.dedup();
    let n = cells.len();
    cells.truncate(LISTED_CELLS);
    (cells, n)
}

/// Check the cells of `faces` on `points` (a layout that passes
/// [`FaceMesh::check`]) with their `geometry`: closed, and every face's
/// triangles and pyramid on the outer side of the cell's centre.
#[must_use]
pub fn check(points: &[Point3], faces: &FaceMesh, geometry: &Geometry) -> CellCheck {
    // Closure: each cell's edges, directed out of the cell, pair up.
    let mut edges: Vec<(u32, u32, u32, bool)> = Vec::with_capacity(2 * faces.faces.nodes().len());
    let mut triangles = Vec::new();
    let mut pyramids = Vec::new();
    for (k, face) in faces.faces.iter().enumerate() {
        let a = average(points, face);
        for (c, out) in sides(faces, k) {
            let cell = c as u32;
            let centre = geometry.cell_centres[c];
            for (i, &u) in face.iter().enumerate() {
                let w = face[(i + 1) % face.len()];
                // The edge as the cell sees it, out of the cell.
                let (from, to) = if out > 0.0 { (u, w) } else { (w, u) };
                edges.push((cell, from.min(to), from.max(to), from < to));
                // Outward triangle (a, from, to): the centre must be below.
                let want = Sign::Negative;
                let got = orient3d(a, points[from as usize], points[to as usize], centre);
                if got != want {
                    triangles.push(cell);
                }
            }
            let pyramid = out
                * dot(
                    geometry.face_areas[k],
                    sub(geometry.face_centres[k], centre),
                );
            if pyramid <= 0.0 || !pyramid.is_finite() {
                pyramids.push(cell);
            }
        }
    }
    edges.sort_unstable();
    let mut open = Vec::new();
    let mut i = 0;
    while i < edges.len() {
        let key = (edges[i].0, edges[i].1, edges[i].2);
        let j = (i..edges.len())
            .find(|&j| (edges[j].0, edges[j].1, edges[j].2) != key)
            .unwrap_or(edges.len());
        if j - i != 2 || edges[i].3 == edges[i + 1].3 {
            open.push(key.0);
        }
        i = j;
    }
    let volume = geometry.cell_volumes.iter().sum();
    let min_volume = geometry
        .cell_volumes
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    CellCheck {
        cells: faces.cells,
        open: listed(open),
        inverted_triangles: listed(triangles),
        inverted_pyramids: listed(pyramids),
        volume,
        min_volume,
    }
}

/// The most rounds [`merge`] runs.
pub const MERGE_ROUNDS: usize = 8;

/// What [`merge`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct Merged {
    /// The dual with its invalid cells merged; each cell's vertex the
    /// lowest of the cells merged into it.
    pub dual: Dual,
    /// Merges made (each joins two cells, so the dual has this many fewer).
    pub merges: usize,
    /// Rounds run.
    pub rounds: usize,
}

/// A face as merging handles it: its nodes, listed out of cell `a`; `b`
/// the cell on its other side, or `None` on the boundary, in patch `patch`.
#[derive(Debug, Clone, PartialEq)]
struct Side {
    nodes: Vec<u32>,
    a: u32,
    b: Option<u32>,
    patch: usize,
}

/// Merge the cells of `dual` whose face pyramids are not positive into
/// neighbours, so that their centres move inside them: in rounds of at
/// most [`MERGE_ROUNDS`], each invalid cell (worst first) joins the
/// neighbour across one of its faces that leaves the best merged cell, if
/// that is better than the cell was; no two merges of a round touch the same
/// cell or each other's neighbours, so each merged cell is the one weighed.
/// The faces a merged cell has with one other cell become one face (their
/// shared edges dropped); a merge whose faces would not form one loop is
/// not made, so no two faces join the same two cells. Merging stops when
/// every cell is valid or no merge helps; cells that stay invalid are left
/// for [`check`] to find. Deterministic.
#[must_use]
pub fn merge(dual: &Dual) -> Merged {
    let points = &dual.points;
    let patch_names: Vec<String> = dual.faces.patches.iter().map(|p| p.name.clone()).collect();
    let mut zone_of: Vec<u32> = Vec::with_capacity(dual.faces.cells);
    let zones: Vec<String> = dual.faces.zones.iter().map(|z| z.name.clone()).collect();
    for (z, zone) in dual.faces.zones.iter().enumerate() {
        zone_of.extend(std::iter::repeat_n(z as u32, zone.count));
    }
    let mut vertices = dual.vertices.clone();
    let mut sides = sides_of(&dual.faces);
    let mut merges = 0;
    let mut rounds = 0;
    let mut faces = dual.faces.clone();
    while rounds < MERGE_ROUNDS {
        let g = geometry(points, &faces);
        let score = scores(&faces, &g);
        let mut bad: Vec<u32> = (0..faces.cells as u32)
            .filter(|&c| score[c as usize] <= 0.0)
            .collect();
        if bad.is_empty() {
            break;
        }
        rounds += 1;
        bad.sort_by(|&x, &y| {
            score[x as usize]
                .total_cmp(&score[y as usize])
                .then(x.cmp(&y))
        });
        let mut of_cell: Vec<Vec<usize>> = vec![Vec::new(); faces.cells];
        for (k, s) in sides.iter().enumerate() {
            of_cell[s.a as usize].push(k);
            if let Some(b) = s.b {
                of_cell[b as usize].push(k);
            }
        }
        let neighbours = |c: u32| -> Vec<u32> {
            let mut n: Vec<u32> = of_cell[c as usize]
                .iter()
                .filter_map(|&k| {
                    let s = &sides[k];
                    s.b.map(|b| if s.a == c { b } else { s.a })
                })
                .collect();
            n.sort_unstable();
            n.dedup();
            n
        };
        let mut locked = vec![false; faces.cells];
        let mut chosen: Vec<(u32, u32, Vec<Side>)> = Vec::new();
        for &c in &bad {
            if locked[c as usize] {
                continue;
            }
            let mut best: Option<(f64, u32, Vec<Side>)> = None;
            for n in neighbours(c) {
                if locked[n as usize] || zone_of[n as usize] != zone_of[c as usize] {
                    continue;
                }
                let Some(merged) = merged_sides(&sides, &of_cell, c, n) else {
                    continue;
                };
                let s = cell_score(points, &merged, c);
                if best.as_ref().is_none_or(|(b, _, _)| s > *b) {
                    best = Some((s, n, merged));
                }
            }
            if let Some((_, n, merged)) = best.filter(|(s, _, _)| *s > score[c as usize]) {
                for m in neighbours(c).into_iter().chain(neighbours(n)).chain([c, n]) {
                    locked[m as usize] = true;
                }
                chosen.push((c, n, merged));
            }
        }
        if chosen.is_empty() {
            break;
        }
        // Apply: the merged cells' sides replace theirs; then renumber.
        let mut drop = vec![false; sides.len()];
        let mut into: Vec<u32> = (0..faces.cells as u32).collect();
        let mut added: Vec<Side> = Vec::new();
        for (c, n, merged) in chosen {
            for &k in of_cell[c as usize].iter().chain(&of_cell[n as usize]) {
                drop[k] = true;
            }
            into[n as usize] = c;
            added.extend(merged);
            merges += 1;
        }
        let mut kept: Vec<Side> = sides
            .into_iter()
            .zip(drop)
            .filter_map(|(s, d)| (!d).then_some(s))
            .collect();
        kept.extend(added);
        // New numbers: cells in order, merged ones gone; each keeps the
        // lowest vertex of the cells merged into it.
        let mut lowest = vertices.clone();
        for c in 0..faces.cells {
            let root = into[c] as usize;
            lowest[root] = lowest[root].min(vertices[c]);
        }
        let mut number = vec![u32::MAX; faces.cells];
        let mut next = 0u32;
        let mut new_vertices = Vec::with_capacity(faces.cells);
        let mut new_zone = Vec::with_capacity(faces.cells);
        for c in 0..faces.cells {
            if into[c] as usize == c {
                number[c] = next;
                next += 1;
                new_vertices.push(lowest[c]);
                new_zone.push(zone_of[c]);
            }
        }
        for c in 0..faces.cells {
            number[c] = number[into[c] as usize];
        }
        for s in &mut kept {
            s.a = number[s.a as usize];
            s.b = s.b.map(|b| number[b as usize]);
        }
        sides = kept;
        vertices = new_vertices;
        zone_of = new_zone;
        faces = build(sides.clone(), &vertices, &zone_of, &zones, &patch_names);
    }
    Merged {
        dual: Dual {
            points: dual.points.clone(),
            faces,
            vertices,
        },
        merges,
        rounds,
    }
}

/// Each cell's worst face: the smallest of its faces' pyramid heights (the
/// distance from the cell's centre to the face's plane through its centre,
/// out of the cell), over the cell's size (the cube root of its volume);
/// not positive for an invalid cell.
fn scores(faces: &FaceMesh, g: &Geometry) -> Vec<f64> {
    let mut worst = vec![f64::INFINITY; faces.cells];
    for k in 0..faces.faces.len() {
        for (c, out) in sides(faces, k) {
            let s = g.face_areas[k];
            let h = out * dot(s, sub(g.face_centres[k], g.cell_centres[c]))
                / norm(s).max(f64::MIN_POSITIVE);
            let size = libm::cbrt(g.cell_volumes[c].max(0.0)).max(f64::MIN_POSITIVE);
            let r = if g.cell_volumes[c] > 0.0 {
                h / size
            } else {
                -1.0
            };
            worst[c] = worst[c].min(if r.is_finite() { r } else { -1.0 });
        }
    }
    worst
}

/// The score of a cell `c` made of the faces `sides` (all with `a == c`),
/// as [`scores`] gives it.
fn cell_score(points: &[Point3], sides: &[Side], c: u32) -> f64 {
    let geo: Vec<(Point3, Point3)> = sides.iter().map(|s| polygon(points, &s.nodes)).collect();
    let n = geo.len().max(1) as f64;
    let estimate = scale(geo.iter().fold([0.0; 3], |e, (cf, _)| add(e, *cf)), 1.0 / n);
    let mut volume = 0.0;
    let mut moment = [0.0; 3];
    for (s, (cf, area)) in sides.iter().zip(&geo) {
        debug_assert_eq!(s.a, c);
        let v = dot(*area, sub(*cf, estimate)) / 3.0;
        volume += v;
        moment = add(
            moment,
            scale(add(scale(*cf, 0.75), scale(estimate, 0.25)), v),
        );
    }
    if volume <= 0.0 || !volume.is_finite() {
        return -1.0;
    }
    let centre = scale(moment, 1.0 / volume);
    let size = libm::cbrt(volume);
    geo.iter()
        .map(|(cf, area)| dot(*area, sub(*cf, centre)) / norm(*area).max(f64::MIN_POSITIVE) / size)
        .fold(f64::INFINITY, f64::min)
}

/// The faces of `faces` as sides: each internal face once, out of its
/// owner; each boundary face with its patch.
fn sides_of(faces: &FaceMesh) -> Vec<Side> {
    let mut patch_of = vec![0usize; faces.faces.len()];
    for (p, patch) in faces.patches.iter().enumerate() {
        for slot in &mut patch_of[patch.start..patch.start + patch.count] {
            *slot = p;
        }
    }
    faces
        .faces
        .iter()
        .enumerate()
        .map(|(k, nodes)| Side {
            nodes: nodes.to_vec(),
            a: faces.owner[k],
            b: faces.neighbour.get(k).copied(),
            patch: patch_of[k],
        })
        .collect()
}

/// The faces of the cell `c` and `n` would make, all listed out of it
/// (`a == c`): their faces with each other dropped, and their faces with
/// each other cell joined into one; `None` if some such faces do not join
/// into one loop.
fn merged_sides(sides: &[Side], of_cell: &[Vec<usize>], c: u32, n: u32) -> Option<Vec<Side>> {
    // Every face of either, out of it; grouped by the cell beyond.
    let mut beyond: Vec<(Option<u32>, usize, Vec<u32>)> = Vec::new();
    for &k in of_cell[c as usize].iter().chain(&of_cell[n as usize]) {
        let s = &sides[k];
        let inside = |x: u32| x == c || x == n;
        if inside(s.a) && s.b.is_some_and(inside) {
            continue;
        }
        if inside(s.a) {
            beyond.push((s.b, s.patch, s.nodes.clone()));
        } else {
            let mut nodes = s.nodes.clone();
            nodes.reverse();
            beyond.push((Some(s.a), s.patch, nodes));
        }
    }
    beyond.sort_by_key(|(b, p, _)| (b.map_or((1, 0), |b| (0, b)), *p));
    let mut out = Vec::new();
    let mut i = 0;
    while i < beyond.len() {
        // Boundary faces stay apart; faces with one other cell are joined.
        let j = if beyond[i].0.is_none() {
            i + 1
        } else {
            (i + 1..beyond.len())
                .find(|&j| beyond[j].0 != beyond[i].0)
                .unwrap_or(beyond.len())
        };
        let group: Vec<&Vec<u32>> = beyond[i..j].iter().map(|(_, _, nodes)| nodes).collect();
        let nodes = if group.len() == 1 {
            group[0].clone()
        } else {
            join(&group)?
        };
        out.push(Side {
            nodes,
            a: c,
            b: beyond[i].0,
            patch: beyond[i].1,
        });
        i = j;
    }
    Some(out)
}

/// Polygons with the same orientation joined into one: the edges two of
/// them share (once each way) dropped, the rest followed as one loop from
/// the first face's first kept edge; `None` unless they form exactly one
/// loop through distinct nodes.
fn join(faces: &[&Vec<u32>]) -> Option<Vec<u32>> {
    let mut edges: Vec<(u32, u32)> = Vec::new();
    for f in faces {
        for (i, &u) in f.iter().enumerate() {
            edges.push((u, f[(i + 1) % f.len()]));
        }
    }
    let mut sorted = edges.clone();
    sorted.sort_unstable();
    let kept: Vec<(u32, u32)> = edges
        .into_iter()
        .filter(|&(u, w)| sorted.binary_search(&(w, u)).is_err())
        .collect();
    let first = *kept.first()?;
    let mut by_start = kept.clone();
    by_start.sort_unstable();
    if by_start.windows(2).any(|w| w[0].0 == w[1].0) {
        return None;
    }
    let mut loop_nodes = vec![first.0];
    let mut at = first.1;
    while at != first.0 {
        if loop_nodes.len() > kept.len() {
            return None;
        }
        loop_nodes.push(at);
        let k = by_start.binary_search_by_key(&at, |e| e.0).ok()?;
        at = by_start[k].1;
    }
    (loop_nodes.len() == kept.len() && loop_nodes.len() >= 3).then_some(loop_nodes)
}

/// A face view from sides: internal faces out of the lower cell, in owner
/// and neighbour order; boundary faces by patch, then cell.
fn build(
    sides: Vec<Side>,
    vertices: &[u32],
    zone_of: &[u32],
    zones: &[String],
    patches: &[String],
) -> FaceMesh {
    let mut internal: Vec<(u32, u32, Vec<u32>)> = Vec::new();
    let mut boundary: Vec<(usize, u32, Vec<u32>)> = Vec::new();
    for s in sides {
        match s.b {
            Some(b) if b < s.a => {
                let mut nodes = s.nodes;
                nodes.reverse();
                internal.push((b, s.a, nodes));
            }
            Some(b) => internal.push((s.a, b, s.nodes)),
            None => boundary.push((s.patch, s.a, s.nodes)),
        }
    }
    internal.sort_by_key(|&(o, n, _)| (o, n));
    boundary.sort_by_key(|&(p, c, _)| (p, c));
    let mut faces = FaceMesh {
        faces: Polygons::default(),
        owner: Vec::new(),
        neighbour: Vec::new(),
        patches: Vec::new(),
        zones: Vec::new(),
        cells: vertices.len(),
        shape: CellShape::Polyhedron,
    };
    for (o, n, nodes) in internal {
        faces.faces.push(&nodes);
        faces.owner.push(o);
        faces.neighbour.push(n);
    }
    let mut b = 0;
    for (p, name) in patches.iter().enumerate() {
        let start = faces.faces.len();
        while b < boundary.len() && boundary[b].0 == p {
            faces.faces.push(&boundary[b].2);
            faces.owner.push(boundary[b].1);
            b += 1;
        }
        if faces.faces.len() > start {
            faces.patches.push(Patch {
                name: name.clone(),
                start,
                count: faces.faces.len() - start,
            });
        }
    }
    let mut start = 0;
    for (z, name) in zones.iter().enumerate() {
        let count = zone_of.iter().filter(|&&x| x as usize == z).count();
        faces.zones.push(Zone {
            name: name.clone(),
            start,
            count,
        });
        start += count;
    }
    faces
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::panic)]
    use super::*;
    use crate::mesh::Block;

    fn one_part(points: Vec<Point3>, tets: Vec<[u32; 4]>) -> VolumeMesh {
        let mut mesh = VolumeMesh {
            points,
            blocks: vec![Block {
                name: "b".into(),
                tets,
            }],
            side_sets: Vec::new(),
        };
        mesh.add_boundary_side_sets();
        mesh
    }

    fn valid(mesh: &VolumeMesh) -> (Dual, CellCheck) {
        let d = dual(mesh).unwrap();
        let g = geometry(&d.points, &d.faces);
        let c = check(&d.points, &d.faces, &g);
        assert!(c.valid(), "{c:?}");
        assert!(
            (c.volume - mesh.volume()).abs() <= 1e-12 * mesh.volume().abs(),
            "{} against {}",
            c.volume,
            mesh.volume()
        );
        (d, c)
    }

    /// One tet: four cells, six faces (one per edge, each through the
    /// tet's centroid, two face centroids and the edge's midpoint), and
    /// twelve boundary quadrilaterals, three per triangle.
    #[test]
    fn one_tet() {
        let mesh = one_part(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            vec![[0, 1, 2, 3]],
        );
        let (d, c) = valid(&mesh);
        assert_eq!(d.vertices, [0, 1, 2, 3]);
        assert_eq!(
            (d.faces.cells, d.faces.internal(), d.faces.faces.len()),
            (4, 6, 18)
        );
        assert!(d.faces.faces.iter().take(6).all(|f| f.len() == 4));
        // 1 tet centroid, 4 face centroids, 6 edge midpoints, 4 vertices.
        assert_eq!(d.points.len(), 15);
        assert_eq!(d.faces.patches.len(), 1);
        assert_eq!(d.faces.patches[0].count, 12);
        // The median dual splits a tet into four equal shares.
        let g = geometry(&d.points, &d.faces);
        for v in &g.cell_volumes {
            assert!((v - 1.0 / 24.0).abs() < 1e-15, "{v}");
        }
        assert_eq!(c.cells, 4);
    }

    /// A unit cube of six tets around its diagonal: eight cells, one per
    /// corner, of the cube's volume together; every internal face's normal
    /// points from its owner's vertex towards its neighbour's.
    #[test]
    fn the_cube() {
        let points: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ]
            })
            .collect();
        // The Kuhn triangulation: each tet 0 → 7 along one path of the cube.
        let mut tets = Vec::new();
        for (x, y) in [(1, 2), (1, 4), (2, 1), (2, 4), (4, 1), (4, 2)] {
            let mut t = [0u32, x, x | y, 7];
            let p = |v: u32| points[v as usize];
            let det = dot(
                sub(p(t[1]), p(t[0])),
                cross(sub(p(t[2]), p(t[0])), sub(p(t[3]), p(t[0]))),
            );
            if det < 0.0 {
                t.swap(1, 2);
            }
            tets.push(t);
        }
        let mesh = one_part(points.clone(), tets);
        let (d, c) = valid(&mesh);
        assert_eq!(d.faces.cells, 8);
        assert!((c.volume - 1.0).abs() < 1e-14);
        let g = geometry(&d.points, &d.faces);
        for k in 0..d.faces.internal() {
            let (o, n) = (d.faces.owner[k], d.faces.neighbour[k]);
            let towards = sub(
                points[d.vertices[n as usize] as usize],
                points[d.vertices[o as usize] as usize],
            );
            assert!(dot(g.face_areas[k], towards) > 0.0, "face {k}");
        }
        // Each corner's cell is an eighth of the cube by symmetry of the
        // boundary, not of the inside: every cell is positive.
        assert!(c.min_volume > 0.0);
    }

    /// Two parts are two zones; a vertex in tets of both is refused.
    #[test]
    fn parts_are_zones_and_may_not_share_vertices() {
        let mut points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        points.extend(points.clone().iter().map(|p| add(*p, [3.0, 0.0, 0.0])));
        let mut mesh = VolumeMesh {
            points,
            blocks: vec![
                Block {
                    name: "a".into(),
                    tets: vec![[0, 1, 2, 3]],
                },
                Block {
                    name: "b".into(),
                    tets: vec![[4, 5, 6, 7]],
                },
            ],
            side_sets: Vec::new(),
        };
        mesh.add_boundary_side_sets();
        let (d, _) = valid(&mesh);
        let zones: Vec<(&str, usize, usize)> = d
            .faces
            .zones
            .iter()
            .map(|z| (z.name.as_str(), z.start, z.count))
            .collect();
        assert_eq!(zones, [("a", 0, 4), ("b", 4, 4)]);
        let names: Vec<&str> = d.faces.patches.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names, ["a_boundary", "b_boundary"]);
        mesh.blocks[1].tets[0][0] = 0;
        assert_eq!(dual(&mesh), Err(MeshProblem::SharedVertex { node: 0 }));
    }

    /// Two tets meeting only at an edge: the edge's vertices have a cell in
    /// each, and the edge a face in each; meeting only at a vertex, it has
    /// two cells.
    #[test]
    fn stars_that_fall_apart_have_a_cell_per_piece() {
        let base = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        // Edge 0-1 shared; the second tet mirrored through the plane y = -z.
        let mut points = base.to_vec();
        points.push([0.0, -1.0, 0.0]);
        points.push([0.0, 0.0, -1.0]);
        let mesh = one_part(points.clone(), vec![[0, 1, 2, 3], [0, 1, 4, 5]]);
        let (d, c) = valid(&mesh);
        assert_eq!(d.vertices, [0, 0, 1, 1, 2, 3, 4, 5]);
        assert_eq!((d.faces.cells, d.faces.internal()), (8, 12));
        assert!(c.volume > 0.0);
        // A shared vertex only.
        let mut points = base.to_vec();
        points.extend([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]);
        let mesh = one_part(points, vec![[0, 1, 2, 3], [0, 5, 4, 6]]);
        let (d, _) = valid(&mesh);
        assert_eq!(d.vertices, [0, 0, 1, 2, 3, 4, 5, 6]);
    }

    /// Joining polygons: two squares sharing an edge make a hexagon; two
    /// that share no edge, or touch at a corner only, do not join.
    #[test]
    fn joining_faces() {
        let a = vec![0, 1, 4, 3];
        let b = vec![1, 2, 5, 4];
        assert_eq!(join(&[&a, &b]), Some(vec![0, 1, 2, 5, 4, 3]));
        let far = vec![6, 7, 8, 9];
        assert_eq!(join(&[&a, &far]), None);
        let corner = vec![4, 10, 11, 12];
        assert_eq!(join(&[&a, &corner]), None);
    }

    /// A valid dual is left as it is; a dual with cells made invalid (its
    /// boundary vertices' cells squashed by flat tets at a wall) has them
    /// merged until valid, with fewer cells and the same volume.
    #[test]
    fn merging_fixes_flat_boundary_cells() {
        let points: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ]
            })
            .collect();
        let mut tets = Vec::new();
        for (x, y) in [(1, 2), (1, 4), (2, 1), (2, 4), (4, 1), (4, 2)] {
            tets.push(positive(&points, [0u32, x, x | y, 7]));
        }
        let mesh = one_part(points.clone(), tets.clone());
        let (d, _) = valid(&mesh);
        let m = merge(&d);
        assert_eq!((m.merges, m.rounds), (0, 0));
        assert_eq!(m.dual, d);

        // A slab: a grid of points 4 × 4 in x and y, at z = 0, 0.02 and 1,
        // each cube of it in Kuhn tets: the cells of the z = 0 vertices are
        // flat, their centres off to the side.
        let n = 4u32;
        let zs: [f64; 3] = [0.0, 0.02, 1.0];
        let mut points = Vec::new();
        for &z in &zs {
            for j in 0..n {
                for i in 0..n {
                    // Shear the thin layer so its cells lean.
                    let shift = if z == 0.0 {
                        0.0
                    } else {
                        0.3 * z.min(0.02) / 0.02
                    };
                    points.push([f64::from(i) + shift * f64::from(j % 2), f64::from(j), z]);
                }
            }
        }
        let id = |i: u32, j: u32, k: u32| i + n * j + n * n * k;
        let mut tets = Vec::new();
        for k in 0..2 {
            for j in 0..n - 1 {
                for i in 0..n - 1 {
                    let corner = |b: u32| id(i + (b & 1), j + ((b >> 1) & 1), k + ((b >> 2) & 1));
                    for (x, y) in [(1, 2), (1, 4), (2, 1), (2, 4), (4, 1), (4, 2)] {
                        tets.push(positive(&points, [0, x, x | y, 7].map(corner)));
                    }
                }
            }
        }
        let mesh = one_part(points, tets);
        let d = dual(&mesh).unwrap();
        let before = check(&d.points, &d.faces, &geometry(&d.points, &d.faces));
        let m = merge(&d);
        let after = check(
            &m.dual.points,
            &m.dual.faces,
            &geometry(&m.dual.points, &m.dual.faces),
        );
        m.dual.faces.check(m.dual.points.len()).unwrap();
        assert_eq!(after.open.1, 0, "{after:?}");
        assert!((after.volume - mesh.volume()).abs() < 1e-9, "{after:?}");
        assert_eq!(m.dual.faces.cells + m.merges, d.faces.cells);
        // The sheared layer makes some cells invalid (3 of 48), and merging
        // leaves fewer.
        assert!(before.inverted_pyramids.1 > 0, "{before:?}");
        assert!(m.merges > 0);
        assert!(
            after.inverted_pyramids.1 < before.inverted_pyramids.1,
            "{before:?} {after:?}"
        );
    }

    fn positive(points: &[Point3], mut t: [u32; 4]) -> [u32; 4] {
        let p = |v: u32| points[v as usize];
        let det = dot(
            sub(p(t[1]), p(t[0])),
            cross(sub(p(t[2]), p(t[0])), sub(p(t[3]), p(t[0]))),
        );
        if det < 0.0 {
            t.swap(1, 2);
        }
        t
    }

    #[test]
    fn even_permutations() {
        assert!(even([0, 1, 2, 3]));
        assert!(!even([1, 0, 2, 3]));
        assert!(even([1, 2, 0, 3]));
        assert!(!even([0, 1, 3, 2]));
    }
}
