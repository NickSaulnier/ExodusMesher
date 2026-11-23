//! How well a part's feature lines and named patches survive meshing
//! (Milestone 18, stage 5): the input's feature curves (`features`, on the
//! repaired surface) against the feature edges of the mesh's boundary, and
//! each OBJ group's area against its side set's.
//!
//! For each input feature edge, points every quarter of the element size
//! are measured to the nearest output feature edge; the largest distance
//! (relative to the size) and the share of feature length within a tenth of
//! the size of one say how sharp edges came through. A conformal mesh of the
//! input keeps them exactly; a remeshed one to rounding; stuffing rounds
//! them off.

use std::collections::BTreeMap;

use mesher_geom::Point3;
use mesher_geom::bvh::closest_point_on_segment;
use mesher_geom::vec3::{cross, norm, sub};
use mesher_tet::remesh::{RemeshParams, features};

use crate::features::{face_groups, remesh_params};
use crate::model::SurfaceModel;
use crate::repair::RepairedPart;
use crate::report::JobConfig;
use crate::volume::MeshedPart;

/// Faces of a positive tet, outward.
const OUTWARD: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

/// How one part's features and patches came through.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Survival {
    /// Total length of the input's feature edges.
    pub input_length: f64,
    /// The largest distance from an input feature point to the output's
    /// feature edges, over the element size (capped at `SEARCH` sizes).
    pub max_distance: f64,
    /// The share of the input's feature length within a tenth of the size
    /// of an output feature edge.
    pub kept: f64,
    /// Per OBJ group: its area on the input and its side set's area.
    pub patches: Vec<(String, f64, f64)>,
    /// The input feature point farthest from the output's feature edges.
    pub worst_point: Option<Point3>,
}

/// How far, in element sizes, the search for an output feature edge goes.
pub const SEARCH: f64 = 4.0;

/// The share of the feature angle an output edge must turn by to count.
pub const OUTPUT_ANGLE: f64 = 0.5;

fn area(a: Point3, b: Point3, c: Point3) -> f64 {
    0.5 * norm(cross(sub(b, a), sub(c, a)))
}

/// Measure part `r`'s mesh `m` at element size `h`.
#[must_use]
pub fn survival(
    model: &SurfaceModel,
    r: &RepairedPart,
    m: &MeshedPart,
    config: &JobConfig,
    h: f64,
) -> Survival {
    let params = remesh_params(r, config);
    let groups_in = face_groups(model, r);
    let input = features(&r.points, &r.faces, &groups_in, &params);

    // The mesh's boundary, outward, with each face's group.
    let mut count: BTreeMap<[u32; 3], ([u32; 3], u32)> = BTreeMap::new();
    for t in &m.tets {
        for f in OUTWARD {
            let face = f.map(|k| t[k]);
            let mut key = face;
            key.sort_unstable();
            count.entry(key).or_insert((face, 0)).1 += 1;
        }
    }
    let mut group_of: BTreeMap<[u32; 3], u32> = BTreeMap::new();
    for (name, faces) in &m.face_groups {
        let g = model
            .groups
            .iter()
            .position(|n| n == name)
            .map_or(u32::MAX, |g| g as u32);
        for f in faces {
            let mut key = *f;
            key.sort_unstable();
            group_of.insert(key, g);
        }
    }
    let (mut boundary, mut groups_out) = (Vec::new(), Vec::new());
    for (key, (face, n)) in &count {
        if *n == 1 {
            boundary.push(*face);
            groups_out.push(group_of.get(key).copied().unwrap_or(u32::MAX));
        }
    }
    // An output edge counts where it turns by half the feature angle: the
    // remesher keeps a feature line, but smooths faceted curved faces beside
    // it, which can take its turn below the angle (a 51° feature between
    // 28° facets comes out at 30° to 38°).
    let output = features(
        &m.points,
        &boundary,
        &groups_out,
        &RemeshParams {
            feature_degrees: params.feature_degrees * OUTPUT_ANGLE,
            ..RemeshParams::default()
        },
    );

    // Output feature edges in a grid of cells of size h.
    let cell = h.max(f64::MIN_POSITIVE);
    let key = |p: Point3| p.map(|x| libm::floor(x / cell) as i64);
    let mut grid: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
    for (s, &(u, v)) in output.segments.iter().enumerate() {
        let (a, b) = (m.points[u as usize], m.points[v as usize]);
        let (ka, kb) = (key(a), key(b));
        for i in ka[0].min(kb[0])..=ka[0].max(kb[0]) {
            for j in ka[1].min(kb[1])..=ka[1].max(kb[1]) {
                for k in ka[2].min(kb[2])..=ka[2].max(kb[2]) {
                    grid.entry([i, j, k]).or_default().push(s);
                }
            }
        }
    }
    let rings = libm::ceil(SEARCH) as i64;
    let nearest = |p: Point3| -> f64 {
        let c = key(p);
        let mut best = SEARCH * h;
        for i in -rings..=rings {
            for j in -rings..=rings {
                for k in -rings..=rings {
                    for &s in grid
                        .get(&[c[0] + i, c[1] + j, c[2] + k])
                        .into_iter()
                        .flatten()
                    {
                        let (u, v) = output.segments[s];
                        let q =
                            closest_point_on_segment(p, m.points[u as usize], m.points[v as usize]);
                        best = best.min(norm(sub(p, q)));
                    }
                }
            }
        }
        best
    };
    let mut s = Survival::default();
    let mut kept = 0.0;
    for &(u, v) in &input.segments {
        let (a, b) = (r.points[u as usize], r.points[v as usize]);
        let length = norm(sub(b, a));
        let n = (libm::ceil(length / (0.25 * h)) as usize).max(1);
        // Samples at the middles of n equal pieces, each standing for its
        // piece's length.
        for i in 0..n {
            let t = (i as f64 + 0.5) / n as f64;
            let p = [0, 1, 2].map(|k| a[k] + t * (b[k] - a[k]));
            let d = nearest(p);
            if s.worst_point.is_none() || d / h > s.max_distance {
                s.worst_point = Some(p);
            }
            s.max_distance = s.max_distance.max(d / h);
            if d <= 0.1 * h {
                kept += length / n as f64;
            }
        }
        s.input_length += length;
    }
    s.kept = if s.input_length > 0.0 {
        kept / s.input_length
    } else {
        1.0
    };

    // Patches: each group's area on the input and in its side set.
    let mut input_area: BTreeMap<u32, f64> = BTreeMap::new();
    for (t, &g) in r.faces.iter().zip(&groups_in) {
        if g != u32::MAX {
            let [a, b, c] = t.map(|v| r.points[v as usize]);
            *input_area.entry(g).or_insert(0.0) += area(a, b, c);
        }
    }
    for (g, a) in input_area {
        let Some(name) = model.groups.get(g as usize) else {
            continue;
        };
        let out: f64 = m
            .face_groups
            .iter()
            .filter(|(n, _)| n == name)
            .flat_map(|(_, faces)| faces)
            .map(|f| {
                let [p, q, w] = f.map(|v| m.points[v as usize]);
                area(p, q, w)
            })
            .sum();
        s.patches.push((name.clone(), a, out));
    }
    s
}
