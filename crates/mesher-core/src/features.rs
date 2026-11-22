//! Feature lines of a part's repaired surface (Milestone 18): the edges the
//! remesher keeps as lines (where the surface turns by more than the feature
//! angle, where OBJ groups meet, and the user's forced edges, less the forced
//! smooth ones), reported per part and written to `features.vtp`.
//!
//! Overrides name edges by parsed vertex indices (`JobConfig::feature_edges`,
//! `smooth_edges`), as report locations do. They apply to every repaired edge
//! whose two vertices came from those model vertices; one that matches no
//! edge of any selected part, or would smooth a group boundary, is a
//! `FEAT_OVERRIDE_IGNORED` warning.

use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::vec3::{norm, sub};
use mesher_tet::remesh::{FeatureCause, Features, RemeshParams, features};

use crate::model::SurfaceModel;
use crate::repair::RepairedPart;
use crate::report::{DefectCode, FeatureSummary, JobConfig, Record, RecordLog, Stage};

/// The feature lines of one part.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PartFeatures {
    pub features: Features,
    /// The feature curve of each segment: segments chained through vertices
    /// on one smooth line, numbered in segment order.
    pub curves: Vec<u32>,
    /// The angle used, degrees.
    pub angle_degrees: f64,
}

impl PartFeatures {
    /// Each segment's cause as `features.vtp` numbers it: 0 the angle, 1 a
    /// group boundary, 2 forced.
    #[must_use]
    pub fn cause_codes(&self) -> Vec<u8> {
        self.features
            .causes
            .iter()
            .map(|c| match c {
                FeatureCause::Angle => 0,
                FeatureCause::Group => 1,
                FeatureCause::Forced => 2,
            })
            .collect()
    }

    /// Counts and total length for the report.
    #[must_use]
    pub fn summary(&self, points: &[mesher_geom::Point3]) -> FeatureSummary {
        let f = &self.features;
        let count = |c: FeatureCause| f.causes.iter().filter(|&&x| x == c).count() as u64;
        FeatureSummary {
            feature_angle_degrees: self.angle_degrees,
            edges: f.segments.len() as u64,
            angle_edges: count(FeatureCause::Angle),
            group_edges: count(FeatureCause::Group),
            forced_edges: count(FeatureCause::Forced),
            smoothed_edges: f.smoothed as u64,
            corners: f.corners.len() as u64,
            curves: self.curves.iter().max().map_or(0, |&c| u64::from(c) + 1),
            length: f
                .segments
                .iter()
                .map(|&(u, v)| norm(sub(points[u as usize], points[v as usize])))
                .sum(),
        }
    }
}

/// Each repaired face's OBJ group (`u32::MAX` for a fill face or none), so
/// group boundaries are feature lines.
#[must_use]
pub fn face_groups(model: &SurfaceModel, r: &RepairedPart) -> Vec<u32> {
    r.face_origin
        .iter()
        .map(|o| {
            o.and_then(|t| model.polygon_of(t))
                .and_then(|p| p.group)
                .unwrap_or(u32::MAX)
        })
        .collect()
}

/// The override edges of `pairs` (model vertex pairs) on the repaired
/// surface, and which of `pairs` matched one.
fn matched(r: &RepairedPart, pairs: &[[u32; 2]]) -> (Vec<(u32, u32)>, BTreeSet<usize>) {
    let mut wanted: BTreeMap<(u32, u32), Vec<usize>> = BTreeMap::new();
    for (i, &[a, b]) in pairs.iter().enumerate() {
        wanted.entry((a.min(b), a.max(b))).or_default().push(i);
    }
    let mut edges = BTreeSet::new();
    let mut hit = BTreeSet::new();
    if wanted.is_empty() {
        return (Vec::new(), hit);
    }
    let origin = |v: u32| r.vertex_origin.get(v as usize).copied().flatten();
    for t in &r.faces {
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            let (Some(a), Some(b)) = (origin(u), origin(v)) else {
                continue;
            };
            if let Some(which) = wanted.get(&(a.min(b), a.max(b))) {
                edges.insert((u.min(v), u.max(v)));
                hit.extend(which.iter().copied());
            }
        }
    }
    (edges.into_iter().collect(), hit)
}

/// The remesher's parameters for a part: the feature angle and the user's
/// overrides on it.
#[must_use]
pub fn remesh_params(r: &RepairedPart, config: &JobConfig) -> RemeshParams {
    RemeshParams {
        feature_degrees: config.remesh_feature_degrees(),
        forced_features: matched(r, &config.feature_edges).0,
        forced_smooth: matched(r, &config.smooth_edges).0,
        ..RemeshParams::default()
    }
}

/// The feature lines of a part's repaired surface.
#[must_use]
pub fn part_features(model: &SurfaceModel, r: &RepairedPart, config: &JobConfig) -> PartFeatures {
    let params = remesh_params(r, config);
    let f = features(&r.points, &r.faces, &face_groups(model, r), &params);
    let curves = curves(r.points.len(), &f);
    PartFeatures {
        features: f,
        curves,
        angle_degrees: params.feature_degrees,
    }
}

/// Chain segments into curves: two segments at a vertex that is not a
/// corner belong to one curve.
fn curves(n: usize, f: &Features) -> Vec<u32> {
    let mut parent: Vec<usize> = (0..f.segments.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let corner: BTreeSet<u32> = f.corners.iter().copied().collect();
    let mut at: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (s, &(u, v)) in f.segments.iter().enumerate() {
        at[u as usize].push(s);
        at[v as usize].push(s);
    }
    for (v, segs) in at.iter().enumerate() {
        if segs.len() == 2 && !corner.contains(&(v as u32)) {
            let (a, b) = (find(&mut parent, segs[0]), find(&mut parent, segs[1]));
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut ids: BTreeMap<usize, u32> = BTreeMap::new();
    (0..f.segments.len())
        .map(|s| {
            let root = find(&mut parent, s);
            let next = ids.len() as u32;
            *ids.entry(root).or_insert(next)
        })
        .collect()
}

/// Record the overrides that took no effect: matching no edge of any
/// selected part, or asking to smooth a group boundary.
pub fn record_ignored_overrides(
    model: &SurfaceModel,
    repaired: &[RepairedPart],
    config: &JobConfig,
    log: &mut RecordLog,
) {
    if config.feature_edges.is_empty() && config.smooth_edges.is_empty() {
        return;
    }
    let mut hit_features = BTreeSet::new();
    let mut hit_smooth = BTreeSet::new();
    let mut boundaries: BTreeSet<(u32, u32)> = BTreeSet::new();
    for r in repaired {
        hit_features.extend(matched(r, &config.feature_edges).1);
        let (edges, hit) = matched(r, &config.smooth_edges);
        hit_smooth.extend(hit);
        if edges.is_empty() {
            continue;
        }
        // Group boundaries among the edges to smooth, as model vertex pairs.
        let f = features(
            &r.points,
            &r.faces,
            &face_groups(model, r),
            &RemeshParams {
                feature_degrees: config.remesh_feature_degrees(),
                ..RemeshParams::default()
            },
        );
        let wanted: BTreeSet<(u32, u32)> = edges.into_iter().collect();
        for (&(u, v), &cause) in f.segments.iter().zip(&f.causes) {
            if cause == FeatureCause::Group && wanted.contains(&(u, v)) {
                let origin = |x: u32| r.vertex_origin.get(x as usize).copied().flatten();
                if let (Some(a), Some(b)) = (origin(u), origin(v)) {
                    boundaries.insert((a.min(b), a.max(b)));
                }
            }
        }
    }
    let ignored = |log: &mut RecordLog, setting: &str, [a, b]: [u32; 2], reason: &str| {
        log.push(
            Record::new(DefectCode::FeatOverrideIgnored, Stage::VolumeMesh)
                .vertices([a, b])
                .measured("setting", setting)
                .measured("reason", reason)
                .hint("the override took no effect; name an edge by the 0-based vertex indices of report locations"),
        );
    };
    for (i, &pair) in config.feature_edges.iter().enumerate() {
        if !hit_features.contains(&i) {
            ignored(
                log,
                "feature_edges",
                pair,
                "no selected part's repaired surface has this edge",
            );
        }
    }
    for (i, &[a, b]) in config.smooth_edges.iter().enumerate() {
        if !hit_smooth.contains(&i) {
            ignored(
                log,
                "smooth_edges",
                [a, b],
                "no selected part's repaired surface has this edge",
            );
        } else if boundaries.contains(&(a.min(b), a.max(b))) {
            ignored(
                log,
                "smooth_edges",
                [a, b],
                "a group boundary stays a feature line",
            );
        }
    }
}
