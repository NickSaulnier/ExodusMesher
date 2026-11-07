//! Splitting the model into parts.
//!
//! - Any `o` statement: one part per object; faces before the first `o` form
//!   the part `default`.
//! - Otherwise any `g` statement: one part per group, likewise.
//! - Otherwise connected components (through shared merged vertices). A
//!   component that is inward-facing and enclosed by another component is a
//!   cavity, not a part of its own: it joins the innermost component around it.
//!
//! Parts are ordered by their smallest triangle index.

use std::collections::BTreeMap;

use mesher_geom::predicates::collinear;
use mesher_geom::vec3::{cross, dot, sub};
use mesher_geom::{Aabb, FastWindingNumber, Point3, TriangleSoup};

use super::InternalError;
use super::uf::DisjointSets;
use crate::model::SurfaceModel;
use crate::report::PartSource;

pub(crate) struct Group {
    pub name: String,
    pub source: PartSource,
    pub triangles: Vec<u32>,
}

pub(crate) fn split(model: &SurfaceModel, rep: &[u32]) -> Result<Vec<Group>, InternalError> {
    let by_name = |pick: &dyn Fn(u32) -> Option<u32>, names: &[String], source| {
        let mut buckets: BTreeMap<Option<u32>, Vec<u32>> = BTreeMap::new();
        for t in 0..model.triangles.len() as u32 {
            buckets.entry(pick(t)).or_default().push(t);
        }
        let mut groups: Vec<Group> = buckets
            .into_iter()
            .map(|(key, triangles)| Group {
                name: key
                    .and_then(|k| names.get(k as usize).cloned())
                    .unwrap_or_else(|| "default".to_owned()),
                source,
                triangles,
            })
            .collect();
        merge_same_names(&mut groups);
        groups.sort_by_key(|g| g.triangles.first().copied());
        groups
    };
    if !model.objects.is_empty() {
        let pick = |t: u32| model.polygon_of(t).and_then(|p| p.object);
        return Ok(by_name(&pick, &model.objects, PartSource::Object));
    }
    if !model.groups.is_empty() {
        let pick = |t: u32| model.polygon_of(t).and_then(|p| p.group);
        return Ok(by_name(&pick, &model.groups, PartSource::Group));
    }
    components(model, rep)
}

/// An object literally named `default` and the unnamed faces share a name;
/// treat them as one part.
fn merge_same_names(groups: &mut Vec<Group>) {
    let mut merged: Vec<Group> = Vec::with_capacity(groups.len());
    for g in groups.drain(..) {
        match merged.iter_mut().find(|m| m.name == g.name) {
            Some(m) => {
                m.triangles.extend(g.triangles);
                m.triangles.sort_unstable();
            }
            None => merged.push(g),
        }
    }
    *groups = merged;
}

fn components(model: &SurfaceModel, rep: &[u32]) -> Result<Vec<Group>, InternalError> {
    let mut sets = DisjointSets::new(model.vertices.len());
    for t in &model.triangles {
        let [a, b, c] = t.vertices.map(|v| rep[v as usize]);
        sets.union(a, b);
        sets.union(b, c);
    }
    let mut by_root: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (i, t) in model.triangles.iter().enumerate() {
        let root = sets.find(rep[t.vertices[0] as usize]);
        by_root.entry(root).or_default().push(i as u32);
    }
    let mut comps: Vec<Vec<u32>> = by_root.into_values().collect();
    comps.sort_by_key(|c| c.first().copied());

    let points = &model.vertices;
    let faces = |c: &[u32]| -> Vec<[u32; 3]> {
        c.iter()
            .map(|&t| {
                model.triangles[t as usize]
                    .vertices
                    .map(|v| rep[v as usize])
            })
            .filter(|f| !degenerate(points, *f))
            .collect()
    };
    let info: Vec<(f64, Aabb, Option<Point3>)> = comps
        .iter()
        .map(|c| {
            let f = faces(c);
            let bounds = f
                .iter()
                .flatten()
                .fold(Aabb::EMPTY, |b, &v| b.with_point(points[v as usize]));
            let probe = f.first().map(|&t| centroid(points, t));
            (signed_volume(points, &f), bounds, probe)
        })
        .collect();

    // Attach each inward-facing enclosed component to the smallest outward
    // component whose winding number marks the probe point as inside.
    let mut owner: Vec<usize> = (0..comps.len()).collect();
    let mut container_faces: BTreeMap<usize, Vec<[u32; 3]>> = BTreeMap::new();
    for i in 0..comps.len() {
        let (vol, bounds, Some(probe)) = info[i] else {
            continue;
        };
        if vol >= 0.0 {
            continue;
        }
        let mut best: Option<(f64, usize)> = None;
        for j in 0..comps.len() {
            let (vj, bj, _) = info[j];
            if j == i || vj <= 0.0 || !contains_box(&bj, &bounds) {
                continue;
            }
            let fj = container_faces.entry(j).or_insert_with(|| faces(&comps[j]));
            let soup = TriangleSoup::new(points, fj).map_err(|e| InternalError(e.to_string()))?;
            if FastWindingNumber::new(soup).winding_number(probe).abs() > 0.5 {
                let size = box_volume(&bj);
                if best.is_none_or(|(s, k)| (size, j) < (s, k)) {
                    best = Some((size, j));
                }
            }
        }
        if let Some((_, j)) = best {
            owner[i] = j;
        }
    }

    let mut groups: Vec<Group> = Vec::new();
    let mut slot: BTreeMap<usize, usize> = BTreeMap::new();
    for (i, c) in comps.into_iter().enumerate() {
        let mut root = owner[i];
        // Owners are outward components, which own themselves.
        while owner[root] != root {
            root = owner[root];
        }
        match slot.get(&root) {
            Some(&g) => {
                groups[g].triangles.extend(c);
                groups[g].triangles.sort_unstable();
            }
            None => {
                slot.insert(root, groups.len());
                groups.push(Group {
                    name: String::new(),
                    source: PartSource::Component,
                    triangles: c,
                });
            }
        }
    }
    groups.sort_by_key(|g| g.triangles.first().copied());
    for (k, g) in groups.iter_mut().enumerate() {
        g.name = format!("component_{}", k + 1);
    }
    Ok(groups)
}

pub(crate) fn degenerate(points: &[Point3], f: [u32; 3]) -> bool {
    f[0] == f[1]
        || f[1] == f[2]
        || f[0] == f[2]
        || collinear(
            points[f[0] as usize],
            points[f[1] as usize],
            points[f[2] as usize],
        )
}

pub(crate) fn centroid(points: &[Point3], f: [u32; 3]) -> Point3 {
    let [a, b, c] = f.map(|v| points[v as usize]);
    [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0)
}

/// Signed volume enclosed by the faces, taken about the first vertex to limit
/// cancellation.
pub(crate) fn signed_volume(points: &[Point3], faces: &[[u32; 3]]) -> f64 {
    let Some(r) = faces.first().map(|f| points[f[0] as usize]) else {
        return 0.0;
    };
    faces
        .iter()
        .map(|f| {
            let [a, b, c] = f.map(|v| sub(points[v as usize], r));
            dot(a, cross(b, c))
        })
        .sum::<f64>()
        / 6.0
}

fn contains_box(outer: &Aabb, inner: &Aabb) -> bool {
    (0..3).all(|k| outer.min[k] <= inner.min[k] && inner.max[k] <= outer.max[k])
}

fn box_volume(b: &Aabb) -> f64 {
    (0..3).map(|k| (b.max[k] - b.min[k]).max(0.0)).product()
}
