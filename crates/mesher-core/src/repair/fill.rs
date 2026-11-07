//! Filling holes: each boundary loop that is a simple cycle and small enough
//! is closed by ear-clipping the loop in its own plane, or, if ear clipping
//! cannot (a strongly non-planar loop), by a fan from a new vertex at the
//! loop's centroid.

use std::collections::BTreeMap;

use mesher_geom::Point3;
use mesher_geom::polygon::{TriangulationMethod, triangulate_polygon};
use mesher_geom::vec3::{cross, norm, sub};

use crate::diagnose::Tolerances;
use crate::diagnose::analysis::BoundaryLoop;

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Fill {
    Filled {
        faces_added: u64,
        vertex_added: bool,
        area: f64,
        fraction: f64,
    },
    Failed {
        reason: &'static str,
        fraction: Option<f64>,
    },
}

pub(crate) fn fill_loop(
    points: &mut Vec<Point3>,
    faces: &mut Vec<[u32; 3]>,
    origin: &mut Vec<Option<u32>>,
    hole: &BoundaryLoop,
    shell_area: f64,
    tol: &Tolerances,
) -> Fill {
    if hole.edges.len() > tol.max_hole_edges {
        return Fill::Failed {
            reason: "too_many_edges",
            fraction: None,
        };
    }
    // The loop must be one simple cycle: every vertex with one edge in and one out.
    let mut next: BTreeMap<u32, u32> = BTreeMap::new();
    let mut incoming: BTreeMap<u32, u32> = BTreeMap::new();
    for &([i, j], _) in &hole.edges {
        if next.insert(i, j).is_some() || incoming.insert(j, i).is_some() {
            return Fill::Failed {
                reason: "pinched_loop",
                fraction: None,
            };
        }
    }
    let Some((&start, _)) = next.iter().next() else {
        return Fill::Failed {
            reason: "empty_loop",
            fraction: None,
        };
    };
    let mut cycle = vec![start];
    let mut at = start;
    while let Some(&to) = next.get(&at) {
        if to == start {
            break;
        }
        if cycle.len() > hole.edges.len() {
            break;
        }
        cycle.push(to);
        at = to;
    }
    if cycle.len() != hole.edges.len() || next.get(&at) != Some(&start) {
        return Fill::Failed {
            reason: "pinched_loop",
            fraction: None,
        };
    }

    // The fill runs against the boundary direction.
    cycle.reverse();
    let positions: Vec<Point3> = cycle.iter().map(|&v| points[v as usize]).collect();
    let result = triangulate_polygon(&positions);
    let (new_faces, centroid): (Vec<[u32; 3]>, Option<Point3>) = match result.method {
        TriangulationMethod::EarClip => (
            result
                .triangles
                .iter()
                .map(|t| t.map(|i| cycle[i]))
                .collect(),
            None,
        ),
        TriangulationMethod::Fan(_) => {
            let n = positions.len() as f64;
            let c = [0, 1, 2].map(|k| positions.iter().map(|p| p[k]).sum::<f64>() / n);
            let apex = points.len() as u32;
            let tris = (0..cycle.len())
                .map(|i| [apex, cycle[i], cycle[(i + 1) % cycle.len()]])
                .collect();
            (tris, Some(c))
        }
    };
    let at_point = |v: u32| {
        if v as usize == points.len() {
            centroid.unwrap_or([0.0; 3])
        } else {
            points[v as usize]
        }
    };
    let area: f64 = new_faces
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(at_point);
            norm(cross(sub(b, a), sub(c, a))) / 2.0
        })
        .sum();
    if area <= 0.0 || !area.is_finite() {
        return Fill::Failed {
            reason: "zero_area",
            fraction: None,
        };
    }
    let fraction = area / (shell_area + area);
    if fraction > tol.max_hole_fraction {
        return Fill::Failed {
            reason: "too_large",
            fraction: Some(fraction),
        };
    }
    if let Some(c) = centroid {
        points.push(c);
        origin.push(None);
    }
    let faces_added = new_faces.len() as u64;
    faces.extend(new_faces);
    Fill::Filled {
        faces_added,
        vertex_added: centroid.is_some(),
        area,
        fraction,
    }
}
