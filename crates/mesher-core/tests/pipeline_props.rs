//! Diagnosis and repair never panic or report an internal error, on random
//! triangle soups dense in coincident, collinear and coplanar configurations,
//! and their output is internally consistent.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::model::{Polygon, SurfaceModel, Triangle};
use mesher_core::repair::repair;
use mesher_core::{JobConfig, RecordLog};
use proptest::prelude::*;

fn model(points: Vec<[f64; 3]>, tris: Vec<[u32; 3]>, objects: usize) -> SurfaceModel {
    let n = tris.len();
    SurfaceModel {
        vertex_lines: (1..=points.len() as u32).collect(),
        vertices: points,
        polygons: tris
            .iter()
            .enumerate()
            .map(|(i, t)| Polygon {
                vertices: t.to_vec(),
                line: i as u32 + 1,
                object: (objects > 0).then(|| (i % objects) as u32),
                group: None,
            })
            .collect(),
        triangles: (0..n)
            .map(|i| Triangle {
                vertices: tris[i],
                polygon: i as u32,
            })
            .collect(),
        objects: (0..objects).map(|k| format!("o{k}")).collect(),
        ..SurfaceModel::default()
    }
}

fn soup() -> impl Strategy<Value = SurfaceModel> {
    (4usize..14, 0usize..3).prop_flat_map(|(n, objects)| {
        (
            prop::collection::vec(prop::array::uniform3((-2i32..=2).prop_map(f64::from)), n),
            prop::collection::vec(prop::array::uniform3(0..n as u32), 1..24),
            Just(objects),
        )
            .prop_map(|(p, t, o)| model(p, t, o))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(600))]

    #[test]
    fn pipeline_is_total_and_consistent(m in soup()) {
        let tol = Tolerances::default();
        let mut log = RecordLog::new();
        let diagnosis = diagnose(&m, &tol, &mut log).unwrap();
        let selected = select_parts(&diagnosis, &JobConfig::default(), &mut log).unwrap();
        let mut overlay = diagnosis.overlay.clone();
        check_part_pairs(&m, &diagnosis, &selected, &mut log, &mut overlay).unwrap();
        let repaired = repair(&m, &diagnosis, &selected, &tol, &mut log).unwrap();

        for r in &repaired {
            let n = r.points.len() as u32;
            prop_assert_eq!(r.vertex_origin.len(), r.points.len());
            prop_assert_eq!(r.face_origin.len(), r.faces.len());
            prop_assert!(r.faces.iter().flatten().all(|&v| v < n));
            prop_assert!(r.faces.iter().all(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]));
            prop_assert_eq!(r.summary.faces_after as usize, r.faces.len());
            prop_assert_eq!(r.summary.vertices_after as usize, r.points.len());
        }
        // Same input, same records: diagnosis and repair are deterministic.
        let mut again = RecordLog::new();
        let d2 = diagnose(&m, &tol, &mut again).unwrap();
        let s2 = select_parts(&d2, &JobConfig::default(), &mut again).unwrap();
        let mut o2 = d2.overlay.clone();
        check_part_pairs(&m, &d2, &s2, &mut again, &mut o2).unwrap();
        repair(&m, &d2, &s2, &tol, &mut again).unwrap();
        prop_assert_eq!(log.records(), again.records());
    }
}
