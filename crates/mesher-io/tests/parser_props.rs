//! Property tests: the parser never panics, and its output is internally
//! consistent, on arbitrary bytes and on OBJ-like token soup.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_core::{DefectCode, RecordLog};
use mesher_io::parse_obj;
use proptest::prelude::*;

fn check_invariants(bytes: &[u8]) {
    let mut log = RecordLog::new();
    let model = parse_obj(bytes, &mut log, &mut |_| {});
    let nv = model.vertices.len() as u32;
    assert_eq!(model.vertex_lines.len(), model.vertices.len());
    assert!(model.vertices.iter().flatten().all(|c| c.is_finite()));
    for t in &model.triangles {
        assert!(t.vertices.iter().all(|&v| v < nv));
        assert!((t.polygon as usize) < model.polygons.len());
    }
    for p in &model.polygons {
        assert!(p.vertices.len() >= 3);
        assert!(p.vertices.iter().all(|&v| v < nv));
        assert!(p.line >= 1 && u64::from(p.line) <= model.lines);
    }
    let expected: usize = model.polygons.iter().map(|p| p.vertices.len() - 2).sum();
    assert_eq!(model.triangles.len(), expected);
    for r in log.records() {
        assert!(
            r.location
                .obj_lines
                .iter()
                .all(|&l| l >= 1 && u64::from(l) <= model.lines)
        );
    }
    let no_geometry = log
        .records()
        .iter()
        .any(|r| r.code == DefectCode::ParseNoGeometry);
    assert_eq!(no_geometry, model.triangles.is_empty());
}

fn token() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("v".to_owned()),
        Just("f".to_owned()),
        Just("o".to_owned()),
        Just("g".to_owned()),
        Just("vn".to_owned()),
        Just("l".to_owned()),
        Just("#".to_owned()),
        Just("\\".to_owned()),
        Just("\n".to_owned()),
        Just("\r\n".to_owned()),
        Just("\r".to_owned()),
        Just("/".to_owned()),
        Just("nan".to_owned()),
        Just("1e400".to_owned()),
        (-12i64..12).prop_map(|i| i.to_string()),
        (-1e3f64..1e3).prop_map(|x| x.to_string()),
        (1i64..4).prop_map(|i| format!("{i}/{i}/{i}")),
        "[a-z]{1,4}",
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn arbitrary_bytes_never_panic(bytes in prop::collection::vec(any::<u8>(), 0..2048)) {
        check_invariants(&bytes);
    }

    #[test]
    fn obj_like_token_soup_never_panics(tokens in prop::collection::vec(token(), 0..300)) {
        check_invariants(tokens.join(" ").as_bytes());
    }

    #[test]
    fn well_formed_fans_parse_cleanly(n in 3usize..30, scale in 1e-6f64..1e6) {
        let mut text = String::new();
        for i in 0..n {
            let a = i as f64 / n as f64 * std::f64::consts::TAU;
            text.push_str(&format!(
                "v {} {} 0\n",
                scale * libm::cos(a),
                scale * libm::sin(a)
            ));
        }
        text.push('f');
        for i in 1..=n {
            text.push_str(&format!(" {i}"));
        }
        text.push('\n');
        let mut log = RecordLog::new();
        let model = parse_obj(text.as_bytes(), &mut log, &mut |_| {});
        prop_assert!(log.records().is_empty());
        prop_assert_eq!(model.triangles.len(), n - 2);
    }
}
