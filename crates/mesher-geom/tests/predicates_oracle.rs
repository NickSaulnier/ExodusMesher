//! orient3d and insphere against exact integer and rational references on
//! degenerate-heavy integer grids, and under exact scalings and translations.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

mod common;

use common::{Q, cross_i, grid_point, to_f64};
use mesher_geom::Sign;
use mesher_geom::predicates::{insphere, orient3d};
use proptest::prelude::*;

fn sign(v: i128) -> Sign {
    match v.signum() {
        1 => Sign::Positive,
        -1 => Sign::Negative,
        _ => Sign::Zero,
    }
}

fn det_i(a: [i64; 3], b: [i64; 3], c: [i64; 3], d: [i64; 3]) -> i128 {
    let n = cross_i(a, b, c);
    (0..3).map(|k| n[k] as i128 * (d[k] - a[k]) as i128).sum()
}

/// Inside / on / outside the circumsphere via the exact rational circumcentre.
fn sphere_side(a: [i64; 3], b: [i64; 3], c: [i64; 3], d: [i64; 3], e: [i64; 3]) -> Sign {
    // Solve 2(q - a)·x = |q|² - |a|² for q in b, c, d by Cramer's rule.
    let rows: Vec<([i128; 3], i128)> = [b, c, d]
        .iter()
        .map(|q| {
            let r = [0, 1, 2].map(|k| 2 * (q[k] - a[k]) as i128);
            let s = (0..3).map(|k| (q[k] * q[k] - a[k] * a[k]) as i128).sum();
            (r, s)
        })
        .collect();
    let det3 = |m: [[i128; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let m = [rows[0].0, rows[1].0, rows[2].0];
    let det = det3(m);
    let centre: Vec<Q> = (0..3)
        .map(|k| {
            let mut mk = m;
            for r in 0..3 {
                mk[r][k] = rows[r].1;
            }
            Q::int(det3(mk)).div(Q::int(det))
        })
        .collect();
    let dist2 = |p: [i64; 3]| {
        (0..3)
            .map(|k| {
                let d = Q::int(p[k] as i128).sub(centre[k]);
                d.mul(d)
            })
            .fold(Q::int(0), Q::add)
    };
    let diff = dist2(a).sub(dist2(e)); // positive when e is inside
    sign(diff.n)
}

proptest! {
    // Flat tetrahedra are skipped by insphere's assumption; allow plenty of rejects.
    #![proptest_config(ProptestConfig {
        cases: 2000,
        max_global_rejects: 1 << 20,
        ..ProptestConfig::default()
    })]

    #[test]
    fn orient3d_matches_integer_determinant(
        a in grid_point(), b in grid_point(), c in grid_point(), d in grid_point()
    ) {
        let expected = sign(det_i(a, b, c, d));
        prop_assert_eq!(orient3d(to_f64(a), to_f64(b), to_f64(c), to_f64(d)), expected);
        // Exact scaling by 2^-40 and translation by 2^10 change nothing.
        let t = |p: [i64; 3]| to_f64(p).map(|x| x * 2f64.powi(-40) + 1024.0);
        prop_assert_eq!(orient3d(t(a), t(b), t(c), t(d)), expected);
    }

    #[test]
    fn insphere_matches_rational_circumcentre(
        a in grid_point(), b in grid_point(), c in grid_point(), d in grid_point(), e in grid_point()
    ) {
        let orientation = sign(det_i(a, b, c, d));
        prop_assume!(orientation != Sign::Zero);
        let side = sphere_side(a, b, c, d, e);
        let expected = if orientation == Sign::Positive { side } else { side.flipped() };
        prop_assert_eq!(insphere(to_f64(a), to_f64(b), to_f64(c), to_f64(d), to_f64(e)), expected);
    }
}
