//! Shared generators and exact references for the property tests.

#![allow(
    dead_code,
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::needless_range_loop
)]

use std::collections::BTreeMap;

use mesher_geom::Point3;
use proptest::prelude::*;

/// Small integer coordinates: dense in coincident, collinear and coplanar cases.
pub fn grid_point() -> impl Strategy<Value = [i64; 3]> {
    prop::array::uniform3(-3i64..=3)
}

pub fn to_f64(p: [i64; 3]) -> Point3 {
    p.map(|c| c as f64)
}

/// Random float triangle soup over `n_points` points.
pub fn float_soup(
    max_points: usize,
    max_tris: usize,
) -> impl Strategy<Value = (Vec<Point3>, Vec<[u32; 3]>)> {
    (3..max_points).prop_flat_map(move |n| {
        (
            prop::collection::vec(prop::array::uniform3(-10.0f64..10.0), n),
            prop::collection::vec(prop::array::uniform3(0..n as u32), 0..max_tris),
        )
    })
}

/// Random soup on the integer grid (degenerate-heavy).
pub fn grid_soup(
    max_points: usize,
    max_tris: usize,
) -> impl Strategy<Value = (Vec<Point3>, Vec<[u32; 3]>)> {
    (3..max_points).prop_flat_map(move |n| {
        (
            prop::collection::vec(grid_point().prop_map(to_f64), n),
            prop::collection::vec(prop::array::uniform3(0..n as u32), 0..max_tris),
        )
    })
}

/// Unit icosphere, outward-oriented.
pub fn icosphere(subdivisions: u32) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut pts: Vec<Point3> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut tris: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..subdivisions {
        let mut mid: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        let mut next = Vec::with_capacity(tris.len() * 4);
        for [a, b, c] in tris {
            let mut m = |i: u32, j: u32| {
                *mid.entry((i.min(j), i.max(j))).or_insert_with(|| {
                    let (p, q) = (pts[i as usize], pts[j as usize]);
                    pts.push([
                        (p[0] + q[0]) / 2.0,
                        (p[1] + q[1]) / 2.0,
                        (p[2] + q[2]) / 2.0,
                    ]);
                    (pts.len() - 1) as u32
                })
            };
            let (ab, bc, ca) = (m(a, b), m(b, c), m(c, a));
            next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        tris = next;
    }
    for p in &mut pts {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        *p = p.map(|c| c / r);
    }
    (pts, tris)
}

/// Closed surface of the cube `[0, n]^3`, unit quads split in two, outward,
/// with vertices shared along edges and corners.
pub fn grid_cube(n: i64) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let mut ids: BTreeMap<[i64; 3], u32> = BTreeMap::new();
    let mut pts = Vec::new();
    let mut tris = Vec::new();
    let x = [1, 0, 0];
    let y = [0, 1, 0];
    let z = [0, 0, 1];
    // (origin, u, v) with u × v the outward normal.
    let faces = [
        ([0, 0, 0], y, x),
        ([0, 0, n], x, y),
        ([0, 0, 0], x, z),
        ([0, n, 0], z, x),
        ([0, 0, 0], z, y),
        ([n, 0, 0], y, z),
    ];
    for (o, u, v) in faces {
        let mut id = |i: i64, j: i64| {
            let p = [0, 1, 2].map(|k| o[k] + i * u[k] + j * v[k]);
            *ids.entry(p).or_insert_with(|| {
                pts.push(to_f64(p));
                (pts.len() - 1) as u32
            })
        };
        for i in 0..n {
            for j in 0..n {
                let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            }
        }
    }
    (pts, tris)
}

/// Exact rationals over i128, always reduced with a positive denominator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Q {
    pub n: i128,
    pub d: i128,
}

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

impl Q {
    pub fn int(n: i128) -> Q {
        Q { n, d: 1 }
    }
    fn new(n: i128, d: i128) -> Q {
        let g = gcd(n, d).max(1) * d.signum();
        Q { n: n / g, d: d / g }
    }
    pub fn add(self, o: Q) -> Q {
        Q::new(self.n * o.d + o.n * self.d, self.d * o.d)
    }
    pub fn sub(self, o: Q) -> Q {
        Q::new(self.n * o.d - o.n * self.d, self.d * o.d)
    }
    pub fn mul(self, o: Q) -> Q {
        Q::new(self.n * o.n, self.d * o.d)
    }
    pub fn div(self, o: Q) -> Q {
        Q::new(self.n * o.d, self.d * o.n)
    }
    pub fn is_zero(self) -> bool {
        self.n == 0
    }
}

/// Every basic feasible solution `(λ0, λ1, λ2, μ0, μ1, μ2)` of
/// `Σλ_i a_i = Σμ_j b_j`, `Σλ = Σμ = 1`, `λ, μ ≥ 0`. The triangles (convex hulls,
/// so degenerate ones are segments or points) intersect iff this is non-empty,
/// and any linear function of `(λ, μ)` attains its extremes over the
/// intersection at one of these (fundamental theorem of linear programming).
pub fn basic_feasible_solutions(a: [[i64; 3]; 3], b: [[i64; 3]; 3]) -> Vec<[Q; 6]> {
    let mut columns = [[0i128; 5]; 6];
    for i in 0..3 {
        columns[i] = [a[i][0] as i128, a[i][1] as i128, a[i][2] as i128, 1, 0];
        columns[i + 3] = [-b[i][0] as i128, -b[i][1] as i128, -b[i][2] as i128, 0, 1];
    }
    let rhs = [0i128, 0, 0, 1, 1];
    let mut out = Vec::new();
    for mask in 1u32..64 {
        let support: Vec<usize> = (0..6).filter(|&c| mask & (1 << c) != 0).collect();
        let k = support.len();
        if k > 5 {
            continue;
        }
        let mut m: Vec<Vec<Q>> = (0..5)
            .map(|r| {
                let mut row: Vec<Q> = support.iter().map(|&c| Q::int(columns[c][r])).collect();
                row.push(Q::int(rhs[r]));
                row
            })
            .collect();
        let mut pivots = Vec::new();
        let mut row = 0;
        for col in 0..k {
            let Some(p) = (row..5).find(|&r| !m[r][col].is_zero()) else {
                continue;
            };
            m.swap(row, p);
            let lead = m[row][col];
            for c in 0..=k {
                m[row][c] = m[row][c].div(lead);
            }
            for r in 0..5 {
                if r != row && !m[r][col].is_zero() {
                    let f = m[r][col];
                    for c in 0..=k {
                        let v = m[row][c].mul(f);
                        m[r][c] = m[r][c].sub(v);
                    }
                }
            }
            pivots.push(col);
            row += 1;
        }
        if row < k || (row..5).any(|r| !m[r][k].is_zero()) {
            continue; // dependent columns or inconsistent
        }
        let mut x = [Q::int(0); 6];
        for (i, &col) in pivots.iter().enumerate() {
            x[support[col]] = m[i][k];
        }
        if x.iter().all(|q| q.n >= 0) {
            out.push(x);
        }
    }
    out
}

/// Exact integer cross product.
pub fn cross_i(a: [i64; 3], b: [i64; 3], c: [i64; 3]) -> [i64; 3] {
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}
