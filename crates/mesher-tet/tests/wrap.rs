//! Wrapping (Milestone 17): on closed, open, sheet-like and
//! self-intersecting inputs the wrap is closed and oriented outward, no two
//! of its triangles cross, its vertices lie on the offset surface, and it
//! encloses the input.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;
use std::time::Instant;

use mesher_geom::vec3::{cross, dot};
use mesher_geom::{Bvh, FastWindingNumber, Point3, TriangleSoup, faces_intersect};
use mesher_tet::{WrapParams, Wrapped, wrap};

type Surface = (Vec<Point3>, Vec<[u32; 3]>);

/// An axis-aligned box, outward.
fn cuboid(lo: Point3, hi: Point3) -> Surface {
    let p: Vec<Point3> = (0..8)
        .map(|k| {
            [
                if k & 1 == 0 { lo[0] } else { hi[0] },
                if k & 2 == 0 { lo[1] } else { hi[1] },
                if k & 4 == 0 { lo[2] } else { hi[2] },
            ]
        })
        .collect();
    let quads = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let f = quads
        .iter()
        .flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]])
        .collect();
    (p, f)
}

fn join(a: Surface, b: Surface) -> Surface {
    let (mut p, mut f) = a;
    let n = p.len() as u32;
    p.extend(b.0);
    f.extend(b.1.into_iter().map(|t| t.map(|v| v + n)));
    (p, f)
}

fn volume(w: &Wrapped) -> f64 {
    w.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| w.points[v as usize]);
            dot(a, cross(b, c)) / 6.0
        })
        .sum()
}

/// Check the wrap's guarantees against the input.
fn check(name: &str, input: &Surface, params: WrapParams) -> Wrapped {
    let started = Instant::now();
    let w = wrap(&input.0, &input.1, params).unwrap();
    let seconds = started.elapsed().as_secs_f64();

    // Closed and oriented: each edge used as often in each direction.
    let mut uses: BTreeMap<(u32, u32), i64> = BTreeMap::new();
    for t in &w.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += if a < b { 1 } else { -1 };
        }
    }
    assert!(uses.values().all(|&u| u == 0), "{name}: not closed");

    // A 2-manifold: two triangles per edge, one fan per vertex.
    let mut count: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    for t in &w.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *count.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    assert!(count.values().all(|&n| n == 2), "{name}: non-manifold edge");
    let mut fans: Vec<Vec<usize>> = vec![Vec::new(); w.points.len()];
    for (i, t) in w.triangles.iter().enumerate() {
        for &v in t {
            fans[v as usize].push(i);
        }
    }
    for (v, around) in fans.iter().enumerate() {
        // Walk the fan from its first triangle across edges at v.
        let mut reached = vec![around[0]];
        let mut k = 0;
        while k < reached.len() {
            let t = w.triangles[reached[k]];
            k += 1;
            for &o in around {
                if reached.contains(&o) {
                    continue;
                }
                let s = w.triangles[o];
                let shared = t.iter().filter(|x| s.contains(x)).count();
                if shared >= 2 {
                    reached.push(o);
                }
            }
        }
        assert_eq!(
            reached.len(),
            around.len(),
            "{name}: non-manifold vertex {v}"
        );
    }
    assert!(
        !w.triangles.is_empty(),
        "{name}: empty: {} carved, {} + {} points, {} stuck",
        w.carved,
        w.edge_points,
        w.projected_points,
        w.stuck
    );

    // No two triangles cross.
    let soup = TriangleSoup::new(&w.points, &w.triangles).unwrap();
    let bvh = Bvh::new(soup);
    let corners = |t: [u32; 3]| t.map(|v| w.points[v as usize]);
    let crossing = bvh
        .self_overlaps()
        .into_iter()
        .filter(|&(i, j)| {
            let (a, b) = (w.triangles[i as usize], w.triangles[j as usize]);
            faces_intersect(a, corners(a), b, corners(b))
        })
        .count();
    assert_eq!(crossing, 0, "{name}: crossing triangles");

    // Vertices on the offset surface.
    let input_soup = TriangleSoup::new(&input.0, &input.1).unwrap();
    let input_bvh = Bvh::new(input_soup);
    // With the winding guard, a vertex farther out is on the cap: where the
    // input's winding number is 1/2.
    let input_fwn = FastWindingNumber::new(input_soup);
    let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
    for &p in &w.points {
        let d = input_bvh.nearest(p).unwrap().distance_squared.sqrt();
        if params.winding_guard && d > 1.001 * params.offset {
            let n = input_fwn.winding_number(p).abs();
            assert!(
                (n - 0.5).abs() < 1e-6,
                "{name}: vertex {p:?} at {d}, winding {n}"
            );
            continue;
        }
        lo = lo.min(d);
        hi = hi.max(d);
    }
    let delta = params.offset;
    assert!(
        lo > 0.5 * delta && hi < 1.001 * delta,
        "{name}: vertex distances {lo} to {hi}, offset {delta}"
    );

    // The input is enclosed: its vertices are inside the wrap.
    let fwn = FastWindingNumber::new(TriangleSoup::new(&w.points, &w.triangles).unwrap());
    for &p in &input.0 {
        let n = fwn.winding_number(p);
        assert!(
            (n - 1.0).abs() < 1e-6,
            "{name}: input point {p:?} has winding {n}"
        );
    }
    let v = volume(&w);
    assert!(v > 0.0, "{name}: volume {v}");
    println!(
        "{name}: {} triangles, {} points ({} on edges, {} projected), {} carved, {} stuck, non-manifold {} edges {} vertices ({} relabelled), volume {v:.4}, distance {lo:.4e} to {hi:.4e}, {seconds:.2} s",
        w.triangles.len(),
        w.points.len(),
        w.edge_points,
        w.projected_points,
        w.carved,
        w.stuck,
        w.nonmanifold_edges,
        w.nonmanifold_vertices,
        w.relabelled
    );
    w
}

const PARAMS: WrapParams = WrapParams {
    alpha: 0.1,
    offset: 0.01,
    max_points: 1_000_000,
    winding_guard: false,
};

#[test]
fn a_cube_is_wrapped_closely() {
    let cube = cuboid([0.0; 3], [1.0; 3]);
    let w = check("cube", &cube, PARAMS);
    // The offset cube with rounded edges: a little more than the cube.
    let v = volume(&w);
    assert!(v > 1.0 && v < 1.0 + 7.0 * PARAMS.offset, "volume {v}");
}

#[test]
fn an_open_box_is_wrapped_as_a_shell() {
    // A box missing its top: the hole is far wider than 2α, so the wrap
    // goes in and wraps both sides of the walls.
    let (p, mut f) = cuboid([0.0; 3], [1.0; 3]);
    f.retain(|t| !t.iter().all(|&v| p[v as usize][2] == 1.0));
    let w = check("open box", &(p, f), PARAMS);
    assert!(volume(&w) < 0.2, "volume {}", volume(&w));
}

#[test]
fn a_single_triangle_is_wrapped() {
    let t = (
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    );
    check("triangle", &t, PARAMS);
}

#[test]
fn overlapping_cubes_are_wrapped_as_their_union() {
    let input = join(
        cuboid([0.0; 3], [1.0; 3]),
        cuboid([0.5, 0.3, 0.2], [1.5, 1.2, 0.9]),
    );
    let w = check("overlapping cubes", &input, PARAMS);
    let union = 1.0 + 0.5 * 0.9 * 0.7 + 0.5 * 0.2 * 0.7;
    let v = volume(&w);
    assert!(v > union && v < union * 1.1, "volume {v}, union {union}");
}

#[test]
fn a_gap_narrower_than_alpha_is_closed() {
    // Two boxes 0.05 apart (less than 2α): one wrap around both.
    let input = join(
        cuboid([0.0; 3], [1.0; 3]),
        cuboid([1.05, 0.0, 0.0], [2.05, 1.0, 1.0]),
    );
    let w = check("narrow gap", &input, PARAMS);
    let v = volume(&w);
    assert!(v > 2.05, "volume {v}");
}

/// Near-touching boxes and random triangle soups give carvings with
/// non-manifold edges and vertices; the wrap removes them all.
#[test]
fn near_contacts_and_soups_are_wrapped_as_manifolds() {
    let mut found = 0;
    let params = WrapParams {
        alpha: 0.05,
        offset: 0.005,
        max_points: 1_000_000,
        winding_guard: false,
    };
    for gap in [0.0, 0.005, 0.0101, 0.02, 0.06] {
        // Along an edge, and at a corner.
        let edge = join(
            cuboid([0.0; 3], [1.0; 3]),
            cuboid([1.0 + gap, 1.0 + gap, 0.0], [2.0, 2.0, 1.0]),
        );
        let corner = join(cuboid([0.0; 3], [1.0; 3]), cuboid([1.0 + gap; 3], [2.0; 3]));
        for (name, input) in [("edge contact", edge), ("corner contact", corner)] {
            let w = check(&format!("{name} {gap}"), &input, params);
            found += w.nonmanifold_edges + w.nonmanifold_vertices;
        }
    }
    let mut seed: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for round in 0..8 {
        let mut points = Vec::new();
        let mut triangles = Vec::new();
        for _ in 0..60 {
            let c = [next(), next(), next()];
            let n = points.len() as u32;
            for _ in 0..3 {
                points.push([0, 1, 2].map(|k| c[k] + 0.2 * (next() - 0.5)));
            }
            triangles.push([n, n + 1, n + 2]);
        }
        let w = check(&format!("soup {round}"), &(points, triangles), params);
        found += w.nonmanifold_edges + w.nonmanifold_vertices;
    }
    assert!(found > 0, "no case exercised the manifold repair");
}

/// With the winding guard, an open box is wrapped as the box: the wrap
/// closes the hole across the winding number's level set.
#[test]
fn the_winding_guard_closes_an_open_box() {
    let (p, mut f) = cuboid([0.0; 3], [1.0; 3]);
    f.retain(|t| !t.iter().all(|&v| p[v as usize][2] == 1.0));
    let params = WrapParams {
        winding_guard: true,
        ..PARAMS
    };
    let w = check("open box, guarded", &(p, f), params);
    let v = volume(&w);
    assert!(v > 0.95 && v < 1.07, "volume {v}");
    assert!(w.cap_points > 0);
}
