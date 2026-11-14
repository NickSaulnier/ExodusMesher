//! Closed surfaces against the solid an input encloses, for measuring
//! (Milestone 17: a wrapped part's mesh against stuffing's):
//!
//! ```text
//! cargo run --release -p mesher-tet --example surface_compare -- input.obj a.obj b.obj ...
//! ```
//!
//! The solid is where the input's generalized winding number is at least 1/2
//! in magnitude. Prints, per candidate (a closed, outward surface, such as a
//! volume mesh's boundary): on a 48³ grid over the input's box, the
//! intersection over union of the solid and the candidate's inside, the
//! fraction of the solid kept and the candidate's extra over the solid; the
//! largest distance from the centres of the input faces bounding the solid
//! (winding number below 1/2 just off one side) to the candidate, and from
//! the candidate's vertices to the input, both over the diagonal. The OBJ is
//! read naively (`v` and polygon `f` lines, fanned).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_geom::vec3::{add, cross, norm, scale, sub};
use mesher_geom::{Aabb, Bvh, FastWindingNumber, Point3, TriangleSoup};

fn read(path: &str) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut points = Vec::new();
    let mut triangles = Vec::new();
    for line in text.lines() {
        let mut w = line.split_whitespace();
        match w.next() {
            Some("v") => {
                let c: Vec<f64> = w.take(3).map(|x| x.parse().unwrap()).collect();
                points.push([c[0], c[1], c[2]]);
            }
            Some("f") => {
                let v: Vec<u32> = w
                    .map(|x| x.split('/').next().unwrap().parse::<u32>().unwrap() - 1)
                    .collect();
                for k in 1..v.len().saturating_sub(1) {
                    triangles.push([v[0], v[k], v[k + 1]]);
                }
            }
            _ => {}
        }
    }
    (points, triangles)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (points, triangles) = read(&args[0]);
    let soup = TriangleSoup::new(&points, &triangles).unwrap();
    let fwn = FastWindingNumber::new(soup);
    let input_bvh = Bvh::new(soup);
    let b = Aabb::of_points(&points);
    let diagonal = b.diagonal();
    let off = 1e-4 * diagonal;
    let bounding: Vec<Point3> = triangles
        .iter()
        .filter_map(|t| {
            let [p, q, r] = t.map(|v| points[v as usize]);
            let n = cross(sub(q, p), sub(r, p));
            let len = norm(n);
            if len.is_nan() || len <= 0.0 {
                return None;
            }
            let m = [0, 1, 2].map(|k| (p[k] + q[k] + r[k]) / 3.0);
            let step = scale(n, off / len);
            let outside = |x: Point3| fwn.winding_number(x).abs() < 0.5;
            (outside(add(m, step)) || outside(sub(m, step))).then_some(m)
        })
        .collect();
    let n = 48;
    let grid: Vec<Point3> = (0..n * n * n)
        .map(|i| {
            let q = [i / (n * n), (i / n) % n, i % n].map(|x| (x as f64 + 0.5) / n as f64);
            [0, 1, 2].map(|a| b.min[a] + q[a] * (b.max[a] - b.min[a]))
        })
        .collect();
    let solid: Vec<bool> = grid
        .iter()
        .map(|&p| fwn.winding_number(p).abs() >= 0.5)
        .collect();
    println!("candidate\tiou\tkept\textra\tlost\tbulge");
    for path in &args[1..] {
        let (cp, ct) = read(path);
        let csoup = TriangleSoup::new(&cp, &ct).unwrap();
        let cfwn = FastWindingNumber::new(csoup);
        let cbvh = Bvh::new(csoup);
        let (mut both, mut s, mut w) = (0usize, 0usize, 0usize);
        for (p, &inside) in grid.iter().zip(&solid) {
            let c = cfwn.winding_number(*p) >= 0.5;
            s += usize::from(inside);
            w += usize::from(c);
            both += usize::from(inside && c);
        }
        let iou = both as f64 / (s + w - both).max(1) as f64;
        let kept = both as f64 / s.max(1) as f64;
        let extra = (w - both) as f64 / s.max(1) as f64;
        let lost = bounding
            .iter()
            .map(|&p| cbvh.nearest(p).unwrap().distance_squared)
            .fold(0.0f64, f64::max)
            .sqrt()
            / diagonal;
        let bulge = cp
            .iter()
            .map(|&p| input_bvh.nearest(p).unwrap().distance_squared)
            .fold(0.0f64, f64::max)
            .sqrt()
            / diagonal;
        println!("{path}\t{iou:.3}\t{kept:.3}\t{extra:.3}\t{lost:.4}\t{bulge:.4}");
    }
}
