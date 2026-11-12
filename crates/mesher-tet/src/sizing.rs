//! Sizing field: the target element size at each point of a part.
//!
//! Sizes come from the surface. Each vertex gets the smaller of
//!
//! - a *curvature* size: across a smooth edge whose faces turn by φ (below
//!   the feature angle; sharper edges are features, not curvature), elements
//!   should turn by at most θ each, so the size is (θ / φ) times the distance
//!   between the two face centroids, but never below half that distance
//!   (finer than the facets would not follow the input any better);
//! - a *thickness* size: a ray cast inward from each face centroid measures
//!   the local thickness of the solid, and an element should be no larger
//!   than it (`thickness` times it).
//!
//! These samples are clamped to `[min_size, max_size]` and spread through a
//! background grid so that the size grows by at most `gradation` per unit
//! distance (Dijkstra over the grid's 26-neighbourhood, which overestimates
//! Euclidean distance by at most a few percent). Queries interpolate the grid
//! trilinearly. The field is deterministic: it depends only on the surface
//! and the parameters.

use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap};

use mesher_geom::bvh::{Aabb, Bvh, Ray, closest_point_on_segment, closest_point_on_triangle};
use mesher_geom::vec3::{add, cross, dot, norm, scale, sub};
use mesher_geom::{Point3, TriangleSoup};
use rayon::prelude::*;

/// How the field is built.
#[derive(Debug, Clone, PartialEq)]
pub struct SizingParams {
    /// The size where nothing asks for smaller (the job's element size).
    pub max_size: f64,
    /// No sample asks for less than this.
    pub min_size: f64,
    /// Largest growth of the size per unit distance.
    pub gradation: f64,
    /// Target turning angle per element on curved surfaces, degrees.
    pub curvature_degrees: f64,
    /// Edges turning by more than this are features, not curvature, degrees.
    pub feature_degrees: f64,
    /// Elements across the local thickness: size = thickness / this.
    pub thickness_cells: f64,
    /// Background grid nodes allowed; the cell grows to stay below it.
    pub max_nodes: usize,
}

impl SizingParams {
    /// Defaults for a job whose element size is `h`.
    #[must_use]
    pub fn for_size(h: f64) -> SizingParams {
        SizingParams {
            max_size: h,
            min_size: h / 4.0,
            gradation: 0.4,
            curvature_degrees: 25.0,
            feature_degrees: 45.0,
            thickness_cells: 1.0,
            max_nodes: 2_000_000,
        }
    }
}

/// A target element size over space.
#[derive(Clone, PartialEq)]
pub struct SizeField {
    origin: Point3,
    cell: f64,
    dims: [usize; 3],
    /// Node values, x fastest; empty for a uniform field.
    values: Vec<f64>,
    max: f64,
    min: f64,
    gradation: f64,
}

impl std::fmt::Debug for SizeField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SizeField")
            .field("min", &self.min)
            .field("max", &self.max)
            .field("dims", &self.dims)
            .field("cell", &self.cell)
            .finish_non_exhaustive()
    }
}

/// A vertex budget for meshing the region a closed, outward surface bounds
/// at size `h`, or by `field` capped at `h`: 20 × (volume / h³ + the
/// surface's area over its local size squared + the input's vertices),
/// clamped to [10,000, 4,000,000]. Elements fill the volume at the maximum
/// size, plus the surface at its local size; where the size is smaller, its
/// graded layer (about h / gradation thick) adds ∫ 1/(h + g d)³ dd =
/// 1/(2 g h²) per unit area over what the maximum size would put there.
#[must_use]
pub fn point_budget(
    points: &[Point3],
    triangles: &[[u32; 3]],
    h: f64,
    field: Option<&SizeField>,
) -> usize {
    let g = field.map_or(0.0, SizeField::gradation);
    let (mut surface, mut volume) = (0.0f64, 0.0f64);
    for t in triangles {
        let [a, b, c] = t.map(|v| points[v as usize]);
        let area = norm(cross(sub(b, a), sub(c, a))) / 2.0;
        let local = field.map_or(h, |f| {
            f.at(scale(add(add(a, b), c), 1.0 / 3.0)).clamp(f.min(), h)
        });
        surface += area / (local * local);
        if g > 0.0 && local < h {
            surface += area * (1.0 / (local * local) - 1.0 / (h * h)) / (2.0 * g);
        }
        // Relative to a vertex: far from the origin the products cancel.
        let o = points[0];
        volume += dot(sub(a, o), cross(sub(b, o), sub(c, o))) / 6.0;
    }
    let estimate = volume.abs() / h.powi(3) + surface + points.len() as f64;
    ((20.0 * estimate) as usize).clamp(10_000, 4_000_000)
}

/// The closest point to `p` on the line through `a` and `b`.
fn closest_point_on_line(p: Point3, a: Point3, b: Point3) -> Point3 {
    let d = sub(b, a);
    let l = dot(d, d);
    if l > 0.0 {
        add(a, scale(d, dot(sub(p, a), d) / l))
    } else {
        a
    }
}

/// A piece of the surface that asks for a smaller size.
enum Seed {
    Segment(Point3, Point3),
    Triangle([Point3; 3]),
}

/// An `f64` that orders by `total_cmp`, for the priority queue.
#[derive(Clone, Copy, PartialEq)]
struct Key(f64);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl SizeField {
    /// The same size everywhere.
    #[must_use]
    pub fn uniform(h: f64) -> SizeField {
        SizeField {
            origin: [0.0; 3],
            cell: 1.0,
            dims: [0; 3],
            values: Vec::new(),
            max: h,
            min: h,
            gradation: 0.0,
        }
    }

    /// The field of a closed, consistently oriented (outward) surface. A
    /// surface that cannot be indexed gives a uniform field.
    #[must_use]
    pub fn from_surface(points: &[Point3], triangles: &[[u32; 3]], p: &SizingParams) -> SizeField {
        let h0 = p.max_size;
        if points.is_empty() || !(h0 > 0.0 && h0.is_finite()) {
            return SizeField::uniform(h0);
        }
        let Ok(soup) = TriangleSoup::new(points, triangles) else {
            return SizeField::uniform(h0);
        };
        let corners = |t: &[u32; 3]| t.map(|v| points[v as usize]);
        // Unit normals (zero for degenerate triangles) and centroids.
        let normals: Vec<Point3> = triangles
            .iter()
            .map(|t| {
                let [a, b, c] = corners(t);
                let n = cross(sub(b, a), sub(c, a));
                let l = norm(n);
                if l > 0.0 { scale(n, 1.0 / l) } else { [0.0; 3] }
            })
            .collect();
        let centroids: Vec<Point3> = triangles
            .iter()
            .map(|t| {
                let [a, b, c] = corners(t);
                scale(add(add(a, b), c), 1.0 / 3.0)
            })
            .collect();
        let clamp = |h: f64| h.clamp(p.min_size, h0);
        // What asks for less than the maximum: curved edges and thin faces.
        let mut seeds: Vec<(Seed, f64)> = Vec::new();

        // Curvature, per smooth manifold edge.
        let mut edges: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (t, tri) in triangles.iter().enumerate() {
            for k in 0..3 {
                let (u, v) = (tri[k], tri[(k + 1) % 3]);
                edges
                    .entry((u.min(v), u.max(v)))
                    .or_default()
                    .push(t as u32);
            }
        }
        let (theta, feature) = (
            p.curvature_degrees.to_radians(),
            p.feature_degrees.to_radians(),
        );
        let mut curved: Vec<((u32, u32), f64)> = Vec::new();
        for (&(u, v), around) in &edges {
            let &[s, t] = around.as_slice() else {
                continue;
            };
            let (ns, nt) = (normals[s as usize], normals[t as usize]);
            if ns == [0.0; 3] || nt == [0.0; 3] {
                continue;
            }
            let phi = libm::acos(dot(ns, nt).clamp(-1.0, 1.0));
            if !(phi > 1e-6 && phi < feature) {
                continue;
            }
            // The arc across the edge, from the middle of one face to the
            // middle of the other: a centroid is a third of the way across
            // its face, so 1.5 times the centroids' distances to the edge.
            // Then φ / w is the curvature and θ w / φ is θ times the radius.
            let (a, b) = (points[u as usize], points[v as usize]);
            let w = [s, t]
                .iter()
                .map(|&f| {
                    norm(sub(
                        centroids[f as usize],
                        closest_point_on_line(centroids[f as usize], a, b),
                    ))
                })
                .sum::<f64>()
                * 1.5;
            let h = clamp((theta / phi * w).max(0.5 * w));
            if h < h0 {
                curved.push(((u, v), h));
            }
        }
        // Sorted, so the seeds (and with them the field) never depend on
        // hash order.
        curved.sort_by_key(|c| c.0);
        for ((u, v), h) in curved {
            seeds.push((Seed::Segment(points[u as usize], points[v as usize]), h));
        }

        // Thickness, by a ray inward from each face centroid. The origin
        // steps off the face by a hair so the face itself is not hit.
        let bvh = Bvh::new(soup);
        let diag = bvh.bounds().diagonal();
        let thickness: Vec<Option<f64>> = (0..triangles.len())
            .into_par_iter()
            .map(|t| {
                let n = normals[t];
                if n == [0.0; 3] {
                    return None;
                }
                let step = 1e-9 * diag;
                let ray = Ray {
                    origin: sub(centroids[t], scale(n, step)),
                    dir: scale(n, -1.0),
                };
                bvh.ray_cast(&ray).map(|hit| hit.t + step)
            })
            .collect();
        for (t, d) in thickness.into_iter().enumerate() {
            if let Some(d) = d {
                let h = clamp(d / p.thickness_cells);
                if h < h0 {
                    seeds.push((Seed::Triangle(corners(&triangles[t])), h));
                }
            }
        }
        if seeds.is_empty() {
            return SizeField::uniform(h0);
        }
        let min = seeds.iter().map(|s| s.1).fold(h0, f64::min);

        // Background grid over the bounds, padded by a cell.
        let bounds: Aabb = bvh.bounds();
        let extent = sub(bounds.max, bounds.min);
        let volume = extent.iter().map(|e| e.max(h0)).product::<f64>();
        let cell = (h0 / 6.0).max(libm::cbrt(volume / p.max_nodes as f64));
        let origin = sub(bounds.min, [cell; 3]);
        let dims = extent.map(|e| (e / cell).ceil() as usize + 3);
        let n = dims[0] * dims[1] * dims[2];
        let mut values = vec![h0; n];
        let index = |i: usize, j: usize, k: usize| i + dims[0] * (j + dims[1] * k);
        let node = |i: usize, j: usize, k: usize| {
            add(origin, [i as f64 * cell, j as f64 * cell, k as f64 * cell])
        };
        let g = p.gradation;
        // Every node within a cell of a seed gets its exact graded value;
        // the propagation below carries it further.
        let mut heap: BinaryHeap<Reverse<(Key, usize)>> = BinaryHeap::new();
        for (seed, h) in &seeds {
            let corners: &[Point3] = match seed {
                Seed::Segment(a, b) => &[*a, *b],
                Seed::Triangle(t) => t,
            };
            let lo = Aabb::of_points(corners);
            let range = |a: usize| {
                let first = ((lo.min[a] - origin[a]) / cell - 1.0).floor().max(0.0) as usize;
                let last =
                    (((lo.max[a] - origin[a]) / cell + 1.0).ceil() as usize).min(dims[a] - 1);
                first..=last
            };
            for k in range(2) {
                for j in range(1) {
                    for i in range(0) {
                        let q = node(i, j, k);
                        let closest = match seed {
                            Seed::Segment(a, b) => closest_point_on_segment(q, *a, *b),
                            Seed::Triangle(t) => closest_point_on_triangle(q, *t),
                        };
                        let d = norm(sub(q, closest));
                        if d > 1.75 * cell {
                            continue;
                        }
                        let v = h + g * d;
                        let x = index(i, j, k);
                        if v < values[x] {
                            values[x] = v;
                            heap.push(Reverse((Key(v), x)));
                        }
                    }
                }
            }
        }
        let steps: Vec<(isize, isize, isize, f64)> = (-1..=1isize)
            .flat_map(|a| (-1..=1isize).flat_map(move |b| (-1..=1isize).map(move |c| (a, b, c))))
            .filter(|&o| o != (0, 0, 0))
            .map(|(a, b, c)| (a, b, c, g * cell * ((a * a + b * b + c * c) as f64).sqrt()))
            .collect();
        while let Some(Reverse((Key(v), x))) = heap.pop() {
            if v > values[x] {
                continue;
            }
            let (i, j, k) = (
                x % dims[0],
                (x / dims[0]) % dims[1],
                x / (dims[0] * dims[1]),
            );
            for &(a, b, c, d) in &steps {
                let (ni, nj, nk) = (i as isize + a, j as isize + b, k as isize + c);
                if ni < 0
                    || nj < 0
                    || nk < 0
                    || ni as usize >= dims[0]
                    || nj as usize >= dims[1]
                    || nk as usize >= dims[2]
                {
                    continue;
                }
                let y = index(ni as usize, nj as usize, nk as usize);
                let w = v + d;
                if w < values[y] {
                    values[y] = w;
                    heap.push(Reverse((Key(w), y)));
                }
            }
        }
        SizeField {
            origin,
            cell,
            dims,
            values,
            max: h0,
            min,
            gradation: g,
        }
    }

    /// The target size at `p`.
    #[must_use]
    pub fn at(&self, p: Point3) -> f64 {
        if self.values.is_empty() {
            return self.max;
        }
        let d = self.dims;
        let mut base = [0usize; 3];
        let mut frac = [0.0f64; 3];
        for a in 0..3 {
            let x = ((p[a] - self.origin[a]) / self.cell).clamp(0.0, (d[a] - 1) as f64);
            let i = (x.floor() as usize).min(d[a] - 2);
            base[a] = i;
            frac[a] = x - i as f64;
        }
        let mut h = 0.0;
        for m in 0..8 {
            let o = [m & 1, (m >> 1) & 1, m >> 2];
            let w: f64 = (0..3)
                .map(|a| if o[a] == 1 { frac[a] } else { 1.0 - frac[a] })
                .product();
            if w > 0.0 {
                let x = base[0] + o[0] + d[0] * (base[1] + o[1] + d[1] * (base[2] + o[2]));
                h += w * self.values[x];
            }
        }
        h.min(self.max)
    }

    /// The largest size (the job's element size).
    #[must_use]
    pub fn max(&self) -> f64 {
        self.max
    }

    /// The smallest size any sample asks for.
    #[must_use]
    pub fn min(&self) -> f64 {
        self.min
    }

    /// Largest growth of the size per unit distance (0 for a uniform field).
    #[must_use]
    pub fn gradation(&self) -> f64 {
        self.gradation
    }

    /// Whether the size is the same everywhere.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.values.is_empty()
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    /// Fixed parameters, so the tests do not move with the defaults.
    fn params(h: f64) -> SizingParams {
        SizingParams {
            min_size: h / 8.0,
            curvature_degrees: 15.0,
            ..SizingParams::for_size(h)
        }
    }

    /// An axis-aligned box as 12 outward triangles.
    fn boxed(lo: Point3, hi: Point3) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let p: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    if i & 1 == 0 { lo[0] } else { hi[0] },
                    if i & 2 == 0 { lo[1] } else { hi[1] },
                    if i & 4 == 0 { lo[2] } else { hi[2] },
                ]
            })
            .collect();
        let f = vec![
            [0, 2, 3],
            [0, 3, 1],
            [4, 5, 7],
            [4, 7, 6],
            [0, 1, 5],
            [0, 5, 4],
            [2, 6, 7],
            [2, 7, 3],
            [0, 4, 6],
            [0, 6, 2],
            [1, 3, 7],
            [1, 7, 5],
        ];
        (p, f)
    }

    #[test]
    fn a_thick_box_is_uniform() {
        let (p, f) = boxed([0.0; 3], [1.0; 3]);
        let field = SizeField::from_surface(&p, &f, &params(0.2));
        assert!(field.is_uniform(), "{field:?}");
        assert_eq!(field.at([0.5; 3]), 0.2);
    }

    #[test]
    fn a_thin_plate_is_graded_to_its_thickness() {
        // 4 × 4 × 0.05: the thickness asks for 0.05, the job for 0.5.
        let (p, f) = boxed([0.0; 3], [4.0, 4.0, 0.05]);
        let params = params(0.5);
        let field = SizeField::from_surface(&p, &f, &params);
        assert!(!field.is_uniform());
        assert!(
            (field.min() - 0.0625).abs() < 1e-12,
            "clamped to h/8: {}",
            field.min()
        );
        let inside = field.at([2.0, 2.0, 0.025]);
        assert!(inside < 0.15, "{inside}");
        // Away from the plate the size grows by at most the gradation.
        let away = field.at([2.0, 2.0, 1.0]);
        assert!(
            away > inside && away <= inside + 1.05 * 0.4 * 1.0 + 0.05,
            "{away}"
        );
        assert!(field.at([2.0, 2.0, 50.0]) <= 0.5);
    }

    #[test]
    fn grading_is_lipschitz() {
        let (p, f) = boxed([0.0; 3], [4.0, 4.0, 0.05]);
        let params = params(0.5);
        let field = SizeField::from_surface(&p, &f, &params);
        let mut x = 0.0;
        while x < 4.0 {
            let (a, b) = ([x, 1.0, 0.5], [x + 0.1, 1.0, 0.5]);
            let d = (field.at(a) - field.at(b)).abs();
            assert!(d <= 0.4 * 0.1 * 1.1 + 1e-12, "{x}: {d}");
            x += 0.1;
        }
    }

    /// A closed cylinder of `n` sides, radius 1, from z = 0 to `height`.
    fn cylinder(n: u32, height: f64) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let mut p = Vec::new();
        for z in [0.0, height] {
            for i in 0..n {
                let a = f64::from(i) / f64::from(n) * std::f64::consts::TAU;
                p.push([libm::cos(a), libm::sin(a), z]);
            }
        }
        p.push([0.0, 0.0, 0.0]);
        p.push([0.0, 0.0, height]);
        let (bottom, top) = (2 * n, 2 * n + 1);
        let mut f = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            f.push([i, j, n + j]);
            f.push([i, n + j, n + i]);
            f.push([bottom, j, i]);
            f.push([top, n + i, n + j]);
        }
        (p, f)
    }

    #[test]
    fn curvature_asks_for_smaller_elements_on_the_curved_side() {
        // A 15° target on a radius-1 cylinder asks for about 0.26 (θ times
        // the radius), however finely the side is tessellated; the flat caps
        // and their sharp rims ask for nothing.
        for sides in [16, 128] {
            let (p, f) = cylinder(sides, 4.0);
            let field = SizeField::from_surface(&p, &f, &params(0.5));
            let side = field.at([1.0, 0.0, 2.0]);
            assert!((0.22..0.32).contains(&side), "{sides} sides: {side}");
            assert!(
                field.at([0.0, 0.0, 2.0]) > side,
                "{sides} sides: the axis is farther"
            );
        }
        // With a larger radius the curvature asks for nothing.
        let (p, f) = cylinder(16, 4.0);
        let p: Vec<Point3> = p.iter().map(|q| [q[0] * 4.0, q[1] * 4.0, q[2]]).collect();
        let field = SizeField::from_surface(&p, &f, &params(0.5));
        assert!(
            field.at([4.0, 0.0, 2.0]) > 0.45,
            "{}",
            field.at([4.0, 0.0, 2.0])
        );
    }

    #[test]
    fn queries_outside_the_grid_are_clamped() {
        let (p, f) = boxed([0.0; 3], [4.0, 4.0, 0.05]);
        let field = SizeField::from_surface(&p, &f, &params(0.5));
        for q in [[-100.0, 2.0, 0.0], [2.0, 200.0, 0.0], [1e9, -1e9, 1e9]] {
            let h = field.at(q);
            assert!(h > 0.0 && h <= 0.5, "{q:?}: {h}");
        }
    }
}
