//! Triangulation of planar-ish 3D polygons (OBJ quads and n-gons).
//!
//! Ear clipping runs in the plane of the polygon's Newell normal, with every
//! convexity and containment test decided by the exact `orient2d` predicate.
//! When ear clipping cannot apply (zero normal, too many vertices, or no ear
//! found because the polygon overlaps itself) the polygon is fan-triangulated and
//! the caller is told why, so it can record the fallback.

use crate::predicates::{Sign, orient2d};

/// Polygons with more vertices than this are fan-triangulated rather than
/// ear-clipped, bounding the quadratic-to-cubic cost of ear clipping on hostile input.
pub const MAX_EAR_CLIP_VERTICES: usize = 256;

/// Why a polygon was fan-triangulated instead of ear-clipped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FanReason {
    /// The Newell normal is zero or non-finite: the polygon has no plane.
    DegenerateNormal,
    /// More than [`MAX_EAR_CLIP_VERTICES`] vertices.
    TooManyVertices,
    /// Ear clipping stalled: the polygon is self-overlapping or degenerate.
    NoEarFound,
}

impl FanReason {
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            FanReason::DegenerateNormal => "degenerate_normal",
            FanReason::TooManyVertices => "too_many_vertices",
            FanReason::NoEarFound => "no_ear_found",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriangulationMethod {
    /// Polygon was a triangle already, or ear clipping succeeded.
    EarClip,
    /// Fan triangulation from the first vertex (for `NoEarFound`, from the
    /// polygon left after the ears that were found).
    Fan(FanReason),
}

/// Triangles as local vertex indices `0..n`, preserving the polygon's winding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Triangulation {
    pub triangles: Vec<[usize; 3]>,
    pub method: TriangulationMethod,
}

/// Triangulate a polygon given by its vertex positions in order.
///
/// Returns `n - 2` triangles for `n >= 3` vertices, and none for fewer. Never panics.
#[must_use]
pub fn triangulate_polygon(points: &[[f64; 3]]) -> Triangulation {
    let n = points.len();
    if n < 3 {
        return Triangulation {
            triangles: Vec::new(),
            method: TriangulationMethod::EarClip,
        };
    }
    if n == 3 {
        return Triangulation {
            triangles: vec![[0, 1, 2]],
            method: TriangulationMethod::EarClip,
        };
    }
    if n > MAX_EAR_CLIP_VERTICES {
        return fan(
            &(0..n).collect::<Vec<_>>(),
            Vec::new(),
            FanReason::TooManyVertices,
        );
    }

    let normal = newell_normal(points);
    let Some(plane) = Projection::from_normal(normal) else {
        return fan(
            &(0..n).collect::<Vec<_>>(),
            Vec::new(),
            FanReason::DegenerateNormal,
        );
    };
    let projected: Vec<[f64; 2]> = points.iter().map(|p| plane.project(*p)).collect();
    let orient = |a: usize, b: usize, c: usize| {
        let s = orient2d(projected[a], projected[b], projected[c]);
        if plane.flip { s.flipped() } else { s }
    };

    let mut ring: Vec<usize> = (0..n).collect();
    let mut triangles = Vec::with_capacity(n - 2);
    while ring.len() > 3 {
        let m = ring.len();
        let ear = (0..m).find(|&i| {
            let (p, c, q) = (ring[(i + m - 1) % m], ring[i], ring[(i + 1) % m]);
            if orient(p, c, q) != Sign::Positive {
                return false;
            }
            let corners = [projected[p], projected[c], projected[q]];
            !ring.iter().any(|&j| {
                j != p
                    && j != c
                    && j != q
                    && !corners.contains(&projected[j])
                    && orient(p, c, j) != Sign::Negative
                    && orient(c, q, j) != Sign::Negative
                    && orient(q, p, j) != Sign::Negative
            })
        });
        match ear {
            Some(i) => {
                triangles.push([ring[(i + m - 1) % m], ring[i], ring[(i + 1) % m]]);
                ring.remove(i);
            }
            None => return fan(&ring, triangles, FanReason::NoEarFound),
        }
    }
    triangles.push([ring[0], ring[1], ring[2]]);
    Triangulation {
        triangles,
        method: TriangulationMethod::EarClip,
    }
}

fn fan(ring: &[usize], mut triangles: Vec<[usize; 3]>, reason: FanReason) -> Triangulation {
    for w in ring.windows(2).skip(1) {
        triangles.push([ring[0], w[0], w[1]]);
    }
    Triangulation {
        triangles,
        method: TriangulationMethod::Fan(reason),
    }
}

fn newell_normal(points: &[[f64; 3]]) -> [f64; 3] {
    let mut n = [0.0; 3];
    for (i, a) in points.iter().enumerate() {
        let b = points[(i + 1) % points.len()];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    n
}

/// Projection onto the coordinate plane most aligned with the polygon. The kept
/// axes are cyclic `(x,y)`, `(y,z)` or `(z,x)`, so a polygon with a positive
/// normal component on the dropped axis stays counter-clockwise; `flip` records
/// a negative component.
struct Projection {
    u: usize,
    v: usize,
    flip: bool,
}

impl Projection {
    fn from_normal(n: [f64; 3]) -> Option<Projection> {
        if !n.iter().all(|c| c.is_finite()) {
            return None;
        }
        let abs = n.map(f64::abs);
        let k = if abs[0] >= abs[1] && abs[0] >= abs[2] {
            0
        } else if abs[1] >= abs[2] {
            1
        } else {
            2
        };
        if abs[k] == 0.0 {
            return None;
        }
        Some(Projection {
            u: (k + 1) % 3,
            v: (k + 2) % 3,
            flip: n[k] < 0.0,
        })
    }

    fn project(&self, p: [f64; 3]) -> [f64; 2] {
        [p[self.u], p[self.v]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn shoelace(points: &[[f64; 2]]) -> f64 {
        let mut a = 0.0;
        for (i, p) in points.iter().enumerate() {
            let q = points[(i + 1) % points.len()];
            a += p[0] * q[1] - q[0] * p[1];
        }
        a / 2.0
    }

    fn tri_area(points: &[[f64; 2]], t: [usize; 3]) -> f64 {
        let [a, b, c] = t.map(|i| points[i]);
        ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2.0
    }

    fn lift(points: &[[f64; 2]]) -> Vec<[f64; 3]> {
        points.iter().map(|p| [p[0], p[1], 0.0]).collect()
    }

    #[test]
    fn quad_gives_two_triangles() {
        let t = triangulate_polygon(&lift(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]));
        assert_eq!(t.method, TriangulationMethod::EarClip);
        assert_eq!(t.triangles.len(), 2);
    }

    #[test]
    fn concave_polygon_is_not_fanned_across_the_notch() {
        // An arrow shape whose fan from vertex 0 would leave the polygon.
        let pts = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [2.0, 1.0], [0.0, 4.0]];
        let t = triangulate_polygon(&lift(&pts));
        assert_eq!(t.method, TriangulationMethod::EarClip);
        let total: f64 = t.triangles.iter().map(|&tr| tri_area(&pts, tr)).sum();
        assert!((total - shoelace(&pts)).abs() < 1e-12);
        assert!(t.triangles.iter().all(|&tr| tri_area(&pts, tr) > 0.0));
    }

    #[test]
    fn clockwise_and_tilted_polygons_keep_winding() {
        // Clockwise in xy: normal points to -z.
        let pts = [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]];
        let t = triangulate_polygon(&lift(&pts));
        assert_eq!(t.method, TriangulationMethod::EarClip);
        assert!(t.triangles.iter().all(|&tr| tri_area(&pts, tr) < 0.0));
        // A quad in the yz plane.
        let yz = [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
        ];
        assert_eq!(
            triangulate_polygon(&yz).method,
            TriangulationMethod::EarClip
        );
    }

    #[test]
    fn collinear_vertices_on_an_edge() {
        let pts = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [0.0, 3.0]];
        let t = triangulate_polygon(&lift(&pts));
        assert_eq!(t.method, TriangulationMethod::EarClip);
        assert_eq!(t.triangles.len(), 3);
    }

    #[test]
    fn degenerate_polygon_falls_back_to_fan() {
        let line = [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0],
        ];
        let t = triangulate_polygon(&line);
        assert_eq!(
            t.method,
            TriangulationMethod::Fan(FanReason::DegenerateNormal)
        );
        assert_eq!(t.triangles.len(), 2);
    }

    #[test]
    fn huge_polygon_falls_back_to_fan() {
        let n = MAX_EAR_CLIP_VERTICES + 1;
        let pts: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let a = i as f64 / n as f64 * std::f64::consts::TAU;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let t = triangulate_polygon(&pts);
        assert_eq!(
            t.method,
            TriangulationMethod::Fan(FanReason::TooManyVertices)
        );
        assert_eq!(t.triangles.len(), n - 2);
    }

    #[test]
    fn self_overlapping_bowtie_falls_back() {
        let pts = [[0.0, 0.0], [2.0, 2.0], [2.0, 0.0], [0.0, 2.0]];
        let t = triangulate_polygon(&lift(&pts));
        assert_eq!(t.triangles.len(), 2);
    }

    fn star(n: usize, radii: &[f64]) -> Vec<[f64; 2]> {
        (0..n)
            .map(|i| {
                let a = i as f64 / n as f64 * std::f64::consts::TAU;
                let r = radii[i % radii.len()];
                [r * a.cos(), r * a.sin()]
            })
            .collect()
    }

    proptest! {
        #[test]
        fn any_polygon_gives_n_minus_2_valid_triangles(
            pts in prop::collection::vec(prop::array::uniform3(-1e3f64..1e3), 3..40)
        ) {
            let t = triangulate_polygon(&pts);
            prop_assert_eq!(t.triangles.len(), pts.len() - 2);
            for tr in &t.triangles {
                prop_assert!(tr.iter().all(|&i| i < pts.len()));
                prop_assert!(tr[0] != tr[1] && tr[1] != tr[2] && tr[0] != tr[2]);
            }
        }

        #[test]
        fn star_shaped_polygons_tile_exactly(
            n in 4usize..60,
            radii in prop::collection::vec(0.2f64..5.0, 1..8),
        ) {
            let pts = star(n, &radii);
            let t = triangulate_polygon(&lift(&pts));
            prop_assert_eq!(t.method, TriangulationMethod::EarClip);
            let total: f64 = t.triangles.iter().map(|&tr| tri_area(&pts, tr)).sum();
            let expected = shoelace(&pts);
            prop_assert!((total - expected).abs() <= 1e-9 * expected.abs().max(1.0));
            // Non-overlapping tiling: no triangle has negative area.
            prop_assert!(t.triangles.iter().all(|&tr| tri_area(&pts, tr) >= 0.0));
        }
    }
}
