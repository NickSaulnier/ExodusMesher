//! Bounding volume hierarchy over a triangle soup: ray casts (picking), nearest
//! points (deviation measurement), and bounding-box overlap pairs (candidates
//! for exact intersection tests).
//!
//! Boxes are exact: a triangle's box is the min and max of its coordinates and a
//! node's box is the union of its children's, so two triangles that touch always
//! have overlapping boxes and [`Bvh::self_overlaps`] never misses a candidate.
//! Ray and nearest-point queries are floating-point measurements; their culling
//! is padded so it never rejects a node a brute-force scan would accept.

use crate::Point3;
use crate::soup::TriangleSoup;
use crate::vec3::{add, cross, distance_squared, dot, scale, sub};

/// Triangles per leaf at most.
const LEAF_SIZE: usize = 4;

/// Axis-aligned bounding box, closed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    pub min: Point3,
    pub max: Point3,
}

impl Aabb {
    /// The empty box: unions with it are the identity.
    pub const EMPTY: Aabb = Aabb {
        min: [f64::INFINITY; 3],
        max: [f64::NEG_INFINITY; 3],
    };

    #[must_use]
    pub fn of_points(points: &[Point3]) -> Aabb {
        points.iter().fold(Aabb::EMPTY, |b, &p| b.with_point(p))
    }

    #[must_use]
    pub fn with_point(self, p: Point3) -> Aabb {
        Aabb {
            min: [0, 1, 2].map(|k| self.min[k].min(p[k])),
            max: [0, 1, 2].map(|k| self.max[k].max(p[k])),
        }
    }

    #[must_use]
    pub fn union(self, other: Aabb) -> Aabb {
        Aabb {
            min: [0, 1, 2].map(|k| self.min[k].min(other.min[k])),
            max: [0, 1, 2].map(|k| self.max[k].max(other.max[k])),
        }
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        (0..3).any(|k| self.min[k] > self.max[k])
    }

    /// Closed overlap: boxes sharing only a face, edge or corner overlap.
    #[must_use]
    pub fn overlaps(&self, other: &Aabb) -> bool {
        (0..3).all(|k| self.min[k] <= other.max[k] && other.min[k] <= self.max[k])
    }

    /// Closed containment.
    #[must_use]
    pub fn contains(&self, p: Point3) -> bool {
        (0..3).all(|k| self.min[k] <= p[k] && p[k] <= self.max[k])
    }

    #[must_use]
    pub fn distance_squared(&self, p: Point3) -> f64 {
        (0..3)
            .map(|k| {
                let d = (self.min[k] - p[k]).max(p[k] - self.max[k]).max(0.0);
                d * d
            })
            .sum()
    }

    /// Length of the diagonal; zero for an empty box.
    #[must_use]
    pub fn diagonal(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            crate::vec3::norm(sub(self.max, self.min))
        }
    }

    fn center(&self) -> Point3 {
        scale(add(self.min, self.max), 0.5)
    }
}

/// A ray `origin + t·dir`, `t ≥ 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    pub origin: Point3,
    pub dir: Point3,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayHit {
    pub triangle: u32,
    /// Ray parameter of the hit.
    pub t: f64,
    /// Barycentric coordinates of the hit point.
    pub barycentric: [f64; 3],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Nearest {
    pub triangle: u32,
    pub point: Point3,
    pub distance_squared: f64,
}

#[derive(Debug, Clone, Copy)]
enum NodeKind {
    /// Triangles `order[start..end]`.
    Leaf {
        start: u32,
        end: u32,
    },
    Inner {
        left: u32,
        right: u32,
    },
}

#[derive(Debug, Clone, Copy)]
struct Node {
    bbox: Aabb,
    kind: NodeKind,
}

#[derive(Debug, Clone)]
pub struct Bvh<'a> {
    soup: TriangleSoup<'a>,
    nodes: Vec<Node>,
    /// Triangle ids in leaf order.
    order: Vec<u32>,
    /// Box of each triangle, by triangle id.
    boxes: Vec<Aabb>,
}

impl<'a> Bvh<'a> {
    /// Build by median splits along the longest axis of the centroid bounds.
    /// Deterministic, and depth is about `log2(n / LEAF_SIZE)` even when every
    /// centroid coincides.
    #[must_use]
    pub fn new(soup: TriangleSoup<'a>) -> Bvh<'a> {
        let boxes: Vec<Aabb> = soup
            .triangles()
            .iter()
            .map(|&t| Aabb::of_points(&soup.corners(t)))
            .collect();
        let mut bvh = Bvh {
            soup,
            nodes: Vec::new(),
            order: (0..boxes.len() as u32).collect(),
            boxes,
        };
        if !bvh.order.is_empty() {
            let centroids: Vec<Point3> = bvh.boxes.iter().map(Aabb::center).collect();
            bvh.build(0, bvh.order.len(), &centroids);
        }
        bvh
    }

    fn build(&mut self, start: usize, end: usize, centroids: &[Point3]) -> u32 {
        let ids = &self.order[start..end];
        let bbox = ids
            .iter()
            .fold(Aabb::EMPTY, |b, &t| b.union(self.boxes[t as usize]));
        let id = self.nodes.len();
        self.nodes.push(Node {
            bbox,
            kind: NodeKind::Leaf {
                start: start as u32,
                end: end as u32,
            },
        });
        if end - start <= LEAF_SIZE {
            return id as u32;
        }
        let cbox = ids
            .iter()
            .fold(Aabb::EMPTY, |b, &t| b.with_point(centroids[t as usize]));
        let extent = sub(cbox.max, cbox.min);
        let axis = (0..3)
            .max_by(|&i, &j| extent[i].total_cmp(&extent[j]).then(j.cmp(&i)))
            .unwrap_or(0);
        let mid = (start + end) / 2;
        self.order[start..end].select_nth_unstable_by(mid - start, |&a, &b| {
            centroids[a as usize][axis]
                .total_cmp(&centroids[b as usize][axis])
                .then(a.cmp(&b))
        });
        let left = self.build(start, mid, centroids);
        let right = self.build(mid, end, centroids);
        self.nodes[id].kind = NodeKind::Inner { left, right };
        id as u32
    }

    #[must_use]
    pub fn soup(&self) -> TriangleSoup<'a> {
        self.soup
    }

    /// Box of all triangles (empty for an empty soup).
    #[must_use]
    pub fn bounds(&self) -> Aabb {
        self.nodes.first().map_or(Aabb::EMPTY, |n| n.bbox)
    }

    /// Closest hit with `t ≥ 0`; ties go to the lower triangle id.
    #[must_use]
    pub fn ray_cast(&self, ray: &Ray) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;
        let mut stack = self.root();
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n as usize];
            let limit = best.map_or(f64::INFINITY, |b| b.t);
            if !ray_box(ray, &node.bbox, limit) {
                continue;
            }
            match node.kind {
                NodeKind::Leaf { start, end } => {
                    for &t in &self.order[start as usize..end as usize] {
                        let corners = self.soup.corners(self.soup.triangles()[t as usize]);
                        if let Some(hit) = ray_triangle(ray, corners, t) {
                            if best.is_none_or(|b| (hit.t, hit.triangle) < (b.t, b.triangle)) {
                                best = Some(hit);
                            }
                        }
                    }
                }
                NodeKind::Inner { left, right } => stack.extend([right, left]),
            }
        }
        best
    }

    /// Closest point on the surface; ties go to the lower triangle id.
    #[must_use]
    pub fn nearest(&self, p: Point3) -> Option<Nearest> {
        // Pad pruning by a few ulps of the problem scale so floating error in a
        // triangle's distance can never make pruning reject the true nearest.
        let scale = self.bounds().diagonal() + crate::vec3::norm(p);
        let pad = 1e-12 * scale;
        let mut best: Option<Nearest> = None;
        let mut stack = self.root();
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n as usize];
            if let Some(b) = best {
                let reach = b.distance_squared.sqrt() + pad;
                if node.bbox.distance_squared(p) > reach * reach {
                    continue;
                }
            }
            match node.kind {
                NodeKind::Leaf { start, end } => {
                    for &t in &self.order[start as usize..end as usize] {
                        let corners = self.soup.corners(self.soup.triangles()[t as usize]);
                        let q = closest_point_on_triangle(p, corners);
                        let cand = Nearest {
                            triangle: t,
                            point: q,
                            distance_squared: distance_squared(p, q),
                        };
                        if best.is_none_or(|b| {
                            (cand.distance_squared, cand.triangle)
                                < (b.distance_squared, b.triangle)
                        }) {
                            best = Some(cand);
                        }
                    }
                }
                NodeKind::Inner { left, right } => {
                    let dl = self.nodes[left as usize].bbox.distance_squared(p);
                    let dr = self.nodes[right as usize].bbox.distance_squared(p);
                    // Visit the nearer child first (pushed last).
                    if dl <= dr {
                        stack.extend([right, left]);
                    } else {
                        stack.extend([left, right]);
                    }
                }
            }
        }
        best
    }

    /// All pairs `(i, j)`, `i < j`, of triangles whose closed boxes overlap,
    /// sorted. Candidates for exact intersection tests.
    #[must_use]
    pub fn self_overlaps(&self) -> Vec<(u32, u32)> {
        let mut out = Vec::new();
        if !self.nodes.is_empty() {
            self.self_pairs(0, &mut out);
        }
        out.sort_unstable();
        out
    }

    /// All pairs `(i, j)` with triangle `i` of `self` and triangle `j` of
    /// `other` whose closed boxes overlap, sorted.
    #[must_use]
    pub fn overlaps_with(&self, other: &Bvh<'_>) -> Vec<(u32, u32)> {
        let mut out = Vec::new();
        if !self.nodes.is_empty() && !other.nodes.is_empty() {
            cross_pairs(self, 0, other, 0, &mut out);
        }
        out.sort_unstable();
        out
    }

    fn self_pairs(&self, n: u32, out: &mut Vec<(u32, u32)>) {
        match self.nodes[n as usize].kind {
            NodeKind::Leaf { start, end } => {
                let ids = &self.order[start as usize..end as usize];
                for (k, &i) in ids.iter().enumerate() {
                    for &j in &ids[k + 1..] {
                        if self.boxes[i as usize].overlaps(&self.boxes[j as usize]) {
                            out.push((i.min(j), i.max(j)));
                        }
                    }
                }
            }
            NodeKind::Inner { left, right } => {
                self.self_pairs(left, out);
                self.self_pairs(right, out);
                let mut pairs = Vec::new();
                cross_pairs(self, left, self, right, &mut pairs);
                out.extend(pairs.into_iter().map(|(i, j)| (i.min(j), i.max(j))));
            }
        }
    }

    fn root(&self) -> Vec<u32> {
        if self.nodes.is_empty() {
            Vec::new()
        } else {
            vec![0]
        }
    }

    /// Node count and the triangle range of each leaf, for the winding-number
    /// hierarchy.
    pub(crate) fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub(crate) fn node_box(&self, n: u32) -> Aabb {
        self.nodes[n as usize].bbox
    }

    /// `Ok(triangles)` for a leaf, `Err((left, right))` for an inner node.
    pub(crate) fn node_children(&self, n: u32) -> Result<&[u32], (u32, u32)> {
        match self.nodes[n as usize].kind {
            NodeKind::Leaf { start, end } => Ok(&self.order[start as usize..end as usize]),
            NodeKind::Inner { left, right } => Err((left, right)),
        }
    }

    /// Number of triangles under node `n`.
    pub(crate) fn node_triangle_count(&self, n: u32) -> usize {
        match self.node_children(n) {
            Ok(ids) => ids.len(),
            Err((l, r)) => self.node_triangle_count(l) + self.node_triangle_count(r),
        }
    }
}

fn cross_pairs(a: &Bvh<'_>, na: u32, b: &Bvh<'_>, nb: u32, out: &mut Vec<(u32, u32)>) {
    let (node_a, node_b) = (&a.nodes[na as usize], &b.nodes[nb as usize]);
    if !node_a.bbox.overlaps(&node_b.bbox) {
        return;
    }
    match (node_a.kind, node_b.kind) {
        (NodeKind::Leaf { start: sa, end: ea }, NodeKind::Leaf { start: sb, end: eb }) => {
            for &i in &a.order[sa as usize..ea as usize] {
                for &j in &b.order[sb as usize..eb as usize] {
                    if a.boxes[i as usize].overlaps(&b.boxes[j as usize]) {
                        out.push((i, j));
                    }
                }
            }
        }
        (NodeKind::Inner { left, right }, NodeKind::Leaf { .. }) => {
            cross_pairs(a, left, b, nb, out);
            cross_pairs(a, right, b, nb, out);
        }
        (_, NodeKind::Inner { left, right }) => {
            cross_pairs(a, na, b, left, out);
            cross_pairs(a, na, b, right, out);
        }
    }
}

/// Slab test against `[0, limit]`, padded so floating error never rejects a box
/// the ray truly enters (the padding only costs extra triangle tests).
fn ray_box(ray: &Ray, bbox: &Aabb, limit: f64) -> bool {
    const PAD: f64 = 1e-9;
    let (mut t0, mut t1) = (0.0f64, limit);
    for k in 0..3 {
        let (o, d) = (ray.origin[k], ray.dir[k]);
        let span = (bbox.max[k] - bbox.min[k]).abs() + o.abs();
        let (lo, hi) = (bbox.min[k] - PAD * span, bbox.max[k] + PAD * span);
        if d == 0.0 {
            if o < lo || o > hi {
                return false;
            }
            continue;
        }
        let (mut near, mut far) = ((lo - o) / d, (hi - o) / d);
        if near > far {
            std::mem::swap(&mut near, &mut far);
        }
        t0 = t0.max(near);
        t1 = t1.min(far * (1.0 + PAD) + PAD);
        if t0 > t1 {
            return false;
        }
    }
    true
}

/// Two-sided Möller–Trumbore intersection; `None` when parallel or behind.
#[must_use]
pub fn ray_triangle(ray: &Ray, [a, b, c]: [Point3; 3], triangle: u32) -> Option<RayHit> {
    let e1 = sub(b, a);
    let e2 = sub(c, a);
    let p = cross(ray.dir, e2);
    let det = dot(e1, p);
    if det == 0.0 || !det.is_finite() {
        return None;
    }
    let inv = 1.0 / det;
    let s = sub(ray.origin, a);
    let u = dot(s, p) * inv;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross(s, e1);
    let v = dot(ray.dir, q) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = dot(e2, q) * inv;
    (t >= 0.0 && t.is_finite()).then_some(RayHit {
        triangle,
        t,
        barycentric: [1.0 - u - v, u, v],
    })
}

/// Closest point on a (possibly degenerate) triangle, after Ericson,
/// *Real-Time Collision Detection* §5.1.5.
#[must_use]
pub fn closest_point_on_triangle(p: Point3, [a, b, c]: [Point3; 3]) -> Point3 {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 && d1 - d3 > 0.0 {
        return add(a, scale(ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 && d2 - d6 > 0.0 {
        return add(a, scale(ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 && (d4 - d3) + (d5 - d6) > 0.0 {
        return add(b, scale(sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6))));
    }
    let denom = va + vb + vc;
    if denom > 0.0 && denom.is_finite() {
        let v = vb / denom;
        let w = vc / denom;
        return add(a, add(scale(ab, v), scale(ac, w)));
    }
    // Degenerate triangle: fall back to the nearest of its three edges.
    [(a, b), (b, c), (c, a)]
        .into_iter()
        .map(|(x, y)| closest_point_on_segment(p, x, y))
        .min_by(|&q, &r| distance_squared(p, q).total_cmp(&distance_squared(p, r)))
        .unwrap_or(a)
}

#[must_use]
pub fn closest_point_on_segment(p: Point3, a: Point3, b: Point3) -> Point3 {
    let ab = sub(b, a);
    let len2 = dot(ab, ab);
    if len2 <= 0.0 || !len2.is_finite() {
        return a;
    }
    let t = (dot(sub(p, a), ab) / len2).clamp(0.0, 1.0);
    add(a, scale(ab, t))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closest_point_regions() {
        let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert_eq!(closest_point_on_triangle([-1.0, -1.0, 5.0], tri), tri[0]);
        assert_eq!(
            closest_point_on_triangle([0.25, 0.25, 3.0], tri),
            [0.25, 0.25, 0.0]
        );
        assert_eq!(
            closest_point_on_triangle([0.5, -2.0, 0.0], tri),
            [0.5, 0.0, 0.0]
        );
        let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        assert_eq!(
            closest_point_on_triangle([1.5, 1.0, 0.0], flat),
            [1.5, 0.0, 0.0]
        );
    }

    #[test]
    fn empty_bvh_queries() {
        let soup = TriangleSoup::new(&[], &[]).unwrap();
        let bvh = Bvh::new(soup);
        let ray = Ray {
            origin: [0.0; 3],
            dir: [1.0, 0.0, 0.0],
        };
        assert!(bvh.ray_cast(&ray).is_none());
        assert!(bvh.nearest([0.0; 3]).is_none());
        assert!(bvh.self_overlaps().is_empty());
    }
}
