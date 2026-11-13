//! Implicit points and exact predicates on them (Milestone 11).
//!
//! An implicit point is defined by how it was constructed, from explicit
//! (f64) points and f64 parameters, and is never rounded:
//!
//! - [`Kind::Lnc`]: a point on a segment, `t·a + (1−t)·b`.
//! - [`Kind::Bpt`]: a point of a triangle's plane by barycentric weights,
//!   `(u·a + v·b + w·c) / (u + v + w)`.
//! - [`Kind::Lpi`]: where the line through `p` and `q` meets the plane
//!   through `r`, `s` and `t`.
//! - [`Kind::Tpi`]: the common point of three planes, each through three
//!   points.
//!
//! Each has homogeneous coordinates `(X, Y, Z, W)`, polynomials in the
//! defining coordinates, and the predicates are signs of determinants of
//! them. A predicate is written once over [`arith::Ring`] and evaluated in
//! stages, stopping at the first that certifies the sign: an f64 filter with a
//! rigorous error bound, interval arithmetic, then exact dyadic arithmetic.
//! Coordinates are first translated by an explicit point of the query, inside
//! the arithmetic, so the filter's bounds scale with the size of the
//! configuration rather than its distance from the origin.
//!
//! Conventions match [`crate::predicates`]: [`orient3d`] is `Positive` for a
//! positive-volume tet, [`insphere`] is `Positive` when `e` is inside the
//! sphere of a positively oriented `abcd`. [`Point::approx`] is for spatial
//! queries and output only, never for a topological decision.

mod arith;

use thiserror::Error;

use crate::Point3;
use crate::predicates::{self, Sign};
use arith::{Exact, Filtered, Interval, Ring, det};

/// What defines an implicit point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Kind {
    Explicit(Point3),
    /// `t·a + (1−t)·b`.
    Lnc {
        a: Point3,
        b: Point3,
        t: f64,
    },
    /// `(u·a + v·b + w·c) / (u + v + w)`.
    Bpt {
        tri: [Point3; 3],
        weights: [f64; 3],
    },
    /// The line through `line[0]` and `line[1]` meets the plane through the
    /// three points of `plane`.
    Lpi {
        line: [Point3; 2],
        plane: [Point3; 3],
    },
    /// The common point of three planes, each through three points.
    Tpi {
        planes: [[Point3; 3]; 3],
    },
}

/// A point, explicit or implicit. Built only by the checked constructors, so
/// every point exists: its homogeneous `W` is exactly nonzero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    kind: Kind,
    /// The exact sign of `W` is negative.
    w_negative: bool,
}

/// Why an implicit point could not be built.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum ImplicitError {
    #[error("a coordinate or parameter is not finite")]
    NonFinite,
    #[error("the barycentric weights sum to zero")]
    ZeroWeight,
    #[error("the line is parallel to the plane, or the plane is degenerate")]
    NoLinePlaneIntersection,
    #[error("the three planes have no single common point")]
    NoThreePlaneIntersection,
}

/// Which stage settled a predicate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Stage {
    /// All points explicit: Shewchuk's adaptive predicates.
    Explicit,
    /// The f64 filter with its error bound.
    Filter,
    /// Interval arithmetic.
    Interval,
    /// Exact arithmetic.
    Exact,
}

fn finite(points: &[Point3]) -> bool {
    points.iter().flatten().all(|c| c.is_finite())
}

impl Point {
    /// An explicit point.
    ///
    /// # Errors
    /// [`ImplicitError::NonFinite`] for a non-finite coordinate.
    pub fn explicit(p: Point3) -> Result<Point, ImplicitError> {
        if !finite(&[p]) {
            return Err(ImplicitError::NonFinite);
        }
        Ok(Point {
            kind: Kind::Explicit(p),
            w_negative: false,
        })
    }

    /// `t·a + (1−t)·b`: on the segment `ab` for `t` in `[0, 1]`.
    ///
    /// # Errors
    /// [`ImplicitError::NonFinite`] for a non-finite input.
    pub fn lnc(a: Point3, b: Point3, t: f64) -> Result<Point, ImplicitError> {
        if !finite(&[a, b]) || !t.is_finite() {
            return Err(ImplicitError::NonFinite);
        }
        Ok(Point {
            kind: Kind::Lnc { a, b, t },
            w_negative: false,
        })
    }

    /// `(u·a + v·b + w·c) / (u + v + w)`, the point of triangle `abc`'s plane
    /// with barycentric weights `u`, `v`, `w`.
    ///
    /// # Errors
    /// [`ImplicitError::NonFinite`], or [`ImplicitError::ZeroWeight`] when
    /// the weights sum to exactly zero.
    pub fn bpt(tri: [Point3; 3], weights: [f64; 3]) -> Result<Point, ImplicitError> {
        if !finite(&tri) || !weights.iter().all(|w| w.is_finite()) {
            return Err(ImplicitError::NonFinite);
        }
        Point::with_w(Kind::Bpt { tri, weights }, ImplicitError::ZeroWeight)
    }

    /// Where the line through `line[0]` and `line[1]` meets the plane through
    /// `plane`.
    ///
    /// # Errors
    /// [`ImplicitError::NonFinite`], or
    /// [`ImplicitError::NoLinePlaneIntersection`] when the line is exactly
    /// parallel to the plane (or lies in it), or the line or the plane is
    /// degenerate.
    pub fn lpi(line: [Point3; 2], plane: [Point3; 3]) -> Result<Point, ImplicitError> {
        if !finite(&line) || !finite(&plane) {
            return Err(ImplicitError::NonFinite);
        }
        Point::with_w(
            Kind::Lpi { line, plane },
            ImplicitError::NoLinePlaneIntersection,
        )
    }

    /// The common point of three planes, each through three points.
    ///
    /// # Errors
    /// [`ImplicitError::NonFinite`], or
    /// [`ImplicitError::NoThreePlaneIntersection`] when the planes' normals
    /// are exactly linearly dependent (or a plane is degenerate).
    pub fn tpi(planes: [[Point3; 3]; 3]) -> Result<Point, ImplicitError> {
        if !planes.iter().all(|p| finite(p)) {
            return Err(ImplicitError::NonFinite);
        }
        Point::with_w(
            Kind::Tpi { planes },
            ImplicitError::NoThreePlaneIntersection,
        )
    }

    /// Settle the exact sign of `W`; a zero `W` means the point does not exist.
    fn with_w(kind: Kind, zero: ImplicitError) -> Result<Point, ImplicitError> {
        let o = kind.origin();
        let (sign, _) = settle(|stage| match stage {
            0 => homog::<Filtered>(&kind, o)[3].sign(),
            1 => homog::<Interval>(&kind, o)[3].sign(),
            _ => homog::<Exact>(&kind, o)[3].sign(),
        });
        match sign {
            Sign::Zero => Err(zero),
            s => Ok(Point {
                kind,
                w_negative: s == Sign::Negative,
            }),
        }
    }

    /// What defines the point.
    #[must_use]
    pub fn kind(&self) -> &Kind {
        &self.kind
    }

    /// Is the point explicit?
    #[must_use]
    pub fn is_explicit(&self) -> bool {
        matches!(self.kind, Kind::Explicit(_))
    }

    /// The f64 position, within a few ulps of the exact one (for spatial
    /// queries and output only, never a topological decision). The f64
    /// evaluation is kept when its error bound puts it within 2^-42 of the
    /// configuration's size; otherwise (a line nearly parallel to its plane,
    /// planes nearly sharing a line) the exact coordinates are divided.
    #[must_use]
    pub fn approx(&self) -> Point3 {
        if let Kind::Explicit(p) = self.kind {
            return p;
        }
        let o = self.kind.origin();
        let h = homog::<Filtered>(&self.kind, o);
        let (w, we) = (h[3].value(), h[3].error());
        if w.is_finite() && we.is_finite() && w.abs() > 2.0 * we {
            let size = self.kind.size();
            let fast = [0, 1, 2].map(|k| {
                let x = h[k].value() / w;
                let err = (h[k].error() + x.abs() * we) / (w.abs() - we);
                (
                    x,
                    err <= libm::ldexp(1.0, -42) * size.max(x.abs()) && err.is_finite(),
                )
            });
            if fast.iter().all(|&(x, ok)| ok && x.is_finite()) {
                return [0, 1, 2].map(|k| o[k] + fast[k].0);
            }
        }
        let h = homog::<Exact>(&self.kind, o);
        [0, 1, 2].map(|k| o[k] + h[k].quotient(&h[3]).unwrap_or(f64::NAN))
    }

    /// A box certainly holding the point: `[lo, hi]` per axis, from interval
    /// arithmetic rounded outwards (unbounded where the intervals are too
    /// wide to say). For spatial filters that must never miss the point.
    #[must_use]
    pub fn bounds(&self) -> [[f64; 2]; 3] {
        if let Kind::Explicit(p) = self.kind {
            return p.map(|c| [c, c]);
        }
        let o = self.kind.origin();
        let h = homog::<Interval>(&self.kind, o);
        [0, 1, 2].map(|k| {
            h[k].offset_quotient(&h[3], o[k])
                .unwrap_or([f64::NEG_INFINITY, f64::INFINITY])
        })
    }
}

impl From<Point3> for Point {
    /// An explicit point, unchecked: predicates on a non-finite point are
    /// meaningless (use [`Point::explicit`] to check).
    fn from(p: Point3) -> Point {
        Point {
            kind: Kind::Explicit(p),
            w_negative: false,
        }
    }
}

impl Kind {
    /// The largest distance (per axis) of a defining point from the origin:
    /// the size of the configuration.
    fn size(&self) -> f64 {
        let o = self.origin();
        let far = |ps: &[Point3]| {
            ps.iter()
                .flat_map(|p| (0..3).map(move |k| (p[k] - o[k]).abs()))
                .fold(0.0, f64::max)
        };
        match *self {
            Kind::Explicit(_) => 0.0,
            Kind::Lnc { a, b, .. } => far(&[a, b]),
            Kind::Bpt { tri, .. } => far(&tri),
            Kind::Lpi { line, plane } => far(&[line[0], line[1], plane[0], plane[1], plane[2]]),
            Kind::Tpi { planes } => far(&planes.concat()),
        }
    }

    /// An explicit point the coordinates are translated by: one of the
    /// defining points.
    fn origin(&self) -> Point3 {
        match *self {
            Kind::Explicit(p) => p,
            Kind::Lnc { b, .. } => b,
            Kind::Bpt { tri, .. } => tri[0],
            Kind::Lpi { line, .. } => line[0],
            Kind::Tpi { planes } => planes[0][0],
        }
    }
}

/// Evaluate a sign in stages: `f(0)` with the filter, `f(1)` with intervals,
/// `f(2)` exactly.
fn settle(f: impl Fn(u8) -> Option<Sign>) -> (Sign, Stage) {
    if let Some(s) = f(0) {
        return (s, Stage::Filter);
    }
    if let Some(s) = f(1) {
        return (s, Stage::Interval);
    }
    (f(2).unwrap_or(Sign::Zero), Stage::Exact)
}

/// `p - o` in the ring, exactly as far as the stage allows.
fn rel<R: Ring>(p: Point3, o: Point3) -> [R; 3] {
    [0, 1, 2].map(|k| R::from_f64(p[k]).sub(&R::from_f64(o[k])))
}

fn sub3<R: Ring>(a: &[R; 3], b: &[R; 3]) -> [R; 3] {
    [0, 1, 2].map(|k| a[k].sub(&b[k]))
}

fn dot3<R: Ring>(a: &[R; 3], b: &[R; 3]) -> R {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn cross3<R: Ring>(a: &[R; 3], b: &[R; 3]) -> [R; 3] {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

/// The normal `(b−a)×(c−a)` of a plane and `n·a`, translated by `o`.
fn plane<R: Ring>(p: &[Point3; 3], o: Point3) -> ([R; 3], R) {
    let [a, b, c] = p.map(|q| rel::<R>(q, o));
    let n = cross3(&sub3(&b, &a), &sub3(&c, &a));
    let d = dot3(&n, &a);
    (n, d)
}

/// Homogeneous coordinates `(X, Y, Z, W)` of the point translated by `o`:
/// the point is `o + (X, Y, Z) / W`.
fn homog<R: Ring>(kind: &Kind, o: Point3) -> [R; 4] {
    let one = || R::from_f64(1.0);
    match *kind {
        Kind::Explicit(p) => {
            let [x, y, z] = rel::<R>(p, o);
            [x, y, z, one()]
        }
        Kind::Lnc { a, b, t } => {
            // b + t·(a − b)
            let (a, b) = (rel::<R>(a, o), rel::<R>(b, o));
            let t = R::from_f64(t);
            let [x, y, z] = [0, 1, 2].map(|k| b[k].add(&t.mul(&a[k].sub(&b[k]))));
            [x, y, z, one()]
        }
        Kind::Bpt { tri, weights } => {
            let p = tri.map(|q| rel::<R>(q, o));
            let w = weights.map(R::from_f64);
            let [x, y, z] = [0, 1, 2].map(|k| {
                w[0].mul(&p[0][k])
                    .add(&w[1].mul(&p[1][k]))
                    .add(&w[2].mul(&p[2][k]))
            });
            [x, y, z, w[0].add(&w[1]).add(&w[2])]
        }
        Kind::Lpi { line, plane: pl } => {
            // p + λ(q − p), λ = n·(r − p) / n·(q − p).
            let p = rel::<R>(line[0], o);
            let q = rel::<R>(line[1], o);
            let (n, d) = plane::<R>(&pl, o);
            let qp = sub3(&q, &p);
            let den = dot3(&n, &qp);
            let num = d.sub(&dot3(&n, &p));
            let [x, y, z] = [0, 1, 2].map(|k| p[k].mul(&den).add(&qp[k].mul(&num)));
            [x, y, z, den]
        }
        Kind::Tpi { planes } => {
            // Cramer's rule on N·x = d, rows of N the planes' normals.
            let rows = planes.map(|p| plane::<R>(&p, o));
            let matrix = |replace: Option<usize>| -> Vec<Vec<R>> {
                rows.iter()
                    .map(|(n, d)| {
                        (0..3)
                            .map(|k| {
                                if replace == Some(k) {
                                    d.clone()
                                } else {
                                    n[k].clone()
                                }
                            })
                            .collect()
                    })
                    .collect()
            };
            let w = det(&matrix(None));
            let [x, y, z] = [0, 1, 2].map(|k| det(&matrix(Some(k))));
            [x, y, z, w]
        }
    }
}

/// The origin of a query: the first point's.
fn origin_of(points: &[&Point]) -> Point3 {
    points.first().map_or([0.0; 3], |p| p.kind.origin())
}

fn w_sign(points: &[&Point]) -> Sign {
    if points.iter().filter(|p| p.w_negative).count() % 2 == 1 {
        Sign::Negative
    } else {
        Sign::Positive
    }
}

fn times(s: Sign, t: Sign) -> Sign {
    match (s, t) {
        (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
        (a, b) if a == b => Sign::Positive,
        _ => Sign::Negative,
    }
}

/// `det` of the rows `(X, Y, Z, W)`; its sign times the signs of the `W`s is
/// the sign of `det[(x, y, z, 1)]`, which is minus `orient3d`.
fn orient3d_det<R: Ring>(p: [&Point; 4], o: Point3) -> R {
    let rows: Vec<Vec<R>> = p.iter().map(|q| homog::<R>(&q.kind, o).to_vec()).collect();
    det(&rows)
}

/// Exact orientation of four points, `Positive` for a positive-volume tet
/// (the convention of [`predicates::orient3d`]), and the stage that settled it.
#[must_use]
pub fn orient3d_staged(a: &Point, b: &Point, c: &Point, d: &Point) -> (Sign, Stage) {
    let p = [a, b, c, d];
    if let [
        Kind::Explicit(a),
        Kind::Explicit(b),
        Kind::Explicit(c),
        Kind::Explicit(d),
    ] = p.map(|q| q.kind)
    {
        return (predicates::orient3d(a, b, c, d), Stage::Explicit);
    }
    let o = origin_of(&p);
    let (s, stage) = settle(|stage| match stage {
        0 => orient3d_det::<Filtered>(p, o).sign(),
        1 => orient3d_det::<Interval>(p, o).sign(),
        _ => orient3d_det::<Exact>(p, o).sign(),
    });
    (times(s, w_sign(&p)).flipped(), stage)
}

/// Exact orientation of four points (see [`orient3d_staged`]).
#[must_use]
pub fn orient3d(a: &Point, b: &Point, c: &Point, d: &Point) -> Sign {
    orient3d_staged(a, b, c, d).0
}

/// `det` of the rows `(X·W, Y·W, Z·W, X² + Y² + Z², W²)`: the lifted rows
/// `(x, y, z, |x|², 1)` scaled by `W² > 0`, so its sign is that of the
/// classic in-sphere determinant, which is minus [`insphere`].
fn insphere_det<R: Ring>(p: [&Point; 5], o: Point3) -> R {
    let rows: Vec<Vec<R>> = p
        .iter()
        .map(|q| {
            let [x, y, z, w] = homog::<R>(&q.kind, o);
            let lift = x.mul(&x).add(&y.mul(&y)).add(&z.mul(&z));
            vec![x.mul(&w), y.mul(&w), z.mul(&w), lift, w.mul(&w)]
        })
        .collect();
    det(&rows)
}

/// Exact in-sphere test: `Positive` when `e` is strictly inside the sphere
/// through `a`, `b`, `c`, `d` and `orient3d(a, b, c, d)` is `Positive`; the
/// sign reverses for a negatively oriented `abcd`; `Zero` when cospherical
/// (as [`predicates::insphere`]).
#[must_use]
pub fn insphere_staged(a: &Point, b: &Point, c: &Point, d: &Point, e: &Point) -> (Sign, Stage) {
    let p = [a, b, c, d, e];
    if let [
        Kind::Explicit(a),
        Kind::Explicit(b),
        Kind::Explicit(c),
        Kind::Explicit(d),
        Kind::Explicit(e),
    ] = p.map(|q| q.kind)
    {
        return (predicates::insphere(a, b, c, d, e), Stage::Explicit);
    }
    let o = origin_of(&p);
    let (s, stage) = settle(|stage| match stage {
        0 => insphere_det::<Filtered>(p, o).sign(),
        1 => insphere_det::<Interval>(p, o).sign(),
        _ => insphere_det::<Exact>(p, o).sign(),
    });
    (s.flipped(), stage)
}

/// Exact in-sphere test (see [`insphere_staged`]).
#[must_use]
pub fn insphere(a: &Point, b: &Point, c: &Point, d: &Point, e: &Point) -> Sign {
    insphere_staged(a, b, c, d, e).0
}

/// The two coordinates kept when dropping `axis`, in the cyclic order of
/// [`predicates::project`].
fn kept(axis: usize) -> [usize; 2] {
    [(axis + 1) % 3, (axis + 2) % 3]
}

fn orient2d_det<R: Ring>(p: [&Point; 3], o: Point3, axis: usize) -> R {
    let [i, j] = kept(axis);
    let rows: Vec<Vec<R>> = p
        .iter()
        .map(|q| {
            let h = homog::<R>(&q.kind, o);
            vec![h[i].clone(), h[j].clone(), h[3].clone()]
        })
        .collect();
    det(&rows)
}

/// Exact orientation of three points projected along `axis` (the sign of
/// component `axis` of the normal `(b−a)×(c−a)`, as
/// [`predicates::orient2d_projected`]), and the stage that settled it.
#[must_use]
pub fn orient2d_projected_staged(a: &Point, b: &Point, c: &Point, axis: usize) -> (Sign, Stage) {
    let p = [a, b, c];
    let axis = axis % 3;
    if let [Kind::Explicit(a), Kind::Explicit(b), Kind::Explicit(c)] = p.map(|q| q.kind) {
        return (
            predicates::orient2d_projected(a, b, c, axis),
            Stage::Explicit,
        );
    }
    let o = origin_of(&p);
    let (s, stage) = settle(|stage| match stage {
        0 => orient2d_det::<Filtered>(p, o, axis).sign(),
        1 => orient2d_det::<Interval>(p, o, axis).sign(),
        _ => orient2d_det::<Exact>(p, o, axis).sign(),
    });
    (times(s, w_sign(&p)), stage)
}

/// Exact projected orientation (see [`orient2d_projected_staged`]).
#[must_use]
pub fn orient2d_projected(a: &Point, b: &Point, c: &Point, axis: usize) -> Sign {
    orient2d_projected_staged(a, b, c, axis).0
}

fn incircle_det<R: Ring>(p: [&Point; 4], o: Point3, axis: usize) -> R {
    let [i, j] = kept(axis);
    let rows: Vec<Vec<R>> = p
        .iter()
        .map(|q| {
            let h = homog::<R>(&q.kind, o);
            let (u, v, w) = (&h[i], &h[j], &h[3]);
            vec![u.mul(w), v.mul(w), u.mul(u).add(&v.mul(v)), w.mul(w)]
        })
        .collect();
    det(&rows)
}

/// Exact in-circle test on the projection along `axis`: `Positive` when the
/// projection of `d` is inside the circle through those of `a`, `b`, `c`
/// turning counter-clockwise (as [`predicates::incircle`] on
/// [`predicates::project`]ed points). A projection keeps orientation but
/// maps circles to ellipses: this is the consistent test for triangulating a
/// facet in its projection, not the circumcircle in 3D (for that, see
/// [`incircle_coplanar`]).
#[must_use]
pub fn incircle_projected_staged(
    a: &Point,
    b: &Point,
    c: &Point,
    d: &Point,
    axis: usize,
) -> (Sign, Stage) {
    let p = [a, b, c, d];
    let axis = axis % 3;
    if let [
        Kind::Explicit(a),
        Kind::Explicit(b),
        Kind::Explicit(c),
        Kind::Explicit(d),
    ] = p.map(|q| q.kind)
    {
        let pr = |q| predicates::project(q, axis);
        return (
            predicates::incircle(pr(a), pr(b), pr(c), pr(d)),
            Stage::Explicit,
        );
    }
    let o = origin_of(&p);
    let (s, stage) = settle(|stage| match stage {
        0 => incircle_det::<Filtered>(p, o, axis).sign(),
        1 => incircle_det::<Interval>(p, o, axis).sign(),
        _ => incircle_det::<Exact>(p, o, axis).sign(),
    });
    (s, stage)
}

/// Exact projected in-circle test (see [`incircle_projected_staged`]).
#[must_use]
pub fn incircle_projected(a: &Point, b: &Point, c: &Point, d: &Point, axis: usize) -> Sign {
    incircle_projected_staged(a, b, c, d, axis).0
}

fn coincident_det<R: Ring>(a: &Point, b: &Point, o: Point3, k: usize) -> R {
    let ha = homog::<R>(&a.kind, o);
    let hb = homog::<R>(&b.kind, o);
    ha[k].mul(&hb[3]).sub(&hb[k].mul(&ha[3]))
}

/// Are `a` and `b` the same point, exactly? (`X_a·W_b = X_b·W_a` for each
/// coordinate.)
#[must_use]
pub fn coincident(a: &Point, b: &Point) -> bool {
    if let (Kind::Explicit(p), Kind::Explicit(q)) = (a.kind, b.kind) {
        return p == q;
    }
    let o = origin_of(&[a, b]);
    (0..3).all(|k| {
        let (s, _) = settle(|stage| match stage {
            0 => coincident_det::<Filtered>(a, b, o, k).sign(),
            1 => coincident_det::<Interval>(a, b, o, k).sign(),
            _ => coincident_det::<Exact>(a, b, o, k).sign(),
        });
        s == Sign::Zero
    })
}

/// Are `a`, `b` and `c` on one line (or coincident), exactly?
#[must_use]
pub fn collinear(a: &Point, b: &Point, c: &Point) -> bool {
    (0..3).all(|axis| orient2d_projected(a, b, c, axis) == Sign::Zero)
}

/// Is `d`, coplanar with `a`, `b`, `c`, inside their circumcircle in 3D?
/// `Positive` inside, `Negative` outside, `Zero` on the circle or when `abc`
/// is collinear. Exact: every sphere through `a`, `b`, `c` meets their plane
/// in their circumcircle, so this is the in-sphere test with an explicit
/// point `q` exactly off the plane (as
/// [`predicates::in_circumcircle_coplanar`]).
#[must_use]
pub fn incircle_coplanar(a: &Point, b: &Point, c: &Point, d: &Point) -> Sign {
    // An axis along which the normal is exactly nonzero.
    let Some(axis) = (0..3).find(|&k| orient2d_projected(a, b, c, k) != Sign::Zero) else {
        return Sign::Zero;
    };
    // `q`: off the plane along `axis`, far enough to be clear of rounding.
    let pa = a.approx();
    let span = [b.approx(), c.approx()]
        .iter()
        .map(|p| (0..3).map(|k| (p[k] - pa[k]).abs()).fold(0.0, f64::max))
        .fold(0.0, f64::max);
    let mut step = if span > 0.0 && span.is_finite() {
        span
    } else {
        1.0
    };
    for _ in 0..64 {
        let mut q = pa;
        q[axis] += step;
        let q = Point::from(q);
        let side = orient3d(a, b, c, &q);
        if side != Sign::Zero {
            return times(insphere(a, b, c, &q, d), side);
        }
        step *= 2.0;
    }
    Sign::Zero
}

/// The sign of `W` for one point (`W` is positive unless the construction
/// made it negative).
fn w_of(p: &Point) -> Sign {
    if p.w_negative {
        Sign::Negative
    } else {
        Sign::Positive
    }
}

/// The numerator of `p − q` in homogeneous coordinates: `X_p·W_q − X_q·W_p`
/// (the difference is this over `W_p·W_q`).
fn diff_num<R: Ring>(p: &[R; 4], q: &[R; 4]) -> [R; 3] {
    [0, 1, 2].map(|k| p[k].mul(&q[3]).sub(&q[k].mul(&p[3])))
}

/// The numerator of `(p − q)·(r − q)`, over `W_p·W_r·W_q²`.
fn dot_det<R: Ring>(p: [&Point; 3], o: Point3) -> R {
    let [hp, hq, hr] = p.map(|x| homog::<R>(&x.kind, o));
    dot3(&diff_num(&hp, &hq), &diff_num(&hr, &hq))
}

/// The sign of the dot product `(p − q)·(r − q)`: `Positive` when the angle
/// at `q` is acute, `Zero` when it is right (or `p` or `r` is `q`), and the
/// stage that settled it.
#[must_use]
pub fn dot_sign_staged(p: &Point, q: &Point, r: &Point) -> (Sign, Stage) {
    let pts = [p, q, r];
    let o = origin_of(&[q]);
    let (s, stage) = settle(|stage| match stage {
        0 => dot_det::<Filtered>(pts, o).sign(),
        1 => dot_det::<Interval>(pts, o).sign(),
        _ => dot_det::<Exact>(pts, o).sign(),
    });
    (times(s, times(w_of(p), w_of(r))), stage)
}

/// The sign of `(p − q)·(r − q)` (see [`dot_sign_staged`]).
#[must_use]
pub fn dot_sign(p: &Point, q: &Point, r: &Point) -> Sign {
    dot_sign_staged(p, q, r).0
}

/// The numerator of `[(b − d)×(c − d)]·[(a − d)×(c − d)]`, over
/// `W_a·W_b·W_c²·W_d⁴`.
fn dihedral_det<R: Ring>(p: [&Point; 4], o: Point3) -> R {
    let [ha, hb, hc, hd] = p.map(|x| homog::<R>(&x.kind, o));
    let (ua, ub, uc) = (diff_num(&ha, &hd), diff_num(&hb, &hd), diff_num(&hc, &hd));
    dot3(&cross3(&ub, &uc), &cross3(&ua, &uc))
}

/// The sign of `[(b − d)×(c − d)]·[(a − d)×(c − d)]`: `Positive` when the
/// triangles `a c d` and `b c d`, sharing the edge `c d`, meet at an acute
/// dihedral angle, `Zero` at a right angle or when either is degenerate.
#[must_use]
pub fn dihedral_sign_staged(a: &Point, b: &Point, c: &Point, d: &Point) -> (Sign, Stage) {
    let pts = [a, b, c, d];
    let o = origin_of(&[d]);
    let (s, stage) = settle(|stage| match stage {
        0 => dihedral_det::<Filtered>(pts, o).sign(),
        1 => dihedral_det::<Interval>(pts, o).sign(),
        _ => dihedral_det::<Exact>(pts, o).sign(),
    });
    (times(s, times(w_of(a), w_of(b))), stage)
}

/// The sign of the dihedral expression (see [`dihedral_sign_staged`]).
#[must_use]
pub fn dihedral_sign(a: &Point, b: &Point, c: &Point, d: &Point) -> Sign {
    dihedral_sign_staged(a, b, c, d).0
}

/// Does the projection of edge `a q` onto the plane of triangle `a b c`
/// meet the triangle's interior? It does when `q` leans over the triangle
/// past both of its sides at `a`: the dihedral angles at `a b` (between
/// `q a b` and `c a b`) and at `a c` (between `q a c` and `b a c`) are both
/// acute. Exact.
#[must_use]
pub fn edge_projects_into_triangle(a: &Point, q: &Point, b: &Point, c: &Point) -> bool {
    dihedral_sign(q, c, a, b) == Sign::Positive && dihedral_sign(q, b, a, c) == Sign::Positive
}
