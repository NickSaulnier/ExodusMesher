//! Sign-returning wrappers around Shewchuk's adaptive predicates (`robust` crate),
//! plus exact helpers built only from them.
//!
//! Conventions (differ from Shewchuk's raw functions, which the wrappers flip):
//! - [`orient3d`] is `Positive` when `abcd` is a positive-volume tetrahedron, i.e.
//!   `d` lies on the side of plane `abc` that the normal `(b-a)×(c-a)` points to.
//! - [`insphere`] is `Positive` when `e` is inside the circumsphere of a
//!   positively oriented `abcd`.

use crate::Point3;

/// The sign of an exact predicate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Sign {
    Negative,
    Zero,
    Positive,
}

impl Sign {
    /// Sign of a predicate value. A NaN (only possible from non-finite input) is
    /// treated as `Zero`, i.e. degenerate, rather than panicking.
    #[must_use]
    pub fn of(value: f64) -> Sign {
        if value > 0.0 {
            Sign::Positive
        } else if value < 0.0 {
            Sign::Negative
        } else {
            Sign::Zero
        }
    }

    #[must_use]
    pub fn flipped(self) -> Sign {
        match self {
            Sign::Negative => Sign::Positive,
            Sign::Zero => Sign::Zero,
            Sign::Positive => Sign::Negative,
        }
    }

    /// True for `Positive` against `Negative` (or the reverse); false if either is `Zero`.
    #[must_use]
    pub fn opposes(self, other: Sign) -> bool {
        matches!(
            (self, other),
            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
        )
    }
}

fn coord2(p: [f64; 2]) -> robust::Coord<f64> {
    robust::Coord { x: p[0], y: p[1] }
}

fn coord3(p: Point3) -> robust::Coord3D<f64> {
    robust::Coord3D {
        x: p[0],
        y: p[1],
        z: p[2],
    }
}

/// Exact 2D orientation: `Positive` when `a`, `b`, `c` turn counter-clockwise.
#[must_use]
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Sign {
    Sign::of(robust::orient2d(coord2(a), coord2(b), coord2(c)))
}

/// Exact in-circle test: `Positive` when `d` is strictly inside the circle
/// through `a`, `b`, `c` given counter-clockwise; the sign reverses for
/// clockwise `a`, `b`, `c`; `Zero` when concyclic.
#[must_use]
pub fn incircle(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> Sign {
    Sign::of(robust::incircle(coord2(a), coord2(b), coord2(c), coord2(d)))
}

/// Is `d` strictly inside the circumcircle of triangle `abc`, all four
/// coplanar in 3D? `false` for a degenerate triangle.
///
/// Exact without projecting (a projection onto a coordinate plane keeps
/// orientations but turns circles into ellipses): every sphere through `a`,
/// `b`, `c` meets their plane in exactly their circumcircle, so `d` is inside
/// the circle exactly when it is inside the sphere through `a`, `b`, `c` and
/// any point `q` off the plane. `q` is `a` moved along an axis on which the
/// plane's normal is exactly nonzero, so it is exactly off the plane.
#[must_use]
pub fn in_circumcircle_coplanar(a: Point3, b: Point3, c: Point3, d: Point3) -> bool {
    let Some(k) = plane_axis(a, b, c) else {
        return false;
    };
    let mut q = a;
    q[k] = if a[k] + 1.0 != a[k] {
        a[k] + 1.0
    } else {
        a[k] * 2.0
    };
    let side = orient3d(a, b, c, q);
    let inside = insphere(a, b, c, q, d);
    side != Sign::Zero && inside != Sign::Zero && inside == side
}

/// Exact 3D orientation: the sign of `det[b-a, c-a, d-a]`, `Positive` for a
/// positive-volume tetrahedron `abcd`.
#[must_use]
pub fn orient3d(a: Point3, b: Point3, c: Point3, d: Point3) -> Sign {
    // Shewchuk's orient3d is positive when d is *below* plane abc: the opposite sign.
    Sign::of(robust::orient3d(coord3(a), coord3(b), coord3(c), coord3(d))).flipped()
}

/// Exact in-sphere test: `Positive` when `e` is strictly inside the circumsphere
/// of `abcd` and `orient3d(a, b, c, d)` is `Positive`. The sign reverses for a
/// negatively oriented `abcd`; `Zero` means cospherical (or a flat `abcd`).
#[must_use]
pub fn insphere(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> Sign {
    // Shewchuk's insphere is positive-inside for his positive orientation, which is
    // our negative one.
    Sign::of(robust::insphere(
        coord3(a),
        coord3(b),
        coord3(c),
        coord3(d),
        coord3(e),
    ))
    .flipped()
}

/// Drop coordinate `axis`, keeping the other two in cyclic order
/// (`axis+1`, `axis+2`), so orientation in the projection equals the sign of
/// the normal's `axis` component.
#[must_use]
pub fn project(p: Point3, axis: usize) -> [f64; 2] {
    [p[(axis + 1) % 3], p[(axis + 2) % 3]]
}

/// Exact sign of component `axis` of the normal `(b-a)×(c-a)`.
#[must_use]
pub fn orient2d_projected(a: Point3, b: Point3, c: Point3, axis: usize) -> Sign {
    orient2d(project(a, axis), project(b, axis), project(c, axis))
}

/// Exactly collinear (including coincident) points.
#[must_use]
pub fn collinear(a: Point3, b: Point3, c: Point3) -> bool {
    (0..3).all(|axis| orient2d_projected(a, b, c, axis) == Sign::Zero)
}

/// A coordinate axis along which the normal of triangle `abc` has an exactly
/// nonzero component, preferring the largest; projecting along it maps the
/// triangle's plane one-to-one, so coplanar orientation tests can run in 2D.
/// `None` for a degenerate (collinear) triangle.
#[must_use]
pub fn plane_axis(a: Point3, b: Point3, c: Point3) -> Option<usize> {
    let n = crate::vec3::cross(crate::vec3::sub(b, a), crate::vec3::sub(c, a));
    let mut axes = [0, 1, 2];
    axes.sort_by(|&i, &j| n[j].abs().total_cmp(&n[i].abs()).then(i.cmp(&j)));
    axes.into_iter()
        .find(|&axis| orient2d_projected(a, b, c, axis) != Sign::Zero)
}

#[cfg(test)]
mod tests {
    use super::*;

    const O: Point3 = [0.0, 0.0, 0.0];
    const X: Point3 = [1.0, 0.0, 0.0];
    const Y: Point3 = [0.0, 1.0, 0.0];
    const Z: Point3 = [0.0, 0.0, 1.0];

    #[test]
    fn orient2d_signs() {
        assert_eq!(orient2d([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]), Sign::Positive);
        assert_eq!(orient2d([0.0, 0.0], [0.0, 1.0], [1.0, 0.0]), Sign::Negative);
        assert_eq!(orient2d([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]), Sign::Zero);
    }

    #[test]
    fn orient2d_is_exact_near_degeneracy() {
        // A point a hair off the line y = x, where naive evaluation rounds to zero.
        let a = [0.5, 0.5];
        let b = [12.0, 12.0];
        let c = [24.0, 24.0 + 2f64.powi(-48)];
        assert_eq!(orient2d(a, b, c), Sign::Positive);
    }

    #[test]
    fn orient3d_positive_volume_convention() {
        assert_eq!(orient3d(O, X, Y, Z), Sign::Positive);
        assert_eq!(orient3d(O, Y, X, Z), Sign::Negative);
        assert_eq!(orient3d(O, X, Y, [3.0, -2.0, 0.0]), Sign::Zero);
    }

    #[test]
    fn insphere_convention() {
        let inside = [0.25, 0.25, 0.25];
        let outside = [5.0, 5.0, 5.0];
        let on = [1.0, 1.0, 0.0]; // on the sphere through O, X, Y, Z
        assert_eq!(insphere(O, X, Y, Z, inside), Sign::Positive);
        assert_eq!(insphere(O, X, Y, Z, outside), Sign::Negative);
        assert_eq!(insphere(O, X, Y, Z, on), Sign::Zero);
        assert_eq!(insphere(O, Y, X, Z, inside), Sign::Negative);
    }

    #[test]
    fn plane_axis_and_collinear() {
        assert_eq!(plane_axis(O, X, Y), Some(2));
        assert_eq!(plane_axis(O, Y, Z), Some(0));
        assert_eq!(plane_axis(O, X, [2.0, 0.0, 0.0]), None);
        assert!(collinear(O, X, [-7.0, 0.0, 0.0]));
        assert!(collinear(O, O, O));
        assert!(!collinear(O, X, Y));
    }

    #[test]
    fn incircle_and_coplanar_circumcircle() {
        let (a, b, c) = ([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]);
        assert_eq!(incircle(a, b, c, [0.4, 0.4]), Sign::Positive);
        assert_eq!(incircle(a, b, c, [2.0, 2.0]), Sign::Negative);
        assert_eq!(incircle(a, b, c, [1.0, 1.0]), Sign::Zero);
        assert_eq!(incircle(a, c, b, [0.4, 0.4]), Sign::Negative);
        // The same in a tilted plane z = x.
        let lift = |p: [f64; 2]| [p[0], p[1], p[0]];
        assert!(in_circumcircle_coplanar(
            lift(a),
            lift(b),
            lift(c),
            lift([0.4, 0.4])
        ));
        assert!(in_circumcircle_coplanar(
            lift(a),
            lift(c),
            lift(b),
            lift([0.4, 0.4])
        ));
        assert!(!in_circumcircle_coplanar(
            lift(a),
            lift(b),
            lift(c),
            lift([1.0, 1.0])
        ));
        // A plane where projecting onto a coordinate plane distorts circles:
        // x + y + z = 3. (4, 1, 3) are on a circle about (1, 1, 1) of radius
        // sqrt(2); (-1, 2, 2) is outside it, (1, 1, 1) inside.
        let (p4, p1, p3) = ([2.0, 0.0, 1.0], [1.0, 2.0, 0.0], [1.0, 0.0, 2.0]);
        assert!(!in_circumcircle_coplanar(p4, p1, p3, [-1.0, 2.0, 2.0]));
        assert!(in_circumcircle_coplanar(p4, p1, p3, [1.0, 1.0, 1.0]));
        assert!(in_circumcircle_coplanar(p4, p3, p1, [1.0, 1.0, 1.0]));
        // Far from the origin, where a + 1 would round.
        let big = |p: Point3| p.map(|x| x + 1e17);
        assert!(!in_circumcircle_coplanar(
            big(p4),
            big(p1),
            big(p3),
            big([-1.0, 2.0, 2.0])
        ));
    }

    #[test]
    fn nan_is_degenerate_not_panic() {
        assert_eq!(Sign::of(f64::NAN), Sign::Zero);
        let nan = [f64::NAN, 0.0, 0.0];
        let _ = orient3d(nan, X, Y, Z);
        let _ = insphere(nan, O, X, Y, Z);
    }
}
