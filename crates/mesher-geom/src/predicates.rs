//! Sign-returning wrappers around Shewchuk's adaptive predicates (`robust` crate).

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
}

/// Exact 2D orientation: `Positive` when `a`, `b`, `c` turn counter-clockwise.
#[must_use]
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Sign {
    let coord = |p: [f64; 2]| robust::Coord { x: p[0], y: p[1] };
    Sign::of(robust::orient2d(coord(a), coord(b), coord(c)))
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn nan_is_degenerate_not_panic() {
        assert_eq!(Sign::of(f64::NAN), Sign::Zero);
    }
}
