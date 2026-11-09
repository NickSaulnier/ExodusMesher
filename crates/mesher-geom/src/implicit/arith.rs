//! The three arithmetic stages behind the implicit predicates. Each predicate
//! is written once, generically over [`Ring`], and evaluated with
//! [`Filtered`] first, then [`Interval`], then [`Exact`], stopping at the
//! first stage that can certify the sign.

use num_bigint::BigInt;
use num_bigint::Sign as BigSign;

use crate::Sign;

/// The arithmetic a predicate needs: exact inputs, ring operations, and the
/// sign of a result when the stage can certify it.
pub(crate) trait Ring: Clone {
    /// An f64 input, held exactly.
    fn from_f64(x: f64) -> Self;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    /// The sign of the exact value, if this stage can prove it.
    fn sign(&self) -> Option<Sign>;
}

/// Unit roundoff of f64 (round to nearest), doubled for slack.
const U2: f64 = f64::EPSILON;
/// Inflation that covers the rounding of the bound's own computation.
const GROW: f64 = 1.0 + 16.0 * f64::EPSILON;
/// Covers the absolute error of an operation whose result underflows.
const ETA: f64 = f64::MIN_POSITIVE;

/// An f64 value with a rigorous bound on its distance from the exact value:
/// `|v - V| <= e`. A bound of zero means the value is exact; results of
/// exact inputs that round without error stay exact, so degenerate inputs
/// (exact zeros) are often settled here too.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Filtered {
    v: f64,
    e: f64,
}

impl Filtered {
    /// The approximate value.
    pub(crate) fn value(&self) -> f64 {
        self.v
    }
}

/// `a + b` and its rounding error, exactly (Knuth's TwoSum).
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

impl Ring for Filtered {
    fn from_f64(x: f64) -> Filtered {
        Filtered { v: x, e: 0.0 }
    }

    fn add(&self, o: &Filtered) -> Filtered {
        if self.e == 0.0 && o.e == 0.0 {
            let (s, err) = two_sum(self.v, o.v);
            if err == 0.0 && s.is_finite() {
                return Filtered { v: s, e: 0.0 };
            }
        }
        let s = self.v + o.v;
        Filtered {
            v: s,
            e: (self.e + o.e + U2 * s.abs() + ETA) * GROW,
        }
    }

    fn sub(&self, o: &Filtered) -> Filtered {
        self.add(&Filtered { v: -o.v, e: o.e })
    }

    fn mul(&self, o: &Filtered) -> Filtered {
        let p = self.v * o.v;
        if self.e == 0.0 && o.e == 0.0 {
            // The product is exact when the fused remainder is zero and
            // nothing underflowed.
            let err = self.v.mul_add(o.v, -p);
            if err == 0.0 && p.is_finite() && (p == 0.0) == (self.v == 0.0 || o.v == 0.0) {
                let tiny = p != 0.0 && p.abs() < f64::MIN_POSITIVE;
                if !tiny {
                    return Filtered { v: p, e: 0.0 };
                }
            }
        }
        Filtered {
            v: p,
            e: (self.v.abs() * o.e + o.v.abs() * self.e + self.e * o.e + U2 * p.abs() + ETA) * GROW,
        }
    }

    fn sign(&self) -> Option<Sign> {
        if !(self.v.is_finite() && self.e.is_finite()) {
            return None;
        }
        if self.e == 0.0 {
            return Some(Sign::of(self.v));
        }
        if self.v > self.e {
            Some(Sign::Positive)
        } else if self.v < -self.e {
            Some(Sign::Negative)
        } else {
            None
        }
    }
}

/// The next f64 up (towards +∞); NaN and +∞ stay.
fn next_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let bits = x.to_bits();
    f64::from_bits(if x > 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(x: f64) -> f64 {
    -next_up(-x)
}

/// A closed interval containing the exact value, with every bound rounded
/// outwards.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    const WHOLE: Interval = Interval {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    fn checked(lo: f64, hi: f64) -> Interval {
        if lo.is_finite() && hi.is_finite() {
            Interval { lo, hi }
        } else {
            Interval::WHOLE
        }
    }
}

impl Ring for Interval {
    fn from_f64(x: f64) -> Interval {
        Interval::checked(x, x)
    }

    fn add(&self, o: &Interval) -> Interval {
        Interval::checked(next_down(self.lo + o.lo), next_up(self.hi + o.hi))
    }

    fn sub(&self, o: &Interval) -> Interval {
        Interval::checked(next_down(self.lo - o.hi), next_up(self.hi - o.lo))
    }

    fn mul(&self, o: &Interval) -> Interval {
        if !(self.lo.is_finite() && o.lo.is_finite() && self.hi.is_finite() && o.hi.is_finite()) {
            return Interval::WHOLE;
        }
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        let lo = p.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Interval::checked(next_down(lo), next_up(hi))
    }

    fn sign(&self) -> Option<Sign> {
        if self.lo > 0.0 {
            Some(Sign::Positive)
        } else if self.hi < 0.0 {
            Some(Sign::Negative)
        } else if self.lo == 0.0 && self.hi == 0.0 {
            Some(Sign::Zero)
        } else {
            None
        }
    }
}

/// An exact dyadic rational `m · 2^e`. Every f64 is one, and sums and
/// products of them stay exact with no overflow or underflow, whatever the
/// degree of the polynomial.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Exact {
    m: BigInt,
    e: i64,
}

impl Exact {
    fn zero() -> Exact {
        Exact {
            m: BigInt::ZERO,
            e: 0,
        }
    }

    fn is_zero(&self) -> bool {
        self.m.sign() == BigSign::NoSign
    }
}

impl Ring for Exact {
    fn from_f64(x: f64) -> Exact {
        if x == 0.0 || !x.is_finite() {
            return Exact::zero();
        }
        let bits = x.to_bits();
        let negative = bits >> 63 == 1;
        let exponent = ((bits >> 52) & 0x7ff) as i64;
        let fraction = bits & ((1u64 << 52) - 1);
        let (mut mantissa, mut e) = if exponent == 0 {
            (fraction, -1074)
        } else {
            (fraction | (1u64 << 52), exponent - 1075)
        };
        let tz = mantissa.trailing_zeros();
        mantissa >>= tz;
        e += i64::from(tz);
        let m = BigInt::from(mantissa);
        Exact {
            m: if negative { -m } else { m },
            e,
        }
    }

    fn add(&self, o: &Exact) -> Exact {
        if self.is_zero() {
            return o.clone();
        }
        if o.is_zero() {
            return self.clone();
        }
        let e = self.e.min(o.e);
        // Exponents stay within a few thousand of each other for f64 inputs
        // and the polynomial degrees used here.
        let shift = |x: &Exact| -> BigInt { x.m.clone() << ((x.e - e) as usize) };
        Exact {
            m: shift(self) + shift(o),
            e,
        }
    }

    fn sub(&self, o: &Exact) -> Exact {
        self.add(&Exact {
            m: -o.m.clone(),
            e: o.e,
        })
    }

    fn mul(&self, o: &Exact) -> Exact {
        if self.is_zero() || o.is_zero() {
            return Exact::zero();
        }
        Exact {
            m: &self.m * &o.m,
            e: self.e + o.e,
        }
    }

    fn sign(&self) -> Option<Sign> {
        Some(match self.m.sign() {
            BigSign::Minus => Sign::Negative,
            BigSign::NoSign => Sign::Zero,
            BigSign::Plus => Sign::Positive,
        })
    }
}

/// Determinant of a 2×2 to 5×5 matrix by Laplace expansion along the first
/// row (rows of equal length).
pub(crate) fn det<R: Ring>(m: &[Vec<R>]) -> R {
    let n = m.len();
    if n == 1 {
        return m[0][0].clone();
    }
    if n == 2 {
        return m[0][0].mul(&m[1][1]).sub(&m[0][1].mul(&m[1][0]));
    }
    let mut total: Option<R> = None;
    for j in 0..n {
        let minor: Vec<Vec<R>> = m[1..]
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|&(k, _)| k != j)
                    .map(|(_, x)| x.clone())
                    .collect()
            })
            .collect();
        let term = m[0][j].mul(&det(&minor));
        total = Some(match total {
            None => term,
            Some(t) if j % 2 == 0 => t.add(&term),
            Some(t) => t.sub(&term),
        });
    }
    total.unwrap_or_else(|| R::from_f64(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exact(x: f64) -> Exact {
        Exact::from_f64(x)
    }

    #[test]
    fn exact_arithmetic_is_exact() {
        // (1 + 2^-60) - 1 = 2^-60, lost in f64, kept here.
        let tiny = 2f64.powi(-60);
        let d = exact(1.0).add(&exact(tiny)).sub(&exact(1.0));
        assert_eq!(d, exact(tiny));
        // Far beyond the f64 range in both directions.
        let (up, down) = (2f64.powi(1000), 2f64.powi(-1000));
        let big = exact(up).mul(&exact(up)).mul(&exact(down));
        assert_eq!(big.sub(&exact(up)).sign(), Some(Sign::Zero));
        let near = exact(1e300)
            .mul(&exact(1e300))
            .sub(&exact(1e300).mul(&exact(1e300)));
        assert_eq!(near.sign(), Some(Sign::Zero));
        let small = exact(5e-324).mul(&exact(5e-324));
        assert_eq!(small.sign(), Some(Sign::Positive));
        assert_eq!(exact(-0.0).sign(), Some(Sign::Zero));
    }

    #[test]
    fn the_filter_bound_holds_and_proves_exact_zeros() {
        let a = Filtered::from_f64(0.1);
        let b = Filtered::from_f64(0.2);
        let s = a.add(&b).sub(&Filtered::from_f64(0.3));
        // 0.1 + 0.2 - 0.3 is 2^-54 or so exactly: the bound must cover it.
        let exact_value = exact(0.1).add(&exact(0.2)).sub(&exact(0.3));
        assert_eq!(s.sign().unwrap_or(Sign::Positive), Sign::Positive);
        assert_eq!(exact_value.sign(), Some(Sign::Positive));
        // Exact zeros of exact inputs are certain at the filter.
        let z = Filtered::from_f64(3.0)
            .mul(&Filtered::from_f64(4.0))
            .sub(&Filtered::from_f64(12.0));
        assert_eq!(z.sign(), Some(Sign::Zero));
        // An underflowing product is not taken as exact.
        let u = Filtered::from_f64(1e-200).mul(&Filtered::from_f64(1e-200));
        assert_eq!(u.sign(), None);
    }

    #[test]
    fn intervals_round_outwards() {
        let third = Interval::from_f64(1.0).mul(&Interval::from_f64(1.0 / 3.0));
        assert!(third.lo < 1.0 / 3.0 && third.hi > 1.0 / 3.0);
        let z = Interval::from_f64(2.0).sub(&Interval::from_f64(2.0));
        assert_eq!(z.sign(), None, "rounded outwards, so not exactly zero");
        assert_eq!(Interval::from_f64(0.0).sign(), Some(Sign::Zero));
        assert_eq!(next_up(0.0), f64::from_bits(1));
        assert_eq!(next_down(1.0), 1.0 - f64::EPSILON / 2.0);
    }

    #[test]
    fn determinants_match_by_hand() {
        let m: Vec<Vec<Exact>> = [[2.0, 0.0, 1.0], [1.0, 3.0, 2.0], [1.0, 1.0, 1.0]]
            .iter()
            .map(|r| r.iter().map(|&x| exact(x)).collect())
            .collect();
        // 2(3-2) - 0 + 1(1-3) = 0
        assert_eq!(det(&m).sign(), Some(Sign::Zero));
    }
}
