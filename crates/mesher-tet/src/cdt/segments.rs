//! Segment recovery: make every edge of the input triangles an edge of the
//! tetrahedralization by inserting exact points on it.
//!
//! Every input triangle is its own facet, so every triangle edge is a
//! segment. A missing subsegment is split by an LNC point of its input
//! segment (`t·a + (1−t)·b`, exactly on the segment). Where exactly one end
//! of the subsegment is an input vertex, the split is on a sphere around that
//! vertex whose radius is a power of two (concentric shells, after Ruppert and
//! Shewchuk): segments meeting there at a small angle split on common
//! spheres, so they do not chase each other forever. Otherwise the split is
//! at the middle. Lengths are approximate (f64 positions): they only choose
//! where to split, never decide topology.
//!
//! Each round first tries to recover every missing subsegment by flips
//! ([`ExactDelaunay::recover_edge`]), then splits what is still missing and
//! tries flips on the two halves. Every subsegment that is an edge is
//! protected: flips never remove it, and split points go in by
//! [`ExactDelaunay::insert_protected`], whose cavity never removes it, so a
//! recovered segment stays recovered. The tetrahedralization is then no
//! longer Delaunay. Rounds repeat until every subsegment is an edge, within a
//! point budget.
//!
//! On the 100-model clean subset this recovers every model, at a median of
//! 125 Steiner points; splitting alone recovered 77 at a median of about
//! 4,000. Needle triangles cost the most: 101548, with angles down to
//! 0.006°, needs 70,000.

use std::collections::BTreeSet;

use mesher_geom::Point3;
use mesher_geom::implicit::Point;

use super::delaunay::{DelaunayError, ExactDelaunay, Inserted};

/// A piece of an input segment between two vertices.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subsegment {
    /// The input segment: its end vertices, lower index first.
    pub segment: (u32, u32),
    /// The vertices at each end.
    pub ends: [u32; 2],
    /// The LNC parameters of the ends on the segment (1 at `segment.0`,
    /// 0 at `segment.1`).
    pub t: [f64; 2],
}

/// Why segment recovery failed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SegmentError {
    Delaunay(DelaunayError),
    /// A non-finite input coordinate.
    NonFinite {
        vertex: u32,
    },
    /// A split point equals an existing vertex: a vertex lies on a segment,
    /// so the input is not a valid PLC.
    VertexOnSegment {
        segment: (u32, u32),
        vertex: u32,
    },
    /// A subsegment is too short to split in f64 parameters.
    TooShort {
        segment: (u32, u32),
    },
    /// More points than the budget.
    Budget {
        points: usize,
        limit: usize,
    },
}

impl std::fmt::Display for SegmentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SegmentError::Delaunay(e) => write!(f, "{e}"),
            SegmentError::NonFinite { vertex } => write!(f, "vertex {vertex} is not finite"),
            SegmentError::VertexOnSegment { segment, vertex } => write!(
                f,
                "vertex {vertex} lies on segment {}-{}",
                segment.0, segment.1
            ),
            SegmentError::TooShort { segment } => write!(
                f,
                "a piece of segment {}-{} is too short to split",
                segment.0, segment.1
            ),
            SegmentError::Budget { points, limit } => {
                write!(f, "{points} points exceed the budget of {limit}")
            }
        }
    }
}

impl From<DelaunayError> for SegmentError {
    fn from(e: DelaunayError) -> SegmentError {
        SegmentError::Delaunay(e)
    }
}

/// A Delaunay tetrahedralization in which every input segment is a chain of
/// edges.
#[derive(Debug, Clone)]
pub struct RecoveredSegments {
    pub dt: ExactDelaunay,
    /// Every subsegment, sorted by segment and parameter.
    pub subsegments: Vec<Subsegment>,
    /// Steiner points inserted on segments.
    pub steiner: usize,
    /// Subsegments recovered by flips.
    pub by_flips: usize,
    /// Segments missing from the first tetrahedralization of the inputs.
    pub initially_missing: usize,
    /// Rounds of splitting.
    pub rounds: usize,
}

/// The unique edges of `triangles`, lower vertex first, sorted.
#[must_use]
pub fn segments_of(triangles: &[[u32; 3]]) -> Vec<(u32, u32)> {
    let set: BTreeSet<(u32, u32)> = triangles
        .iter()
        .flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])])
        .filter(|(u, v)| u != v)
        .map(|(u, v)| (u.min(v), u.max(v)))
        .collect();
    set.into_iter().collect()
}

/// Flip steps tried per missing subsegment before splitting it.
const FLIP_STEPS: usize = 64;

fn distance(a: Point3, b: Point3) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt()
}

/// Tetrahedralize `points` and recover every edge of `triangles`. `max_points`
/// bounds the vertex count (0: 20 × the inputs + 10,000).
///
/// # Errors
/// See [`SegmentError`].
pub fn recover_segments(
    points: &[Point3],
    triangles: &[[u32; 3]],
    max_points: usize,
    seed: u64,
) -> Result<RecoveredSegments, SegmentError> {
    let explicit: Vec<Point> = points
        .iter()
        .enumerate()
        .map(|(i, &p)| Point::explicit(p).map_err(|_| SegmentError::NonFinite { vertex: i as u32 }))
        .collect::<Result<_, _>>()?;
    let mut dt = ExactDelaunay::new(explicit, seed)?;
    if let Some(&(v, _)) = dt.duplicates().first() {
        // A duplicate input vertex would make segments ambiguous.
        return Err(SegmentError::VertexOnSegment {
            segment: (v, v),
            vertex: v,
        });
    }
    let limit = if max_points == 0 {
        20 * points.len() + 10_000
    } else {
        max_points
    };
    // Shells: radii are powers of two of one unit for the whole input.
    let (lo, hi) = points.iter().fold(
        ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
        |(lo, hi), p| {
            (
                [0, 1, 2].map(|k| lo[k].min(p[k])),
                [0, 1, 2].map(|k| hi[k].max(p[k])),
            )
        },
    );
    let diagonal = distance(lo, hi);
    let unit = if diagonal > 0.0 && diagonal.is_finite() {
        libm::exp2(libm::log2(diagonal).floor())
    } else {
        1.0
    };
    let inputs = points.len() as u32;
    let mut subs: Vec<Subsegment> = segments_of(triangles)
        .into_iter()
        .map(|s| Subsegment {
            segment: s,
            ends: [s.0, s.1],
            t: [1.0, 0.0],
        })
        .collect();
    // Whether each subsegment is known to be an edge: recovered ones are
    // protected, so they stay edges and need no checking again.
    let mut done: Vec<bool> = vec![false; subs.len()];
    let mut steiner = 0usize;
    let mut by_flips = 0usize;
    let mut initially_missing = None;
    let mut rounds = 0usize;
    let key = |u: u32, v: u32| (u.min(v), u.max(v));
    loop {
        let missing: Vec<usize> = {
            let _s = crate::prof::scope("cdt segments: find missing");
            let mut missing = Vec::new();
            for k in 0..subs.len() {
                if !done[k] {
                    if dt.has_edge(subs[k].ends[0], subs[k].ends[1]) {
                        done[k] = true;
                    } else {
                        missing.push(k);
                    }
                }
            }
            missing
        };
        initially_missing.get_or_insert(missing.len());
        if missing.is_empty() {
            break;
        }
        rounds += 1;
        // Every subsegment that is an edge is protected: no flip or
        // insertion below removes it.
        let mut protected: BTreeSet<(u32, u32)> = subs
            .iter()
            .zip(&done)
            .filter(|&(_, &d)| d)
            .map(|(s, _)| key(s.ends[0], s.ends[1]))
            .collect();
        for &k in &missing {
            let [u, v] = subs[k].ends;
            let _s = crate::prof::scope("cdt segments: flips");
            if dt.recover_edge(u, v, &|x, y| protected.contains(&key(x, y)), FLIP_STEPS) {
                by_flips += 1;
                protected.insert(key(u, v));
            }
        }
        // Split what flips could not recover, from the back, so earlier
        // indices stay valid.
        for &k in missing.iter().rev() {
            let s = subs[k];
            if dt.has_edge(s.ends[0], s.ends[1]) {
                continue;
            }
            if dt.points().len() >= limit {
                return Err(SegmentError::Budget {
                    points: dt.points().len(),
                    limit,
                });
            }
            let (a, b) = (points[s.segment.0 as usize], points[s.segment.1 as usize]);
            let length = distance(a, b);
            let t = split_parameter(s, length, unit, inputs);
            let (t_hi, t_lo) = (s.t[0].max(s.t[1]), s.t[0].min(s.t[1]));
            if !(t > t_lo && t < t_hi) {
                return Err(SegmentError::TooShort { segment: s.segment });
            }
            let p =
                Point::lnc(a, b, t).map_err(|_| SegmentError::TooShort { segment: s.segment })?;
            let _s = crate::prof::scope("cdt segments: protected insertion");
            let m = match dt.insert_protected(p, &|x, y| protected.contains(&key(x, y)))? {
                Inserted::New(m) => m,
                Inserted::Duplicate(v) => {
                    return Err(SegmentError::VertexOnSegment {
                        segment: s.segment,
                        vertex: v,
                    });
                }
            };
            steiner += 1;
            subs[k] = Subsegment {
                ends: [s.ends[0], m],
                t: [s.t[0], t],
                ..s
            };
            subs.insert(
                k + 1,
                Subsegment {
                    ends: [m, s.ends[1]],
                    t: [t, s.t[1]],
                    ..s
                },
            );
            done[k] = false;
            done.insert(k + 1, false);
            for (u, v) in [(s.ends[0], m), (m, s.ends[1])] {
                let _s = crate::prof::scope("cdt segments: flips after a split");
                if dt.has_edge(u, v)
                    || dt.recover_edge(u, v, &|x, y| protected.contains(&key(x, y)), FLIP_STEPS)
                {
                    protected.insert(key(u, v));
                }
            }
        }
    }
    Ok(RecoveredSegments {
        dt,
        subsegments: subs,
        steiner,
        by_flips,
        initially_missing: initially_missing.unwrap_or(0),
        rounds,
    })
}

/// Where to split a subsegment, as an LNC parameter: on a power-of-two
/// sphere around its end that is an input vertex if exactly one is, else at
/// the middle.
fn split_parameter(s: Subsegment, length: f64, unit: f64, inputs: u32) -> f64 {
    let middle = 0.5 * (s.t[0] + s.t[1]);
    let at_input = [s.ends[0] < inputs, s.ends[1] < inputs];
    if at_input[0] == at_input[1] || !(length > 0.0 && length.is_finite()) {
        return middle;
    }
    // Distance from the input end, a power of two of the unit nearest half
    // the subsegment's length.
    let piece = (s.t[0] - s.t[1]).abs() * length;
    let d = unit * libm::exp2(libm::log2(piece / (2.0 * unit)).round());
    if !(d > 0.0 && d < piece) {
        return middle;
    }
    // t falls from 1 at the segment's first end to 0 at its second.
    let (t_end, toward) = if at_input[0] {
        (s.t[0], s.t[1])
    } else {
        (s.t[1], s.t[0])
    };
    let step = d / length;
    if toward < t_end {
        t_end - step
    } else {
        t_end + step
    }
}
