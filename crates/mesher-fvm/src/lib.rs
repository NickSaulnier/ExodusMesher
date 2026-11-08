//! Cell-centred finite-volume Poisson solver, the downstream check that a mesh
//! is usable by a finite-volume code.
//!
//! Solves `-∇²φ = f` with Dirichlet values on the boundary, one unknown per
//! tetrahedron at its centroid. Across a face with area vector `S` and
//! centroid-to-centroid vector `d`, the diffusive flux `S·∇φ` is split as in
//! OpenFOAM's over-relaxed corrected Laplacian: an orthogonal part
//! `Δ = (S·S)/(S·d) d`, treated implicitly as `|Δ|/|d| (φN − φP)`, plus an
//! explicit non-orthogonal correction `(S − Δ)·∇φ_f` from least-squares cell
//! gradients, iterated as a deferred correction. With [`Scheme::TwoPoint`] the
//! correction is dropped: exact on orthogonal meshes, inconsistent otherwise.
//! Bad cells cost accuracy either way, which is the point of the check. Each
//! linear solve uses Jacobi-preconditioned conjugate gradients.

use mesher_geom::Point3;
use mesher_geom::vec3::{add, cross, dot, norm, scale, sub};

/// Faces of a positive tet, outward.
const OUTWARD: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scheme {
    /// Orthogonal part only.
    TwoPoint,
    /// With the non-orthogonal correction (the default).
    Corrected,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PoissonResult {
    pub cells: usize,
    /// Volume-weighted RMS error at cell centroids.
    pub l2_error: f64,
    pub max_error: f64,
    /// Conjugate-gradient iterations, summed over correction passes.
    pub iterations: usize,
    pub correction_passes: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FvmError {
    Empty,
    /// A cell has zero volume or a face zero area.
    Degenerate {
        cell: usize,
    },
    NotConverged {
        iterations: usize,
    },
}

impl std::fmt::Display for FvmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FvmError::Empty => write!(f, "no cells"),
            FvmError::Degenerate { cell } => write!(f, "cell {cell} is degenerate"),
            FvmError::NotConverged { iterations } => {
                write!(
                    f,
                    "conjugate gradients did not converge in {iterations} iterations"
                )
            }
        }
    }
}

/// A sparse symmetric matrix in compressed rows.
struct Csr {
    start: Vec<usize>,
    col: Vec<usize>,
    val: Vec<f64>,
}

impl Csr {
    fn from_triplets(n: usize, mut t: Vec<(usize, usize, f64)>) -> Csr {
        t.sort_by_key(|a| (a.0, a.1));
        let mut start = vec![0usize; n + 1];
        let mut col = Vec::with_capacity(t.len());
        let mut val: Vec<f64> = Vec::with_capacity(t.len());
        let mut last: Option<(usize, usize)> = None;
        for (r, c, v) in t {
            if last == Some((r, c)) {
                if let Some(x) = val.last_mut() {
                    *x += v;
                }
                continue;
            }
            last = Some((r, c));
            col.push(c);
            val.push(v);
            start[r + 1] += 1;
        }
        for r in 0..n {
            start[r + 1] += start[r];
        }
        Csr { start, col, val }
    }

    fn mul(&self, x: &[f64], y: &mut [f64]) {
        for (r, out) in y.iter_mut().enumerate() {
            *out = (self.start[r]..self.start[r + 1])
                .map(|k| self.val[k] * x[self.col[k]])
                .sum();
        }
    }

    fn diagonal(&self) -> Vec<f64> {
        (0..self.start.len() - 1)
            .map(|r| {
                (self.start[r]..self.start[r + 1])
                    .find(|&k| self.col[k] == r)
                    .map_or(1.0, |k| self.val[k])
            })
            .collect()
    }

    /// Jacobi-preconditioned conjugate gradients from `x`.
    fn solve(&self, b: &[f64], x: &mut [f64]) -> Result<usize, FvmError> {
        let n = b.len();
        let inv: Vec<f64> = self.diagonal().iter().map(|&d| 1.0 / d).collect();
        let mut ax = vec![0.0; n];
        self.mul(x, &mut ax);
        let mut r: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
        let mut z: Vec<f64> = r.iter().zip(&inv).map(|(a, b)| a * b).collect();
        let mut p = z.clone();
        let mut ap = vec![0.0f64; n];
        let bnorm = b
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE);
        let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let max_iterations = 20 * n + 100;
        let mut it = 0;
        while r.iter().map(|v| v * v).sum::<f64>().sqrt() / bnorm > 1e-12 {
            if it >= max_iterations {
                return Err(FvmError::NotConverged { iterations: it });
            }
            self.mul(&p, &mut ap);
            let alpha = rz / p.iter().zip(&ap).map(|(a, b)| a * b).sum::<f64>();
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
                z[i] = r[i] * inv[i];
            }
            let next: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
            let beta = next / rz;
            rz = next;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
            it += 1;
        }
        Ok(it)
    }
}

struct Face {
    owner: usize,
    /// `None` on the boundary.
    neighbour: Option<usize>,
    /// Area vector, outward from the owner.
    s: Point3,
    centroid: Point3,
}

/// Solve with the default (corrected) scheme.
pub fn solve_manufactured(
    points: &[Point3],
    tets: &[[u32; 4]],
    exact: &dyn Fn(Point3) -> f64,
    source: &dyn Fn(Point3) -> f64,
) -> Result<PoissonResult, FvmError> {
    solve_with(points, tets, exact, source, Scheme::Corrected)
}

/// Solve `-∇²φ = source` with `φ = exact` on the boundary and compare with
/// `exact` at the cell centroids.
pub fn solve_with(
    points: &[Point3],
    tets: &[[u32; 4]],
    exact: &dyn Fn(Point3) -> f64,
    source: &dyn Fn(Point3) -> f64,
    scheme: Scheme,
) -> Result<PoissonResult, FvmError> {
    let n = tets.len();
    if n == 0 {
        return Err(FvmError::Empty);
    }
    let corners = |t: &[u32; 4]| t.map(|v| points[v as usize]);
    let mean = |p: &[Point3]| {
        let k = p.len() as f64;
        [0, 1, 2].map(|i| p.iter().map(|q| q[i]).sum::<f64>() / k)
    };
    let centroids: Vec<Point3> = tets.iter().map(|t| mean(&corners(t))).collect();
    let volumes: Vec<f64> = tets
        .iter()
        .map(|t| {
            let p = corners(t);
            dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))) / 6.0
        })
        .collect();
    if let Some(cell) = volumes.iter().position(|&v| v <= 0.0 || !v.is_finite()) {
        return Err(FvmError::Degenerate { cell });
    }

    // Faces, interior in discovery order, then boundary sorted.
    let mut open: std::collections::HashMap<[u32; 3], (usize, usize)> =
        std::collections::HashMap::with_capacity(2 * n);
    let mut raw: Vec<(usize, Option<usize>, usize)> = Vec::with_capacity(2 * n);
    for (i, t) in tets.iter().enumerate() {
        for (k, f) in OUTWARD.iter().enumerate() {
            let mut key = f.map(|j| t[j]);
            key.sort_unstable();
            match open.remove(&key) {
                Some((j, kj)) => raw.push((j, Some(i), kj)),
                None => {
                    open.insert(key, (i, k));
                }
            }
        }
    }
    let mut boundary: Vec<(usize, usize)> = open.into_values().collect();
    boundary.sort_unstable();
    raw.extend(boundary.into_iter().map(|(i, k)| (i, None, k)));
    let mut faces = Vec::with_capacity(raw.len());
    for (owner, neighbour, k) in raw {
        let f = OUTWARD[k].map(|m| points[tets[owner][m] as usize]);
        let s = scale(cross(sub(f[1], f[0]), sub(f[2], f[0])), 0.5);
        if norm(s) <= 0.0 {
            return Err(FvmError::Degenerate { cell: owner });
        }
        faces.push(Face {
            owner,
            neighbour,
            s,
            centroid: mean(&f),
        });
    }

    // Orthogonal parts (implicit) and corrections (explicit).
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(7 * n);
    let mut rhs: Vec<f64> = (0..n).map(|i| source(centroids[i]) * volumes[i]).collect();
    let mut correction: Vec<Point3> = Vec::with_capacity(faces.len());
    for f in &faces {
        let other = f.neighbour.map_or(f.centroid, |q| centroids[q]);
        let d = sub(other, centroids[f.owner]);
        let sd = dot(f.s, d).max(1e-6 * norm(f.s) * norm(d));
        let delta = scale(d, dot(f.s, f.s) / sd);
        let a = norm(delta) / norm(d);
        correction.push(sub(f.s, delta));
        match f.neighbour {
            Some(q) => {
                let p = f.owner;
                triplets.extend([(p, p, a), (q, q, a), (p, q, -a), (q, p, -a)]);
            }
            None => {
                triplets.push((f.owner, f.owner, a));
                rhs[f.owner] += a * exact(f.centroid);
            }
        }
    }
    let m = Csr::from_triplets(n, triplets);

    let mut x = vec![0.0f64; n];
    let mut iterations = m.solve(&rhs, &mut x)?;
    let mut passes = 0;
    if scheme == Scheme::Corrected {
        for _ in 0..60 {
            passes += 1;
            let grad = gradients(&faces, &centroids, &x, exact);
            let mut b = rhs.clone();
            for (f, k) in faces.iter().zip(&correction) {
                let g = match f.neighbour {
                    Some(q) => scale(add(grad[f.owner], grad[q]), 0.5),
                    None => grad[f.owner],
                };
                let flux = dot(*k, g);
                b[f.owner] += flux;
                if let Some(q) = f.neighbour {
                    b[q] -= flux;
                }
            }
            let before = x.clone();
            iterations += m.solve(&b, &mut x)?;
            let change = x
                .iter()
                .zip(&before)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            let size = x.iter().map(|v| v.abs()).fold(0.0, f64::max).max(1e-300);
            if change <= 1e-10 * size {
                break;
            }
        }
    }

    let total: f64 = volumes.iter().sum();
    let mut sq = 0.0;
    let mut max_error = 0.0f64;
    for i in 0..n {
        let e = (x[i] - exact(centroids[i])).abs();
        sq += volumes[i] * e * e;
        max_error = max_error.max(e);
    }
    Ok(PoissonResult {
        cells: n,
        l2_error: (sq / total).sqrt(),
        max_error,
        iterations,
        correction_passes: passes,
    })
}

/// Least-squares cell gradients from face neighbours and boundary values.
fn gradients(
    faces: &[Face],
    centroids: &[Point3],
    x: &[f64],
    exact: &dyn Fn(Point3) -> f64,
) -> Vec<Point3> {
    let n = centroids.len();
    let mut ata = vec![[0.0f64; 6]; n];
    let mut atb = vec![[0.0f64; 3]; n];
    let mut add_pair = |i: usize, d: Point3, dv: f64| {
        let w = 1.0 / dot(d, d).max(f64::MIN_POSITIVE);
        let m = &mut ata[i];
        m[0] += w * d[0] * d[0];
        m[1] += w * d[0] * d[1];
        m[2] += w * d[0] * d[2];
        m[3] += w * d[1] * d[1];
        m[4] += w * d[1] * d[2];
        m[5] += w * d[2] * d[2];
        for k in 0..3 {
            atb[i][k] += w * d[k] * dv;
        }
    };
    for f in faces {
        let p = f.owner;
        match f.neighbour {
            Some(q) => {
                let d = sub(centroids[q], centroids[p]);
                add_pair(p, d, x[q] - x[p]);
                add_pair(q, scale(d, -1.0), x[p] - x[q]);
            }
            None => {
                let d = sub(f.centroid, centroids[p]);
                add_pair(p, d, exact(f.centroid) - x[p]);
            }
        }
    }
    (0..n)
        .map(|i| {
            let [a, b, c, d, e, f] = ata[i];
            let m = [[a, b, c], [b, d, e], [c, e, f]];
            solve3(m, atb[i]).unwrap_or([0.0; 3])
        })
        .collect()
}

fn solve3(m: [[f64; 3]; 3], r: [f64; 3]) -> Option<[f64; 3]> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    if d.abs() <= 1e-300 {
        return None;
    }
    Some([0, 1, 2].map(|k| {
        let mut mk = m;
        for row in 0..3 {
            mk[row][k] = r[row];
        }
        det(mk) / d
    }))
}
