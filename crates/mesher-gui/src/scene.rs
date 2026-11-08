//! A loaded model: the OBJ read and diagnosed with the same code the worker
//! runs, so the parts, records and defect overlay the GUI shows are exactly
//! those the worker will mesh. Diagnosis only; meshing happens in the worker.

use std::path::{Path, PathBuf};

use mesher_core::diagnose::{Tolerances, diagnose};
use mesher_core::model::SurfaceModel;
use mesher_core::overlay::Overlay;
use mesher_core::{Record, RecordLog, Severity};

/// One part, as the worker numbers it (position + 1 is its index).
#[derive(Debug, Clone, PartialEq)]
pub struct PartInfo {
    pub name: String,
    /// Model triangle indices.
    pub triangles: Vec<u32>,
    /// `[min, max]`, zero when empty.
    pub bounds: [[f64; 3]; 2],
    /// Records of each severity that name this part: info, warning, error, fatal.
    pub badges: [usize; 4],
}

#[derive(Debug, Clone)]
pub struct Scene {
    pub path: PathBuf,
    pub model: SurfaceModel,
    pub parts: Vec<PartInfo>,
    /// Every record of reading and diagnosis, in order.
    pub records: Vec<Record>,
    pub overlay: Overlay,
    /// Part position of each model triangle (`u32::MAX` if in none).
    pub triangle_part: Vec<u32>,
    /// Bounds of the parts, `[min, max]`.
    pub bounds: [[f64; 3]; 2],
}

/// Position of a severity in [`PartInfo::badges`].
#[must_use]
pub fn severity_slot(s: Severity) -> usize {
    match s {
        Severity::Info => 0,
        Severity::Warning => 1,
        Severity::Error => 2,
        Severity::Fatal => 3,
    }
}

impl Scene {
    /// Read and diagnose an OBJ file. Never fails: problems are records (a
    /// file that cannot be read has no parts and a fatal record).
    #[must_use]
    pub fn load(path: &Path) -> Scene {
        let mut log = RecordLog::new();
        let model = mesher_io::read_obj(path, &mut log, &mut |_| {});
        let mut scene = Scene {
            path: path.to_path_buf(),
            model: model.clone().unwrap_or_default(),
            parts: Vec::new(),
            records: Vec::new(),
            overlay: Overlay::default(),
            triangle_part: Vec::new(),
            bounds: [[0.0; 3]; 2],
        };
        if let Some(model) = model.as_ref().filter(|_| !log.has_fatal()) {
            if let Ok(d) = diagnose(model, &Tolerances::default(), &mut log) {
                scene.triangle_part = vec![u32::MAX; model.triangles.len()];
                for (i, p) in d.parts.iter().enumerate() {
                    for &t in &p.triangles {
                        if let Some(slot) = scene.triangle_part.get_mut(t as usize) {
                            *slot = i as u32;
                        }
                    }
                    let bounds = if p.bounds.is_empty() {
                        [[0.0; 3]; 2]
                    } else {
                        [p.bounds.min, p.bounds.max]
                    };
                    scene.parts.push(PartInfo {
                        name: p.name.clone(),
                        triangles: p.triangles.clone(),
                        bounds,
                        badges: [0; 4],
                    });
                }
                scene.overlay = d.overlay;
            }
        }
        scene.records = log.records().to_vec();
        for r in &scene.records {
            if let Some(p) = scene
                .parts
                .iter_mut()
                .find(|p| r.part.as_deref() == Some(p.name.as_str()))
            {
                p.badges[severity_slot(r.severity)] += 1;
            }
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in scene.parts.iter().filter(|p| !p.triangles.is_empty()) {
            for k in 0..3 {
                lo[k] = lo[k].min(p.bounds[0][k]);
                hi[k] = hi[k].max(p.bounds[1][k]);
            }
        }
        if lo[0].is_finite() {
            scene.bounds = [lo, hi];
        }
        scene
    }

    /// The part the ray hits first (Möller–Trumbore against every triangle;
    /// a click is rare enough that no acceleration structure is needed).
    #[must_use]
    pub fn pick(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<usize> {
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let cross = |a: [f64; 3], b: [f64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let mut best: Option<(f64, usize)> = None;
        for (t, tri) in self.model.triangles.iter().enumerate() {
            let part = self.triangle_part.get(t).copied().unwrap_or(u32::MAX);
            if part == u32::MAX {
                continue;
            }
            let [a, b, c] = tri.vertices.map(|v| {
                self.model
                    .vertices
                    .get(v as usize)
                    .copied()
                    .unwrap_or([0.0; 3])
            });
            let (e1, e2) = (sub(b, a), sub(c, a));
            let pv = cross(dir, e2);
            let det = dot(e1, pv);
            if det.abs() < 1e-300 {
                continue;
            }
            let tv = sub(origin, a);
            let u = dot(tv, pv) / det;
            if !(0.0..=1.0).contains(&u) {
                continue;
            }
            let qv = cross(tv, e1);
            let v = dot(dir, qv) / det;
            if v < 0.0 || u + v > 1.0 {
                continue;
            }
            let dist = dot(e2, qv) / det;
            if dist > 0.0 && best.is_none_or(|(d, _)| dist < d) {
                best = Some((dist, part as usize));
            }
        }
        best.map(|(_, p)| p)
    }
}
