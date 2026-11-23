//! Feature and patch survival on real models (Milestone 18, stage 5).
//!
//!     cargo run --release -p mesher-synthetic --example survival -- <obj files>
//!
//! Each file goes through the pipeline at the default settings (in this
//! process, without the worker's time limit); for each meshed part it prints
//! the kernel, whether the surface was remeshed, whether the mesh was
//! accepted, its tets and largest neighbour size ratio, the input's feature length,
//! the largest distance from an input feature point to the mesh's feature
//! edges (in element sizes, capped at `survival::SEARCH`), the share of
//! feature length within a tenth of a size, the worst patch area error, the
//! worst point, and the smallest dihedral angle and largest non-orthogonality.
//! With `SURVIVAL_PATCHES` set, each patch's input and side-set areas follow
//! (in square element sizes), on lines starting `#patch`.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::print_stdout)]

use std::path::Path;

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::repair::repair;
use mesher_core::survival::survival;
use mesher_core::volume::mesh_parts;
use mesher_core::{JobConfig, RecordLog};
use mesher_io::read_obj;

fn main() {
    let config = JobConfig::default();
    println!(
        "model\tpart\tkernel\tremeshed\taccepted\ttets\tmax_size_ratio\tfeature_length\tmax_distance\tkept\tpatches\tworst_patch_error\tworst_point\tmin_dihedral\tmax_non_orthogonality"
    );
    for arg in std::env::args().skip(1) {
        let path = Path::new(&arg);
        let id = path
            .file_stem()
            .map_or(arg.clone(), |s| s.to_string_lossy().into_owned());
        let mut log = RecordLog::new();
        let Some(obj) = read_obj(path, &mut log, &mut |_| {}) else {
            println!("{id}\t-\tunreadable");
            continue;
        };
        let model = obj;
        let Ok(diagnosis) = diagnose(&model, &Tolerances::default(), &mut log) else {
            println!("{id}\t-\tdiagnosis failed");
            continue;
        };
        let Some(selected) = select_parts(&diagnosis, &config, &mut log) else {
            println!("{id}\t-\tno parts");
            continue;
        };
        let mut overlay = diagnosis.overlay.clone();
        if check_part_pairs(&model, &diagnosis, &selected, &mut log, &mut overlay).is_err() {
            println!("{id}\t-\tpart check failed");
            continue;
        }
        let Ok(repaired) = repair(
            &model,
            &diagnosis,
            &selected,
            &Tolerances::default(),
            &mut log,
        ) else {
            println!("{id}\t-\trepair failed");
            continue;
        };
        let Ok(outcome) = mesh_parts(&model, &repaired, &config, &mut log) else {
            println!("{id}\t-\tmeshing failed");
            continue;
        };
        if outcome.meshed.is_empty() {
            println!("{id}\t-\tno mesh");
        }
        for m in &outcome.meshed {
            let Some((_, summary)) = outcome.summaries.iter().find(|(i, _)| *i == m.part) else {
                continue;
            };
            let Some(r) = outcome
                .wrapped
                .iter()
                .find(|w| w.part == m.part)
                .or_else(|| repaired.iter().find(|r| r.part == m.part))
            else {
                continue;
            };
            let s = survival(&model, r, m, &config, summary.element_size);
            let patch = s
                .patches
                .iter()
                .map(|(_, a, b)| if *a > 0.0 { (a - b).abs() / a } else { 0.0 })
                .fold(0.0, f64::max);
            println!(
                "{id}\t{}\t{:?}\t{}\t{}\t{}\t{:.3}\t{:.4}\t{:.4}\t{:.4}\t{}\t{:.4}\t{:?}\t{:.2}\t{:.2}",
                m.name,
                m.kernel,
                summary.remeshed,
                summary.accepted,
                summary.tets,
                summary.max_size_ratio,
                s.input_length,
                s.max_distance,
                s.kept,
                s.patches.len(),
                patch,
                s.worst_point,
                summary.min_dihedral_degrees,
                summary.max_non_orthogonality_degrees
            );
            if std::env::var_os("SURVIVAL_PATCHES").is_some() {
                let h2 = summary.element_size * summary.element_size;
                for (name, a, b) in &s.patches {
                    println!(
                        "#patch\t{id}\t{}\t{name}\t{:.4}\t{:.4}",
                        m.name,
                        a / h2,
                        b / h2
                    );
                }
            }
        }
    }
}
