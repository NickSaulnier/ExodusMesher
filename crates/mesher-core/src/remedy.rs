//! Suggested remedies (Milestone 20 stage 5): for a record, the settings
//! changes worth trying next, in order, each with a one-line reason, in a
//! form a program can apply ([`Remedy::apply`], `Report.next_job()` in
//! Python).
//!
//! The table is grounded in what the pipeline already does: a remedy is a
//! step the fallback ladder did not take under the job's settings (the
//! exact kernel when the conforming one was chosen, remeshing when it was
//! turned off, the fallbacks when `strict` forbade them) or a size that
//! moves a budget or a distance the record measured the right way. A code
//! with no step known to help gets none: input problems (parsing, names,
//! units, capacity), repairs that already ran, kernel bugs, warnings about
//! output that is still valid. `no_repair` is never undone: diagnosis only
//! is a choice, not a failure.
//!
//! Remedies go on records that left a part without a mesh (errors and
//! fatal records), and on the conformal kernel's failures that a fallback
//! (stuffing, the stuffing-only retry, the wrap) covered, where the remedy
//! is a way to a conformal mesh. Neither another kernel nor remeshing is
//! suggested for a part whose repair split a non-manifold edge or vertex:
//! its surface then touches itself at coincident vertices, which no
//! Delaunay kernel can mesh (the synthetic suite's cases of it all fall to
//! stuffing whatever the kernel).

use std::collections::BTreeSet;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::report::{Action, DefectCode, JobConfig, KernelChoice, Outcome, Record, Severity};

/// The setting a remedy changes: a `config` field, by name.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Setting {
    /// `parts`.
    Parts,
    /// `strict`.
    Strict,
    /// `no_remesh`.
    NoRemesh,
    /// `kernel`.
    Kernel,
    /// `element_size`.
    ElementSize,
    /// `max_deviation`.
    MaxDeviation,
    /// `wrap_alpha`.
    WrapAlpha,
}

impl Setting {
    /// Every setting a remedy can change.
    pub const ALL: [Setting; 7] = [
        Setting::Parts,
        Setting::Strict,
        Setting::NoRemesh,
        Setting::Kernel,
        Setting::ElementSize,
        Setting::MaxDeviation,
        Setting::WrapAlpha,
    ];

    /// The [`JobConfig`] field's name.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Setting::Parts => "parts",
            Setting::Strict => "strict",
            Setting::NoRemesh => "no_remesh",
            Setting::Kernel => "kernel",
            Setting::ElementSize => "element_size",
            Setting::MaxDeviation => "max_deviation",
            Setting::WrapAlpha => "wrap_alpha",
        }
    }
}

/// One suggested next step: set `setting` to `value`, or multiply it by
/// `scale` (exactly one of the two).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Remedy {
    /// The setting to change.
    pub setting: Setting,
    /// Its new value, as in `config`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<Value>,
    /// A factor to multiply its value by (a numeric setting).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scale: Option<f64>,
    /// Why it should help, in one line.
    pub reason: String,
}

/// The finest element size a remedy suggests (relative to the diagonal).
pub const MIN_SUGGESTED_SIZE: f64 = 1.0 / 512.0;
/// The coarsest element size a remedy suggests.
pub const MAX_SUGGESTED_SIZE: f64 = 0.25;
/// How far above a measured distance a suggested `max_deviation` goes.
pub const DEVIATION_MARGIN: f64 = 1.25;

impl Remedy {
    fn set(setting: Setting, value: impl Into<Value>, reason: &str) -> Remedy {
        Remedy {
            setting,
            value: Some(value.into()),
            scale: None,
            reason: reason.to_owned(),
        }
    }

    fn scale(setting: Setting, scale: f64, reason: &str) -> Remedy {
        Remedy {
            setting,
            value: None,
            scale: Some(scale),
            reason: reason.to_owned(),
        }
    }

    /// `config` with this remedy applied.
    ///
    /// # Errors
    /// The remedy is not exactly one of a value and a scale, a scale is
    /// applied to a setting that is not a number, or the result is not a
    /// valid configuration (an unknown value, a size or distance that is not
    /// positive).
    pub fn apply(&self, config: &JobConfig) -> Result<JobConfig, String> {
        let name = self.setting.as_str();
        let mut v = serde_json::to_value(config).map_err(|e| e.to_string())?;
        let fields = v
            .as_object_mut()
            .ok_or_else(|| "the configuration is not an object".to_owned())?;
        let new = match (&self.value, self.scale) {
            (Some(value), None) => value.clone(),
            (None, Some(scale)) => {
                let old = fields
                    .get(name)
                    .and_then(Value::as_f64)
                    .ok_or_else(|| format!("`{name}` is not a number to scale"))?;
                serde_json::Number::from_f64(old * scale)
                    .map(Value::Number)
                    .ok_or_else(|| format!("`{name}` × {scale} is not finite"))?
            }
            _ => return Err("a remedy has exactly one of `value` and `scale`".to_owned()),
        };
        fields.insert(name.to_owned(), new);
        let out: JobConfig = serde_json::from_value(v).map_err(|e| format!("`{name}`: {e}"))?;
        for (what, x) in [
            ("element_size", out.element_size),
            ("max_deviation", out.max_deviation),
            ("wrap_alpha", out.wrap_alpha),
        ] {
            if !(x > 0.0 && x.is_finite()) {
                return Err(format!("`{what}` must be a positive number, not {x}"));
            }
        }
        Ok(out)
    }
}

/// Fill in the remedies of every record of a job run with `config`.
pub fn fill(records: &mut [Record], config: &JobConfig) {
    let coincident = coincident_parts(records);
    for r in records.iter_mut() {
        r.remedies = suggest(
            r,
            config,
            r.part.as_ref().is_some_and(|p| coincident.contains(p)),
        );
    }
}

/// The parts whose repair split a non-manifold edge or vertex, leaving
/// coincident vertices.
#[must_use]
pub fn coincident_parts(records: &[Record]) -> BTreeSet<String> {
    records
        .iter()
        .filter(|r| {
            matches!(r.action, Some(Action::SplitEdge | Action::DuplicateVertex))
                && r.outcome == Some(Outcome::Applied)
        })
        .filter_map(|r| r.part.clone())
        .collect()
}

/// The remedies for `record` in a job run with `config`, best first; empty
/// where none is known or none is needed. `coincident`: the record's part
/// has coincident vertices ([`coincident_parts`]).
#[must_use]
pub fn suggest(record: &Record, config: &JobConfig, coincident: bool) -> Vec<Remedy> {
    let unresolved = record.severity >= Severity::Error;
    let covered = record.severity == Severity::Warning
        && matches!(record.action, Some(Action::Fallback | Action::Wrap))
        && record.outcome == Some(Outcome::Applied);
    if !unresolved && !covered {
        return Vec::new();
    }
    let mut out = Vec::new();
    let exact_kernel = |out: &mut Vec<Remedy>| {
        if config.kernel == KernelChoice::Conformal && !coincident {
            out.push(Remedy::set(
                Setting::Kernel,
                "exact",
                "the exact constrained kernel recovers surfaces without refining them",
            ));
        }
    };
    let coarser = |out: &mut Vec<Remedy>, reason: &str| {
        if config.element_size * 2.0 <= MAX_SUGGESTED_SIZE {
            out.push(Remedy::scale(Setting::ElementSize, 2.0, reason));
        }
    };
    let measured = |key: &str| record.measured.get(key);
    let exact_kernel_record = measured("kernel").and_then(Value::as_str) == Some("exact");
    match record.code {
        DefectCode::PartEmptySelection if unresolved && !config.parts.is_empty() => {
            out.push(Remedy::set(
                Setting::Parts,
                Value::Array(Vec::new()),
                "mesh every part: the selection names a part the input does not have",
            ));
        }
        DefectCode::MeshStrictUnmet if unresolved => {
            exact_kernel(&mut out);
            if config.no_remesh && !coincident {
                out.push(Remedy::set(
                    Setting::NoRemesh,
                    false,
                    "remeshing the surface retries the conformal kernel on a better surface",
                ));
            }
            if config.strict {
                out.push(Remedy::set(
                    Setting::Strict,
                    false,
                    "allow the stuffing fallback, whose boundary only approximates the surface",
                ));
            }
        }
        DefectCode::MeshRecoveryFailed if !exact_kernel_record && !coincident => {
            if config.no_remesh {
                out.push(Remedy::set(
                    Setting::NoRemesh,
                    false,
                    "remeshing the surface retries recovery on a better surface",
                ));
            }
            exact_kernel(&mut out);
        }
        DefectCode::MeshBudgetExceeded if !exact_kernel_record => {
            coarser(&mut out, "a coarser element size needs fewer points");
        }
        DefectCode::MeshWorkerFailed => {
            let reason = measured("reason").and_then(Value::as_str);
            if matches!(reason, Some("timeout" | "out_of_memory")) {
                coarser(
                    &mut out,
                    "a coarser element size needs less time and memory",
                );
            }
        }
        DefectCode::ExoCapacityExceeded if unresolved => {
            coarser(&mut out, "a coarser element size makes fewer elements");
        }
        DefectCode::VolBoundaryDeviation if unresolved => {
            if config.element_size * 0.5 >= MIN_SUGGESTED_SIZE {
                out.push(Remedy::scale(
                    Setting::ElementSize,
                    0.5,
                    "a finer element size follows the surface more closely",
                ));
            }
            if let Some(d) = measured("hausdorff_relative").and_then(Value::as_f64) {
                accept(
                    &mut out,
                    d,
                    config,
                    "accept the boundary at the distance measured",
                );
            }
        }
        _ => {}
    }
    // A wrap that came out too far from the surface (repair rung 3), on
    // whatever record sent the part to it.
    if unresolved
        && record.action == Some(Action::Wrap)
        && record.outcome == Some(Outcome::Failed)
        && let Some(d) = measured("deviation").and_then(Value::as_f64)
    {
        out.push(Remedy::scale(
            Setting::WrapAlpha,
            0.5,
            "a smaller wrap closes fewer openings and stays closer to the surface",
        ));
        accept(
            &mut out,
            d,
            config,
            "accept the wrap at the distance measured",
        );
    }
    out
}

/// Raise `max_deviation` to just above a measured relative distance `d`.
fn accept(out: &mut Vec<Remedy>, d: f64, config: &JobConfig, reason: &str) {
    let value = d * DEVIATION_MARGIN;
    if value > config.max_deviation && value.is_finite() {
        out.push(Remedy::set(Setting::MaxDeviation, value, reason));
    }
}

/// The next job's settings: the first remedy that applies and changes
/// `config`, of the worst record (fatal first, then errors, ...; in report
/// order within a severity) at `min` severity or above.
#[derive(Debug, Clone, PartialEq)]
pub struct Next {
    /// The record's index in the report.
    pub record: usize,
    /// The remedy's index in the record's `remedies`.
    pub remedy: usize,
    pub config: JobConfig,
}

/// The next job's settings ([`Next`]); `None` when no record at `min`
/// severity or above has a remedy that applies.
#[must_use]
pub fn next_config(records: &[Record], config: &JobConfig, min: Severity) -> Option<Next> {
    // Fatal records first, then errors, ...: the worst problem is the one
    // to act on.
    let mut order: Vec<usize> = (0..records.len())
        .filter(|&i| records[i].severity >= min)
        .collect();
    order.sort_by_key(|&i| std::cmp::Reverse(records[i].severity));
    order.into_iter().find_map(|i| {
        records[i]
            .remedies
            .iter()
            .enumerate()
            .filter_map(|(k, r)| r.apply(config).ok().map(|c| (k, c)))
            .find(|(_, c)| c != config)
            .map(|(k, c)| Next {
                record: i,
                remedy: k,
                config: c,
            })
    })
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use crate::report::Stage;

    fn error(code: DefectCode) -> Record {
        Record::new(code, Stage::VolumeMesh).severity(Severity::Error)
    }

    fn settings(rs: &[Remedy]) -> Vec<&'static str> {
        rs.iter().map(|r| r.setting.as_str()).collect()
    }

    #[test]
    fn every_setting_is_a_config_field() {
        let schema = serde_json::to_value(schemars::schema_for!(JobConfig)).unwrap();
        let fields = schema["properties"].as_object().unwrap();
        for s in Setting::ALL {
            assert!(fields.contains_key(s.as_str()), "{}", s.as_str());
            assert_eq!(serde_json::to_value(s).unwrap(), s.as_str());
        }
    }

    /// Every remedy the table gives, for every code under a spread of
    /// settings, applies to a valid configuration that differs.
    #[test]
    fn every_suggested_remedy_applies() {
        let configs = [
            JobConfig::default(),
            JobConfig {
                parts: vec!["a".into()],
                strict: true,
                no_remesh: true,
                ..JobConfig::default()
            },
            JobConfig {
                kernel: KernelChoice::Exact,
                element_size: 0.2,
                ..JobConfig::default()
            },
            JobConfig {
                element_size: 0.003,
                max_deviation: 0.5,
                ..JobConfig::default()
            },
        ];
        let mut seen = 0;
        for config in &configs {
            for &code in DefectCode::ALL {
                for severity in [Severity::Warning, Severity::Error, Severity::Fatal] {
                    for (action, outcome) in [
                        (None, None),
                        (Some(Action::Fallback), Some(Outcome::Applied)),
                        (Some(Action::Wrap), Some(Outcome::Failed)),
                    ] {
                        let mut r = Record::new(code, Stage::VolumeMesh)
                            .severity(severity)
                            .measured("reason", "timeout")
                            .measured("hausdorff_relative", 0.08)
                            .measured("deviation", 0.07);
                        r.action = action;
                        r.outcome = outcome;
                        for remedy in suggest(&r, config, false) {
                            let applied = remedy.apply(config).unwrap();
                            assert_ne!(&applied, config, "{code:?}: {remedy:?}");
                            // The remedy survives the report's JSON.
                            let text = serde_json::to_string(&remedy).unwrap();
                            let back: Remedy = serde_json::from_str(&text).unwrap();
                            assert_eq!(back, remedy);
                            seen += 1;
                        }
                    }
                }
            }
        }
        assert!(seen > 20, "{seen}");
    }

    #[test]
    fn remedies_follow_the_settings_and_the_record() {
        let d = JobConfig::default();
        // Conformal failures: the exact kernel, unless it was the one.
        assert_eq!(
            settings(&suggest(&error(DefectCode::MeshRecoveryFailed), &d, false)),
            ["kernel"]
        );
        let exact = JobConfig {
            kernel: KernelChoice::Exact,
            ..d.clone()
        };
        assert!(suggest(&error(DefectCode::MeshRecoveryFailed), &exact, false).is_empty());
        let no_remesh = JobConfig {
            no_remesh: true,
            ..d.clone()
        };
        assert_eq!(
            settings(&suggest(
                &error(DefectCode::MeshRecoveryFailed),
                &no_remesh,
                false
            )),
            ["no_remesh", "kernel"]
        );
        // Covered by stuffing: still a way to a conformal mesh.
        let mut covered = error(DefectCode::MeshRecoveryFailed).severity(Severity::Warning);
        covered.action = Some(Action::Fallback);
        covered.outcome = Some(Outcome::Applied);
        assert_eq!(settings(&suggest(&covered, &d, false)), ["kernel"]);
        // Remeshed into a mesh: nothing to do.
        let mut remeshed = covered.clone();
        remeshed.action = Some(Action::Remesh);
        assert!(suggest(&remeshed, &d, false).is_empty());
        // Strict: the exact kernel, then the fallbacks.
        let strict = JobConfig {
            strict: true,
            ..d.clone()
        };
        assert_eq!(
            settings(&suggest(
                &error(DefectCode::MeshStrictUnmet),
                &strict,
                false
            )),
            ["kernel", "strict"]
        );
        // Budgets: coarser, except for the exact kernel's, which does not
        // refine.
        let budget = error(DefectCode::MeshBudgetExceeded);
        assert_eq!(settings(&suggest(&budget, &d, false)), ["element_size"]);
        assert!(suggest(&budget.clone().measured("kernel", "exact"), &d, false).is_empty());
        let coarse = JobConfig {
            element_size: 0.2,
            ..d.clone()
        };
        assert!(suggest(&budget, &coarse, false).is_empty());
        // The worker: time and memory, never a crash.
        let worker = |reason: &str| {
            Record::new(DefectCode::MeshWorkerFailed, Stage::Worker)
                .severity(Severity::Error)
                .measured("reason", reason)
        };
        assert_eq!(
            settings(&suggest(&worker("timeout"), &d, false)),
            ["element_size"]
        );
        assert!(suggest(&worker("died"), &d, false).is_empty());
        // Deviation: finer, then accepting what was measured.
        let dev = error(DefectCode::VolBoundaryDeviation).measured("hausdorff_relative", 0.08);
        let rs = suggest(&dev, &d, false);
        assert_eq!(settings(&rs), ["element_size", "max_deviation"]);
        assert_eq!(rs[0].scale, Some(0.5));
        assert_eq!(rs[1].value, Some(Value::from(0.08 * DEVIATION_MARGIN)));
        // A wrap too far from the surface.
        let mut wrap = error(DefectCode::TopBoundaryEdge).measured("deviation", 0.07);
        wrap.action = Some(Action::Wrap);
        wrap.outcome = Some(Outcome::Failed);
        assert_eq!(
            settings(&suggest(&wrap, &d, false)),
            ["wrap_alpha", "max_deviation"]
        );
        // A selection naming a missing part.
        let missing = Record::new(DefectCode::PartEmptySelection, Stage::Select);
        let picked = JobConfig {
            parts: vec!["nope".into()],
            ..d.clone()
        };
        assert_eq!(settings(&suggest(&missing, &picked, false)), ["parts"]);
        assert!(suggest(&missing, &d, false).is_empty());
        // A part touching itself: no other kernel, no remeshing.
        assert!(suggest(&error(DefectCode::MeshRecoveryFailed), &no_remesh, true).is_empty());
        assert_eq!(
            settings(&suggest(&error(DefectCode::MeshStrictUnmet), &strict, true)),
            ["strict"]
        );
        let mut split = Record::new(DefectCode::TopNonmanifoldEdge, Stage::Repair).part("a");
        split.action = Some(Action::SplitEdge);
        split.outcome = Some(Outcome::Applied);
        let mut records = vec![
            split,
            error(DefectCode::MeshRecoveryFailed).part("a"),
            error(DefectCode::MeshRecoveryFailed).part("b"),
        ];
        fill(&mut records, &d);
        assert!(records[1].remedies.is_empty());
        assert_eq!(settings(&records[2].remedies), ["kernel"]);
        // Input problems and valid output: none.
        for code in [
            DefectCode::ParseReadFailed,
            DefectCode::TopOpenSheet,
            DefectCode::PartOverlap,
            DefectCode::VolInvertedTet,
            DefectCode::OutWriteFailed,
        ] {
            assert!(suggest(&error(code), &d, false).is_empty(), "{code:?}");
        }
        assert!(
            suggest(
                &error(DefectCode::MeshRecoveryFailed).severity(Severity::Info),
                &d,
                false
            )
            .is_empty()
        );
    }

    #[test]
    fn the_next_config_acts_on_the_worst_record_first() {
        let d = JobConfig {
            parts: vec!["nope".into()],
            ..JobConfig::default()
        };
        let mut records = vec![
            error(DefectCode::MeshBudgetExceeded),
            Record::new(DefectCode::PartEmptySelection, Stage::Select),
        ];
        for r in &mut records {
            r.remedies = suggest(r, &d, false);
        }
        let n = next_config(&records, &d, Severity::Error).unwrap();
        assert_eq!((n.record, n.remedy), (1, 0));
        assert!(n.config.parts.is_empty());
        // Below the severity asked for: nothing.
        assert!(next_config(&records[..1], &d, Severity::Fatal).is_none());
        let n = next_config(&records[..1], &d, Severity::Error).unwrap();
        assert_eq!(n.record, 0);
        assert_eq!(n.config.element_size, d.element_size * 2.0);
        // A remedy that changes nothing, or does not apply, is passed over.
        let mut r = error(DefectCode::MeshRecoveryFailed);
        r.remedies = vec![
            Remedy::set(Setting::Kernel, "conformal", "no change"),
            Remedy::set(Setting::Kernel, "fastest", "not a kernel"),
            Remedy::scale(Setting::Kernel, 2.0, "not a number"),
            Remedy::scale(Setting::ElementSize, -1.0, "not positive"),
        ];
        assert!(next_config(std::slice::from_ref(&r), &d, Severity::Error).is_none());
        r.remedies
            .push(Remedy::scale(Setting::ElementSize, 0.5, "finer"));
        let n = next_config(&[r], &d, Severity::Error).unwrap();
        assert_eq!(n.remedy, 4);
    }
}
