//! What the face-based writers (OpenFOAM, Fluent) share: the names each
//! format accepts, and the points a face view uses.

use std::collections::BTreeSet;

use mesher_core::mesh::FaceMesh;
use mesher_core::{Action, DefectCode, Outcome, Record, RecordLog, Severity, Stage};

/// A format's rule for names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NameRule {
    /// OpenFOAM words, as written here: ASCII letters, digits and `_`, `.`,
    /// `-`, starting with a letter or `_` (a subset of what OpenFOAM
    /// accepts, free of its punctuation and of characters a shell or file
    /// system treats specially).
    Foam,
    /// Fluent zone names (Scheme symbols), as written here: lowercase ASCII
    /// letters, digits and `_`, `.`, `+`, `-`, starting with a letter or `_`.
    /// Capitals are lowered.
    Fluent,
}

impl NameRule {
    /// The format's name in `OUT_NAME_CHANGED` records.
    pub(crate) fn format(self) -> &'static str {
        match self {
            NameRule::Foam => "openfoam",
            NameRule::Fluent => "fluent",
        }
    }

    /// `name` under the rule: other characters become `_`, and a name that
    /// is empty or does not start with a letter or `_` gets a leading `_`.
    pub(crate) fn apply(self, name: &str) -> String {
        let keep = |c: char| match self {
            NameRule::Foam => c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'),
            NameRule::Fluent => {
                c.is_ascii_lowercase() || c.is_ascii_digit() || matches!(c, '_' | '.' | '+' | '-')
            }
        };
        let mut w: String = name
            .chars()
            .map(|c| match self {
                NameRule::Fluent => c.to_ascii_lowercase(),
                NameRule::Foam => c,
            })
            .map(|c| if keep(c) { c } else { '_' })
            .collect();
        if !w.starts_with(|c: char| c.is_ascii_alphabetic() || c == '_') {
            w.insert(0, '_');
        }
        w
    }
}

/// Names for `items` (each an entity kind and its name) under `rule`,
/// distinct, each recorded as `OUT_NAME_CHANGED` where it differs (a
/// warning, or information when only the letter case changed). A name
/// that already follows the rule and is not repeated keeps it (the first of
/// a repeated name keeps it); the others take the first of `w`, `w_2`,
/// `w_3`, ... that is free.
pub(crate) fn names(rule: NameRule, items: &[(&str, &str)], log: &mut RecordLog) -> Vec<String> {
    let mut used: BTreeSet<String> = BTreeSet::new();
    let mut out: Vec<Option<String>> = vec![None; items.len()];
    for (i, &(_, name)) in items.iter().enumerate() {
        if rule.apply(name) == name && used.insert(name.to_owned()) {
            out[i] = Some(name.to_owned());
        }
    }
    for (i, &(entity, name)) in items.iter().enumerate() {
        if out[i].is_some() {
            continue;
        }
        let base = rule.apply(name);
        let mut w = base.clone();
        let mut k = 2u32;
        while used.contains(&w) {
            w = format!("{base}_{k}");
            k += 1;
        }
        used.insert(w.clone());
        // Only the letter case changed (Fluent's lowercase names): nothing
        // is lost, so it is recorded as information.
        let severity = if w == name.to_ascii_lowercase() {
            Severity::Info
        } else {
            Severity::Warning
        };
        log.push(
            Record::new(DefectCode::OutNameChanged, Stage::Emit)
                .severity(severity)
                .measured("format", rule.format())
                .measured("entity", entity)
                .measured("name", name)
                .measured("written_as", w.clone())
                .action(Action::Rename, Outcome::Applied),
        );
        out[i] = Some(w);
    }
    out.into_iter().map(Option::unwrap_or_default).collect()
}

/// The points `f`'s faces use, in order: each mesh point's new number
/// (`u32::MAX` where unused) and the mesh points kept.
pub(crate) fn used_points(points: usize, f: &FaceMesh) -> (Vec<u32>, Vec<usize>) {
    let mut new = vec![u32::MAX; points];
    for v in f.faces.iter().flatten() {
        if let Some(n) = new.get_mut(*v as usize) {
            *n = 0;
        }
    }
    let mut kept = Vec::new();
    for (v, n) in new.iter_mut().enumerate() {
        if *n == 0 {
            *n = kept.len() as u32;
            kept.push(v);
        }
    }
    (new, kept)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    #[test]
    fn fluent_names_are_lowercase_symbols() {
        let mut log = RecordLog::new();
        let n = names(
            NameRule::Fluent,
            &[
                ("interior", "interior"),
                ("cell_zone", "Body"),
                ("patch", "body"),
                ("patch", "interior"),
                ("patch", "a+b.c-d"),
                ("patch", "9 lives"),
                ("patch", "Wall"),
            ],
            &mut log,
        );
        assert_eq!(
            n,
            [
                "interior",
                "body_2",
                "body",
                "interior_2",
                "a+b.c-d",
                "_9_lives",
                "wall"
            ]
        );
        assert_eq!(log.records().len(), 4);
        // Lowering `Wall` alone loses nothing: information; `Body`, which also
        // took a suffix, is a warning.
        let severity = |name: &str| {
            log.records()
                .iter()
                .find(|r| r.measured["name"].as_str() == Some(name))
                .unwrap()
                .severity
        };
        assert_eq!(severity("Wall"), Severity::Info);
        assert_eq!(severity("Body"), Severity::Warning);
        assert_eq!(severity("9 lives"), Severity::Warning);
        assert!(
            log.records()
                .iter()
                .all(|r| r.measured["format"].as_str() == Some("fluent"))
        );
    }
}
