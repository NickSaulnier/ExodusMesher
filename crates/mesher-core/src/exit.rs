//! Process exit codes (see `PLAN.md` > Error reporting design > Exit codes).

/// Outcome of a run as a process exit code. Scripts branch on these without
/// parsing JSON.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ExitStatus {
    /// No warnings, errors or fatal records.
    Clean,
    /// Succeeded, but at least one warning (a repair) was recorded.
    Repaired,
    /// Succeeded via the stuffing fallback.
    Fallback,
    /// At least one error or fatal record.
    Failed,
    /// A bug: the worker panicked or broke the protocol.
    InternalBug,
}

impl ExitStatus {
    #[must_use]
    pub fn code(self) -> i32 {
        match self {
            ExitStatus::Clean => 0,
            ExitStatus::Repaired => 10,
            ExitStatus::Fallback => 20,
            ExitStatus::Failed => 30,
            ExitStatus::InternalBug => 70,
        }
    }

    #[must_use]
    pub fn from_code(code: i32) -> Option<ExitStatus> {
        [
            ExitStatus::Clean,
            ExitStatus::Repaired,
            ExitStatus::Fallback,
            ExitStatus::Failed,
            ExitStatus::InternalBug,
        ]
        .into_iter()
        .find(|s| s.code() == code)
    }

    /// The exit-code rule: errors or fatals fail; otherwise a fallback attempt
    /// gives `Fallback`; otherwise any warning gives `Repaired`; otherwise `Clean`.
    #[must_use]
    pub fn from_counts(warnings: u64, errors: u64, fatal: u64, used_fallback: bool) -> ExitStatus {
        if errors > 0 || fatal > 0 {
            ExitStatus::Failed
        } else if used_fallback {
            ExitStatus::Fallback
        } else if warnings > 0 {
            ExitStatus::Repaired
        } else {
            ExitStatus::Clean
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rule_precedence() {
        assert_eq!(ExitStatus::from_counts(0, 0, 0, false), ExitStatus::Clean);
        assert_eq!(
            ExitStatus::from_counts(3, 0, 0, false),
            ExitStatus::Repaired
        );
        assert_eq!(ExitStatus::from_counts(3, 0, 0, true), ExitStatus::Fallback);
        assert_eq!(ExitStatus::from_counts(3, 1, 0, true), ExitStatus::Failed);
        assert_eq!(ExitStatus::from_counts(0, 0, 1, false), ExitStatus::Failed);
    }

    #[test]
    fn codes_round_trip() {
        for code in [0, 10, 20, 30, 70] {
            assert_eq!(
                ExitStatus::from_code(code).map(ExitStatus::code),
                Some(code)
            );
        }
        assert_eq!(ExitStatus::from_code(1), None);
    }
}
