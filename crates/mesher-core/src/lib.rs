//! Report types, exit codes and the worker protocol shared by every ExodusMesher
//! crate. Later milestones add parts, diagnosis, repair, gates and sizing here.

pub mod exit;
pub mod protocol;
pub mod report;

pub use exit::ExitStatus;
pub use report::{
    Action, Attempt, AttemptOutcome, DefectCode, InputInfo, InputStats, JobConfig, Kernel,
    Location, Outcome, Record, RecordLog, Report, Severity, Stage, Summary,
};
