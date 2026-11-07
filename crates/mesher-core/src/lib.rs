//! Report types, exit codes and the worker protocol shared by every ExodusMesher
//! crate. Later milestones add parts, diagnosis, repair, gates and sizing here.

pub mod diagnose;
pub mod exit;
pub mod mesh;
pub mod model;
pub mod overlay;
pub mod protocol;
pub mod repair;
pub mod report;
pub mod volume;

pub use exit::ExitStatus;
pub use report::{
    Action, Attempt, AttemptOutcome, DefectCode, InputInfo, InputStats, JobConfig, Kernel,
    Location, MeshSummary, Outcome, PartSource, PartSummary, Record, RecordLog, RepairSummary,
    Report, Severity, Stage, Summary,
};
