//! Report types, exit codes and the worker protocol shared by every ExodusMesher
//! crate. Later milestones add parts, diagnosis, repair, gates and sizing here.

pub mod diagnose;
pub mod exit;
pub mod mesh;
pub mod model;
pub mod overlay;
pub mod project;
pub mod protocol;
pub mod repair;
pub mod report;
pub mod volume;

pub use exit::ExitStatus;
pub use report::{
    Action, Attempt, AttemptOutcome, DefectCode, InputInfo, InputStats, JobConfig, Kernel,
    KernelChoice, Location, MeshSummary, Outcome, PartSource, PartSummary, QualityMeasure, Record,
    RecordLog, RepairSummary, Report, Severity, Stage, Summary, WorstCell,
};
