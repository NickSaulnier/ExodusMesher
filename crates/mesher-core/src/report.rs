//! The structured report: defect codes, records, attempts and the run summary.
//!
//! Every finding is a [`Record`] with a stable [`DefectCode`]. The code list here
//! and the taxonomy table in `PLAN.md` must stay identical; a test enforces it.

use std::collections::BTreeMap;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::exit::ExitStatus;

/// Version of the `report.json` layout. Bump on any breaking change.
pub const SCHEMA_VERSION: u32 = 1;

/// At most this many detailed records are kept per defect code; the summary
/// still counts every occurrence.
pub const MAX_RECORDS_PER_CODE: u64 = 100;

/// At most this many entries per location list.
pub const MAX_LOCATION_ITEMS: usize = 256;

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    /// Noteworthy, no action needed.
    Info,
    /// Repaired automatically; output still valid.
    Warning,
    /// Blocks a stage; a fallback may still succeed.
    Error,
    /// No usable output possible.
    Fatal,
}

macro_rules! defect_codes {
    ($( $variant:ident = $name:literal, $severity:ident, $doc:literal; )*) => {
        /// Stable defect codes. Serialized as `SCREAMING_SNAKE_CASE` strings.
        #[derive(
            Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash,
            Serialize, Deserialize, JsonSchema,
        )]
        pub enum DefectCode {
            $(
                #[doc = $doc]
                #[serde(rename = $name)]
                $variant,
            )*
        }

        impl DefectCode {
            /// Every code, in taxonomy order.
            pub const ALL: &'static [DefectCode] = &[$(DefectCode::$variant,)*];

            /// The code as it appears in reports and `PLAN.md`.
            #[must_use]
            pub fn as_str(self) -> &'static str {
                match self { $(DefectCode::$variant => $name,)* }
            }

            /// Severity when the default action applies.
            #[must_use]
            pub fn default_severity(self) -> Severity {
                match self { $(DefectCode::$variant => Severity::$severity,)* }
            }

            /// One-line meaning, as in the taxonomy table.
            #[must_use]
            pub fn description(self) -> &'static str {
                match self { $(DefectCode::$variant => $doc,)* }
            }
        }
    };
}

defect_codes! {
    ParseBadLine = "PARSE_BAD_LINE", Warning, "Malformed or non-numeric OBJ line";
    ParseBadIndex = "PARSE_BAD_INDEX", Warning, "Face references a missing vertex";
    ParseReadFailed = "PARSE_READ_FAILED", Fatal, "Input file cannot be opened or read";
    ParseNoGeometry = "PARSE_NO_GEOMETRY", Fatal, "No faces survive parsing";
    ParseUnsupportedElement = "PARSE_UNSUPPORTED_ELEMENT", Info, "Lines, points, curves or surfaces in the OBJ file";
    ParsePolygonFallback = "PARSE_POLYGON_FALLBACK", Warning, "Polygon could not be ear-clipped and was fan-triangulated";
    PartEmptySelection = "PART_EMPTY_SELECTION", Fatal, "No parts selected, or a named part does not exist";
    PartOverlap = "PART_OVERLAP", Error, "Two selected parts intersect";
    PartSharedInterface = "PART_SHARED_INTERFACE", Info, "Two selected parts touch along a surface";
    GeoDuplicateVertex = "GEO_DUPLICATE_VERTEX", Warning, "Vertices closer than tolerance";
    GeoDegenerateFace = "GEO_DEGENERATE_FACE", Warning, "Zero-area or collinear face";
    GeoDuplicateFace = "GEO_DUPLICATE_FACE", Warning, "Same face listed twice";
    TopBoundaryEdge = "TOP_BOUNDARY_EDGE", Warning, "Hole: edge with one adjacent face";
    TopNonmanifoldEdge = "TOP_NONMANIFOLD_EDGE", Warning, "Edge shared by three or more faces";
    TopNonmanifoldVertex = "TOP_NONMANIFOLD_VERTEX", Warning, "Vertex joining separate fans";
    TopInconsistentOrientation = "TOP_INCONSISTENT_ORIENTATION", Warning, "Neighbouring faces disagree on winding";
    TopInvertedShell = "TOP_INVERTED_SHELL", Warning, "Shell with negative signed volume";
    TopOpenSheet = "TOP_OPEN_SHEET", Fatal, "Surface encloses no volume";
    IntSelfIntersection = "INT_SELF_INTERSECTION", Error, "Faces of one shell intersect";
    IntShellOverlap = "INT_SHELL_OVERLAP", Error, "Two shells of one part intersect";
    QualSliverFace = "QUAL_SLIVER_FACE", Warning, "Minimum angle below threshold";
    QualSuspectScale = "QUAL_SUSPECT_SCALE", Info, "Bounding box suggests wrong units";
    MeshRecoveryFailed = "MESH_RECOVERY_FAILED", Error, "Constrained Delaunay could not recover a facet";
    MeshBudgetExceeded = "MESH_BUDGET_EXCEEDED", Error, "Steiner-point or element budget hit";
    GeoCapacityExceeded = "GEO_CAPACITY_EXCEEDED", Fatal, "Part too large for 32-bit vertex and face indices";
    RepairGuardExceeded = "REPAIR_GUARD_EXCEEDED", Error, "A repair changed the shell count or enclosed volume beyond its limit";
    MeshWorkerFailed = "MESH_WORKER_FAILED", Error, "Worker panicked, timed out or ran out of memory";
    VolInvertedTet = "VOL_INVERTED_TET", Error, "Tetrahedron with non-positive volume";
    VolLowQuality = "VOL_LOW_QUALITY", Warning, "Non-orthogonality, skewness or volume ratio out of bounds";
    VolBoundaryDeviation = "VOL_BOUNDARY_DEVIATION", Error, "Hausdorff distance to input over tolerance";
    ExoNameTruncated = "EXO_NAME_TRUNCATED", Warning, "Part or set name over the Exodus limit";
    ExoCapacityExceeded = "EXO_CAPACITY_EXCEEDED", Error, "Mesh too large for the 64-bit-offset format";
    OutWriteFailed = "OUT_WRITE_FAILED", Fatal, "An output file could not be written";
}

/// Pipeline stage that produced a record.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Stage {
    Parse,
    Diagnose,
    Select,
    Repair,
    SurfaceGate,
    Sizing,
    VolumeMesh,
    VolumeGate,
    Emit,
    /// The worker process itself (crash, timeout, protocol failure).
    Worker,
}

/// What was done about a finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Action {
    SkipLine,
    DropFace,
    FanTriangulate,
    Remesh,
    Merge,
    Remove,
    FillHole,
    SplitEdge,
    DuplicateVertex,
    Reorient,
    Flip,
    Truncate,
    Fallback,
    Abort,
}

/// Result of the action.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Outcome {
    /// The action was carried out.
    Applied,
    /// The action was attempted and did not succeed.
    Failed,
    /// Left for a later stage to decide.
    Deferred,
}

/// Where a finding is. `vertices` are 0-based indices into the parsed vertex
/// list and `faces` are 0-based indices into the triangulated face list (a
/// polygon with n corners becomes n - 2 consecutive triangles); `obj_lines` are
/// 1-based line numbers in the input file. Each list holds at most
/// [`MAX_LOCATION_ITEMS`] entries; a truncated list has its full length in the
/// record's `measured` map (`vertices_total`, `faces_total`, `obj_lines_total`).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Location {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub vertices: Vec<u32>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub faces: Vec<u32>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub obj_lines: Vec<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub point: Option<[f64; 3]>,
}

impl Location {
    fn is_empty(&self) -> bool {
        self.vertices.is_empty()
            && self.faces.is_empty()
            && self.obj_lines.is_empty()
            && self.point.is_none()
    }
}

/// One finding.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Record {
    pub code: DefectCode,
    pub severity: Severity,
    pub stage: Stage,
    /// Part the finding belongs to, once parts exist.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub part: Option<String>,
    #[serde(default, skip_serializing_if = "Location::is_empty")]
    pub location: Location,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub measured: BTreeMap<String, Value>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub threshold: BTreeMap<String, Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub action: Option<Action>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub outcome: Option<Outcome>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hint: Option<String>,
}

impl Record {
    /// A record with the code's default severity and nothing else filled in.
    #[must_use]
    pub fn new(code: DefectCode, stage: Stage) -> Record {
        Record {
            code,
            severity: code.default_severity(),
            stage,
            part: None,
            location: Location::default(),
            measured: BTreeMap::new(),
            threshold: BTreeMap::new(),
            action: None,
            outcome: None,
            hint: None,
        }
    }

    #[must_use]
    pub fn severity(mut self, severity: Severity) -> Record {
        self.severity = severity;
        self
    }

    #[must_use]
    pub fn line(mut self, obj_line: u32) -> Record {
        self.location.obj_lines.push(obj_line);
        self
    }

    #[must_use]
    pub fn vertices(mut self, vertices: impl IntoIterator<Item = u32>) -> Record {
        self.location.vertices.extend(vertices);
        self
    }

    #[must_use]
    pub fn face(mut self, face: u32) -> Record {
        self.location.faces.push(face);
        self
    }

    /// Set the face list (sorted, deduplicated, capped at [`MAX_LOCATION_ITEMS`]).
    #[must_use]
    pub fn faces(mut self, faces: impl IntoIterator<Item = u32>) -> Record {
        let (list, total) = capped(faces);
        if let Some(total) = total {
            self = self.measured("faces_total", total);
        }
        self.location.faces = list;
        self
    }

    /// Set the vertex list (sorted, deduplicated, capped).
    #[must_use]
    pub fn vertex_set(mut self, vertices: impl IntoIterator<Item = u32>) -> Record {
        let (list, total) = capped(vertices);
        if let Some(total) = total {
            self = self.measured("vertices_total", total);
        }
        self.location.vertices = list;
        self
    }

    /// Set the input-line list (sorted, deduplicated, capped; zeros dropped).
    #[must_use]
    pub fn lines(mut self, lines: impl IntoIterator<Item = u32>) -> Record {
        let (list, total) = capped(lines.into_iter().filter(|&l| l > 0));
        if let Some(total) = total {
            self = self.measured("obj_lines_total", total);
        }
        self.location.obj_lines = list;
        self
    }

    #[must_use]
    pub fn point(mut self, point: [f64; 3]) -> Record {
        self.location.point = Some(point);
        self
    }

    #[must_use]
    pub fn part(mut self, part: &str) -> Record {
        self.part = Some(part.to_owned());
        self
    }

    #[must_use]
    pub fn measured(mut self, key: &str, value: impl Into<Value>) -> Record {
        self.measured.insert(key.to_owned(), value.into());
        self
    }

    #[must_use]
    pub fn threshold(mut self, key: &str, value: impl Into<Value>) -> Record {
        self.threshold.insert(key.to_owned(), value.into());
        self
    }

    #[must_use]
    pub fn action(mut self, action: Action, outcome: Outcome) -> Record {
        self.action = Some(action);
        self.outcome = Some(outcome);
        self
    }

    #[must_use]
    pub fn hint(mut self, hint: impl Into<String>) -> Record {
        self.hint = Some(hint.into());
        self
    }
}

/// Sort, deduplicate and cap a list; the second value is the full length when
/// it was truncated.
fn capped(items: impl IntoIterator<Item = u32>) -> (Vec<u32>, Option<u64>) {
    let mut list: Vec<u32> = items.into_iter().collect();
    list.sort_unstable();
    list.dedup();
    let total = list.len();
    if total > MAX_LOCATION_ITEMS {
        list.truncate(MAX_LOCATION_ITEMS);
        (list, Some(total as u64))
    } else {
        (list, None)
    }
}

/// Collects records, keeping at most [`MAX_RECORDS_PER_CODE`] per code while
/// counting every one.
#[derive(Debug, Clone, Default)]
pub struct RecordLog {
    records: Vec<Record>,
    by_code: BTreeMap<DefectCode, CodeCount>,
    severity_counts: BTreeMap<Severity, u64>,
}

impl RecordLog {
    #[must_use]
    pub fn new() -> RecordLog {
        RecordLog::default()
    }

    /// Count the record and keep it if its code is under the cap. Returns the
    /// kept record, so callers can stream it.
    pub fn push(&mut self, record: Record) -> Option<&Record> {
        *self.severity_counts.entry(record.severity).or_default() += 1;
        let count = self.by_code.entry(record.code).or_default();
        count.total += 1;
        if count.reported < MAX_RECORDS_PER_CODE {
            count.reported += 1;
            self.records.push(record);
            self.records.last()
        } else {
            None
        }
    }

    #[must_use]
    pub fn records(&self) -> &[Record] {
        &self.records
    }

    #[must_use]
    pub fn count(&self, severity: Severity) -> u64 {
        self.severity_counts.get(&severity).copied().unwrap_or(0)
    }

    #[must_use]
    pub fn has_fatal(&self) -> bool {
        self.count(Severity::Fatal) > 0
    }
}

/// Per-code totals in the summary.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct CodeCount {
    /// Every occurrence found.
    pub total: u64,
    /// Occurrences with a detailed record in `records`.
    pub reported: u64,
}

/// Counts and the resulting exit code.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Summary {
    pub exit_code: i32,
    pub info: u64,
    pub warnings: u64,
    pub errors: u64,
    pub fatal: u64,
    pub by_code: BTreeMap<DefectCode, CodeCount>,
}

/// User-facing job settings (the policy flags). Recorded in the report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields, default)]
pub struct JobConfig {
    /// Names or indices of the parts to mesh; empty means all.
    pub parts: Vec<String>,
    /// Diagnose and report only.
    pub no_repair: bool,
    /// Conformal paths only; never use the stuffing fallback.
    pub strict: bool,
    /// Largest acceptable fallback deviation, relative to the bounding-box diagonal.
    pub max_deviation: f64,
}

impl Default for JobConfig {
    fn default() -> JobConfig {
        JobConfig {
            parts: Vec::new(),
            no_repair: false,
            strict: false,
            max_deviation: 0.01,
        }
    }
}

/// Volume kernel used by a meshing attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Kernel {
    Delaunay,
    Stuffing,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum AttemptOutcome {
    Succeeded,
    Failed,
}

/// One meshing attempt, in the order tried.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Attempt {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub part: Option<String>,
    pub kernel: Kernel,
    pub outcome: AttemptOutcome,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub stdout: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub stderr: String,
}

/// What the parser found in the input.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct InputStats {
    pub bytes: u64,
    pub lines: u64,
    pub vertices: u64,
    /// Faces as written in the file (polygons of any size).
    pub polygons: u64,
    /// Triangles after triangulating the polygons.
    pub triangles: u64,
    /// `o` names in order of first appearance.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub objects: Vec<String>,
    /// `g` names in order of first appearance.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub groups: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct InputInfo {
    /// The input path exactly as given.
    pub path: String,
    /// Absent when the file could not be read.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stats: Option<InputStats>,
}

/// How a part was formed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum PartSource {
    /// An `o` statement (faces before any `o` form the part `default`).
    Object,
    /// A `g` statement, in a file without `o` statements.
    Group,
    /// A connected component, in a file with neither; enclosed cavities join
    /// the component that contains them.
    Component,
}

/// One part, as the GUI part tree shows it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct PartSummary {
    /// 1-based; `--parts` accepts this or the name.
    pub index: u32,
    pub name: String,
    pub source: PartSource,
    pub triangles: u64,
    /// Distinct vertices after merging duplicates.
    pub vertices: u64,
    /// Edge-connected surface pieces.
    pub shells: u64,
    /// Bounding box `[min, max]`.
    pub bounds: [[f64; 3]; 2],
    pub selected: bool,
    /// Present when the part was repaired.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub repair: Option<RepairSummary>,
}

/// Face and vertex counts through repair, and the surface gate's verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct RepairSummary {
    /// Triangles in the part as parsed.
    pub faces_before: u64,
    pub faces_after: u64,
    /// Vertices entering repair (after merging duplicates).
    pub vertices_before: u64,
    pub vertices_after: u64,
    pub gate_passed: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct ToolInfo {
    pub name: String,
    pub version: String,
}

/// The complete `report.json`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Report {
    pub schema_version: u32,
    pub tool: ToolInfo,
    pub input: InputInfo,
    pub config: JobConfig,
    /// Stages that ran, in order.
    pub stages_run: Vec<Stage>,
    /// Parts found by diagnosis, in index order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub parts: Vec<PartSummary>,
    pub records: Vec<Record>,
    pub attempts: Vec<Attempt>,
    pub summary: Summary,
}

impl Report {
    /// An empty report for a job on `input_path`.
    #[must_use]
    pub fn new(input_path: &str, config: JobConfig) -> Report {
        Report {
            schema_version: SCHEMA_VERSION,
            tool: ToolInfo {
                name: "exodus-mesher".to_owned(),
                version: env!("CARGO_PKG_VERSION").to_owned(),
            },
            input: InputInfo {
                path: input_path.to_owned(),
                stats: None,
            },
            config,
            stages_run: Vec::new(),
            parts: Vec::new(),
            records: Vec::new(),
            attempts: Vec::new(),
            summary: Summary::default(),
        }
    }

    /// Move the log's records into the report and compute the summary and exit code.
    pub fn finish(&mut self, log: RecordLog) {
        let used_fallback = self
            .attempts
            .iter()
            .any(|a| a.kernel == Kernel::Stuffing && a.outcome == AttemptOutcome::Succeeded);
        let (info, warnings, errors, fatal) = (
            log.count(Severity::Info),
            log.count(Severity::Warning),
            log.count(Severity::Error),
            log.count(Severity::Fatal),
        );
        self.summary = Summary {
            exit_code: ExitStatus::from_counts(warnings, errors, fatal, used_fallback).code(),
            info,
            warnings,
            errors,
            fatal,
            by_code: log.by_code,
        };
        self.records = log.records;
    }

    /// Serialize as pretty JSON with a trailing newline. Deterministic: same
    /// report, same bytes.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = serde_json::to_string_pretty(self).unwrap_or_default();
        s.push('\n');
        s
    }

    /// The JSON schema of `report.json`, pretty-printed with a trailing newline.
    #[must_use]
    pub fn json_schema() -> String {
        let schema = schemars::schema_for!(Report);
        let mut s = serde_json::to_string_pretty(&schema).unwrap_or_default();
        s.push('\n');
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_serialize_as_their_names() {
        for &code in DefectCode::ALL {
            let json = serde_json::to_string(&code).unwrap();
            assert_eq!(json, format!("\"{}\"", code.as_str()));
            assert_eq!(serde_json::from_str::<DefectCode>(&json).unwrap(), code);
        }
    }

    #[test]
    fn log_caps_detail_but_counts_everything() {
        let mut log = RecordLog::new();
        for line in 0..250 {
            log.push(Record::new(DefectCode::ParseBadLine, Stage::Parse).line(line));
        }
        log.push(Record::new(DefectCode::ParseNoGeometry, Stage::Parse));
        let mut report = Report::new("x.obj", JobConfig::default());
        report.finish(log);
        assert_eq!(report.records.len() as u64, MAX_RECORDS_PER_CODE + 1);
        let bad = report.summary.by_code[&DefectCode::ParseBadLine];
        assert_eq!(bad.total, 250);
        assert_eq!(bad.reported, MAX_RECORDS_PER_CODE);
        assert_eq!(report.summary.warnings, 250);
        assert_eq!(report.summary.fatal, 1);
        assert_eq!(report.summary.exit_code, 30);
    }

    #[test]
    fn report_round_trips() {
        let mut log = RecordLog::new();
        log.push(
            Record::new(DefectCode::TopNonmanifoldEdge, Stage::Diagnose)
                .severity(Severity::Error)
                .vertices([1042, 1043])
                .face(2201)
                .line(1088)
                .measured("incident_faces", 3)
                .threshold("max_incident_faces", 2)
                .action(Action::SplitEdge, Outcome::Failed)
                .hint("two touching parts"),
        );
        let mut report = Report::new("in.obj", JobConfig::default());
        report.finish(log);
        let json = report.to_json();
        let back: Report = serde_json::from_str(&json).unwrap();
        assert_eq!(back, report);
        assert_eq!(back.to_json(), json);
    }
}
