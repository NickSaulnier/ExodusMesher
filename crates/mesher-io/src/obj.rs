//! Tolerant OBJ reader with line provenance.
//!
//! Accepts any byte sequence. Every line it cannot use becomes a record
//! (`PARSE_BAD_LINE`, `PARSE_BAD_INDEX`, `PARSE_UNSUPPORTED_ELEMENT`), never a
//! panic and never a silent skip. Normals, texture coordinates, materials and
//! smoothing groups are read past without records, since they carry no geometry.
//!
//! Index conventions: OBJ vertex numbers are 1-based and may be negative
//! (relative to the vertices defined so far). Positive indices may refer forward;
//! they are checked against the final vertex count. A malformed `v` line still
//! occupies its vertex number, so later indices keep their meaning; any face
//! using such a vertex is dropped with `PARSE_BAD_INDEX`.

use std::collections::HashMap;
use std::path::Path;

use mesher_core::{Action, DefectCode, InputStats, Outcome, Record, RecordLog, Severity, Stage};
use mesher_geom::polygon::{
    FanReason, MAX_EAR_CLIP_VERTICES, TriangulationMethod, triangulate_polygon,
};

/// A source face as written in the file, with vertex indices into
/// [`ObjModel::vertices`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polygon {
    pub vertices: Vec<u32>,
    /// 1-based OBJ line (the first line of a continued statement).
    pub line: u32,
    /// Index into [`ObjModel::objects`] of the active `o` statement.
    pub object: Option<u32>,
    /// Index into [`ObjModel::groups`] of the active `g` statement.
    pub group: Option<u32>,
}

/// A triangle from triangulating a polygon; winding follows the polygon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Triangle {
    pub vertices: [u32; 3],
    /// Index into [`ObjModel::polygons`].
    pub polygon: u32,
}

/// The parsed geometry plus provenance.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ObjModel {
    /// Well-formed vertices, in file order.
    pub vertices: Vec<[f64; 3]>,
    /// 1-based OBJ line of each vertex.
    pub vertex_lines: Vec<u32>,
    pub polygons: Vec<Polygon>,
    pub triangles: Vec<Triangle>,
    /// `o` names in order of first appearance.
    pub objects: Vec<String>,
    /// `g` names in order of first appearance.
    pub groups: Vec<String>,
    pub bytes: u64,
    pub lines: u64,
}

impl ObjModel {
    #[must_use]
    pub fn stats(&self) -> InputStats {
        InputStats {
            bytes: self.bytes,
            lines: self.lines,
            vertices: self.vertices.len() as u64,
            polygons: self.polygons.len() as u64,
            triangles: self.triangles.len() as u64,
            objects: self.objects.clone(),
            groups: self.groups.clone(),
        }
    }
}

/// Read and parse an OBJ file. On an I/O error, records `PARSE_READ_FAILED` and
/// returns `None`.
pub fn read_obj(
    path: &Path,
    log: &mut RecordLog,
    progress: &mut dyn FnMut(f64),
) -> Option<ObjModel> {
    match std::fs::read(path) {
        Ok(bytes) => Some(parse_obj(&bytes, log, progress)),
        Err(err) => {
            log.push(
                Record::new(DefectCode::ParseReadFailed, Stage::Parse)
                    .measured("io_error", err.kind().to_string())
                    .action(Action::Abort, Outcome::Applied)
                    .hint(format!("could not read {}: {err}", path.display())),
            );
            None
        }
    }
}

/// Statements that describe lines, points, free-form curves and surfaces, or
/// that would run external commands. None of them are meshed.
const UNSUPPORTED: &[&[u8]] = &[
    b"l", b"p", b"vp", b"curv", b"curv2", b"surf", b"cstype", b"deg", b"bmat", b"step", b"parm",
    b"trim", b"hole", b"scrv", b"sp", b"end", b"con", b"call", b"csh",
];

/// Statements with no geometric content, read past without a record.
const IGNORED: &[&[u8]] = &[
    b"vt",
    b"vn",
    b"s",
    b"usemtl",
    b"mtllib",
    b"mg",
    b"lod",
    b"bevel",
    b"c_interp",
    b"d_interp",
    b"shadow_obj",
    b"trace_obj",
    b"maplib",
    b"usemap",
];

/// Progress is reported every this many physical lines.
const PROGRESS_INTERVAL: u64 = 1 << 16;

/// Parse OBJ bytes. Never panics; every problem becomes a record in `log`.
/// `progress` receives the fraction of bytes consumed.
pub fn parse_obj(bytes: &[u8], log: &mut RecordLog, progress: &mut dyn FnMut(f64)) -> ObjModel {
    let mut parser = Parser::default();
    let mut pending: Vec<u8> = Vec::new();
    let mut pending_start: Option<u64> = None;
    let mut lines = 0;

    for (number, raw, end) in PhysicalLines::new(bytes) {
        lines = number;
        if number % PROGRESS_INTERVAL == 0 {
            progress(end as f64 / bytes.len().max(1) as f64);
        }
        let content = strip_comment(raw);
        if let Some(body) = content.trim_ascii_end().strip_suffix(b"\\") {
            pending_start.get_or_insert(number);
            pending.extend_from_slice(body);
            pending.push(b' ');
            continue;
        }
        match pending_start.take() {
            Some(start) => {
                pending.extend_from_slice(content);
                parser.statement(start, &pending, log);
                pending.clear();
            }
            None => parser.statement(number, content, log),
        }
    }
    if let Some(start) = pending_start {
        parser.statement(start, &pending, log);
    }
    progress(1.0);

    let mut model = parser.finish(log);
    model.bytes = bytes.len() as u64;
    model.lines = lines;
    model
}

/// Physical lines split on LF, CRLF or lone CR; yields (1-based number, line,
/// end offset).
struct PhysicalLines<'a> {
    bytes: &'a [u8],
    pos: usize,
    number: u64,
}

impl<'a> PhysicalLines<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        PhysicalLines {
            bytes,
            pos: 0,
            number: 0,
        }
    }
}

impl<'a> Iterator for PhysicalLines<'a> {
    type Item = (u64, &'a [u8], usize);

    fn next(&mut self) -> Option<Self::Item> {
        let rest = self.bytes.get(self.pos..)?;
        if rest.is_empty() {
            return None;
        }
        self.number += 1;
        let (len, terminator) = match rest.iter().position(|&b| b == b'\n' || b == b'\r') {
            Some(i) if rest[i] == b'\r' && rest.get(i + 1) == Some(&b'\n') => (i, 2),
            Some(i) => (i, 1),
            None => (rest.len(), 0),
        };
        self.pos += len + terminator;
        Some((self.number, &rest[..len], self.pos))
    }
}

fn strip_comment(line: &[u8]) -> &[u8] {
    match line.iter().position(|&b| b == b'#') {
        Some(i) => &line[..i],
        None => line,
    }
}

fn line_u32(line: u64) -> u32 {
    u32::try_from(line).unwrap_or(u32::MAX)
}

/// A face before index validation. Vertex numbers are 0-based into the raw
/// (pre-compaction) vertex list; positive OBJ indices are unchecked until the end.
struct RawFace {
    line: u32,
    indices: Vec<RawIndex>,
    object: Option<u32>,
    group: Option<u32>,
}

#[derive(Clone, Copy)]
enum RawIndex {
    /// Resolved 0-based vertex number.
    Resolved(u64),
    /// Positive OBJ index (1-based), checked against the final count.
    Forward(u64),
}

#[derive(Default)]
struct Parser {
    /// Every `v` statement, `None` where the line was malformed.
    raw_vertices: Vec<Option<[f64; 3]>>,
    raw_vertex_lines: Vec<u32>,
    faces: Vec<RawFace>,
    objects: Vec<String>,
    object_ids: HashMap<String, u32>,
    groups: Vec<String>,
    group_ids: HashMap<String, u32>,
    object: Option<u32>,
    group: Option<u32>,
}

fn bad_line(line: u32, hint: impl Into<String>) -> Record {
    Record::new(DefectCode::ParseBadLine, Stage::Parse)
        .line(line)
        .action(Action::SkipLine, Outcome::Applied)
        .hint(hint)
}

fn intern(name: String, names: &mut Vec<String>, ids: &mut HashMap<String, u32>) -> Option<u32> {
    if let Some(&id) = ids.get(&name) {
        return Some(id);
    }
    let id = u32::try_from(names.len()).ok()?;
    ids.insert(name.clone(), id);
    names.push(name);
    Some(id)
}

impl Parser {
    fn statement(&mut self, line: u64, content: &[u8], log: &mut RecordLog) {
        let line = line_u32(line);
        let mut tokens = content
            .split(|b| b.is_ascii_whitespace())
            .filter(|t| !t.is_empty());
        let Some(keyword) = tokens.next() else {
            return;
        };
        match keyword {
            b"v" => self.vertex(line, tokens, log),
            b"f" | b"fo" => self.face(line, tokens, log),
            b"o" => {
                self.object =
                    name(tokens).and_then(|n| intern(n, &mut self.objects, &mut self.object_ids));
            }
            b"g" => {
                self.group =
                    name(tokens).and_then(|n| intern(n, &mut self.groups, &mut self.group_ids));
            }
            k if IGNORED.contains(&k) => {}
            k if UNSUPPORTED.contains(&k) => {
                log.push(
                    Record::new(DefectCode::ParseUnsupportedElement, Stage::Parse)
                        .line(line)
                        .measured("statement", String::from_utf8_lossy(k).into_owned())
                        .action(Action::SkipLine, Outcome::Applied)
                        .hint("only polygonal faces are meshed"),
                );
            }
            k => {
                let shown: String = String::from_utf8_lossy(k).chars().take(32).collect();
                log.push(bad_line(line, format!("unknown statement `{shown}`")));
            }
        }
    }

    fn vertex<'a>(
        &mut self,
        line: u32,
        tokens: impl Iterator<Item = &'a [u8]>,
        log: &mut RecordLog,
    ) {
        let mut values = Vec::with_capacity(3);
        let mut problem = None;
        for t in tokens {
            match std::str::from_utf8(t)
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
            {
                Some(v) if v.is_finite() => values.push(v),
                Some(_) => problem = Some("non-finite coordinate"),
                None => problem = Some("non-numeric coordinate"),
            }
            if problem.is_some() || values.len() > 7 {
                break;
            }
        }
        let problem = problem.or(match values.len() {
            0..=2 => Some("vertex needs 3 coordinates"),
            3..=7 => None,
            _ => Some("too many values on vertex line"),
        });
        self.raw_vertex_lines.push(line);
        match problem {
            None => self
                .raw_vertices
                .push(Some([values[0], values[1], values[2]])),
            Some(hint) => {
                self.raw_vertices.push(None);
                log.push(bad_line(line, hint));
            }
        }
    }

    fn face<'a>(&mut self, line: u32, tokens: impl Iterator<Item = &'a [u8]>, log: &mut RecordLog) {
        let defined = self.raw_vertices.len() as u64;
        let mut indices = Vec::new();
        for token in tokens {
            let mut parts = token.split(|&b| b == b'/');
            let vertex = parts.next().and_then(parse_int);
            // Texture and normal indices are ignored but must be well-formed.
            let rest_ok = parts.all(|p| p.is_empty() || parse_int(p).is_some());
            let Some(vertex) = vertex.filter(|_| rest_ok) else {
                log.push(bad_line(line, "malformed face index"));
                return;
            };
            let resolved = match vertex {
                0 => None,
                v if v > 0 => Some(RawIndex::Forward(v.unsigned_abs())),
                v => defined
                    .checked_sub(v.unsigned_abs())
                    .map(RawIndex::Resolved),
            };
            let Some(resolved) = resolved else {
                log.push(
                    Record::new(DefectCode::ParseBadIndex, Stage::Parse)
                        .line(line)
                        .measured("index", vertex)
                        .measured("vertices_defined", defined)
                        .action(Action::DropFace, Outcome::Applied)
                        .hint("face index is zero or reaches before the first vertex"),
                );
                return;
            };
            indices.push(resolved);
        }
        if indices.len() < 3 {
            log.push(bad_line(line, "face needs at least 3 vertices"));
            return;
        }
        self.faces.push(RawFace {
            line,
            indices,
            object: self.object,
            group: self.group,
        });
    }

    fn finish(self, log: &mut RecordLog) -> ObjModel {
        let mut model = ObjModel {
            objects: self.objects,
            groups: self.groups,
            ..ObjModel::default()
        };

        // Compact out malformed vertices, keeping a map from raw number to index.
        let mut remap: Vec<Option<u32>> = Vec::with_capacity(self.raw_vertices.len());
        for (v, &line) in self.raw_vertices.iter().zip(&self.raw_vertex_lines) {
            let id = v.and_then(|_| u32::try_from(model.vertices.len()).ok());
            if let (Some(p), Some(_)) = (v, id) {
                model.vertices.push(*p);
                model.vertex_lines.push(line);
            }
            remap.push(id);
        }

        let total = self.raw_vertices.len() as u64;
        'faces: for face in self.faces {
            let mut vertices = Vec::with_capacity(face.indices.len());
            for index in &face.indices {
                let raw = match *index {
                    RawIndex::Resolved(r) => Some(r),
                    RawIndex::Forward(i) => Some(i - 1).filter(|&r| r < total),
                };
                let mapped = raw
                    .and_then(|r| usize::try_from(r).ok())
                    .and_then(|r| remap.get(r).copied().flatten());
                match mapped {
                    Some(id) => vertices.push(id),
                    None => {
                        let hint = if raw.is_some() {
                            "face uses a vertex whose line was malformed"
                        } else {
                            "face index is past the last vertex"
                        };
                        let shown = match *index {
                            RawIndex::Forward(i) => i.to_string(),
                            RawIndex::Resolved(r) => format!("vertex #{}", r + 1),
                        };
                        log.push(
                            Record::new(DefectCode::ParseBadIndex, Stage::Parse)
                                .line(face.line)
                                .measured("index", shown)
                                .measured("vertices_defined", total)
                                .action(Action::DropFace, Outcome::Applied)
                                .hint(hint),
                        );
                        continue 'faces;
                    }
                }
            }
            let Ok(polygon_id) = u32::try_from(model.polygons.len()) else {
                log.push(bad_line(face.line, "too many faces"));
                continue;
            };
            triangulate(&mut model, polygon_id, &vertices, face.line, log);
            model.polygons.push(Polygon {
                vertices,
                line: face.line,
                object: face.object,
                group: face.group,
            });
        }

        if model.triangles.is_empty() {
            log.push(
                Record::new(DefectCode::ParseNoGeometry, Stage::Parse)
                    .severity(Severity::Fatal)
                    .measured("vertices", model.vertices.len() as u64)
                    .measured("polygons", model.polygons.len() as u64)
                    .action(Action::Abort, Outcome::Applied)
                    .hint("the file contains no usable faces"),
            );
        }
        model
    }
}

fn triangulate(
    model: &mut ObjModel,
    polygon: u32,
    vertices: &[u32],
    line: u32,
    log: &mut RecordLog,
) {
    let points: Vec<[f64; 3]> = vertices
        .iter()
        .map(|&v| model.vertices[v as usize])
        .collect();
    let result = triangulate_polygon(&points);
    if let TriangulationMethod::Fan(reason) = result.method {
        let mut record = Record::new(DefectCode::ParsePolygonFallback, Stage::Parse)
            .line(line)
            .face(polygon)
            .measured("polygon_vertices", vertices.len() as u64)
            .measured("reason", reason.as_str())
            .action(Action::FanTriangulate, Outcome::Applied);
        if reason == FanReason::TooManyVertices {
            record = record.threshold("max_ear_clip_vertices", MAX_EAR_CLIP_VERTICES as u64);
        }
        log.push(record);
    }
    model
        .triangles
        .extend(result.triangles.iter().map(|t| Triangle {
            vertices: t.map(|i| vertices[i]),
            polygon,
        }));
}

fn parse_int(token: &[u8]) -> Option<i64> {
    std::str::from_utf8(token).ok()?.parse().ok()
}

/// The statement's name: the remaining tokens joined by single spaces, with
/// invalid UTF-8 replaced. `None` for an unnamed `o`/`g`.
fn name<'a>(tokens: impl Iterator<Item = &'a [u8]>) -> Option<String> {
    let joined = tokens
        .map(|t| String::from_utf8_lossy(t).into_owned())
        .collect::<Vec<_>>()
        .join(" ");
    (!joined.is_empty()).then_some(joined)
}

#[cfg(test)]
mod tests {
    use super::*;
    use mesher_core::{DefectCode, Severity};

    fn parse(text: &str) -> (ObjModel, RecordLog) {
        let mut log = RecordLog::new();
        let model = parse_obj(text.as_bytes(), &mut log, &mut |_| {});
        (model, log)
    }

    fn codes(log: &RecordLog) -> Vec<(DefectCode, Vec<u32>)> {
        log.records()
            .iter()
            .map(|r| (r.code, r.location.obj_lines.clone()))
            .collect()
    }

    const CUBE: &str = "\
# unit cube, quads
o cube
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
v 0 0 1
v 1 0 1
v 1 1 1
v 0 1 1
vn 0 0 1
vt 0 0
s off
usemtl steel
f 1 4 3 2
f 5 6 7 8
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
";

    #[test]
    fn clean_cube() {
        let (m, log) = parse(CUBE);
        assert!(log.records().is_empty(), "{:?}", codes(&log));
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.polygons.len(), 6);
        assert_eq!(m.triangles.len(), 12);
        assert_eq!(m.objects, vec!["cube"]);
        assert_eq!(m.vertex_lines[0], 3);
        assert_eq!(m.polygons[0].line, 15);
        assert!(m.polygons.iter().all(|p| p.object == Some(0)));
    }

    #[test]
    fn face_tokens_with_texture_and_normal_indices() {
        let (m, log) = parse("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1/1/1 2//2 3/3\n");
        assert!(log.records().is_empty());
        assert_eq!(m.triangles[0].vertices, [0, 1, 2]);
    }

    #[test]
    fn negative_and_forward_indices() {
        let (m, log) = parse("f 1 2 3\nv 0 0 0\nv 1 0 0\nv 0 1 0\nf -3 -2 -1\n");
        assert!(log.records().is_empty(), "{:?}", codes(&log));
        assert_eq!(m.triangles.len(), 2);
        assert_eq!(m.triangles[1].vertices, [0, 1, 2]);
    }

    #[test]
    fn bad_indices_drop_the_face() {
        let (m, log) = parse("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 9\nf 0 1 2\nf -4 1 2\nf 1 2 3\n");
        assert_eq!(m.polygons.len(), 1);
        assert_eq!(
            codes(&log),
            vec![
                (DefectCode::ParseBadIndex, vec![5]),
                (DefectCode::ParseBadIndex, vec![6]),
                (DefectCode::ParseBadIndex, vec![4]),
            ]
        );
    }

    #[test]
    fn malformed_vertex_keeps_numbering() {
        // Vertex 2 is malformed; vertex 4 must still be the fourth `v` line.
        let (m, log) = parse("v 0 0 0\nv 1 nope 0\nv 1 0 0\nv 0 1 0\nf 1 3 4\nf 1 2 3\n");
        assert_eq!(m.vertices.len(), 3);
        assert_eq!(m.triangles.len(), 1);
        assert_eq!(m.triangles[0].vertices, [0, 1, 2]);
        assert_eq!(
            codes(&log),
            vec![
                (DefectCode::ParseBadLine, vec![2]),
                (DefectCode::ParseBadIndex, vec![6])
            ]
        );
    }

    #[test]
    fn bad_lines_are_reported_not_fatal() {
        let text = "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 1 2\nv nan 0 0\nv 1e999 0 0\nf 1 2\nf 1 x 3\nwat 1 2\nf 1 2 3\n";
        let (m, log) = parse(text);
        assert_eq!(m.triangles.len(), 1);
        let lines: Vec<u32> = log
            .records()
            .iter()
            .flat_map(|r| r.location.obj_lines.clone())
            .collect();
        assert_eq!(lines, vec![4, 5, 6, 7, 8, 9]);
        assert!(
            log.records()
                .iter()
                .all(|r| r.code == DefectCode::ParseBadLine)
        );
        assert!(!log.has_fatal());
    }

    #[test]
    fn unsupported_elements_are_info() {
        let (_, log) = parse("v 0 0 0\nv 1 0 0\nv 0 1 0\nl 1 2\ncurv 0 1 1 2\nf 1 2 3\n");
        assert_eq!(log.records().len(), 2);
        assert!(
            log.records()
                .iter()
                .all(|r| r.code == DefectCode::ParseUnsupportedElement
                    && r.severity == Severity::Info)
        );
    }

    #[test]
    fn line_endings_continuations_and_comments() {
        let text = "v 0 0 0\r\nv 1 0 0\rv 0 1 0 # trailing comment\nf 1 \\\n 2 \\\n 3\n# done";
        let (m, log) = parse(text);
        assert!(log.records().is_empty(), "{:?}", codes(&log));
        assert_eq!(m.vertex_lines, vec![1, 2, 3]);
        assert_eq!(m.polygons[0].line, 4);
        assert_eq!(m.lines, 7);
    }

    #[test]
    fn groups_and_objects_are_tracked() {
        let text =
            "v 0 0 0\nv 1 0 0\nv 0 1 0\no A\ng left side\nf 1 2 3\ng\nf 1 2 3\no B\nf 1 2 3\n";
        let (m, _) = parse(text);
        assert_eq!(m.objects, vec!["A", "B"]);
        assert_eq!(m.groups, vec!["left side"]);
        let tags: Vec<_> = m.polygons.iter().map(|p| (p.object, p.group)).collect();
        assert_eq!(
            tags,
            vec![(Some(0), Some(0)), (Some(0), None), (Some(1), None)]
        );
    }

    #[test]
    fn empty_and_garbage_input_is_fatal_no_geometry() {
        for input in [
            &b""[..],
            b"\x00\xff\xfe\x01",
            b"# nothing\n\n",
            b"v 1 2 3\n",
        ] {
            let mut log = RecordLog::new();
            let m = parse_obj(input, &mut log, &mut |_| {});
            assert!(m.triangles.is_empty());
            assert!(log.has_fatal());
            assert!(
                log.records()
                    .iter()
                    .any(|r| r.code == DefectCode::ParseNoGeometry)
            );
        }
    }

    #[test]
    fn invalid_utf8_names_survive() {
        let mut log = RecordLog::new();
        let m = parse_obj(
            b"o caf\xe9\nv 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n",
            &mut log,
            &mut |_| {},
        );
        assert!(log.records().is_empty());
        assert_eq!(m.objects.len(), 1);
    }

    #[test]
    fn degenerate_polygon_records_fallback() {
        let (m, log) = parse("v 0 0 0\nv 1 0 0\nv 2 0 0\nv 3 0 0\nf 1 2 3 4\n");
        assert_eq!(m.triangles.len(), 2);
        assert_eq!(
            codes(&log),
            vec![(DefectCode::ParsePolygonFallback, vec![5])]
        );
        assert_eq!(log.records()[0].location.faces, vec![0]);
    }
}
