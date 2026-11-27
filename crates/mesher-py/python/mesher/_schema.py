"""The report's types, generated from ``schema/report.schema.json`` by
``crates/mesher-py/tests/report_module.rs``: do not edit; run
``UPDATE_REPORT_MODULE=1 cargo test -p mesher-py --test report_module``.

Each object in the schema is a frozen dataclass (built from the report's
JSON by its ``_parse``), each set of string constants a ``str`` enum. A field
the worker may leave out is its schema default when it has one, else empty
for a list or a map, else ``None``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Mapping, Optional, Tuple

__all__ = [
    "Action",
    "Attempt",
    "AttemptOutcome",
    "CodeCount",
    "DefectCode",
    "FeatureSummary",
    "InputInfo",
    "InputStats",
    "JobConfig",
    "Kernel",
    "KernelChoice",
    "Location",
    "MeshSummary",
    "Outcome",
    "OutputFormat",
    "PartSource",
    "PartSummary",
    "QualityMeasure",
    "Record",
    "Remedy",
    "RepairSummary",
    "Report",
    "Setting",
    "Severity",
    "SizeControl",
    "Stage",
    "Summary",
    "ToolInfo",
    "WorstCell",
]


class Action(str, Enum):
    """What was done about a finding."""

    SKIP_LINE = "skip_line"
    DROP_FACE = "drop_face"
    FAN_TRIANGULATE = "fan_triangulate"
    REMESH = "remesh"
    MERGE = "merge"
    REMOVE = "remove"
    FILL_HOLE = "fill_hole"
    SPLIT_EDGE = "split_edge"
    DUPLICATE_VERTEX = "duplicate_vertex"
    REORIENT = "reorient"
    FLIP = "flip"
    TRUNCATE = "truncate"
    FALLBACK = "fallback"
    ABORT = "abort"
    #: A name an output format cannot hold, written as one it can.
    RENAME = "rename"
    #: A tet rounding turned over, fixed by moving a constructed vertex
    #: within its rounding or by a flip.
    UNTANGLE = "untangle"
    #: Repair rung 2: a surface that intersects itself or whose shells
    #: overlap, cut along its intersections and replaced by the boundary of
    #: the solid it encloses.
    RESOLVE = "resolve"
    #: Repair rung 3: a surface the gate still rejects (open, non-manifold,
    #: intersecting where rung 2 could not resolve it), replaced by a wrap:
    #: a closed, manifold, intersection-free surface at a small offset
    #: around it.
    WRAP = "wrap"


@dataclass(frozen=True)
class Attempt:
    """One meshing attempt, in the order tried."""

    kernel: Kernel

    outcome: AttemptOutcome

    part: Optional[str] = None

    remeshed: Optional[bool] = None
    """The attempt meshed the remeshed surface, not the repaired one."""

    stderr: Optional[str] = None

    stdout: Optional[str] = None

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> Attempt:
        return cls(
            kernel=Kernel(v["kernel"]),
            outcome=AttemptOutcome(v["outcome"]),
            part=v.get("part", None),
            remeshed=v.get("remeshed", None),
            stderr=v.get("stderr", None),
            stdout=v.get("stdout", None),
        )


class AttemptOutcome(str, Enum):
    """AttemptOutcome"""

    SUCCEEDED = "succeeded"
    FAILED = "failed"


@dataclass(frozen=True)
class CodeCount:
    """Per-code totals in the summary."""

    reported: int
    """Occurrences with a detailed record in `records`."""

    total: int
    """Every occurrence found."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> CodeCount:
        return cls(
            reported=v["reported"],
            total=v["total"],
        )


class DefectCode(str, Enum):
    """Stable defect codes. Serialized as `SCREAMING_SNAKE_CASE` strings."""

    #: Malformed or non-numeric OBJ line, or a .meshb vertex that is not finite
    PARSE_BAD_LINE = "PARSE_BAD_LINE"
    #: Face references a missing vertex
    PARSE_BAD_INDEX = "PARSE_BAD_INDEX"
    #: Input file cannot be opened or read
    PARSE_READ_FAILED = "PARSE_READ_FAILED"
    #: Project file cannot be read or is not a valid project
    PARSE_PROJECT_INVALID = "PARSE_PROJECT_INVALID"
    #: No faces survive parsing
    PARSE_NO_GEOMETRY = "PARSE_NO_GEOMETRY"
    #: A binary input file (.meshb) is malformed or truncated
    PARSE_FORMAT_INVALID = "PARSE_FORMAT_INVALID"
    #: Lines, points, curves or surfaces in the OBJ file, or elements other than triangles in a .meshb file
    PARSE_UNSUPPORTED_ELEMENT = "PARSE_UNSUPPORTED_ELEMENT"
    #: Polygon could not be ear-clipped and was fan-triangulated
    PARSE_POLYGON_FALLBACK = "PARSE_POLYGON_FALLBACK"
    #: No parts selected, or a named part does not exist
    PART_EMPTY_SELECTION = "PART_EMPTY_SELECTION"
    #: Two selected parts intersect
    PART_OVERLAP = "PART_OVERLAP"
    #: Two selected parts touch along a surface
    PART_SHARED_INTERFACE = "PART_SHARED_INTERFACE"
    #: Vertices closer than tolerance
    GEO_DUPLICATE_VERTEX = "GEO_DUPLICATE_VERTEX"
    #: Zero-area or collinear face
    GEO_DEGENERATE_FACE = "GEO_DEGENERATE_FACE"
    #: Same face listed twice
    GEO_DUPLICATE_FACE = "GEO_DUPLICATE_FACE"
    #: Hole: edge with one adjacent face
    TOP_BOUNDARY_EDGE = "TOP_BOUNDARY_EDGE"
    #: Edge shared by three or more faces
    TOP_NONMANIFOLD_EDGE = "TOP_NONMANIFOLD_EDGE"
    #: Vertex joining separate fans
    TOP_NONMANIFOLD_VERTEX = "TOP_NONMANIFOLD_VERTEX"
    #: Neighbouring faces disagree on winding
    TOP_INCONSISTENT_ORIENTATION = "TOP_INCONSISTENT_ORIENTATION"
    #: Shell with negative signed volume
    TOP_INVERTED_SHELL = "TOP_INVERTED_SHELL"
    #: Surface encloses no volume
    TOP_OPEN_SHEET = "TOP_OPEN_SHEET"
    #: Faces of one shell intersect
    INT_SELF_INTERSECTION = "INT_SELF_INTERSECTION"
    #: Two shells of one part intersect
    INT_SHELL_OVERLAP = "INT_SHELL_OVERLAP"
    #: Minimum angle below threshold
    QUAL_SLIVER_FACE = "QUAL_SLIVER_FACE"
    #: Bounding box suggests wrong units
    QUAL_SUSPECT_SCALE = "QUAL_SUSPECT_SCALE"
    #: A feature or smooth edge override took no effect
    FEAT_OVERRIDE_IGNORED = "FEAT_OVERRIDE_IGNORED"
    #: A size control took no effect
    SIZE_CONTROL_IGNORED = "SIZE_CONTROL_IGNORED"
    #: Constrained Delaunay could not recover a facet
    MESH_RECOVERY_FAILED = "MESH_RECOVERY_FAILED"
    #: Steiner-point or element budget hit
    MESH_BUDGET_EXCEEDED = "MESH_BUDGET_EXCEEDED"
    #: Part too large for 32-bit vertex and face indices
    GEO_CAPACITY_EXCEEDED = "GEO_CAPACITY_EXCEEDED"
    #: A repair changed the shell count or enclosed volume beyond its limit
    REPAIR_GUARD_EXCEEDED = "REPAIR_GUARD_EXCEEDED"
    #: --strict was given but no conformal mesh could be produced
    MESH_STRICT_UNMET = "MESH_STRICT_UNMET"
    #: Worker panicked, timed out or ran out of memory
    MESH_WORKER_FAILED = "MESH_WORKER_FAILED"
    #: A tet that rounding to f64 turned over or flattened was untangled
    MESH_ROUNDING_REPAIRED = "MESH_ROUNDING_REPAIRED"
    #: Tets that rounding to f64 turned over or flattened could not be untangled
    MESH_ROUNDING_FAILED = "MESH_ROUNDING_FAILED"
    #: Tetrahedron with non-positive volume
    VOL_INVERTED_TET = "VOL_INVERTED_TET"
    #: Mesh face held by more than two tets, or by two on the same side
    VOL_NONMANIFOLD = "VOL_NONMANIFOLD"
    #: Non-orthogonality, skewness or volume ratio out of bounds
    VOL_LOW_QUALITY = "VOL_LOW_QUALITY"
    #: Hausdorff distance to input over tolerance
    VOL_BOUNDARY_DEVIATION = "VOL_BOUNDARY_DEVIATION"
    #: Neighbouring cells' sizes differ beyond the gradation limit
    VOL_GRADATION = "VOL_GRADATION"
    #: Part or set name over the Exodus limit
    EXO_NAME_TRUNCATED = "EXO_NAME_TRUNCATED"
    #: Mesh too large for the 64-bit-offset format
    EXO_CAPACITY_EXCEEDED = "EXO_CAPACITY_EXCEEDED"
    #: A part or patch name was changed to one the output format accepts
    OUT_NAME_CHANGED = "OUT_NAME_CHANGED"
    #: An output file could not be written
    OUT_WRITE_FAILED = "OUT_WRITE_FAILED"


@dataclass(frozen=True)
class FeatureSummary:
    """The feature lines of a part's surface: edges where it turns by more than
    the feature angle, where OBJ groups meet, and those the user forced, as
    the remesher keeps them. `features.vtp` holds the lines.
    """

    angle_edges: int
    """Of them, edges turning by more than the angle."""

    corners: int
    """Feature vertices not on one smooth line (they never move)."""

    curves: int
    """Chains of feature edges between corners (or closed loops)."""

    edges: int
    """Feature edges in all."""

    feature_angle_degrees: float
    """The angle used, degrees."""

    forced_edges: int
    """Edges forced by `feature_edges`."""

    group_edges: int
    """Edges between faces of different OBJ groups (fill faces are a group
    of their own).
    """

    length: float
    """Total length of the feature edges, in model units."""

    smoothed_edges: int
    """Edges `smooth_edges` turned smooth."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> FeatureSummary:
        return cls(
            angle_edges=v["angle_edges"],
            corners=v["corners"],
            curves=v["curves"],
            edges=v["edges"],
            feature_angle_degrees=float(v["feature_angle_degrees"]),
            forced_edges=v["forced_edges"],
            group_edges=v["group_edges"],
            length=float(v["length"]),
            smoothed_edges=v["smoothed_edges"],
        )


@dataclass(frozen=True)
class InputInfo:
    """InputInfo"""

    path: str
    """The input path exactly as given."""

    stats: Optional[InputStats] = None
    """Absent when the file could not be read."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> InputInfo:
        return cls(
            path=v["path"],
            stats=(None if v["stats"] is None else InputStats._parse(v["stats"])) if "stats" in v else None,
        )


@dataclass(frozen=True)
class InputStats:
    """What the parser found in the input."""

    bytes: int

    lines: int

    polygons: int
    """Faces as written in the file (polygons of any size)."""

    triangles: int
    """Triangles after triangulating the polygons."""

    vertices: int

    groups: Tuple[str, ...] = ()
    """`g` names in order of first appearance."""

    objects: Tuple[str, ...] = ()
    """`o` names in order of first appearance."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> InputStats:
        return cls(
            bytes=v["bytes"],
            lines=v["lines"],
            polygons=v["polygons"],
            triangles=v["triangles"],
            vertices=v["vertices"],
            groups=(tuple(v["groups"])) if "groups" in v else (),
            objects=(tuple(v["objects"])) if "objects" in v else (),
        )


@dataclass(frozen=True)
class JobConfig:
    """User-facing job settings (the policy flags). Recorded in the report."""

    element_size: float = 0.03
    """Target element size (the stuffing lattice spacing), relative to the
    part's bounding-box diagonal.
    """

    feature_angle: Optional[float] = None
    """The feature angle, degrees: edges where the surface turns by more
    than this are feature lines, which remeshing keeps and the report
    lists. Unset: 40° for remeshing and the report, 45° for the size
    field's curvature term; set, both use it.
    """

    feature_edges: Tuple[Tuple[int, int], ...] = ()
    """Edges made feature lines whatever the angle, each a pair of 0-based
    indices into the parsed vertex list (as in report locations).
    """

    formats: Tuple[OutputFormat, ...] = ()
    """The mesh files written: by default Exodus II, Gmsh and VTU."""

    kernel: Optional[KernelChoice] = None
    """The volume kernel tried first."""

    max_deviation: float = 0.05
    """Largest acceptable fallback deviation, relative to the part's bounding-box diagonal."""

    no_remesh: Optional[bool] = None
    """Never remesh a part's surface (by default it is remeshed when its
    conformal mesh fails or its cells are beyond the quality limits).
    """

    no_repair: bool = False
    """Diagnose and report only."""

    parts: Tuple[str, ...] = ()
    """Names or indices of the parts to mesh; empty means all."""

    sizes: Tuple[SizeControl, ...] = ()
    """Element sizes per part and per patch (OBJ group), overriding
    `element_size` there.
    """

    smooth_edges: Tuple[Tuple[int, int], ...] = ()
    """Edges that are not feature lines although they turn by more than the
    feature angle, as pairs like `feature_edges` (a group boundary stays
    a feature line).
    """

    strict: bool = False
    """Conformal paths only; never use the stuffing fallback."""

    stuffing_only: Optional[bool] = None
    """Skip the conformal kernel and mesh with stuffing only. Set by the
    client when a worker running the conformal kernel timed out or ran out
    of memory, for its one retry (`MESH_WORKER_FAILED`, fallback).
    """

    uniform_size: Optional[bool] = None
    """Mesh at `element_size` everywhere, without the sizing field (which
    grades the size down to a quarter of it where the surface is curved
    or the part is thin).
    """

    wrap_alpha: float = 1.0
    """The wrap's α (repair rung 3: openings narrower than about 2α are
    closed), relative to the element size.
    """

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> JobConfig:
        return cls(
            element_size=(float(v["element_size"])) if "element_size" in v else 0.03,
            feature_angle=(None if v["feature_angle"] is None else float(v["feature_angle"])) if "feature_angle" in v else None,
            feature_edges=(tuple(tuple(e0) for e0 in v["feature_edges"])) if "feature_edges" in v else (),
            formats=(tuple(OutputFormat(e0) for e0 in v["formats"])) if "formats" in v else (),
            kernel=(KernelChoice(v["kernel"])) if "kernel" in v else None,
            max_deviation=(float(v["max_deviation"])) if "max_deviation" in v else 0.05,
            no_remesh=v.get("no_remesh", None),
            no_repair=v.get("no_repair", False),
            parts=(tuple(v["parts"])) if "parts" in v else (),
            sizes=(tuple(SizeControl._parse(e0) for e0 in v["sizes"])) if "sizes" in v else (),
            smooth_edges=(tuple(tuple(e0) for e0 in v["smooth_edges"])) if "smooth_edges" in v else (),
            strict=v.get("strict", False),
            stuffing_only=v.get("stuffing_only", None),
            uniform_size=v.get("uniform_size", None),
            wrap_alpha=(float(v["wrap_alpha"])) if "wrap_alpha" in v else 1.0,
        )


class Kernel(str, Enum):
    """Volume kernel used by a meshing attempt."""

    DELAUNAY = "delaunay"
    STUFFING = "stuffing"
    #: The exact constrained Delaunay kernel, rounded to f64.
    EXACT = "exact"


class KernelChoice(str, Enum):
    """The volume kernel tried first; the rest of the fallback ladder follows it."""

    #: The conforming Delaunay kernel, refined to the element size.
    CONFORMAL = "conformal"
    #: The exact constrained Delaunay kernel (`--kernel exact`): the input
    #: surface tetrahedralized exactly, then rounded to f64, without
    #: refinement; the conforming kernel and stuffing follow if it fails.
    EXACT = "exact"


@dataclass(frozen=True)
class Location:
    """Where a finding is. `vertices` are 0-based indices into the parsed vertex
    list and `faces` are 0-based indices into the triangulated face list (a
    polygon with n corners becomes n - 2 consecutive triangles); `obj_lines` are
    1-based line numbers in the input file. Each list holds at most
    [`MAX_LOCATION_ITEMS`] entries; a truncated list has its full length in the
    record's `measured` map (`vertices_total`, `faces_total`, `obj_lines_total`).
    """

    faces: Tuple[int, ...] = ()

    obj_lines: Tuple[int, ...] = ()

    point: Optional[Tuple[float, float, float]] = None

    vertices: Tuple[int, ...] = ()

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> Location:
        return cls(
            faces=(tuple(v["faces"])) if "faces" in v else (),
            obj_lines=(tuple(v["obj_lines"])) if "obj_lines" in v else (),
            point=(None if v["point"] is None else tuple(float(e0) for e0 in v["point"])) if "point" in v else None,
            vertices=(tuple(v["vertices"])) if "vertices" in v else (),
        )


@dataclass(frozen=True)
class MeshSummary:
    """The volume mesh of one part and its quality."""

    accepted: bool
    """Whether the mesh passed the volume gate and was written."""

    element_size: float
    """Element size used (absolute): the largest target size."""

    hausdorff: float
    """Two-sided Hausdorff distance between the input surface and the mesh
    boundary (sampled at vertices and face centroids), absolute and
    relative to the part's bounding-box diagonal.
    """

    hausdorff_relative: float

    kernel: Kernel

    max_dihedral_degrees: float

    max_size_ratio: float
    """Largest ratio between face neighbours' sizes (each cell's longest
    edge), against `VOL_GRADATION`'s limit.
    """

    min_dihedral_degrees: float

    nodes: int

    surface_volume: float
    """Signed volume enclosed by the repaired surface (for a surface that
    failed the gate, overlapping shells count their overlap twice).
    """

    tets: int

    volume: float

    dihedral_histogram: Tuple[int, ...] = ()
    """Every dihedral angle (six per tet) in 10° bins from [0°, 10°) to [170°, 180°]."""

    max_non_orthogonality_degrees: float = 0.0
    """Finite-volume quality on interior faces: angle between face normal and
    centroid line, its mean, skewness, and neighbour volume ratio.
    """

    max_skewness: float = 0.0

    max_volume_ratio: float = 0.0

    mean_non_orthogonality_degrees: float = 0.0

    min_element_size: Optional[float] = None
    """Smallest target size (absolute), when the size was graded: the
    sizing field's for a conformal mesh, the surface cell size for graded
    stuffing.
    """

    remeshed: Optional[bool] = None
    """The mesh conforms to the remeshed surface, not the repaired one (the
    conformal mesh of the repaired surface failed); its deviation from
    the repaired surface is `hausdorff`.
    """

    steiner_points: int = 0
    """Points the kernel added (Delaunay: boundary and interior Steiner points)."""

    worst_cells: Tuple[WorstCell, ...] = ()
    """Where the worst cells are: the smallest dihedral angle, and the
    largest non-orthogonality and skewness over interior faces and,
    separately, over the faces of wall-adjacent cells.
    """

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> MeshSummary:
        return cls(
            accepted=v["accepted"],
            element_size=float(v["element_size"]),
            hausdorff=float(v["hausdorff"]),
            hausdorff_relative=float(v["hausdorff_relative"]),
            kernel=Kernel(v["kernel"]),
            max_dihedral_degrees=float(v["max_dihedral_degrees"]),
            max_size_ratio=float(v["max_size_ratio"]),
            min_dihedral_degrees=float(v["min_dihedral_degrees"]),
            nodes=v["nodes"],
            surface_volume=float(v["surface_volume"]),
            tets=v["tets"],
            volume=float(v["volume"]),
            dihedral_histogram=(tuple(v["dihedral_histogram"])) if "dihedral_histogram" in v else (),
            max_non_orthogonality_degrees=(float(v["max_non_orthogonality_degrees"])) if "max_non_orthogonality_degrees" in v else 0.0,
            max_skewness=(float(v["max_skewness"])) if "max_skewness" in v else 0.0,
            max_volume_ratio=(float(v["max_volume_ratio"])) if "max_volume_ratio" in v else 0.0,
            mean_non_orthogonality_degrees=(float(v["mean_non_orthogonality_degrees"])) if "mean_non_orthogonality_degrees" in v else 0.0,
            min_element_size=(None if v["min_element_size"] is None else float(v["min_element_size"])) if "min_element_size" in v else None,
            remeshed=v.get("remeshed", None),
            steiner_points=v.get("steiner_points", 0),
            worst_cells=(tuple(WorstCell._parse(e0) for e0 in v["worst_cells"])) if "worst_cells" in v else (),
        )


class Outcome(str, Enum):
    """Result of the action."""

    #: The action was carried out.
    APPLIED = "applied"
    #: The action was attempted and did not succeed.
    FAILED = "failed"
    #: Left for a later stage to decide.
    DEFERRED = "deferred"


class OutputFormat(str, Enum):
    """A mesh file format the worker writes."""

    #: Exodus II (`mesh.exo`): a `TETRA4` block per part, side sets.
    EXODUS = "exodus"
    #: Gmsh MSH 4.1 (`mesh.msh`).
    GMSH = "gmsh"
    #: VTK unstructured grid (`mesh.vtu`).
    VTU = "vtu"
    #: An OpenFOAM case (`openfoam/`): an ASCII polyMesh with a patch per
    #: boundary patch and a cell zone per part. Not written by default.
    OPENFOAM = "openfoam"
    #: A Fluent mesh file (`mesh_fluent.msh`, ASCII): a cell zone per part,
    #: an interior face zone, a wall face zone per boundary patch. Not
    #: written by default.
    FLUENT = "fluent"
    #: A libMeshb binary file (`mesh.meshb`): tetrahedra referenced by part,
    #: boundary triangles by patch, with `mesh_refs.txt` naming each
    #: reference. Not written by default.
    MESHB = "meshb"


class PartSource(str, Enum):
    """How a part was formed."""

    #: An `o` statement (faces before any `o` form the part `default`).
    OBJECT = "object"
    #: A `g` statement, in a file without `o` statements.
    GROUP = "group"
    #: A connected component, in a file with neither; enclosed cavities join
    #: the component that contains them.
    COMPONENT = "component"


@dataclass(frozen=True)
class PartSummary:
    """One part, as the GUI part tree shows it."""

    bounds: Tuple[Tuple[float, float, float], Tuple[float, float, float]]
    """Bounding box `[min, max]`."""

    index: int
    """1-based; `--parts` accepts this or the name."""

    name: str

    selected: bool

    shells: int
    """Edge-connected surface pieces."""

    source: PartSource

    triangles: int

    vertices: int
    """Distinct vertices after merging duplicates."""

    features: Optional[FeatureSummary] = None
    """The feature lines of the repaired surface (the final one: the wrap,
    for a wrapped part), present when the part was repaired.
    """

    mesh: Optional[MeshSummary] = None
    """Present when a volume mesh was attempted."""

    repair: Optional[RepairSummary] = None
    """Present when the part was repaired."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> PartSummary:
        return cls(
            bounds=tuple(tuple(float(e1) for e1 in e0) for e0 in v["bounds"]),
            index=v["index"],
            name=v["name"],
            selected=v["selected"],
            shells=v["shells"],
            source=PartSource(v["source"]),
            triangles=v["triangles"],
            vertices=v["vertices"],
            features=(None if v["features"] is None else FeatureSummary._parse(v["features"])) if "features" in v else None,
            mesh=(None if v["mesh"] is None else MeshSummary._parse(v["mesh"])) if "mesh" in v else None,
            repair=(None if v["repair"] is None else RepairSummary._parse(v["repair"])) if "repair" in v else None,
        )


class QualityMeasure(str, Enum):
    """A quality measure whose worst value is located."""

    #: A cell's smallest dihedral angle, in degrees (the worst is smallest).
    MIN_DIHEDRAL = "min_dihedral"
    #: A face's angle between its normal and the line joining the cell
    #: centroids, in degrees (on the boundary, the line from the cell's
    #: centroid to the face's).
    NON_ORTHOGONALITY = "non_orthogonality"
    #: A face's distance from its centroid to where the centroid line
    #: crosses it, over the centroid distance.
    SKEWNESS = "skewness"
    #: The ratio of the sizes (longest edges) of a face's two cells.
    SIZE_RATIO = "size_ratio"


@dataclass(frozen=True)
class Record:
    """One finding."""

    code: DefectCode

    severity: Severity

    stage: Stage

    action: Optional[Action] = None

    hint: Optional[str] = None

    location: Optional[Location] = None

    measured: Dict[str, Any] = field(default_factory=dict)

    outcome: Optional[Outcome] = None

    part: Optional[str] = None
    """Part the finding belongs to, once parts exist."""

    remedies: Tuple[Remedy, ...] = ()
    """Settings changes to try next, best first (see the report reference);
    filled in when the report is finished.
    """

    threshold: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> Record:
        return cls(
            code=DefectCode(v["code"]),
            severity=Severity(v["severity"]),
            stage=Stage(v["stage"]),
            action=(None if v["action"] is None else Action(v["action"])) if "action" in v else None,
            hint=v.get("hint", None),
            location=(Location._parse(v["location"])) if "location" in v else None,
            measured=(dict(v["measured"])) if "measured" in v else {},
            outcome=(None if v["outcome"] is None else Outcome(v["outcome"])) if "outcome" in v else None,
            part=v.get("part", None),
            remedies=(tuple(Remedy._parse(e0) for e0 in v["remedies"])) if "remedies" in v else (),
            threshold=(dict(v["threshold"])) if "threshold" in v else {},
        )


@dataclass(frozen=True)
class Remedy:
    """One suggested next step: set `setting` to `value`, or multiply it by
    `scale` (exactly one of the two).
    """

    reason: str
    """Why it should help, in one line."""

    setting: Setting
    """The setting to change."""

    scale: Optional[float] = None
    """A factor to multiply its value by (a numeric setting)."""

    value: Optional[Any] = None
    """Its new value, as in `config`."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> Remedy:
        return cls(
            reason=v["reason"],
            setting=Setting(v["setting"]),
            scale=(None if v["scale"] is None else float(v["scale"])) if "scale" in v else None,
            value=v.get("value", None),
        )


@dataclass(frozen=True)
class RepairSummary:
    """Face and vertex counts through repair, and the surface gate's verdict."""

    faces_after: int

    faces_before: int
    """Triangles in the part as parsed."""

    gate_passed: bool

    vertices_after: int

    vertices_before: int
    """Vertices entering repair (after merging duplicates)."""

    components_after: int = 0

    components_before: int = 0
    """Connected pieces of the surface (by shared vertices) and its Euler
    characteristic (vertices − edges + faces; 2 per sphere-like shell),
    entering repair and after it.
    """

    euler_after: int = 0

    euler_before: int = 0

    rung: int = 1
    """The highest repair rung that changed the surface: 1 for local fixes
    (or none), 2 when its intersections were resolved.
    """

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> RepairSummary:
        return cls(
            faces_after=v["faces_after"],
            faces_before=v["faces_before"],
            gate_passed=v["gate_passed"],
            vertices_after=v["vertices_after"],
            vertices_before=v["vertices_before"],
            components_after=v.get("components_after", 0),
            components_before=v.get("components_before", 0),
            euler_after=v.get("euler_after", 0),
            euler_before=v.get("euler_before", 0),
            rung=v.get("rung", 1),
        )


class Setting(str, Enum):
    """The setting a remedy changes: a `config` field, by name."""

    #: `parts`.
    PARTS = "parts"
    #: `strict`.
    STRICT = "strict"
    #: `no_remesh`.
    NO_REMESH = "no_remesh"
    #: `kernel`.
    KERNEL = "kernel"
    #: `element_size`.
    ELEMENT_SIZE = "element_size"
    #: `max_deviation`.
    MAX_DEVIATION = "max_deviation"
    #: `wrap_alpha`.
    WRAP_ALPHA = "wrap_alpha"


class Severity(str, Enum):
    """Severity"""

    #: Noteworthy, no action needed.
    INFO = "info"
    #: Repaired automatically; output still valid.
    WARNING = "warning"
    #: Blocks a stage; a fallback may still succeed.
    ERROR = "error"
    #: No usable output possible.
    FATAL = "fatal"


@dataclass(frozen=True)
class SizeControl:
    """One size control (`[[settings.sizes]]` in the project file)."""

    element_size: float
    """The element size, relative to the part's bounding-box diagonal. For a
    whole part it replaces `element_size` (the last such control listed
    wins); for a patch it is the size on its faces, graded from there,
    and only refines: it cannot exceed the part's size.
    """

    group: Optional[str] = None
    """The OBJ group (patch); absent: the whole part."""

    part: Optional[str] = None
    """The part, by name or 1-based index; absent: every part."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> SizeControl:
        return cls(
            element_size=float(v["element_size"]),
            group=v.get("group", None),
            part=v.get("part", None),
        )


class Stage(str, Enum):
    """Pipeline stage that produced a record."""

    PARSE = "parse"
    DIAGNOSE = "diagnose"
    SELECT = "select"
    REPAIR = "repair"
    SURFACE_GATE = "surface_gate"
    SIZING = "sizing"
    VOLUME_MESH = "volume_mesh"
    VOLUME_GATE = "volume_gate"
    EMIT = "emit"
    #: The worker process itself (crash, timeout, protocol failure).
    WORKER = "worker"


@dataclass(frozen=True)
class Summary:
    """Counts and the resulting exit code."""

    by_code: Dict[str, CodeCount]

    errors: int

    exit_code: int

    fatal: int

    info: int

    warnings: int

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> Summary:
        return cls(
            by_code={k0: CodeCount._parse(e0) for k0, e0 in v["by_code"].items()},
            errors=v["errors"],
            exit_code=v["exit_code"],
            fatal=v["fatal"],
            info=v["info"],
            warnings=v["warnings"],
        )


@dataclass(frozen=True)
class ToolInfo:
    """ToolInfo"""

    name: str

    version: str

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> ToolInfo:
        return cls(
            name=v["name"],
            version=v["version"],
        )


@dataclass(frozen=True)
class WorstCell:
    """The worst cell or face for one measure."""

    element: int
    """The cell: its 1-based element number in the part's element block.
    For a face measure, the first of the face's two cells.
    """

    measure: QualityMeasure

    point: Tuple[float, float, float]
    """The cell's centroid, or the face's for a face measure."""

    value: float

    neighbour: Optional[int] = None
    """For an interior face, the cell across it."""

    wall: Optional[bool] = None
    """Over the wall-adjacent cells only (cells with a boundary face; their
    boundary faces count for non-orthogonality).
    """

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> WorstCell:
        return cls(
            element=v["element"],
            measure=QualityMeasure(v["measure"]),
            point=tuple(float(e0) for e0 in v["point"]),
            value=float(v["value"]),
            neighbour=v.get("neighbour", None),
            wall=v.get("wall", None),
        )


@dataclass(frozen=True)
class Report:
    """The complete `report.json`."""

    attempts: Tuple[Attempt, ...]

    config: JobConfig

    input: InputInfo

    records: Tuple[Record, ...]

    schema_version: int

    stages_run: Tuple[Stage, ...]
    """Stages that ran, in order."""

    summary: Summary

    tool: ToolInfo

    parts: Tuple[PartSummary, ...] = ()
    """Parts found by diagnosis, in index order."""

    @classmethod
    def _parse(cls, v: Mapping[str, Any]) -> Report:
        return cls(
            attempts=tuple(Attempt._parse(e0) for e0 in v["attempts"]),
            config=JobConfig._parse(v["config"]),
            input=InputInfo._parse(v["input"]),
            records=tuple(Record._parse(e0) for e0 in v["records"]),
            schema_version=v["schema_version"],
            stages_run=tuple(Stage(e0) for e0 in v["stages_run"]),
            summary=Summary._parse(v["summary"]),
            tool=ToolInfo._parse(v["tool"]),
            parts=(tuple(PartSummary._parse(e0) for e0 in v["parts"])) if "parts" in v else (),
        )
