# The report (`report.json`)

Every run writes `report.json` to its output directory, whatever happens: a clean mesh, a repaired one, a fallback, or a failure with no mesh at all. It records every finding, every repair, every meshing attempt and the quality of every mesh written, so a script can decide what to do without reading the logs. The layout is versioned (`schema_version`, currently 1), and its JSON schema is [`schema/report.schema.json`](../schema/report.schema.json), also printed by `mesher schema`.

The same input and settings give a byte-identical report.

## Top level

| Field | Contents |
| --- | --- |
| `schema_version` | Layout version; bumped on any breaking change |
| `tool` | `name` and `version` of the mesher |
| `input` | `path` as given, and `stats` (absent if the file could not be read): `bytes`, `lines`, `vertices`, `polygons` (faces as written), `triangles` (after triangulating polygons), `objects` and `groups` (names in order of appearance; a `.meshb` input has no objects or lines, and its groups are its triangles' references, `ref<k>`, which name patches, not parts) |
| `config` | The job's settings; see [Settings](#settings) |
| `stages_run` | Stages that ran, in order: `parse`, `diagnose`, `select`, `repair`, `surface_gate`, `volume_mesh`, `volume_gate`, `emit`; empty when the worker failed before reporting (its record has stage `worker`) |
| `parts` | One entry per part found, in index order; see [Parts](#parts) |
| `records` | The findings; see [Records](#records) |
| `attempts` | Every meshing attempt, in the order tried; see [Attempts](#attempts) |
| `summary` | Counts by severity and by code, and the exit code; see [Summary and exit code](#summary-and-exit-code) |

## Settings

`config` holds the settings the job ran with (the CLI flags or the project file's `[settings]`):

| Field | Meaning |
| --- | --- |
| `parts` | Names or 1-based indices of the parts to mesh; empty means all |
| `no_repair` | Diagnose and report only |
| `strict` | Conformal meshes only; never use the stuffing fallback |
| `max_deviation` | Largest acceptable distance between the mesh boundary and the surface, relative to the part's bounding-box diagonal (default 0.05) |
| `element_size` | Target element size, relative to the part's bounding-box diagonal (default 0.03): the largest size anywhere |
| `uniform_size` | Present (true) when the size was not graded |
| `no_remesh` | Present (true) when remeshing a part's surface was turned off |
| `wrap_alpha` | The wrap's α for repair rung 3, relative to `element_size` (default 1): openings narrower than about 2α are closed |
| `stuffing_only` | Present (true) when the client retried a worker that timed out or ran out of memory with the stuffing fallback only |
| `feature_angle` | Present when set (`--feature-angle`): the feature angle in degrees, between 0 and 180. Edges where the surface turns by more than it are feature lines, which remeshing keeps (their vertices move only along them, corners never) and `features` lists. Unset, remeshing and `features` use 40° and the size field's curvature term 45°; set, all use it |
| `feature_edges` | Present when set (project file only): edges made feature lines whatever the angle, each `[a, b]`, 0-based indices into the parsed vertex list as in record locations; every repaired edge between vertices that came from `a` and `b` counts. One that matches no edge of a selected part is `FEAT_OVERRIDE_IGNORED` |
| `sizes` | Present when set (project file only, `[[settings.sizes]]`): size controls, each with `element_size` (relative to the part's bounding-box diagonal, like `element_size`) and `part` (a name or 1-based index; absent: every part) and/or `group` (an OBJ group, the patch). A control with a part and no group sets that part's element size (the last listed wins); one with a group sets the size on that group's faces, graded away from them (it only refines: the floor of a quarter of the size moves down to it, and it cannot exceed the part's size). Under `uniform_size` only the patch sizes grade the field. One that takes no effect is `SIZE_CONTROL_IGNORED` |
| `formats` | Present when not the default (`--formats`, or the project file): the mesh files written, of `exodus` (`mesh.exo`), `gmsh` (`mesh.msh`), `vtu` (`mesh.vtu`), `openfoam` (the case directory `openfoam/`), `fluent` (`mesh_fluent.msh`) and `meshb` (`mesh.meshb`, `mesh_refs.txt`). Default `exodus`, `gmsh` and `vtu` |
| `smooth_edges` | Present when set (project file only): edges, as for `feature_edges`, that are not feature lines although they turn by more than the angle. A group boundary stays a feature line (`FEAT_OVERRIDE_IGNORED`) |

## Parts

Each part is an OBJ object (`o`), or in a file without objects an OBJ group (`g`), or in a file with neither a connected component (`source`: `object`, `group` or `component`).

| Field | Meaning |
| --- | --- |
| `index` | 1-based; `--parts` accepts it or the name |
| `name`, `source` | As above |
| `triangles`, `vertices` | After triangulating and merging duplicate vertices |
| `shells` | Edge-connected surface pieces |
| `bounds` | Bounding box `[min, max]` |
| `selected` | Whether the job meshed it |
| `repair` | `faces_before`, `faces_after`, `vertices_before`, `vertices_after`, and `gate_passed`: whether the repaired surface is a valid closed solid; `rung`: the highest repair rung that changed the surface (1 for local fixes or none, 2 when its intersections were resolved, 3 when it was wrapped as the last resort); `components_before`, `components_after`, `euler_before`, `euler_after`: the surface's connected pieces and Euler characteristic (vertices − edges + faces, 2 per sphere-like shell, 0 for a ring) entering repair and after it |
| `mesh` | Present when a volume mesh was attempted; see [Mesh quality](#mesh-quality) |
| `features` | Present when the part was repaired: the feature lines of its final surface (the wrap, for a wrapped part), as the remesher keeps them. `feature_angle_degrees` (the angle used); `edges` in all, of them `angle_edges` (turning by more than the angle, or with a degenerate face), `group_edges` (between faces of different OBJ groups; fill faces are a group of their own) and `forced_edges` (`feature_edges`); `smoothed_edges` (`smooth_edges` that took effect); `corners` (feature vertices not on one smooth line: never moved); `curves` (chains of feature edges between corners, or closed loops); `length` (total, in model units). The lines are in `features.vtp` |

## Records

A record is one finding:

```json
{
  "code": "TOP_NONMANIFOLD_EDGE",
  "severity": "error",
  "stage": "diagnose",
  "part": "body",
  "location": {
    "vertices": [1042, 1043],
    "faces": [2201, 2209, 3310],
    "obj_lines": [1088, 1096, 2197],
    "point": [0.412, -1.07, 3.3]
  },
  "measured": { "incident_faces": 3 },
  "threshold": { "max_incident_faces": 2 },
  "action": "split_edge",
  "outcome": "applied",
  "hint": "Three faces share this edge; likely two touching parts exported as one object."
}
```

| Field | Meaning |
| --- | --- |
| `code` | A stable defect code; see [Defect codes](#defect-codes) |
| `severity` | `info` (noteworthy, nothing to do), `warning` (repaired or worked around; the output is still valid), `error` (a stage failed; a fallback may still have succeeded), `fatal` (no usable output) |
| `stage` | The stage that found it |
| `part` | The part, once parts exist |
| `location` | Where it is: `vertices` (0-based indices into the parsed vertex list), `faces` (0-based indices into the triangulated face list; a polygon of n corners becomes n − 2 consecutive triangles), `obj_lines` (1-based lines of the input file), `point` (a position, when one stands for the finding). Each list holds at most 256 entries; a longer one keeps its full length in `measured` as `vertices_total`, `faces_total` or `obj_lines_total` |
| `measured`, `threshold` | What was measured, and the limit it was held against, by name |
| `action`, `outcome` | What was done (`skip_line`, `drop_face`, `fan_triangulate`, `remesh`, `merge`, `remove`, `fill_hole`, `split_edge`, `duplicate_vertex`, `reorient`, `flip`, `truncate`, `rename`, `fallback`, `abort`, `untangle`, `resolve`, `wrap`) and whether it was `applied`, `failed` or `deferred` to a later stage |
| `hint` | A sentence on what it means or what to try |

At most 100 records are kept per code and stage (a repair's and the surface gate's records stay when diagnosis has used its share); `summary.by_code` still counts every occurrence. A finding is never repaired without a record.

The same records are in `defects.vtp` (open it in ParaView beside the mesh) and are drawn by the GUI, which shows each record's OBJ lines.

### Records from repair rung 2

A closed part whose surface intersects itself or whose shells overlap is cut along its intersections and replaced by the boundary of the solid it encloses (exact arithmetic, then rounded to f64 and checked; nested shells inside the solid and interfaces between solids touching face to face go). One `INT_SELF_INTERSECTION` or `INT_SHELL_OVERLAP` record per code found, stage `repair`, severity `info`, action `resolve` / `applied`, with `rung` (2), `intersecting_pairs` and `coplanar_pairs` (pairs of input triangles that meet), `faces_before`, `faces_after`, `pieces_dropped` (pieces of the cut surface inside or outside the solid), `components_before`, `components_after`, `euler_before`, `euler_after`, `rounding_rounds` (times the rounded surface had to be resolved again); the first also `faces_removed`, `faces_added` and `vertices_added` (net). Diagnosis's findings of those codes become warnings with `resolve` / `applied`. When it cannot be done (the part is open, or rounding cannot make it clean), the records have `resolve` / `failed` and `reason`, and the part goes to the fallback. The part's `repair` summary carries `rung` and the counts too, since the per-code cap can leave a record out.

### Records from repair rung 3

The last resort: a part whose repaired surface failed the surface gate (a hole too large to fill, an edge whose faces do not pair up, an open part that intersects itself, rounding that could not be made clean), or passed it but did not mesh conformally (not under `strict`, which meshes such a surface as it is or fails), and whose stuffing fallback was rejected by the volume gate (or not allowed, under `strict`) is replaced by a wrap: a closed, oriented, manifold surface with no intersecting triangles, at a small offset δ around the surface, that closes openings narrower than about 2α (α = `wrap_alpha` × `element_size` × the part's diagonal, δ = α / 30). For an open part, the wrap closes holes across the level set where the surface's generalized winding number is 1/2 (the winding guard), instead of wrapping both sides of its walls. The wrap is kept when the faces bounding the solid the surface encloses (those with the winding number below 1/2 just off one side) are all within `max_deviation` of it. One record per gate check the surface failed (`TOP_BOUNDARY_EDGE`, `TOP_NONMANIFOLD_EDGE`, `TOP_NONMANIFOLD_VERTEX`, `TOP_INCONSISTENT_ORIENTATION`, `TOP_OPEN_SHEET`, `TOP_INVERTED_SHELL`, `INT_SELF_INTERSECTION`, `INT_SHELL_OVERLAP`; for a surface that passed the gate, one `MESH_RECOVERY_FAILED`), stage `repair`, severity `info`, action `wrap` / `applied`, with `rung` (3), `alpha` and `offset` (relative to the diagonal), `winding_guard`, `deviation` (the largest distance from those faces to the wrap, relative to the diagonal: features closed over), `faces_before`, `faces_after`, `components_before`, `components_after`, `euler_before`, `euler_after`, `cap_points` (wrap vertices on the winding guard's level set) and `cells_relabelled` (to make it manifold); the first also `faces_removed`, `faces_added` and `vertices_added` (net). The part's errors of those codes become warnings with `wrap` / `applied` (a repair that failed keeps its own action). A wrap triangle's `faces` location is the model triangle nearest its centre, so side sets keep their groups; its vertices have no model vertex. The wrapped part then goes through the volume pipeline from the start (the conformal kernel, then stuffing on the wrap), and its `repair` summary, and `repaired.obj`, are the wrap's; once it meshes, the records that sent the part to the last resort (the gate's, the rejected stuffing's) become warnings with `wrap` / `applied`. These records come after the rejected stuffing's. When the wrap fails or is too far from the surface, the records have `wrap` / `failed` and `reason`, and the part has no mesh.

### Records from meshing

These carry the measurements most often needed to act on them:

| Code | `measured` (and `threshold`) |
| --- | --- |
| `MESH_RECOVERY_FAILED` | `missing_faces`, or for a recovery stopped early because it grew instead of converging, `stalled` (`missing_grew`: the first recovery's missing subsegments and subfacets passed 8 × the input triangles + 1,024; `points_grew`: its new points passed 8 × the input points + 4,096; `cascade`: a later recovery's missing count passed 8 × its count at the start + 1,024) with `missing_pieces` against `max_missing_pieces` or `new_points` against `max_new_points`; after a remesh retry, `remeshed_triangles`, `remesh_rung` (the rung whose surface was kept, below) and, when remeshing merged features too small to mesh (corners closer than a tenth of the local size along a feature edge, and the two sides of strips that narrow), `merged_small_features`, and when other remeshings were tried, `other_remeshings`: each such rung in order, `name: why it was not kept` (it did not mesh, or its `worst cell` was worse; a surface that meshes with a worst cell under 0.25 does not end the ladder: the rungs below are tried until one meshes with a worst cell of at least 0.25, and the best kept; a checked rung runs only when the same rung unchecked did not mesh) (the rungs: `snapped`, which snaps narrow strips shut and collapses edges shorter than a tenth of the local size whatever lengths they leave; `snapped, checked`, with every remeshing operation checked against making the surface cross itself; `free collapses` and `free collapses, checked`, without snapping; `plain`, without either) (it succeeded: the record is then a warning with `remesh` / `applied`) or `remesh_retry` (why the retry did not help; for a remeshed surface that failed the surface gate, each failed check with its count, such as `3 self-intersecting triangle pairs`) |
| `MESH_BUDGET_EXCEEDED` | `points` against `max_points`, `work_steps` against `max_work_steps` (the constrained kernel), or `lattice_vertices` against `max_lattice_vertices` (stuffing) |
| `MESH_WORKER_FAILED` | `reason`: `timeout` (with `timeout_seconds`), `out_of_memory` (with `memory_limit_bytes`), `cancelled`, `spawn_failed` (with `io_error`), `died` (with `exit_status`) or `unreadable_result` (the worker finished but its result line did not parse, with `parse_error`; a bug, exit 70); `stderr_tail`. After a timeout or out-of-memory, the client retries with stuffing only and the record becomes a warning with `fallback` / `applied` if that produced a mesh |
| `MESH_ROUNDING_REPAIRED` | One record per fix of the exact kernel's rounding (`--kernel exact`): `fix` is `vertex_move` (a constructed point moved within 64 ulps per coordinate of its exact position, by `moved`; two records when two points move together), `flip_2_3` or `edge_removal`; `location.point` is where. Action `untangle` / `applied` |
| `MESH_ROUNDING_FAILED` | `tets_left` still turned over or flat in floating point after untangling, `turned_over_before`, `fixes`; `location.point` is the first one. The conformal kernel and stuffing follow |
| `VOL_LOW_QUALITY` | One record per measure beyond its limit: `min_dihedral_degrees` and `max_dihedral_degrees` against the kernel's bounds (10° and 170° for the constrained Delaunay kernel; 9.3171° and 161.6432°, the bounds Labelle and Shewchuk prove for its warping parameters, for stuffing); `max_non_orthogonality_degrees` against 70°; `max_skewness` against 2. Each is located at the worst cell: `element` (its 1-based number in the part's element block), for a face measure `neighbour` (the cell across the face), and the location's `point` (the cell's centroid, or the face's). An `info` record with `remesh/applied` instead says the surface was remeshed for quality: `input_worst_cell` and `remeshed_worst_cell` (1 is at the limits, smaller is worse), `remeshed_triangles`, and `remesh_rung`, `merged_small_features` and `other_remeshings` as for `MESH_RECOVERY_FAILED`. A warning with `fallback/applied` says every conformal mesh was far beyond the limits (worst cell under 0.25) and the stuffing mesh, whose cells were better, was kept: `conformal_worst_cell`, `stuffing_worst_cell` and `volume_error` (it must be at most 1 %), and what remeshing gave though not kept: `remeshed_worst_cell`, `remeshed_triangles`, `remesh_rung`, `merged_small_features` and `other_remeshings` as for `MESH_RECOVERY_FAILED`, or `remesh_retry` when no rung meshed; its boundary only approximates the surface (the mesh's `hausdorff`) |
| `VOL_BOUNDARY_DEVIATION` | `hausdorff` and `hausdorff_relative` against `max_deviation_relative` |
| `VOL_INVERTED_TET` | `tets` not positive by the exact predicate (a kernel bug), or `flat_in_float`: tets positive by the exact predicate but with no positive volume when computed in floating point, as a solver computes it (flat to rounding). Either rejects the candidate mesh, and the next step of the fallback ladder runs |
| `VOL_NONMANIFOLD` | `crowded_faces` (faces held by more than two tets) and `same_side_faces` (faces whose two tets lie on the same side, so they overlap), with `kernel`. With every tet positive and the boundary on the surface, a mesh with neither covers each point of the part once. A kernel bug: it rejects the candidate mesh, and the next step of the fallback ladder runs |

## Defect codes

Severity is the code's default; records can be upgraded or downgraded where the outcome warrants it. A conformal failure that a fallback recovered is a warning. So are a part's errors from diagnosis, repair and the surface gate (a self-intersection, an overlap, a hole too large to fill, a non-manifold edge that could not be split), and the volume gate's rejections of earlier candidate meshes, once a later step (the remesh retry or the stuffing fallback) has meshed the part and its mesh passed the volume gate: a finding left for later becomes `remesh` or `fallback` / `applied`, and a failed repair keeps its `failed` outcome. `REPAIR_GUARD_EXCEEDED` stays an error.

| Code | Default severity | Meaning | Default action |
| --- | --- | --- | --- |
| `PARSE_BAD_LINE` | warning | Malformed or non-numeric OBJ line, or a .meshb vertex that is not finite | Skip line, warning (`.meshb`: `vertex`, 1-based) |
| `PARSE_BAD_INDEX` | warning | Face references a missing vertex | Drop face, warning |
| `PARSE_READ_FAILED` | fatal | Input file cannot be opened or read | Fatal |
| `PARSE_PROJECT_INVALID` | fatal | Project file cannot be read or is not a valid project | Fatal |
| `PARSE_NO_GEOMETRY` | fatal | No faces survive parsing | Fatal |
| `PARSE_FORMAT_INVALID` | error | A binary input file (.meshb) is malformed or truncated | Abort reading at `offset` (bytes) with a `reason` (byte-order code, version, truncation, a bad next position); what was read before is kept |
| `PARSE_UNSUPPORTED_ELEMENT` | info | Lines, points, curves or surfaces in the OBJ file, or elements other than triangles in a .meshb file | Skip line, info (`.meshb`: one per keyword, `statement` its name and `count`) |
| `PARSE_POLYGON_FALLBACK` | warning | Polygon could not be ear-clipped and was fan-triangulated | Fan-triangulate, warning |
| `PART_EMPTY_SELECTION` | fatal | No parts selected, or a named part does not exist | Fatal |
| `PART_OVERLAP` | error | Two selected parts intersect | Error |
| `PART_SHARED_INTERFACE` | info | Two selected parts touch along a surface | Info; meshed independently in v1 |
| `GEO_DUPLICATE_VERTEX` | warning | Vertices closer than tolerance | Merge, warning |
| `GEO_DEGENERATE_FACE` | warning | Zero-area or collinear face | Remove, warning |
| `GEO_DUPLICATE_FACE` | warning | Same face listed twice | Remove, warning |
| `TOP_BOUNDARY_EDGE` | warning | Hole: edge with one adjacent face | Fill if small, else error, triggers fallback (stuffing, then the wrap of repair rung 3) |
| `TOP_NONMANIFOLD_EDGE` | warning | Edge shared by three or more faces | Split, else error, triggers fallback (stuffing, then the wrap of repair rung 3) |
| `TOP_NONMANIFOLD_VERTEX` | warning | Vertex joining separate fans | Duplicate vertex, warning |
| `TOP_INCONSISTENT_ORIENTATION` | warning | Neighbouring faces disagree on winding | Reorient, warning |
| `TOP_INVERTED_SHELL` | warning | Shell with negative signed volume | Flip, warning |
| `TOP_OPEN_SHEET` | fatal | Surface encloses no volume | Fatal for that part |
| `INT_SELF_INTERSECTION` | error | Faces of one shell intersect | Resolve (repair rung 2) for a closed part, warning; else error, triggers fallback (stuffing, then the wrap of repair rung 3) |
| `INT_SHELL_OVERLAP` | error | Two shells of one part intersect | Resolve (repair rung 2) for a closed part, warning; else error, triggers fallback (stuffing, then the wrap of repair rung 3) |
| `QUAL_SLIVER_FACE` | warning | Minimum angle below threshold | Remesh, warning |
| `QUAL_SUSPECT_SCALE` | info | Bounding box suggests wrong units | Info; tolerances are relative |
| `FEAT_OVERRIDE_IGNORED` | warning | A feature or smooth edge override took no effect | Warning: the override is ignored |
| `SIZE_CONTROL_IGNORED` | warning | A size control took no effect | Warning: the control is ignored |
| `MESH_RECOVERY_FAILED` | error | Constrained Delaunay could not recover a facet | Error, triggers fallback (stuffing, then the wrap of repair rung 3) |
| `MESH_BUDGET_EXCEEDED` | error | Steiner-point or element budget hit | Error, triggers fallback |
| `GEO_CAPACITY_EXCEEDED` | fatal | Part too large for 32-bit vertex and face indices | Fatal |
| `REPAIR_GUARD_EXCEEDED` | error | A repair changed the shell count or enclosed volume beyond its limit | Error |
| `MESH_STRICT_UNMET` | error | --strict was given but no conformal mesh could be produced | Error |
| `MESH_WORKER_FAILED` | error | Worker panicked, timed out or ran out of memory | Error, triggers fallback |
| `MESH_ROUNDING_REPAIRED` | warning | A tet that rounding to f64 turned over or flattened was untangled | Untangle, warning |
| `MESH_ROUNDING_FAILED` | error | Tets that rounding to f64 turned over or flattened could not be untangled | Error, triggers fallback |
| `VOL_INVERTED_TET` | error | Tetrahedron with non-positive volume | Error |
| `VOL_NONMANIFOLD` | error | Mesh face held by more than two tets, or by two on the same side | Error |
| `VOL_LOW_QUALITY` | warning | Non-orthogonality, skewness or volume ratio out of bounds | Warning or error by threshold |
| `VOL_BOUNDARY_DEVIATION` | error | Hausdorff distance to input over tolerance | Error |
| `VOL_GRADATION` | warning | Neighbouring cells' sizes differ beyond the gradation limit | Warning, located at the face with the largest ratio of its cells' longest edges (`max_size_ratio` against 4; `element`, `neighbour`, `point`). Conformal meshes are graded to the limit before the gate, so it is rare: where a cell could not be split, or grading ran out of budget |
| `EXO_NAME_TRUNCATED` | warning | Part or set name over the Exodus limit | Truncate, warning |
| `EXO_CAPACITY_EXCEEDED` | error | Mesh too large for the 64-bit-offset format | Error |
| `OUT_NAME_CHANGED` | warning | A part or patch name was changed to one the output format accepts | Rename, warning: `format` (`openfoam` or `fluent`), `entity` (`patch` or `cell_zone`), `name`, `written_as`. OpenFOAM names here are ASCII letters, digits, `_`, `.` and `-`, starting with a letter or `_`; Fluent zone names lowercase letters, digits, `_`, `.`, `+` and `-`, starting with a letter or `_` (capitals are lowered, and the interior zone is `interior`). Other characters become `_`, and a repeated name takes `_2`, `_3`, ... |
| `OUT_WRITE_FAILED` | fatal | An output file could not be written | Fatal |

`VOL_LOW_QUALITY` is raised, a warning, for each of a mesh's dihedral angles outside its kernel's bounds, its worst non-orthogonality over 70° and its worst skewness over 2 (the targets the cell improvement aims at). The volume ratio is reported in the mesh summary but not gated.

## Attempts

Each meshing attempt, in the order tried: `part`, `kernel` (`delaunay`, the constrained Delaunay kernel, or `stuffing`, isosurface stuffing), `outcome` (`succeeded` or `failed`), `remeshed` (present when the attempt meshed the remeshed surface), and for a failed worker run its `stderr`. Within one worker run a part goes through at most three: the conformal mesh of the repaired surface, the conformal mesh of the remeshed surface, and stuffing. When every conformal mesh is far beyond the quality limits, stuffing is tried before they are gated: kept, the conformal attempts are listed as `failed` (meshed, not kept); not kept, the stuffing attempt is listed as `failed` and the conformal ones follow. When the client retried a timed-out worker with stuffing only, the failed run comes first, as a `delaunay` attempt with no part.

## Mesh quality

`parts[].mesh` describes the final attempt's mesh of the part:

| Field | Meaning |
| --- | --- |
| `kernel`, `accepted` | The kernel, and whether the mesh passed the volume gate (and was written) |
| `tets`, `nodes`, `steiner_points` | Counts; Steiner points are those the constrained kernel added |
| `element_size` | The largest target size (absolute) |
| `min_element_size` | The smallest target size, when it was graded: the sizing field's for a conformal mesh, the surface cell size for graded stuffing |
| `remeshed` | Present when the mesh conforms to the remeshed surface rather than the repaired one |
| `volume`, `surface_volume` | The mesh's volume, and the volume the repaired surface encloses |
| `min_dihedral_degrees`, `max_dihedral_degrees` | Extreme dihedral angles |
| `dihedral_histogram` | Every dihedral angle (six per tet) in 10° bins, [0°, 10°) to [170°, 180°] |
| `max_non_orthogonality_degrees`, `mean_non_orthogonality_degrees` | On interior faces, the angle between the face normal and the line joining the two cell centroids; what a finite-volume solver's non-orthogonal correction has to make up |
| `max_skewness` | On interior faces, the distance from the face centroid to where the line joining the two cell centroids crosses the face, over the distance between the cell centroids |
| `max_volume_ratio` | Largest volume ratio between neighbouring cells, over pairs whose floating-point volumes are both positive |
| `max_size_ratio` | Largest ratio between neighbouring cells' sizes (each cell's longest edge): how abruptly the size changes. Above 4 it is `VOL_GRADATION` (stuffing's graded lattice stays under 2.9) |
| `hausdorff`, `hausdorff_relative` | Two-sided distance between the input surface and the mesh boundary (sampled at vertices and face centroids), absolute and relative to the part's diagonal; zero to rounding for a conformal mesh of the repaired surface |
| `worst_cells` | Where the worst cells are, one entry per measure: `measure` (`min_dihedral`, `non_orthogonality`, `skewness` or `size_ratio`), `wall` (present for the worst over the wall-adjacent cells, those with a boundary face; their boundary faces count for non-orthogonality, by the line from the cell's centroid to the face's), `value`, `element` (the cell's 1-based number in the part's element block; for a face, the first of its two cells), `neighbour` (for an interior face, the cell across it) and `point` (the cell's centroid, or the face's). Entries in order: smallest dihedral angle, then non-orthogonality and skewness over every interior face, then both over the wall-adjacent cells. Ties go to the lowest element number |

Every number is finite: JSON has no infinity or NaN. A value that could not be measured (the deviation of an empty mesh, for one) is written as -1.

## Summary and exit code

`summary` counts records by severity (`info`, `warnings`, `errors`, `fatal`) and by code (`by_code`: `total` occurrences and `reported` records), and gives `exit_code`, which `mesher` also exits with:

| Code | Meaning | Rule |
| --- | --- | --- |
| 0 | Clean | No warnings, errors or fatal records |
| 10 | Succeeded after repairs | At least one warning, no errors |
| 20 | Succeeded via the stuffing fallback | A part was meshed by stuffing, no errors |
| 30 | Failed | At least one error or fatal record; the report says why |
| 70 | Internal bug | The worker panicked, could not start, or broke its protocol |

## Other outputs

Beside `report.json`, a run writes:

- `mesh.exo`: Exodus II (netCDF-3 64-bit offset). One `TETRA4` element block per meshed part, named after it; a side set `<part>_boundary` per block, and one per OBJ group (on every kernel: a stuffed mesh's boundary faces take the group of a nearby input triangle facing their way, so near a group boundary they may go to either side; fill faces are in `<part>_boundary` only; a group in several parts is one side set spanning their blocks), using Exodus's tetrahedron face numbering; QA records with the tool and version. Names longer than 32 characters are truncated (`EXO_NAME_TRUNCATED`).
- `mesh.msh` (Gmsh 4) and `mesh.vtu` (VTK), the same mesh.
- `openfoam/`, with `--formats openfoam`: an OpenFOAM case, ready for `checkMesh`. `constant/polyMesh/` holds the mesh in ASCII (`points`, `faces`, `owner`, `neighbour`, `boundary`, `cellZones`): a cell zone per part, and a patch (type `patch`) per OBJ group, each boundary face in exactly one (a face in a group's side set goes to that group's patch, the rest of a part's boundary to `<part>_boundary`; empty patches are left out). `system/` holds a minimal `controlDict`, `fvSchemes` and `fvSolution` for a solver's case to replace. Names that are not OpenFOAM words are changed (`OUT_NAME_CHANGED`).
- `mesh_fluent.msh`, with `--formats fluent`: a Fluent mesh file (ASCII): a cell zone of tetrahedra per part (type `fluid`), the internal faces as the zone `interior`, and a `wall` face zone per patch, the same patches as the OpenFOAM case. Names that are not lowercase Fluent symbols are changed (`OUT_NAME_CHANGED`).
- `mesh.meshb`, with `--formats meshb`: a libMeshb binary file (version 2, or 3 beyond 2 GiB): the points, the boundary triangles out of the mesh with each patch's 1-based number as reference, and the tetrahedra with each part's. `mesh_refs.txt` names each reference (`tetrahedra <ref> <part>`, `triangles <ref> <patch>`), as the format holds no names.
- `repaired.obj`: the repaired surface.
- `defects.vtp`: the findings' faces, edges and vertices, for ParaView.
- `features.vtp`: the feature lines of every repaired part (as `features` counts them), one line cell per feature edge, with cell arrays `part` (1-based), `cause` (0 the angle, 1 a group boundary, 2 forced) and `curve` (within the part); written whenever a part was repaired.
- `project.toml`, when the GUI ran the job: the input, output directory and settings, for `mesher run project.toml`.
