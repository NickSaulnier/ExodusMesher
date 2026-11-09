# The report (`report.json`)

Every run writes `report.json` to its output directory, whatever happens: a clean mesh, a repaired one, a fallback, or a failure with no mesh at all. It records every finding, every repair, every meshing attempt and the quality of every mesh written, so a script can decide what to do without reading the logs. The layout is versioned (`schema_version`, currently 1), and its JSON schema is [`schema/report.schema.json`](../schema/report.schema.json), also printed by `mesher schema`.

The same input and settings give a byte-identical report.

## Top level

| Field | Contents |
| --- | --- |
| `schema_version` | Layout version; bumped on any breaking change |
| `tool` | `name` and `version` of the mesher |
| `input` | `path` as given, and `stats` (absent if the file could not be read): `bytes`, `lines`, `vertices`, `polygons` (faces as written), `triangles` (after triangulating polygons), `objects` and `groups` (names in order of appearance) |
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
| `no_remesh` | Present (true) when the remesh-and-retry step was turned off |
| `stuffing_only` | Present (true) when the client retried a worker that timed out or ran out of memory with the stuffing fallback only |

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
| `repair` | `faces_before`, `faces_after`, `vertices_before`, `vertices_after`, and `gate_passed`: whether the repaired surface is a valid closed solid |
| `mesh` | Present when a volume mesh was attempted; see [Mesh quality](#mesh-quality) |

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
| `action`, `outcome` | What was done (`skip_line`, `drop_face`, `fan_triangulate`, `remesh`, `merge`, `remove`, `fill_hole`, `split_edge`, `duplicate_vertex`, `reorient`, `flip`, `truncate`, `fallback`, `abort`) and whether it was `applied`, `failed` or `deferred` to a later stage |
| `hint` | A sentence on what it means or what to try |

At most 100 records are kept per code; `summary.by_code` still counts every occurrence. A finding is never repaired without a record.

The same records are in `defects.vtp` (open it in ParaView beside the mesh) and are drawn by the GUI, which shows each record's OBJ lines.

### Records from meshing

These carry the measurements most often needed to act on them:

| Code | `measured` (and `threshold`) |
| --- | --- |
| `MESH_RECOVERY_FAILED` | `missing_faces`; after a remesh retry, `remeshed_triangles` (it succeeded: the record is then a warning with `remesh` / `applied`) or `remesh_retry` (why the retry did not help) |
| `MESH_BUDGET_EXCEEDED` | `points` against `max_points`, `work_steps` against `max_work_steps` (the constrained kernel), or `lattice_vertices` against `max_lattice_vertices` (stuffing) |
| `MESH_WORKER_FAILED` | `reason`: `timeout` (with `timeout_seconds`), `out_of_memory` (with `memory_limit_bytes`), `cancelled`, `spawn_failed` (with `io_error`) or `died` (with `exit_status`); `stderr_tail`. After a timeout or out-of-memory, the client retries with stuffing only and the record becomes a warning with `fallback` / `applied` if that produced a mesh |
| `VOL_LOW_QUALITY` | `min_dihedral_degrees` and `max_dihedral_degrees` against the kernel's bounds (10° and 170° for the constrained Delaunay kernel; 10.7° and 164.8°, the proven bounds, for stuffing), and `max_non_orthogonality_degrees` |
| `VOL_BOUNDARY_DEVIATION` | `hausdorff` and `hausdorff_relative` against `max_deviation_relative` |

## Defect codes

Severity is the code's default; records can be upgraded or downgraded where the outcome warrants it (a conformal failure that a fallback recovered is a warning).

| Code | Default severity | Meaning | Default action |
| --- | --- | --- | --- |
| `PARSE_BAD_LINE` | warning | Malformed or non-numeric OBJ line | Skip line, warning |
| `PARSE_BAD_INDEX` | warning | Face references a missing vertex | Drop face, warning |
| `PARSE_READ_FAILED` | fatal | Input file cannot be opened or read | Fatal |
| `PARSE_PROJECT_INVALID` | fatal | Project file cannot be read or is not a valid project | Fatal |
| `PARSE_NO_GEOMETRY` | fatal | No faces survive parsing | Fatal |
| `PARSE_UNSUPPORTED_ELEMENT` | info | Lines, points, curves or surfaces in the OBJ file | Skip line, info |
| `PARSE_POLYGON_FALLBACK` | warning | Polygon could not be ear-clipped and was fan-triangulated | Fan-triangulate, warning |
| `PART_EMPTY_SELECTION` | fatal | No parts selected, or a named part does not exist | Fatal |
| `PART_OVERLAP` | error | Two selected parts intersect | Error |
| `PART_SHARED_INTERFACE` | info | Two selected parts touch along a surface | Info; meshed independently in v1 |
| `GEO_DUPLICATE_VERTEX` | warning | Vertices closer than tolerance | Merge, warning |
| `GEO_DEGENERATE_FACE` | warning | Zero-area or collinear face | Remove, warning |
| `GEO_DUPLICATE_FACE` | warning | Same face listed twice | Remove, warning |
| `TOP_BOUNDARY_EDGE` | warning | Hole: edge with one adjacent face | Fill if small, else error |
| `TOP_NONMANIFOLD_EDGE` | warning | Edge shared by three or more faces | Split, else error |
| `TOP_NONMANIFOLD_VERTEX` | warning | Vertex joining separate fans | Duplicate vertex, warning |
| `TOP_INCONSISTENT_ORIENTATION` | warning | Neighbouring faces disagree on winding | Reorient, warning |
| `TOP_INVERTED_SHELL` | warning | Shell with negative signed volume | Flip, warning |
| `TOP_OPEN_SHEET` | fatal | Surface encloses no volume | Fatal for that part |
| `INT_SELF_INTERSECTION` | error | Faces of one shell intersect | Error, triggers fallback |
| `INT_SHELL_OVERLAP` | error | Two shells of one part intersect | Error, triggers fallback |
| `QUAL_SLIVER_FACE` | warning | Minimum angle below threshold | Remesh, warning |
| `QUAL_SUSPECT_SCALE` | info | Bounding box suggests wrong units | Info; tolerances are relative |
| `MESH_RECOVERY_FAILED` | error | Constrained Delaunay could not recover a facet | Error, triggers fallback |
| `MESH_BUDGET_EXCEEDED` | error | Steiner-point or element budget hit | Error, triggers fallback |
| `GEO_CAPACITY_EXCEEDED` | fatal | Part too large for 32-bit vertex and face indices | Fatal |
| `REPAIR_GUARD_EXCEEDED` | error | A repair changed the shell count or enclosed volume beyond its limit | Error |
| `MESH_STRICT_UNMET` | error | --strict was given but no conformal mesh could be produced | Error |
| `MESH_WORKER_FAILED` | error | Worker panicked, timed out or ran out of memory | Error, triggers fallback |
| `VOL_INVERTED_TET` | error | Tetrahedron with non-positive volume | Error |
| `VOL_LOW_QUALITY` | warning | Non-orthogonality, skewness or volume ratio out of bounds | Warning or error by threshold |
| `VOL_BOUNDARY_DEVIATION` | error | Hausdorff distance to input over tolerance | Error |
| `EXO_NAME_TRUNCATED` | warning | Part or set name over the Exodus limit | Truncate, warning |
| `EXO_CAPACITY_EXCEEDED` | error | Mesh too large for the 64-bit-offset format | Error |
| `OUT_WRITE_FAILED` | fatal | An output file could not be written | Fatal |

`VOL_LOW_QUALITY` is raised when a mesh's dihedral angles fall outside its kernel's bounds; the other finite-volume metrics are reported in the mesh summary but not gated.

## Attempts

Each meshing attempt, in the order tried: `part`, `kernel` (`delaunay`, the constrained Delaunay kernel, or `stuffing`, isosurface stuffing), `outcome` (`succeeded` or `failed`), `remeshed` (present when the attempt meshed the remeshed surface), and for a failed worker run its `stderr`. Within one worker run a part goes through at most three: the conformal mesh of the repaired surface, the conformal mesh of the remeshed surface, and stuffing. When the client retried a timed-out worker with stuffing only, the failed run comes first, as a `delaunay` attempt with no part.

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
| `max_volume_ratio` | Largest volume ratio between neighbouring cells |
| `hausdorff`, `hausdorff_relative` | Two-sided distance between the input surface and the mesh boundary (sampled at vertices and face centroids), absolute and relative to the part's diagonal; zero to rounding for a conformal mesh of the repaired surface |

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

- `mesh.exo`: Exodus II (netCDF-3 64-bit offset). One `TETRA4` element block per meshed part, named after it; a side set `<part>_boundary` per block, and one per OBJ group on a conformal mesh, using Exodus's tetrahedron face numbering; QA records with the tool and version. Names longer than 32 characters are truncated (`EXO_NAME_TRUNCATED`).
- `mesh.msh` (Gmsh 4) and `mesh.vtu` (VTK), the same mesh.
- `repaired.obj`: the repaired surface.
- `defects.vtp`: the findings' faces, edges and vertices, for ParaView.
- `project.toml`, when the GUI ran the job: the input, output directory and settings, for `mesher run project.toml`.
