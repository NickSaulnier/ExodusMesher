"""Round trip of our .meshb files through meshio (MIT), an independent reader
and writer of libMeshb.

For each case directory under the given root (as the synthetic suite writes
them with SUITE_FOAM_DIR=<root>) holding `mesh.meshb` beside its OpenFOAM
case:

1. meshio reads `mesh.meshb`; its points must equal the case's
   `constant/polyMesh/points` exactly (the same points in the same order), its
   tetrahedra count the case's cells, per reference each cell zone's in order,
   and its triangles the boundary faces, per reference each patch's in order.
2. meshio writes what it read as `meshio.meshb`, which the Rust test
   `tests/synthetic/tests/meshb_meshio.rs` reads with our reader and compares
   with `mesh.meshb` (run it with MESHB_MESHIO_DIR=<root>).

Prints one line per case and exits 1 if any case differs. Needs numpy and
meshio (pinned in CI: meshio==5.3.5). Run as
`python -I tests/meshb/meshio_roundtrip.py <root>`.
"""

import collections
import pathlib
import re
import sys

import meshio
import numpy as np


def body(path):
    """A FoamFile's data after its header."""
    text = path.read_text()
    return text.split("}\n\n", 1)[1]


def foam_points(case):
    lines = body(case / "constant/polyMesh/points").split("\n")
    n = int(lines[0])
    rows = [line.strip("()").split() for line in lines[2 : 2 + n]]
    return np.array(rows, dtype=float).reshape(-1, 3)


def foam_counts(case, name, key):
    """The `key` count of each entry of a FoamFile list of dictionaries."""
    return [int(v) for v in re.findall(rf"{key}\s+(\d+);", body(case / name))]


def zone_sizes(case):
    """Each cell zone's size, from its `cellLabels List<label> N`."""
    return [int(v) for v in re.findall(r"cellLabels List<label> (\d+)", body(case / "constant/polyMesh/cellZones"))]


def by_ref(mesh, cell_type):
    counts = collections.Counter()
    for block, refs in zip(mesh.cells, mesh.cell_data.get("medit:ref", [])):
        if block.type == cell_type:
            counts.update(int(r) for r in refs)
    return [counts[r] for r in range(1, max(counts, default=0) + 1)]


def check(case):
    mesh = meshio.read(case / "mesh.meshb")
    problems = []
    points = foam_points(case)
    if mesh.points.shape != points.shape or not np.array_equal(mesh.points, points):
        problems.append(f"points differ ({mesh.points.shape} against {points.shape})")
    zones = zone_sizes(case)
    tets = by_ref(mesh, "tetra")
    if tets != zones[: len(tets)] or sum(tets) != sum(zones):
        problems.append(f"tetrahedra per reference {tets}, cell zones {zones}")
    patches = foam_counts(case, "constant/polyMesh/boundary", "nFaces")
    triangles = by_ref(mesh, "triangle")
    if triangles != patches:
        problems.append(f"triangles per reference {triangles}, patches {patches}")
    meshio.write(case / "meshio.meshb", mesh)
    return problems, len(mesh.points), sum(tets), sum(triangles)


def main():
    root = pathlib.Path(sys.argv[1])
    cases = sorted(p.parent for p in root.glob("*/mesh.meshb"))
    failed = 0
    for case in cases:
        try:
            problems, points, tets, triangles = check(case)
        except Exception as e:  # any reader failure is a failed case
            problems, points, tets, triangles = [f"{type(e).__name__}: {e}"], 0, 0, 0
        result = "FAILED" if problems else "ok"
        failed += bool(problems)
        print(f"{result}\t{case.name}\tpoints {points}\ttetrahedra {tets}\ttriangles {triangles}\t{'; '.join(problems)}")
    print(f"meshio: {len(cases) - failed} of {len(cases)} .meshb files read back")
    sys.exit(1 if failed or not cases else 0)


if __name__ == "__main__":
    main()
