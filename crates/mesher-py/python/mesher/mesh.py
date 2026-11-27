"""A run's mesh, read back from its ``.meshb`` file.

The worker writes ``mesh.meshb`` (with ``mesh_refs.txt`` naming its
references) when ``"meshb"`` is among the run's formats. ``read_mesh``
reads it: points, tetrahedra (0-based) with each one's part, and boundary
triangles with each one's patch, as numpy arrays when numpy is installed,
else tuples.
"""

from __future__ import annotations

import os
import pathlib
from dataclasses import dataclass
from typing import Any, Optional, Tuple, Union

from ._mesher import _read_mesh

__all__ = ["Mesh", "read_mesh"]


@dataclass(frozen=True, eq=False)
class Mesh:
    """A volume mesh: ``points`` (n × 3 floats), ``tets`` (m × 4 indices into
    ``points``, 0-based, as the worker wrote them), ``tet_parts`` (per
    tetrahedron, its index in ``parts``), ``triangles`` (k × 3: each
    boundary face once, ordered so its normal points out of its cell) and
    ``triangle_patches`` (per triangle, its index in ``patches``).

    Arrays are read-only numpy arrays (float64, int64) when the mesh was
    read with numpy, else tuples of tuples (and tuples of ints for the
    per-element indices). ``parts`` and ``patches`` are names, in order of
    the file's references: parts are the report's parts that meshed;
    patches are the boundary's side sets that hold faces (each OBJ group,
    and ``<part>_boundary`` for the faces in none). Faces between two parts
    are internal, so in no patch.
    """

    points: Any
    tets: Any
    tet_parts: Any
    parts: Tuple[str, ...]
    triangles: Any
    triangle_patches: Any
    patches: Tuple[str, ...]
    #: The ``.meshb`` file read.
    path: pathlib.Path

    def part_tets(self, name: str) -> Any:
        """The tetrahedra of the part named ``name``; ``KeyError`` if the mesh
        has no such part."""
        return _select(self.tets, self.tet_parts, _index(self.parts, name))

    def patch_triangles(self, name: str) -> Any:
        """The boundary triangles of the patch named ``name``; ``KeyError`` if
        the mesh has no such patch."""
        return _select(self.triangles, self.triangle_patches, _index(self.patches, name))


def _index(names: Tuple[str, ...], name: str) -> int:
    try:
        return names.index(name)
    except ValueError:
        raise KeyError(name) from None


def _select(rows: Any, groups: Any, i: int) -> Any:
    if isinstance(rows, tuple):
        return tuple(r for r, g in zip(rows, groups) if g == i)
    return rows[groups == i]


def _numpy(wanted: Optional[bool]) -> Any:
    """numpy, if it is to be used: ``None`` uses it when it is installed."""
    if wanted is False:
        return None
    try:
        import numpy
    except ImportError:
        if wanted:
            raise
        return None
    return numpy


def read_mesh(path: Union[str, "os.PathLike[str]"], *, numpy: Optional[bool] = None) -> Mesh:
    """The mesh at ``path``: a run's output directory (its ``mesh.meshb``) or
    a ``.meshb`` file (named by ``<stem>_refs.txt`` beside it; references
    it does not name are ``ref<k>``).

    ``numpy``: ``None`` (default) gives numpy arrays when numpy is
    installed, else tuples; ``True`` requires numpy; ``False`` gives tuples.
    Raises ``FileNotFoundError`` (or another ``OSError``) when the file
    cannot be read and ``ValueError`` when it is not a complete mesh.
    """
    path = pathlib.Path(path)
    raw = _read_mesh(os.fspath(path))
    np = _numpy(numpy)
    if np is not None:

        def rows(data: bytes, dtype: Any, k: int) -> Any:
            return np.frombuffer(data, dtype=dtype).reshape(-1, k)

        def flat(data: bytes) -> Any:
            return np.frombuffer(data, dtype=np.int64)

        points = rows(raw["points"], np.float64, 3)
        tets = rows(raw["tets"], np.int64, 4)
        triangles = rows(raw["triangles"], np.int64, 3)
        tet_parts = flat(raw["tet_parts"])
        triangle_patches = flat(raw["triangle_patches"])
    else:

        def rows(data: bytes, code: str, k: int) -> Any:  # type: ignore[misc]
            values = memoryview(data).cast(code)
            return tuple(zip(*[iter(values)] * k))

        def flat(data: bytes) -> Any:  # type: ignore[misc]
            return tuple(memoryview(data).cast("q"))

        points = rows(raw["points"], "d", 3)
        tets = rows(raw["tets"], "q", 4)
        triangles = rows(raw["triangles"], "q", 3)
        tet_parts = flat(raw["tet_parts"])
        triangle_patches = flat(raw["triangle_patches"])
    file = path / "mesh.meshb" if path.is_dir() else path
    return Mesh(
        points=points,
        tets=tets,
        tet_parts=tet_parts,
        parts=tuple(raw["parts"]),
        triangles=triangles,
        triangle_patches=triangle_patches,
        patches=tuple(raw["patches"]),
        path=file,
    )
