"""ExodusMesher from Python.

    import mesher
    report = mesher.run("part.obj", element_size=0.02)
    print(report.exit_code, report.parts[0].mesh.tets)
    mesh = mesher.run("part.obj", formats=["exodus", "meshb"]).mesh()
    print(mesh.points.shape, mesh.parts, mesh.patches)

Jobs run in a separate ``mesher-worker`` process (shipped in this package,
or ``$MESHER_WORKER``), under a timeout and a memory cap, so a crash or a
runaway job never takes down the interpreter. Every run ends in a report,
the same ``report.json`` the command line writes: ``run`` and
``run_project`` return it typed (``Report``), ``run_job`` and
``run_request`` as text. ``Report.mesh()`` and ``read_mesh`` read the
mesh back from a run's ``.meshb``. Records carry ``remedies``, settings
changes to try next; ``Report.next_job()`` applies the first that fits.
"""

from __future__ import annotations

import json
import os
import pathlib
import sys
from enum import Enum
from typing import Any, Callable, Iterable, NamedTuple, Optional, Union

from . import mesh as _mesh
from . import report as _report
from ._mesher import __version__, _build_job, _run, _run_request
from .mesh import Mesh, read_mesh
from .report import *  # noqa: F401,F403 (the report's types)
from .report import CellChoice, KernelChoice, OutputFormat, Report

__all__ = [
    "RawRun",
    "build_job",
    "run",
    "run_job",
    "run_project",
    "run_request",
    "worker_path",
    "__version__",
] + _mesh.__all__ + _report.__all__

PathLike = Union[str, "os.PathLike[str]"]
Progress = Callable[[str, float], None]


class RawRun(NamedTuple):
    """A finished run: the report as JSON text, and the path of ``report.json``
    (``None`` if it could not be written)."""

    report_json: str
    report_path: Optional[str]

    def report(self) -> Report:
        """The report, typed."""
        return Report.from_json(self.report_json, self.report_path)


def worker_path() -> str:
    """The worker binary: ``$MESHER_WORKER``, else the one in this package."""
    env = os.environ.get("MESHER_WORKER")
    if env:
        return env
    exe = "mesher-worker.exe" if sys.platform == "win32" else "mesher-worker"
    return str(pathlib.Path(__file__).parent / "bin" / exe)


def run_job(
    job: dict[str, Any],
    *,
    timeout: float = 3600.0,
    memory_limit: Optional[int] = 8192,
    worker: Optional[str] = None,
    progress: Optional[Callable[[str, float], None]] = None,
) -> RawRun:
    """Run a job, given as the worker protocol's ``Job`` (``input``,
    ``output_dir``, ``config``), and wait for its report.

    ``timeout`` is in seconds and ``memory_limit`` in MiB (``None``: no cap),
    as the command line's ``--timeout`` and ``--memory-limit``.
    ``progress(stage, fraction)`` is called as the worker reports progress.
    Ctrl-C, or an exception from ``progress``, cancels the worker and is
    raised once its report is written.
    """
    memory = None if memory_limit is None else int(memory_limit) * (1 << 20)
    text, path = _run(json.dumps(job), worker or worker_path(), float(timeout), memory, progress)
    return RawRun(text, path)


def build_job(request: dict[str, Any]) -> Optional[dict[str, Any]]:
    """The job a request builds, exactly as ``mesher run`` builds it: the
    request names an ``input`` (an OBJ or ``.meshb`` file, or a project file)
    and the command line's options (``out``, ``parts``, ``element_size``,
    ``kernel``, ``formats``, ``cells``, ...). ``None`` for a project file that cannot be
    used. Raises ``ValueError`` for an unknown option or one out of range."""
    text = _build_job(json.dumps(_paths(request)))
    return None if text is None else json.loads(text)


def run_request(
    request: dict[str, Any],
    *,
    timeout: float = 3600.0,
    memory_limit: Optional[int] = 8192,
    worker: Optional[str] = None,
    progress: Optional[Callable[[str, float], None]] = None,
) -> RawRun:
    """Run a request as ``mesher run`` does (see ``build_job``): the same job,
    so the same report, as the command line given the same options; a
    project file that cannot be used gives a ``PARSE_PROJECT_INVALID``
    report. Raises ``ValueError`` before anything runs for an unknown option
    or one out of range. Limits and cancelling as for ``run_job``."""
    memory = None if memory_limit is None else int(memory_limit) * (1 << 20)
    text, path = _run_request(
        json.dumps(_paths(request)), worker or worker_path(), float(timeout), memory, progress
    )
    return RawRun(text, path)


def _paths(request: dict[str, Any]) -> dict[str, Any]:
    """Path-like values as strings, so ``pathlib.Path`` arguments work."""
    return {k: os.fspath(v) if isinstance(v, os.PathLike) else v for k, v in request.items()}


def run(
    input: PathLike,
    out: Optional[PathLike] = None,
    *,
    parts: Iterable[Union[str, int]] = (),
    element_size: Optional[float] = None,
    max_deviation: Optional[float] = None,
    kernel: Optional[Union[KernelChoice, str]] = None,
    formats: Iterable[Union[OutputFormat, str]] = (),
    cells: Optional[Union[CellChoice, str]] = None,
    strict: bool = False,
    no_repair: bool = False,
    no_remesh: bool = False,
    uniform_size: bool = False,
    wrap_alpha: Optional[float] = None,
    feature_angle: Optional[float] = None,
    timeout: float = 3600.0,
    memory_limit: Optional[int] = 8192,
    worker: Optional[str] = None,
    progress: Optional[Progress] = None,
) -> Report:
    """Mesh ``input`` (an OBJ or ``.meshb`` file, or a project file) as
    ``mesher run`` does with the same options, and return its report.

    ``out`` is the output directory (default: the project's, else
    ``<input stem>_out`` in the current directory). ``parts`` are names or
    1-based indices (default: all). ``element_size``, ``max_deviation`` and
    ``wrap_alpha`` are relative to each part's bounding-box diagonal (or,
    for the alpha, the element size), ``feature_angle`` in degrees; unset,
    each keeps the project's value or the default. ``formats`` are written
    sorted, each once (default: Exodus, Gmsh and VTU). ``cells`` is
    ``"tet"`` or ``"dual"`` (experimental: the polyhedral dual of each
    part's tet mesh, written to Exodus, VTU, OpenFOAM and Fluent; default:
    the project's, else tets). Raises
    ``ValueError`` before anything runs for an option out of range.

    ``timeout`` is in seconds and ``memory_limit`` in MiB (``None``: no cap).
    ``progress(stage, fraction)`` is called as the worker reports progress.
    Ctrl-C, or an exception from ``progress``, cancels the worker and is
    raised once its report is written.
    """
    request = {
        "input": input,
        "out": out,
        "parts": [str(p) for p in parts],
        "element_size": element_size,
        "max_deviation": max_deviation,
        "kernel": _value(kernel),
        "formats": [_value(f) for f in formats],
        "cells": _value(cells),
        "strict": strict,
        "no_repair": no_repair,
        "no_remesh": no_remesh,
        "uniform_size": uniform_size,
        "wrap_alpha": wrap_alpha,
        "feature_angle": feature_angle,
    }
    return run_request(
        request, timeout=timeout, memory_limit=memory_limit, worker=worker, progress=progress
    ).report()


def run_project(
    path: PathLike,
    *,
    timeout: float = 3600.0,
    memory_limit: Optional[int] = 8192,
    worker: Optional[str] = None,
    progress: Optional[Progress] = None,
) -> Report:
    """Run a project file (``.toml``: its input, parts, output directory
    and settings) and return its report; one that cannot be used gives a
    ``PARSE_PROJECT_INVALID`` report. To override its settings, pass it to
    ``run`` with the options. Raises ``ValueError`` for a path that is not a
    ``.toml`` file. Limits and cancelling as for ``run``."""
    if pathlib.Path(path).suffix.lower() != ".toml":
        raise ValueError(f"not a project file (.toml): {os.fspath(path)}")
    return run_request(
        {"input": path}, timeout=timeout, memory_limit=memory_limit, worker=worker, progress=progress
    ).report()


def _value(x: Any) -> Any:
    """An enum as its value, so ``KernelChoice.EXACT`` and ``"exact"`` both work."""
    return x.value if isinstance(x, Enum) else x
