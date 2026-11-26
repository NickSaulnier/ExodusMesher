"""ExodusMesher from Python.

Jobs run in a separate ``mesher-worker`` process (shipped in this package,
or ``$MESHER_WORKER``), under a timeout and a memory cap, so a crash or a
runaway job never takes down the interpreter. Every run ends in a report,
the same ``report.json`` the command line writes.
"""

from __future__ import annotations

import json
import os
import pathlib
import sys
from typing import Any, Callable, NamedTuple, Optional

from ._mesher import __version__, _build_job, _run, _run_request

__all__ = ["RawRun", "build_job", "run_job", "run_request", "worker_path", "__version__"]


class RawRun(NamedTuple):
    """A finished run: the report as JSON text, and the path of ``report.json``
    (``None`` if it could not be written)."""

    report_json: str
    report_path: Optional[str]


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
    ``kernel``, ``formats``, ...). ``None`` for a project file that cannot be
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
