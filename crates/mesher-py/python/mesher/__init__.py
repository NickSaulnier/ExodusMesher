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

from ._mesher import __version__, _run

__all__ = ["RawRun", "run_job", "worker_path", "__version__"]


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
