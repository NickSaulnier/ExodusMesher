"""A run's report, typed.

The dataclasses and enums are generated from the report schema
(``_schema``); ``Report`` adds what the schema cannot say: the worker's own
text, views over the records, the run's mesh, and the job to run next.
"""

from __future__ import annotations

import json
import os
import pathlib
from typing import Any, Dict, NamedTuple, Optional, Tuple, Union

from ._mesher import _next_job
from .mesh import Mesh, read_mesh

from . import _schema
from ._schema import (
    Action,
    Attempt,
    AttemptOutcome,
    CodeCount,
    DefectCode,
    FeatureSummary,
    InputInfo,
    InputStats,
    JobConfig,
    Kernel,
    KernelChoice,
    Location,
    MeshSummary,
    Outcome,
    OutputFormat,
    PartSource,
    PartSummary,
    QualityMeasure,
    Record,
    Remedy,
    RepairSummary,
    Setting,
    Severity,
    SizeControl,
    Stage,
    Summary,
    ToolInfo,
    WorstCell,
)

__all__ = list(_schema.__all__) + ["NextJob", "read_report"]


class NextJob(NamedTuple):
    """The job to run next after a report (``Report.next_job``): ``job`` for
    ``mesher.run_job``, made by applying ``remedy`` of ``record``."""

    job: Dict[str, Any]
    record: Record
    remedy: Remedy


class Report(_schema.Report):
    """A run's report: ``report.json``, typed.

    Build one with ``Report.from_json`` (or ``read_report``); ``mesher.run``
    and ``mesher.run_project`` return one. ``to_json()`` gives back the
    worker's text exactly, byte for byte.
    """

    @classmethod
    def from_json(cls, text: str, path: Optional[Union[str, os.PathLike]] = None) -> Report:
        """The report in ``text`` (the worker's JSON); ``path`` is where it
        was written, if anywhere. Raises ``KeyError`` or ``ValueError`` for
        JSON that is not a report this version reads."""
        report = cls._parse(json.loads(text))
        # The dataclass is frozen; the text and path are not fields, so they
        # do not take part in equality.
        object.__setattr__(report, "_text", text)
        object.__setattr__(report, "_path", None if path is None else pathlib.Path(path))
        return report

    def to_json(self) -> str:
        """The worker's JSON text, exactly as ``report.json`` holds it."""
        text = getattr(self, "_text", None)
        if text is None:
            raise ValueError("this report was not read from the worker's JSON")
        return text

    @property
    def path(self) -> Optional[pathlib.Path]:
        """Where ``report.json`` was written (``None`` if it could not be, or
        the report was not read from a run)."""
        return getattr(self, "_path", None)

    @property
    def exit_code(self) -> int:
        """The exit code ``mesher run`` gives for this report: 0 clean, 10
        succeeded after repairs (warnings), 20 succeeded by the stuffing
        fallback, 30 failed (an error or a fatal record), 70 an internal bug
        (see the report reference)."""
        return self.summary.exit_code

    @property
    def records_by_code(self) -> Dict[DefectCode, Tuple[Record, ...]]:
        """The records grouped by code, each group in report order, codes in
        order of first appearance."""
        groups: Dict[DefectCode, list] = {}
        for record in self.records:
            groups.setdefault(record.code, []).append(record)
        return {code: tuple(records) for code, records in groups.items()}

    def part(self, name: str) -> PartSummary:
        """The part named ``name``; ``KeyError`` if there is none."""
        for part in self.parts:
            if part.name == name:
                return part
        raise KeyError(name)

    def mesh(self, *, numpy: Optional[bool] = None) -> Mesh:
        """The run's mesh, read from the ``mesh.meshb`` beside ``report.json``
        (see ``read_mesh``). The run must have written it: ``"meshb"`` among
        its formats (``mesher.run(..., formats=["exodus", "meshb"])``, or the
        project's ``formats``). Raises ``ValueError`` when it did not, or the
        report was not read from a run, and ``FileNotFoundError`` when no
        mesh was written (no part meshed)."""
        if self.path is None:
            raise ValueError("this report was not read from a run's report.json")
        if OutputFormat.MESHB not in self.config.formats:
            raise ValueError(
                "the run wrote no .meshb: run it with \"meshb\" among its formats"
            )
        return read_mesh(self.path.parent, numpy=numpy)

    def next_job(
        self,
        out: Optional[Union[str, os.PathLike]] = None,
        *,
        severity: Union[Severity, str] = Severity.ERROR,
    ) -> Optional[NextJob]:
        """The job to run next: this report's input with the settings of the
        first remedy that applies and changes them, of the worst record at
        ``severity`` or above (fatal records first, then errors, in report
        order). ``None`` when there is none: the run needs nothing, or
        nothing known would help.

        The job writes to ``out`` (default: this run's output directory, so
        its files are replaced). Pass ``severity=Severity.WARNING`` to also
        act on failures a fallback covered (a part meshed by stuffing rather
        than conformally). Run it with ``mesher.run_job(next.job)``.
        Raises ``ValueError`` for a report not read from the worker's JSON,
        or without a path when ``out`` is not given."""
        if out is None:
            if self.path is None:
                raise ValueError("this report has no output directory: pass `out`")
            out = self.path.parent
        level = severity.value if isinstance(severity, Severity) else Severity(severity).value
        found = _next_job(self.to_json(), os.fspath(out), level)
        if found is None:
            return None
        text, record, remedy = found
        chosen = self.records[record]
        return NextJob(json.loads(text), chosen, chosen.remedies[remedy])


def read_report(path: Union[str, os.PathLike]) -> Report:
    """The report in a ``report.json`` file."""
    path = pathlib.Path(path)
    return Report.from_json(path.read_text(encoding="utf-8"), path)
