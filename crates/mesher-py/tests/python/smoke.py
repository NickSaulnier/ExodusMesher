"""The package's check, run against an installed wheel (no test framework
needed):

    python -I crates/mesher-py/tests/python/smoke.py <output dir>

Meshes the cube fixture through the packaged worker, checks the report and
`report.json` agree, that progress arrives, and that a bad job raises before
anything runs (stage 1); that requests build the CLI's jobs (stage 2); and
`mesher.run`, `run_project` and the typed report (stage 3).
"""

import json
import pathlib
import sys

import mesher

ROOT = pathlib.Path(__file__).resolve().parents[4]


def main():
    out = pathlib.Path(sys.argv[1])
    cube = ROOT / "tests" / "fixtures" / "obj" / "cube.obj"
    assert pathlib.Path(mesher.worker_path()).is_file(), mesher.worker_path()

    stages = []
    run = mesher.run_job(
        {"input": str(cube), "output_dir": str(out), "config": {}},
        timeout=600,
        progress=lambda stage, fraction: stages.append((stage, fraction)),
    )
    report = json.loads(run.report_json)
    assert report["summary"]["exit_code"] == 0, report["summary"]
    assert run.report_path is not None
    assert pathlib.Path(run.report_path).read_text() == run.report_json
    tets = report["parts"][0]["mesh"]["tets"]
    assert tets > 0
    assert (out / "mesh.exo").is_file()
    assert stages and stages[0][0] == "parse", stages[:3]

    # An exception from the progress callback cancels the worker; the report
    # is still written, and the exception raised after it.
    class Stop(Exception):
        pass

    def stop(stage, fraction):
        raise Stop(stage)

    cancelled = out / "cancelled"
    try:
        mesher.run_job({"input": str(cube), "output_dir": str(cancelled), "config": {}}, progress=stop)
    except Stop:
        pass
    else:
        raise AssertionError("the callback's exception was not raised")
    after = json.loads((cancelled / "report.json").read_text())
    reasons = [r["measured"].get("reason") for r in after["records"] if r["code"] == "MESH_WORKER_FAILED"]
    assert reasons == ["cancelled"], reasons

    # Requests, as `mesher run` takes them.
    job = mesher.build_job({"input": cube, "element_size": 0.05, "formats": ["vtu", "exodus"]})
    assert job["output_dir"] == "cube_out", job
    assert job["config"]["element_size"] == 0.05
    assert job["config"]["formats"] == ["exodus", "vtu"]
    for bad in [{"input": cube, "element_size": 0}, {"input": cube, "colour": 1}, {"out": "x"}]:
        try:
            mesher.build_job(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{bad} built a job")
    broken = out / "broken.toml"
    broken.write_text("not a project")
    assert mesher.build_job({"input": broken}) is None
    failed = mesher.run_request({"input": broken, "out": out / "broken_out"})
    codes = [r["code"] for r in json.loads(failed.report_json)["records"]]
    assert codes == ["PARSE_PROJECT_INVALID"], codes
    assert pathlib.Path(failed.report_path).read_text() == failed.report_json
    by_request = mesher.run_request({"input": cube, "out": out / "by_request"}, timeout=600)
    assert by_request.report_json == run.report_json, "a request and its job differ"

    for bad in [{"input": 3}, {"output_dir": "x"}]:
        try:
            mesher.run_job(bad)
        except ValueError as e:
            assert "not a valid job" in str(e), e
        else:
            raise AssertionError(f"{bad} ran")

    typed(out, cube, run)

    print(f"mesher {mesher.__version__}: cube meshed, {tets} tets, {len(stages)} progress messages")


def typed(out, cube, run):
    """Stage 3: `mesher.run` and `run_project`, and the typed report."""
    report = mesher.run(cube, out / "typed", timeout=600)
    assert isinstance(report, mesher.Report)
    assert report.to_json() == run.report_json, "run and run_job differ"
    assert report.path == out / "typed" / "report.json", report.path
    assert report.path.read_text() == report.to_json()
    assert report.exit_code == 0 and report.summary.exit_code == 0
    assert report.schema_version >= 1
    part = report.parts[0]
    assert report.part(part.name) is part
    mesh = part.mesh
    assert isinstance(mesh, mesher.MeshSummary) and mesh.accepted
    assert mesh.kernel is mesher.Kernel.DELAUNAY, mesh.kernel
    assert mesh.tets == json.loads(run.report_json)["parts"][0]["mesh"]["tets"]
    assert isinstance(part.bounds, tuple) and len(part.bounds) == 2 and len(part.bounds[0]) == 3
    assert all(isinstance(s, mesher.Stage) for s in report.stages_run)
    assert report.stages_run[0] is mesher.Stage.PARSE
    assert report.records_by_code == {}, report.records_by_code
    assert mesher.read_report(report.path) == report
    try:
        report.schema_version = 0
    except Exception as e:
        assert type(e).__name__ == "FrozenInstanceError", e
    else:
        raise AssertionError("a report could be changed")

    # Options as keywords, enums or their values.
    opts = mesher.run(
        cube,
        out / "opts",
        element_size=0.2,
        formats=[mesher.OutputFormat.VTU, "exodus"],
        kernel="conformal",
        timeout=600,
    )
    assert opts.config.formats == (mesher.OutputFormat.EXODUS, mesher.OutputFormat.VTU)
    assert opts.config.element_size == 0.2
    assert (out / "opts" / "mesh.vtu").is_file() and not (out / "opts" / "mesh.msh").exists()
    assert opts.parts[0].mesh.tets < report.parts[0].mesh.tets
    for bad in [dict(element_size=0), dict(kernel="fastest"), dict(formats=["stl"])]:
        try:
            mesher.run(cube, out / "bad", **bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{bad} ran")

    # A project, and one that cannot be used.
    project = out / "cube.toml"
    project.write_text(f"input = '{cube}'\noutput = 'project_out'\n\n[settings]\nelement_size = 0.2\n")
    from_project = mesher.run_project(project, timeout=600)
    assert from_project.config.element_size == 0.2
    assert from_project.path == out / "project_out" / "report.json", from_project.path
    assert from_project.parts[0].mesh.tets == opts.parts[0].mesh.tets
    broken = mesher.run_project(out / "broken.toml")
    assert list(broken.records_by_code) == [mesher.DefectCode.PARSE_PROJECT_INVALID]
    record = broken.records_by_code[mesher.DefectCode.PARSE_PROJECT_INVALID][0]
    assert record.severity is mesher.Severity.FATAL, record.severity
    assert broken.exit_code == 30
    try:
        mesher.run_project(cube)
    except ValueError:
        pass
    else:
        raise AssertionError("an OBJ ran as a project")
    # The raw runners' reports, typed.
    assert run.report() == report


if __name__ == "__main__":
    main()
