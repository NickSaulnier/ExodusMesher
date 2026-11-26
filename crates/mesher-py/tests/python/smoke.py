"""Stage 1's check, run against an installed wheel (no test framework needed):

    python -I crates/mesher-py/tests/python/smoke.py <output dir>

Meshes the cube fixture through the packaged worker, checks the report and
`report.json` agree, that progress arrives, and that a bad job raises before
anything runs.
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

    print(f"mesher {mesher.__version__}: cube meshed, {tets} tets, {len(stages)} progress messages")


if __name__ == "__main__":
    main()
