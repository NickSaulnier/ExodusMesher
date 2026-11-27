"""Milestone 20's exit run: mesh each input the synthetic suite wrote (with
`SUITE_INPUT_DIR=<inputs>`) with the CLI and with `mesher.run`, and require
the two runs to write the same files, byte for byte, `report.json` among
them, and to end with the same exit code.

    python -I crates/mesher-py/tests/python/compare_cli.py <inputs> <mesher> [--out DIR]

`<mesher>` is the CLI (`target/release/mesher`, its worker beside it); the
Python runs use the installed package's worker. Each run's outputs go to
`<out>/<run>/cli` and `<out>/<run>/python` (default `<out>`: `<inputs>`).
Exits 1 on any difference, listing each.
"""

import argparse
import json
import pathlib
import subprocess
import sys
import time

import mesher


def files(root):
    """Every file under `root`, by its path relative to it."""
    return {p.relative_to(root).as_posix(): p for p in sorted(root.rglob("*")) if p.is_file()}


def compare(label, cli_dir, py_dir, cli_code, report):
    problems = []
    if report.exit_code != cli_code:
        problems.append(f"{label}: exit code {cli_code} from the CLI, {report.exit_code} from Python")
    cli, py = files(cli_dir), files(py_dir)
    if "report.json" not in cli:
        problems.append(f"{label}: the CLI wrote no report.json")
    for name in sorted(cli.keys() | py.keys()):
        if name not in py:
            problems.append(f"{label}: {name} written by the CLI only")
        elif name not in cli:
            problems.append(f"{label}: {name} written by Python only")
        elif cli[name].read_bytes() != py[name].read_bytes():
            problems.append(f"{label}: {name} differs")
    if "report.json" in py and report.to_json().encode() != py["report.json"].read_bytes():
        problems.append(f"{label}: Report.to_json() is not the report.json written")
    return problems


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", type=pathlib.Path)
    parser.add_argument("cli", type=pathlib.Path)
    parser.add_argument("--out", type=pathlib.Path)
    args = parser.parse_args()
    out_root = (args.out or args.inputs).resolve()

    runs = sorted(p.parent for p in args.inputs.resolve().glob("*/options.json"))
    if not runs:
        sys.exit(f"no runs under {args.inputs} (set SUITE_INPUT_DIR for the synthetic suite)")
    problems, codes = [], {}
    started = time.monotonic()
    for run in runs:
        options = json.loads((run / "options.json").read_text(encoding="utf-8"))
        label, parts, size = options["label"], options["parts"], options["element_size"]
        # The same input path for both: the report records it as given.
        source = str(run / "input.obj")
        cli_dir, py_dir = out_root / run.name / "cli", out_root / run.name / "python"
        for d in (cli_dir, py_dir):
            if d.exists() and any(d.iterdir()):
                sys.exit(f"{d} is not empty: give an empty --out")
        command = [str(args.cli), "run", source, "--out", str(cli_dir), "--element-size", repr(size)]
        if parts:
            assert not any("," in p for p in parts), parts
            command += ["--parts", ",".join(parts)]
        cli_code = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
        report = mesher.run(source, py_dir, parts=parts, element_size=size)
        found = compare(label, cli_dir, py_dir, cli_code, report)
        problems += found
        codes[report.exit_code] = codes.get(report.exit_code, 0) + 1
        print(f"{'differs' if found else 'same':>8}  exit {report.exit_code:>2}  {label}", flush=True)
    elapsed = time.monotonic() - started
    by_code = ", ".join(f"{n} exit {c}" for c, n in sorted(codes.items()))
    print(f"{len(runs)} runs in {elapsed:.0f} s ({by_code}); {len(problems)} differences")
    for p in problems:
        print(p)
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
