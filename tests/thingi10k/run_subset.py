"""Mesh a set of Thingi10k models with the `mesher` CLI in parallel and
summarise the outcomes: exit codes, kernel (conformal or stuffing), and the
records that sent a model to the fallback.

    python tests/thingi10k/run_subset.py [ids.txt] [--jobs N] [--timeout S]
        [--out DIR] [--arg=--flag ...]

The models are read from target/thingi10k/obj/<id>.obj (welded OBJ written by
the oracle harness) and meshed into target/thingi10k/run/<id>/ (or
DIR/<id>/); each `--arg` is passed on to `mesher run`. Needs a
release build of `mesher` and `mesher-worker`. The ids file defaults to
tests/thingi10k/subset.txt.
"""

import argparse
import collections
import concurrent.futures
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EXE = os.path.join(ROOT, "target", "release", "mesher" + (".exe" if os.name == "nt" else ""))


def run(model, out, timeout, extra):
    start = time.monotonic()
    proc = subprocess.run(
        [EXE, "run", model, "--out", out, "--timeout", str(timeout), *extra],
        capture_output=True,
        text=True,
    )
    seconds = time.monotonic() - start
    kernels, codes = [], []
    try:
        with open(os.path.join(out, "report.json"), encoding="utf-8") as f:
            report = json.load(f)
        for part in report.get("parts", []):
            mesh = part.get("mesh") or {}
            if mesh.get("kernel"):
                kernels.append(mesh["kernel"])
        for r in report.get("records", []):
            if r["code"].startswith("MESH_"):
                codes.append(r["code"])
    except (OSError, ValueError):
        pass
    return proc.returncode, seconds, kernels, codes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="?", default=os.path.join(ROOT, "tests", "thingi10k", "subset.txt"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--timeout", type=float, default=600)
    ap.add_argument("--out", default=None)
    ap.add_argument("--arg", action="append", default=[])
    args = ap.parse_args()
    with open(args.ids, encoding="utf-8") as f:
        ids = [line.split()[0] for line in f if line.strip() and not line.startswith("#")]
    data = os.path.join(ROOT, "target", "thingi10k")
    results = {}
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = {
            pool.submit(
                run,
                os.path.join(data, "obj", f"{i}.obj"),
                os.path.join(args.out or os.path.join(data, "run"), i),
                args.timeout,
                args.arg,
            ): i
            for i in ids
        }
        for fut in concurrent.futures.as_completed(futures):
            results[futures[fut]] = fut.result()
    exits = collections.Counter(r[0] for r in results.values())
    conformal = sorted(i for i, r in results.items() if r[2] and all(k == "delaunay" for k in r[2]))
    nomesh = sorted(i for i, r in results.items() if not r[2])
    codes = collections.Counter(c for r in results.values() for c in set(r[3]))
    total = sum(r[1] for r in results.values())
    print(f"models {len(results)}  conformal {len(conformal)}  no mesh {len(nomesh)}")
    print("exit codes", dict(sorted(exits.items())))
    print("mesh records", dict(sorted(codes.items())))
    print(f"total time {total:.0f} s, slowest:")
    for i, r in sorted(results.items(), key=lambda kv: -kv[1][1])[:5]:
        print(f"  {i}: {r[1]:.0f} s exit {r[0]} {r[2]} {sorted(set(r[3]))}")
    print("no mesh:", " ".join(nomesh))
    with open(os.path.join(data, "run", "summary.json"), "w", encoding="utf-8") as f:
        json.dump({i: {"exit": r[0], "seconds": round(r[1], 2), "kernels": r[2], "codes": r[3]}
                   for i, r in sorted(results.items())}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
