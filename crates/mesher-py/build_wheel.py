"""Build the `exodusmesher` wheel: the release worker first, copied into the
package (`python/mesher/bin/`), then the wheel with maturin.

    python -I crates/mesher-py/build_wheel.py [--out DIR] [-j N]

Needs maturin on PATH (or importable). Cargo runs with `-j N` (default 4),
as the dev machine cannot take two all-core builds at once.
"""

import argparse
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(ROOT / "target" / "wheels"))
    parser.add_argument("-j", "--jobs", default="4")
    args = parser.parse_args()

    subprocess.run(
        ["cargo", "build", "--release", "-j", args.jobs, "-p", "mesher-worker"],
        cwd=ROOT,
        check=True,
    )
    exe = "mesher-worker.exe" if sys.platform == "win32" else "mesher-worker"
    bin_dir = HERE / "python" / "mesher" / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "target" / "release" / exe, bin_dir / exe)

    maturin = shutil.which("maturin")
    command = [maturin] if maturin else [sys.executable, "-m", "maturin"]
    subprocess.run(
        command
        + ["build", "--release", "-j", args.jobs, "--manifest-path", str(HERE / "Cargo.toml"), "--out", args.out],
        cwd=ROOT,
        check=True,
    )


if __name__ == "__main__":
    main()
