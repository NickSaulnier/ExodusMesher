"""Download Thingi10k meshes and metadata for the benchmarks.

    python tests/thingi10k/fetch.py [ids.txt ...] [--jobs N]

Metadata goes to target/thingi10k/metadata/ and meshes to
target/thingi10k/raw/<id>.stl (or $THINGI10K_DIR), from the dataset's
Hugging Face mirror. A few models exist only as OBJ; those are saved as
raw/<id>.obj. Files already present are kept. The meshes carry their
authors' Thingiverse licences and are never committed.
"""

import argparse
import concurrent.futures
import os
import sys
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = "https://huggingface.co/datasets/Thingi10K/Thingi10K/resolve/main"
METADATA = ["geometry_data.csv", "input_summary.csv", "contextual_data.csv"]


def fetch(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return "kept"
    tmp = path + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    os.replace(tmp, path)
    return "fetched"


def fetch_mesh(raw, i):
    """The model's STL, or its OBJ when the dataset has no STL for it."""
    for ext in ("stl", "obj"):
        if os.path.exists(os.path.join(raw, f"{i}.{ext}")):
            return "kept"
    try:
        return fetch(f"{BASE}/raw_meshes/{i}.stl", os.path.join(raw, f"{i}.stl"))
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
    return fetch(f"{BASE}/raw_meshes/{i}.obj", os.path.join(raw, f"{i}.obj"))


def ids_from(path):
    with open(path, encoding="utf-8") as f:
        return [line.split()[0] for line in f if line.strip() and not line.startswith("#")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lists", nargs="*", default=[os.path.join(ROOT, "tests", "thingi10k", "subset.txt")])
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    data = os.environ.get("THINGI10K_DIR") or os.path.join(ROOT, "target", "thingi10k")
    os.makedirs(os.path.join(data, "metadata"), exist_ok=True)
    os.makedirs(os.path.join(data, "raw"), exist_ok=True)
    for name in METADATA:
        fetch(f"{BASE}/metadata/{name}", os.path.join(data, "metadata", name))
    ids = sorted({i for path in args.lists for i in ids_from(path)}, key=int)
    failed = []
    counts = {"kept": 0, "fetched": 0}
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = {
            pool.submit(fetch_mesh, os.path.join(data, "raw"), i): i
            for i in ids
        }
        for fut in concurrent.futures.as_completed(futures):
            try:
                counts[fut.result()] += 1
            except Exception as e:  # noqa: BLE001 - report and carry on
                failed.append(f"{futures[fut]}: {e}")
    print(f"{len(ids)} meshes: {counts['fetched']} fetched, {counts['kept']} already present, {len(failed)} failed")
    for line in failed:
        print("  " + line)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
