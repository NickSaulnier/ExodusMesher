"""Choose the stratified Thingi10k benchmark sample (Milestone 10).

    python tests/thingi10k/sample.py > tests/thingi10k/sample1000.txt

Reads the dataset's metadata (target/thingi10k/metadata/geometry_data.csv,
from https://huggingface.co/datasets/Thingi10K/Thingi10K) and picks 1,000
models of at most 100,000 faces. Each model gets one defect class, the first
that applies:

    self_intersecting  faces intersect
    non_manifold       an edge or vertex is non-manifold
    open               boundary edges (holes)
    multi_component    several connected components, otherwise clean
    clean              none of the above

Classes get fixed quotas; within a class, models are spread over three size
bins (under 1,000 faces, under 10,000, up to 100,000) in proportion to the
dataset, and picked with a fixed seed, so the sample is reproducible.
"""

import csv
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
META = os.path.join(ROOT, "target", "thingi10k", "metadata", "geometry_data.csv")
QUOTAS = {
    "clean": 350,
    "self_intersecting": 300,
    "multi_component": 150,
    "non_manifold": 120,
    "open": 80,
}
MAX_FACES = 100_000
SEED = 20261108


def defect_class(r):
    def yes(k):
        return r[k].strip().lower() in ("true", "1")

    def count(k):
        return int(float(r[k] or 0))

    if count("num_self_intersections") > 0:
        return "self_intersecting"
    if not yes("edge_manifold") or not yes("vertex_manifold"):
        return "non_manifold"
    if count("num_boundary_edges") > 0:
        return "open"
    if count("num_connected_components") > 1:
        return "multi_component"
    return "clean"


def size_bin(faces):
    return 0 if faces < 1_000 else 1 if faces < 10_000 else 2


def main():
    with open(META, encoding="utf-8", errors="replace") as f:
        rows = list(csv.DictReader(f))
    strata = {}
    for r in rows:
        faces = int(r["num_faces"])
        if faces > MAX_FACES:
            continue
        key = (defect_class(r), size_bin(faces))
        strata.setdefault(key, []).append((int(r["file_id"]), faces))
    rng = random.Random(SEED)
    chosen = []
    for cls, quota in QUOTAS.items():
        bins = [sorted(strata.get((cls, b), [])) for b in range(3)]
        total = sum(len(b) for b in bins)
        # Largest-remainder apportionment of the quota over the bins.
        shares = [quota * len(b) / total for b in bins]
        take = [int(s) for s in shares]
        for b in sorted(range(3), key=lambda b: shares[b] - take[b], reverse=True):
            if sum(take) == quota:
                break
            take[b] += 1
        for b, n in zip(bins, take):
            chosen += [(i, cls, faces) for i, faces in rng.sample(b, min(n, len(b)))]
    chosen.sort()
    out = sys.stdout
    out.write("# Stratified Thingi10k sample for the Milestone 10 benchmark; see sample.py.\n")
    out.write("# id class faces\n")
    for i, cls, faces in chosen:
        out.write(f"{i} {cls} {faces}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
