"""Prepare the ABC sample for the feature and patch survival checks
(Milestone 18).

    python tests/abc/prepare.py [ids.txt] [--archives DIR] [--out DIR]

Needs chunk 0000 of the ABC dataset (Koch et al., 2019,
https://deep-geometry.github.io/abc-dataset/): abc_0000_obj_v00.7z and
abc_0000_feat_v00.7z, in target/abc/ by default, and 7-Zip (`7z` on the
path, or C:/Program Files/7-Zip/7z.exe). The models are their creators'
(Onshape public documents); they are never committed, and only aggregate
results by model id are published.

For each id (default: sample300.txt, 300 models drawn with seed 18 from the
chunk's models that have a features file and an OBJ of at most 20 MB), the
OBJ and its features YAML are extracted to <out>/raw/<id>/, and
<out>/obj/<id>.obj is written: the same vertices and triangles, one OBJ
object per connected body (without objects the parser makes each group a
part), and in it one group per CAD surface (`g p<k>_<type>`), so each surface
becomes a named patch. <out>/sharp/<id>.txt lists the edges of the curves
marked sharp, as 0-based vertex pairs. The features file's indices are
0-based (checked on the sample: every face in exactly one surface, every
sharp edge a mesh edge between two surfaces but one).

Then: cargo run --release -p mesher-synthetic --example survival -- <out>/obj/*.obj
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def seven_zip():
    exe = shutil.which("7z")
    if exe:
        return exe
    exe = r"C:\Program Files\7-Zip\7z.exe"
    if os.path.exists(exe):
        return exe
    sys.exit("7-Zip not found")


def listing(exe, archive):
    """Map model id to its path in the archive."""
    out = subprocess.run([exe, "l", "-slt", archive], capture_output=True,
                         text=True, check=True).stdout
    paths = {}
    for line in out.splitlines():
        if line.startswith("Path = "):
            p = line[7:].strip()
            if p.endswith((".obj", ".yml")):
                paths[p.replace("\\", "/").split("/")[0]] = p
    return paths


def extract(exe, archive, members, dest):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("\n".join(members) + "\n")
        names = f.name
    try:
        subprocess.run([exe, "x", archive, f"-o{dest}", f"-i@{names}", "-bd", "-y"],
                       check=True, stdout=subprocess.DEVNULL)
    finally:
        os.unlink(names)


def convert(obj, yml, out_obj, out_sharp):
    verts, faces = [], []
    with open(obj, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("v "):
                verts.append(line)
            elif line.startswith("f "):
                faces.append([int(t.split("/")[0]) for t in line.split()[1:4]])
    with open(yml, encoding="utf-8") as f:
        feat = yaml.load(f, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader)) or {}
    surfaces = feat.get("surfaces") or []
    curves = feat.get("curves") or []

    owner = [-1] * len(faces)
    for k, s in enumerate(surfaces):
        for i in s.get("face_indices", []):
            if 0 <= i < len(faces) and owner[i] < 0:
                owner[i] = k

    parent = list(range(len(verts)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, c in faces:
        for u, v in ((a, b), (b, c)):
            ru, rv = find(u - 1), find(v - 1)
            if ru != rv:
                parent[max(ru, rv)] = min(ru, rv)
    roots = {}
    body = [roots.setdefault(find(f[0] - 1), len(roots)) for f in faces]

    with open(out_obj, "w", encoding="utf-8", newline="\n") as w:
        w.writelines(verts)
        cur = None
        for i in sorted(range(len(faces)), key=lambda i: (body[i], owner[i], i)):
            if (body[i], owner[i]) != cur:
                if cur is None or body[i] != cur[0]:
                    w.write(f"o b{body[i]}\n")
                cur = (body[i], owner[i])
                if owner[i] < 0:
                    w.write("g none\n")
                else:
                    w.write(f"g p{owner[i]}_{surfaces[owner[i]].get('type', 'x')}\n")
            a, b, c = faces[i]
            w.write(f"f {a} {b} {c}\n")

    with open(out_sharp, "w", encoding="utf-8", newline="\n") as w:
        for c in curves:
            if c.get("sharp"):
                vi = c.get("vert_indices", [])
                for u, v in zip(vi, vi[1:]):
                    if u != v:
                        w.write(f"{u} {v}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="?", default=os.path.join(HERE, "sample300.txt"))
    ap.add_argument("--archives", default=os.path.join(ROOT, "target", "abc"))
    ap.add_argument("--out", default=os.path.join(ROOT, "target", "abc"))
    args = ap.parse_args()
    with open(args.ids, encoding="utf-8") as f:
        ids = [l.strip() for l in f if l.strip()]
    exe = seven_zip()
    raw = os.path.join(args.out, "raw")
    for kind in ("obj", "feat"):
        archive = os.path.join(args.archives, f"abc_0000_{kind}_v00.7z")
        paths = listing(exe, archive)
        missing = [i for i in ids if i not in paths]
        if missing:
            sys.exit(f"{len(missing)} ids not in {archive}, e.g. {missing[0]}")
        todo = [paths[i] for i in ids
                if not os.path.exists(os.path.join(raw, paths[i]))]
        if todo:
            extract(exe, archive, todo, raw)
    for d in ("obj", "sharp"):
        os.makedirs(os.path.join(args.out, d), exist_ok=True)
    for i in ids:
        files = os.listdir(os.path.join(raw, i))
        obj = next(f for f in files if f.endswith(".obj"))
        yml = next(f for f in files if f.endswith(".yml"))
        convert(os.path.join(raw, i, obj), os.path.join(raw, i, yml),
                os.path.join(args.out, "obj", f"{i}.obj"),
                os.path.join(args.out, "sharp", f"{i}.txt"))
    print(f"{len(ids)} models in {os.path.join(args.out, 'obj')}")


if __name__ == "__main__":
    main()
