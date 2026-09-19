"""Stage the calibrated maize reconstructions into the Sorghum VOXEL_PATH
layout the experiment loaders expect:
  <root>/reconstructed/<plant>/voxels.txt
  <root>/skeletons/<plant>/optim_skeleton.txt
  <root>/skeletons/<plant>/angles.txt   (written from traits.csv: h, theta, phi)
  <root>/dataset/<plant>/               (empty placeholder)
Uses hardlinks for the big files (no duplication). Only plants WITH traits are
staged (so require_skeleton=True passes)."""
import argparse
import os
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, required=True, help="Existing reconstruction directories")
parser.add_argument("--output-root", type=Path, required=True, help="New staged voxel root")
args = parser.parse_args()
SRC = args.source.resolve()
ROOT = args.output_root.resolve()
if not SRC.is_dir():
    parser.error(f"Source directory does not exist: {SRC}")
if ROOT == SRC or SRC in ROOT.parents:
    parser.error("Output root must not be inside the source directory")
for sub in ("reconstructed", "skeletons", "dataset"):
    (ROOT / sub).mkdir(parents=True, exist_ok=True)


def link(src, dst):
    if dst.exists():
        return
    try:
        os.link(src, dst)
    except OSError:
        import shutil; shutil.copy(src, dst)


n = 0
for d in SRC.iterdir():
    if not d.is_dir():
        continue
    vt, st, tr = d / "voxels.txt", d / "optim_skeleton.txt", d / "traits.csv"
    if not (vt.exists() and st.exists() and tr.exists()):
        continue
    try:
        df = pd.read_csv(tr)
        if df.empty:
            continue
    except Exception:
        continue
    rec = ROOT / "reconstructed" / d.name
    sk = ROOT / "skeletons" / d.name
    rec.mkdir(parents=True, exist_ok=True)
    sk.mkdir(parents=True, exist_ok=True)
    (ROOT / "dataset" / d.name).mkdir(parents=True, exist_ok=True)
    link(vt, rec / "voxels.txt")
    link(st, sk / "optim_skeleton.txt")
    # angles.txt: one tab-separated line, interleaved height, theta, phi per leaf
    vals = []
    for _, r in df.iterrows():
        vals += [f"{r['height_norm']:.6f}", f"{r['theta']:.4f}", f"{r['phi']:.4f}"]
    (sk / "angles.txt").write_text("\t".join(vals))
    n += 1

print(f"staged {n} maize plants into {ROOT}")
