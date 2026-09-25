"""Cache true voxel bounding-box height per plant (vertical extent of the raw
voxel grid along the up-axis), for use as the geometric-fidelity normalizer.

Usage: voxel_heights.py <species> <voxel_root> <e4_csv> <out_csv> [limit]
Height_m = (max-min voxel index along axis 2) * 0.002 m.
"""
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

species, voxel_root, e4_csv, out_csv = sys.argv[1:5]
limit = int(sys.argv[5]) if len(sys.argv) > 5 else 0

REPO = Path(__file__).resolve().parents[2]
os.environ["PHENOFRAME_VOXEL_PATH"] = voxel_root
sys.path.insert(0, str(REPO / "experiments"))
sys.path.insert(0, str(REPO))
from experiments import data_loaders as dl  # noqa: E402

VOXEL_SIZE_M = 0.002
VAXIS = 2

plants = pd.read_csv(e4_csv)["plant_id"].tolist()
if limit:
    plants = plants[:limit]

t0 = time.time()
rows, fail = [], 0
for i, pid in enumerate(plants):
    try:
        v = dl.load_voxel_grid(pid)
        h = (int(v[:, VAXIS].max()) - int(v[:, VAXIS].min())) * VOXEL_SIZE_M
        rows.append((pid, h))
    except Exception:
        fail += 1
    if limit == 0 and (i + 1) % 250 == 0:
        print(f"  {i+1}/{len(plants)}  ({time.time()-t0:.0f}s)", flush=True)

df = pd.DataFrame(rows, columns=["plant_id", "vox_height_m"])
df.to_csv(out_csv, index=False)
dt = time.time() - t0
print(f"{species}: {len(df)} heights, {fail} failed, {dt:.1f}s "
      f"({1000*dt/max(len(plants),1):.0f} ms/plant)  "
      f"median height {df.vox_height_m.median():.2f} m")
