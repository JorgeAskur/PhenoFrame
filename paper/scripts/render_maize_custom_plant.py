"""Render a CHOSEN maize plant in four formats (voxel/skeleton/spline/procedural),
generating its meshes on the fly. Mirrors sorghum Fig. 1.

Usage:
    python render_maize_custom_plant.py --voxel-root <root> [--plant-id <id>]

If plant_id is omitted, picks the most-developed clean-fit plant
(max leaf count, tie-broken by lowest median RMS).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from _data_paths import GENERATED, MAIZE, REPO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--voxel-root", type=Path, required=True)
parser.add_argument("--plant-id", help="Plant identifier; default is the most developed clean fit")
parser.add_argument("--out-name", default="sec2_2_maize_plant_developed")
args = parser.parse_args()
os.environ["PHENOFRAME_VOXEL_PATH"] = str(args.voxel_root.resolve())
sys.path.insert(0, str(REPO / "experiments"))
sys.path.insert(0, str(REPO))

from experiments import data_loaders as dl  # noqa: E402
from phenoframe.skeleton_to_descriptor import (  # noqa: E402
    segment_skeleton, stem_aligned_rotation, voxel_to_world,
    skeleton_to_procedural_xml, skeleton_to_override_xml, skeleton_to_point_cloud_obj,
    DEFAULT_VOXEL_SIZE_M, DEFAULT_VERTICAL_AXIS,
)
from phenoframe.wrapper import from_xml_to_obj  # noqa: E402

OUT = MAIZE
FIG = GENERATED / "figures" / "maize"
FIG.mkdir(parents=True, exist_ok=True)
TMP = GENERATED / "visual_comparison_custom"
TMP.mkdir(exist_ok=True, parents=True)

SPLINE_POINTS = 8
VIS_STEM_RADIUS = 0.005
VIS_LEAF_WIDTH = 0.02

fits = pd.read_csv(OUT / "e5_fitted_params.csv")

# --- choose plant -----------------------------------------------------------
if args.plant_id:
    plant_id = args.plant_id
else:
    g = fits.groupby("plant_id").agg(n=("leaf_index", "max"),
                                     rms=("rms_error_cm", "median"))
    g["n"] = g["n"].astype(int) + 1
    # most developed with a still-clean fit (<=1.1 cm), tie-break by rms
    cand = g[g.rms <= 1.1].sort_values(["n", "rms"], ascending=[False, True])
    plant_id = cand.index[0]
out_name = args.out_name

rms_val = fits[fits.plant_id == plant_id]["rms_error_cm"].median()
n_leaves = int(fits[fits.plant_id == plant_id]["leaf_index"].max()) + 1
import re as _re
_m = _re.search(r"_\d+-\d+-([A-Za-z0-9]+)-(\d+)_", plant_id)
geno = f"{_m.group(1)} rep {_m.group(2)}" if _m else "?"
print(f"plant: {plant_id}")
print(f"  genotype~{geno}  leaves={n_leaves}  median RMS={rms_val:.2f} cm")

# --- generate the three meshes ----------------------------------------------
sk = dl.load_skeleton(plant_id)
p = str(TMP / "dev")
skeleton_to_point_cloud_obj(sk, p + "_skeleton.obj")
skeleton_to_procedural_xml(sk, p + "_procedural.xml", stem_radius=VIS_STEM_RADIUS,
                           spline_points=SPLINE_POINTS, leaf_width=VIS_LEAF_WIDTH)
assert from_xml_to_obj(p + "_procedural.xml", p + "_procedural.obj"), "proc OBJ failed"
skeleton_to_override_xml(sk, p + "_override.xml", stem_radius=VIS_STEM_RADIUS,
                         leaf_width=VIS_LEAF_WIDTH)
assert from_xml_to_obj(p + "_override.xml", p + "_override.obj"), "override OBJ failed"


def load_obj(path):
    verts, faces = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                q = line.split(); verts.append((float(q[1]), float(q[2]), float(q[3])))
            elif line.startswith("f "):
                idx = [int(q.split("/")[0]) - 1 for q in line.split()[1:]]
                for j in range(1, len(idx) - 1):
                    faces.append((idx[0], idx[j], idx[j + 1]))
    return np.array(verts, np.float32), np.array(faces, np.int32)


voxel_idx = dl.load_voxel_grid(plant_id)
skel_verts, _ = load_obj(p + "_skeleton.obj")
ovr_verts, ovr_faces = load_obj(p + "_override.obj")
proc_verts, proc_faces = load_obj(p + "_procedural.obj")
print(f"  voxels={len(voxel_idx):,} skel={len(skel_verts):,} "
      f"ovr={len(ovr_verts):,}v proc={len(proc_verts):,}v")

skel_int = dl.load_skeleton(plant_id)
seg = segment_skeleton(skel_int, vertical_axis=DEFAULT_VERTICAL_AXIS)
R, centroid = stem_aligned_rotation(skel_int, seg.stem_indices, DEFAULT_VERTICAL_AXIS)
voxel_m = voxel_to_world(voxel_idx, R, centroid, DEFAULT_VOXEL_SIZE_M)

ALL = np.concatenate([voxel_m, skel_verts, ovr_verts, proc_verts], axis=0)
bb_min, bb_max = ALL.min(axis=0), ALL.max(axis=0)
ctr = 0.5 * (bb_min + bb_max)
extent = (bb_max - bb_min).max() * 0.55


def set_axes(ax):
    ax.set_xlim(ctr[0] - extent, ctr[0] + extent)
    ax.set_ylim(ctr[2] - extent, ctr[2] + extent)
    ax.set_zlim(bb_min[1], bb_min[1] + 2 * extent)
    ax.set_box_aspect((1, 1, 1.4))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False; pane.set_edgecolor((1, 1, 1, 0))
    ax.grid(False); ax.view_init(elev=14, azim=-60)


def pts(ax, P, c, s, a=1.0):
    ax.scatter(P[:, 0], P[:, 2], P[:, 1], s=s, c=c, alpha=a, marker=".",
               edgecolors="none", rasterized=True)


def mesh(ax, V, F, c, a=0.95):
    ax.add_collection3d(Poly3DCollection(V[:, [0, 2, 1]][F], facecolor=c,
                                         edgecolor="none", linewidths=0.05, alpha=a))


LAB = ["(a) voxel grid", "(b) skeleton", "(c) skeleton descriptor",
       "(d) procedural descriptor"]
fig = plt.figure(figsize=(16, 5.5))
axes = [fig.add_subplot(1, 4, k + 1, projection="3d") for k in range(4)]
pts(axes[0], voxel_m, "#3a7fb0", 0.6, 0.5)
pts(axes[1], skel_verts, "#222222", 3.5, 0.9)
mesh(axes[2], ovr_verts, ovr_faces, "#7aaa66")
mesh(axes[3], proc_verts, proc_faces, "#cc6677")
for ax, lab in zip(axes, LAB):
    set_axes(ax)
    ax.text2D(0.5, -0.05, lab, transform=ax.transAxes, ha="center", va="top", fontsize=12)
fig.suptitle(f"One maize plant in four representations  "
             f"(median-tier fit, RMS = {rms_val:.2f} cm)",
             y=0.98, fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
outp = FIG / f"{out_name}.png"
fig.savefig(outp, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved: {outp}")
