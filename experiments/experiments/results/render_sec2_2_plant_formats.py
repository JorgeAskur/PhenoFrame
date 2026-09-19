"""Render one plant in four formats for paper §2.2.

Panels: (a) raw voxel grid, (b) skeleton point cloud, (c) full-spline mesh,
(d) procedural-fit mesh — all from the same plant, same viewing angle.

Picks the median-tier example (50th-percentile fit RMS) used by E5's
visual_comparison output. The skeleton/procedural/override OBJ files are
already in `outputs/visual_comparison/`; only the raw voxel grid is loaded
fresh from the Zenodo archive.

Run:
    python experiments/experiments/results/render_sec2_2_plant_formats.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from experiments import data_loaders as dl  # noqa: E402
from phenoframe.skeleton_to_descriptor import (  # noqa: E402
    segment_skeleton,
    stem_aligned_rotation,
    voxel_to_world,
    DEFAULT_VOXEL_SIZE_M,
    DEFAULT_VERTICAL_AXIS,
)

OUT = Path(__file__).resolve().parent / 'outputs'
FIG = Path(__file__).resolve().parent / 'figures'
VIS = OUT / 'visual_comparison'
FIG.mkdir(exist_ok=True, parents=True)

TIER = 'median'  # 'good' | 'median' | 'poor'

# ---------------------------------------------------------------------------
# Pick the plant for this tier (re-derives the E5 selection deterministically)
# ---------------------------------------------------------------------------
fits = pd.read_csv(OUT / 'e5_fitted_params.csv')
plant_rms = fits.groupby('plant_id')['rms_error_cm'].median().sort_values()
percentile = {'good': 0.10, 'median': 0.50, 'poor': 0.90}[TIER]
plant_id = plant_rms.index[int(len(plant_rms) * percentile)]
rms_val = plant_rms[plant_id]
print(f'{TIER}-tier plant: {plant_id}')
print(f'  median per-leaf RMS = {rms_val:.2f} cm')


# ---------------------------------------------------------------------------
# OBJ parser (vertices + triangular faces; ignores normals/UVs/materials)
# ---------------------------------------------------------------------------
def load_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (vertices Nx3, faces Mx3 as int indices into vertices). Faces
    are triangulated; OBJs with quads are split via fan."""
    verts: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    with open(path) as f:
        for line in f:
            if line.startswith('v '):
                parts = line.split()
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('f '):
                # face indices may include /uv/normal: 'f 1/2/3 4/5/6 ...'
                idx = [int(p.split('/')[0]) - 1 for p in line.split()[1:]]
                # fan-triangulate
                for j in range(1, len(idx) - 1):
                    faces.append((idx[0], idx[j], idx[j + 1]))
    return np.array(verts, dtype=np.float32), np.array(faces, dtype=np.int32)


# ---------------------------------------------------------------------------
# Load all four representations of the same plant
# ---------------------------------------------------------------------------
print('\nLoading data...')

# (a) Voxel grid: integer indices into a 512^3 grid, from Zenodo archive
voxel_idx = dl.load_voxel_grid(plant_id)  # (N, 3) int
print(f'  voxels:    {len(voxel_idx):,} occupied cells')

# (b) Skeleton: also integer indices (same 512^3 grid).
skel_verts, _ = load_obj(VIS / f'{TIER}_skeleton.obj')
print(f'  skeleton:  {len(skel_verts):,} points')

# (c) Override (spline) mesh: full B-spline through skeleton positions
ovr_verts, ovr_faces = load_obj(VIS / f'{TIER}_override.obj')
print(f'  override:  {len(ovr_verts):,} verts, {len(ovr_faces):,} faces')

# (d) Procedural mesh: 3-parameter fit
proc_verts, proc_faces = load_obj(VIS / f'{TIER}_procedural.obj')
print(f'  procedural:{len(proc_verts):,} verts, {len(proc_faces):,} faces')


# ---------------------------------------------------------------------------
# Common-frame transform for the raw voxel grid.
#
# The voxel grid comes in 512^3 integer indices. The skeleton-derived
# meshes are in metric Y-up world coordinates after a stem-aligned
# rotation (see phenoframe.skeleton_to_descriptor.skeleton_to_point_cloud_obj).
# We apply the exact same rotation + scale to the raw voxel grid so it
# overlays in the same frame.
# ---------------------------------------------------------------------------

skel_int = dl.load_skeleton(plant_id)  # int-grid skeleton (same plant as voxel grid)
seg = segment_skeleton(skel_int, vertical_axis=DEFAULT_VERTICAL_AXIS)
R, centroid = stem_aligned_rotation(skel_int, seg.stem_indices, DEFAULT_VERTICAL_AXIS)
voxel_m = voxel_to_world(voxel_idx, R, centroid, DEFAULT_VOXEL_SIZE_M)

sk_min, sk_max = skel_verts.min(axis=0), skel_verts.max(axis=0)
vx_min, vx_max = voxel_m.min(axis=0),    voxel_m.max(axis=0)
print(f'\nskeleton bbox (m):  {sk_min.round(3)} -> {sk_max.round(3)}')
print(f'voxel    bbox (m):  {vx_min.round(3)} -> {vx_max.round(3)}')


# ---------------------------------------------------------------------------
# Pick a consistent viewing angle + axes box for all four panels
# ---------------------------------------------------------------------------
ALL_PTS = np.concatenate([voxel_m, skel_verts, ovr_verts, proc_verts], axis=0)
bb_min = ALL_PTS.min(axis=0)
bb_max = ALL_PTS.max(axis=0)
ctr = 0.5 * (bb_min + bb_max)
extent = (bb_max - bb_min).max() * 0.55  # pad so the plant fits comfortably

# Y is "up" in PhenoFrame; we'll match that with matplotlib's z being "up" on screen
# by swapping axes when plotting (x->X, z->Y_screen-depth, y->Z_screen-vertical).
def set_axes(ax):
    ax.set_xlim(ctr[0] - extent, ctr[0] + extent)
    ax.set_ylim(ctr[2] - extent, ctr[2] + extent)
    ax.set_zlim(bb_min[1], bb_min[1] + 2 * extent)
    ax.set_box_aspect((1, 1, 1.4))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    # Hide pane edges for a clean paper look
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor((1, 1, 1, 0))
    ax.grid(False)
    ax.view_init(elev=14, azim=-60)


def plot_points(ax, P, color, size, alpha=1.0):
    # P is Nx3 in (x, y_world_up, z) PhenoFrame frame; remap to mpl (x, z, y_up)
    ax.scatter(P[:, 0], P[:, 2], P[:, 1],
               s=size, c=color, alpha=alpha, marker='.',
               edgecolors='none', rasterized=True)


def plot_mesh(ax, V, F, color, alpha=0.9, edgecolor=None):
    # Same axis remap as plot_points
    Vp = V[:, [0, 2, 1]]
    tris = Vp[F]
    poly = Poly3DCollection(tris, facecolor=color, edgecolor=edgecolor or 'none',
                            linewidths=0.05, alpha=alpha)
    ax.add_collection3d(poly)


# ---------------------------------------------------------------------------
# Render four panels
# ---------------------------------------------------------------------------
LABELS = ['(a) voxel grid', '(b) skeleton', '(c) full-spline descriptor', '(d) procedural descriptor']
fig = plt.figure(figsize=(16, 5.5))
axes = [fig.add_subplot(1, 4, k + 1, projection='3d') for k in range(4)]

# (a) Voxel grid — sparse cloud
plot_points(axes[0], voxel_m, color='#3a7fb0', size=0.6, alpha=0.5)

# (b) Skeleton — heavier centerline points, dark
plot_points(axes[1], skel_verts, color='#222222', size=3.5, alpha=0.9)

# (c) Full-spline mesh — green leaves
plot_mesh(axes[2], ovr_verts, ovr_faces, color='#7aaa66', alpha=0.95)

# (d) Procedural mesh — red leaves
plot_mesh(axes[3], proc_verts, proc_faces, color='#cc6677', alpha=0.95)

for ax, label in zip(axes, LABELS):
    set_axes(ax)
    # Use ax.text2D below the axes box so the label can't collide with the 3D axis lines
    ax.text2D(0.5, -0.05, label,
              transform=ax.transAxes, ha='center', va='top', fontsize=12)

fig.suptitle(f'One plant in four representations  '
             f'(median-tier fit, RMS = {rms_val:.2f} cm)',
             y=0.98, fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out_path = FIG / 'sec2_2_plant_four_formats.png'
fig.savefig(out_path, dpi=180, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved: {out_path}')
