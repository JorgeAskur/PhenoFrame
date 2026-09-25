"""Convert a 3D voxel skeleton into a PhenoFrame descriptor / spline representation.

Pipeline:

1. Segment the skeleton into stem and per-leaf voxel sets (delegates to
   :mod:`phenoframe.skeleton_traits`).
2. Compute the stem's principal directions via PCA; rotate the whole
   skeleton so the stem axis aligns with PhenoFrame's world Y axis (the
   convention used throughout the descriptor format and trait extractor).
3. Per leaf, order voxels by graph distance from the leaf-stem junction
   and downsample to a small set of control points for the leaf's center
   spline.
4. Either return the leaf splines as world-coordinate point lists (which
   can be fed directly to :func:`phenoframe.traits._leaf_traits_from_center_spline`)
   or write a full PhenoFrame descriptor XML with ``useCtrlOverrides=1`` and
   ``<ctrlCenter>`` elements for each leaf.

The spline representation uses world Y-up meters: voxel ``k`` axis maps
to world ``Y``; voxel ``i`` and ``j`` map to world ``X`` and ``Z`` after
the stem-alignment rotation. This matches PhenoFrame's documented coordinate
convention so that the values returned by the standard PhenoFrame trait
extractor are directly comparable to Mathieu/Jensina's ``θ`` (= 90° −
``inclination_deg``) and ``φ`` (= ``azimuth_deg``).
"""

from __future__ import annotations

from collections import deque
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .skeleton_traits import (
    _bfs_distances_within,
    _build_neighbors,
    _orient_stem_axis,
    _principal_directions,
    DEFAULT_VERTICAL_AXIS,
    segment_skeleton,
)


DEFAULT_VOXEL_SIZE_M = 0.002    # 2 mm/voxel — matches the Sorghum dataset's ~1m plant in 512^3
DEFAULT_N_LEAF_CTRL = 4         # number of control points per leaf spline (matches PhenoFrame default)
DEFAULT_FIT_LEAF_CTRL = 12      # higher resolution for inverse fitting (captures S-curves)


# ---------------------------------------------------------------------------
# Coordinate transform
# ---------------------------------------------------------------------------

def stem_aligned_rotation(
    voxels: np.ndarray,
    stem_indices: np.ndarray,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(R, centroid)`` such that ``(voxels - centroid) @ R.T`` puts the stem along world +Y.

    R is a 3×3 rotation matrix taking voxel-space coordinates to PhenoFrame world
    coordinates: stem PCA primary direction → +Y; PCA secondary → +X; the
    third axis (right-handed) → +Z. The stem-voxel centroid is subtracted
    so the plant is centered on the world origin's horizontal plane.
    """
    if len(stem_indices) < 2:
        raise ValueError(f"Need ≥ 2 stem voxels for PCA; got {len(stem_indices)}")
    stem_pts = voxels[stem_indices].astype(float)
    v1, v2, _ = _principal_directions(stem_pts)
    v1 = _orient_stem_axis(v1, vertical_axis)

    # Right-handed orthonormal basis (v2, v1, v2 × v1):
    #   v1 → +Y, v2 → +X, v2 × v1 → +Z
    v3 = np.cross(v2, v1)
    R = np.array([v2, v1, v3])  # rows: world-X-in-voxel, world-Y-in-voxel, world-Z-in-voxel
    centroid = stem_pts.mean(axis=0)
    return R, centroid


def voxel_to_world(
    voxels_subset: np.ndarray,
    R: np.ndarray,
    centroid: np.ndarray,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
) -> np.ndarray:
    """Apply the rotation + scale to a subset of voxel-space points → world meters."""
    return ((voxels_subset.astype(float) - centroid) @ R.T) * voxel_size_m


def _orient_v2_by_first_leaf(
    R: np.ndarray,
    centroid: np.ndarray,
    voxels: np.ndarray,
    seg,
) -> np.ndarray:
    """Orient v₂ (world X-axis) deterministically using the first leaf's direction.

    PCA defines v₂ only up to sign, which propagates as ±180° azimuthal
    flips between plants.  This function resolves the ambiguity by
    requiring the **lowest leaf's centroid** to have a non-negative X
    component in world space.  The choice is arbitrary but plant-intrinsic
    and deterministic, eliminating the random sign flips that previously
    created ±180° corner clusters in the φ scatter.
    """
    if not seg.leaf_indices or len(seg.leaf_indices[0]) == 0:
        return R

    # Centroid of the lowest leaf's voxels → robust direction estimate
    first_leaf_voxels = seg.leaf_indices[0]
    first_junction = seg.junction_indices[0]
    leaf_mean = voxels[first_leaf_voxels].astype(float).mean(axis=0)
    junction_pt = voxels[first_junction].astype(float)
    direction_world = (leaf_mean - junction_pt) @ R.T   # scale doesn't matter for sign

    if direction_world[0] < 0:
        R = R.copy()
        R[0] *= -1   # flip world-X  (v₂)
        R[2] *= -1   # flip world-Z  (v₃ = v₂ × v₁) to keep right-handed
    return R


# ---------------------------------------------------------------------------
# Per-leaf spline generation
# ---------------------------------------------------------------------------

def order_leaf_voxels(
    voxels: np.ndarray,
    leaf_voxel_indices: np.ndarray,
    junction_idx: int,
    neighbors: Optional[list[list[int]]] = None,
) -> np.ndarray:
    """Return the leaf voxel indices ordered by graph distance from the junction (ascending)."""
    if neighbors is None:
        neighbors = _build_neighbors(voxels)
    leaf_set = set(int(i) for i in leaf_voxel_indices)
    dist = _bfs_distances_within(junction_idx, neighbors, leaf_set | {junction_idx})
    leaf_dists = sorted([(dist.get(int(v), 10**9), int(v)) for v in leaf_voxel_indices])
    return np.array([v for _, v in leaf_dists], dtype=int)


def downsample_polyline(points: np.ndarray, n: int) -> np.ndarray:
    """Subsample ``points`` to ``n`` evenly-spaced indices (preserving first and last)."""
    if len(points) <= n:
        return points
    idx = np.linspace(0, len(points) - 1, n).round().astype(int)
    # np.linspace can produce duplicate indices when n is close to len(points); de-dup
    idx = np.unique(idx)
    return points[idx]


def skeleton_to_leaf_splines(
    voxels: np.ndarray,
    n_leaf_ctrl: int = DEFAULT_N_LEAF_CTRL,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
    include_junction: bool = True,
) -> tuple[list[list[tuple[float, float, float]]], list[int]]:
    """Segment a skeleton, build per-leaf center splines in PhenoFrame world coords.

    Parameters
    ----------
    voxels:
        ``(N, 3)`` integer voxel indices.
    n_leaf_ctrl:
        Number of control points per leaf spline. Lower values give a longer
        baseline tangent (closer to Mathieu's "first 6 cm" recipe); higher
        values track curvature more closely.
    vertical_axis:
        Voxel axis aligned with the stem (default 2 = ``k``).
    voxel_size_m:
        Physical size of one voxel edge in meters.
    include_junction:
        If True, the junction stem voxel is prepended as control point 0
        (PhenoFrame convention: ``center[0]`` is the leaf-stem connection point).

    Returns
    -------
    leaf_splines, leaf_voxel_counts
        ``leaf_splines[i]`` is a list of ``(x, y, z)`` tuples in PhenoFrame world
        meters. ``leaf_voxel_counts[i]`` is the original number of skeleton
        voxels making up the leaf (before downsampling), useful for QC.
    """
    seg = segment_skeleton(voxels, vertical_axis=vertical_axis)
    R, centroid = stem_aligned_rotation(voxels, seg.stem_indices, vertical_axis)
    neighbors = _build_neighbors(voxels)

    leaf_splines: list[list[tuple[float, float, float]]] = []
    leaf_voxel_counts: list[int] = []
    for leaf_voxel_indices, junction_idx in zip(seg.leaf_indices, seg.junction_indices):
        if len(leaf_voxel_indices) == 0:
            continue
        ordered = order_leaf_voxels(voxels, leaf_voxel_indices, junction_idx, neighbors=neighbors)
        # Optionally prepend junction so center[0] is the stem-leaf connection
        if include_junction:
            full_voxels = np.concatenate([[junction_idx], ordered])
        else:
            full_voxels = ordered
        if len(full_voxels) < 2:
            continue
        sampled = downsample_polyline(voxels[full_voxels], n_leaf_ctrl + (1 if include_junction else 0))
        world = voxel_to_world(sampled, R, centroid, voxel_size_m)
        leaf_splines.append([tuple(p) for p in world])
        leaf_voxel_counts.append(int(len(leaf_voxel_indices)))

    return leaf_splines, leaf_voxel_counts


# ---------------------------------------------------------------------------
# Trait extraction via PhenoFrame
# ---------------------------------------------------------------------------

def extract_traits_via_phenoframe(
    voxels: np.ndarray,
    n_leaf_ctrl: int = DEFAULT_N_LEAF_CTRL,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
) -> pd.DataFrame:
    """End-to-end: skeleton voxels → PhenoFrame spline → PhenoFrame trait extractor.

    Returns a DataFrame with one row per leaf, sorted bottom-to-top by junction
    height. Columns:

    - ``leaf_index`` (int)
    - ``leaf_length_m`` (float, meters)
    - ``azimuth_deg`` (PhenoFrame convention, ``[0, 360)``)
    - ``inclination_deg`` (PhenoFrame, angle of base tangent above horizontal)
    - ``theta`` — Mathieu/Jensina's ``θ`` ``= 90° − inclination_deg``
    - ``phi`` — Mathieu/Jensina's ``φ``, mapped from ``azimuth_deg`` to ``(-180°, 180°]``
    - ``leaf_voxel_count`` — number of skeleton voxels (pre-downsampling)
    - ``connection_y_m`` — Y coordinate of the connection point (≈ junction height)
    """
    from .traits import _leaf_traits_from_center_spline   # local import to avoid cycle at module load

    leaf_splines, vox_counts = skeleton_to_leaf_splines(
        voxels,
        n_leaf_ctrl=n_leaf_ctrl,
        vertical_axis=vertical_axis,
        voxel_size_m=voxel_size_m,
        include_junction=True,
    )

    rows = []
    for li, (spline, n_voxels) in enumerate(zip(leaf_splines, vox_counts)):
        if len(spline) < 2:
            continue
        traits = _leaf_traits_from_center_spline(spline)
        azimuth = traits["leaf_angle"]["azimuth_deg"]                  # [0, 360)
        # PhenoFrame's azimuth is atan2(tz, tx) — rotation around +Y is one chirality;
        # Mathieu/Jensina's φ uses the opposite rotational direction (empirically
        # verified by 6× drop in median |Δφ| against gold when sign is negated).
        # So Mathieu's φ = − PhenoFrame's azimuth, then wrap to (-180°, 180°].
        phi = ((180.0 - azimuth) % 360.0) - 180.0
        inclination = traits["leaf_angle"]["inclination_deg"]
        rows.append({
            "leaf_index": li,
            "leaf_length_m": traits["leaf_length"],
            "azimuth_deg": azimuth,
            "inclination_deg": inclination,
            "theta": 90.0 - inclination,                                # Mathieu's convention
            "phi": phi,                                                  # Mathieu's convention
            "leaf_voxel_count": n_voxels,
            "connection_y_m": traits["connection_point"]["y"],
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Inverse fitting: world spline -> procedural parameters
# ---------------------------------------------------------------------------

DEFAULT_FIT_SPLINE_POINTS = 8
DEFAULT_STEM_RADIUS = 0.015     # metres; typical sorghum stem radius
_N_RESAMPLE = 50                # arc-length resampling density for fitting
_N_STARTS = 5                   # multi-start restarts for optimization
_BASE_TANGENT_RESAMPLE = 20     # points used to estimate noisy leaf-base tangents
_BASE_TANGENT_FRACTION = 0.25   # first fraction of the leaf used for base PCA


def _resample_polyline(pts: np.ndarray, n: int) -> np.ndarray:
    """Resample a polyline to *n* points spaced evenly by arc length."""
    diffs = np.diff(pts, axis=0)
    seg_lens = np.linalg.norm(diffs, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_lens)])
    total = cum[-1]
    if total < 1e-12:
        return np.tile(pts[0], (n, 1))
    targets = np.linspace(0.0, total, n)
    out = np.empty((n, pts.shape[1]), dtype=float)
    j = 0
    for i, t in enumerate(targets):
        while j < len(cum) - 2 and cum[j + 1] < t:
            j += 1
        span = cum[j + 1] - cum[j]
        frac = (t - cum[j]) / span if span > 1e-12 else 0.0
        out[i] = pts[j] + frac * diffs[j]
    return out


def _estimate_base_tangent(
    leaf_pts: np.ndarray,
    junction: np.ndarray,
    n_resample: int = _BASE_TANGENT_RESAMPLE,
    fraction: float = _BASE_TANGENT_FRACTION,
) -> np.ndarray:
    """Estimate the leaf-base tangent from the first part of a leaf polyline.

    The first skeleton/control-point step is often noisy near the stem.  This
    helper fits a line to the basal portion of the leaf, then orients the line
    from base to tip.  It falls back to the junction-to-first-leaf vector for
    very short or degenerate inputs.
    """
    leaf_pts = np.asarray(leaf_pts, dtype=float)
    junction = np.asarray(junction, dtype=float)
    fallback = leaf_pts[0] - junction
    if len(leaf_pts) < 3:
        return fallback

    resampled = _resample_polyline(leaf_pts, n_resample)
    n_base = max(4, int(round(n_resample * fraction)))
    n_base = min(n_base, len(resampled))
    base = resampled[:n_base]
    centered = base - base.mean(axis=0)

    try:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        tangent = vh[0]
    except np.linalg.LinAlgError:
        return fallback

    forward = base[-1] - base[0]
    if np.dot(tangent, forward) < 0.0:
        tangent = -tangent
    if np.linalg.norm(tangent) < 1e-12:
        return fallback
    return tangent


def _wrap_angle_deg(angle: float) -> float:
    """Normalize an angle to ``[0, 360)`` degrees."""
    return float(angle % 360.0)


def _estimate_stem_base_y(
    voxels: np.ndarray,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
) -> float:
    """Estimate the basal stem Y coordinate in the PhenoFrame-aligned frame."""
    seg = segment_skeleton(voxels, vertical_axis=vertical_axis)
    R, centroid = stem_aligned_rotation(voxels, seg.stem_indices, vertical_axis)
    stem_world = voxel_to_world(voxels[seg.stem_indices], R, centroid, voxel_size_m)
    return float(np.min(stem_world[:, 1]))


def _leaf_local_basis(
    stem_dir: np.ndarray,
    azimuth_rad: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the (blade_x, blade_y, blade_z) orthonormal basis.

    Matches the frame used by ``_world_leaf_center`` in :mod:`phenoframe.traits`:
    blade_y = stem tangent, blade_z = radial direction, blade_x = tangential.
    """
    blade_y = stem_dir / np.linalg.norm(stem_dir)
    rhat = np.array([np.cos(azimuth_rad), 0.0, np.sin(azimuth_rad)])
    blade_x = np.cross(blade_y, rhat)
    if np.linalg.norm(blade_x) < 1e-6:
        rhat = np.array([0.0, 0.0, 1.0])
        blade_x = np.cross(blade_y, rhat)
    blade_x = blade_x / np.linalg.norm(blade_x)
    blade_z = np.cross(blade_x, blade_y)
    blade_z = blade_z / np.linalg.norm(blade_z)
    return blade_x, blade_y, blade_z


def _world_to_local(
    world_pts: np.ndarray,
    connection: np.ndarray,
    blade_x: np.ndarray,
    blade_y: np.ndarray,
    blade_z: np.ndarray,
) -> np.ndarray:
    """Project world-space points into the leaf-local frame.

    Returns ``(N, 3)`` array of ``(local_x, local_y, local_z)`` where
    local_y is radial and local_z is along the stem (matching the frame
    that ``_build_leaf_center_local`` generates into).
    """
    d = world_pts - connection
    return np.column_stack([d @ blade_x, d @ blade_z, d @ blade_y])


def _forward_leaf_local(
    leaf_angle: float,
    droopiness: float,
    leaf_length: float,
    spline_points: int,
    stem_radius: float,
    leaf_curl: float = 0.0,
) -> np.ndarray:
    """Thin wrapper: generate local-frame control points from procedural params."""
    import math
    from .traits import LeafDesc as _LeafDesc, _build_leaf_center_local

    ld = _LeafDesc(
        leaf_length=leaf_length,
        leaf_angle=leaf_angle,
        droopiness=droopiness,
        leaf_curl=leaf_curl,
        spline_points=spline_points,
    )
    pts = _build_leaf_center_local(ld, stem_radius)
    return np.array(pts, dtype=float)


def fit_procedural_params(
    world_spline: Sequence[tuple[float, float, float]],
    stem_direction: Optional[np.ndarray] = None,
    stem_radius: float = DEFAULT_STEM_RADIUS,
    spline_points: int = DEFAULT_FIT_SPLINE_POINTS,
) -> dict:
    """Fit procedural descriptor parameters from a world-space leaf center spline.

    Given a leaf center spline in world coordinates (as produced by
    :func:`skeleton_to_leaf_splines`), find the ``leaf_angle`` and
    ``droopiness`` values that best reproduce the observed spline shape
    when run through the forward procedural model
    (:func:`phenoframe.traits._build_leaf_center_local`).

    ``leaf_curl`` is always returned as 0.0 because curl does not affect
    the center spline in the authoritative C++ model — it only modifies
    the left/right edge splines.

    Uses multi-start L-BFGS-B with geometry-informed initial guesses to
    avoid local minima.  The base azimuth is estimated from a line fit to
    the basal leaf segment and held fixed so the fitted descriptor preserves
    the biologically measured leaf orientation.

    ``leaf_length`` is set from the observed arc length (exact); the
    azimuth is extracted from the base tangent.

    Parameters
    ----------
    world_spline:
        List of ``(x, y, z)`` tuples in PhenoFrame world metres (Y-up).
        ``world_spline[0]`` is the stem-leaf connection point.
    stem_direction:
        Unit vector of the stem at this leaf's junction (default ``+Y``).
    stem_radius:
        Stem radius at the leaf exit, metres (controls the offset of the
        first control point from the stem axis).
    spline_points:
        Number of procedural control points to use in the forward model.

    Returns
    -------
    dict with keys:
        ``leaf_length``, ``leaf_angle``, ``droopiness``, ``leaf_curl``,
        ``azimuth_deg``, ``azimuth_delta_deg`` (currently 0), ``inclination_deg``,
        ``distance`` (internode — set to 0, caller should compute from
        junction heights),
        ``spline_points``, ``residual`` (sum of squared distances after
        fitting), ``residual_per_point`` (mean squared distance per
        resampled point).
    """
    from scipy.optimize import minimize as _minimize

    pts = np.asarray(world_spline, dtype=float)
    if len(pts) < 3:
        raise ValueError("Need >= 3 spline points for fitting (junction + >=2 leaf)")

    # pts[0] = junction (stem centre); pts[1:] = leaf points from sheath exit
    junction = pts[0]
    leaf_pts = pts[1:]

    # -- robust base tangent -> initial azimuth + inclination --
    tangent = _estimate_base_tangent(leaf_pts, junction)
    tx, ty, tz = tangent
    azimuth0_deg = _wrap_angle_deg(float(np.degrees(np.arctan2(tz, tx))))
    inclination_deg = float(np.degrees(np.arctan2(ty, np.sqrt(tx * tx + tz * tz))))

    # -- leaf length from arc length of the leaf portion only --
    diffs = np.diff(leaf_pts, axis=0)
    leaf_length = float(np.sum(np.linalg.norm(diffs, axis=1)))
    if leaf_length < 1e-6:
        raise ValueError("Degenerate spline (zero arc length)")

    # -- build local frame --
    if stem_direction is None:
        stem_direction = np.array([0.0, 1.0, 0.0])
    blade_x, blade_y, blade_z = _leaf_local_basis(stem_direction, np.radians(azimuth0_deg))

    # -- initial world -> local frame for start guesses only --
    local_leaf0 = _world_to_local(leaf_pts, junction, blade_x, blade_y, blade_z)

    obs_resampled = _resample_polyline(local_leaf0, _N_RESAMPLE)

    # -- cost function: compare leaf-only local points with fixed basal azimuth --
    # NOTE: leaf_curl does NOT affect the center spline in the authoritative C++
    # forward model (it only modifies left/right edge splines via
    # LeafSectionEmitter).  We therefore fit only leaf_angle and droopiness
    # here and report leaf_curl as 0.0 (the XML value is preserved when the
    # descriptor is loaded, but it has no effect on the center-spline shape).
    def _cost(params: np.ndarray) -> float:
        la, droop = params
        gen = _forward_leaf_local(la, droop, leaf_length, spline_points, stem_radius, 0.0)
        gen_r = _resample_polyline(gen, _N_RESAMPLE)
        return float(np.sum((obs_resampled - gen_r) ** 2))

    # -- geometry-informed initial guesses --
    la0 = max(1.0, min(89.0, abs(inclination_deg)))
    # Estimate droop direction from the tip's local height relative to base
    tip_z = local_leaf0[-1, 2]
    base_z = local_leaf0[0, 2]
    droop_sign = 1.0 if tip_z > base_z else -1.0

    # Droopiness bounds: the forward model now uses droopiness directly (no
    # *5 multiplier). Typical XML values range from about -105 to +15 degrees.
    starts = [
        [la0, 0.0],                                          # neutral
        [la0, droop_sign * 25.0],                            # moderate droop
        [la0, -droop_sign * 25.0],                           # opposite droop
        [la0, droop_sign * 80.0],                            # heavy droop
        [max(1.0, la0 * 0.7), droop_sign * 50.0],           # varied angle
    ]
    bounds = [(0.5, 89.5), (-300.0, 300.0)]

    # -- multi-start optimization --
    best = None
    for x0 in starts:
        result = _minimize(
            _cost,
            x0=x0,
            method="L-BFGS-B",
            bounds=bounds,
        )
        if best is None or result.fun < best.fun:
            best = result

    return {
        "leaf_length": leaf_length,
        "leaf_angle": float(best.x[0]),
        "droopiness": float(best.x[1]),
        "leaf_curl": 0.0,
        "azimuth_deg": azimuth0_deg,
        "azimuth_delta_deg": 0.0,
        "inclination_deg": inclination_deg,
        "distance": 0.0,
        "spline_points": spline_points,
        "residual": float(best.fun),
        "residual_per_point": float(best.fun / _N_RESAMPLE),
    }


def fit_skeleton_procedural(
    voxels: np.ndarray,
    n_leaf_ctrl: int = DEFAULT_FIT_LEAF_CTRL,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
    stem_radius: float = DEFAULT_STEM_RADIUS,
    spline_points: int = DEFAULT_FIT_SPLINE_POINTS,
) -> pd.DataFrame:
    """End-to-end: skeleton voxels -> per-leaf procedural parameters.

    Returns a DataFrame with one row per leaf sorted bottom-to-top,
    containing the fitted procedural parameters (``leaf_angle``,
    ``droopiness``, ``leaf_length``) alongside the extracted/refined
    azimuth, internode distance, and fit residual. ``leaf_curl`` is
    always 0.0 (curl does not affect the center spline).

    Uses multi-start L-BFGS-B with geometry-informed initial guesses.
    Defaults to ``n_leaf_ctrl=12`` (instead of the trait pipeline's 4)
    to give the optimizer enough shape detail.

    Parameters
    ----------
    voxels:
        ``(N, 3)`` integer voxel indices.
    n_leaf_ctrl:
        Control points for skeleton -> world spline conversion.
        Defaults to 12 for fitting; the trait pipeline uses 4.
    vertical_axis:
        Voxel axis aligned with the stem (default 2 = ``k``).
    voxel_size_m:
        Physical size of one voxel edge in metres.
    stem_radius:
        Assumed stem radius for the procedural forward model.
    spline_points:
        Number of procedural control points in the forward model.

    Returns
    -------
    DataFrame with columns: ``leaf_index``, ``leaf_length``,
    ``leaf_angle``, ``droopiness``, ``leaf_curl``, ``azimuth_deg``,
    ``azimuth_delta_deg``, ``inclination_deg``, ``distance``, ``spline_points``,
    ``connection_y_m``, ``leaf_voxel_count``, ``residual``,
    ``residual_per_point``.
    """
    leaf_splines, vox_counts = skeleton_to_leaf_splines(
        voxels,
        n_leaf_ctrl=n_leaf_ctrl,
        vertical_axis=vertical_axis,
        voxel_size_m=voxel_size_m,
        include_junction=True,
    )
    stem_base_y = _estimate_stem_base_y(voxels, vertical_axis, voxel_size_m)

    rows: list[dict] = []
    prev_connection_y = stem_base_y
    for li, (spline, n_voxels) in enumerate(zip(leaf_splines, vox_counts)):
        if len(spline) < 2:
            continue
        try:
            params = fit_procedural_params(
                spline,
                stem_radius=stem_radius,
                spline_points=spline_points,
            )
        except (ValueError, RuntimeError):
            continue

        connection_y = spline[0][1]
        params["distance"] = max(0.0, connection_y - prev_connection_y)
        prev_connection_y = connection_y

        params["leaf_index"] = li
        params["connection_y_m"] = connection_y
        params["leaf_voxel_count"] = n_voxels
        rows.append(params)

    return pd.DataFrame(rows)


def skeleton_to_procedural_xml(
    voxels: np.ndarray,
    output_path: Path | str,
    n_leaf_ctrl: int = DEFAULT_FIT_LEAF_CTRL,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
    stem_radius: float = DEFAULT_STEM_RADIUS,
    spline_points: int = DEFAULT_FIT_SPLINE_POINTS,
    species_name: str = "Sorghum_Fitted",
    phyllotaxy_deg: float = 180.0,
    leaf_width: float = 0.08,
) -> pd.DataFrame:
    """Skeleton -> procedural descriptor XML (``useCtrlOverrides=0``).

    Writes a full PhenoFrame descriptor with fitted procedural parameters
    for each leaf. The descriptor can be loaded by the C++ engine or
    the Python trait pipeline to regenerate geometry or extract traits.

    Returns the same per-leaf DataFrame as :func:`fit_skeleton_procedural`.
    """
    import xml.etree.ElementTree as ET

    output_path = Path(output_path)
    df = fit_skeleton_procedural(
        voxels,
        n_leaf_ctrl=n_leaf_ctrl,
        vertical_axis=vertical_axis,
        voxel_size_m=voxel_size_m,
        stem_radius=stem_radius,
        spline_points=spline_points,
    )
    if df.empty:
        raise ValueError("No leaves could be fitted from this skeleton")

    root = ET.Element("plant")
    species = ET.SubElement(root, "species", name=species_name)
    tiller = ET.SubElement(species, "Tiller",
                           type="main",
                           radius=f"{stem_radius:.6f}",
                           alpha="0",
                           beta=f"{phyllotaxy_deg:.1f}",
                           stemShrink="0.001",
                           randomSeed="1337")
    leaves_elem = ET.SubElement(tiller, "leaves", number=str(len(df)))

    for _, row in df.iterrows():
        ET.SubElement(leaves_elem, "leaf",
                      id=str(int(row["leaf_index"])),
                      distance=f"{row['distance']:.6f}",
                      leafLength=f"{row['leaf_length']:.6f}",
                      leafWidth=f"{leaf_width:.6f}",
                      leafAngle=f"{row['leaf_angle']:.4f}",
                      droopiness=f"{row['droopiness']:.4f}",
                      leafCurl=f"{row['leaf_curl']:.4f}",
                      leafAzimuthDeg=f"{row['azimuth_deg']:.4f}",
                      stemInclinationDeg="0",
                      splinePoints=str(int(row["spline_points"])),
                      widthTaper="2.5",
                      leafTwist="0")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree = ET.ElementTree(root)
    try:
        ET.indent(tree, space="  ")
    except AttributeError:
        pass
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    return df


def skeleton_to_override_xml(
    voxels: np.ndarray,
    output_path: Path | str,
    n_leaf_ctrl: int = DEFAULT_N_LEAF_CTRL,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
    stem_radius: float = DEFAULT_STEM_RADIUS,
    species_name: str = "Sorghum_Override",
    phyllotaxy_deg: float = 180.0,
    leaf_width: float = 0.08,
) -> list[list[tuple[float, float, float]]]:
    """Skeleton -> spline-override descriptor XML (``useCtrlOverrides=1``).

    Writes a full PhenoFrame descriptor where each leaf carries its center
    spline as ``<ctrlCenter>`` control points in the leaf-local frame
    (the same frame that ``_build_leaf_center_local`` generates into).
    The C++ engine and Python trait pipeline transform these local points
    to world space using the stem position + azimuth rotation.

    Internode ``distance`` values and ``leafAzimuthDeg`` are set from the
    skeleton so that leaves are positioned and oriented correctly.

    Returns the list of world-space leaf splines.
    """
    import xml.etree.ElementTree as ET

    output_path = Path(output_path)
    leaf_splines, vox_counts = skeleton_to_leaf_splines(
        voxels,
        n_leaf_ctrl=n_leaf_ctrl,
        vertical_axis=vertical_axis,
        voxel_size_m=voxel_size_m,
        include_junction=True,
    )
    if not leaf_splines:
        raise ValueError("No leaves could be extracted from this skeleton")
    stem_base_y = _estimate_stem_base_y(voxels, vertical_axis, voxel_size_m)

    root = ET.Element("plant")
    species = ET.SubElement(root, "species", name=species_name)
    tiller = ET.SubElement(species, "Tiller",
                           type="main",
                           radius=f"{stem_radius:.6f}",
                           alpha="0",
                           beta=f"{phyllotaxy_deg:.1f}",
                           stemShrink="0.001",
                           randomSeed="1337")
    leaves_elem = ET.SubElement(tiller, "leaves", number=str(len(leaf_splines)))

    stem_dir = np.array([0.0, 1.0, 0.0])
    prev_connection_y = stem_base_y

    for li, (spline, n_voxels) in enumerate(zip(leaf_splines, vox_counts)):
        pts = np.asarray(spline, dtype=float)
        connection = pts[0]
        leaf_pts = pts[1:]

        diffs = np.diff(leaf_pts, axis=0)
        arc_length = float(np.sum(np.linalg.norm(diffs, axis=1)))
        if len(leaf_pts) >= 1:
            arc_length += float(np.linalg.norm(leaf_pts[0] - connection))

        connection_y = float(connection[1])
        distance = max(0.0, connection_y - prev_connection_y)
        prev_connection_y = connection_y

        tangent = _estimate_base_tangent(leaf_pts, connection) if len(leaf_pts) > 0 else stem_dir
        tx, _, tz = tangent
        azimuth_deg = float(np.degrees(np.arctan2(tz, tx))) % 360.0
        azimuth_rad = np.radians(azimuth_deg)

        blade_x, blade_y, blade_z = _leaf_local_basis(stem_dir, azimuth_rad)
        local_pts = _world_to_local(pts, connection, blade_x, blade_y, blade_z)

        leaf_elem = ET.SubElement(leaves_elem, "leaf",
                                  id=str(li),
                                  distance=f"{distance:.6f}",
                                  leafLength=f"{arc_length:.6f}",
                                  leafWidth=f"{leaf_width:.6f}",
                                  leafAngle="45",
                                  droopiness="0",
                                  leafAzimuthDeg=f"{azimuth_deg:.4f}",
                                  stemInclinationDeg="0",
                                  splinePoints=str(len(spline)),
                                  widthTaper="2.5",
                                  leafTwist="0",
                                  leafCurl="0",
                                  useCtrlOverrides="1")

        ctrl_attrs = {"count": str(len(local_pts))}
        for i, (lx, ly, lz) in enumerate(local_pts):
            ctrl_attrs[f"p{i}x"] = f"{lx:.6f}"
            ctrl_attrs[f"p{i}y"] = f"{ly:.6f}"
            ctrl_attrs[f"p{i}z"] = f"{lz:.6f}"
        ET.SubElement(leaf_elem, "ctrlCenter", **ctrl_attrs)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree = ET.ElementTree(root)
    try:
        ET.indent(tree, space="  ")
    except AttributeError:
        pass
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    return leaf_splines


def skeleton_to_point_cloud_obj(
    voxels: np.ndarray,
    output_path: Path | str,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    voxel_size_m: float = DEFAULT_VOXEL_SIZE_M,
) -> None:
    """Export a voxel skeleton as a point-cloud OBJ (vertices only, no faces).

    The skeleton is rotated to PhenoFrame's Y-up world frame using the same
    stem-alignment transform as the descriptor pipeline, so the point
    cloud is spatially registered with OBJ meshes produced from
    descriptors.
    """
    output_path = Path(output_path)
    seg = segment_skeleton(voxels, vertical_axis=vertical_axis)
    R, centroid = stem_aligned_rotation(voxels, seg.stem_indices, vertical_axis)
    world = voxel_to_world(voxels, R, centroid, voxel_size_m)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("# Skeleton point cloud (PhenoFrame Y-up world frame)\n")
        for x, y, z in world:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
