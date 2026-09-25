"""Per-leaf trait extraction from a 3D voxel skeleton.

Reproduces the recipe from Tross et al. (2021) PeerJ 9:e12628 and
Davis et al. (2025) Plant Phenomics 7:100023, both of which use the
Gaillard et al. voxel-carving + skeletonization pipeline:

1. Treat skeleton voxels as nodes in a 26-connected graph.
2. Identify the ground voxel (lowest position along the vertical axis).
3. For every other endpoint (degree-1 node), BFS the shortest path back
   to ground. Voxels visited by ≥ 2 paths are **stem**; voxels visited
   by exactly 1 path belong to the corresponding **leaf**.
4. PCA on stem voxels gives the stem's primary direction ``v₁,T`` and
   secondary direction ``v₂,T``.
5. Per leaf, take the first ``n`` voxels after the leaf-stem junction
   (default 20 ≈ 4 cm at the dataset's voxel resolution); PCA gives the
   leaf principal direction ``vL``.
6. Compute the polar angle ``θ = ∠(vL, v₁,T) ∈ [0°, 180°]`` and the
   azimuthal angle ``φ ∈ (-180°, 180°]`` measured in the plane
   perpendicular to ``v₁,T`` between ``v₂,T`` and the projection of ``vL``.

Coordinate convention: voxel indices are ``(i, j, k)`` integers in a
``[0, grid_size)³`` cube. The default vertical axis is ``k`` (axis 2),
matching the Sorghum Zenodo dataset orientation.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd


# Number of skeleton voxels at the leaf base used for the leaf-direction PCA.
# Tross 2021 used 20 voxels (≈ 4 cm at the Zenodo dataset's resolution).
DEFAULT_LEAF_PCA_N = 20

# Vertical axis index (0=i, 1=j, 2=k). For the Sorghum Zenodo dataset, k is up.
DEFAULT_VERTICAL_AXIS = 2


# ---------------------------------------------------------------------------
# Skeleton graph
# ---------------------------------------------------------------------------

def _build_neighbors(voxels: np.ndarray) -> list[list[int]]:
    """For each voxel, return the list of indices of its 26-connected neighbors."""
    voxel_set = {tuple(v): i for i, v in enumerate(voxels)}
    neighbors: list[list[int]] = [[] for _ in range(len(voxels))]
    offsets = [
        (di, dj, dk)
        for di in (-1, 0, 1)
        for dj in (-1, 0, 1)
        for dk in (-1, 0, 1)
        if (di, dj, dk) != (0, 0, 0)
    ]
    for idx, v in enumerate(voxels):
        for di, dj, dk in offsets:
            n = (v[0] + di, v[1] + dj, v[2] + dk)
            j = voxel_set.get(n)
            if j is not None:
                neighbors[idx].append(j)
    return neighbors


def _connected_components(neighbors: list[list[int]]) -> list[list[int]]:
    """Return the connected components of the skeleton graph."""
    n = len(neighbors)
    visited = [False] * n
    comps: list[list[int]] = []
    for start in range(n):
        if visited[start]:
            continue
        comp: list[int] = []
        queue = deque([start])
        visited[start] = True
        while queue:
            u = queue.popleft()
            comp.append(u)
            for v in neighbors[u]:
                if not visited[v]:
                    visited[v] = True
                    queue.append(v)
        comps.append(comp)
    return comps


def _bridge_components(
    voxels: np.ndarray,
    neighbors: list[list[int]],
    max_gap: float = 10.0,
) -> list[list[int]]:
    """Add virtual edges between nearby components to repair small skeleton gaps.

    Skeletons produced by voxel thinning (e.g. the Sorghum Zenodo dataset) are
    often disconnected — the stem and parts of the leaves end up in separate
    components. Gaillard et al. apply a joining step inside their pipeline; we
    approximate it here by repeatedly bridging the two closest components
    until either everything is connected or the next gap exceeds ``max_gap``
    (in voxel units, e.g. ~2mm at the dataset's resolution).

    Modifies ``neighbors`` in place and also returns it.
    """
    while True:
        comps = _connected_components(neighbors)
        if len(comps) <= 1:
            return neighbors
        # Find the closest pair across components — O(N²) brute force over voxel pairs;
        # fine for a few thousand voxels. For larger skeletons we'd switch to a KD-tree.
        comps.sort(key=len, reverse=True)
        best_dist: float = float("inf")
        best_pair: Optional[tuple[int, int]] = None
        # Compare only the largest component against the rest, then the next largest, etc.
        # This is O(C * N) where C = #components and N = #voxels — usually only a handful of comps.
        merged_membership = np.empty(len(voxels), dtype=int)
        for ci, comp in enumerate(comps):
            for v in comp:
                merged_membership[v] = ci
        # For each non-largest component, find its closest voxel in any other component.
        for ci, comp in enumerate(comps[1:], start=1):
            comp_pts = voxels[comp].astype(float)
            for vi, v_idx in enumerate(comp):
                v_pos = voxels[v_idx].astype(float)
                # Compare against all voxels not in this component
                others_mask = merged_membership != ci
                others = np.where(others_mask)[0]
                d2 = ((voxels[others].astype(float) - v_pos) ** 2).sum(axis=1)
                j_local = int(np.argmin(d2))
                d = float(np.sqrt(d2[j_local]))
                if d < best_dist:
                    best_dist = d
                    best_pair = (int(v_idx), int(others[j_local]))
        if best_pair is None or best_dist > max_gap:
            return neighbors
        a, b = best_pair
        neighbors[a].append(b)
        neighbors[b].append(a)


def _bfs_paths_from(source: int, neighbors: list[list[int]], targets: Sequence[int]) -> dict[int, list[int]]:
    """BFS from ``source``; return ``{target: path[from source to target]}`` for each target."""
    parents: dict[int, int] = {source: source}
    queue = deque([source])
    target_set = set(targets)
    found: dict[int, list[int]] = {}
    while queue and target_set:
        u = queue.popleft()
        if u in target_set:
            target_set.discard(u)
            # Reconstruct path
            path = [u]
            while parents[path[-1]] != path[-1]:
                path.append(parents[path[-1]])
            found[u] = path[::-1]
            if not target_set:
                break
        for v in neighbors[u]:
            if v not in parents:
                parents[v] = u
                queue.append(v)
    return found


def _bfs_distances_within(source: int, neighbors: list[list[int]], allowed: set[int]) -> dict[int, int]:
    """BFS from ``source`` restricted to the subgraph induced by ``allowed``."""
    if source not in allowed:
        return {}
    dist = {source: 0}
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v in neighbors[u]:
            if v in allowed and v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


# ---------------------------------------------------------------------------
# Segmentation: stem vs leaves
# ---------------------------------------------------------------------------

@dataclass
class SkeletonSegmentation:
    """Result of segmenting a skeleton into stem and per-leaf voxel sets."""

    stem_indices: np.ndarray            # (N_stem,) int indices into the input voxel array
    leaf_indices: list[np.ndarray]      # one (N_leaf_i,) int array per leaf, ordered by ascending junction height
    junction_indices: list[int]         # one int per leaf — the stem voxel adjacent to the leaf base
    ground_index: int                   # the index of the voxel taken as ground


def segment_skeleton(
    voxels: np.ndarray,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    max_bridge_gap: float = 10.0,
) -> SkeletonSegmentation:
    """Segment a skeleton into a stem and per-leaf voxel sets.

    Parameters
    ----------
    voxels:
        ``(N, 3)`` array of integer voxel indices.
    vertical_axis:
        Which axis (0, 1, or 2) is vertical. Default is 2 (= k axis), matching
        the Sorghum Zenodo dataset.
    max_bridge_gap:
        Maximum voxel-unit distance allowed when bridging disconnected
        skeleton components. Components further apart than this are left
        disconnected (and their voxels are excluded from analysis).

    Returns
    -------
    SkeletonSegmentation
    """
    if voxels.ndim != 2 or voxels.shape[1] != 3:
        raise ValueError(f"voxels must be (N, 3); got shape {voxels.shape}")

    neighbors = _build_neighbors(voxels)
    neighbors = _bridge_components(voxels, neighbors, max_gap=max_bridge_gap)
    degrees = np.array([len(n) for n in neighbors])

    # The ground voxel is the lowest endpoint (degree-1 node) along the vertical axis,
    # falling back to the lowest voxel overall if no degree-1 node is present.
    candidate_endpoints = np.where(degrees == 1)[0]
    if len(candidate_endpoints) == 0:
        candidate_endpoints = np.arange(len(voxels))
    ground = candidate_endpoints[np.argmin(voxels[candidate_endpoints, vertical_axis])]

    # All other endpoints are leaf tips
    tips = [int(i) for i in candidate_endpoints if i != ground]
    if not tips:
        raise ValueError("Skeleton has no leaf tips (only one endpoint found)")

    # BFS from ground to every tip; collect paths. Tips that are still
    # unreachable (because they lie in a stranded component too far away to
    # bridge) are quietly dropped.
    paths = _bfs_paths_from(ground, neighbors, tips)
    if len(paths) < len(tips):
        # Keep only the tips we could actually reach
        tips = [t for t in tips if t in paths]
        if not tips:
            raise ValueError("Skeleton has no leaf tips reachable from ground after bridging")

    # A voxel is "stem" iff it lies on ≥ 2 paths
    visit_count = np.zeros(len(voxels), dtype=int)
    for path in paths.values():
        for v in path:
            visit_count[v] += 1
    stem_mask = visit_count >= 2
    stem_indices = np.where(stem_mask)[0]

    # Each tip's leaf voxels = path voxels not in the stem
    # Use BFS within the leaf voxels (from the junction) for per-leaf ordering
    leaf_voxel_sets: list[np.ndarray] = []
    junction_indices: list[int] = []
    for tip in tips:
        path = paths[tip]
        # Walk from tip backward until we re-enter the stem; everything until
        # then is the leaf, and the stem voxel adjacent to the leaf base is the junction.
        leaf_path: list[int] = []
        junction = None
        for v in reversed(path):
            if stem_mask[v]:
                junction = v
                break
            leaf_path.append(v)
        if junction is None:
            # The whole path is non-stem (only 1 tip in the skeleton — pathological)
            junction = ground
        leaf_voxel_sets.append(np.array(leaf_path[::-1], dtype=int))
        junction_indices.append(int(junction))

    # Sort leaves bottom-to-top by junction height
    junction_heights = [voxels[j, vertical_axis] for j in junction_indices]
    order = np.argsort(junction_heights)
    leaf_voxel_sets = [leaf_voxel_sets[i] for i in order]
    junction_indices = [junction_indices[i] for i in order]

    return SkeletonSegmentation(
        stem_indices=stem_indices,
        leaf_indices=leaf_voxel_sets,
        junction_indices=junction_indices,
        ground_index=int(ground),
    )


# ---------------------------------------------------------------------------
# PCA + trait math
# ---------------------------------------------------------------------------

def _principal_directions(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the three principal directions of ``points`` (rows = samples), ordered by descending variance.

    Eigenvector signs are made deterministic with sklearn's convention: the
    component with the largest absolute value is forced to be positive.
    Without this, ``np.linalg.eigh`` can return either sign per call, which
    propagates as ±180° flips in the azimuthal angle.
    """
    if len(points) < 2:
        raise ValueError("Need at least 2 points for PCA")
    centered = points - points.mean(axis=0)
    cov = centered.T @ centered  # (3, 3)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvecs = eigvecs[:, order]
    # sklearn-style sign convention: pick the sign so the largest |component| is positive
    for i in range(eigvecs.shape[1]):
        max_abs = np.argmax(np.abs(eigvecs[:, i]))
        if eigvecs[max_abs, i] < 0:
            eigvecs[:, i] *= -1
    return eigvecs[:, 0], eigvecs[:, 1], eigvecs[:, 2]


def _orient_stem_axis(v1: np.ndarray, vertical_axis: int) -> np.ndarray:
    """Flip ``v1`` so that its component along the vertical axis is positive."""
    if v1[vertical_axis] < 0:
        return -v1
    return v1


def _orient_leaf_axis(vL: np.ndarray, leaf_voxels: np.ndarray, junction_voxel: np.ndarray) -> np.ndarray:
    """Flip ``vL`` so it points away from the junction (toward the leaf centroid)."""
    direction = leaf_voxels.mean(axis=0) - junction_voxel
    if vL @ direction < 0:
        return -vL
    return vL


def _polar_angle_deg(v_leaf: np.ndarray, v_stem: np.ndarray) -> float:
    """θ = angle between leaf and stem principal directions, in [0°, 180°]."""
    cos = np.clip(v_leaf @ v_stem, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos)))


def _azimuthal_angle_deg(v_leaf: np.ndarray, v1: np.ndarray, v2: np.ndarray) -> float:
    """φ ∈ (-180°, 180°] = signed angle from v2 to (vL projected onto plane perpendicular to v1) around v1.

    Computed via ``atan2(vL_perp · (v1 × v2), vL_perp · v2)``.
    """
    v_perp = v_leaf - (v_leaf @ v1) * v1
    norm = np.linalg.norm(v_perp)
    if norm < 1e-12:
        return 0.0
    v_perp = v_perp / norm
    cross = np.cross(v1, v2)
    return float(np.degrees(np.arctan2(v_perp @ cross, v_perp @ v2)))


# ---------------------------------------------------------------------------
# Top-level extractor
# ---------------------------------------------------------------------------

def extract_traits_from_skeleton(
    voxels: np.ndarray,
    leaf_pca_n: int = DEFAULT_LEAF_PCA_N,
    vertical_axis: int = DEFAULT_VERTICAL_AXIS,
    return_segmentation: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, SkeletonSegmentation]:
    """Extract per-leaf θ and φ from a 3D voxel skeleton.

    Parameters
    ----------
    voxels:
        ``(N, 3)`` array of integer voxel indices on a regular grid.
    leaf_pca_n:
        Number of skeleton voxels at the leaf base used for the leaf-direction
        PCA. Tross/Mathieu used 20 (≈ 4 cm); Jensina also computed 30 / 40 / 60
        for the 6cm/8cm/12cm CSV variants.
    vertical_axis:
        0, 1, or 2 — which axis is up. Default 2 (k-axis), matching the
        Sorghum Zenodo dataset.
    return_segmentation:
        If True, also return the underlying ``SkeletonSegmentation``.

    Returns
    -------
    DataFrame
        One row per leaf, ordered bottom-to-top by junction height. Columns:

        - ``leaf_index`` — 0-indexed
        - ``junction_i``, ``junction_j``, ``junction_k`` — voxel coordinates
          of the leaf-stem junction
        - ``height_norm`` — junction height normalized to ``[-1, 1]``
          (offset by skeleton centroid, scaled by half-extent along vertical axis)
        - ``theta`` — leaf insertion angle in degrees (∈ ``[0, 180]``)
        - ``phi`` — leaf azimuth in degrees (∈ ``(-180, 180]``)
        - ``leaf_voxel_count`` — number of skeleton voxels assigned to this leaf
    """
    seg = segment_skeleton(voxels, vertical_axis=vertical_axis)

    # Stem PCA → v1 (axial), v2 (in horizontal plane)
    if len(seg.stem_indices) < 2:
        raise ValueError(f"Stem has only {len(seg.stem_indices)} voxels; cannot do PCA")
    stem_pts = voxels[seg.stem_indices].astype(float)
    v1, v2, _ = _principal_directions(stem_pts)
    v1 = _orient_stem_axis(v1, vertical_axis)

    # Normalize for height: junction height relative to skeleton extent along vertical
    z_all = voxels[:, vertical_axis].astype(float)
    z_min, z_max = z_all.min(), z_all.max()
    z_center = 0.5 * (z_min + z_max)
    z_half = 0.5 * (z_max - z_min) if z_max > z_min else 1.0

    rows = []
    for li, (leaf_voxel_indices, junction_idx) in enumerate(zip(seg.leaf_indices, seg.junction_indices)):
        if len(leaf_voxel_indices) < 2:
            # Leaf too small for PCA; record placeholders and continue
            jv = voxels[junction_idx]
            rows.append(
                {
                    "leaf_index": li,
                    "junction_i": int(jv[0]),
                    "junction_j": int(jv[1]),
                    "junction_k": int(jv[2]),
                    "height_norm": float((jv[vertical_axis] - z_center) / z_half),
                    "theta": float("nan"),
                    "phi": float("nan"),
                    "leaf_voxel_count": int(len(leaf_voxel_indices)),
                }
            )
            continue

        # Order leaf voxels by graph distance from the junction (BFS within leaf-only subgraph
        # rooted at the junction's neighbor inside the leaf).
        # First collect all leaf voxels reachable from the junction via leaf voxels:
        leaf_set = set(int(i) for i in leaf_voxel_indices)
        # The BFS source is the junction itself (so the leaf extends outward),
        # but we only walk within {leaf voxels ∪ junction} to keep distances meaningful.
        neighbors_local = _build_neighbors(voxels[list(leaf_set | {junction_idx})])
        # Re-do BFS in original index space for clarity:
        all_neighbors = _build_neighbors(voxels)
        dist = _bfs_distances_within(junction_idx, all_neighbors, leaf_set | {junction_idx})
        # Sort leaf voxels by distance from junction; truncate to first leaf_pca_n
        leaf_dists = sorted([(dist.get(int(v), 10**9), int(v)) for v in leaf_voxel_indices])
        first_n = [v for _, v in leaf_dists[:leaf_pca_n]]

        if len(first_n) < 2:
            theta = phi = float("nan")
        else:
            base_pts = voxels[first_n].astype(float)
            vL, _, _ = _principal_directions(base_pts)
            vL = _orient_leaf_axis(vL, base_pts, voxels[junction_idx].astype(float))
            theta = _polar_angle_deg(vL, v1)
            phi = _azimuthal_angle_deg(vL, v1, v2)

        jv = voxels[junction_idx]
        rows.append(
            {
                "leaf_index": li,
                "junction_i": int(jv[0]),
                "junction_j": int(jv[1]),
                "junction_k": int(jv[2]),
                "height_norm": float((jv[vertical_axis] - z_center) / z_half),
                "theta": theta,
                "phi": phi,
                "leaf_voxel_count": int(len(leaf_voxel_indices)),
            }
        )

    df = pd.DataFrame(rows)
    if return_segmentation:
        return df, seg
    return df
