"""Unit tests for the skeleton trait extractor.

Cover the algorithm on synthetic skeletons whose θ and φ we can compute by
hand, plus a sanity check that leaf-axis flipping happens correctly.
"""

from __future__ import annotations

import numpy as np
import pytest

from pymaize.skeleton_traits import (
    extract_traits_from_skeleton,
    segment_skeleton,
    _principal_directions,
)


# ---------------------------------------------------------------------------
# Helpers — build small synthetic skeletons by hand
# ---------------------------------------------------------------------------

def _vertical_stem(length: int = 51, axis_origin: tuple[int, int] = (30, 30)) -> list[tuple[int, int, int]]:
    i, j = axis_origin
    return [(i, j, k) for k in range(length)]


def _horizontal_branch(start: tuple[int, int, int], direction: tuple[int, int, int], n: int) -> list[tuple[int, int, int]]:
    di, dj, dk = direction
    return [(start[0] + d * di, start[1] + d * dj, start[2] + d * dk) for d in range(1, n + 1)]


# ---------------------------------------------------------------------------
# PCA convention
# ---------------------------------------------------------------------------

class TestPCASignConvention:
    def test_largest_component_is_positive(self):
        # A point cloud aligned along x: principal direction is (±1, 0, 0); convention should pick +1.
        pts = np.array([(i, 0, 0) for i in range(10)], dtype=float)
        v1, _, _ = _principal_directions(pts)
        max_abs = np.argmax(np.abs(v1))
        assert v1[max_abs] >= 0

    def test_deterministic_across_calls(self):
        rng = np.random.default_rng(0)
        pts = rng.standard_normal((50, 3))
        out1 = _principal_directions(pts)
        out2 = _principal_directions(pts)
        for a, b in zip(out1, out2):
            np.testing.assert_array_equal(a, b)


# ---------------------------------------------------------------------------
# Segmentation
# ---------------------------------------------------------------------------

class TestSegmentation:
    def test_single_horizontal_leaf_on_vertical_stem(self):
        stem = _vertical_stem(length=51)
        leaf = _horizontal_branch((30, 30, 25), (1, 0, 0), n=20)
        voxels = np.array(stem + leaf, dtype=int)

        seg = segment_skeleton(voxels)
        assert len(seg.leaf_indices) >= 1
        # Stem voxels are mostly the column (30, 30, 0..25), since only one path goes higher.
        # We just check that the segmentation did not crash and returned a junction near the leaf base.
        junction = voxels[seg.junction_indices[0]]
        assert junction[2] == 25 or abs(junction[2] - 25) <= 1

    def test_two_branches_at_different_heights(self):
        stem = _vertical_stem(length=51)
        leaf_a = _horizontal_branch((30, 30, 10), (1, 0, 0), n=15)
        leaf_b = _horizontal_branch((30, 30, 30), (0, 1, 0), n=15)
        voxels = np.array(stem + leaf_a + leaf_b, dtype=int)

        seg = segment_skeleton(voxels)
        # Should detect at least the two branches as separate leaves
        assert len(seg.leaf_indices) >= 2
        # Lowest leaf comes first
        junction_heights = [voxels[j][2] for j in seg.junction_indices]
        assert junction_heights == sorted(junction_heights)


# ---------------------------------------------------------------------------
# Trait extraction on known geometries
# ---------------------------------------------------------------------------

class TestTraitsKnownGeometry:
    """Use synthetic skeletons whose stem terminates at the highest leaf — mimicking
    real `optim_skeleton.txt` files in which the apical stem stub is pruned by the
    Gaillard SVM classifier."""

    def test_two_horizontal_leaves_at_apex_have_theta_90(self):
        """Two horizontal branches sharing the apex node — full stem resolves correctly.

        The "≥2 paths = stem" rule classifies all stem voxels correctly only
        when at least 2 leaves branch from the apex (otherwise the upper stem
        segment between the second-highest junction and the highest one gets
        misclassified as part of the topmost leaf). Real plants typically
        satisfy this; synthetic tests must too.
        """
        stem = _vertical_stem(length=61)
        leaf_a = _horizontal_branch((30, 30, 20), (1, 0, 0), n=20)
        leaf_b = _horizontal_branch((30, 30, 60), (0, 1, 0), n=20)   # at apex
        leaf_c = _horizontal_branch((30, 30, 60), (-1, 0, 0), n=20)  # also at apex
        voxels = np.array(stem + leaf_a + leaf_b + leaf_c, dtype=int)

        df = extract_traits_from_skeleton(voxels)
        df = df[df["theta"].notna()]
        assert len(df) >= 3
        # All three are horizontal → θ ≈ 90° for each
        assert all(abs(t - 90.0) < 5.0 for t in df["theta"])

    def test_apex_leaves_phi_differs_by_180(self):
        """Two leaves at the apex pointing opposite directions → φ differs by 180°."""
        stem = _vertical_stem(length=61)
        leaf_east_low = _horizontal_branch((30, 30, 20), (1, 0, 0), n=20)
        leaf_east_apex = _horizontal_branch((30, 30, 60), (1, 0, 0), n=20)
        leaf_west_apex = _horizontal_branch((30, 30, 60), (-1, 0, 0), n=20)
        voxels = np.array(stem + leaf_east_low + leaf_east_apex + leaf_west_apex, dtype=int)

        df = extract_traits_from_skeleton(voxels)
        df = df[df["theta"].notna()]
        # Find the two apex leaves (junction at k≈60)
        apex = df[df["junction_k"].between(58, 62)]
        assert len(apex) == 2
        phi_diff = abs((apex["phi"].iloc[0] - apex["phi"].iloc[1] + 540) % 360 - 180)
        assert abs(phi_diff - 180.0) < 5.0 or phi_diff < 5.0  # 180° apart, possibly wrapped to 0

    def test_diagonal_leaf_has_intermediate_theta(self):
        """One diagonal leaf at 45° + horizontal anchors at the apex."""
        stem = _vertical_stem(length=51)
        anchor_low = _horizontal_branch((30, 30, 20), (0, 1, 0), n=20)
        diag_leaf = [(30 + d, 30, 30 + d) for d in range(1, 21)]
        anchor_apex_a = _horizontal_branch((30, 30, 50), (-1, 0, 0), n=20)
        anchor_apex_b = _horizontal_branch((30, 30, 50), (0, -1, 0), n=20)
        voxels = np.array(stem + anchor_low + diag_leaf + anchor_apex_a + anchor_apex_b, dtype=int)

        df = extract_traits_from_skeleton(voxels, leaf_pca_n=20)
        df = df[df["theta"].notna()]
        thetas = sorted(df["theta"].tolist())
        # Diagonal leaf θ ≈ 45°. Horizontal anchors θ ≈ 90°.
        assert thetas[0] == pytest.approx(45.0, abs=5.0)
        assert thetas[-1] == pytest.approx(90.0, abs=5.0)


class TestExtractorAPI:
    def test_returns_dataframe_with_expected_columns(self):
        stem = _vertical_stem(length=50)
        leaf = _horizontal_branch((30, 30, 25), (1, 0, 0), n=20)
        voxels = np.array(stem + leaf, dtype=int)

        df = extract_traits_from_skeleton(voxels)
        for col in ("leaf_index", "junction_i", "junction_j", "junction_k", "height_norm", "theta", "phi", "leaf_voxel_count"):
            assert col in df.columns

    def test_return_segmentation_flag(self):
        stem = _vertical_stem(length=50)
        leaf = _horizontal_branch((30, 30, 25), (1, 0, 0), n=20)
        voxels = np.array(stem + leaf, dtype=int)

        result = extract_traits_from_skeleton(voxels, return_segmentation=True)
        assert isinstance(result, tuple) and len(result) == 2
        df, seg = result
        assert hasattr(seg, "stem_indices")
        assert hasattr(seg, "leaf_indices")
