"""Tests for the inverse procedural parameter fitting pipeline.

Verifies that fit_procedural_params recovers known parameters from
synthetically generated leaf splines, and that the full skeleton ->
procedural descriptor pipeline works end-to-end.
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from phenoframe.traits import LeafDesc, _build_leaf_center_local
from phenoframe.skeleton_to_descriptor import (
    fit_procedural_params,
    _resample_polyline,
    _estimate_base_tangent,
    _leaf_local_basis,
    _world_to_local,
    _forward_leaf_local,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_world_spline(
    leaf_angle: float,
    droopiness: float,
    leaf_length: float = 0.4,
    spline_points: int = 8,
    stem_radius: float = 0.015,
    azimuth_deg: float = 90.0,
    connection: tuple[float, float, float] = (0.0, 0.3, 0.0),
) -> list[tuple[float, float, float]]:
    """Generate a junction-prepended world-space leaf spline from known params.

    Mimics the output of ``skeleton_to_leaf_splines(include_junction=True)``:
    point 0 is the junction (stem centre), followed by leaf points starting
    at the sheath exit.
    """
    ld = LeafDesc(
        leaf_length=leaf_length,
        leaf_angle=leaf_angle,
        droopiness=droopiness,
        spline_points=spline_points,
    )
    local_pts = _build_leaf_center_local(ld, stem_radius)

    # Build local -> world transform (matching _world_leaf_center)
    stem_dir = np.array([0.0, 1.0, 0.0])
    az_rad = math.radians(azimuth_deg)
    rhat = np.array([math.cos(az_rad), 0.0, math.sin(az_rad)])

    blade_y = stem_dir
    blade_x = np.cross(blade_y, rhat)
    blade_x = blade_x / np.linalg.norm(blade_x)
    blade_z = np.cross(blade_x, blade_y)

    conn = np.array(connection)
    # Prepend junction (stem centre = connection point, before sheath exit)
    world = [tuple(conn)]
    for p in local_pts:
        w = conn + blade_x * p[0] + blade_y * p[2] + blade_z * p[1]
        world.append(tuple(w))
    return world


# ---------------------------------------------------------------------------
# Resampling
# ---------------------------------------------------------------------------

class TestResamplePolyline:
    def test_preserves_endpoints(self):
        pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]], dtype=float)
        resampled = _resample_polyline(pts, 10)
        np.testing.assert_allclose(resampled[0], [0, 0, 0], atol=1e-10)
        np.testing.assert_allclose(resampled[-1], [1, 1, 0], atol=1e-10)

    def test_output_length(self):
        pts = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float)
        for n in [5, 20, 100]:
            assert len(_resample_polyline(pts, n)) == n


class TestBaseTangentEstimate:
    def test_uses_basal_segment_not_single_noisy_step(self):
        pts = np.array([
            [0.00, 0.00, 0.00],
            [0.02, 0.03, 0.08],  # deliberately noisy lateral first step
            [0.00, 0.06, 0.16],
            [0.00, 0.09, 0.24],
            [0.00, 0.12, 0.32],
            [0.00, 0.15, 0.40],
        ])
        tangent = _estimate_base_tangent(pts, np.array([0.0, 0.0, -0.01]))
        tangent = tangent / np.linalg.norm(tangent)
        first_step = pts[1] - pts[0]
        first_step = first_step / np.linalg.norm(first_step)

        expected = np.array([0.0, 0.03, 0.08])
        expected = expected / np.linalg.norm(expected)
        assert abs(tangent[0]) < abs(first_step[0])
        assert abs(np.dot(tangent, expected)) > 0.95


# ---------------------------------------------------------------------------
# Local frame
# ---------------------------------------------------------------------------

class TestLocalFrame:
    def test_roundtrip_identity(self):
        stem = np.array([0.0, 1.0, 0.0])
        az = math.radians(45.0)
        bx, by, bz = _leaf_local_basis(stem, az)

        pts = np.array([[0.1, 0.5, 0.2], [-0.1, 0.3, 0.7]], dtype=float)
        conn = np.array([0.0, 0.3, 0.0])
        local = _world_to_local(pts, conn, bx, by, bz)

        # Reconstruct world from local: w = conn + bx*lx + bz*ly + by*lz
        for i in range(len(pts)):
            w = conn + bx * local[i, 0] + bz * local[i, 1] + by * local[i, 2]
            np.testing.assert_allclose(w, pts[i], atol=1e-10)


# ---------------------------------------------------------------------------
# Forward model wrapper
# ---------------------------------------------------------------------------

class TestForwardLeafLocal:
    def test_first_point_at_stem_radius(self):
        pts = _forward_leaf_local(45.0, 0.0, 0.5, 8, 0.015)
        assert pts[0, 0] == pytest.approx(0.0, abs=1e-10)
        assert pts[0, 1] == pytest.approx(0.015, abs=1e-10)
        assert pts[0, 2] == pytest.approx(0.0, abs=1e-10)

    def test_arc_length_matches_leaf_length(self):
        L = 0.45
        pts = _forward_leaf_local(35.0, -5.0, L, 8, 0.015)
        diffs = np.diff(pts, axis=0)
        arc = np.sum(np.linalg.norm(diffs, axis=1))
        assert arc == pytest.approx(L, rel=0.05)


# ---------------------------------------------------------------------------
# Inverse fitting round-trips
# ---------------------------------------------------------------------------

class TestFitProceduralParams:
    @pytest.mark.parametrize("leaf_angle,droopiness", [
        (30.0, 0.0),
        (45.0, -10.0),
        (60.0, 5.0),
        (15.0, -20.0),
        (75.0, 15.0),
    ])
    def test_recovers_known_params(self, leaf_angle, droopiness):
        world = _make_world_spline(
            leaf_angle=leaf_angle,
            droopiness=droopiness,
            leaf_length=0.4,
            spline_points=8,
            azimuth_deg=120.0,
        )
        result = fit_procedural_params(world, spline_points=8)

        assert result["leaf_angle"] == pytest.approx(leaf_angle, abs=2.0)
        assert result["droopiness"] == pytest.approx(droopiness, abs=2.0)
        assert result["leaf_length"] == pytest.approx(0.4, rel=0.05)

    def test_different_azimuths_give_same_shape_params(self):
        results = []
        for az in [0.0, 90.0, 180.0, 270.0]:
            world = _make_world_spline(leaf_angle=40.0, droopiness=-5.0, azimuth_deg=az)
            r = fit_procedural_params(world, spline_points=8)
            results.append(r)

        for r in results[1:]:
            assert r["leaf_angle"] == pytest.approx(results[0]["leaf_angle"], abs=2.0)
            assert r["droopiness"] == pytest.approx(results[0]["droopiness"], abs=2.0)

    def test_residual_is_near_zero_for_exact_input(self):
        world = _make_world_spline(leaf_angle=45.0, droopiness=0.0)
        result = fit_procedural_params(world, spline_points=8)
        assert result["residual_per_point"] < 1e-6

    def test_azimuth_extracted_correctly(self):
        for target_az in [0.0, 45.0, 135.0, 270.0]:
            world = _make_world_spline(leaf_angle=40.0, droopiness=0.0, azimuth_deg=target_az)
            result = fit_procedural_params(world, spline_points=8)
            az_diff = abs((result["azimuth_deg"] - target_az + 180) % 360 - 180)
            assert az_diff < 5.0
            assert "azimuth_delta_deg" in result
            assert abs(result["azimuth_delta_deg"]) < 1e-3

    def test_rejects_degenerate_input(self):
        with pytest.raises(ValueError, match="Need >= 3"):
            fit_procedural_params([(0, 0, 0), (0, 0.1, 0)])
