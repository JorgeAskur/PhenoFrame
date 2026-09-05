"""Unit tests for the pure-Python trait computation.

These tests target two layers:

1. The core trait math (`_leaf_traits_from_center_spline`) against synthetic
   center splines whose analytic length / orientation we can compute directly.
2. The full descriptor pipeline (`compute_traits_from_descriptor`) against
   small XML files written into a tmp_path with control-point overrides that
   pin the leaf geometry to known shapes.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from phenosuite import compute_traits_from_descriptor
from phenosuite.traits import _leaf_traits_from_center_spline


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_descriptor_with_override(
    xml_path: Path,
    ctrl_points: list[tuple[float, float, float]],
    *,
    leaf_length: float = 0.5,
    leaf_angle: float = 0.0,
    droopiness: float = 0.0,
    distance: float = 0.0,
    radius: float = 0.0,
) -> None:
    """Write a single-tiller, single-leaf descriptor with explicit ctrlCenter overrides.

    The override path makes leaf geometry deterministic regardless of azimuth
    noise — only the rigid azimuthal rotation around the stem axis varies, and
    that preserves the Y component of every tangent (so inclination and arc
    length stay invariant).
    """
    n = len(ctrl_points)
    ctrl_attrs = " ".join(
        f'p{i}x="{p[0]}" p{i}y="{p[1]}" p{i}z="{p[2]}"'
        for i, p in enumerate(ctrl_points)
    )
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<plant>
    <species name="Maize_Procedural">
        <Tiller type="main" radius="{radius}" alpha="0" beta="0" randomSeed="0" stemShrink="0">
            <leaves number="1">
                <leaf id="0" distance="{distance}" leafLength="{leaf_length}" leafWidth="0.05"
                      leafAngle="{leaf_angle}" droopiness="{droopiness}"
                      stemInclinationDeg="0" splinePoints="{max(n, 4)}"
                      useCtrlOverrides="1">
                    <ctrlCenter count="{n}" {ctrl_attrs}/>
                </leaf>
            </leaves>
        </Tiller>
    </species>
</plant>
"""
    xml_path.write_text(xml, encoding="utf-8")


# ---------------------------------------------------------------------------
# Direct math: _leaf_traits_from_center_spline
# ---------------------------------------------------------------------------

class TestLeafTraitsMath:
    def test_horizontal_two_point_spline(self) -> None:
        """Length equals Euclidean distance, inclination ~ 0 deg."""
        center = [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0)]
        traits = _leaf_traits_from_center_spline(center)

        assert traits["leaf_length"] == pytest.approx(0.5, abs=1e-12)
        assert traits["leaf_angle"]["inclination_deg"] == pytest.approx(0.0, abs=1e-9)
        assert traits["connection_point"] == {"x": 0.0, "y": 0.0, "z": 0.0}
        assert traits["tip_position"] == {"x": 0.5, "y": 0.0, "z": 0.0}

    def test_vertical_spline_inclination_90(self) -> None:
        """Pure-Y tangent gives inclination = 90 deg."""
        center = [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0)]
        traits = _leaf_traits_from_center_spline(center)

        assert traits["leaf_length"] == pytest.approx(0.5, abs=1e-12)
        assert traits["leaf_angle"]["inclination_deg"] == pytest.approx(90.0, abs=1e-9)
        assert traits["tip_position"]["y"] == pytest.approx(0.5, abs=1e-12)

    def test_polyline_length_sums_segments(self) -> None:
        """Length of an L-shaped polyline is the sum of leg lengths."""
        center = [(0.0, 0.0, 0.0), (0.3, 0.0, 0.0), (0.3, 0.4, 0.0)]
        traits = _leaf_traits_from_center_spline(center)
        assert traits["leaf_length"] == pytest.approx(0.7, abs=1e-12)

    def test_quarter_circle_arc_length(self) -> None:
        """Polyline length of a finely sampled quarter-circle approaches R*pi/2."""
        radius = 0.4
        n_samples = 200
        thetas = [i * (math.pi / 2.0) / (n_samples - 1) for i in range(n_samples)]
        center = [(radius * math.cos(t), radius * math.sin(t), 0.0) for t in thetas]

        traits = _leaf_traits_from_center_spline(center)
        analytic = radius * (math.pi / 2.0)
        # Polyline under-estimates by O(1/N^2). With N=200 the error is ~1e-5.
        assert traits["leaf_length"] == pytest.approx(analytic, rel=1e-4)

    def test_azimuth_normalized_to_0_360(self) -> None:
        """atan2-based azimuth is wrapped into [0, 360)."""
        # Tangent pointing toward -X gives atan2(0, -1) = 180 deg.
        center = [(0.0, 0.0, 0.0), (-0.3, 0.0, 0.0)]
        traits = _leaf_traits_from_center_spline(center)
        assert 0.0 <= traits["leaf_angle"]["azimuth_deg"] < 360.0
        assert traits["leaf_angle"]["azimuth_deg"] == pytest.approx(180.0, abs=1e-9)

    def test_empty_spline_raises(self) -> None:
        with pytest.raises(RuntimeError):
            _leaf_traits_from_center_spline([])


# ---------------------------------------------------------------------------
# Full pipeline: compute_traits_from_descriptor
# ---------------------------------------------------------------------------

class TestDescriptorPipeline:
    def test_horizontal_leaf_via_override(self, tmp_path: Path) -> None:
        """A leaf whose local control points lie along +Y (radial) is horizontal in world space."""
        xml = tmp_path / "horizontal.xml"
        # Local frame: px=0, py varies along radial (becomes horizontal in world),
        # pz=0 keeps the leaf out of the stem axis.
        _write_descriptor_with_override(
            xml,
            ctrl_points=[(0.0, 0.0, 0.0), (0.0, 0.25, 0.0), (0.0, 0.5, 0.0), (0.0, 0.75, 0.0)],
        )
        traits = compute_traits_from_descriptor(xml)

        assert len(traits) == 1
        t = traits[0]
        assert t["leaf_length"] == pytest.approx(0.75, abs=1e-9)
        assert t["leaf_angle"]["inclination_deg"] == pytest.approx(0.0, abs=1e-6)
        assert 0.0 <= t["leaf_angle"]["azimuth_deg"] < 360.0

    def test_vertical_leaf_via_override(self, tmp_path: Path) -> None:
        """A leaf whose local control points lie along +Z (stem-aligned) is vertical in world space."""
        xml = tmp_path / "vertical.xml"
        # pz varies → world Y varies → 90-degree inclination tangent.
        _write_descriptor_with_override(
            xml,
            ctrl_points=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.2), (0.0, 0.0, 0.4), (0.0, 0.0, 0.6)],
        )
        traits = compute_traits_from_descriptor(xml)

        assert len(traits) == 1
        t = traits[0]
        assert t["leaf_length"] == pytest.approx(0.6, abs=1e-9)
        assert t["leaf_angle"]["inclination_deg"] == pytest.approx(90.0, abs=1e-6)

    def test_quarter_circle_arc_via_override(self, tmp_path: Path) -> None:
        """Polyline length of a sampled arc descriptor matches analytic length."""
        xml = tmp_path / "arc.xml"
        radius = 0.5
        n = 32
        # Arc in local YZ plane: py = R*cos(theta), pz = R*sin(theta), theta in [0, pi/2].
        # Rotation around the stem axis preserves arc length.
        ctrl = []
        for i in range(n):
            theta = (math.pi / 2.0) * i / (n - 1)
            ctrl.append((0.0, radius * math.cos(theta), radius * math.sin(theta)))
        _write_descriptor_with_override(xml, ctrl_points=ctrl)
        traits = compute_traits_from_descriptor(xml)

        analytic = radius * (math.pi / 2.0)
        # n=32 sample polyline of a quarter-circle is within ~0.2% of analytic.
        assert traits[0]["leaf_length"] == pytest.approx(analytic, rel=5e-3)

    def test_sample_descriptor_yields_eight_leaves(self, sample_plant_xml: Path) -> None:
        """plants/plant_0.xml is a single 8-leaf tiller; verify counts and basic sanity."""
        traits = compute_traits_from_descriptor(sample_plant_xml)
        assert len(traits) == 8
        assert all(t["tiller_index"] == 0 for t in traits)
        assert [t["leaf_id"] for t in traits] == list(range(8))
        for t in traits:
            assert t["leaf_length"] > 0.0
            assert 0.0 <= t["leaf_angle"]["azimuth_deg"] < 360.0
            assert -90.0 <= t["leaf_angle"]["inclination_deg"] <= 90.0
