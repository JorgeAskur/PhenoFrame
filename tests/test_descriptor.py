"""Round-trip and cross-pipeline consistency tests for plant descriptors.

These tests verify three things:

1. Parsing a descriptor XML twice yields identical traits (parser determinism).
2. Building a plant via the C++ wrapper, saving its XML, reloading it, and
   regenerating it produces the same spline traits within tight tolerance
   (round-trip preserves geometry).
3. The pure-Python trait pipeline and the C++ wrapper agree on the leaf-stem
   connection point for the bundled sample plant within a loose tolerance
   (the two pipelines build leaf splines independently from the same
   descriptor; they should land within a few millimeters of each other).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from phenosuite import Leaf, Maize, Tiller, compute_traits_from_descriptor

from .conftest import requires_wrapper


def _build_three_leaf_plant(gen: Maize) -> None:
    gen.reset()
    gen.set_species("maize")
    gen.set_global_settings(
        density_v=80.0,
        density_u=80.0,
        stem_density_scale=0.15,
        stem_row_override=0,
        stem_ribbon_spacing=0.15,
    )
    t = gen.add_tiller(Tiller(random_seed=42).to_ctype())
    for i, leaf in enumerate([
        Leaf(id=0, distance=0.18, leaf_length=0.40, leaf_width=0.07, leaf_angle=30.0, droopiness=-15.0),
        Leaf(id=1, distance=0.22, leaf_length=0.55, leaf_width=0.09, leaf_angle=42.0, droopiness=-10.0),
        Leaf(id=2, distance=0.30, leaf_length=0.50, leaf_width=0.10, leaf_angle=50.0, droopiness=-5.0),
    ]):
        gen.add_leaf(t, leaf.to_ctype())
    gen.rebuild()


class TestPythonParserDeterminism:
    def test_parsing_same_xml_twice_is_identical(self, sample_plant_xml: Path) -> None:
        first = compute_traits_from_descriptor(sample_plant_xml)
        second = compute_traits_from_descriptor(sample_plant_xml)

        assert len(first) == len(second)
        for a, b in zip(first, second):
            assert a["leaf_length"] == pytest.approx(b["leaf_length"], abs=1e-12)
            assert a["leaf_angle"]["azimuth_deg"] == pytest.approx(b["leaf_angle"]["azimuth_deg"], abs=1e-12)
            assert a["leaf_angle"]["inclination_deg"] == pytest.approx(b["leaf_angle"]["inclination_deg"], abs=1e-12)
            for axis in ("x", "y", "z"):
                assert a["connection_point"][axis] == pytest.approx(b["connection_point"][axis], abs=1e-12)
                assert a["tip_position"][axis] == pytest.approx(b["tip_position"][axis], abs=1e-12)


@requires_wrapper
@pytest.mark.wrapper
class TestWrapperXmlRoundTrip:
    def test_save_then_load_preserves_spline_traits(self, tmp_path: Path) -> None:
        """Round-trip via the C++ wrapper: build → save_xml → load_xml → traits match."""
        with Maize() as gen:
            _build_three_leaf_plant(gen)
            originals = [gen.get_leaf_spline_traits(i) for i in range(gen.leaf_spline_trait_count())]
            xml_path = tmp_path / "round_trip.xml"
            assert gen.save_xml(str(xml_path)) is True

        with Maize() as reloaded:
            assert reloaded.load_xml(str(xml_path)) is True
            assert reloaded.leaf_spline_trait_count() == len(originals)
            for i, original in enumerate(originals):
                rt = reloaded.get_leaf_spline_traits(i)
                assert rt["leaf_length"] == pytest.approx(original["leaf_length"], rel=1e-3, abs=1e-4)
                for key in ("connection_point", "tip_position"):
                    for axis in ("x", "y", "z"):
                        assert rt[key][axis] == pytest.approx(original[key][axis], rel=1e-3, abs=1e-4)

    def test_round_tripped_xml_is_parseable_by_python_pipeline(self, tmp_path: Path) -> None:
        """The XML the C++ engine writes must remain consumable by the Python parser."""
        with Maize() as gen:
            _build_three_leaf_plant(gen)
            xml_path = tmp_path / "for_python.xml"
            assert gen.save_xml(str(xml_path)) is True

        traits = compute_traits_from_descriptor(xml_path)
        assert len(traits) == 3
        assert all(t["leaf_length"] > 0 for t in traits)


@requires_wrapper
@pytest.mark.wrapper
class TestCrossPipelineConsistency:
    def test_connection_points_agree_within_tolerance(self, sample_plant_xml: Path) -> None:
        """Python and C++ trait pipelines should locate leaf-stem joints within ~1cm."""
        py_traits = compute_traits_from_descriptor(sample_plant_xml)

        with Maize() as gen:
            assert gen.load_xml(str(sample_plant_xml)) is True
            cpp_traits = [gen.get_leaf_spline_traits(i) for i in range(gen.leaf_spline_trait_count())]

        assert len(py_traits) == len(cpp_traits)
        for py_t, cpp_t in zip(py_traits, cpp_traits):
            for axis in ("x", "y", "z"):
                # Both pipelines build the stem from the same parameters; connection
                # points are simple stem positions and should agree within a few mm.
                assert py_t["connection_point"][axis] == pytest.approx(
                    cpp_t["connection_point"][axis], abs=0.01
                )
