"""Smoke tests for the C++ procedural model bridge.

These exercise the ctypes wrapper end-to-end: load the shared library, build
a minimal plant, write outputs to disk, and verify the basic invariants hold
(non-empty OBJ file, correct tiller / leaf counts, spline traits returned).
Tests skip cleanly when the compiled library is unavailable.
"""

from __future__ import annotations

import ctypes
from pathlib import Path

import pytest

from phenosuite import Leaf, Maize, Tiller
from phenosuite.wrapper import EXPECTED_C_API_VERSION, LeafDesc, LeafSplineTraits, TillerDesc

from .conftest import requires_wrapper


pytestmark = [requires_wrapper, pytest.mark.wrapper]


def _build_one_leaf_plant(gen: Maize) -> int:
    gen.reset()
    gen.set_species("maize")
    gen.set_global_settings(
        density_v=60.0,
        density_u=60.0,
        stem_density_scale=0.15,
        stem_row_override=0,
        stem_ribbon_spacing=0.15,
    )
    tiller_index = gen.add_tiller(Tiller().to_ctype())
    gen.add_leaf(
        tiller_index,
        Leaf(id=0, distance=0.2, leaf_length=0.4, leaf_width=0.08, leaf_angle=40.0).to_ctype(),
    )
    return tiller_index


class TestWrapperBasics:
    def test_create_and_destroy(self) -> None:
        with Maize() as gen:
            assert gen.tiller_count() == 0

    def test_set_species_returns_true(self) -> None:
        with Maize() as gen:
            gen.reset()
            assert gen.set_species("maize") is True

    def test_add_tiller_and_leaf_counts(self) -> None:
        with Maize() as gen:
            t = _build_one_leaf_plant(gen)
            assert gen.tiller_count() == 1
            assert gen.leaf_count(t) == 1

    def test_rebuild_produces_triangles(self) -> None:
        with Maize() as gen:
            _build_one_leaf_plant(gen)
            assert gen.rebuild() is True
            assert gen.triangle_count() > 0

    def test_leaf_spline_traits_after_rebuild(self) -> None:
        with Maize() as gen:
            _build_one_leaf_plant(gen)
            gen.rebuild()
            assert gen.leaf_spline_trait_count() == 1
            traits = gen.get_leaf_spline_traits(0)
            # Sanity: tip is somewhere away from connection along blade.
            cp = traits["connection_point"]
            tip = traits["tip_position"]
            assert (tip["x"] - cp["x"]) ** 2 + (tip["y"] - cp["y"]) ** 2 + (tip["z"] - cp["z"]) ** 2 > 0.0
            assert traits["leaf_length"] > 0.0


class TestObjExport:
    def test_save_obj_writes_non_empty_file(self, tmp_path: Path) -> None:
        with Maize() as gen:
            _build_one_leaf_plant(gen)
            gen.rebuild()
            out = tmp_path / "plant"
            # Pass empty texture paths so we don't depend on bundled textures.
            assert gen.save_obj(str(out), leaf_texture_path="", stem_texture_path="") is True

        obj_file = tmp_path / "plant.obj"
        assert obj_file.exists()
        assert obj_file.stat().st_size > 0
        text = obj_file.read_text(encoding="utf-8", errors="ignore")
        # An OBJ for a real mesh always contains vertex (`v`) and face (`f`) records.
        assert "\nv " in text or text.startswith("v ")
        assert "\nf " in text or text.startswith("f ")

class TestApiVersionCheck:
    def test_loaded_library_reports_expected_major_version(self) -> None:
        """The shipped DLL should match the wrapper's expected major version."""
        with Maize() as gen:
            major = gen._lib.maize_c_api_version_major()
            assert major == EXPECTED_C_API_VERSION[0]

    def test_struct_sizes_match_python_ctypes(self) -> None:
        """C-side and Python ctypes struct sizes must agree, byte-for-byte."""
        with Maize() as gen:
            assert gen._lib.maize_leaf_desc_size() == ctypes.sizeof(LeafDesc)
            assert gen._lib.maize_tiller_desc_size() == ctypes.sizeof(TillerDesc)
            assert gen._lib.maize_leaf_spline_traits_size() == ctypes.sizeof(LeafSplineTraits)


class TestObjExportRoundTrip:
    def test_save_xml_round_trip_loads(self, tmp_path: Path) -> None:
        with Maize() as gen:
            _build_one_leaf_plant(gen)
            xml_out = tmp_path / "saved.xml"
            assert gen.save_xml(str(xml_out)) is True
            assert xml_out.exists()
            assert xml_out.stat().st_size > 0

        with Maize() as reloaded:
            assert reloaded.load_xml(str(xml_out)) is True
            assert reloaded.tiller_count() == 1
            assert reloaded.leaf_count(0) == 1
