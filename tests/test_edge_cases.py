"""Edge-case and error-path tests for descriptor parsing and trait computation.

Covers:
- Missing XML file → FileNotFoundError
- Malformed XML → ParseError
- Structurally invalid descriptor (missing <species>) → ValueError
- Empty descriptor (zero tillers / zero leaves)
- Single-leaf plant
- splinePoints clamped to the [4, 40] range
- Control-point override path used when useCtrlOverrides=1
"""

from __future__ import annotations

import warnings
from pathlib import Path
from xml.etree.ElementTree import ParseError

import pytest

from phenoframe import compute_traits_from_descriptor
from phenoframe.traits import _parse_descriptor


def _write(xml_path: Path, body: str) -> Path:
    xml_path.write_text(body, encoding="utf-8")
    return xml_path


class TestMissingAndMalformed:
    def test_missing_file_raises_file_not_found(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.xml"
        with pytest.raises(FileNotFoundError):
            compute_traits_from_descriptor(missing)

    def test_malformed_xml_raises_parse_error(self, tmp_path: Path) -> None:
        bad = _write(tmp_path / "bad.xml", "<plant><Tiller><leaves><leaf></plant>")
        with pytest.raises(ParseError):
            compute_traits_from_descriptor(bad)

    def test_missing_tiller_elements_raises_value_error(self, tmp_path: Path) -> None:
        no_tillers = _write(
            tmp_path / "no_tillers.xml",
            "<?xml version='1.0'?><plant><not_species/></plant>",
        )
        with pytest.raises(ValueError, match="no <Tiller> elements found"):
            compute_traits_from_descriptor(no_tillers)


class TestEmptyDescriptors:
    def test_zero_tillers_in_legacy_species_raises(self, tmp_path: Path) -> None:
        empty = _write(
            tmp_path / "empty.xml",
            "<?xml version='1.0'?><plant><species name='Maize_Procedural'/></plant>",
        )
        with pytest.raises(ValueError, match="no <Tiller> elements found"):
            compute_traits_from_descriptor(empty)

    def test_tiller_with_no_leaves_yields_no_traits(self, tmp_path: Path) -> None:
        no_leaves = _write(
            tmp_path / "no_leaves.xml",
            """<?xml version='1.0'?>
            <plant><species name='Maize_Procedural'>
              <Tiller type='main' radius='0.02' alpha='0' beta='180' randomSeed='1' stemShrink='0'/>
            </species></plant>""",
        )
        assert compute_traits_from_descriptor(no_leaves) == []


class TestSingleLeafPlant:
    def test_single_leaf_plant_produces_one_trait_record(self, tmp_path: Path) -> None:
        xml = _write(
            tmp_path / "single_leaf.xml",
            """<?xml version='1.0'?>
            <plant><species name='Maize_Procedural'>
              <Tiller type='main' radius='0.02' alpha='0' beta='180' randomSeed='7' stemShrink='0'>
                <leaves number='1'>
                  <leaf id='0' distance='0.20' leafLength='0.45' leafWidth='0.07'
                        leafAngle='35' droopiness='-10' stemInclinationDeg='0' splinePoints='8'/>
                </leaves>
              </Tiller>
            </species></plant>""",
        )
        traits = compute_traits_from_descriptor(xml)
        assert len(traits) == 1
        t = traits[0]
        assert t["leaf_id"] == 0
        assert t["tiller_index"] == 0
        assert t["leaf_length"] > 0.0


class TestSplinePointsClamping:
    @pytest.mark.parametrize(
        ("requested", "expected", "should_warn"),
        [
            (0, 4, True),    # below the minimum
            (1, 4, True),
            (3, 4, True),
            (4, 4, False),   # boundary: in-range
            (10, 10, False),
            (40, 40, False), # boundary: in-range
            (50, 40, True),  # above the maximum
            (1000, 40, True),
        ],
    )
    def test_spline_points_clamped_to_valid_range(
        self, tmp_path: Path, requested: int, expected: int, should_warn: bool
    ) -> None:
        xml = _write(
            tmp_path / f"spline_{requested}.xml",
            f"""<?xml version='1.0'?>
            <plant><species name='Maize_Procedural'>
              <Tiller type='main' radius='0.02' alpha='0' beta='180' randomSeed='1' stemShrink='0'>
                <leaves number='1'>
                  <leaf id='0' distance='0.2' leafLength='0.5' leafWidth='0.05'
                        leafAngle='30' droopiness='0' stemInclinationDeg='0'
                        splinePoints='{requested}'/>
                </leaves>
              </Tiller>
            </species></plant>""",
        )

        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            tillers = _parse_descriptor(xml)

        assert tillers[0].leaves[0].spline_points == expected

        clamp_warnings = [w for w in captured if "splinePoints" in str(w.message)]
        if should_warn:
            assert len(clamp_warnings) == 1
            assert str(requested) in str(clamp_warnings[0].message)
            assert str(expected) in str(clamp_warnings[0].message)
        else:
            assert clamp_warnings == []


class TestCtrlOverridePath:
    def test_ctrl_override_replaces_procedural_spline(self, tmp_path: Path) -> None:
        """With useCtrlOverrides=1 and an explicit ctrlCenter, leafLength/leafAngle/droopiness are ignored
        and the supplied control points define the entire blade."""
        # Local ctrl points along radial (+Y in local) → horizontal leaf in world.
        ctrl_pts = [(0.0, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.6, 0.0), (0.0, 0.9, 0.0)]
        attrs = " ".join(f"p{i}x='{p[0]}' p{i}y='{p[1]}' p{i}z='{p[2]}'" for i, p in enumerate(ctrl_pts))
        xml = _write(
            tmp_path / "override.xml",
            f"""<?xml version='1.0'?>
            <plant><species name='Maize_Procedural'>
              <Tiller type='main' radius='0' alpha='0' beta='0' randomSeed='0' stemShrink='0'>
                <leaves number='1'>
                  <leaf id='0' distance='0' leafLength='9.99' leafWidth='0.05'
                        leafAngle='80' droopiness='-50' stemInclinationDeg='0' splinePoints='4'
                        useCtrlOverrides='1'>
                    <ctrlCenter count='4' {attrs}/>
                  </leaf>
                </leaves>
              </Tiller>
            </species></plant>""",
        )

        # Override path is taken: descriptor flag set and ctrl_center_override populated.
        tillers = _parse_descriptor(xml)
        leaf = tillers[0].leaves[0]
        assert leaf.use_ctrl_overrides is True
        assert len(leaf.ctrl_center_override) == 4

        # The trait pipeline ignores the bogus leafLength=9.99 and uses the override polyline.
        traits = compute_traits_from_descriptor(xml)
        assert traits[0]["leaf_length"] == pytest.approx(0.9, abs=1e-9)
        assert traits[0]["leaf_angle"]["inclination_deg"] == pytest.approx(0.0, abs=1e-6)

class TestSchemaShippedWithPackage:
    """The descriptor XSD is part of the distribution and should be parseable."""

    def test_xsd_ships_with_package_and_is_well_formed(self) -> None:
        import xml.etree.ElementTree as ET

        import phenoframe as pkg

        xsd_path = Path(pkg.__file__).resolve().parent / "schemas" / "descriptor.xsd"
        assert xsd_path.exists(), f"XSD missing from package: {xsd_path}"

        # Well-formed XML and rooted at the XSD schema element.
        tree = ET.parse(xsd_path)
        root = tree.getroot()
        assert root.tag.endswith("schema")


class TestCtrlOverridePathExtra:
    def test_override_flag_without_points_falls_back_to_procedural(self, tmp_path: Path) -> None:
        """useCtrlOverrides=1 but missing <ctrlCenter> should not flip the flag on."""
        xml = _write(
            tmp_path / "override_no_points.xml",
            """<?xml version='1.0'?>
            <plant><species name='Maize_Procedural'>
              <Tiller type='main' radius='0.02' alpha='0' beta='180' randomSeed='1' stemShrink='0'>
                <leaves number='1'>
                  <leaf id='0' distance='0.2' leafLength='0.5' leafWidth='0.05'
                        leafAngle='30' droopiness='0' stemInclinationDeg='0'
                        splinePoints='8' useCtrlOverrides='1'/>
                </leaves>
              </Tiller>
            </species></plant>""",
        )
        tillers = _parse_descriptor(xml)
        leaf = tillers[0].leaves[0]
        assert leaf.use_ctrl_overrides is False
        assert leaf.ctrl_center_override == []
        # And computation still succeeds via the procedural path.
        assert compute_traits_from_descriptor(xml)[0]["leaf_length"] > 0.0
