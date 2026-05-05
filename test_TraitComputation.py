from __future__ import annotations

import argparse
from pathlib import Path

from Traits.trait_lib import compute_traits_from_descriptor, write_traits_xml


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute spline-based maize leaf traits from descriptor XML.")
    parser.add_argument("--xml", default="plant_0.xml", help="Path to descriptor XML")
    parser.add_argument("--out-xml", default="traits.xml", help="Output traits XML path")
    args = parser.parse_args()

    xml_path = Path(args.xml)
    if not xml_path.exists():
        fallback = Path("plants") / args.xml
        if fallback.exists():
            xml_path = fallback
        else:
            raise FileNotFoundError(f"XML file not found: {args.xml}")

    traits = compute_traits_from_descriptor(xml_path)
    out_xml = Path(args.out_xml)
    write_traits_xml(out_xml, xml_path.resolve(), traits)

    print(f"Computed traits for {len(traits)} leaves.")
    print(f"Saved XML: {out_xml}")
    for t in traits:
        ang = t["leaf_angle"]
        tip = t["tip_position"]
        print(
            f"Leaf {t['leaf_index']:02d}: "
            f"len={t['leaf_length']:.4f}, "
            f"az={ang['azimuth_deg']:.2f} deg, "
            f"inc={ang['inclination_deg']:.2f} deg, "
            f"tip=({tip['x']:.4f}, {tip['y']:.4f}, {tip['z']:.4f})"
        )


if __name__ == "__main__":
    main()
