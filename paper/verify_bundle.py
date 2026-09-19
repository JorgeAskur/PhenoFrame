"""Check the tracked paper inputs against the study populations."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PAPER = Path(__file__).resolve().parent
SOURCE = PAPER / "source_data"
EXPECTED = {
    "sorghum/e1b_per_leaf_pymaize_vs_gold.csv": (2476, 324),
    "sorghum/e4_geometric_fidelity.csv": (324, 324),
    "sorghum/e5_fitted_params.csv": (2479, 324),
    "sorghum/e5_roundtrip_vs_gold.csv": (2469, 324),
    "maize/e1b_per_leaf_pymaize_vs_gold.csv": (24537, 3318),
    "maize/e4_geometric_fidelity.csv": (3317, 3317),
    "maize/e5_fitted_params.csv": (24646, 3317),
    "maize/e5_roundtrip_vs_gold.csv": (24453, 3317),
    "generation/sec2_5_synthetic_traits.csv": (2205, 324),
    "generation/sec2_5_holdout_split.csv": (324, 324),
}


def main() -> None:
    for relative, (rows, plants) in EXPECTED.items():
        data = pd.read_csv(SOURCE / relative)
        actual = (len(data), data["plant_id"].nunique())
        if actual != (rows, plants):
            raise AssertionError(f"{relative}: expected {(rows, plants)}, found {actual}")
        print(f"OK {relative}: {rows:,} rows; {plants:,} plant IDs")

    xmls = list((PAPER / "synthetic_descriptors").glob("plant_*.xml"))
    if len(xmls) != 500:
        raise AssertionError(f"Expected 500 synthetic descriptors, found {len(xmls)}")
    print("OK synthetic_descriptors: 500 XMLs")

    for relative, rows in (
        ("summary/R123_table3_definitive.csv", 78),
        ("summary/B1_table2_by_leafset.csv", 12),
    ):
        actual = len(pd.read_csv(SOURCE / relative))
        if actual != rows:
            raise AssertionError(f"{relative}: expected {rows} rows, found {actual}")
        print(f"OK {relative}: {rows} rows")

    print("Bundle checks passed. Raw image/voxel and GWAS inputs remain external.")


if __name__ == "__main__":
    main()
