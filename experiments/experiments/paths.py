"""Filesystem paths to external datasets.

Both datasets live outside the repo (we do not redistribute them). Override
the defaults via environment variables when running on a different machine:

    PHENOFRAME_JENSINA_PATH  -> root of jdavis-132/phyllotaxy clone
    PHENOFRAME_VOXEL_PATH    -> root of the Sorghum voxel-reconstruction archive
                              (the directory containing dataset/, reconstructed/, skeletons/)
"""

from __future__ import annotations

import os
from pathlib import Path


# Default fallbacks point at a sibling `datasets/` folder; override with the
# environment variables below (recommended) to use datasets stored elsewhere.
_DATA_ROOT = Path(os.environ.get("PHENOFRAME_DATA_ROOT", os.environ.get("PYMAIZE_DATA_ROOT", Path.home() / "phenoframe_datasets")))

JENSINA_REPO = Path(
    os.environ.get(
        "PHENOFRAME_JENSINA_PATH",
        os.environ.get("PYMAIZE_JENSINA_PATH", str(_DATA_ROOT / "phyllotaxy")),
    )
)
JENSINA_DATA = JENSINA_REPO / "Data"

VOXEL_ROOT = Path(
    os.environ.get(
        "PHENOFRAME_VOXEL_PATH",
        os.environ.get("PYMAIZE_VOXEL_PATH", str(_DATA_ROOT / "Sorghum")),
    )
)
VOXEL_DATASET = VOXEL_ROOT / "dataset"
VOXEL_RECONSTRUCTED = VOXEL_ROOT / "reconstructed"
VOXEL_SKELETONS = VOXEL_ROOT / "skeletons"


def assert_data_present() -> None:
    """Raise a clear error if the configured dataset paths don't exist."""
    missing = []
    for label, p in [
        ("JENSINA_DATA", JENSINA_DATA),
        ("VOXEL_RECONSTRUCTED", VOXEL_RECONSTRUCTED),
        ("VOXEL_SKELETONS", VOXEL_SKELETONS),
    ]:
        if not p.exists():
            missing.append(f"  {label} -> {p}")
    if missing:
        raise FileNotFoundError(
            "External datasets are missing on this machine. Either:\n"
            "  - install them at the default location, or\n"
            "  - set PHENOFRAME_JENSINA_PATH / PHENOFRAME_VOXEL_PATH env vars.\n\n"
            "Missing paths:\n" + "\n".join(missing)
        )
