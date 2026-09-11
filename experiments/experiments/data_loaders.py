"""Loaders for the external datasets used by the validation experiments.

Two data sources:

1. Jensina's phyllotaxy repo (`jdavis-132/phyllotaxy`) — population-scale
   trait CSVs derived from voxel reconstructions of 366 sorghum plants
   (panel of 2018) plus 10-plant manual ground-truth set (panel of 2024).

2. The Sorghum voxel-reconstruction archive (Zenodo 10.5281/zenodo.4426620)
   — per-plant raw voxels, segmented skeletons, gold-standard `(height, θ, φ)`
   per leaf, and camera calibrations for 351 plants from the 2018 panel.

All loaders return either pandas DataFrames or numpy arrays. Plant IDs in the
voxel set are full directory names like
`4-9-18_Schnable_49-281-JS39-65_2018-04-11_12-09-35_9968800`.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .paths import (
    JENSINA_DATA,
    VOXEL_DATASET,
    VOXEL_RECONSTRUCTED,
    VOXEL_SKELETONS,
)


VOXEL_LENGTHS = ("4cm", "6cm", "8cm", "12cm")


# ---------------------------------------------------------------------------
# Jensina's CSV loaders
# ---------------------------------------------------------------------------

def load_jensina_phi_csv(voxel_len: Optional[str] = "4cm") -> pd.DataFrame:
    """Load `processedData.csv`: per-plant phyllotaxic angles φᵢ for i=0..14.

    Each row is one (plant, timepoint, voxel-length) combination. Pass
    ``voxel_len=None`` to keep all four voxel-length variants in one frame.
    """
    df = pd.read_csv(JENSINA_DATA / "processedData.csv").dropna(how="all")
    if voxel_len is not None:
        df = df[df["voxel_len"] == voxel_len].copy()
    return df.reset_index(drop=True)


def load_jensina_theta_csv(voxel_len: Optional[str] = "4cm") -> pd.DataFrame:
    """Load `processedData_Thetas.csv`: per-plant leaf-insertion angles θᵢ for i=0..14."""
    df = pd.read_csv(JENSINA_DATA / "processedData_Thetas.csv").dropna(how="all")
    if voxel_len is not None:
        df = df[df["voxel_len"] == voxel_len].copy()
    return df.reset_index(drop=True)


def load_jensina_raw_angles(voxel_len: str = "4cm") -> pd.DataFrame:
    """Load `angles_three_days_{voxel_len}.csv`: raw per-leaf (height, θ, φ).

    Each row corresponds to one reconstruction (plant × timepoint). Up to 15
    leaves per plant, each with three columns: `height`, `theta`, `phi`.
    The header has two rows; we read the second.
    """
    if voxel_len not in VOXEL_LENGTHS:
        raise ValueError(f"voxel_len must be one of {VOXEL_LENGTHS}, got {voxel_len!r}")
    path = JENSINA_DATA / f"angles_three_days_{voxel_len}.csv"
    df = pd.read_csv(path, header=1, encoding="utf-8-sig")
    df = df.dropna(axis=1, how="all")
    return df


def load_jensina_manual() -> pd.DataFrame:
    """Load the 10-plant 2024 manual ground-truth set (`leafLevelTraits.csv`).

    Plant IDs here are simple integers (187..196) and are NOT in the
    Zenodo voxel set — those plants were grown in 2024 specifically for
    manual validation. Useful for measuring inter-rater agreement; not
    usable as direct ground truth against our voxel-fit descriptors.
    """
    return pd.read_csv(JENSINA_DATA / "leafLevelTraits.csv")


def load_jensina_genotypes() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (genotypes_df, JS-to-PI mapping)."""
    geno = pd.read_csv(JENSINA_DATA / "genotypes.csv")
    js_pi = pd.read_csv(JENSINA_DATA / "JSname_PIname.csv")
    return geno, js_pi


# ---------------------------------------------------------------------------
# Voxel-set loaders
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def list_voxel_plants(require_skeleton: bool = False) -> tuple[str, ...]:
    """Return the sorted list of plant directory names with voxel reconstructions.

    With ``require_skeleton=True``, restricts to plants that also have a
    successful skeleton fit (`optim_skeleton.txt` and `angles.txt`). About
    19/351 plants in the Zenodo set have voxels but no usable skeleton —
    excluding them is the right default for any experiment that depends on
    the per-leaf trait extraction.
    """
    if not VOXEL_RECONSTRUCTED.exists():
        return ()
    plants = sorted(p.name for p in VOXEL_RECONSTRUCTED.iterdir() if p.is_dir())
    if require_skeleton:
        plants = [
            p
            for p in plants
            if (VOXEL_SKELETONS / p / "optim_skeleton.txt").exists()
            and (VOXEL_SKELETONS / p / "angles.txt").exists()
        ]
    return tuple(plants)


def _read_voxel_indices(path: Path) -> np.ndarray:
    """Read a `.txt` file of voxel indices (first line = count, then `i j k` lines)."""
    with open(path) as f:
        n = int(f.readline().strip())
        # Use loadtxt for the rest; could also use np.fromfile but loadtxt handles whitespace cleanly.
        arr = np.loadtxt(f, dtype=int)
    if arr.ndim == 1:
        arr = arr.reshape(1, 3)
    if arr.shape[0] != n:
        raise ValueError(f"Expected {n} voxels in {path}, got {arr.shape[0]}")
    return arr


def load_voxel_grid(plant_id: str) -> np.ndarray:
    """Load `reconstructed/<plant_id>/voxels.txt` as an (N, 3) int array of 512³ indices."""
    return _read_voxel_indices(VOXEL_RECONSTRUCTED / plant_id / "voxels.txt")


def load_skeleton(plant_id: str) -> np.ndarray:
    """Load `skeletons/<plant_id>/optim_skeleton.txt` as an (N, 3) int array of 512³ indices."""
    return _read_voxel_indices(VOXEL_SKELETONS / plant_id / "optim_skeleton.txt")


def load_angles_txt(plant_id: str) -> pd.DataFrame:
    """Load gold-standard `(height_norm, θ, φ)` per leaf from `skeletons/<plant_id>/angles.txt`.

    Format: a single tab-separated line of `n_leaves * 3` floats interleaved
    as height, theta, phi. Returns one row per leaf.
    """
    path = VOXEL_SKELETONS / plant_id / "angles.txt"
    raw = path.read_text().strip().split("\t")
    raw = [float(x) for x in raw if x.strip()]
    if len(raw) % 3 != 0:
        raise ValueError(f"angles.txt for {plant_id} has {len(raw)} values; expected a multiple of 3")
    n_leaves = len(raw) // 3
    return pd.DataFrame(
        {
            "leaf_index": list(range(n_leaves)),
            "height_norm": raw[0::3][:n_leaves],
            "theta": raw[1::3][:n_leaves],
            "phi": raw[2::3][:n_leaves],
        }
    )


def load_reconstruction_error(plant_id: str) -> float:
    """Return the Dice coefficient between reprojected voxels and segmented input views.

    Returns NaN if the file is missing (a few plants have voxels but no error log).
    """
    path = VOXEL_RECONSTRUCTED / plant_id / "error.txt"
    if not path.exists():
        return float("nan")
    return float(path.read_text().strip())


def load_skeleton_error(plant_id: str) -> float:
    """Return the skeleton-fit error reported by the Gaillard pipeline.

    Returns NaN if the file is missing (about 19/351 plants in the Zenodo set
    have voxels but no successful skeleton fit).
    """
    path = VOXEL_SKELETONS / plant_id / "error.txt"
    if not path.exists():
        return float("nan")
    return float(path.read_text().strip())


def has_skeleton(plant_id: str) -> bool:
    """True iff the plant has a usable skeleton + per-leaf angle extraction."""
    return (VOXEL_SKELETONS / plant_id / "optim_skeleton.txt").exists() and (
        VOXEL_SKELETONS / plant_id / "angles.txt"
    ).exists()


# ---------------------------------------------------------------------------
# Plant-ID parsing and cross-referencing
# ---------------------------------------------------------------------------

_PLANT_ID_RE = re.compile(
    r"^\d+-\d+-\d+_Schnable_\d+-(?P<plant_num>\d+)-(?P<js_id>.+?)"
    r"_(?P<img_date>\d{4}-\d{2}-\d{2})_"
)


def parse_plant_id(plant_id: str) -> dict:
    """Parse metadata fields out of a voxel-set directory name.

    Returns a dict with `plant_num` (int), `js_id` (str), `img_date` (str)
    or all-None values if the format doesn't match.
    """
    m = _PLANT_ID_RE.match(plant_id)
    if not m:
        return {"plant_num": None, "js_id": None, "img_date": None}
    d = m.groupdict()
    d["plant_num"] = int(d["plant_num"])
    return d


def voxel_plant_metadata() -> pd.DataFrame:
    """Return a DataFrame summarizing all voxel-set plants: name, parsed metadata, errors, skeleton presence."""
    rows = []
    for plant_id in list_voxel_plants():
        meta = parse_plant_id(plant_id)
        rows.append(
            {
                "plant_id": plant_id,
                **meta,
                "recon_error": load_reconstruction_error(plant_id),
                "skel_error": load_skeleton_error(plant_id),
                "has_skeleton": has_skeleton(plant_id),
            }
        )
    return pd.DataFrame(rows)
