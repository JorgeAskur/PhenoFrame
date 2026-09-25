"""Export a per-leaf phyllotaxis (phi azimuth, plus theta) CSV for the SORGHUM
voxel dataset, with the angle from three sources:

  1. voxel plants     -> gold standard from skeletons/<plant>/angles.txt (the
                         Gaillard/Mathieu PCA extraction on the voxel skeleton).
  2. descriptor spline -> spline-based descriptor traits via
                         extract_traits_via_phenoframe (the E1b path, notebook 03).
  3. fitted descriptor -> 3-parameter procedural descriptor via
                         skeleton_to_procedural_xml + compute_traits_from_descriptor
                         (the E5 round-trip path, notebook 05).

phi convention matches the notebooks exactly: phi = ((180 - azimuth) % 360) - 180,
with a per-plant +/-180 residual sign alignment against gold (resolving the v2
PCA eigenvector ambiguity), applied independently for the spline and fitted
sources. This is the clean Zenodo sorghum skeleton set, so NO clean_segment
patch (that fix is maize-only).
"""
import argparse
import os, sys, tempfile
from pathlib import Path
import numpy as np
import pandas as pd

from _data_paths import GENERATED, PAPER, REPO

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--voxel-root", type=Path, required=True)
parser.add_argument("--output", type=Path, default=GENERATED / "sorghum_phyllotaxis_angles.csv")
args = parser.parse_args()
os.environ["PHENOFRAME_VOXEL_PATH"] = str(args.voxel_root.resolve())
os.environ.setdefault("PHENOFRAME_JENSINA_PATH", str(PAPER / "reference" / "phyllotaxy"))
sys.path.insert(0, str(REPO / "experiments"))

from experiments import paths, data_loaders as dl
from phenoframe.skeleton_to_descriptor import (
    extract_traits_via_phenoframe, skeleton_to_procedural_xml)
from phenoframe.traits import compute_traits_from_descriptor

paths.assert_data_present()

N_LEAF_CTRL = 4        # spline control points (E1b)
SPLINE_POINTS = 8      # procedural forward-model control points (E5)
STEM_RADIUS = 0.015    # m, typical sorghum stem (E5)
OUT = args.output


def wrap180(x):
    return ((np.asarray(x, float) + 180) % 360) - 180


def procedural_traits(sk, tmp, pid):
    """E5 path: procedural descriptor -> per-leaf (theta, phi)."""
    xml = tmp / f"{pid}.xml"
    skeleton_to_procedural_xml(sk, xml, stem_radius=STEM_RADIUS, spline_points=SPLINE_POINTS)
    out = []
    for t in compute_traits_from_descriptor(xml):
        az = t["leaf_angle"]["azimuth_deg"]
        inc = t["leaf_angle"]["inclination_deg"]
        out.append((90.0 - inc, ((180.0 - az) % 360.0) - 180.0))
    return out  # list of (theta, phi) bottom-to-top


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    plants = dl.list_voxel_plants(require_skeleton=True)
    print(f"Processing {len(plants)} sorghum plants...")
    records, failures = [], []
    tmp = Path(tempfile.mkdtemp())
    for i, pid in enumerate(plants):
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(plants)}")
        try:
            gold = dl.load_angles_txt(pid)          # height_norm, theta, phi
            sk = dl.load_skeleton(pid)
            spline = extract_traits_via_phenoframe(sk, n_leaf_ctrl=N_LEAF_CTRL)  # theta, phi cols
            fitted = procedural_traits(sk, tmp, pid)
        except Exception as e:
            failures.append((pid, type(e).__name__, str(e)))
            continue
        n = min(len(gold), len(spline), len(fitted))
        for li in range(n):
            records.append({
                "plant_id": pid,
                "leaf_index": li,
                "height_norm": gold.iloc[li]["height_norm"],
                "theta_voxel": gold.iloc[li]["theta"],
                "phi_voxel": gold.iloc[li]["phi"],
                "theta_spline": spline.iloc[li]["theta"],
                "phi_spline": spline.iloc[li]["phi"],
                "theta_fitted": fitted[li][0],
                "phi_fitted": fitted[li][1],
            })

    df = pd.DataFrame(records)
    print(f"\n{len(df)} leaves across {df['plant_id'].nunique()} plants. "
          f"Failures: {len(failures)}.")

    # Per-plant +/-180 phi sign alignment against gold, independently for each
    # descriptor source (same resolution of v2 ambiguity as E1b/E5 notebooks).
    for src in ("phi_spline", "phi_fitted"):
        flip = set()
        for pid, grp in df.groupby("plant_id"):
            g, o = grp["phi_voxel"].values, grp[src].values
            if (wrap180(o + 180 - g) ** 2).sum() < (wrap180(o - g) ** 2).sum():
                flip.add(pid)
        m = df["plant_id"].isin(flip)
        df.loc[m, src] = wrap180(df.loc[m, src].values + 180)
        print(f"  {src}: +180 flip applied to {len(flip)} plants")

    # signed phi deltas for convenience
    df["dphi_spline"] = wrap180(df["phi_spline"].values - df["phi_voxel"].values)
    df["dphi_fitted"] = wrap180(df["phi_fitted"].values - df["phi_voxel"].values)

    df.to_csv(OUT, index=False)
    print(f"\nSaved {OUT}  ({len(df)} rows)")
    print(f"  median |dphi| spline vs voxel: {df['dphi_spline'].abs().median():.1f} deg")
    print(f"  median |dphi| fitted vs voxel: {df['dphi_fitted'].abs().median():.1f} deg")
    if failures:
        pd.DataFrame(failures, columns=["plant_id", "err", "msg"]).to_csv(
            OUT.with_name("sorghum_phyllotaxis_failures.csv"), index=False)
        print(f"  {len(failures)} failures logged")


if __name__ == "__main__":
    main()
