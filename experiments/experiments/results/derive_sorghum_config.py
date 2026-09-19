"""Derive sorghum_generator_config.xml from E5 fitted parameter statistics.

The C++ MaizeGenerator (`MaizeProceduralModel/MaizeGenerator.cpp`) consumes
an XML config of the form described in §4.2.2 of the paper:

  <PerParameter mean= std= />   for 14 plant parameters
  <Curve bottom=  middle=  top= />  for per-rank canopy modulation
  <traitLoadingsMatrix>         4-latent loadings (Vigor/Slenderness/Posture/Texture)

This script reads `outputs/e5_fitted_params.csv`, fits the means/stds and
canopy curves for the architectural parameters PhenoFrame can fit
(leafLength, leafAngle, droopiness, internode distance -> plant height),
and writes a sorghum-flavored config. Visual / texture / non-fit
architectural parameters (leafTwist, leafCurl, surface noise) inherit the
maize-config values so the rendering still looks like a grass plant.

The trait loadings matrix (the 4-latent structure) is preserved from the
maize config -- those loadings encode biology (vigor -> large leaves +
upright) that is roughly species-invariant.

Run:
    python experiments/experiments/results/derive_sorghum_config.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from phenoframe import MaizeGenerator

# Outputs from E1b/E5 live alongside this script (the results/ notebook dir).
OUT = Path(__file__).resolve().parent / "outputs"

# Source maize config to inherit visual/texture defaults + loadings matrix.
MAIZE_CFG = REPO / "paper" / "configs" / "maize_generator_config.xml"
SORGHUM_CFG_OUT = REPO / "paper" / "generated" / "sorghum_generator_config.xml"

# ---------------------------------------------------------------------------
# 1. Load E5 fitted parameters (same filter as E5 / E9)
# ---------------------------------------------------------------------------
raw = pd.read_csv(OUT / "e5_fitted_params.csv")
rms_cap = raw["rms_error_cm"].quantile(0.90)
fits = raw[(raw["leaf_voxel_count"] >= 20) & (raw["rms_error_cm"] <= rms_cap)].copy()
print(f"E5 filtered population: {len(fits):,} leaves across {fits['plant_id'].nunique()} plants")

# ---------------------------------------------------------------------------
# 2. Compute per-plant aggregates we need for the config
# ---------------------------------------------------------------------------
per_plant = fits.groupby("plant_id").agg(
    n_leaves=("leaf_index", "count"),
    plant_height=("distance", "sum"),
)
n_leaves_mean = float(per_plant["n_leaves"].mean())
n_leaves_std = float(per_plant["n_leaves"].std())
height_mean = float(per_plant["plant_height"].mean())
height_std = float(per_plant["plant_height"].std())
print(f"Plants: n_leaves = {n_leaves_mean:.2f} +/- {n_leaves_std:.2f}")
print(f"        height   = {height_mean:.3f} +/- {height_std:.3f} m")

# ---------------------------------------------------------------------------
# 3. Per-node canopy curves
#
# MaizeGenerator samples (base_value, std) for each parameter and then
# modulates it across the canopy via either a multiplicative scale curve
# (leafLength, leafWidth, widthTaper, waveAmp, waveFreq) or an additive
# offset curve (leafAngle, droopiness, leafTwist, leafCurl,
# stemInclination). Curves interpolate among bottom (t=0), middle (t=0.5),
# top (t=1) values where t = leaf_index / (n_leaves - 1).
#
# We compute per-rank means on a normalized t axis, then sample at t=0,
# 0.5, 1.0 to extract bottom/middle/top.
# ---------------------------------------------------------------------------
# Map each leaf to its normalized canopy position t
fits = fits.merge(per_plant["n_leaves"].rename("plant_n_leaves"),
                  left_on="plant_id", right_index=True)
fits["t"] = fits["leaf_index"] / np.maximum(fits["plant_n_leaves"] - 1, 1)

def per_t_curve(values, param, n_bins=11):
    """Bin by t and return mean (bin centers, bin means)."""
    bins = np.linspace(0, 1, n_bins + 1)
    fits["_bin"] = pd.cut(fits["t"], bins, include_lowest=True, labels=False)
    return fits.groupby("_bin")[param].mean().reindex(range(n_bins))

def sample_curve_at(values, t_query):
    """Linear-interpolate the binned curve at a target t in [0, 1]."""
    centers = (np.arange(len(values)) + 0.5) / len(values)
    return float(np.interp(t_query, centers, values))

def fit_curve(param, multiplicative=True):
    """Return (base_mean, base_std, bottom, middle, top).

    For multiplicative curves: base = empirical middle value; scale = empirical/middle.
    For additive curves:      base = empirical middle value; offset = empirical - middle.
    """
    binned = per_t_curve(fits, param).values
    bottom = sample_curve_at(binned, 0.0)
    middle = sample_curve_at(binned, 0.5)
    top    = sample_curve_at(binned, 1.0)
    # Use pooled std across all leaves -- representative of within-plant variation
    base_std = float(fits[param].std())
    if multiplicative:
        if middle == 0:
            return middle, base_std, 1.0, 1.0, 1.0
        return middle, base_std, bottom / middle, 1.0, top / middle
    else:
        return middle, base_std, bottom - middle, 0.0, top - middle

leaf_length_base, leaf_length_std, ll_bot, ll_mid, ll_top = fit_curve("leaf_length", multiplicative=True)
droopiness_base,  droopiness_std,  dr_bot, dr_mid, dr_top = fit_curve("droopiness",  multiplicative=False)

# ---------------------------------------------------------------------------
# leafAngle from gold-standard theta (NOT E5 procedural-fit angles)
#
# Using E5 leaf_angle directly inherits the procedural-fit bias documented
# in §2.3 (median |d theta| ~ 5 deg between procedural fit and gold). Pull
# the calibration target from angles.txt instead -- it is the same source
# Davis 2025 uses for its theta. The descriptor leafAngle that produces a
# measured theta of theta_gold is approximately 90 - theta_gold, with a
# small droopiness correction (~2-3 deg) that we omit since the generator's
# droopiness is sampled, not fixed.
# ---------------------------------------------------------------------------
e1b_for_angle = pd.read_csv(OUT / "e1b_per_leaf_pymaize_vs_gold.csv")
gold = e1b_for_angle[["plant_id", "leaf_index", "theta_gold"]].dropna().copy()
gold = gold.merge(per_plant["n_leaves"].rename("plant_n_leaves"),
                  left_on="plant_id", right_index=True)
gold["t"] = gold["leaf_index"] / np.maximum(gold["plant_n_leaves"] - 1, 1)
# Convert gold theta to descriptor leafAngle: leafAngle = 90 - theta
gold["leafAngle_descriptor"] = 90.0 - gold["theta_gold"]

def per_t_curve_on(df, param, n_bins=11):
    bins = np.linspace(0, 1, n_bins + 1)
    df["_bin"] = pd.cut(df["t"], bins, include_lowest=True, labels=False)
    return df.groupby("_bin")[param].mean().reindex(range(n_bins))

binned_la = per_t_curve_on(gold, "leafAngle_descriptor").values
la_bot_gold = sample_curve_at(binned_la, 0.0)
la_mid_gold = sample_curve_at(binned_la, 0.5)
la_top_gold = sample_curve_at(binned_la, 1.0)
leaf_angle_base = la_mid_gold
leaf_angle_std  = float(gold["leafAngle_descriptor"].std())
la_bot = la_bot_gold - la_mid_gold
la_top = la_top_gold - la_mid_gold
print(f"\nleafAngle calibrated from gold-standard theta (not E5 fits):")
print(f"  bottom/middle/top descriptor leafAngle = {la_bot_gold:.1f}/{la_mid_gold:.1f}/{la_top_gold:.1f} deg")

print(f"\nleafLengthBase = {leaf_length_base:.4f} +/- {leaf_length_std:.4f}  "
      f"scale curve bottom/mid/top = {ll_bot:.3f}/1.000/{ll_top:.3f}")
print(f"leafAngleBase  = {leaf_angle_base:.2f} +/- {leaf_angle_std:.2f}  "
      f"offset curve bottom/mid/top = {la_bot:+.2f}/0.00/{la_top:+.2f}")
print(f"droopinessBase = {droopiness_base:.2f} +/- {droopiness_std:.2f}  "
      f"offset curve bottom/mid/top = {dr_bot:+.2f}/0.00/{dr_top:+.2f}")

# ---------------------------------------------------------------------------
# 4. Auto-calibrate tillerAzimuthNoise from real |Delta phi - 180| median
# ---------------------------------------------------------------------------
e1b = pd.read_csv(OUT / "e1b_per_leaf_pymaize_vs_gold.csv")
_p = e1b[["plant_id", "leaf_index", "phi_gold"]].dropna().copy()
_p = _p.sort_values(["plant_id", "leaf_index"])
_p["phi_norm"] = _p["phi_gold"] % 360.0
_p["phi_next"] = _p.groupby("plant_id")["phi_norm"].shift(-1)
_diff_dev = ((_p["phi_next"] - _p["phi_norm"]) % 360.0 - 180.0).abs().dropna()
# MaizeGenerator code: noise = hash_to_float_minus1_1(seed) * tillerAzimuthNoise
# hash_to_float_minus1_1 returns uniform[-1, 1], so noise ~ U(-A, +A) with A = tillerAzimuthNoise.
# phi_diff_dev = |noise_{i+1} - noise_i|, where each noise is U(-A, A).
# noise_diff = noise_{i+1} - noise_i has triangular distribution on [-2A, 2A]
# with std = A * sqrt(2/3). median(|triangular|) = 2A * (1 - 1/sqrt(2)) ~ 0.586 * A.
# Empirically calibrate amplitude A from observed median:
tiller_azimuth_noise = float(_diff_dev.median() / 0.5858)
print(f"\nempirical median |Phi_diff - 180 deg| = {float(_diff_dev.median()):.1f}  ->  "
      f"tillerAzimuthNoise (uniform amplitude) = {tiller_azimuth_noise:.1f}")

# ---------------------------------------------------------------------------
# 5. Whorl-compression parameters from distances
# ---------------------------------------------------------------------------
# Compare lower-half mean vs upper-half mean of distance per plant to get
# the compression factor. Pick whorlStart at 0.55 (the t above which the
# means drop sharply in sorghum).
binned_dist = per_t_curve(fits, "distance").values
lower_mean = float(np.nanmean(binned_dist[:int(len(binned_dist) * 0.55)]))
upper_mean = float(np.nanmean(binned_dist[int(len(binned_dist) * 0.55):]))
whorl_min_weight = float(np.clip(upper_mean / max(lower_mean, 1e-6), 0.10, 0.95))
print(f"Lower-half mean distance = {lower_mean:.3f} m  Upper-half mean = {upper_mean:.3f} m")
print(f"=> whorlMinInternodeWeight = {whorl_min_weight:.3f}")

# ---------------------------------------------------------------------------
# 6. Hand the canopy specs to phenoframe.MaizeGenerator to write the XML.
#    The wrapper takes care of the element-name conventions
#    (leafLengthBase -> leafLengthScaleCurve, leafAngleBase ->
#    leafAngleOffsetCurve, droopinessBase -> droopinessOffsetCurve)
#    plus pretty-print indentation. We inherit the trait-loadings matrix
#    and visual defaults from the maize template unchanged.
# ---------------------------------------------------------------------------
canopy_params = {
    # Fitted from E5 + gold
    "leafLengthBase":  dict(mean=leaf_length_base, std=leaf_length_std,
                            bottom_scale=ll_bot,   middle_scale=1.0,   top_scale=ll_top),
    "leafAngleBase":   dict(mean=leaf_angle_base,  std=leaf_angle_std,
                            bottom_offset=la_bot,  middle_offset=0.0,  top_offset=la_top),
    "droopinessBase":  dict(mean=droopiness_base,  std=droopiness_std,
                            bottom_offset=dr_bot,  middle_offset=0.0,  top_offset=dr_top),
    # Plant-level Gaussians (no curve)
    "height":          dict(mean=height_mean,     std=height_std),
    "numLeaves":       dict(mean=n_leaves_mean,   std=n_leaves_std),
    # Sorghum-typical defaults for non-E5 knobs
    "stemRadius":      dict(mean=0.012,           std=0.003),
    "leafWidthBase":   dict(mean=0.040,           std=0.012),     # ~4 cm blade
    "widthTaper":      dict(mean=2.2,             std=0.20),      # more pointed than maize
    "leafCurl":        dict(mean=0.05,            std=0.04),
    "leafTwist":       dict(mean=0.0,             std=25.0),
}

# Per-leaf jitter scale: 0.60 is the empirical compromise between the maize
# default (0.687, which over-disperses theta) and tighter jitter (0.50,
# which under-disperses the heavy phi tail by reducing near-vertical leaf
# flips). 0.60 keeps the gold-calibrated theta tight without collapsing the
# phi tail.

written = MaizeGenerator().derive_config(
    canopy_params=canopy_params,
    base_config_path=MAIZE_CFG,
    output_config_path=SORGHUM_CFG_OUT,
    tiller_azimuth_noise=tiller_azimuth_noise,
    leaf_jitter_scale=0.60,
    extra_scalars={"whorlMinInternodeWeight": whorl_min_weight},
    root_attrs={"speciesName": "Sorghum", "phenotypeId": "SAP_E5"},
)
print(f"\nWrote sorghum config: {written}")
print(f"  inherits trait-loadings matrix + visual defaults from maize config")
