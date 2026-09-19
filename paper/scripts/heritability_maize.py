"""E8 broad-sense heritability, adapted from notebook 06 for the GIC 20-genotype
maize dataset.

Differences from the sorghum notebook:
  * Genotype is parsed directly from the maize plant_id (`_85-<num>-<GENO>-<rep>_`)
    instead of the JS->PI Jensina lookup. No reference genotype is dropped.
  * The GIC data scans each physical plant (genotype x rep) at many lifecycle
    timepoints. Treating each scan as a biological replicate would pseudo-
    replicate, so we collapse to ONE value per physical plant: the median of
    each per-scan trait across all scans that have >=5 leaves (a developed lower
    canopy). That yields 80 plants = 20 genotypes x 4 reps.

Everything else follows Davis 2025: lower-5 leaves, consecutive-leaf azimuth
deviation Phi_dev_i = |((phi_{i+1}-phi_i) mod 360) - 180|, and the lme4-equiv
random-intercept REML model with H2 (Eq. 3, n=2) = sigma2_G / (sigma2_G + sigma2_e/2).
"""
import re, warnings
from pathlib import Path
from _data_paths import GENERATED, MAIZE
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

INPUT = MAIZE
OUT = GENERATED / "legacy_heritability"
FIG = GENERATED / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
LOWER_LEAVES = list(range(5))
MIN_REPS = 2
MIN_LEAVES_PER_SCAN = 5     # a scan must have a developed lower canopy to count
SRC = ["gold", "e1b", "e5"]
_RX = re.compile(r"_\d+-(\d+)-([A-Za-z0-9]+)-(\d+)_(\d{4}-\d{2}-\d{2})_")


def norm360(x):
    return np.mod(x, 360.0)


def norm180(x):
    return np.abs(norm360(x) - 180.0)


def parse_geno(pid):
    m = _RX.search(pid)
    return (None, None) if not m else (m.group(2), f"{m.group(2)}-{m.group(3)}")


def per_scan_features(leaf):
    """One row per scan (plant_id) with Davis lower-5-leaf features for each src."""
    leaf = leaf[leaf["leaf_index"].isin(LOWER_LEAVES)].copy()
    for s in SRC:
        leaf[f"phi_{s}"] = norm360(leaf[f"phi_{s}"])
    wide = leaf.pivot_table(
        index="plant_id", columns="leaf_index",
        values=[f"phi_{s}" for s in SRC] + [f"theta_{s}" for s in SRC],
        aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    # require a developed lower canopy: all 5 lower theta_gold present
    need = [f"theta_gold_{i}" for i in LOWER_LEAVES]
    wide = wide[wide[need].notna().sum(axis=1) >= MIN_LEAVES_PER_SCAN].copy()
    for s in SRC:
        dev = []
        for i in range(1, 5):
            c = f"Phi_dev_{s}_{i}"
            wide[c] = norm180(wide[f"phi_{s}_{i}"] - wide[f"phi_{s}_{i-1}"])
            dev.append(c)
        wide[f"med_phi_dev_{s}"] = wide[dev].median(axis=1, skipna=True)
        wide[f"am_phi_dev_{s}"] = wide[dev].mean(axis=1, skipna=True)
        tcol = [f"theta_{s}_{i}" for i in LOWER_LEAVES]
        wide[f"med_theta_{s}"] = wide[tcol].median(axis=1, skipna=True)
        wide[f"am_theta_{s}"] = wide[tcol].mean(axis=1, skipna=True)
    return wide


def calc_heritability(df, trait, group="geno"):
    sub = df[[group, trait]].dropna().rename(columns={trait: "y"})
    n_per = sub.groupby(group).size()
    sub = sub[sub[group].isin(n_per[n_per > 1].index)].copy()
    out = {"trait": trait, "n_plants": len(sub), "n_genotypes": sub[group].nunique(),
           "sigma2_G": np.nan, "sigma2_e": np.nan, "H2": np.nan, "H2_n2": np.nan}
    if out["n_genotypes"] < 2:
        return out
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sg2 = se2 = None
        for method in ("bfgs", "powell", "lbfgs"):
            try:
                mf = smf.mixedlm("y ~ 1", sub, groups=sub[group]).fit(method=method, reml=True)
                if mf.converged:
                    sg2, se2 = float(mf.cov_re.iloc[0, 0]), float(mf.scale); break
            except Exception:
                continue
        if sg2 is None:
            return out
    out["sigma2_G"], out["sigma2_e"] = sg2, se2
    out["H2"] = sg2 / (sg2 + se2) if sg2 + se2 > 0 else 0.0
    out["H2_n2"] = sg2 / (sg2 + se2 / 2.0) if sg2 + se2 / 2.0 > 0 else 0.0
    return out


def trait_list_for(src):
    out = [("per-node theta", f"theta_{i}", f"theta_{src}_{i}") for i in LOWER_LEAVES]
    out += [("per-node Phi_dev", f"Phi_dev_{i}", f"Phi_dev_{src}_{i}") for i in range(1, 5)]
    out += [("per-plant theta", "med_theta", f"med_theta_{src}"),
            ("per-plant theta", "am_theta", f"am_theta_{src}"),
            ("per-plant Phi_dev", "med_phi_dev", f"med_phi_dev_{src}"),
            ("per-plant Phi_dev", "am_phi_dev", f"am_phi_dev_{src}")]
    return out


def main():
    e1b = pd.read_csv(INPUT / "e1b_per_leaf_pymaize_vs_gold.csv")
    e5 = pd.read_csv(INPUT / "e5_roundtrip_vs_gold.csv")
    leaf = (e1b[["plant_id", "leaf_index", "theta_gold", "phi_gold",
                 "theta_ours", "phi_ours"]]
            .rename(columns={"theta_ours": "theta_e1b", "phi_ours": "phi_e1b"})
            .merge(e5[["plant_id", "leaf_index", "theta_rt", "phi_rt"]]
                   .rename(columns={"theta_rt": "theta_e5", "phi_rt": "phi_e5"}),
                   on=["plant_id", "leaf_index"], how="inner"))
    print(f"merged per-leaf rows: {len(leaf):,} across {leaf.plant_id.nunique()} scans")

    scans = per_scan_features(leaf)
    scans["geno"], scans["physical"] = zip(*scans["plant_id"].map(parse_geno))
    scans = scans[scans["geno"].notna()].copy()
    print(f"scans with >={MIN_LEAVES_PER_SCAN} leaves: {len(scans)} "
          f"across {scans.physical.nunique()} physical plants, {scans.geno.nunique()} genotypes")

    # collapse lifecycle timepoints -> one row per physical plant (median trait)
    feat_cols = [c for c in scans.columns if c.startswith(("theta_", "phi_", "Phi_dev_",
                 "med_theta_", "am_theta_", "med_phi_dev_", "am_phi_dev_"))]
    plant = scans.groupby(["physical", "geno"])[feat_cols].median().reset_index()
    reps = plant.groupby("geno").size()
    keep = set(reps[reps >= MIN_REPS].index)
    plant = plant[plant.geno.isin(keep)].copy()
    print(f"per-physical-plant rows: {len(plant)} | genotypes(reps>={MIN_REPS}): {plant.geno.nunique()} "
          f"| reps/geno min/med/max = {reps.min()}/{int(reps.median())}/{reps.max()}")

    rows = []
    for src in SRC:
        for cat, label, col in trait_list_for(src):
            r = calc_heritability(plant, col)
            r.update(source=src, category=cat, label=label)
            rows.append(r)
    results = pd.DataFrame(rows)[["source", "category", "label", "trait",
                                  "n_plants", "n_genotypes", "sigma2_G", "sigma2_e", "H2", "H2_n2"]]
    wide_h2 = (results.pivot_table(index=["category", "label"], columns="source",
               values="H2_n2", aggfunc="first")[["gold", "e1b", "e5"]].reset_index())
    wide_h2["extractor_loss"] = wide_h2["gold"] - wide_h2["e1b"]
    wide_h2["compression_loss"] = wide_h2["e1b"] - wide_h2["e5"]

    results.to_csv(OUT / "e8_heritability_maize.csv", index=False)
    wide_h2.to_csv(OUT / "e8_h2_by_trait_maize.csv", index=False)

    # summary
    print("\n=== H2 (Eq.3, n=2) by source ===")
    for src in SRC:
        sub = results[results.source == src]
        print(f"  {src:4s}: median H2={sub.H2_n2.median():.3f}  mean={sub.H2_n2.mean():.3f}  "
              f">=0.20: {(sub.H2_n2>=0.20).sum()}/{len(sub)}  >=0.40: {(sub.H2_n2>=0.40).sum()}  "
              f"max={sub.H2_n2.max():.3f}")
    print("\n=== headline table (H2_n2) ===")
    print(wide_h2.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # figure: grouped bars per trait
    n_plants, n_geno = plant.physical.nunique(), plant.geno.nunique()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    panels = [("per-node theta", r"Per-node leaf angle $\theta_i$"),
              ("per-node Phi_dev", r"Per-node phyllotaxic deviation $\Phi_i$"),
              ("per-plant", "Per-plant aggregates")]
    for ax, (cm, title) in zip(axes, panels):
        sub = wide_h2[wide_h2.category.str.startswith(cm)].reset_index(drop=True)
        x = np.arange(len(sub)); w = 0.27
        ax.bar(x - w, sub.gold, w, label="Voxel", color="#999999", edgecolor="black", lw=0.5)
        ax.bar(x, sub.e1b, w, label="Skeleton", color="#4477aa", edgecolor="black", lw=0.5)
        ax.bar(x + w, sub.e5, w, label="Procedural", color="#cc6677", edgecolor="black", lw=0.5)
        ax.axhline(0.20, color="black", ls="--", lw=0.8, alpha=0.6)
        ax.set_xticks(x); ax.set_xticklabels(sub.label, rotation=20, ha="right", fontsize=9)
        ax.set_ylim(0, 1); ax.grid(axis="y", alpha=0.3); ax.set_title(title, fontsize=11)
        if ax is axes[0]:
            ax.set_ylabel(r"Broad-sense $H^2$ (Eq. 3, $n=2$)"); ax.legend(fontsize=9)
    fig.suptitle(f"Maize heritability across pipelines (n={n_plants} plants, "
                 f"{n_geno} genotypes, lower 5 leaves)", y=1.02, fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "e8_h2_by_trait_maize.png", dpi=150, bbox_inches="tight")
    print(f"\nsaved {OUT/'e8_heritability_maize.csv'}")
    print(f"saved {OUT/'e8_h2_by_trait_maize.csv'}")
    print(f"saved {FIG/'e8_h2_by_trait_maize.png'}")


if __name__ == "__main__":
    main()
