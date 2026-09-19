"""Confidence intervals for the PhenoFrame paper (issue #5).

Part A  R^2 bootstrap over plants for Tables 1 & 2 (theta/phi, S and P vs voxel V).
Part B  Chamfer/Hausdorff median + IQR for Table 2.
Part C  H^2 bootstrap over genotypes for Table 3, plus the fraction of resamples
        where procedural median H^2 beats voxel.

The H^2 point estimator in the paper is REML (statsmodels MixedLM). For 2,000
bootstrap resamples that is too slow, so the bootstrap uses a one-way random-
effects method-of-moments (ANOVA) estimator, validated here to reproduce the
published REML medians before it is trusted for the CI.
"""
import re, sys, os, warnings
from pathlib import Path
from _data_paths import GENERATED, MAIZE, PAPER, REPO, SORGHUM
import numpy as np
import pandas as pd

warnings.simplefilter("ignore")
os.environ.setdefault("PHENOFRAME_JENSINA_PATH", str(PAPER / "reference" / "phyllotaxy"))
RNG = np.random.default_rng(20260826)
NB = 2000
MZ = MAIZE
PY = SORGHUM
SOR_BK = SORGHUM
SOR_E4 = SOR_BK / "e4_geometric_fidelity.csv"
LOWER5 = list(range(5))
SRC = ["gold", "e1b", "e5"]          # voxel, skeleton(spline), procedural


def norm360(x): return np.mod(x, 360.0)
def norm180(x): return np.abs(norm360(x) - 180.0)


# ---------------------------------------------------------------- MoM H^2
def mom_h2(y, g):
    """One-way random-effects method-of-moments H^2 (Davis Eq.3, n=2).
    y: values, g: integer group codes. Returns H2_n2."""
    m = np.isfinite(y)
    y, g = y[m], g[m]
    if y.size < 3:
        return np.nan
    groups, gidx = np.unique(g, return_inverse=True)
    k = groups.size
    if k < 2:
        return np.nan
    ni = np.bincount(gidx)
    keep = ni >= 1
    N = y.size
    grand = y.mean()
    gsum = np.bincount(gidx, weights=y)
    gmean = gsum / ni
    ssb = np.sum(ni * (gmean - grand) ** 2)
    ssw = np.sum((y - gmean[gidx]) ** 2)
    dfb, dfw = k - 1, N - k
    if dfw <= 0 or dfb <= 0:
        return np.nan
    msb, msw = ssb / dfb, ssw / dfw
    n0 = (N - np.sum(ni ** 2) / N) / (k - 1)
    sg2 = max((msb - msw) / n0, 0.0)
    se2 = msw
    den = sg2 + se2 / 2.0
    return sg2 / den if den > 0 else 0.0


# ---------------------------------------------------------------- plant tables
def wide_from_leaf(leaf, group_col):
    """leaf has plant_id, group_col, leaf_index, theta_/phi_{src}. -> per-plant wide."""
    leaf = leaf[leaf.leaf_index.isin(LOWER5)].copy()
    for s in SRC:
        leaf[f"phi_{s}"] = norm360(leaf[f"phi_{s}"])
    wide = leaf.pivot_table(index=["plant_id", group_col], columns="leaf_index",
                            values=[f"phi_{s}" for s in SRC] + [f"theta_{s}" for s in SRC],
                            aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    for s in SRC:
        dev = []
        for i in range(1, 5):
            c = f"Phi_dev_{s}_{i}"
            wide[c] = norm180(wide[f"phi_{s}_{i}"] - wide[f"phi_{s}_{i-1}"])
            dev.append(c)
        wide[f"med_phi_dev_{s}"] = wide[dev].median(axis=1, skipna=True)
        wide[f"am_phi_dev_{s}"] = wide[dev].mean(axis=1, skipna=True)
        tc = [f"theta_{s}_{i}" for i in LOWER5]
        wide[f"med_theta_{s}"] = wide[tc].median(axis=1, skipna=True)
        wide[f"am_theta_{s}"] = wide[tc].mean(axis=1, skipna=True)
    return wide


def trait_cols(src):
    c = [f"theta_{src}_{i}" for i in LOWER5]
    c += [f"Phi_dev_{src}_{i}" for i in range(1, 5)]
    c += [f"med_theta_{src}", f"am_theta_{src}", f"med_phi_dev_{src}", f"am_phi_dev_{src}"]
    return c  # 13 traits


# ---------------------------------------------------------------- MAIZE data
_RXM = re.compile(r"_\d+-(\d+)-([A-Za-z0-9]+)-(\d+)_(\d{4}-\d{2}-\d{2})_")
def maize_geno(pid):
    m = _RXM.search(pid)
    return (m.group(2), f"{m.group(2)}-{m.group(3)}") if m else (None, None)


def load_maize():
    e1b = pd.read_csv(MZ / "e1b_per_leaf_pymaize_vs_gold.csv")
    e5 = pd.read_csv(MZ / "e5_roundtrip_vs_gold.csv")
    leaf = (e1b[["plant_id", "leaf_index", "theta_gold", "phi_gold",
                 "theta_ours", "phi_ours"]]
            .rename(columns={"theta_ours": "theta_e1b", "phi_ours": "phi_e1b"})
            .merge(e5[["plant_id", "leaf_index", "theta_rt", "phi_rt"]]
                   .rename(columns={"theta_rt": "theta_e5", "phi_rt": "phi_e5"}),
                   on=["plant_id", "leaf_index"], how="inner"))
    gp = leaf.plant_id.map(maize_geno)
    leaf["geno"] = [a for a, b in gp]
    leaf["physical"] = [b for a, b in gp]
    leaf = leaf[leaf.geno.notna()].copy()
    return leaf


def maize_plant_table():
    leaf = load_maize()
    # collapse lifecycle scans -> one median value per physical plant, >=5-leaf scans
    scans = wide_from_leaf(leaf.rename(columns={"physical": "physical_"}), "geno")
    # need physical plant id: re-derive from plant_id
    scans["physical"] = scans.plant_id.map(lambda p: maize_geno(p)[1])
    need = [f"theta_gold_{i}" for i in LOWER5]
    scans = scans[scans[need].notna().sum(axis=1) >= 5].copy()
    feat = [c for c in scans.columns if c.startswith(("theta_", "phi_", "Phi_dev_",
            "med_theta_", "am_theta_", "med_phi_dev_", "am_phi_dev_"))]
    plant = scans.groupby(["physical", "geno"])[feat].median().reset_index()
    reps = plant.groupby("geno").size()
    plant = plant[plant.geno.isin(reps[reps >= 2].index)].copy()
    return plant


# ---------------------------------------------------------------- SORGHUM data
REFERENCE_GENO = "PI656058"
def load_sorghum_leaf():
    """Per-leaf sorghum table with PI genotype (mirrors notebook 06)."""
    sys.path.insert(0, str(REPO / "experiments"))
    from experiments import data_loaders as dl
    e1b = pd.read_csv(PY / "e1b_per_leaf_pymaize_vs_gold.csv")
    e5 = pd.read_csv(SOR_BK / "e5_roundtrip_vs_gold.csv")   # 324-plant backup
    leaf = (e1b[["plant_id", "leaf_index", "theta_gold", "phi_gold",
                 "theta_ours", "phi_ours"]]
            .rename(columns={"theta_ours": "theta_e1b", "phi_ours": "phi_e1b"})
            .merge(e5[["plant_id", "leaf_index", "theta_rt", "phi_rt"]]
                   .rename(columns={"theta_rt": "theta_e5", "phi_rt": "phi_e5"}),
                   on=["plant_id", "leaf_index"], how="inner"))
    def norm_js(s):
        if not isinstance(s, str):
            return None
        m = re.match(r"(?i)js(\d+)", s)
        return f"JS{m.group(1)}" if m else None
    geno, js_pi = dl.load_jensina_genotypes()
    js_pi["JS_norm"] = js_pi["JS_ID"].apply(norm_js)
    js_to_pi = dict(zip(js_pi["JS_norm"], js_pi["PI#"].str.replace("_", "").str.upper()))
    leaf["geno"] = leaf.plant_id.apply(
        lambda pid: js_to_pi.get(norm_js(dl.parse_plant_id(pid)["js_id"])))
    leaf = leaf[leaf.geno.notna() & (leaf.geno != REFERENCE_GENO)].copy()
    return leaf


def sorghum_plant_table():
    leaf = load_sorghum_leaf()
    wide = wide_from_leaf(leaf, "geno")            # one row per plant (324 -> replicated)
    reps = wide.groupby("geno").size()
    wide = wide[wide.geno.isin(reps[reps >= 2].index)].copy()
    return wide


# ---------------------------------------------------------------- R^2 bootstrap
def _r2(x, y):
    if x.size < 3:
        return np.nan
    r = np.corrcoef(x, y)[0, 1]
    return r * r


def r2_ci(d, gcol, mcol, cluster_col, nb=NB):
    """Clustered bootstrap of R^2=r^2 between measured (mcol) and gold (gcol).
    `d` is already filtered to the analysis window. Resamples clusters
    (physical plants) with replacement via a numpy fast path."""
    d = d[[cluster_col, gcol, mcol]].dropna()
    codes, clusters = pd.factorize(d[cluster_col])
    g = d[gcol].to_numpy(); m = d[mcol].to_numpy()
    # per-cluster row arrays
    order = np.argsort(codes, kind="stable")
    codes_s, g_s, m_s = codes[order], g[order], m[order]
    bounds = np.searchsorted(codes_s, np.arange(clusters.size + 1))
    gx = [g_s[bounds[i]:bounds[i + 1]] for i in range(clusters.size)]
    mx = [m_s[bounds[i]:bounds[i + 1]] for i in range(clusters.size)]
    point = _r2(g, m)
    boot = np.empty(nb)
    K = clusters.size
    for b in range(nb):
        pick = RNG.integers(0, K, K)
        gb = np.concatenate([gx[i] for i in pick])
        mb = np.concatenate([mx[i] for i in pick])
        boot[b] = _r2(gb, mb)
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return dict(point=point, lo=lo, hi=hi, n_plants=K, n_leaves=len(d))


# ---------------------------------------------------------------- H^2 bootstrap
def h2_bootstrap(plant, nb=NB):
    """Resample genotypes with replacement; per resample median H^2 across 13 traits/source.
    Returns per-source (point, lo, hi) of median-H^2 and P(proc median > voxel median)."""
    genos = plant.geno.unique()
    # precompute per-source trait columns
    tcols = {s: trait_cols(s) for s in SRC}
    # point
    def median_h2(df, src):
        gcode = pd.factorize(df.geno)[0]
        h2 = np.array([mom_h2(df[c].to_numpy(), gcode) for c in tcols[src]], float)
        return np.nanmedian(h2)
    point = {s: median_h2(plant, s) for s in SRC}
    by_g = {g: plant[plant.geno == g] for g in genos}
    boot = {s: np.empty(nb) for s in SRC}
    proc_beats_vox = 0
    for b in range(nb):
        pick = RNG.choice(genos, genos.size, replace=True)
        # relabel to keep resampled genotypes distinct
        parts = []
        for j, g in enumerate(pick):
            p = by_g[g].copy()
            p["geno"] = f"{g}__{j}"
            parts.append(p)
        rs = pd.concat(parts, ignore_index=True)
        mh = {s: median_h2(rs, s) for s in SRC}
        for s in SRC:
            boot[s][b] = mh[s]
        if mh["e5"] > mh["gold"]:
            proc_beats_vox += 1
    out = {}
    for s in SRC:
        lo, hi = np.nanpercentile(boot[s], [2.5, 97.5])
        out[s] = dict(point=point[s], lo=lo, hi=hi)
    out["frac_proc_beats_voxel"] = proc_beats_vox / nb
    return out


# ---------------------------------------------------------------- Chamfer/Hausdorff
def geom_iqr(e4_path, species):
    d = pd.read_csv(e4_path)
    rows = []
    for pfx, col in [("vo", "Skeleton"), ("vp", "Procedural")]:
        for metric, mname in [("cd", "Chamfer"), ("hd", "Hausdorff")]:
            x = d[f"{pfx}_{metric}_sym_cm"].dropna()
            rows.append(dict(species=species, rep=col, metric=mname, n=len(x),
                             median=x.median(), q1=x.quantile(.25), q3=x.quantile(.75),
                             iqr=x.quantile(.75) - x.quantile(.25)))
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()

    mplant = maize_plant_table()
    print(f"maize plant table: {len(mplant)} plants, {mplant.geno.nunique()} genotypes")
    splant = sorghum_plant_table()
    print(f"sorghum plant table: {len(splant)} plants, {splant.geno.nunique()} genotypes")

    print("\n=== validate MoM H2 vs published REML ===")
    print("  (maize target .77/.82/.79 ; sorghum target .32/.29/.48)")
    for tag, pl in [("maize", mplant), ("sorghum", splant)]:
        for src, name in zip(SRC, ["voxel", "skeleton", "procedural"]):
            gcode = pd.factorize(pl.geno)[0]
            h2 = np.array([mom_h2(pl[c].to_numpy(), gcode) for c in trait_cols(src)], float)
            print(f"  {tag:8s} {name:11s}: median={np.nanmedian(h2):.3f}  "
                  f"mean={np.nanmean(h2):.3f}  >=0.20:{(h2>=0.20).sum()}/13")

    if not args.full:
        sys.exit(0)

    # ---- Part B: geometry IQR ----
    print("\n=== Chamfer / Hausdorff median + IQR (cm) ===")
    grows = geom_iqr(SOR_E4, "sorghum") + geom_iqr(MZ / "e4_geometric_fidelity.csv", "maize")
    gdf = pd.DataFrame(grows)
    print(gdf.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    GENERATED.mkdir(exist_ok=True)
    gdf.to_csv(GENERATED / "ci_geometry_iqr.csv", index=False)

    # ---- Part A: R^2 bootstrap ----
    # Window that reproduces the published tables: lower-4 leaves; maize uses NO
    # keep-filter, sorghum uses the per-leaf keep-filter (matches Tables 1 & 2).
    print("\n=== R^2 bootstrap over plants (2000x, lower-4 leaves) ===")

    def src_frame(e1b_path, e5_path, cluster_fn, use_keep):
        e1b = pd.read_csv(e1b_path)
        e5 = pd.read_csv(e5_path)
        frames = {}
        # skeleton (S) from e1b: theta_ours/phi_ours vs gold
        sk = e1b[["plant_id", "leaf_index", "theta_gold", "phi_gold",
                  "theta_ours", "phi_ours"] + (["keep"] if "keep" in e1b.columns else [])].copy()
        sk = sk[sk.leaf_index < 4]
        if use_keep and "keep" in sk.columns:
            sk = sk[sk.keep]
        sk["cluster"] = sk.plant_id.map(cluster_fn)
        frames["Skeleton"] = (sk, "theta_ours", "phi_ours")
        # procedural (P) from e5: theta_rt/phi_rt vs gold
        pr = e5[["plant_id", "leaf_index", "theta_gold", "phi_gold",
                 "theta_rt", "phi_rt"] + (["keep"] if "keep" in e5.columns else [])].copy()
        pr = pr[pr.leaf_index < 4]
        if use_keep and "keep" in pr.columns:
            pr = pr[pr.keep]
        pr["cluster"] = pr.plant_id.map(cluster_fn)
        frames["Procedural"] = (pr, "theta_rt", "phi_rt")
        return frames

    maize_frames = src_frame(MZ / "e1b_per_leaf_pymaize_vs_gold.csv",
                             MZ / "e5_roundtrip_vs_gold.csv",
                             lambda p: maize_geno(p)[1], use_keep=False)
    sorg_frames = src_frame(PY / "e1b_per_leaf_pymaize_vs_gold.csv",
                            SOR_BK / "e5_roundtrip_vs_gold.csv",
                            lambda p: p, use_keep=True)
    r2rows = []
    for tag, frames in [("maize", maize_frames), ("sorghum", sorg_frames)]:
        for rep, (df, tcol, pcol) in frames.items():
            for trait, mcol, gcol in [("R2(theta)", tcol, "theta_gold"),
                                      ("R2(phi)", pcol, "phi_gold")]:
                r = r2_ci(df, gcol, mcol, "cluster")
                r.update(species=tag, rep=rep, trait=trait)
                r2rows.append(r)
                print(f"  {tag:8s} {rep:10s} {trait}: {r['point']:.3f} "
                      f"[{r['lo']:.3f}, {r['hi']:.3f}]  n_plants={r['n_plants']} n_leaves={r['n_leaves']}")
    pd.DataFrame(r2rows).to_csv(GENERATED / "ci_r2_bootstrap.csv", index=False)

    # ---- Part C: H^2 bootstrap over genotypes ----
    print("\n=== H2 bootstrap over genotypes (2000x) ===")
    for tag, pl in [("maize", mplant), ("sorghum", splant)]:
        res = h2_bootstrap(pl)
        for src, name in zip(SRC, ["voxel", "skeleton", "procedural"]):
            r = res[src]
            print(f"  {tag:8s} {name:11s} median H2 = {r['point']:.3f} "
                  f"[{r['lo']:.3f}, {r['hi']:.3f}]")
        print(f"  {tag:8s} P(procedural median H2 > voxel median H2) = {res['frac_proc_beats_voxel']:.3f}")
