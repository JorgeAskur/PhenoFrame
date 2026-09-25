"""Rebuild the PAPER-matching sorghum figures (324 plants) from the
outputs_sorghum_backup data, using the verbatim notebook plotting code.

The live experiments/.../outputs was overwritten by a 2832-plant run, so the
live e4/e5 figures no longer match the manuscript (which reports 324 plants,
2479 leaves). The matching E1 figure is preserved under source_data/figures/.
"""
import shutil
from _data_paths import GENERATED, PAPER, SORGHUM
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BK = SORGHUM
LIVE_FIG = PAPER / "source_data" / "figures" / "sorghum"
OUT = GENERATED / "figures" / "sorghum"
OUT.mkdir(exist_ok=True, parents=True)

# ---- e5 fitted procedural parameters (notebook-05 verbatim) ---------------
fits = pd.read_csv(BK / "e5_fitted_params.csv")
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
cmap = plt.get_cmap("viridis")
max_node = min(9, int(fits["leaf_index"].max()))
for ax, col, label, color in [
    (axes[0, 0], "leaf_angle",  "Fitted leaf_angle (deg)", "#c25b56"),
    (axes[0, 1], "droopiness",  "Fitted droopiness",       "#3672b5"),
    (axes[0, 2], "leaf_length", "Fitted leaf_length (m)",  "#3a8c5e"),
]:
    ax.hist(fits[col], bins=50, color=color, alpha=0.85, edgecolor="white")
    ax.axvline(fits[col].median(), color="black", ls="--", lw=1,
               label=f"median = {fits[col].median():.2f}")
    ax.set_xlabel(label); ax.set_ylabel("Count"); ax.legend(fontsize=8)
for ax, col, label in [
    (axes[1, 0], "leaf_angle",  "leaf_angle (deg)"),
    (axes[1, 1], "droopiness",  "droopiness"),
    (axes[1, 2], "leaf_length", "leaf_length (m)"),
]:
    nodes = sorted(fits["leaf_index"].unique())
    data = [fits[fits["leaf_index"] == n][col].dropna().values for n in nodes]
    bp = ax.boxplot(data, positions=nodes, widths=0.6, patch_artist=True,
                    medianprops=dict(color="black", lw=1.5),
                    flierprops=dict(marker=".", markersize=2, alpha=0.3))
    for patch, n in zip(bp["boxes"], nodes):
        patch.set_facecolor(cmap(n / max(max_node, 1))); patch.set_alpha(0.7)
    ax.set_xlabel("Leaf node (bottom → top)"); ax.set_ylabel(label)
fig.suptitle(f"Fitted procedural parameters, N = {len(fits)} leaves across "
             f"{fits['plant_id'].nunique()} plants")
plt.tight_layout()
plt.savefig(OUT / "e5_fitted_parameter_distributions.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"e5: {fits.plant_id.nunique()} plants, {len(fits)} leaves")

# ---- e4 chamfer/hausdorff (notebook-04 verbatim) --------------------------
df = pd.read_csv(BK / "e4_geometric_fidelity.csv")
pairs = [("vp_cd_sym_cm", "Voxel ↔ Procedural", "#c25b56"),
         ("vo_cd_sym_cm", "Voxel ↔ Skeleton",   "#3672b5"),
         ("op_cd_sym_cm", "Skeleton ↔ Procedural", "#e8a838")]
fig, axes = plt.subplots(1, 3, figsize=(17, 4.5))
ax = axes[0]
all_cd = pd.concat([df[c].dropna() for c, _, _ in pairs])
bins = np.linspace(0, all_cd.quantile(0.99) * 1.05, 40)
for col, label, color in pairs:
    s = df[col].dropna()
    ax.hist(s, bins=bins, alpha=0.55, color=color, label=f"{label} (med {s.median():.2f})", edgecolor="white")
    ax.axvline(s.median(), color=color, ls="--", lw=1.2)
ax.set_xlabel("Symmetric Chamfer distance (cm)"); ax.set_ylabel("Count")
ax.set_title("CD distributions"); ax.legend(fontsize=8)
ax = axes[1]
v = df.dropna(subset=["vp_cd_sym_cm", "vo_cd_sym_cm"])
ax.scatter(v["vo_cd_sym_cm"], v["vp_cd_sym_cm"], s=10, alpha=0.4, color="#3a8c5e")
lim = max(v["vp_cd_sym_cm"].max(), v["vo_cd_sym_cm"].max()) * 1.05
ax.plot([0, lim], [0, lim], "k--", lw=0.7, alpha=0.5)
ax.set_xlabel("Voxel ↔ Skeleton CD (cm)"); ax.set_ylabel("Voxel ↔ Procedural CD (cm)")
ax.set_title("Volume-based CD: Procedural vs Skeleton")
ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
ax = axes[2]
hd_pairs = [("vp_hd_sym_cm", "Voxel ↔ Procedural", "#c25b56"),
            ("vo_hd_sym_cm", "Voxel ↔ Skeleton",   "#3672b5"),
            ("op_hd_sym_cm", "Skeleton ↔ Procedural", "#e8a838")]
all_hd = pd.concat([df[c].dropna() for c, _, _ in hd_pairs])
bins_h = np.linspace(0, all_hd.quantile(0.99) * 1.05, 40)
for col, label, color in hd_pairs:
    s = df[col].dropna()
    ax.hist(s, bins=bins_h, alpha=0.55, color=color, label=f"{label} (med {s.median():.2f})", edgecolor="white")
    ax.axvline(s.median(), color=color, ls="--", lw=1.2)
ax.set_xlabel("Symmetric Hausdorff distance (cm)"); ax.set_ylabel("Count")
ax.set_title("HD distributions"); ax.legend(fontsize=8)
fig.suptitle(f"Three-way geometric fidelity, N = {len(df)} plants", fontsize=12)
plt.tight_layout()
plt.savefig(OUT / "e4_chamfer_distributions.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"e4: {len(df)} plants")

# ---- e8 heritability (notebook-06 verbatim, legend moved to middle panel
#      to match the maize figure) ------------------------------------------
wide_h2 = pd.read_csv(SORGHUM / "e8_h2_by_trait.csv")
n_plants, n_geno, n_lower = 128, 44, 5


def nice_label(raw):
    if raw.startswith("theta_"):
        return r"$\theta_{" + raw.split("_")[1] + "}$"
    if raw.startswith("Phi_dev_"):
        return r"$\Phi_{" + raw.split("_")[2] + "}$"
    return {"med_theta": r"median $\theta$", "am_theta": r"mean $\theta$",
            "med_phi_dev": r"median $\Phi$", "am_phi_dev": r"mean $\Phi$"}.get(raw, raw)


fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
panels = [("per-node theta",   r"Per-node leaf angle $\theta_i$"),
          ("per-node Phi_dev", r"Per-node phyllotaxic deviation $\Phi_i$"),
          ("per-plant",        "Per-plant aggregates")]
for ax, (cat_match, title) in zip(axes, panels):
    sub = wide_h2[wide_h2["category"].str.startswith(cat_match)].copy().reset_index(drop=True)
    labels = [nice_label(l) for l in sub["label"].tolist()]
    x = np.arange(len(labels)); width = 0.27
    ax.bar(x - width, sub["gold"].values, width, label="Voxel",
           color="#999999", edgecolor="black", linewidth=0.5)
    ax.bar(x,         sub["e1b"].values,  width, label="Skeleton",
           color="#4477aa", edgecolor="black", linewidth=0.5)
    ax.bar(x + width, sub["e5"].values,   width, label="Procedural",
           color="#cc6677", edgecolor="black", linewidth=0.5)
    ax.axhline(0.20, color="black", linestyle="--", linewidth=0.8, alpha=0.6,
               label=r"Davis 2025 $H^2=0.20$ threshold" if ax is axes[1] else None)
    ax.set_xticks(x)
    if cat_match == "per-plant":
        ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=11)
    else:
        ax.set_xticklabels(labels, rotation=0, fontsize=13)
    ax.tick_params(axis="x", pad=2)
    ax.set_ylabel(r"Broad-sense heritability $H^2$  (Davis Eq. 3, $n=2$)")
    ax.set_title(title, fontsize=11); ax.set_ylim(0, 1.0); ax.grid(axis="y", alpha=0.3)
    if ax is axes[1]:
        ax.legend(loc="upper right", fontsize=9)
fig.suptitle(f"Heritability across pipelines  (n={n_plants} plants, {n_geno} genotypes, "
             f"lower {n_lower} leaves)", y=1.02, fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "e8_h2_by_trait.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("e8: regenerated with legend on middle panel")

# ---- copy the already-paper-matching §2.2 figure across -------------------
shutil.copy(LIVE_FIG / "e1_delta_by_node.png", OUT / "e1_delta_by_node.png")
print("copied e1_delta_by_node.png")

print(f"\nPaper-matching sorghum figures in {OUT}")
