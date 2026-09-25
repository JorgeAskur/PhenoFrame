"""Regenerate the two maize figures that had diverged from the notebook/paper
style, using the EXACT notebook plotting code so they match the sorghum figures.

  * e5_fitted_parameter_distributions.png  -- verbatim notebook-05 code, with
    the boxplot node range capped to ranks with >=20 samples (maize has up to
    22 nodes; the tail would otherwise collide and show singleton boxes).
  * e8_h2_by_trait_maize.png               -- verbatim notebook-06 code
    (nice_label math ticks, same colors/legend/threshold line).
"""
from pathlib import Path
from _data_paths import GENERATED, MAIZE
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = MAIZE
FIG = GENERATED / "figures" / "maize"
FIG.mkdir(parents=True, exist_ok=True)
FIG_DIR = FIG

# =====================================================================
# (1) e5 fitted parameter distributions -- notebook-05 style, capped nodes
# =====================================================================
fits = pd.read_csv(OUT / "e5_fitted_params.csv")
NODE_MIN_N = 20
cnt = fits.groupby("leaf_index").size()
box_nodes = sorted(int(n) for n, c in cnt.items() if c >= NODE_MIN_N)

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
    ax.set_xlabel(label)
    ax.set_ylabel("Count")
    ax.legend(fontsize=8)

for ax, col, label in [
    (axes[1, 0], "leaf_angle",  "leaf_angle (deg)"),
    (axes[1, 1], "droopiness",  "droopiness"),
    (axes[1, 2], "leaf_length", "leaf_length (m)"),
]:
    data = [fits[fits["leaf_index"] == n][col].dropna().values for n in box_nodes]
    bp = ax.boxplot(data, positions=box_nodes, widths=0.6, patch_artist=True,
                    medianprops=dict(color="black", lw=1.5),
                    flierprops=dict(marker=".", markersize=2, alpha=0.3))
    for patch, n in zip(bp["boxes"], box_nodes):
        patch.set_facecolor(cmap(n / max(max_node, 1)))
        patch.set_alpha(0.7)
    ax.set_xlabel("Leaf node (bottom → top)")
    ax.set_ylabel(label)
    ax.set_xticks(box_nodes)
    ax.set_xticklabels([str(n) for n in box_nodes])

fig.suptitle(f"Fitted procedural parameters, N = {len(fits)} leaves across "
             f"{fits['plant_id'].nunique()} plants")
plt.tight_layout()
plt.savefig(FIG_DIR / "e5_fitted_parameter_distributions.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("saved e5_fitted_parameter_distributions.png (nodes", box_nodes[0], "..", box_nodes[-1], ")")

# =====================================================================
# (2) e8 heritability -- verbatim notebook-06 code
# =====================================================================
wide_h2 = pd.read_csv(OUT / "e8_h2_by_trait_maize.csv")
n_plants, n_geno, LOWER_LEAVES = 80, 20, range(5)


def nice_label(raw):
    if raw.startswith("theta_"):
        return r"$\theta_{" + raw.split("_")[1] + "}$"
    if raw.startswith("Phi_dev_"):
        return r"$\Phi_{" + raw.split("_")[2] + "}$"
    aggregates = {
        "med_theta":   r"median $\theta$",
        "am_theta":    r"mean $\theta$",
        "med_phi_dev": r"median $\Phi$",
        "am_phi_dev":  r"mean $\Phi$",
    }
    return aggregates.get(raw, raw)


fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
panels = [
    ("per-node theta",   r"Per-node leaf angle $\theta_i$"),
    ("per-node Phi_dev", r"Per-node phyllotaxic deviation $\Phi_i$"),
    ("per-plant",        "Per-plant aggregates"),
]
for ax, (cat_match, title) in zip(axes, panels):
    sub = wide_h2[wide_h2["category"].str.startswith(cat_match)].copy().reset_index(drop=True)
    labels = [nice_label(l) for l in sub["label"].tolist()]
    x = np.arange(len(labels))
    width = 0.27
    ax.bar(x - width, sub["gold"].values, width, label="Voxel",
           color="#999999", edgecolor="black", linewidth=0.5)
    ax.bar(x,         sub["e1b"].values,  width, label="Skeleton",
           color="#4477aa", edgecolor="black", linewidth=0.5)
    ax.bar(x + width, sub["e5"].values,   width, label="Procedural",
           color="#cc6677", edgecolor="black", linewidth=0.5)
    # Maize theta bars are near-ceiling; the middle panel leaves room for the legend.
    ax.axhline(0.20, color="black", linestyle="--", linewidth=0.8, alpha=0.6,
               label=r"Davis 2025 $H^2=0.20$ threshold" if ax is axes[1] else None)
    ax.set_xticks(x)
    if cat_match == "per-plant":
        ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=11)
    else:
        ax.set_xticklabels(labels, rotation=0, fontsize=13)
    ax.tick_params(axis="x", pad=2)
    ax.set_ylabel(r"Broad-sense heritability $H^2$  (Davis Eq. 3, $n=2$)")
    ax.set_title(title, fontsize=11)
    ax.set_ylim(0, 1.0)
    ax.grid(axis="y", alpha=0.3)
    if ax is axes[1]:
        ax.legend(loc="upper right", fontsize=9)
fig.suptitle(f"Heritability across pipelines  (n={n_plants} plants, {n_geno} genotypes, "
             f"lower {len(LOWER_LEAVES)} leaves)", y=1.02, fontsize=12)
fig.tight_layout()
fig.savefig(FIG / "e8_h2_by_trait_maize.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("saved e8_h2_by_trait_maize.png (notebook-06 style)")
