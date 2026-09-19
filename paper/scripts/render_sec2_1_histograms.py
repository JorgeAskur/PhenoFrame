"""Section 2.1 figure (delta-histogram version): descriptor self-consistency.
For each of the four reported quantities, the distribution of the residual
(measured - declared) over 500 synthetic descriptors / 6,213 leaves. A tight
spike at 0 shows the trait pipeline returns the values written into the
descriptor; the leaf-angle panel is deliberately wider (droopiness bends the
midrib tangent, an expected effect, not a recovery error).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAPER = Path(__file__).resolve().parents[1]
CSV = PAPER / "source_data" / "sorghum" / "descriptor_vs_measured_all_leaves.csv"
OUT = PAPER / "generated" / "figures"
OUT.mkdir(exist_ok=True, parents=True)

d = pd.read_csv(CSV).sort_values(["plant", "leaf_id"])
d["xml_cum_distance"] = d.groupby("plant")["xml_distance"].cumsum()
az = d[d["xml_angle"].abs() <= 85].copy()
n_flip = len(d) - len(az)


def wrap180(x):
    return (x + 180.0) % 360.0 - 180.0


# (title, frame, declared, measured, scale, unit, color, wrap, symbol)
panels = [
    ("Leaf length",      d,  "xml_length",       "py_length",      1000.0, "mm",  "#3a8c5e", False, "L"),
    ("Leaf azimuth",     az, "xml_azimuth",      "py_azimuth",     1.0,    "deg", "#3672b5", True,  r"\phi"),
    ("Internode distance", d, "xml_cum_distance", "py_conn_y",     1000.0, "mm",  "#e8a838", False, "d"),
    ("Leaf angle",       d,  "xml_angle",        "py_inclination", 1.0,    "deg", "#c25b56", False, r"\theta"),
]

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, (name, dd, xc, yc, scale, unit, color, wrap, sym) in zip(axes.ravel(), panels):
    delta = (dd[yc].to_numpy() - dd[xc].to_numpy())
    if wrap:
        delta = wrap180(delta)
    delta = delta[np.isfinite(delta)] * scale
    med = np.median(delta)
    mad = np.median(np.abs(delta))                       # typical error magnitude
    rng = (np.nanmax(dd[xc]) - np.nanmin(dd[xc])) * scale  # trait range, same units
    relpct = 100.0 * mad / rng if rng else float("nan")
    lim = np.percentile(np.abs(delta), 99.5)
    lim = lim if lim > 0 else (np.abs(delta).max() or 1.0)
    bins = np.linspace(-lim, lim, 61)
    ax.hist(np.clip(delta, -lim, lim), bins=bins, color=color, alpha=0.85, edgecolor="white")
    ax.axvline(med, color="black", ls="--", lw=1,
               label=f"median |Δ| = {mad:.3g} {unit}\n({relpct:.2g}% of trait range)")
    ax.set_xlabel(rf"$\Delta {sym}$ ({unit})")
    ax.set_ylabel("Count")
    ax.set_title(name, fontsize=12)
    ax.set_xlim(-lim, lim)
    ax.legend(fontsize=8)

fig.suptitle(f"Measured − declared difference (N = 500 descriptors, {len(d):,} leaves)",
             fontsize=13, y=0.99)
fig.tight_layout(rect=(0, 0, 1, 0.97))
out = OUT / "sec2_1_self_consistency_hist.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
plt.close(fig)
print("saved", out)
for name, dd, xc, yc, scale, unit, color, wrap, note in panels:
    delta = dd[yc].to_numpy() - dd[xc].to_numpy()
    if wrap:
        delta = wrap180(delta)
    delta = delta[np.isfinite(delta)] * scale
    print(f"  {name:16s} bias={np.mean(delta):+.3g}{unit} MAE={np.mean(np.abs(delta)):.3g}{unit} "
          f"max={np.max(np.abs(delta)):.3g}{unit}")
