"""Stack the sorghum (top) and maize (bottom) version of each experiment
figure into a single paper-ready panel, with a bold species row label."""
from pathlib import Path
from _data_paths import GENERATED, PAPER
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

SORG = PAPER / "source_data" / "figures" / "sorghum"
MAIZE = PAPER / "source_data" / "figures" / "maize"
OUT = GENERATED / "figures" / "combined"
OUT.mkdir(exist_ok=True, parents=True)

# (output name, sorghum file, maize file)
PAIRS = [
    ("combined_sec2_2_four_formats",
     SORG / "sec2_2_plant_four_formats.png",
     MAIZE / "sec2_2_maize_plant_developed.png"),
    ("combined_sec2_2_delta_by_node",
     SORG / "e1_delta_by_node.png",
     MAIZE / "e1_delta_by_node.png"),
    ("combined_sec2_3_chamfer_hausdorff",
     SORG / "e4_chamfer_distributions.png",
     MAIZE / "e4_chamfer_distributions.png"),
    ("combined_sec2_3_fitted_params",
     SORG / "e5_fitted_parameter_distributions.png",
     MAIZE / "e5_fitted_parameter_distributions.png"),
    ("combined_sec2_4_heritability",
     SORG / "e8_h2_by_trait.png",
     MAIZE / "e8_h2_by_trait_maize.png"),
]

W = 13.0  # common figure width (inches)

for name, sorg, maize in PAIRS:
    if not sorg.exists() or not maize.exists():
        print(f"SKIP {name}: missing {'sorghum' if not sorg.exists() else 'maize'} file")
        continue
    im1, im2 = mpimg.imread(sorg), mpimg.imread(maize)
    a1 = im1.shape[0] / im1.shape[1]
    a2 = im2.shape[0] / im2.shape[1]
    h1, h2 = W * a1, W * a2
    fig = plt.figure(figsize=(W, h1 + h2 + 0.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[h1, h2], hspace=0.04)
    for row, (im, lab) in enumerate([(im1, "Sorghum"), (im2, "Maize")]):
        ax = fig.add_subplot(gs[row])
        ax.imshow(im)
        ax.axis("off")
        ax.text(-0.015, 0.5, lab, transform=ax.transAxes, rotation=90,
                va="center", ha="right", fontsize=16, fontweight="bold")
    out = OUT / f"{name}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out.name}")

print(f"\nAll combined figures in {OUT}")
