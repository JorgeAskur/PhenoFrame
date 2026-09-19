# PhenoFrame paper reproduction

This directory contains compact inputs, analysis scripts, generator
configurations, and 500 descriptor XMLs used for the manuscript
*PhenoFrame: Phenotyping Plant Descriptor and Software Suite*. Paths in the
portable table scripts are relative to this directory. Generated outputs go
under paper/generated/ and are not tracked.
See [DATA_PROVENANCE.md](DATA_PROVENANCE.md) for source details and the
remaining external dependencies.

## Environment

From the repository root, with Python 3.10 or newer:

    python -m pip install -e ".[experiments,test]"
    python -m pytest -q
    python paper/verify_bundle.py

The phenoframe package contains the descriptor parser, skeleton conversion,
trait pipeline, and C++ mesh/trait API wrapper. The Windows DLL is bundled.
On Linux or macOS, build the C++ library with MaizeModel/build_api.sh.

## Results from the included CSVs

    python paper/scripts/round3_B1.py
    python paper/scripts/R123.py
    python paper/scripts/render_sec2_1_histograms.py

round3_B1.py reproduces the paired lower-four-leaf Table 1/2 trait values
and whole-plant geometry medians. Select its "(b) both-keep" rows; its
heritability rows are a sensitivity analysis, not the final Table 3.
R123.py computes Table 3 with REML, shared leaf quality filters, a common
divergence-pair filter, and maize timepoints collapsed to 80 physical plants.
The saved reference results are
source_data/summary/B1_table2_by_leafset.csv and
source_data/summary/R123_table3_definitive.csv.
R89.py reruns the 1,000-resample genotype bootstrap for Table 3 and formats
generation supplement tables; it takes substantially longer.

The source_data/sorghum/ E1b CSV comes from the 324-plant run. Its E4/E5
CSVs came from the preserved outputs_sorghum_backup run; the live notebook
outputs had subsequently been overwritten by a larger exploratory run.
source_data/maize/ contains the GIC 20-genotype results from 3,317
reconstructed scans. source_data/generation/ contains measured traits and
comparison tables from the 324-plant sorghum generator run.
source_data/summary/ retains uncertainty and sensitivity tables used during
manuscript preparation. The published Davis genotype mapping is under
reference/phyllotaxy/Data/.
The paper-matching figure PNGs are under source_data/figures/; run
paper/scripts/combine_species_figures.py to assemble the five cross-species
panels in paper/generated/figures/combined/. The sorghum and maize E4/E5
plot scripts can regenerate their corresponding panels from the bundled CSVs.
The E1b maize source CSV contains 3,318 scans; the paired E4/E5 population
contains 3,317 scans. Do not call either count unique plants.

## Regenerating descriptors and notebooks

The stochastic C++ generator is maintained in a separate repository. Clone
the URL and check out the exact revision in generator.lock, then build its
Maize.exe using that repository's README. Set PHENOFRAME_GENERATOR_EXE to
the executable path. The three configuration XMLs used here are frozen in
configs/.

    python paper/generate_synthetic.py --count 324

This writes the main synthetic descriptor population to
paper/generated/sorghum/. The tracked synthetic_descriptors/ directory
holds the separate 500-plant self-consistency population. To rerun the
sorghum notebooks, obtain the Gaillard voxel archive
(https://doi.org/10.5281/zenodo.4426620) and the Davis phyllotaxy data
(https://github.com/jdavis-132/phyllotaxy), set PHENOFRAME_VOXEL_PATH and
PHENOFRAME_JENSINA_PATH, then run
experiments/experiments/results/_run_all.py. The notebook kernel uses
experiments/experiments/results/ as its working directory.
The tracked source CSVs and figures suffice to regenerate the paper's
summary tables and combined figures without downloading either voxel archive.
Recomputing those CSVs from skeletons requires the external archives.

For maize, point the runner at a root containing reconstructed/ and
skeletons/ for the GIC scans:

    python paper/scripts/run_gic_notebook.py --voxel-root C:\path\to\gic_voxel_root 03_sec2.2_phenoframe_traits_vs_gold.ipynb 04_sec2.3_geometric_fidelity.ipynb 05_sec2.3_procedural_round_trip.ipynb

The GIC data comprise 80 physical plants (20 genotypes x 4 replicates)
imaged repeatedly; 3,317 scans entered the geometric analysis. Raw maize
imagery and voxel grids are not included in Git. They need an archival
location and access instructions before an independent researcher can
rebuild the descriptor CSVs from images. The original sorghum archive is
roughly 34 GB and is likewise linked rather than copied.
The local GIC download is a flat scan directory. Create two directory links
named `reconstructed` and `skeletons` under a new voxel root, both pointing
to that flat directory; see DATA_PROVENANCE.md. This avoids duplicating the
large files.

For one-plant four-format renders, use
paper/scripts/render_sorghum_paper_plant.py --voxel-root <sorghum-root> or
paper/scripts/render_maize_custom_plant.py --voxel-root <maize-root>. For the
per-leaf sorghum angle export, use
paper/scripts/sorghum_phyllotaxis_csv.py --voxel-root <sorghum-root>.

## Remaining external material

The manuscript's genome-wide PANICLE analysis is distinct from the tracked
targeted-marker notebook 10_phase1_gwas_marker_replication.ipynb. The full
GWAS scripts, combined marker archive, and maize data DOI are not available
in this checkout. They must be deposited and linked before the Code and Data
availability statements can claim complete end-to-end reproduction.
