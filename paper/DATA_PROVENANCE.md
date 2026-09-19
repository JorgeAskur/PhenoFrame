# Data provenance and scope

The files in `source_data/` are compact analysis inputs and paper figures,
not copies of the raw image series or complete voxel archives. Run
`python paper/verify_bundle.py` to check the expected row and population
counts. `source_data/sorghum/` contains the 324-plant paper analysis; its
E4/E5 CSVs came from the preserved `outputs_sorghum_backup` run, not the
later exploratory outputs in `experiments/experiments/results/outputs/`.

| Input | Tracked here | Original source | Role |
| --- | --- | --- | --- |
| Sorghum E1b/E4/E5/E8 CSVs | Yes | Gaillard voxel reconstruction and Davis reference traits | Figures and Tables 1-3 |
| GIC maize E1b/E4/E5/E8 CSVs | Yes | 20-genotype GIC lifecycle voxel and skeleton scans | Figures and Tables 1-3 |
| Generation CSVs/configs | Yes | Sorghum-fitted stochastic generator | Section 2.6 comparisons |
| Self-consistency XML descriptors | Yes, 500 | ProceduralModelMaize generator | Section 2.2 validation |
| Davis phyllotaxy mapping CSVs | Yes, under `reference/` | `jdavis-132/phyllotaxy` repository | Genotype mapping |
| Full sorghum voxel archive | No | https://doi.org/10.5281/zenodo.4426620 | Rerunning voxel notebooks |
| GIC maize voxel/skeleton scans | External deposit | Compact archive prepared by `paper/create_gic_deposit.py`; DOI pending | Rerunning maize notebooks |
| Full PANICLE GWAS scripts and marker matrix | No | Not located in this checkout | Genome-wide analysis |

Maize is longitudinal: 80 physical plants from 20 genotypes and four
replicates were scanned repeatedly. The paired E4/E5 population has 3,317
scan IDs; the unpaired E1b export has 3,318. A scan ID is not an independent
plant. Table 3 collapses timepoints to physical plants before estimating
heritability. Sorghum E4/E5 use 324 plant IDs; the replicated heritability
subset has 128 plants from 44 genotypes.

The flat maize voxel tree holds a directory per scan, containing
`voxels.txt`, `optim_skeleton.txt`, and `angles.txt`. The loader expects
`<voxel-root>/reconstructed/<scan>/voxels.txt` and
`<voxel-root>/skeletons/<scan>/{optim_skeleton.txt,angles.txt}`. Create two
directory links to the same flat tree rather than copying the dataset. On
Windows, use PowerShell `New-Item -ItemType Junction`; on POSIX, use symlinks.
The `paper/scripts/run_gic_notebook.py` runner takes that staged root with
`--voxel-root`. The Gaillard archive has the loader's native split layout.

The stochastic executable is not vendored; `generator.lock` pins its
separate repository and revision. The copied Davis mapping CSVs should be
checked against the upstream redistribution terms before public release.
The compact GIC deposit excludes regenerated OBJ meshes and contains the voxel,
skeleton, angle, and status text files needed by the PhenoFrame experiments.
The verified archive is 569,712,966 bytes with SHA-256
`88f6dc5a40a598f17323db4baab3796d3d4b2e65344c4738cf4b5624722f8b1a`.
The missing genome-wide GWAS inputs/scripts remain a separate reproducibility
dependency. Until those are available, the included tables and figures are
reproducible from tracked analysis data, but the full paper is not end-to-end
reproducible.
