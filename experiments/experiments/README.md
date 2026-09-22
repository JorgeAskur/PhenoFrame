# PhenoFrame validation experiments

This directory holds the validation experiments that back the PhenoFrame paper:
trait extractor concordance, descriptor-format fidelity, procedural round-trip,
heritability preservation, and population-level distribution match.

For the submission-ready input bundle, exact paper sample sizes, and runnable
table scripts, see `../../paper/README.md`. The detailed counts below describe
historical notebook runs and may use an earlier leaf filter or manuscript
section numbering.

> **Scope and section mapping.** The notebooks here cover the **sorghum**
> (out-of-calibration) experiments on the Gaillard/Davis voxel data. The
> **maize** (in-calibration, GIC 20-genotype) experiments reported in the paper
> have their compact result CSVs and runner in `../../paper/`. The `§2.x` tags
> below follow an earlier draft's numbering; the
> submitted paper uses §2.2 (self-consistency), §2.3 (real architecture),
> §2.4 (procedural geometry), §2.5 (heritability), §2.6 (generation), and
> §2.7 (GWAS). Add one to each `§2.x` tag below to match the submitted section.

## Layout

```
experiments/
  paths.py                # Filesystem paths to external datasets (env-overridable)
  data_loaders.py         # Pandas/numpy loaders for every external file
  README.md               # This file
  EXPERIMENT_RESULTS.md   # Detailed results report for all experiments
  results/
    00_data_exploration.ipynb                    # Prereq: explore both datasets (E0)
    01_sec2.1_python_cpp_validation.ipynb        # §2.1 / E3: Python vs C++ pipeline parity
    02_sec2.1_descriptor_trait_validation.ipynb  # §2.1 / E6: XML descriptor vs measured trait
    03_sec2.2_phenoframe_traits_vs_gold.ipynb       # §2.2 / E1b: PhenoFrame traits vs gold standard
    04_sec2.3_geometric_fidelity.ipynb           # §2.3 / E4: Chamfer/Hausdorff mesh vs voxel
    05_sec2.3_procedural_round_trip.ipynb        # §2.3 / E5: procedural-fit trait round-trip
    06_sec2.4_heritability.ipynb                 # §2.4 / E8: broad-sense H² (gold/E1b/E5)
    07_sec2.5_distribution_match.ipynb           # §2.5 / E9: synthetic-plant population match (parameters)
    08_sec2.5_data_driven_generation.ipynb       # §2.5 / E9b: 4-latent-factor generator -> C++ engine -> measured trait match
    09_sec2.5_holdout_validation.ipynb           # §2.5 / E9c: 80/20 held-out validation with measurement-gap baseline
    10_phase1_gwas_marker_replication.ipynb      # §2.6 / E10: targeted GWAS marker-effect replication at Davis 2025 hits
    99_appendix_e1_trait_concordance.ipynb       # Appendix / E1: PCA re-implementation check
    derive_sorghum_config.py                     # Fits sorghum generator config from E5 population (calls phenoframe.MaizeGenerator)
    render_sec2_2_plant_formats.py               # 3D-viz: one plant in voxel / skeleton / mesh formats
    figures/                                      # Generated PNGs (gitignored)
    outputs/                                      # Generated CSVs + OBJs (gitignored)
```

**C++ generator integration.** Sections §2.5 (notebooks 08, 09) and §2.6 (notebook 10) invoke the C++ stochastic generator from the sibling `MaizeProceduralModel` repository. Calls go through `phenoframe.MaizeGenerator`, which discovers the `Maize.exe` binary, derives a fitted config XML, and runs `Maize.exe --headless` to produce a batch of descriptor XMLs. See the main repo `README.md` for the wrapper's API; see `derive_sorghum_config.py` for the project-specific calibration pipeline.

## Experiment status

### Completed

| Paper section | ID | Notebook | What it tests | Key result |
|----------------|----|----------|----------------|-------------|
| (prereq) | E0 | `00_data_exploration.ipynb` | Dataset sanity checks, plant/leaf counts, data availability | 332 plants with usable skeletons, 351 total |
| §2.1 | E3 | `01_sec2.1_python_cpp_validation.ipynb` | Python forward model vs C++ engine parity across 500 plants (6,213 leaves) | Length MAE = 0.0001 mm, azimuth MAE = 0.0001°, inclination MAE = 0.0015°. Pipelines are functionally identical. |
| §2.1 | E6 | `02_sec2.1_descriptor_trait_validation.ipynb` | XML descriptor parameter values vs measured traits from both pipelines | leafLength MAE = 0.10 mm (r = 1.0), azimuth MAE = 0.39° excluding near-vertical leaves, distance vs height r = 0.999999. Full round-trip validation. |
| §2.2 | E1b | `03_sec2.2_phenoframe_traits_vs_gold.ipynb` | PhenoFrame descriptor format + trait pipeline vs gold-standard θ and φ | r(θ) = 0.69, R²(θ) = 0.48; r(φ) = 0.73, R²(φ) = 0.53 (filtered lower 4 leaves). PhenoFrame R²(φ) = 0.53 exceeds the Davis et al. 2025 cross-timepoint repeatability (R² = 0.41) and method-vs-manual agreement (R² = 0.48) for φ, approaching the inter-observer manual ceiling (R² = 0.55). Median \|error\| not reported in Davis 2025 (R² only). |
| §2.3 | E4 | `04_sec2.3_geometric_fidelity.ipynb` | 3D geometric fidelity: Chamfer/Hausdorff distance between descriptor meshes and voxel skeletons | Median CD vs voxel: 3.82 cm (procedural) vs 3.58 cm (override) — virtually identical. The 3-parameter compression is geometrically lossless at the whole-plant level. Median Hausdorff vs voxel: 20.03 cm (procedural) vs 18.04 cm (override) — ~2 cm tip penalty from the procedural form. |
| §2.3 | E5 | `05_sec2.3_procedural_round_trip.ipynb` | Procedural descriptor (3 params/leaf) round-trip fidelity | r(θ) = 0.54, R²(θ) = 0.29; r(φ) = 0.73, R²(φ) = 0.53 (filtered). φ preserved exactly; θ median error doubles (~5° → ~10°) from procedural compression. |
| §2.4 | E8 | `06_sec2.4_heritability.ipynb` | Broad-sense heritability (Davis 2025 lme4 model) on θ and φ-derivative traits across three pipelines: gold-standard voxel PCA, PhenoFrame override extractor, and PhenoFrame procedural fit. | gold: 8/13 traits clear H²≥0.20; E1b: 7/13; **E5: 10/13**. Procedural compression preserves or enhances heritability vs both gold and E1b. Top trait `am_theta_e5` H²=0.607. N=44 replicated genotypes. |
| §2.5 | E9 | `07_sec2.5_distribution_match.ipynb` | Population distribution match — fit two generative models (per-node marginal vs per-parameter multivariate Gaussian) to E5 fitted parameters; sample 324 synthetic plants; compare per-node KS, plant-level aggregates, and within-plant cross-node correlations. | Both models pass marginal validation (~67% of param-node pairs at p≥0.05). Only the joint model (B) preserves the real within-plant cross-node correlations (~0.40 for leaf_length and leaf_angle, ~0 for Model A). Plant-aggregate KS: B passes 3 of 4 at p>0.5, A fails on `mean_leaf_angle` due to within-plant variance collapse. |
| §2.5 | E9b | `08_sec2.5_data_driven_generation.ipynb` | End-to-end generative validation: 4-latent-factor stochastic generator (Vigor / Slenderness / Posture / Texture) → descriptor XMLs → **C++ Maize engine** rendering + trait extraction → compare measured θ and consecutive-leaf phyllotaxic deviation against Davis 2025 gold-standard distributions over 324 real plants. | Bulk θ distribution matches well (plant median 30.6° real vs 27.5° synth, KS=0.17; the 3° offset is the known E5 procedural θ bias). Bulk Φ_diff_dev also matches after auto-calibrating PHI_JITTER_SD against the empirical median; the residual secondary mode near 180° in real Φ is the v₂ PCA sign-ambiguity artifact already documented in §2.2 (E1b), not a true biological feature. |
| §2.5 | E9c | `09_sec2.5_holdout_validation.ipynb` | 80/20 held-out validation of E9b: derive config from train-set only, KS- and W₁-compare synthetic plants to the held-out test set, with a measurement-gap baseline (gold vs PhenoFrame on the same test plants). | Plant-aggregate W₁ (deg): synth-vs-gold matches the measurement-gap floor for both φ aggregates (5.41° vs floor 5.93° on mean \|Δφ\|; 8.26° vs 6.23° on median). θ residual is 3-5° above the floor — real generator-vs-population mismatch, not measurement-method confounding. |
| §2.6 | E10 | `10_phase1_gwas_marker_replication.ipynb` | Targeted GWAS marker-effect replication: single-marker linear regression of per-genotype `med` phenotype on minor-allele dosage at Davis 2025's three top phyllotaxy-associated markers (Chr05:12,109,370, Chr05:65,733,791, Chr06:41,390,777) for three phenotype sources (gold, spline, procedural). | All 9/9 (pipeline × marker) effects reproduce Davis's positive direction (p = 0.002 under random-sign null). Procedural pipeline reaches conventional significance (p < 0.05) at Chr05:65,733,791 and Chr06:41,390,777, matching or exceeding the gold-standard control (significant at 1 of 3). All three pipelines miss the published top hit at Chr05:12,109,370 — single-timepoint power loss vs Davis's three-timepoint design. |
| appendix | E1 | `99_appendix_e1_trait_concordance.ipynb` | Re-implementation of Mathieu's PCA trait pipeline; code-correctness check | Appendix-quality; confirms our voxel-PCA code matches the gold standard |

### Not yet started

| ID | What it will test | Dependencies | Notes |
|----|-------------------|--------------|-------|
| E7 | **Human ground truth** — Compare pipeline traits against manual measurements on the same plants. | Manual measurement data (not yet available) | Currently blocked on data collection. (Was originally labeled E6 in earlier drafts; renamed to avoid collision with the completed E6 descriptor-validation notebook.) |

## PhenoFrame source changes for experiments

The following functions were added to `phenoframe/skeleton_to_descriptor.py` to support the experiment pipeline:

- **`extract_traits_via_phenoframe()`** — Converts a voxel skeleton to a PhenoFrame spline representation and extracts θ/φ traits. Used by E1b.
- **`fit_skeleton_procedural()`** — Inverse-fits `(leaf_angle, droopiness, leaf_length)` from a voxel skeleton. Used by E5.
- **`skeleton_to_procedural_xml()`** — Writes a procedural descriptor XML from a skeleton. Parameterized `leaf_width` for visual comparison.
- **`skeleton_to_override_xml()`** — Writes a spline-override descriptor XML (`useCtrlOverrides=1`) from a skeleton. Computes internode distances, per-leaf azimuth, and transforms spline points to the leaf-local coordinate frame.
- **`skeleton_to_point_cloud_obj()`** — Exports a skeleton as a vertex-only OBJ file for visual comparison.

A new module `phenoframe/generator.py` (`phenoframe.MaizeGenerator`) was added to wrap the C++ stochastic generator from the sibling `MaizeProceduralModel` project. It is used by:

- **`derive_sorghum_config.py`** — Calls `MaizeGenerator.derive_config(...)` to write `MaizeProceduralModel/sorghum_generator_config.xml` from E5 + gold-standard population statistics. The wrapper replaces ~65 lines of inline ET-mutation code with a single call and handles the element-name conventions (`leafLengthBase` → `leafLengthScaleCurve`, etc.) automatically.
- **`08_sec2.5_data_driven_generation.ipynb` (§2.5 / E9b)** — Consumes the descriptor XMLs produced by `MaizeGenerator.generate(...)` running the C++ engine on the fitted sorghum config.
- **`09_sec2.5_holdout_validation.ipynb` (§2.5 / E9c)** — Calls both `MaizeGenerator.derive_config(...)` (train-only config) and `MaizeGenerator.generate(...)` (held-out synthetic plants) inline, replacing earlier `subprocess.run([Maize.exe, ...])` invocations. The strict-error semantics of the wrapper surfaced three element-name typos in the original notebook that had been silently no-op'ing (see EXPERIMENT_RESULTS.md § E9c for details).

The wrapper itself is independent of any sorghum-specific assumptions and is documented in the main repo `README.md`.

## Figures and outputs (E1b)

```
figures/
  e1b_scatter_phenoframe_vs_gold.png   # θ and φ scatter: PhenoFrame vs gold
  e1b_bland_altman_by_node.png       # Per-node Δθ and Δφ histograms
  e1b_scatter_filtered.png           # Filtered scatter (voxel count ≥ 20)

outputs/
  e1b_per_leaf_phenoframe_vs_gold.csv   # Per-leaf comparison (all, with keep flag)
  e1b_per_leaf_filtered.csv          # Filtered subset only
  e1b_per_node_summary.csv           # Per-node r, R², RMSE, median |Δ|
  e1b_failures.csv                   # Plants that failed processing
```

## Figures and outputs (E5)

```
figures/
  e5_fit_residuals.png               # RMS error histogram, vs leaf length, by node
  e5_parameter_distributions.png     # Marginal + by-node boxplots for 3 parameters
  e5_scatter_roundtrip_vs_gold.png   # θ and φ scatter: procedural round-trip vs gold
  e5_bland_altman_by_node.png        # Per-node Δθ and Δφ histograms
  e5_scatter_filtered.png            # Filtered scatter (voxel ≥ 20, RMS ≤ 90th pct)

outputs/
  e5_fitted_params.csv               # Per-leaf fitted procedural parameters
  e5_roundtrip_vs_gold.csv           # Per-leaf round-trip comparison (with keep flag)
  e5_roundtrip_filtered.csv          # Filtered subset only
  e5_per_node_summary.csv            # Per-node r, R², RMSE, median |Δ|
  visual_comparison/                 # 3 tiers × 3 files each:
    {good,median,poor}_skeleton.obj  #   Vertex-only point cloud of skeleton
    {good,median,poor}_procedural.obj #  Mesh from fitted procedural descriptor
    {good,median,poor}_override.obj  #   Mesh from spline-override descriptor
```

## Figures and outputs (E4)

```
figures/
  e4_chamfer_distributions.png       # CD and HD histograms + proc vs ovr scatter
  e4_cd_vs_fit_residual.png          # CD vs E5 fit RMS (both descriptors)
  e4_cd_vs_skeleton_size.png         # CD vs skeleton voxel count

outputs/
  e4_geometric_fidelity.csv          # Per-plant CD/HD metrics for both descriptors
```

## Key findings so far

1. **PhenoFrame trait extraction (E1b)** matches the gold standard above the published noise floor for φ. Davis et al. 2025 report a cross-timepoint repeatability of R²(φ) = 0.41, a 3D-vs-manual agreement of R²(φ) = 0.48, and an inter-observer manual ceiling of R²(φ) = 0.55, all for lower-canopy phyllotaxic angles. PhenoFrame achieves R²(φ) = 0.53 on the same skeletons, sitting above the cross-timepoint and method-vs-manual baselines and just below the inter-observer ceiling. R²(θ) = 0.48 is reported standalone (no comparable published baseline exists in Davis 2025 or Tross 2021). Median |error| is not reported in either source paper. The ±180° corner clusters in the φ scatter are v₂ sign ambiguity artifacts, not pipeline errors.

2. **Procedural compression (E5)** preserves azimuth exactly (R²(φ) unchanged) but halves the explained θ variance (R²(θ): 0.48 → 0.29). The ~5° additional median θ error is acceptable for population-level phenotyping but not for individual-leaf precision.

3. **Geometric fidelity (E4)** shows that procedural and override descriptors produce meshes with virtually identical Chamfer distance to the voxel skeleton (median 3.82 cm procedural vs 3.58 cm override). The 3-parameter procedural compression is geometrically lossless at the whole-plant level. The dominant CD error source is the mesh-vs-voxel modality gap (the override-vs-procedural mesh CD is only 0.65 cm, an order of magnitude smaller than either-vs-voxel). Hausdorff distance reveals a ~2 cm procedural penalty at leaf tips (20.03 cm vs 18.04 cm).

4. **Fitted parameter gradients** (leaf_angle increases with node, droopiness magnitude decreases, leaf_length follows a bell curve) are biologically coherent with published sorghum architecture data, supporting their use as heritable traits in E8.

5. **Outlier filtering** (voxel count ≥ 20 for E1b; + RMS ≤ 90th percentile for E5) improves tail statistics without censoring the bulk distribution. E1b drops 2% of leaves; E5 drops 12%.

6. **Python↔C++ pipeline parity (E3)** is confirmed across all 6,213 leaves from 500 plants: length MAE = 0.0001 mm, azimuth MAE = 0.0001°, inclination MAE = 0.0015°. The pure-Python forward model is a faithful reimplementation.

7. **Descriptor round-trip consistency (E6)** validates that XML parameters mean what they claim: `leafLength` matches measured arc length to 0.10 mm, `leafAzimuthDeg` matches measured azimuth to 0.39° (for non-vertical leaves), and cumulative `distance` matches connection-point height to 0.18 mm. The `leafAngle`↔inclination gap is explained by droopiness (r = −0.23). Near-vertical leaves (|leafAngle| > 85°, 6.8% of the dataset) show a 180° azimuth ambiguity that is a geometric property, not a bug.

8. **Heritability preservation (E8).** Running Davis et al. 2025's exact lme4 model (`trait ~ 1 + (1|PI_num)`, n≥2 biological replicates per genotype, reference genotype `PI_656058` excluded) on three trait sources for the same 128 plants and 44 genotypes shows that the 3-parameter procedural compression **preserves or enhances** heritable variation. The procedural-fit pipeline (E5) clears the Davis 2025 published H²≥0.20 threshold for 10 of 13 traits, vs. 8 for the gold-standard voxel-PCA pipeline and 7 for the override-spline extractor. The procedural fit acts as a regularizer that suppresses imaging noise (it recovers `theta_0` H² from 0.11 at E1b to 0.52 at E5, matching gold; and extracts upper-canopy `theta_4` heritability of 0.23 where the gold pipeline finds none). This validates the procedural representation as a viable phenotyping target for downstream GWAS / generative work.

9. **Population distribution match (E9).** Two parametric generative models (independent per-node Gaussians, and per-parameter multivariate Gaussian across the lower 7 nodes) fit to E5's procedural parameters can sample synthetic plant populations that statistically match the real 324-plant SAP subset. Both models reproduce the per-node marginal distributions (~67% of param-node pairs pass a two-sample KS test at p≥0.05) and per-node mean/std curves. The joint multivariate model additionally preserves within-plant cross-node correlations (real ~0.40, joint model 0.41, marginal model ~0) — without this, plant-level aggregates like `mean_leaf_angle` lose 40% of their real population variance. The fitted means + covariance matrices in `outputs/e9_joint_params.csv` are the minimal portable spec for generating arbitrarily many synthetic plants downstream.

10. **End-to-end generative validation (E9b, E9c).** The C++ MaizeGenerator, fitted to the real sorghum population via `derive_sorghum_config.py` and run through `phenoframe.MaizeGenerator`, produces synthetic plants whose plant-level aggregate trait distributions reproduce the bulk of the real population's θ and consecutive-leaf phyllotaxic deviation distributions. Under an 80/20 held-out validation (E9c), Wasserstein-1 distances between synthetic and held-out real plant-aggregate \|Δφ\| distributions match the measurement-gap floor (synth-vs-gold 5.41° vs floor 5.93° on `mean |Δφ|`; 8.26° vs 6.23° on `median |Δφ|`), indicating the φ generator is at the measurement-method floor — its residual is fully explained by the published-pipeline-vs-PhenoFrame disagreement quantified in §2.2. The θ generator has a 3-5° genuine residual above the floor.

11. **GWAS marker-effect replication (E10).** Single-marker linear regression at Davis et al. 2025's three top phyllotaxy-associated markers (Chr05:12,109,370, Chr05:65,733,791, Chr06:41,390,777) across three phenotype sources (gold, spline, procedural) recovers Davis's published positive effect direction at all 9/9 pipeline × marker pairs (p = 0.002 under the null of random sign). The procedural pipeline reaches conventional significance (p < 0.05) at Chr05:65,733,791 (β = +11.1°, p = 0.017) and Chr06:41,390,777 (β = +10.4°, p = 0.019), matching or exceeding the gold-standard positive control which is significant at one of three markers. All three pipelines miss the published top hit at Chr05:12,109,370 — consistent with single-timepoint power loss vs Davis's three-timepoint design. The procedural compression therefore preserves not only heritable variance (§2.4) but the locus-level genetic effects driving that variance.

See [`EXPERIMENT_RESULTS.md`](EXPERIMENT_RESULTS.md) for detailed per-experiment results with tables and figure descriptions.

## Figures and outputs (E3)

```
figures/
  e3_py_vs_cpp_scatter.png           # Python vs C++ scatter for length, azimuth, inclination
  e3_error_distributions.png         # Error distribution histograms
  e3_error_vs_leaf_position.png      # Error vs leaf index

outputs/
  e3_python_cpp_comparison.csv       # Per-leaf Python vs C++ trait values
```

## Figures and outputs (E6)

```
figures/
  length_declared_vs_measured.png    # XML leafLength vs measured arc length scatter
  length_error_histogram.png         # Length error distribution (both pipelines)
  azimuth_declared_vs_measured.png   # XML azimuth vs measured, colored by near-vertical
  angle_vs_inclination.png           # leafAngle vs inclination, colored by droopiness
  droopiness_vs_angle_discrepancy.png # Droopiness vs angle discrepancy
  distance_vs_height.png             # Cumulative distance vs connection-point Y height
  cpp_vs_python_agreement.png        # C++ vs Python cross-validation (3 traits)

outputs/
  descriptor_vs_measured_all_leaves.csv  # Per-leaf XML vs measured values (6,213 rows)
  validation_summary.csv                 # Summary statistics table
```

## Figures and outputs (E8)

```
figures/
  e8_h2_by_trait.png         # Grouped bar chart: gold / E1b / E5 per trait, with H²=0.20 threshold
  e8_h2_vs_gold_scatter.png  # PhenoFrame H² vs gold H² scatter (extractor- and compression-loss view)

outputs/
  e8_heritability.csv        # Long form: 39 rows (13 traits × 3 sources) with σ²_G, σ²_e, H², H²_n2
  e8_h2_by_trait.csv         # Wide form: one row per trait, columns gold/e1b/e5 + extractor/compression loss
  e8_summary.csv             # Per-source counts: traits ≥0.20, ≥0.40, median, mean, max, top trait
```

## Figures and outputs (E9)

```
figures/
  e9_marginal_histograms.png       # Per-node overlaid histograms: real / Model A / Model B
  e9_per_node_meanstd.png          # Per-node mean ± std curves for all 4 parameters
  e9_plant_aggregates.png          # Plant-level aggregate distributions (height, mean angle, etc.)
  e9_cross_node_correlations.png   # Cross-node Pearson r — the key plot showing B preserves what A breaks

outputs/
  e9_marginal_params.csv           # Model A spec: per-node Gaussian mean + std + n
  e9_joint_params.csv              # Model B spec: per-parameter mean vector + N_FIT×N_FIT covariance
  e9_synthetic_modelA.csv          # 324 synthetic plants from Model A (one row per leaf)
  e9_synthetic_modelB.csv          # 324 synthetic plants from Model B
  e9_ks_per_node.csv               # Per-node KS test results (real vs each synthetic model)
  e9_aggregate_compare.csv         # Plant-aggregate KS comparison
  e9_cross_node_correlations.csv   # Cross-node correlation comparison table
  e9_summary.csv                   # Per-model KS pass counts
```

## Figures and outputs (§2.5 / E9b — C++ generator + §2.5 / E9c — held-out)

```
figures/
  sec2_5_per_node_distributions.png      # Real vs synth per-node θ + |Δφ| histograms (main text fig)
  sec2_5_plant_aggregates.png            # Plant-aggregate histograms with KS p-values
  sec2_5_holdout_aggregate_ecdfs.png     # Held-out ECDFs: gold / PhenoFrame / synth, with W₁ vs floor
  sec2_5_holdout_per_node_w1.png         # Per-node W₁ trajectories (synth vs gold vs floor)

outputs/
  sec2_5_synthetic_traits.csv            # Per-leaf measured traits from the C++ generator (E9b)
  sec2_5_ks_per_node.csv                 # Per-node KS comparisons (E9b)
  sec2_5_aggregate_compare.csv           # Plant-aggregate KS comparison (E9b)
  sec2_5_plant_log.csv                   # Per-plant latents + n_leaves (E9b)
  sec2_5_holdout_split.csv               # Train/test PI partition for the held-out experiment (E9c)
  sec2_5_holdout_synthetic_traits.csv    # Held-out synth per-leaf traits (E9c)
  sec2_5_holdout_aggregate_compare.csv   # W₁ + KS, three comparisons per aggregate (E9c)
  sec2_5_holdout_per_node_w1.csv         # Per-node W₁ table (E9c)

External (in sibling MaizeProceduralModel/):
  sorghum_generator_config.xml           # Fitted §2.5 main config (produced by derive_sorghum_config.py)
  sorghum_holdout_config.xml             # Train-only fitted config (produced by notebook 09)
  plant_export_sorghum/                  # 324 descriptor XMLs + OBJs from the §2.5 main run
  plant_export_holdout/                  # 324 descriptor XMLs from the held-out run
```

## Figures and outputs (§2.6 / E10 — GWAS marker replication)

```
figures/
  sec2_6_phase1_marker_replication.png   # 3×3 grid: dosage 0/1/2 boxplots per (pipeline, marker)
  sec2_6_phase1_davis_fig6_replica.png   # Davis Fig 6 C/D/E replica (major vs minor homozygotes)

outputs/
  sec2_6_phase1_gwas_marker_replication.csv  # Per (pipeline, marker) effect β₁, p, R²
  sec2_6_phase1_per_genotype_phenotype.csv   # Per-genotype `med` across the three pipelines
```

## External datasets

The two datasets used here are NOT redistributed in this repo. Get them
yourself once and point the loaders at them via environment variables.

### 1. Jensina's phyllotaxy data

Clone or download:

```
git clone https://github.com/jdavis-132/phyllotaxy.git
```

Then point the loader at it via the environment variable (recommended), or drop
it at the default path resolved by `experiments/paths.py`:

```
export PHENOFRAME_JENSINA_PATH=/path/to/phyllotaxy
```

### 2. Sorghum voxel reconstructions

Download from Zenodo: <https://doi.org/10.5281/zenodo.4426620>
(Gaillard et al., 2021. *Voxel Carving Based 3D Reconstruction of Sorghum*.)

Point the loader at the extracted archive via the environment variable
(recommended), or use the default path resolved by `experiments/paths.py`:

```
export PHENOFRAME_VOXEL_PATH=/path/to/Sorghum
```

The directory must contain `dataset/`, `reconstructed/`, and `skeletons/`.

## Quick check that everything is in place

```python
from experiments import paths, data_loaders as dl

paths.assert_data_present()
print(f"Voxel plants:    {len(dl.list_voxel_plants())}")
print(f"phi CSV plants:  {dl.load_jensina_phi_csv()['plant_name'].nunique()}")
```

If the above prints `Voxel plants: 351` and `phi CSV plants: 977`, you're set.

## Install

The notebooks need the `[experiments]` extras (pandas, matplotlib, jupyter,
openpyxl, scipy, numpy):

```
pip install -e ".[experiments]"
```

## Citation context

These experiments validate PhenoFrame against the published outputs of:

- **Jensina M. Davis et al. (2025)** *Plant Phenomics 7, 100023*
  — phyllotaxy φ traits, heritability, GWAS hits.
- **Tross et al. (2021)** *PeerJ 9:e12628* — leaf insertion angle θ traits.

Both share the same imaging facility, the same 2018 sorghum association
panel, and the same Gaillard et al. voxel-carving + skeletonization
pipeline (`cropsinsilico/SorghumVoxelCarving`).
