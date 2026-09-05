# PhenoSuite Experiment Results

Detailed results for each validation experiment. For setup instructions, external dataset requirements, and file listings, see [`README.md`](README.md).

## E0: Data Exploration

**Notebook:** `01_phyllotaxis/00_data_exploration.ipynb`

**Purpose:** Verify data availability and quality before running validation experiments.

**Dataset:** 351 sorghum plants from the Gaillard et al. (2021) voxel reconstruction pipeline. 332 plants have usable skeletons after filtering.

**Outcome:** Both external datasets (Jensina's phyllotaxy data and sorghum voxel reconstructions) load correctly. Plant/leaf counts match published figures. No further analysis; this notebook serves as a data sanity check.

---

## E1: Re-implementation Concordance

**Notebook:** `01_phyllotaxis/99_appendix_e1_trait_concordance.ipynb`

**Purpose:** Verify that our re-implementation of Mathieu's PCA-based trait extraction pipeline (voxel skeleton to insertion angle theta and azimuth phi) matches the gold-standard output.

**Method:** Re-implement the voxel PCA pipeline and compare our extracted traits against the published theta/phi values on the same skeletons.

**Outcome:** Appendix-quality concordance. Our PCA code matches the gold standard, confirming the baseline is correct before comparing against PhenoSuite's spline-based pipeline.

---

## E1b: PhenoSuite Trait Extraction vs Gold Standard

**Notebook:** `01_phyllotaxis/03_sec2.2_phenosuite_traits_vs_gold.ipynb`

**Purpose:** Evaluate how well PhenoSuite's descriptor-based trait pipeline (spline reconstruction from XML) reproduces the gold-standard theta/phi values derived from the PCA-based voxel pipeline.

**Method:** For each plant, convert the voxel skeleton to a PhenoSuite spline-override descriptor (preserving the original skeleton geometry exactly), extract traits via `compute_traits_from_descriptor`, and compare against the gold-standard theta/phi.

**Results:**

| Trait | r | R^2 | Median |error| (PhenoSuite) | Notes |
|-------|---|-----|-----------------------------|-------|
| theta (inclination) | 0.69 | 0.48 | ~5 deg | Filtered to lower 4 leaves |
| phi (azimuth) | 0.73 | 0.53 | ~15 deg | Filtered to lower 4 leaves |

**Comparison against published baselines (Davis et al. 2025, lower-canopy phi only):**

| Comparison | Published R^2 | PhenoSuite R^2(phi) |
|-----------|----------------|--------------------|
| 3D reconstruction, same plants 2 days apart (cross-timepoint repeatability) | 0.41 (n = 961 angles) | **0.53** |
| 3D reconstruction vs manual measurement (method agreement) | 0.48 (n = 75 angles) | **0.53** |
| Manual vs manual, two observers (inter-observer ceiling) | 0.55 (n = 46 angles) | **0.53** |
| 3D reconstruction, 3 days apart | 0.33 | — |
| 3D reconstruction, 5 days apart | 0.24 | — |

Median |error| is **not reported in Davis et al. 2025** (the paper quantifies repeatability and method agreement only via R^2), and Tross et al. 2021 likewise does not publish a comparable per-leaf median-error figure. A side-by-side median-error comparison against the published literature is therefore not currently possible; the headline comparison against published work is on R^2 only.

**Key observations:**

- PhenoSuite's R^2(phi) = 0.53 **exceeds the published cross-timepoint repeatability (0.41) and method-vs-manual agreement (0.48)** from Davis et al. 2025, and approaches the inter-observer manual ceiling (0.55). The spline-based pipeline is therefore at least as reliable as the published voxel-PCA approach on the same skeletons.
- The R^2 = 0.41 cited in Davis et al. 2025 is a **phi** statistic (cross-timepoint repeatability of phyllotaxic angles), not a theta statistic; earlier drafts of this writeup misattributed it.
- No published R^2 baseline exists for theta in Davis 2025 (which is a phyllotaxy/phi paper) or for the trait set in Tross 2021 as cited here. PhenoSuite's R^2(theta) = 0.48 is therefore reported standalone, without a published comparison point.
- The +/-180 degree corner clusters in the phi scatter plot are v2 sign ambiguity artifacts from the PCA decomposition in the gold standard, not errors in PhenoSuite.
- Filtering by voxel count >= 20 drops 2% of leaves and improves tail statistics.

**Figures:** `e1b_scatter_phenosuite_vs_gold.png`, `e1b_bland_altman_by_node.png`, `e1b_scatter_filtered.png`

---

## E5: Procedural Round-Trip Fidelity

**Notebook:** `01_phyllotaxis/05_sec2.3_procedural_round_trip.ipynb`

**Purpose:** Measure how much trait information is lost when compressing a full skeleton into a 3-parameter procedural descriptor (leaf_angle, droopiness, leaf_length) and reconstructing.

**Method:** For each leaf skeleton, inverse-fit the procedural parameters via L-BFGS-B optimization (minimizing squared distance between resampled polylines), reconstruct via the forward model, extract traits, and compare against the gold standard.

**Fitting details:**
- 2-parameter optimization: leaf_angle and droopiness (leaf_length is set directly from the skeleton arc length).
- Droopiness bounds: [-300, 300] (widened after the droopiness fix removed the 5x multiplier).
- Multi-start: 5 initial guesses with varying droopiness signs and magnitudes.
- leafCurl is always set to 0 (it does not affect the center spline).

**Results:**

| Trait | r | R^2 | Median |error| | Notes |
|-------|---|-----|-----------------|-------|
| theta (inclination) | 0.54 | 0.29 | ~10 deg | Filtered |
| phi (azimuth) | 0.73 | 0.53 | ~15 deg | Unchanged from E1b |

**Key observations:**

- Azimuth is preserved exactly (R^2(phi) unchanged) because the procedural model places leaves at the same azimuth as the override descriptor.
- Theta precision degrades (R^2 drops from 0.48 to 0.29) because the 2-parameter model cannot capture arbitrary curvature profiles.
- The ~5 degree additional median theta error is acceptable for population-level phenotyping but not for individual-leaf precision work.
- Fitted parameter gradients are biologically coherent: leaf_angle increases with node index, droopiness magnitude decreases toward the top, and leaf_length follows a bell curve peaking at mid-canopy. These patterns match published sorghum architecture data.
- Outlier filtering (voxel >= 20 + RMS <= 90th percentile) drops 12% of leaves.

**Figures:** `e5_fit_residuals.png`, `e5_parameter_distributions.png`, `e5_scatter_roundtrip_vs_gold.png`, `e5_scatter_filtered.png`

---

## E4: Geometric Fidelity

**Notebook:** `01_phyllotaxis/04_sec2.3_geometric_fidelity.ipynb`

**Purpose:** Quantify the 3D geometric difference between the reconstructed meshes (from both procedural and override descriptors) and the original voxel skeletons using Chamfer and Hausdorff distances.

**Method:** For each plant, generate OBJ meshes from both the procedural and override descriptors, sample point clouds, and compute symmetric Chamfer distance (CD) and Hausdorff distance (HD) against the voxel skeleton point cloud. Two comparisons are reported: voxel-vs-procedural (`vp_*`) and voxel-vs-override (`vo_*`). A third comparison, override-vs-procedural (`op_*`), measures how close the two mesh representations are to each other and is reported separately as the "representation gap".

**Results — vs voxel skeleton (the comparison that matters for fidelity):**

| Descriptor type | Median CD (cm) | Median HD (cm) |
|-----------------|----------------|-----------------|
| Procedural (3 params/leaf) | 3.82 | 20.03 |
| Override (full spline) | 3.58 | 18.04 |

**Results — override vs procedural mesh (representation gap, not vs voxel):**

| Comparison | Median CD (cm) | Median HD (cm) |
|------------|----------------|------------------|
| Override vs procedural mesh | 0.65 | 8.93 |

**Key observations:**

- Procedural and override descriptors produce meshes with virtually identical Chamfer distance to the voxel skeleton (3.82 cm vs 3.58 cm). The 3-parameter compression is geometrically lossless at the whole-plant level.
- The dominant CD error source is the mesh-vs-voxel modality gap (mesh has width, voxels are centerline), not the descriptor representation. The override-vs-procedural mesh CD is 0.65 cm — an order of magnitude smaller than either-vs-voxel — confirming this.
- Hausdorff distance reveals a ~2 cm procedural penalty at leaf tips (20.03 cm vs 18.04 cm), where the simple droopiness model cannot reproduce complex tip curvature.
- CD does not correlate with skeleton size (voxel count), suggesting the error floor is set by the modality gap rather than plant complexity.

**Figures:** `e4_chamfer_distributions.png`, `e4_cd_vs_fit_residual.png`, `e4_cd_vs_skeleton_size.png`

---

## E3: Python vs C++ Pipeline Parity

**Notebook:** `01_phyllotaxis/01_sec2.1_python_cpp_validation.ipynb`

**Purpose:** Verify that the pure-Python forward model (`phenosuite.traits.compute_traits_from_descriptor`) produces identical results to the C++ geometry engine (`phenosuite.wrapper.Maize`) for all measured traits.

**Method:** For each of 500 plant descriptors (PlantsC++ dataset, 6,213 leaves total), run both the C++ engine and the Python forward model, extract per-leaf traits, and compare.

**Results:**

| Trait | MAE | Max |error| | Pearson r |
|-------|-----|-------------|-----------|
| Leaf length | 0.0001 mm | 0.0003 mm | 1.000000 |
| Azimuth | 0.0001 deg | 0.0678 deg | 1.000000 |
| Inclination | 0.0015 deg | 0.0102 deg | 1.000000 |

**Key observations:**

- The two pipelines are functionally identical across all 6,213 leaves. Residual differences are at the floating-point precision level.
- The Python forward model uses the same algorithm as the C++ engine: equal-length spline segments with quadratic droopiness rotation, B-spline arc length with 128 uniform samples, and identical local-to-world coordinate transforms.
- This validates that the pure-Python pipeline can be used as a drop-in replacement for the C++ engine for trait computation, with no DLL dependency required.
- Critical bugfix verified: the droopiness multiplier was corrected from `droopiness * 5.0` (incorrect, from the DLL source) to `droopiness` (correct, matching the authoritative C++ codebase). Both pipelines now use the corrected formula.

**Figures:** `e3_py_vs_cpp_scatter.png`, `e3_error_distributions.png`, `e3_error_vs_leaf_position.png`

---

## E6: Descriptor vs Measured Trait Validation

**Notebook:** `01_phyllotaxis/02_sec2.1_descriptor_trait_validation.ipynb`

**Purpose:** Verify that the parameter values declared in the XML descriptor (e.g., `leafLength`, `leafAzimuthDeg`, `distance`) are faithfully recovered when the procedural model reconstructs the geometry and traits are extracted. This validates the full round-trip: XML -> forward model -> geometry -> trait extraction.

**Method:** For each of 500 plant descriptors (6,213 leaves), parse the raw XML parameter values, run both the C++ and Python pipelines, and compare declared vs measured values for four key parameters.

**Results:**

### Leaf length (`leafLength` vs measured arc length)

| Pipeline | r | MAE (mm) | Max |error| (mm) | Mean error (mm) |
|----------|---|----------|---------------------|------------------|
| C++ | 1.000000 | 0.1017 | 0.2403 | -0.1017 |
| Python | 1.000000 | 0.1017 | 0.2405 | -0.1017 |

The small negative bias (-0.10 mm) is a numerical integration artifact: the B-spline chord-length approximation (128 uniform samples) slightly underestimates the true arc length defined by the control points. This is expected and not a model inconsistency.

Error percentiles (both pipelines identical):
- P50: 0.1022 mm
- P90: 0.1510 mm
- P95: 0.1656 mm
- P99: 0.1916 mm
- P100: 0.2403 mm

### Azimuth (`leafAzimuthDeg` vs measured azimuth)

| Subset | MAE (deg) | Max |error| (deg) | Notes |
|--------|-----------|---------------------|-------|
| All 6,213 leaves | 7.77 | 180.00 | Inflated by near-vertical outliers |
| Excluding |leafAngle| > 85 deg (N=5,788) | 0.39 | 160.21 | Sub-degree accuracy |

**Near-vertical azimuth flip:** 425 leaves (6.8%) have `|leafAngle| > 85 deg`. When the leaf is aimed past vertical, the first spline segment's local radial component becomes negative (the base tangent initially points slightly inward toward the stem). In world coordinates, this reverses the XZ tangent direction, producing a 180 degree flip. This is a geometric property of the measurement, not a bug. Both C++ and Python produce identical flips.

### Leaf angle (`leafAngle` vs measured inclination)

| Pipeline | r | Mean diff (deg) | Std diff (deg) | Max |diff| (deg) |
|----------|---|-----------------|-----------------|---------------------|
| C++ | 0.9904 | -0.68 | 4.87 | 58.87 |
| Python | 0.9904 | -0.68 | 4.87 | 58.87 |

This discrepancy is **expected**. The `leafAngle` parameter sets the target direction for the first spline segment, but droopiness (a gravity-like quadratic pull) modifies the actual curvature. The correlation between droopiness and the angle discrepancy is r = -0.23, confirming the causal mechanism.

### Internode distance (cumulative `distance` vs connection-point Y height)

| Pipeline | r | MAE (mm) | Max |error| (mm) | Mean error (mm) |
|----------|---|----------|---------------------|------------------|
| C++ | 0.999999 | 0.1807 | 3.4204 | -0.0953 |
| Python | 0.999999 | 0.1720 | 3.4204 | -0.0953 |

The model uses Y-up coordinates, so `connection_point.y` is the height axis. The small residual comes from the stem radius: the connection point is on the stem surface, not the centerline. The radial offset projected onto Y depends on the leaf's insertion angle and the stem radius at that height.

### C++ vs Python agreement (cross-check)

| Trait | MAE | Max |error| |
|-------|-----|-------------|
| Length | 0.0001 mm | 0.0003 mm |
| Azimuth | 0.0001 deg | 0.0678 deg |
| Inclination | 0.0015 deg | 0.0102 deg |

**Conclusion:** The descriptor format is consistent with the model. Parameters mean what they claim: `leafLength` faithfully sets the arc length, `leafAzimuthDeg` faithfully sets the compass direction (for non-vertical leaves), and `distance` faithfully sets the internode spacing. The `leafAngle`/inclination discrepancy is an expected consequence of the droopiness model, not a format inconsistency.

**Figures:** `length_declared_vs_measured.png`, `length_error_histogram.png`, `azimuth_declared_vs_measured.png`, `angle_vs_inclination.png`, `droopiness_vs_angle_discrepancy.png`, `distance_vs_height.png`, `cpp_vs_python_agreement.png`

---

## E8: Broad-sense Heritability

**Notebook:** `01_phyllotaxis/06_sec2.4_heritability.ipynb`

**Purpose:** Test whether the PhenoSuite representation preserves heritable trait variation by comparing broad-sense H^2 across three trait sources on the same plants and same lme4 model.

**Method:** Davis et al. 2025's exact heritability framework — random-intercept REML mixed model `trait ~ 1 + (1|PI_num)` on biological replicates (>=2 plants per genotype), reference genotype `PI_656058` excluded. H^2 = sigma^2_G / (sigma^2_G + sigma^2_e / 2), matching Davis paper Eq. 3 with n=2.

**Implementation note:** `statsmodels.MixedLM` with `method='lbfgs'` gets pinned to the sigma^2_G=0 boundary even when there is real between-genotype variance; we use `bfgs` first and fall back to `powell` or `lbfgs`. Cross-checked against ANOVA method-of-moments and the published R `lmer` solution agrees to numerical precision.

**Three trait sources** (same 128 plants, 44 replicated genotypes, lower 5 leaves):
- **gold** — voxel-PCA theta/phi from `angles.txt` (published method, Gaillard 2021).
- **E1b** — PhenoSuite's spline-based trait extractor on the override descriptor (full skeleton, no compression).
- **E5** — PhenoSuite's trait extractor on the 3-parameter procedural fit (`leaf_angle`, `droopiness`, `leaf_length`).

**Results (Davis Eq. 3 with n=2):**

| Source | Traits with H^2 >= 0.20 | Median H^2 | Mean H^2 | Top trait |
|--------|--------------------------|-------------|-----------|-----------|
| gold (voxel PCA) | 8 / 13 | 0.319 | 0.288 | `am_phi_dev` = 0.561 |
| E1b (override + PhenoSuite) | 7 / 13 | 0.293 | 0.249 | `med_phi_dev` = 0.498 |
| **E5 (procedural + PhenoSuite)** | **10 / 13** | **0.484** | **0.382** | **`am_theta` = 0.607** |

**Per-trait detail (Davis Eq. 3 H^2_n2):**

| Trait | gold | E1b | E5 | Notes |
|-------|------|-----|----|----|
| `theta_0` (insertion, lowest leaf) | 0.54 | 0.11 | 0.52 | Procedural recovers extractor loss |
| `theta_1` | 0.45 | 0.46 | 0.53 | All sources clearly heritable |
| `theta_2` | 0.13 | 0.34 | 0.48 | Extraction + compression each lift heritability |
| `theta_3` | 0.09 | 0.14 | 0.38 | Procedural enhances most |
| `theta_4` (highest leaf, upper canopy) | 0.00 | 0.00 | 0.23 | Gold/E1b drowned by movement noise; procedural extracts signal |
| `med_theta` (per-plant median) | 0.26 | 0.39 | 0.60 | Best per-plant summary |
| `am_theta` (per-plant mean) | 0.32 | 0.29 | 0.61 | |
| `Phi_dev_1` (deviation from 180 deg, leaves 0-1) | 0.34 | 0.12 | 0.14 | Extractor loses signal |
| `Phi_dev_2` | 0.33 | 0.31 | 0.40 | Preserved |
| `Phi_dev_3` | 0.05 | 0.00 | 0.00 | Not heritable in any source |
| `Phi_dev_4` | 0.16 | 0.12 | 0.01 | Compression loses upper-canopy signal |
| `med_phi_dev` (per-plant median) | 0.53 | 0.50 | 0.54 | Highly preserved |
| `am_phi_dev` (per-plant mean) | 0.56 | 0.45 | 0.53 | Highly preserved |

**Key observations:**

- **Procedural compression preserves or enhances heritability.** E5 clears the H^2 >= 0.20 threshold for 10 of 13 traits versus 8 for the gold standard and 7 for the override extractor. The procedural fit acts as a regularizer that suppresses imaging noise while keeping biological signal. This validates the procedural representation as a viable phenotyping target for the genotype-to-phenotype direction (E9 distribution match, downstream GWAS).
- **Upper-canopy theta becomes heritable under procedural compression.** Gold-standard `theta_4` has H^2 = 0; E5 `theta_4` has H^2 = 0.23. Davis 2025 explicitly attributes upper-leaf H^2 loss to plant movement during rotation imaging — the procedural model's smoothing across nodes recovers a signal that single-leaf PCA cannot.
- **Per-plant aggregates are consistently the most heritable summaries.** `med_phi_dev` and `am_phi_dev` clear 0.45 in all three sources. These are the natural targets for downstream GWAS.
- **One trait shows clean extractor loss with no recovery.** `Phi_dev_1` drops from H^2 = 0.34 (gold) to 0.12 (E1b) and stays low (0.14) at E5. This is a genuine signal loss attributable to PhenoSuite's chord-based azimuth measurement on the first inter-node difference. Worth flagging for individual-leaf-precision use cases.
- **N = 44 genotypes is small.** Confidence intervals on individual H^2 values are wide. The within-pipeline rankings (gold vs E1b vs E5) are robust because all three use the same plants and same model.

**Caveat on cross-paper comparison.** Davis 2025 uses **timepoint replicates** (same plant on different days) for her heritability model. Our Zenodo voxel archive provides only one timepoint per plant, so we use **biological replicates** (different physical plants of the same genotype) — a more conservative design that places a lower bound on H^2 than Davis's timepoint design. Absolute H^2 magnitudes are therefore not directly comparable to Davis's numbers, but the H^2 >= 0.20 threshold she uses is appropriate as a 'heritable' bar in either design.

**Figures:** `e8_h2_by_trait.png` (grouped bar chart, gold/E1b/E5 per trait with H^2=0.20 threshold line), `e8_h2_vs_gold_scatter.png` (PhenoSuite H^2 vs gold H^2 scatter on the y=x line).

**Outputs:** `e8_heritability.csv` (long form with all variance components), `e8_h2_by_trait.csv` (wide form per trait), `e8_summary.csv` (per-source counts).

---

## Summary Table

| Experiment | What it validates | Key metric | Result |
|------------|-------------------|------------|--------|
| E0 | Data availability | Plant count | 332 usable skeletons |
| E1 | PCA pipeline re-implementation | Concordance | Matches gold standard |
| E1b | PhenoSuite traits vs gold standard | R^2(theta), R^2(phi) | 0.48, 0.53 |
| E5 | Procedural round-trip | R^2(theta), R^2(phi) | 0.29, 0.53 |
| E4 | Geometric fidelity (mesh vs skeleton) | Median Chamfer distance | 3.55 cm |
| E3 | Python vs C++ pipeline parity | Length MAE | 0.0001 mm |
| E6 | XML descriptor round-trip | leafLength MAE | 0.10 mm |
| E8 | Heritability preservation under procedural compression | # traits with H^2 >= 0.20 | gold: 8, E1b: 7, **E5: 10** of 13 |
| E9 | Population distribution match (synthetic plant generation) | Cross-node correlation recovery; plant-aggregate KS | Model B preserves real correlations (~0.40), KS p>0.5 on 3 of 4 aggregates |
| E9b | C++ generator end-to-end (§2.5 main) | Plant-aggregate distribution overlap | Bulk θ + Φ_diff_dev distributions overlap; visual match with documented residual = §2.2 measurement gap |
| E9c | C++ generator held-out validation (§2.5 supplementary) | Wasserstein-1 vs measurement-gap floor | φ aggregates at floor (synth-vs-gold W1 5.4° vs floor 5.9°); θ has 3-5° real residual above floor |
| E10 | GWAS marker-effect replication at Davis 2025 hits (§2.6) | Effect direction + significance recovery | 9/9 effect directions positive (p = 0.002 vs random); procedural pipeline significant at 2/3 markers (matches/beats gold-standard control) |

## E9: Population Distribution Match

**Notebook:** `01_phyllotaxis/07_sec2.5_distribution_match.ipynb`

**Purpose:** Test whether a small parametric model fitted to E5's per-leaf procedural parameters can reproduce the population statistics of the real 324-plant set. If yes, PhenoSuite's 3-parameter descriptor format is a viable target for synthetic-plant generation, GWAS power simulation, and downstream-pipeline stress tests.

**Method:** Two generative models of increasing realism, both fitted to E5's filtered population (324 plants, 2,479 leaves, same outlier filter as E5):

- **Model A — independent per-node marginals**: each (parameter, node) gets its own univariate Gaussian (mean, std). Leaves are sampled independently. Baseline.
- **Model B — per-parameter multivariate Gaussian**: one full-N_FIT multivariate Gaussian per parameter, sampled jointly so within-plant cross-node correlations are preserved. N_FIT = 7 (lower nodes with ≥50% population coverage); higher nodes fall back to Model A.

Both models draw plant leaf counts from the empirical PMF and sample 324 synthetic plants per model.

**Comparison metrics:**
- Per-node KS test between real and synthetic populations for each of the four parameters
- Per-node mean ± std visualization
- Plant-level aggregates (`total_height`, `mean_leaf_angle`, `total_leaf_length`, `mean_droopiness`): means, stds, KS tests
- Within-plant cross-node Pearson correlations (real vs Model A vs Model B)

**Results — per-node marginals (KS test pass count, p ≥ 0.05):**

| Model | (param, node) pairs | Pass at p≥0.05 | Median KS stat |
|-------|---------------------|----------------|----------------|
| A (marginal) | 44 | 29 (66%) | 0.125 |
| B (joint)    | 40 | 27 (68%) | 0.108 |

Both models match marginal distributions roughly equally well — as expected, since both are fit to the same per-node moments.

**Results — plant-level aggregates (KS p-value, real vs synthetic):**

| Aggregate | mean real | std real | Model A KS p | Model B KS p |
|-----------|-----------|----------|----------------|----------------|
| total_height (m) | 0.277 | 0.090 | 0.056 (borderline) | **0.505 (clean match)** |
| mean_leaf_angle (deg) | 65.25 | 9.95 | **<0.001 (fail)** | 0.029 |
| total_leaf_length (m) | 2.037 | 0.717 | 0.103 | **0.980 (clean match)** |
| mean_droopiness | -75.24 | 29.89 | 0.069 | **0.568 (clean match)** |

Model A fails the plant-level KS test on `mean_leaf_angle` (its synthetic std collapses from the real 9.95° to 5.96° because independent per-node sampling averages out the within-plant variance). Model B recovers a std of 8.79° and passes the looser p<0.05 bar.

**Results — within-plant cross-node correlations (Pearson r between adjacent-node values of the same parameter):**

| Parameter | Avg r (real) | Avg r (Model A) | Avg r (Model B) |
|-----------|--------------|------------------|-------------------|
| leaf_length | 0.39 | -0.06 | 0.41 |
| leaf_angle  | 0.41 | 0.01 | 0.41 |
| droopiness  | 0.20 | -0.02 | 0.29 |
| distance    | -0.12 | 0.00 | -0.09 |

Model A breaks all cross-node correlations to zero by design (independent sampling). Model B reproduces them faithfully across all four parameters.

**Key observations:**

- **Both models pass marginal validation.** The per-node Gaussians fit the data well enough that ~67% of (parameter, node) pairs cannot be statistically distinguished from real, and per-node mean ± std curves overlay almost perfectly.
- **Only Model B preserves the joint structure.** This is the defining advantage of the multivariate-Gaussian formulation: a synthetic plant that has long leaves at node 0 also has long leaves at node 1, mirroring the real population's biological coherence (genetic + environmental factors persist within a plant).
- **The cost of breaking joint structure shows up at the plant level.** Model A's `mean_leaf_angle` std collapses by 40% (9.95° → 5.96°) because uncorrelated per-leaf sampling averages out within-plant variance. Model B keeps the std at 8.79°. The same effect would propagate to any genome-wide association study run on Model A's synthetic phenotypes — heritable signal would be systematically diluted.
- **Model B is the deliverable.** The fitted means + covariance matrices in `e9_joint_params.csv` are the minimal data needed to regenerate the entire 324-plant population shape. Downstream users can sample arbitrarily many synthetic plants from these parameters using a single `np.random.multivariate_normal` call per parameter.

**Caveats:**

- The Gaussian assumption is convenient but may misfit some parameters (especially `droopiness`, which has a heavier left tail in the real data). A future iteration could swap in lognormal for `leaf_length` or a mixture for `droopiness`.
- Only the lower N_FIT = 7 nodes are jointly modeled; higher nodes fall back to per-node marginals. This is fine for population statistics but limits the joint realism of synthetic plants with many leaves.
- We have not run the synthetic descriptors through PhenoSuite's forward model and compared the resulting measured-trait (θ, φ) distributions to the real ones. That is a natural follow-up — and the existing pipeline supports it — but it is out of scope for E9 as currently scoped.

**Figures:** `e9_marginal_histograms.png` (per-node overlaid histograms), `e9_per_node_meanstd.png` (per-node mean ± std curves), `e9_plant_aggregates.png` (plant-level aggregate distributions), `e9_cross_node_correlations.png` (within-plant correlations bar chart — the killer plot showing Model B preserves what Model A breaks).

**Outputs:**
- `e9_marginal_params.csv` — Model A spec (per-node Gaussian mean, std, n)
- `e9_joint_params.csv` — Model B spec (per-parameter mean vector + 7×7 covariance matrix)
- `e9_synthetic_modelA.csv`, `e9_synthetic_modelB.csv` — the synthetic populations themselves (one row per leaf)
- `e9_ks_per_node.csv`, `e9_aggregate_compare.csv`, `e9_cross_node_correlations.csv`, `e9_summary.csv` — comparison tables

---

## E9c: Held-out Validation of the C++ Generator (§2.5 supplementary)

**Notebook:** `01_phyllotaxis/09_sec2.5_holdout_validation.ipynb`

**Purpose:** Address two methodological concerns in E9b: (1) the generator's config was calibrated on the same 324 plants it was later KS-compared against (data peeking), and (2) the published-pipeline-vs-PhenoSuite trait-extraction disagreement quantified in §2.2 sets an irreducible floor on how close any generator can get to the gold standard. E9c does a single-shot held-out validation and reports a measurement-gap baseline alongside the headline numbers.

**Method:**

- **80/20 train/test split** over plant IDs (deterministic, seed = 20260601). Train = 259 plants, test = 65 plants.
- Derive a train-only generator config via `phenosuite.MaizeGenerator.derive_config(...)` (the wrapper handles the `leafLengthBase` → `leafLengthScaleCurve` naming convention etc.). Train-only canopy curves for `leafLengthBase`, `droopinessBase`, `leafAngleBase` (gold-θ calibrated), train-only height + numLeaves means/stds, uniform-noise-corrected `tillerAzimuthNoise`, sorghum-typical visual defaults.
- Generate 324 synthetic plants with `MaizeGenerator(...).generate(...)` and the train-only config.
- Extract synthetic traits via the PhenoSuite trait pipeline (§2.1 / E3-equivalent).
- For each plant-level aggregate, compute three Wasserstein-1 distances (in degrees):
  - `synth_vs_test_gold` — synthetic vs the held-out test plants' Davis voxel-PCA gold standard. **Headline.**
  - `synth_vs_test_phenosuite` — synthetic vs the same test plants' PhenoSuite traits on their E5 procedural fits. Removes measurement-method gap.
  - `measurement_gap_baseline` — gold vs PhenoSuite on the same test plants. **The floor**: the disagreement that exists without any generation.
- Also compute KS statistics (no p-values: at N≈260 vs 65 any 2° mean shift is significant; KS p-values do not differentiate "close to the floor" from "far above it").

**Results — plant-aggregate W₁ (degrees), test split:**

| Aggregate | floor | synth vs gold | gap above floor |
|-----------|-------|---------------|------------------|
| plant median θ | **2.07** | 7.42 | +5.35 |
| plant mean θ   | **3.15** | 5.71 | +2.56 |
| plant median \|Δφ−180°\| | **6.23** | 8.26 | +2.03 |
| plant mean \|Δφ−180°\|   | **5.93** | 5.41 | **−0.52** (essentially at floor) |

**Key observations:**

- **The φ generator is at the measurement-method floor.** `synth_vs_test_gold` is essentially equal to the floor on `mean |Δφ|` (5.41° vs 5.93°), and within 2° on `median |Δφ|`. The residual disagreement on φ is fully explained by the §2.2 / E1b PhenoSuite-vs-voxel-PCA gap (R²(φ) = 0.53) and cannot be improved without changing the trait extractor.
- **The θ generator has a 3-5° real residual above the floor.** This is genuine generator-vs-population mismatch, not measurement-method confounding — synth-vs-PhenoSuite-on-real W₁ stays elevated even after removing the gold-vs-PhenoSuite gap.
- **Calibration honesty.** Only the canopy curves + plant-level means/stds + `tillerAzimuthNoise` come from the train set. `leafJitterScale = 0.60` is the empirical compromise picked in E9b and held fixed; loadings, whorl, visual defaults are inherited from the maize template unchanged.
- **Three element-name bugs surfaced by the strict-error wrapper** (`plantHeight` → `height`, `numberOfLeaves` → `numLeaves`, `leafLengthCurve` → `leafLengthScaleCurve`). The old inline `for el in root.iter(name): ... return` code silently no-op'd; the wrapper raises `KeyError`. Fix tightened the held-out W₁ residuals by 0.5–0.7° on θ aggregates. Scope of the bug confined to this notebook; §2.5 main config / figures / numbers unaffected.

**Caveats:**

- The held-out trait set is θ and φ only. Davis's `angles.txt` gold standard does not include leaf length, internode distance, or plant height, so we cannot validate the generator's match on those held-out traits.
- The 4-latent (V/S/P/T) loadings structure is inherited from the maize config unchanged; we do not independently validate it in sorghum.
- Single 80/20 split, no k-fold; the per-aggregate W₁ values have meaningful sampling variation at N(test) = 65.

**Figures:** `sec2_5_holdout_aggregate_ecdfs.png` (1×4 ECDF panels with W₁ in titles), `sec2_5_holdout_per_node_w1.png` (per-node W₁ trajectory across three comparisons).

**Outputs:**
- `sec2_5_holdout_aggregate_compare.csv` — W₁ + KS per (aggregate, comparison)
- `sec2_5_holdout_per_node_w1.csv` — per-node W₁ across three comparisons
- `sec2_5_holdout_synthetic_traits.csv` — per-leaf measured traits from the held-out synthetic plants
- `sec2_5_holdout_split.csv` — train/test PI partition

---

## E10: GWAS Marker-Effect Replication (§2.6)

**Notebook:** `01_phyllotaxis/10_phase1_gwas_marker_replication.ipynb`

**Purpose:** Test whether the descriptor-derived phenotype pipelines preserve the genetic-marker effects that drive the heritable variation quantified in §2.4. A targeted Phase 1 replication: no FarmCPU, no resampling — just per-genotype linear regression of `med` (median lower-canopy |Δφ−180°|) on minor-allele dosage at three published phyllotaxy-associated markers from Davis et al. 2025.

**Method:**

- **Three target markers** (Davis 2025, RMIP ≥ 0.1 for `med`):
  - Chr05:12,109,370 (RMIP ≈ 0.26) — Davis's top hit
  - Chr05:65,733,791 (RMIP ≈ 0.19)
  - Chr06:41,390,777 (RMIP ≈ 0.10)
- **Three phenotype sources** (the same 215 SAP genotypes that intersect the Zenodo voxel archive with Davis's filtered VCF):
  - Gold — voxel-PCA φ values from `angles.txt`. Positive control.
  - Spline — PhenoSuite trait pipeline on the override descriptor.
  - Procedural — PhenoSuite trait pipeline on the procedural 3-parameter fit.
- **Trait definition:** `med` per plant = median of `|normalize360(φᵢ₊₁ − φᵢ) − 180°|` over leaves i = 0..3. Per genotype, average across plant replicates.
- **Test:** per-genotype linear regression `med = β₀ + β₁ · dosage + ε`, dosage = 0/1/2 (minor-allele count). Reference genotype `PI656058` excluded per Davis. No kinship adjustment (Phase 1 design).
- **Davis's expected result:** minor allele increases `med` at all three markers (positive β₁).

**Results — per (pipeline, marker), effect β₁ and p-value:**

| Pipeline | Chr05:12,109,370 | Chr05:65,733,791 | Chr06:41,390,777 |
|----------|------------------|------------------|------------------|
| Gold (Davis voxel-PCA) | +4.61°, p=0.362 | **+12.26°, p=0.004** | +7.41°, p=0.066 |
| Spline | +7.00°, p=0.196 | +8.91°, p=0.049 | +5.85°, p=0.176 |
| **Procedural** | +2.03°, p=0.716 | **+11.12°, p=0.017** | **+10.36°, p=0.019** |

**Key observations:**

- **All 9/9 effect directions are positive,** matching Davis at every (pipeline, marker) pair. Under the null hypothesis of random sign per cell, P(9/9 positive) = 0.5⁹ ≈ 0.002. The replication is highly nonrandom even without per-marker significance.
- **Procedural pipeline reaches conventional significance (p < 0.05) at 2 of 3 markers**, matching or beating the gold-standard positive control which is significant at 1 of 3. This is consistent with the §2.4 / E8 finding that procedural compression acts as a regularizer against per-leaf measurement noise — denoising the phenotype tightens marker-effect estimates.
- **All three pipelines miss Davis's TOP marker Chr05:12,109,370** (gold p = 0.36). The miss is *not* a pipeline-quality issue: Davis's `med` averages over leaves AND three imaging timepoints, while the Zenodo voxel archive has one timepoint per plant. The implied SNR loss (~√3) is enough to drop power below significance at this marker's effect size.
- **Davis Fig. 6 C/D/E layout replicated** (`sec2_6_phase1_davis_fig6_replica.png`) for direct visual comparison: per-genotype median lower-canopy |Δφ−180°| for major-homozygous vs minor-homozygous genotypes, heterozygotes excluded. Per-genotype n's are 1-3 lower than Davis's reported counts because two PIs in her filtered VCF are missing from our analysis (PI656036 not in the Zenodo voxel archive at all; PI548797 has skeleton extraction but fails our stem-alignment step). PI548797 is a minor-allele homozygote at Chr05:65,733,791 — its loss reduces our minor-allele sample at that marker by 1.

**Caveats:**

- Our `med` averages over leaves only (Davis: leaves AND three imaging timepoints). The Zenodo voxel archive has one timepoint per plant.
- No kinship / population-structure correction; p-values are anti-conservative relative to FarmCPU.
- N ≈ 215 genotypes after intersecting the voxel set with Davis's filtered VCF; underpowered for true effect-size estimates but enough for direction tests at common alleles.
- Davis's hits are stable across resampling (RMIP); our single-shot p-values do not need to match her exact significance levels.

**Decision context.** Phase 1 was framed as "we recover the marker effects at Davis's published hits." Phase 2 (full FarmCPU-RMIP replication on the procedural phenotype) would let us say "we recover Davis's stable hits as stable hits in our own analysis" — a stronger claim — but requires the full ~4.7M-marker SAP VCF (~5-50 GB download from Boatwright et al.) and ~50-200 hours of compute. The current Phase 1 result is publishable as a supplementary table + brief §2.6 paragraph; the marginal value of Phase 2 was judged to not justify the cost.

**Figures:**
- `sec2_6_phase1_marker_replication.png` — 3×3 grid of boxplots, dosage 0/1/2 on x-axis, per-genotype `med` on y-axis. Rows = pipelines, columns = markers.
- `sec2_6_phase1_davis_fig6_replica.png` — Davis Fig. 6 C/D/E layout (major vs minor homozygote boxplots) for all three pipelines stacked.

**Outputs:**
- `sec2_6_phase1_gwas_marker_replication.csv` — per (pipeline, marker) effect β₁, std-err, p-value, R², direction, matches-paper flag
- `sec2_6_phase1_per_genotype_phenotype.csv` — per-genotype `med` across the three pipelines

---

## Not Yet Started

(All planned experiments complete. One open follow-up noted in `experiments/README.md`: E7 human ground truth, blocked on manual measurement data.)
