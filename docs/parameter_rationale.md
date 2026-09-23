# Parameter rationale and analysis notes

This note records **why** each parameter in `config.R` has the value it has,
so that reviewers can judge analyst degrees of freedom. Values marked
【verify】 must be confirmed against the final manuscript text / original runs.

## 2.1 Single-cell processing (Seurat 4.3.0)

| Parameter | Value | Rationale |
|---|---|---|
| cells retained | 200–5,000 detected genes | standard scRNA-seq QC window: removes empty droplets and likely doublets |
| mitochondrial cut-off | < 20% | removes stressed/dying cells |
| myeloid resolution | 0.8 | produces stable myeloid sub-clusters (M1-like, M2-like, monocyte-like) without over-fragmenting |
| M2 definition | C1QA, C1QB, C1QC, CD163, MERTK, IL10 | canonical M2/TAM markers; consistent across both datasets |
| discovery cells QC-passed | 14,327 | reported in the manuscript |

## 2.2 Gene sets

* **63 candidate TFs**: Catalogue of Human Transcription Factors, filtered to
  detectable expression in ≥ 5% of myeloid cells with mean normalised
  expression > 0.1.
* **13-gene M2 signature**: C1QA, C1QB, C1QC, CD163, MERTK, IL10, TGFB1, FN1,
  COL6A3, DPP4, FCGR2B, IL1R2, IL1RN — cross-checked against TRRUST v2.

## 2.3 GRN inference (GENIE3 1.20.0)

| Parameter | Value | Rationale |
|---|---|---|
| trees per TF | 1,000 | accuracy/runtime trade-off; GENIE3 is robust to dropout |
| feature subsampling | 0.7 × √(#genes) | default GENIE3 `mtry` scaled up for high-dimensional scRNA-seq |
| edge pruning | \|Pearson *r*\| > 0.3 | removes weakly supported TF–target pairs returned by the ensemble |
| target universe | 13 M2 genes + co-expressed neighbours (\|*r*\| > 0.3) | keeps the network interpretable around the programme of interest |
| resulting network | 9,847 directed edges, 63 TFs | reported in the manuscript and Figure S1 |

## 2.4 In silico perturbation (RWR)

* **Seed = every candidate TF individually** (n = 63). This differs from
  correlation-based prioritisation, which cannot express functional impact.
* **Impact score** `score(TF) = (1/|M2|) · Σ_{g∈M2} π_g(TF)`, normalised by
  |M2| = 13 so that hub TFs with more targets are not inflated.
* **Restart probability r = 0.7**, selected by grid search over 0.1–0.9 as the
  value maximising stability of the top-25 ranking (Spearman ρ > 0.91 for
  r = 0.65–0.75; Figure S2, `06_robustness_restart_probability.R`).
* **Convergence**: ‖Δπ‖ < 1 × 10⁻⁶, up to 1,000 iterations.
* **"Target influenced"** = predicted change in stationary probability > 0.05
  (≈ 90th percentile of the null distribution).

## 2.4 / 2.7 Statistics

* **Null model**: 1,000 permutations with **TF identity labels shuffled across
  the network edges** → preserves network topology and target-gene degree
  distribution while destroying TF–target associations.
* **Standardized effect size** = (observed − mean_null) / sd_null. (Earlier
  drafts called this a "z-score"; renamed to avoid implying a normal
  approximation.)
* **Fold change** = observed / mean_null.
* **Empirical p** = (1 + #{null ≥ observed}) / (1 + n_perm).
* **Multiple testing**: Benjamini–Hochberg *q* < 0.05 → 25 significant TFs.

## 2.5 Pseudotime

* monocle3 1.3.1 on the myeloid subset; root = monocyte-like pole (maximal
  FCN1/S100A8/S100A9/VCAN/LYZ score).
* TF dynamics smoothed with cubic splines, 3 degrees of freedom.
* Limitation acknowledged in the manuscript: pseudotime assumes a
  unidirectional monocyte → M2 trajectory and cannot capture TAM plasticity.

## Removed analysis

The earlier 5-fold cross-validation (MAE) was **dropped**: it is not a
standard validation for RWR-based prioritisation and its fold definition was
ambiguous. It is replaced by the restart-probability stability analysis
(Figure S2). If a link-prediction-style validation is required later,
implement it as hidden-edge hold-out with AUC.

## 【verify】 before submission

1. Confirm the final numbers in `results/` reproduce the manuscript values
   (25 significant TFs; NR4A2 = 11.77, ATF3 = 9.87, JUN = 6.65; 9,847 edges).
2. Confirm the patient/sample description of GSE140228 and GSE149614 against
   the GEO records (Section 2.1 of the manuscript).
3. If the original runs used different parameter values, update `config.R`
   rather than editing scripts — everything is centralised there.
