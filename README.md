# hcc-m2-tf-in-silico-perturbation

**In silico perturbation of single-cell gene regulatory networks to identify transcription factors governing M2 macrophage polarization in hepatocellular carcinoma.**

Analysis code accompanying the manuscript:

> *Identifying Transcription Factors Governing M2 Polarization in Hepatocellular Carcinoma via In Silico Perturbation of Single-Cell Gene Regulatory Networks* (submitted to **Life**, MDPI)

---

## Overview

Most computational TF-prioritisation strategies rank regulators by **static expression correlation** rather than by their **predicted functional impact** on a gene programme. This repository implements a framework that:

1. processes two public single-cell RNA-seq atlases of human HCC (`GSE140228`, `GSE149614`) with **Seurat**;
2. defines 63 candidate TFs and a 13-gene M2 signature;
3. infers a directed gene regulatory network (GRN) with **GENIE3** (1,000 trees per TF, feature subsampling 0.7) and post-hoc pruning at |Pearson *r*| > 0.3;
4. uses **every TF individually as a seed node** in a random walk with restart (RWR, restart probability *r* = 0.7) and scores its regulatory impact on the M2 programme as

   ```
   score(TF) = (1 / |M2|) * Σ_{g ∈ M2} π_g(TF)
   ```

   where `π_g(TF)` is the stationary probability of gene *g* when *TF* is the seed and |M2| = 13;
5. calibrates the score against **1,000 label permutations** (TF identity shuffled across edges, preserving topology and target degree distribution), reporting standardized effect sizes `(observed − mean_null) / sd_null`, fold-changes and Benjamini–Hochberg *q* values;
6. checks robustness of the top-25 ranking across restart probabilities *r* = 0.1–0.9 (Spearman ρ);
7. reconstructs a monocyte → macrophage **monocle3** pseudotime trajectory and TF–M2 correlation structure;
8. evaluates the inferred network by **5-fold cross-validation** under a
   hidden-edge link-prediction protocol, reporting the AUC (Methods 2.8).

Result: **25 TFs** significantly influence the M2 programme, led by **NR4A2** (standardized effect size = 11.77), **ATF3** (9.87) and **JUN** (6.65), with **PPARG** placed downstream of 14 of them.

---

## Repository structure

```
hcc-m2-tf-in-silico-perturbation/
├── config.R                          # all parameters, gene sets, paths (single source of truth)
├── run_all.R                         # runs 01 → 08 in order
├── 01_data_preparation.R             # QC, myeloid subsetting, M2 scoring      (Methods 2.1)
├── 02_gene_sets.R                    # 63 candidate TFs + 13-gene M2 signature (Methods 2.2)
├── 03_grn_inference.R                # GENIE3 + |r| > 0.3 pruning               (Methods 2.3)
├── 04_insilico_perturbation.R        # RWR per TF, score(TF)                    (Methods 2.4)
├── 05_statistics_permutation.R       # 1,000 permutations, effect sizes, BH q   (Methods 2.4/2.7)
├── 06_robustness_restart_probability.R  # ranking stability across r (Figure S2)
├── 07_pseudotime_monocle3.R          # pseudotime + TF dynamics                 (Methods 2.5)
├── 08_figures.R                      # Figures 1–4, S1, S2
├── data/
│   ├── raw/                          # GEO downloads (NOT tracked; see data/README.md)
│   └── processed/                    # intermediate .rds objects (NOT tracked)
├── results/                          # Table S1 … S5 (.csv, tracked)
├── figures/                          # generated figures (.pdf, tracked)
├── env/                              # sessionInfo() records
└── docs/                             # analysis notes / parameter rationale
```

---

## Requirements

* **R ≥ 4.4** (developed on R 4.5.3)
* CRAN / Bioconductor packages:

```r
install.packages(c("Seurat", "GENIE3", "Matrix", "dplyr", "ggplot2",
                   "pheatmap", "tidyr", "here"))
BiocManager::install(c("monocle3"))
```

Exact versions used are recorded in [`env/sessionInfo.txt`](env/sessionInfo.txt):
Seurat 4.3.0 · GENIE3 1.20.0 · monocle3 1.3.1 · pySCENIC 0.12.1.

---

## Data

Raw single-cell data are **not redistributed** here; both datasets are public.

| Dataset | Role | Source | Accession |
|---|---|---|---|
| `GSE140228` | discovery (CD45+ immune cells, 5 sites) | Zhang Q. *et al.*, *Cell* 2019, 179:829–845 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228) |
| `GSE149614` | validation (10 HCC patients, 4 sites) | Lu Y. *et al.*, *Nat. Commun.* 2022, 13:4077 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614) |

Additionally, the human TF list (`TF_names_v_1.01.txt`) is downloaded from
<https://humantfs.ccbr.utoronto.ca/> (Lambert *et al.*, *Cell* 2018).

Place the unpacked 10x matrices in `data/raw/GSE140228/` and `data/raw/GSE149614/`, then:

```r
Rscript run_all.R
```

or point to another location:

```r
Rscript 01_data_preparation.R --discovery_dir=/path/to/GSE140228 --validation_dir=/path/to/GSE149614
```

## Reproducibility notes

* Fixed seed `SEED <- 20260923L` (`config.R`); permutation block seeds are `SEED + b`.
* All parameters (QC thresholds, gene sets, tree count, pruning threshold, restart probability, permutation count) live in `config.R` — no magic numbers in the analysis scripts.
* Convergence criterion for RWR: `‖Δπ‖ < 1 × 10⁻⁶`, `max_iter = 1000`.
* Network cross-validation (Methods 2.8): 5 folds, 20% of the TF–target edges
  withheld per fold, GENIE3 re-run on the remaining edges, withheld edges scored
  against an equal number of random non-edges; performance reported as AUC.
* Runtimes: GENIE3 is the bottleneck (≈ hours on 8 cores for ~14k myeloid cells); set `GRN$n_trees` lower for a quick smoke test.

---

## Citation

If you use this code, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)).

## License

Code released under the [MIT License](LICENSE). Result tables and figures are released under CC-BY-4.0.

## Contact

Corresponding author: Mengqing Zhou (School of Biological Engineering, Henan University of Technology, Zhengzhou 450001, China) — tiankong168@haut.edu.cn; ORCID: 0009-0008-8091-9266.
