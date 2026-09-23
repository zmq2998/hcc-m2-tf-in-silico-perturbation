# =====================================================================
#  config.R  —  Central configuration for the analysis pipeline
#  Project: TFs governing M2 polarization in HCC (in silico perturbation)
#  Repo   : hcc-m2-tf-in-silico-perturbation
#
#  All paths are RELATIVE to the repository root (see `here()` below).
#  Do NOT hard-code absolute paths -- override with environment variables
#  when needed:
#      Rscript 01_data_preparation.R --input_dir=/path/to/data
# =====================================================================

## ---- repository root ------------------------------------------------
here <- function(...) {
  root <- tryCatch({
    if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
  }, error = function(e) getwd())
  file.path(root, ...)
}

DIRS <- list(
  data_raw    = here("data", "raw"),        # place GEO downloads here (git-ignored)
  data_proc   = here("data", "processed"),  # intermediate .rds objects
  results     = here("results"),            # tables
  figures     = here("figures"),            # output figures
  logs        = here("logs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

## ---- random seed (reproducibility) ----------------------------------
SEED <- 20260923L
set.seed(SEED)

## ---- datasets (public, no raw data redistributed) -------------------
GEO <- list(
  discovery = "GSE140228",   # Zhang Q. et al. Cell 2019, 179:829-845  (CD45+ immune cells, 5 sites)
  validation = "GSE149614"   # Lu Y. et al. Nat Commun 2022, 13:4077   (10 HCC patients, 4 sites)
)

## ---- QC parameters (Methods 2.1) ------------------------------------
QC <- list(
  min_genes_per_cell = 200,
  max_genes_per_cell = 5000,
  max_pct_mito       = 20,
  myeloid_resolution = 0.8
)
MYELOID_MARKERS <- c("CD14", "CD68", "CD163")

## ---- gene sets (Methods 2.2) ----------------------------------------
M2_SIGNATURE <- c(
  "C1QA", "C1QB", "C1QC", "CD163", "MERTK", "IL10", "TGFB1",
  "FN1", "COL6A3", "DPP4", "FCGR2B", "IL1R2", "IL1RN"
)                                          # n = 13, |M2| below
N_M2 <- length(M2_SIGNATURE)

TF_MIN_PCT_CELLS <- 0.05                   # TF detected in >= 5% of myeloid cells
TF_MIN_MEAN_EXPR <- 0.1

## ---- GRN inference (Methods 2.3) ------------------------------------
GRN <- list(
  n_trees        = 1000,
  mtry_ratio     = 0.7,
  prune_abs_r    = 0.3                     # post-hoc |r| pruning threshold
)

## ---- RWR in silico perturbation (Methods 2.4) -----------------------
RWR <- list(
  restart      = 0.7,                      # r, selected by grid search
  r_grid       = seq(0.1, 0.9, by = 0.05),
  tol          = 1e-6,                     # convergence: ||delta|| < tol
  max_iter     = 1000,
  effect_delta = 0.05                      # threshold for "target influenced"
)

## ---- permutation testing (Methods 2.4 / 2.7) ------------------------
PERM <- list(
  n_perm = 1000,
  alpha  = 0.05,
  p_adjust = "BH"
)

## ---- scoring function ----------------------------------------------
## score(TF) = (1 / |M2|) * sum_{g in M2} pi_g(TF)
score_tf <- function(pi_vector, m2_genes, n_m2 = length(m2_genes)) {
  sum(pi_vector[m2_genes], na.rm = TRUE) / n_m2
}

## ---- reproducibility ------------------------------------------------
session_info_to_file <- function(file = here("env", "sessionInfo.txt")) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), con = file)
  invisible(file)
}
