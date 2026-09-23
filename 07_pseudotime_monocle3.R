# =====================================================================
#  07_pseudotime_monocle3.R
#  Monocyte -> macrophage pseudotime trajectory + TF dynamics
#  (Methods 2.5, Figure 3A)
#
#  Out: data/processed/pseudotime.rds
#       results/Table_S4_TF_pseudotime_dynamics.csv
# =====================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(dplyr)
})

set.seed(SEED)
gse <- GEO$discovery
seu <- readRDS(file.path(DIRS$data_proc, paste0(gse, "_myeloid.rds")))

## ---- Seurat -> CellDataSet ------------------------------------------
cds <- as.cell_data_set(seu)
cds <- cluster_cells(cds, resolution = QC$myeloid_resolution)
cds <- learn_graph(cds, use_partition = TRUE)

## ---- root = monocyte-like pole (highest classical monocyte score) ---
mono_markers <- intersect(c("FCN1", "S100A8", "S100A9", "VCAN", "LYZ"),
                          rownames(seu))
mono_score <- Matrix::colMeans(GetAssayData(seu, slot = "data")[mono_markers, , drop = FALSE])
cds <- order_cells(cds, root_cells = names(which.max(mono_score)))
saveRDS(cds, file.path(DIRS$data_proc, "pseudotime.rds"))

## ---- TF dynamics along pseudotime (cubic spline, 3 df) ---------------
pt <- pseudotime(cds)
tfs <- readLines(here("data", "processed", "candidate_tfs.txt"))
top_sig <- read.csv(here("results", "Table_S1_significant_TFs.csv"))$TF
use_tfs <- head(top_sig, 10)

expr_tf <- as.matrix(t(GetAssayData(seu, slot = "data")[use_tfs, , drop = FALSE]))
spline_curves <- apply(expr_tf, 2, function(y) {
  ok <- is.finite(pt) & is.finite(y)
  predict(smooth.spline(pt[ok], y[ok], df = 3), x = seq(min(pt[ok]), max(pt[ok]), length.out = 100))$y
})
out <- data.frame(pseudotime = seq(0, 1, length.out = 100), spline_curves)
write.csv(out, here("results", "Table_S4_TF_pseudotime_dynamics.csv"), row.names = FALSE)

## ---- TF vs M2 marker correlations (Figure 3B) ------------------------
m2_present <- intersect(M2_SIGNATURE, rownames(seu))
corr <- outer(use_tfs, m2_present, Vectorize(function(a, b)
  suppressWarnings(cor(as.numeric(GetAssayData(seu, slot = "data")[a, ]),
                       as.numeric(GetAssayData(seu, slot = "data")[b, ])))))
dimnames(corr) <- list(use_tfs, m2_present)
write.csv(corr, here("results", "Table_S5_TF_M2_correlations.csv"))

session_info_to_file(here("env", "sessionInfo_07.txt"))
message("[07] Done.")
