# =====================================================================
#  03_grn_inference.R
#  Gene regulatory network inference with GENIE3 (Methods 2.3)
#
#  - 1,000 randomised trees per TF, feature subsampling ratio 0.7
#  - post-hoc pruning of edges with |Pearson r| <= 0.3
#  - target universe: 13 M2 genes + co-expressed neighbours
#  Out: data/processed/grn_edges.rds        (data.frame: TF, target, weight)
#       data/processed/grn_adjacency.rds    (sparse, TFs x targets)
#       results/Table_S1_network_summary.csv (node degrees; -> Figure S1)
# =====================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat); library(GENIE3); library(Matrix); library(dplyr)
})

set.seed(SEED)
gse   <- GEO$discovery
seu   <- readRDS(file.path(DIRS$data_proc, paste0(gse, "_myeloid.rds")))
expr  <- GetAssayData(seu, slot = "data")

tfs       <- readLines(here("data", "processed", "candidate_tfs.txt"))
targets   <- readLines(here("data", "processed", "grn_target_genes.txt"))
targets   <- setdiff(unique(c(targets, M2_SIGNATURE)), tfs)   # TFs are regulators, not targets
genes_use <- unique(c(tfs, targets))
expr_use  <- as.matrix(t(expr[genes_use, , drop = FALSE]))    # cells x genes
storage.mode(expr_use) <- "double"

## ---- GENIE3 ---------------------------------------------------------
message("[03] Running GENIE3: ", length(tfs), " regulators x ",
        ncol(expr_use), " genes, ", GRN$n_trees, " trees")
weightMat <- GENIE3(expr_use,
                    regulators = tfs,
                    targets    = colnames(expr_use),
                    nTrees     = GRN$n_trees,
                    mtry       = max(1, floor(sqrt(ncol(expr_use)) * GRN$mtry_ratio)),
                    verbose    = TRUE)
## weightMat is a list (one matrix per target) -> long edge table
edges <- lapply(names(weightMat), function(tg) {
  w <- weightMat[[tg]]
  data.frame(regulator = rownames(w),
             target    = tg,
             weight    = as.numeric(w[, 1]),
             stringsAsFactors = FALSE)
})
edges <- bind_rows(edges)

## ---- |r| post-hoc pruning -------------------------------------------
r_mat <- cor(expr_use)
edges$abs_r <- mapply(function(a, b) abs(r_mat[a, b]), edges$regulator, edges$target)
edges <- edges[!is.na(edges$abs_r) & edges$abs_r > GRN$prune_abs_r, ]

## ---- adjacency ------------------------------------------------------
adj <- xtabs(weight ~ regulator + target, data = edges)   # sparse-friendly
adj <- as(adj, "sparseMatrix")

message("[03] Edges retained: ", nrow(edges),
        "  linking ", length(unique(edges$regulator)), " TFs to ",
        length(unique(edges$target)), " targets  (expected 9,847)")

saveRDS(edges, file.path(DIRS$data_proc, "grn_edges.rds"))
saveRDS(adj,   file.path(DIRS$data_proc, "grn_adjacency.rds"))

## ---- node degree summary (Figure S1) --------------------------------
deg <- data.frame(
  node      = union(rownames(adj), colnames(adj)),
  out_degree = as.numeric(rowSums(adj > 0))[union(rownames(adj), colnames(adj))],
  in_degree  = as.numeric(colSums(adj > 0))[union(rownames(adj), colnames(adj))]
)
deg$out_degree[is.na(deg$out_degree)] <- 0
deg$in_degree[is.na(deg$in_degree)]   <- 0
deg$class <- ifelse(deg$node %in% tfs, "TF", "target")
write.csv(deg, here("results", "Table_S1_network_summary.csv"), row.names = FALSE)

session_info_to_file(here("env", "sessionInfo_03.txt"))
message("[03] Done.")
