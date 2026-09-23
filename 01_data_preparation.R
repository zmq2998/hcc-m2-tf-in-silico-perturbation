# =====================================================================
#  01_data_preparation.R
#  Single-cell RNA-seq processing + myeloid / M2 annotation
#  (Methods 2.1) -- Seurat v4 workflow
#
#  Input : 10x-format count matrices downloaded from GEO
#          GSE140228 (discovery), GSE149614 (validation)
#  Output: data/processed/<GSE>_myeloid.rds  (Seurat object, myeloid subset)
#          data/processed/<GSE>_m2_summary.csv
# =====================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}

## ---------------------------------------------------------------------
## Helper: QC + myeloid subsetting for one dataset
## ---------------------------------------------------------------------
process_dataset <- function(gse, sample_dir) {
  message("[01] Processing ", gse, " from ", sample_dir)

  counts <- Read10X(data.dir = sample_dir)
  seu <- CreateSeuratObject(counts = counts, project = gse,
                            min.cells = 3, min.features = 200)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  seu <- subset(seu,
                subset = nFeature_RNA >= QC$min_genes_per_cell &
                         nFeature_RNA <= QC$max_genes_per_cell &
                         percent.mt   <  QC$max_pct_mito)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.8, verbose = FALSE)

  ## ---- myeloid annotation (canonical markers) ------------------------
  FeaturePlot(seu, features = MYELOID_MARKERS)          # QC visual (optional)
  expr <- GetAssayData(seu, slot = "data")
  myeloid_score <- colSums(expr[intersect(MYELOID_MARKERS, rownames(expr)), , drop = FALSE])
  seu$myeloid_score <- myeloid_score
  seu$is_myeloid    <- myeloid_score > 0

  seu_mye <- subset(seu, subset = is_myeloid)
  seu_mye <- FindVariableFeatures(seu_mye, nfeatures = 2000, verbose = FALSE)
  seu_mye <- ScaleData(seu_mye, verbose = FALSE)
  seu_mye <- RunPCA(seu_mye, npcs = 30, verbose = FALSE)
  seu_mye <- FindNeighbors(seu_mye, dims = 1:30, verbose = FALSE)
  seu_mye <- FindClusters(seu_mye, resolution = QC$myeloid_resolution, verbose = FALSE)
  seu_mye <- RunUMAP(seu_mye, dims = 1:30, verbose = FALSE)

  ## ---- M2 score (13-gene signature) ----------------------------------
  present  <- intersect(M2_SIGNATURE, rownames(seu_mye))
  seu_mye  <- AddModuleScore(seu_mye, features = list(present), name = "M2_score")

  m2_summary <- seu_mye@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(n_cells      = dplyr::n(),
              mean_M2_score = mean(M2_score1, na.rm = TRUE),
              mean_CD163    = mean(expm1(GetAssayData(seu_mye, slot = "data")["CD163", ])),
              .groups = "drop")
  write.csv(m2_summary,
            file.path(DIRS$data_proc, paste0(gse, "_m2_summary.csv")),
            row.names = FALSE)
  saveRDS(seu_mye, file.path(DIRS$data_proc, paste0(gse, "_myeloid.rds")))
  seu_mye
}

## ---------------------------------------------------------------------
invisible(lapply(names(GEO), function(k) {
  gse        <- GEO[[k]]
  sample_dir <- get_arg(paste0(k, "_dir"), file.path(DIRS$data_raw, gse))
  if (!dir.exists(sample_dir)) {
    warning("Directory not found for ", gse, ": ", sample_dir, " -- skipping. ",
            "Download it from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse)
    return(invisible(NULL))
  }
  process_dataset(gse, sample_dir)
}))

session_info_to_file(here("env", "sessionInfo_01.txt"))
message("[01] Done.")
