# =====================================================================
#  02_gene_sets.R
#  Candidate TF list + M2 signature definition (Methods 2.2)
#
#  TFs  : Catalogue of Human Transcription Factors (Lambert et al. 2018)
#         filtered to those detected in >= 5% of myeloid cells
#  M2   : 13-gene signature (C1QA, C1QB, C1QC, CD163, MERTK, IL10,
#         TGFB1, FN1, COL6A3, DPP4, FCGR2B, IL1R2, IL1RN)
#  Out  : data/processed/candidate_tfs.txt   (63 TFs, one per line)
# =====================================================================

source("config.R")
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

gse        <- GEO$discovery
seu        <- readRDS(file.path(DIRS$data_proc, paste0(gse, "_myeloid.rds")))
expr       <- GetAssayData(seu, slot = "data")

## ---- TF universe -----------------------------------------------------
## Download the human TF list once (not redistributed here):
##   https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt
tf_file <- file.path(DIRS$data_raw, "TF_names_v_1.01.txt")
if (!file.exists(tf_file)) {
  stop("Please download the human TF list to ", tf_file,
       " (https://humantfs.ccbr.utoronto.ca/).")
}
tf_all <- readLines(tf_file)
tf_all <- unique(tf_all[nzchar(tf_all)])

## ---- detection filter (>= 5% of cells, mean normalized expr > 0.1) ---
pct_detected <- Matrix::rowMeans(expr > 0)
mean_expr    <- Matrix::rowMeans(expr)
keep <- intersect(tf_all, rownames(expr))
keep <- keep[pct_detected[keep] >= TF_MIN_PCT_CELLS &
             mean_expr[keep]    >  TF_MIN_MEAN_EXPR]

candidate_tfs <- sort(keep)
message("[02] Candidate TFs detected: ", length(candidate_tfs), " (expected 63)")

writeLines(candidate_tfs, here("data", "processed", "candidate_tfs.txt"))

## ---- M2 signature sanity check --------------------------------------
missing_m2 <- setdiff(M2_SIGNATURE, rownames(expr))
if (length(missing_m2))
  warning("M2 signature genes absent from the matrix: ",
          paste(missing_m2, collapse = ", "))

## ---- co-expressed neighbours of the M2 programme (Methods 2.3) ------
## genes with |r| > 0.3 against >= 1 M2 signature gene, computed on
## pseudo-bulk-free single-cell expression (Spearman for robustness)
m2_present <- intersect(M2_SIGNATURE, rownames(expr))
cors <- suppressWarnings(
  apply(expr, 1, function(x) max(abs(cor(as.numeric(x),
                                        as.numeric(Matrix::colMeans(expr[m2_present, , drop = FALSE])))))))
neighbours <- names(which(cors > GRN$prune_abs_r))
writeLines(neighbours, here("data", "processed", "grn_target_genes.txt"))
message("[02] GRN target gene universe: ", length(neighbours),
        " (13 M2 genes + co-expressed neighbours)")

session_info_to_file(here("env", "sessionInfo_02.txt"))
message("[02] Done.")
