# =====================================================================
#  04_insilico_perturbation.R
#  In silico TF perturbation by random walk with restart (Methods 2.4)
#
#  For EVERY candidate TF (n = 63) the TF is used as an individual seed
#  node; the walker propagates over the directed GRN and the stationary
#  probability mass landing on the 13 M2 signature genes defines the
#  regulatory impact score:
#
#      score(TF) = (1 / |M2|) * sum_{g in M2} pi_g(TF)
#
#  Out: results/Table_S1_TF_impact_scores.csv  (63 rows, raw scores)
#       data/processed/rwr_stationary.rds      (matrix pi: genes x TFs)
# =====================================================================

source("config.R")
suppressPackageStartupMessages({ library(Matrix); library(dplyr) })

## ---------------------------------------------------------------------
## RWR core: stationary distribution for one seed node
##   pi = r * P^T pi + (1 - r) * e_seed     (column-stochastic P)
## ---------------------------------------------------------------------
rwr_stationary <- function(adj, seed,
                           restart = RWR$restart,
                           tol = RWR$tol, max_iter = RWR$max_iter) {
  ## column-stochastic transition matrix (out-degree normalised)
  out_deg <- Matrix::colSums(adj)
  out_deg[out_deg == 0] <- 1
  P <- adj %*% Matrix::Diagonal(x = 1 / out_deg)   # genes x genes (column-stochastic)

  n <- nrow(P)
  e <- rep(0, n); names(e) <- rownames(P); e[seed] <- 1
  pi <- e
  for (i in seq_len(max_iter)) {
    pi_new <- restart * as.numeric(P %*% pi) + (1 - restart) * e
    if (sum(abs(pi_new - pi)) < tol) { pi <- pi_new; break }
    pi <- pi_new
  }
  pi / sum(pi)
}

## ---------------------------------------------------------------------
adj <- readRDS(file.path(DIRS$data_proc, "grn_adjacency.rds"))
tfs <- readLines(here("data", "processed", "candidate_tfs.txt"))
adj <- as(adj, "dgCMatrix")

## square the network so that TFs can be reached as walkers (Methods: the
## regulator is the seed, propagation runs over the inferred directed net)
nodes_all <- union(rownames(adj), colnames(adj))
A <- Matrix::Matrix(0, nrow = length(nodes_all), ncol = length(nodes_all),
                    dimnames = list(nodes_all, nodes_all), sparse = TRUE)
A[rownames(adj), colnames(adj)] <- adj

m2_present <- intersect(M2_SIGNATURE, nodes_all)
stopifnot(length(m2_present) > 0)

pi_list <- sapply(tfs, function(tf) {
  if (!tf %in% nodes_all) return(rep(NA, length(nodes_all)))
  rwr_stationary(A, seed = tf)
})
rownames(pi_list) <- nodes_all
saveRDS(pi_list, file.path(DIRS$data_proc, "rwr_stationary.rds"))

impact <- data.frame(
  TF          = tfs,
  impact_raw  = apply(pi_list[m2_present, , drop = FALSE], 2, function(x)
                      sum(x, na.rm = TRUE) / N_M2),
  row.names   = NULL
) %>% arrange(desc(impact_raw))

## predicted per-target effects (used by Figure 2A/2B)
target_effect <- t(pi_list[m2_present, , drop = FALSE])
target_effect <- as.data.frame(target_effect)
write.csv(data.frame(TF = rownames(target_effect), target_effect, row.names = NULL),
          here("results", "Table_S2_per_target_effects.csv"), row.names = FALSE)

## proportion of M2 targets with |delta pi| > 0.05
prop_influenced <- apply(target_effect, 1, function(x)
  mean(abs(x - mean(x, na.rm = TRUE)) > RWR$effect_delta, na.rm = TRUE))
impact$prop_targets_influenced <- prop_influenced[impact$TF]

write.csv(impact, here("results", "Table_S1_TF_impact_scores_raw.csv"), row.names = FALSE)
message("[04] Impact scores computed for ", nrow(impact), " TFs. Top 5: ",
        paste(head(impact$TF, 5), collapse = ", "))

session_info_to_file(here("env", "sessionInfo_04.txt"))
message("[04] Done.")
