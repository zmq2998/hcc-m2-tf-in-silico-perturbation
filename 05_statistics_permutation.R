# =====================================================================
#  05_statistics_permutation.R
#  Permutation testing + effect sizes (Methods 2.4 / 2.7)
#
#  Null distribution: 1,000 permutations in which TF IDENTITY LABELS are
#  shuffled across the network edges.  This preserves network topology and
#  the target-gene degree distribution while destroying the TF-target
#  associations.
#
#  standardized effect size = (observed - mean_null) / sd_null
#  fold change               = observed / mean_null
#  empirical p               = (1 + #{null >= obs}) / (1 + n_perm)
#  q values                  = Benjamini-Hochberg
#
#  Out: results/Table_S1_significant_TFs.csv  (the 25 significant TFs)
# =====================================================================

source("config.R")
suppressPackageStartupMessages({ library(Matrix); library(dplyr) })

source("04_insilico_perturbation.R")   # provides rwr_stationary() and objects

set.seed(SEED)

## ---- observed scores -------------------------------------------------
observed <- impact$impact_raw
names(observed) <- impact$TF

## ---- permutation: shuffle TF labels across edges ---------------------
perm_matrix <- matrix(NA_real_, nrow = PERM$n_perm, ncol = length(tfs),
                      dimnames = list(NULL, tfs))

edges <- readRDS(file.path(DIRS$data_proc, "grn_edges.rds"))

for (b in seq_len(PERM$n_perm)) {
  ## shuffle the regulator labels (topology + target degrees preserved)
  set.seed(SEED + b)
  shuffled <- edges
  shuffled$regulator <- sample(shuffled$regulator)

  A_b <- Matrix::Matrix(0, nrow = length(nodes_all), ncol = length(nodes_all),
                        dimnames = list(nodes_all, nodes_all), sparse = TRUE)
  ab <- xtabs(weight ~ regulator + target, data = shuffled)
  ab <- as(ab, "sparseMatrix")
  A_b[rownames(ab), colnames(ab)] <- ab

  pi_b <- sapply(tfs, function(tf) {
    if (!tf %in% nodes_all) return(rep(NA_real_, length(nodes_all)))
    rwr_stationary(A_b, seed = tf)
  })
  rownames(pi_b) <- nodes_all
  perm_matrix[b, ] <- apply(pi_b[m2_present, , drop = FALSE], 2,
                            function(x) sum(x, na.rm = TRUE) / N_M2)
  if (b %% 100 == 0) message("[05] permutation ", b, "/", PERM$n_perm)
}

saveRDS(perm_matrix, file.path(DIRS$data_proc, "permutation_scores.rds"))

## ---- summarise -------------------------------------------------------
res <- data.frame(
  TF              = tfs,
  impact_score    = as.numeric(observed[tfs]),
  mean_null       = colMeans(perm_matrix, na.rm = TRUE),
  sd_null         = apply(perm_matrix, 2, sd, na.rm = TRUE)
)
res$standardized_effect_size <- (res$impact_score - res$mean_null) / res$sd_null
res$fold_change              <- res$impact_score / res$mean_null
res$p_empirical <- mapply(function(tf, obs) {
  n <- perm_matrix[, tf]
  (1 + sum(n >= obs, na.rm = TRUE)) / (1 + PERM$n_perm)
}, res$TF, res$impact_score)
res$q_value <- p.adjust(res$p_empirical, method = PERM$p_adjust)
res <- res %>% arrange(desc(standardized_effect_size))

sig <- subset(res, q_value < PERM$alpha)
message("[05] Significant TFs (BH q < ", PERM$alpha, "): ", nrow(sig),
        " (expected 25); top: ", paste(head(sig$TF, 3), collapse = ", "))

write.csv(res, here("results", "Table_S1_all_TFs_statistics.csv"), row.names = FALSE)
write.csv(sig, here("results", "Table_S1_significant_TFs.csv"),      row.names = FALSE)

session_info_to_file(here("env", "sessionInfo_05.txt"))
message("[05] Done.")
