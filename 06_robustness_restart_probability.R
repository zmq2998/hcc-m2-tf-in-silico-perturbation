# =====================================================================
#  06_robustness_restart_probability.R
#  Stability of the top-25 ranking across restart probabilities r
#  (Methods 2.4, Figure S2)  -- complements the 5-fold network
#  cross-validation (Methods 2.8); it does not replace it
#
#  Out: results/Table_S3_restart_stability.csv
#       figures/Figure_S2_restart_stability.pdf
# =====================================================================

source("config.R")
suppressPackageStartupMessages({ library(Matrix); library(dplyr); library(ggplot2) })

adj     <- as(readRDS(file.path(DIRS$data_proc, "grn_adjacency.rds")), "dgCMatrix")
tfs     <- readLines(here("data", "processed", "candidate_tfs.txt"))
pi_ref  <- readRDS(file.path(DIRS$data_proc, "rwr_stationary.rds"))
nodes_all <- rownames(pi_ref)
m2_present <- intersect(M2_SIGNATURE, nodes_all)

A <- Matrix::Matrix(0, nrow = length(nodes_all), ncol = length(nodes_all),
                    dimnames = list(nodes_all, nodes_all), sparse = TRUE)
A[rownames(adj), colnames(adj)] <- adj

reference <- apply(pi_ref[m2_present, , drop = FALSE], 2,
                   function(x) sum(x, na.rm = TRUE) / N_M2)
ref_rank  <- rank(-reference)

rows <- lapply(RWR$r_grid, function(r) {
  pi_r <- sapply(tfs, function(tf) rwr_stationary(A, tf, restart = r))
  rownames(pi_r) <- nodes_all
  sc <- apply(pi_r[m2_present, , drop = FALSE], 2,
              function(x) sum(x, na.rm = TRUE) / N_M2)
  data.frame(restart = r,
             spearman_rho = suppressWarnings(cor(sc, reference, method = "spearman")),
             top25_overlap = length(intersect(names(sort(-sc))[1:25],
                                              names(sort(-reference))[1:25])) / 25,
             row.names = NULL)
})
stab <- bind_rows(rows)
write.csv(stab, here("results", "Table_S3_restart_stability.csv"), row.names = FALSE)

p <- ggplot(stab, aes(restart, spearman_rho)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.91, linetype = "dashed") +
  labs(x = "Restart probability r", y = "Spearman rho vs. reference ranking",
       title = "Stability of the top-25 TF ranking across r") +
  theme_bw()
ggsave(here("figures", "Figure_S2_restart_stability.pdf"), p, width = 5, height = 3.5)

message("[06] r = 0.7 selected: rho = ",
        round(stab$spearman_rho[which.min(abs(stab$restart - RWR$restart))], 3))
session_info_to_file(here("env", "sessionInfo_06.txt"))
message("[06] Done.")
