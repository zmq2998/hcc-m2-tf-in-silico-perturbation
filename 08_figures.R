# =====================================================================
#  08_figures.R
#  Reproduce main-text Figures 1-4 and Supplementary Figures S1-S2
#  from the objects produced by scripts 01-07
#
#  Out: figures/Figure1_*.pdf ... figures/Figure_S1_*.pdf
# =====================================================================

source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(Matrix); library(pheatmap)
})

theme_set(theme_bw(base_size = 11))

## ================= Figure 1A: impact vs expression correlation ========
sig_all <- read.csv(here("results", "Table_S1_all_TFs_statistics.csv"))
seu     <- readRDS(file.path(DIRS$data_proc, paste0(GEO$discovery, "_myeloid.rds")))
expr    <- GetAssayData(seu, slot = "data")

m2_present <- intersect(M2_SIGNATURE, rownames(expr))
corr_score <- sapply(sig_all$TF, function(tf)
  mean(sapply(m2_present, function(g)
    suppressWarnings(cor(as.numeric(expr[tf, ]), as.numeric(expr[g, ]))))))

df1 <- data.frame(TF = sig_all$TF,
                  impact = sig_all$standardized_effect_size,
                  correlation = corr_score,
                  significant = sig_all$q_value < PERM$alpha)

p1a <- ggplot(df1, aes(correlation, impact, colour = significant)) +
  geom_point(alpha = .8) +
  geom_text(data = subset(df1, significant & rank(-impact) <= 10),
            aes(label = TF), vjust = -0.7, size = 3, show.legend = FALSE) +
  labs(x = "Expression-based correlation score",
       y = "RWR regulatory impact (standardized effect size)",
       colour = "q < 0.05") +
  theme(legend.position = "bottom")
ggsave(here("figures", "Figure1A_impact_vs_correlation.pdf"), p1a, width = 5, height = 4.2)

## ================= Figure 1B: top TF expression heatmap ===============
top_tfs <- head(read.csv(here("results", "Table_S1_significant_TFs.csv"))$TF, 10)
mat <- as.matrix(expr[top_tfs, ])
mat <- t(scale(t(mat)))
pheatmap(mat, cluster_rows = TRUE, show_colnames = FALSE,
         filename = here("figures", "Figure1B_topTF_expression_heatmap.pdf"),
         width = 6, height = 3.5)

## ================= Figure 2A/2B: perturbation effects =================
eff <- read.csv(here("results", "Table_S2_per_target_effects.csv"), row.names = 1)
prop <- data.frame(TF = rownames(eff),
                   prop = apply(eff, 1, function(x)
                     mean(abs(x - mean(x, na.rm = TRUE)) > RWR$effect_delta, na.rm = TRUE)))
p2a <- ggplot(prop, aes(reorder(TF, prop), prop)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Proportion of M2 targets influenced (|delta pi| > 0.05)")
ggsave(here("figures", "Figure2A_proportion_targets.pdf"), p2a, width = 5, height = 5)

per_target <- data.frame(TF = rownames(eff), mean_delta = rowMeans(eff, na.rm = TRUE))
p2b <- ggplot(per_target, aes(reorder(TF, mean_delta), mean_delta)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = median(per_target$mean_delta), linetype = "dashed") +
  labs(x = NULL, y = "Mean predicted effect per target (delta pi)")
ggsave(here("figures", "Figure2B_per_target_effect.pdf"), p2b, width = 5, height = 5)

## ================= Figure 2C: cross-regulatory network ================
## visualised with igraph / ggraph in 09_network_figure.R (optional)

## ================= Figure 3A/3B =======================================
pt_tab <- read.csv(here("results", "Table_S4_TF_pseudotime_dynamics.csv"))
pt_long <- tidyr::pivot_longer(pt_tab, -pseudotime, names_to = "TF", values_to = "expr")
p3a <- ggplot(pt_long, aes(pseudotime, expr, colour = TF)) +
  geom_line(linewidth = .8) +
  labs(x = "Pseudotime (monocyte -> macrophage)", y = "Smoothed expression") +
  theme(legend.position = "right")
ggsave(here("figures", "Figure3A_TF_pseudotime_curves.pdf"), p3a, width = 6, height = 4)

corr_tab <- as.matrix(read.csv(here("results", "Table_S5_TF_M2_correlations.csv"), row.names = 1))
pheatmap(corr_tab, cluster_rows = TRUE, cluster_cols = TRUE,
         filename = here("figures", "Figure3B_TF_M2_correlation_heatmap.pdf"),
         width = 6, height = 3.5)

## ================= Figure S1: degree distribution =====================
deg <- read.csv(here("results", "Table_S1_network_summary.csv"))
pS1 <- ggplot(deg, aes(out_degree, fill = class)) +
  geom_histogram(position = "identity", alpha = .7, bins = 30) +
  labs(x = "Out-degree", y = "Number of nodes", fill = "Node class")
ggsave(here("figures", "Figure_S1_degree_distribution.pdf"), pS1, width = 5, height = 3.5)

## Figure 4 (schematic model) is drawn in a vector editor from the
## hierarchy derived in Table_S1_significant_TFs.csv -> figures/Figure4_model.pdf

session_info_to_file(here("env", "sessionInfo_08.txt"))
message("[08] Figures written to ", DIRS$figures)
