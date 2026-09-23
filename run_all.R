# =====================================================================
#  run_all.R  —  execute the complete pipeline in order
#      Rscript run_all.R
# =====================================================================

steps <- c(
  "01_data_preparation.R",
  "02_gene_sets.R",
  "03_grn_inference.R",
  "04_insilico_perturbation.R",
  "05_statistics_permutation.R",
  "06_robustness_restart_probability.R",
  "07_pseudotime_monocle3.R",
  "08_figures.R"
)

t0 <- Sys.time()
for (s in steps) {
  message("\n===== ", s, " =====")
  source(s, echo = FALSE)
}
message("\nAll steps finished in ",
        round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

## full environment record
dir.create(here("env"), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), here("env", "sessionInfo.txt"))
