# Main analyses plus the original passage-model display requested for review.
# Set TRUE to reproduce the explicitly exploratory PCA and ASCA views.
run_exploratory <- FALSE
scripts <- c(
  "01_load_data_and_PCA.R",
  "02_post_differentiation_limma.R",
  "03_systematic_lipid_annotation.R",
  "04_global_heatmap.R",
  "05_selected_PCA_and_day_PLSDA.R",
  "06_original_passage_PLSDA.R"
)
if (run_exploratory) scripts <- c(scripts,
  "07_exploratory_PCA_axes.R", "08_exploratory_ASCA_scores.R")
script_directory <- getwd()
if (!all(file.exists(scripts))) stop("Run from RStudio or Code availability")
for (script in scripts) {
  setwd(script_directory)
  source(script, local = .GlobalEnv)
}
cat("Lipidomics analysis completed successfully.\n")
