scripts <- c(
  "01_load_data_and_PCA.R",
  "02_post_differentiation_limma.R",
  "03_systematic_lipid_annotation.R",
  "04_global_heatmap.R",
  "05_selected_PCA_and_three_group_PLSDA.R"
)

if (!all(file.exists(scripts))) stop("Open C2C12_Lipidomics_analysis.Rproj")
for (script in scripts) source(script, local = .GlobalEnv)

cat("Lipidomics analysis completed successfully.\n")
