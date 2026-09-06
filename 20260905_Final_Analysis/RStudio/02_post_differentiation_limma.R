# Post-differentiation passage contrasts in log2-normalized abundance.
# D0 is excluded; thresholds apply to the equal-weight average contrast.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("data.table", "limma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# Load Code 1 objects when needed.
code1_objects <- c(
  "qi_data", "analysis_log2", "analysis_feature_ids",
  "sample_metadata", "feature_qc_table", "analysis_root"
)
if (!all(vapply(code1_objects, exists, logical(1), inherits = TRUE))) {
  code1_candidates <- c(
    "01_load_data_and_PCA.R",
    file.path("RStudio", "01_load_data_and_PCA.R")
  )
  code1_file <- code1_candidates[file.exists(code1_candidates)][1]
  if (is.na(code1_file)) stop("Run Code 1 first or open this script from the RStudio project")
  source(code1_file)
}

table_directory <- file.path(analysis_root, "Tables")
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

model_metadata <- copy(sample_metadata)
model_metadata <- model_metadata[match(colnames(analysis_log2), biological_sample)]
if (!identical(model_metadata$biological_sample, colnames(analysis_log2))) {
  stop("Sample metadata and abundance columns are not aligned")
}

design_matrix <- model.matrix(~0 + condition, data = model_metadata)
colnames(design_matrix) <- levels(model_metadata$condition)

contrast_matrix <- makeContrasts(
  Middle_vs_Early_D1 = Middle_D1 - Early_D1,
  Late_vs_Early_D1 = Late_D1 - Early_D1,
  Middle_vs_Early_D3 = Middle_D3 - Early_D3,
  Late_vs_Early_D3 = Late_D3 - Early_D3,
  Middle_vs_Early_D6 = Middle_D6 - Early_D6,
  Late_vs_Early_D6 = Late_D6 - Early_D6,
  Middle_vs_Early_D9 = Middle_D9 - Early_D9,
  Late_vs_Early_D9 = Late_D9 - Early_D9,
  levels = design_matrix
)

base_fit <- lmFit(analysis_log2, design_matrix)
eight_contrast_fit <- eBayes(
  contrasts.fit(base_fit, contrast_matrix),
  trend = TRUE,
  robust = TRUE
)

# Joint test of the eight post-differentiation contrasts.
joint_results <- as.data.table(
  topTable(
    eight_contrast_fit,
    coef = seq_len(ncol(contrast_matrix)),
    number = Inf,
    sort.by = "none"
  ),
  keep.rownames = "met_id"
)
setnames(
  joint_results,
  c("AveExpr", "F", "P.Value", "adj.P.Val"),
  c("average_log2_abundance", "joint_F", "joint_p_value", "joint_fdr")
)
joint_results[, (colnames(contrast_matrix)) := NULL]

# Equal-weight average of the eight post-differentiation contrasts.
average_contrast <- matrix(
  rowMeans(contrast_matrix),
  ncol = 1,
  dimnames = list(rownames(contrast_matrix), "average_post_effect")
)
average_fit <- eBayes(
  contrasts.fit(base_fit, average_contrast),
  trend = TRUE,
  robust = TRUE
)
average_results <- as.data.table(
  topTable(average_fit, coef = 1, number = Inf, sort.by = "none"),
  keep.rownames = "met_id"
)
average_results <- average_results[, .(
  met_id,
  average_post_log2FC = logFC,
  average_post_t = t,
  average_post_p_value = P.Value,
  average_post_fdr = adj.P.Val
)]

individual_effects <- as.data.table(
  eight_contrast_fit$coefficients,
  keep.rownames = "met_id"
)
condition_means <- sapply(levels(model_metadata$condition), function(condition_name) {
  columns <- model_metadata[condition == condition_name, biological_sample]
  rowMeans(analysis_log2[, columns, drop = FALSE])
})
condition_mean_table <- as.data.table(condition_means, keep.rownames = "met_id")

tested_results <- Reduce(
  function(x, y) merge(x, y, by = "met_id", all = TRUE, sort = FALSE),
  list(joint_results, average_results, individual_effects, condition_mean_table)
)

effect_columns <- colnames(contrast_matrix)
middle_effect_columns <- paste0("Middle_vs_Early_D", c(1, 3, 6, 9))
late_effect_columns <- paste0("Late_vs_Early_D", c(1, 3, 6, 9))
tested_results[, mean_Middle_vs_Early_log2FC := rowMeans(.SD), .SDcols = middle_effect_columns]
tested_results[, mean_Late_vs_Early_log2FC := rowMeans(.SD), .SDcols = late_effect_columns]
tested_results[, average_fold_change_ratio := 2^average_post_log2FC]
tested_results[, average_1.5_fold_pass :=
  average_fold_change_ratio >= 1.5 |
    average_fold_change_ratio <= 1 / 1.5]

tested_results[, positive_contrasts := rowSums(.SD > 0), .SDcols = effect_columns]
tested_results[, negative_contrasts := rowSums(.SD < 0), .SDcols = effect_columns]
tested_results[, trajectory_direction := fcase(
  positive_contrasts == length(effect_columns), "higher",
  negative_contrasts == length(effect_columns), "lower",
  default = "mixed"
)]
tested_results[, all_eight_same_direction := trajectory_direction != "mixed"]
tested_results[, joint_selected := joint_fdr < 0.05]
tested_results[, post_selected :=
  average_post_fdr < 0.05 &
    average_1.5_fold_pass &
    all_eight_same_direction]
tested_results[, absolute_average_post_log2FC := abs(average_post_log2FC)]
tested_results[, post_rank := frank(average_post_fdr, ties.method = "min", na.last = "keep")]

feature_results <- merge(
  feature_qc_table,
  tested_results,
  by = "met_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(
  feature_results,
  -post_selected,
  average_post_fdr,
  -absolute_average_post_log2FC,
  met_id
)
selected_results <- feature_results[post_selected == TRUE]

analysis_summary <- data.table(
  metric = c(
    "input_features", "analysis_filter_features", "post_differentiation_contrasts",
    "joint_FDR_cutoff", "joint_selected_features", "average_effect_FDR_cutoff",
    "average_fold_change_cutoff", "direction_requirement", "post_selected_features",
    "sustained_lower_in_middle_late_vs_early_features",
    "sustained_higher_in_middle_late_vs_early_features"
  ),
  value = as.character(c(
    nrow(feature_qc_table), length(analysis_feature_ids), ncol(contrast_matrix),
    0.05, sum(feature_results$joint_selected, na.rm = TRUE), 0.05,
    1.5, "all eight contrasts must have the same sign", nrow(selected_results),
    sum(selected_results$trajectory_direction == "lower"),
    sum(selected_results$trajectory_direction == "higher")
  ))
)

fwrite(feature_results, file.path(table_directory, "02_post_differentiation_feature_results.csv"), bom = TRUE)
fwrite(selected_results, file.path(table_directory, "02_selected_consistent_direction_features.csv"), bom = TRUE)
fwrite(analysis_summary, file.path(table_directory, "02_differential_analysis_summary.csv"), bom = TRUE)

cat("Code 2 completed successfully.\n")
cat("Features tested:", length(analysis_feature_ids), "\n")
cat("Joint F-test selected features:", sum(feature_results$joint_selected, na.rm = TRUE), "\n")
cat("Final consistent-direction features:", nrow(selected_results), "\n")
cat("  Lower in Middle/Late vs Early:", sum(selected_results$trajectory_direction == "lower"), "\n")
cat("  Higher in Middle/Late vs Early:", sum(selected_results$trajectory_direction == "higher"), "\n")
