# Selected-feature PCA and differentiation-day-supervised PLS-DA.

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260905)

if (!requireNamespace("mixOmics", quietly = TRUE)) {
  stop("Missing R package: mixOmics")
}

code2_objects <- c(
  "analysis_abundance", "analysis_log2", "analysis_pareto",
  "sample_metadata", "selected_results", "analysis_root"
)
if (!all(vapply(code2_objects, exists, logical(1), inherits = TRUE))) {
  code2_candidates <- c(
    "02_post_differentiation_limma.R",
    file.path("RStudio", "02_post_differentiation_limma.R")
  )
  code2_file <- code2_candidates[file.exists(code2_candidates)][1]
  if (is.na(code2_file)) stop("Run Code 2 first or open this script from the RStudio project")
  source(code2_file)
}

figure_directory <- file.path(analysis_root, "Figure")
table_directory <- file.path(analysis_root, "Tables")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

passage_colors <- c(
  "Early Passage" = "#3B6FB6",
  "Middle Passage" = "#E69F00",
  "Late Passage" = "#B24C63"
)
day_shapes <- c(D0 = 21, D1 = 22, D3 = 24, D6 = 25, D9 = 8)

selected_ids <- selected_results$met_id
selected_pareto <- pareto_scale(
  analysis_log2[selected_ids, , drop = FALSE]
)
selected_pca_model <- prcomp(
  t(selected_pareto),
  center = FALSE,
  scale. = FALSE
)
selected_pca_variance <- (
  100 * selected_pca_model$sdev^2 / sum(selected_pca_model$sdev^2)
)
selected_pca_scores <- as.data.table(
  selected_pca_model$x[, 1:5, drop = FALSE],
  keep.rownames = "biological_sample"
)
selected_pca_scores <- sample_metadata[
  selected_pca_scores,
  on = "biological_sample"
]

selected_pca_figure <- ggplot(
  selected_pca_scores,
  aes(PC1, PC2, color = passage, shape = day)
) +
  geom_point(size = 3.2, stroke = 0.9, fill = "white") +
  scale_color_manual(values = passage_colors) +
  scale_shape_manual(values = day_shapes) +
  labs(
    x = sprintf("PC1 (%.1f%%)", selected_pca_variance[1]),
    y = sprintf("PC2 (%.1f%%)", selected_pca_variance[2]),
    color = "Passage",
    shape = "Day"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold")
  )

selected_pca_by_day <- ggplot(
  selected_pca_scores,
  aes(PC1, PC2, color = passage)
) +
  geom_point(shape = 16, size = 2.8) +
  facet_wrap(~day, nrow = 1) +
  scale_color_manual(values = passage_colors) +
  labs(
    x = sprintf("PC1 (%.1f%%)", selected_pca_variance[1]),
    y = sprintf("PC2 (%.1f%%)", selected_pca_variance[2]),
    color = "Passage"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold")
  )

prepare_split <- function(train_index, test_index) {
  train <- analysis_abundance[, train_index, drop = FALSE]
  test <- analysis_abundance[, test_index, drop = FALSE]
  train_min <- apply(train, 1, function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive)) min(positive) else NA_real_
  })
  keep <- is.finite(train_min)
  train <- train[keep, , drop = FALSE]
  test <- test[keep, , drop = FALSE]
  train_min <- train_min[keep]

  for (i in seq_len(nrow(train))) {
    train[i, !is.finite(train[i, ]) | train[i, ] <= 0] <- train_min[i] / 2
    test[i, !is.finite(test[i, ]) | test[i, ] <= 0] <- train_min[i] / 2
  }

  train <- log2(train)
  test <- log2(test)
  feature_mean <- rowMeans(train)
  feature_sd <- apply(train, 1, sd)
  keep <- is.finite(feature_sd) & feature_sd > 0
  train <- sweep(train[keep, , drop = FALSE], 1, feature_mean[keep], "-")
  test <- sweep(test[keep, , drop = FALSE], 1, feature_mean[keep], "-")
  train <- sweep(train, 1, sqrt(feature_sd[keep]), "/")
  test <- sweep(test, 1, sqrt(feature_sd[keep]), "/")
  list(train = t(train), test = t(test))
}


# Same 485 QC-retained features and 45 samples; only supervised target changes.
stopifnot(identical(colnames(analysis_pareto), sample_metadata$biological_sample),
          nrow(analysis_pareto) == 485L, nrow(sample_metadata) == 45L)
plsda_group <- droplevels(sample_metadata$day)
stopifnot(identical(levels(plsda_group), c("D0", "D1", "D3", "D6", "D9")))
plsda_model <- mixOmics::plsda(t(analysis_pareto), plsda_group, ncomp = 2, scale = FALSE)
cat(sprintf("Day-target model fitted: component 1 %.4f%%; component 2 %.4f%% X variance\n",
            100 * plsda_model$prop_expl_var$X[1], 100 * plsda_model$prop_expl_var$X[2]))

make_folds <- function() {
  folds <- integer(nrow(sample_metadata))
  strata <- interaction(sample_metadata$passage, sample_metadata$day, drop = TRUE)
  for (stratum in levels(strata)) folds[strata == stratum] <- sample(1:3)
  folds
}
# Preprocessing depends on the fold, not on the permuted response labels.
prepare_cv <- function(folds) {
  lapply(1:3, function(fold) {
    test_index <- which(folds == fold)
    train_index <- which(folds != fold)
    list(train_index = train_index, test_index = test_index,
         matrices = prepare_split(train_index, test_index))
  })
}
cross_validated_predictions <- function(labels, prepared) {
  predictions <- matrix(NA_character_, nrow(sample_metadata), 2L)
  for (split in prepared) {
    model <- mixOmics::plsda(split$matrices$train,
                            droplevels(labels[split$train_index]),
                            ncomp = 2, scale = FALSE)
    predictions[split$test_index, ] <- predict(model, split$matrices$test,
      dist = "max.dist")$class$max.dist[, 1:2, drop = FALSE]
  }
  stopifnot(!anyNA(predictions))
  predictions
}
classification_metrics <- function(predictions, labels) {
  rbindlist(lapply(1:2, function(component) {
    predicted <- predictions[, component]
    sensitivity <- vapply(levels(labels), function(group) {
      mean(predicted[labels == group] == group)
    }, numeric(1))
    result <- data.table(component = component, accuracy = mean(predicted == labels),
                         balanced_accuracy = mean(sensitivity))
    for (group in levels(labels)) result[[paste0(group, "_sensitivity")]] <- sensitivity[group]
    result
  }))
}
set.seed(20260905)
cv_results <- rbindlist(lapply(seq_len(100L), function(repetition) {
  result <- classification_metrics(
    cross_validated_predictions(plsda_group, prepare_cv(make_folds())), plsda_group)
  result[, repetition := repetition]
  if (repetition %% 25L == 0L) cat("CV repetitions completed:", repetition, "/100\n")
  result
}))
cv_summary <- cv_results[, lapply(.SD, mean), by = component,
  .SDcols = c("accuracy", "balanced_accuracy", paste0(levels(plsda_group), "_sensitivity"))]
fixed_folds <- make_folds()
fixed_prepared <- prepare_cv(fixed_folds)
fixed_observed <- classification_metrics(
  cross_validated_predictions(plsda_group, fixed_prepared), plsda_group
)[component == 2, balanced_accuracy]

# Test day-label information while retaining passage membership.
# Full-class counts are preserved within each passage; CV folds stay fixed.
permutation_results <- rbindlist(lapply(seq_len(999L), function(permutation) {
  permuted_group <- plsda_group
  for (passage_level in levels(sample_metadata$passage)) {
    index <- which(sample_metadata$passage == passage_level)
    permuted_group[index] <- sample(permuted_group[index])
  }
  result <- data.table(permutation = permutation, balanced_accuracy =
    classification_metrics(cross_validated_predictions(permuted_group, fixed_prepared),
                           permuted_group)[component == 2, balanced_accuracy])
  if (permutation %% 100L == 0L) cat("Permutations completed:", permutation, "/999\n")
  result
}))
permutation_p <- (1 + sum(permutation_results$balanced_accuracy >= fixed_observed)) / 1000

plsda_scores <- copy(sample_metadata)
plsda_scores[, pls_component_1 := plsda_model$variates$X[, 1]]
plsda_scores[, pls_component_2 := plsda_model$variates$X[, 2]]
axis1 <- sprintf("PLS component 1 (%.1f%% X variance)", 100 * plsda_model$prop_expl_var$X[1])
axis2 <- sprintf("PLS component 2 (%.1f%% X variance)", 100 * plsda_model$prop_expl_var$X[2])
plsda_figure <- ggplot(plsda_scores,
  aes(pls_component_1, pls_component_2, color = passage, shape = day)) +
  # Descriptive passage ellipses use the same 95% t setting as the original plot.
  # The supervised response remains differentiation day; scores are unchanged.
  stat_ellipse(aes(x = pls_component_1, y = pls_component_2,
                   color = passage, group = passage),
               inherit.aes = FALSE, type = "t", level = 0.95,
               linewidth = 0.65, show.legend = FALSE) +
  geom_point(size = 3.2, stroke = 0.9, fill = "white") +
  scale_color_manual(values = passage_colors) +
  scale_shape_manual(values = day_shapes) +
  labs(title = "PLS-DA: differentiation day", x = axis1, y = axis2,
       color = "Passage", shape = "Day") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"))
plsda_by_day <- ggplot(plsda_scores, aes(pls_component_1, pls_component_2, color = passage)) +
  geom_point(shape = 16, size = 2.8) + facet_wrap(~day, nrow = 1) +
  scale_color_manual(values = passage_colors) +
  labs(x = axis1, y = axis2, color = "Passage") + theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"))
save_figure <- function(name, plot, width, height) {
  ggsave(file.path(figure_directory, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(figure_directory, paste0(name, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
}
save_figure("05_selected_feature_PCA", selected_pca_figure, 7.2, 5.2)
save_figure("05_selected_feature_PCA_by_day", selected_pca_by_day, 11.5, 3.5)
save_figure("05_day_PLSDA", plsda_figure, 7.2, 5.2)
save_figure("05_day_PLSDA_by_day", plsda_by_day, 11.5, 3.5)

summary_values <- list(
 selected_PCA_features = length(selected_ids),
 selected_PCA_PC1_percent = selected_pca_variance[1],
 selected_PCA_PC2_percent = selected_pca_variance[2],
 PLSDA_target = "Differentiation day: D0 vs D1 vs D3 vs D6 vs D9",
 PLSDA_features = nrow(analysis_pareto), PLSDA_samples = nrow(sample_metadata),
 PLSDA_components = 2,
 PLSDA_component_1_X_variance_percent = 100 * plsda_model$prop_expl_var$X[1],
 PLSDA_component_2_X_variance_percent = 100 * plsda_model$prop_expl_var$X[2],
 cross_validation = "3-fold balanced within passage x day; preprocessing estimated in training folds",
 cross_validation_repetitions = 100,
 component_1_mean_balanced_accuracy = cv_summary[component == 1, balanced_accuracy],
 component_2_mean_accuracy = cv_summary[component == 2, accuracy],
 component_2_mean_balanced_accuracy = cv_summary[component == 2, balanced_accuracy],
 component_2_sd_balanced_accuracy = cv_results[component == 2, sd(balanced_accuracy)],
 fixed_fold_observed_balanced_accuracy = fixed_observed,
 permutations = 999,
 permutation_statistic = "Day labels permuted within passage; fixed-fold balanced accuracy",
 permutation_p_value = permutation_p,
 seed = 20260905,
 mixOmics_version = as.character(packageVersion("mixOmics")),
 input_file = normalizePath(input_file, winslash = "/"),
 input_md5 = unname(tools::md5sum(input_file)),
 input_scope = "All 485 QC-retained features; not the 60 differential features",
 detection_filter_CV_scope = "485-feature detection filter fixed from complete input; imputation/log2/Pareto fitted in training folds"
)
for (group in levels(plsda_group)) {
 summary_values[[paste0("component_2_", group, "_sensitivity")]] <-
   cv_summary[component == 2][[paste0(group, "_sensitivity")]]
}
validation_summary <- data.table(metric = names(summary_values),
                                 value = vapply(summary_values, as.character, character(1)))
fwrite(selected_pca_scores, file.path(table_directory, "05_selected_PCA_scores.csv"), bom = TRUE)
fwrite(plsda_scores, file.path(table_directory, "05_day_PLSDA_scores.csv"), bom = TRUE)
fwrite(cv_results, file.path(table_directory, "05_day_PLSDA_cross_validation_results.csv"), bom = TRUE)
fwrite(validation_summary, file.path(table_directory, "05_multivariate_analysis_summary.csv"), bom = TRUE)
fwrite(permutation_results, file.path(table_directory, "05_day_PLSDA_permutation_results.csv"), bom = TRUE)
for (ext in c("png", "pdf")) file.copy(
 file.path(figure_directory, paste0("05_day_PLSDA.", ext)),
 file.path(analysis_root, "Critical", paste0("05_day_PLSDA.", ext)), overwrite = TRUE)
cat("Code 5 day-target analysis completed.\n")
print(cv_summary)
cat("Permutation p:", permutation_p, "\n")
