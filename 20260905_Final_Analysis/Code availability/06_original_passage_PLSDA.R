# Reproduce the original passage-supervised PLS-DA and its requested display.
# Horizontal = negative original component 2; vertical = original component 1.
# Validation targets passage, with passage labels permuted within Day.
options(stringsAsFactors = FALSE, scipen = 999)
if (!requireNamespace("mixOmics", quietly = TRUE)) stop("Missing mixOmics")
needed <- c("analysis_abundance", "analysis_pareto", "sample_metadata", "analysis_root")
if (!all(vapply(needed, exists, logical(1), inherits = TRUE))) {
  candidates <- c("01_load_data_and_PCA.R", file.path("RStudio", "01_load_data_and_PCA.R"),
                  file.path("Code availability", "01_load_data_and_PCA.R"))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f)) stop("Run Code 01 from the analysis project")
  source(f)
}
set.seed(20260905)
figure_directory <- file.path(analysis_root, "Figure")
table_directory <- file.path(analysis_root, "Tables")
stopifnot(identical(colnames(analysis_pareto), sample_metadata$biological_sample),
          nrow(analysis_pareto) == 485L, nrow(sample_metadata) == 45L)
plsda_group <- droplevels(sample_metadata$passage)
plsda_model <- mixOmics::plsda(t(analysis_pareto), plsda_group, ncomp = 2, scale = FALSE)
stopifnot(identical(rownames(plsda_model$variates$X), sample_metadata$biological_sample))
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

make_folds <- function() {
  folds <- integer(nrow(sample_metadata))
  strata <- interaction(
    sample_metadata$passage,
    sample_metadata$day,
    drop = TRUE
  )
  for (stratum in levels(strata)) folds[strata == stratum] <- sample(1:3)
  folds
}

cross_validated_predictions <- function(labels, folds) {
  predictions <- matrix(NA_character_, nrow(sample_metadata), 2L)
  for (fold in 1:3) {
    test_index <- which(folds == fold)
    train_index <- which(folds != fold)
    split <- prepare_split(train_index, test_index)
    model <- mixOmics::plsda(
      split$train,
      droplevels(labels[train_index]),
      ncomp = 2,
      scale = FALSE
    )
    predictions[test_index, ] <- predict(
      model,
      split$test,
      dist = "max.dist"
    )$class$max.dist[, 1:2, drop = FALSE]
  }
  predictions
}

classification_metrics <- function(predictions, labels) {
  rbindlist(lapply(1:2, function(component) {
    predicted <- predictions[, component]
    sensitivity <- vapply(
      levels(labels),
      function(group) mean(predicted[labels == group] == group),
      numeric(1)
    )
    data.table(
      component = component,
      accuracy = mean(predicted == labels),
      balanced_accuracy = mean(sensitivity),
      early_sensitivity = sensitivity["Early Passage"],
      middle_sensitivity = sensitivity["Middle Passage"],
      late_sensitivity = sensitivity["Late Passage"]
    )
  }))
}

cv_results <- rbindlist(lapply(seq_len(100L), function(repetition) {
  result <- classification_metrics(
    cross_validated_predictions(plsda_group, make_folds()),
    plsda_group
  )
  result[, repetition := repetition]
  if (repetition %% 25L == 0L) cat("Passage CV:", repetition, "/100\n")
  result
}))

cv_summary <- cv_results[, .(
  mean_accuracy = mean(accuracy),
  sd_accuracy = sd(accuracy),
  mean_balanced_accuracy = mean(balanced_accuracy),
  sd_balanced_accuracy = sd(balanced_accuracy),
  mean_early_sensitivity = mean(early_sensitivity),
  mean_middle_sensitivity = mean(middle_sensitivity),
  mean_late_sensitivity = mean(late_sensitivity)
), by = component]

fixed_folds <- make_folds()
fixed_observed <- classification_metrics(
  cross_validated_predictions(plsda_group, fixed_folds),
  plsda_group
)[component == 2, balanced_accuracy]

permutation_results <- rbindlist(lapply(seq_len(999L), function(permutation) {
  permuted_group <- plsda_group
  for (day_level in levels(sample_metadata$day)) {
    index <- which(sample_metadata$day == day_level)
    permuted_group[index] <- sample(permuted_group[index])
  }
  if (permutation %% 100L == 0L) cat("Passage permutations:", permutation, "/999\n")
  data.table(
    permutation = permutation,
    balanced_accuracy = classification_metrics(
      cross_validated_predictions(permuted_group, fixed_folds),
      permuted_group
    )[component == 2, balanced_accuracy]
  )
}))

permutation_p <- (
  1 + sum(permutation_results$balanced_accuracy >= fixed_observed)
) / 1000

scores <- copy(sample_metadata)
scores[, pls_component_1 := plsda_model$variates$X[, 1]]
scores[, pls_component_2 := plsda_model$variates$X[, 2]]
scores[, horizontal := -pls_component_2]
scores[, vertical := pls_component_1]
passage_colors <- c("Early Passage" = "#3B6FB6", "Middle Passage" = "#E69F00",
                    "Late Passage" = "#B24C63")
day_shapes <- c(D0 = 21, D1 = 22, D3 = 24, D6 = 25, D9 = 8)
plot <- ggplot(scores, aes(horizontal, vertical, color = passage, shape = day)) +
  stat_ellipse(aes(horizontal, vertical, color = passage, group = passage),
               inherit.aes = FALSE, type = "t", level = 0.95,
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 3.2, stroke = 0.9, fill = "white") +
  scale_color_manual(values = passage_colors) +
  scale_shape_manual(values = day_shapes) +
  labs(title = "Original passage PLS-DA",
       x = sprintf("PLS component 2 (%.1f%% X variance)", 100 * plsda_model$prop_expl_var$X[2]),
       y = sprintf("PLS component 1 (%.1f%% X variance)", 100 * plsda_model$prop_expl_var$X[1]),
       color = "Passage", shape = "Day") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))
stem <- "Original_PLSDA_day_horizontal_passage_vertical"
ggsave(file.path(figure_directory, paste0(stem, ".png")), plot,
       width = 7.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(figure_directory, paste0(stem, ".svg")), plot, width = 7.5, height = 5.5)
for (ext in c("png", "svg")) file.copy(file.path(figure_directory, paste0(stem, ".", ext)),
  file.path(analysis_root, "Critical", paste0(stem, ".", ext)), overwrite = TRUE)
summary_values <- list(
  PLSDA_target = "Early Passage vs Middle Passage vs Late Passage",
  PLSDA_features = nrow(analysis_pareto), PLSDA_samples = nrow(sample_metadata),
  PLSDA_components = 2,
  PLSDA_component_1_X_variance_percent = 100 * plsda_model$prop_expl_var$X[1],
  PLSDA_component_2_X_variance_percent = 100 * plsda_model$prop_expl_var$X[2],
  displayed_horizontal = "negative original component 2; sign convention only",
  displayed_vertical = "original component 1",
  cross_validation = "balanced 3-fold within passage x day; train-fold preprocessing",
  cross_validation_repetitions = 100,
  component_2_mean_balanced_accuracy = cv_summary[component == 2, mean_balanced_accuracy],
  component_2_early_sensitivity = cv_summary[component == 2, mean_early_sensitivity],
  component_2_middle_sensitivity = cv_summary[component == 2, mean_middle_sensitivity],
  component_2_late_sensitivity = cv_summary[component == 2, mean_late_sensitivity],
  fixed_fold_observed_balanced_accuracy = fixed_observed,
  permutations = 999,
  permutation_statistic = "Passage labels permuted within Day; fixed-fold balanced accuracy",
  permutation_p_value = permutation_p,
  seed = 20260905, mixOmics_version = as.character(packageVersion("mixOmics")),
  input_md5 = unname(tools::md5sum(input_file)),
  detection_filter_CV_scope = "485-feature filter fixed from complete input; imputation/Pareto fitted in training folds")
summary_table <- data.table(metric = names(summary_values),
  value = vapply(summary_values, as.character, character(1)))
fwrite(scores, file.path(table_directory, "06_original_passage_PLSDA_scores.csv"), bom = TRUE)
fwrite(cv_results, file.path(table_directory, "06_original_passage_PLSDA_cross_validation_results.csv"), bom = TRUE)
fwrite(permutation_results, file.path(table_directory, "06_original_passage_PLSDA_permutation_results.csv"), bom = TRUE)
fwrite(summary_table, file.path(table_directory, "06_original_passage_PLSDA_summary.csv"), bom = TRUE)
cat("Original passage-target model and display reproduced.\n")
print(cv_summary)
