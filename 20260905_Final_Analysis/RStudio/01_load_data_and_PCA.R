# Load normalized lipidomics abundance, perform QC filtering and PCA.

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260905)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

working_directory <- getwd()
if (basename(working_directory) %in% c("RStudio", "Code availability")) {
  analysis_root <- dirname(working_directory)
} else if (
  dir.exists(file.path(working_directory, "Figure")) &&
    dir.exists(file.path(working_directory, "Tables"))
) {
  analysis_root <- working_directory
} else {
  stop("Open C2C12_Lipidomics_analysis.Rproj")
}

input_file <- file.path(
  dirname(analysis_root),
  "Lipidomics_normalized_abundance_positive.csv"
)
if (!file.exists(input_file)) stop("Input file was not found: ", input_file)

figure_directory <- file.path(analysis_root, "Figure")
table_directory <- file.path(analysis_root, "Tables")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

read_qi_export <- function(path) {
  x <- read.csv(
    path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = character(), fileEncoding = "UTF-8-BOM"
  )
  block_header <- as.character(x[1, , drop = TRUE])
  headers <- as.character(x[3, , drop = TRUE])
  data_rows <- x[-c(1:3), , drop = FALSE]
  normalized_start <- which(block_header == "Normalised abundance")
  raw_start <- which(block_header == "Raw abundance")
  annotation_start <- which(headers == "Accepted Compound ID")
  if (length(normalized_start) != 1L || length(raw_start) != 1L || length(annotation_start) != 1L) {
    stop("QI abundance blocks could not be identified")
  }
  normalized_columns <- normalized_start:(raw_start - 1L)
  raw_columns <- raw_start:(annotation_start - 1L)
  sample_names <- headers[normalized_columns]
  if (!identical(sample_names, headers[raw_columns])) {
    stop("Normalized and raw abundance sample columns differ")
  }
  numeric_column <- function(name) {
    index <- which(headers == name)
    if (length(index)) suppressWarnings(as.numeric(data_rows[[index[1]]])) else rep(NA_real_, nrow(data_rows))
  }
  character_column <- function(name) {
    index <- which(headers == name)
    if (length(index)) as.character(data_rows[[index[1]]]) else rep(NA_character_, nrow(data_rows))
  }
  feature_info <- data.table(
    met_id = sprintf("MET%03d", seq_len(nrow(data_rows))),
    feature = character_column("Compound"),
    neutral_mass = numeric_column("Neutral mass (Da)"),
    mz = numeric_column("m/z"),
    rt_min = numeric_column("Retention time (min)"),
    accepted_compound_id = character_column("Accepted Compound ID"),
    accepted_description = character_column("Accepted Description"),
    adducts = character_column("Adducts"),
    formula = character_column("Formula"),
    score = suppressWarnings(as.numeric(character_column("Score"))),
    fragmentation_score = suppressWarnings(as.numeric(character_column("Fragmentation Score")))
  )
  normalized_matrix <- sapply(
    normalized_columns,
    function(index) suppressWarnings(as.numeric(data_rows[[index]]))
  )
  rownames(normalized_matrix) <- feature_info$met_id
  colnames(normalized_matrix) <- sample_names
  list(feature_info = feature_info, normalized = normalized_matrix)
}

parse_samples <- function(sample_names) {
  pattern <- "_(P11|P22|P33)_D(0|1|3|6|9)_([0-9]{2})_([0-9]+)$"
  matches <- regexec(pattern, sample_names)
  parts <- regmatches(sample_names, matches)
  is_biological <- lengths(parts) == 5L
  is_qc <- grepl("_QC_[0-9]+$", sample_names)
  if (any(!is_biological & !is_qc)) stop("Unexpected LC-MS sample names")
  metadata <- data.table(
    input_injection = sample_names,
    sample_type = fifelse(is_qc, "QC", "Biological"),
    passage = NA_character_, day_numeric = NA_integer_,
    biological_replicate = NA_integer_, injection_replicate = NA_integer_
  )
  passage_lookup <- c(
    P11 = "Early Passage",
    P22 = "Middle Passage",
    P33 = "Late Passage"
  )
  metadata[is_biological, passage := passage_lookup[
    vapply(parts[is_biological], `[[`, character(1), 2)
  ]]
  metadata[is_biological, day_numeric := as.integer(vapply(parts[is_biological], `[[`, character(1), 3))]
  metadata[is_biological, biological_replicate := as.integer(vapply(parts[is_biological], `[[`, character(1), 4))]
  metadata[is_biological, injection_replicate := as.integer(vapply(parts[is_biological], `[[`, character(1), 5))]
  metadata[, biological_sample := fifelse(
    sample_type == "Biological",
    sprintf("%s_D%d_R%02d", gsub(" ", "_", passage), day_numeric, biological_replicate),
    NA_character_
  )]
  metadata[, day := fifelse(sample_type == "Biological", paste0("D", day_numeric), NA_character_)]
  metadata[, condition := fifelse(
    sample_type == "Biological",
    paste0(sub(" Passage$", "", passage), "_", day),
    NA_character_
  )]
  metadata[sample_type == "Biological", injection := sprintf(
    "%s_I%02d",
    biological_sample,
    injection_replicate
  )]
  metadata[sample_type == "QC", injection := sprintf("Pooled_QC_%02d", seq_len(.N))]
  metadata
}

collapse_injections <- function(matrix, injection_metadata) {
  biological <- injection_metadata[sample_type == "Biological"]
  metadata <- unique(biological[, .(
    biological_sample, passage, day_numeric, day, condition, biological_replicate
  )])
  metadata[, passage_order := match(passage, c("Early Passage", "Middle Passage", "Late Passage"))]
  setorder(metadata, passage_order, day_numeric, biological_replicate)
  metadata[, passage_order := NULL]
  collapsed <- sapply(metadata$biological_sample, function(sample_id) {
    columns <- biological[biological_sample == sample_id, injection]
    rowMeans(matrix[, columns, drop = FALSE], na.rm = TRUE)
  })
  rownames(collapsed) <- rownames(matrix)
  colnames(collapsed) <- metadata$biological_sample
  metadata[, passage := factor(
    passage,
    levels = c("Early Passage", "Middle Passage", "Late Passage")
  )]
  metadata[, day := factor(day, levels = c("D0", "D1", "D3", "D6", "D9"))]
  metadata[, condition := factor(
    condition,
    levels = as.vector(outer(c("Early", "Middle", "Late"), levels(day), paste, sep = "_"))
  )]
  list(matrix = collapsed, metadata = metadata)
}

half_min_log2 <- function(matrix) {
  answer <- matrix
  for (i in seq_len(nrow(answer))) {
    positive <- answer[i, is.finite(answer[i, ]) & answer[i, ] > 0]
    if (!length(positive)) {
      answer[i, ] <- NA_real_
    } else {
      answer[i, !is.finite(answer[i, ]) | answer[i, ] <= 0] <- min(positive) / 2
      answer[i, ] <- log2(answer[i, ])
    }
  }
  answer
}

pareto_scale <- function(matrix) {
  feature_means <- rowMeans(matrix, na.rm = TRUE)
  feature_sd <- apply(matrix, 1, sd, na.rm = TRUE)
  keep <- is.finite(feature_sd) & feature_sd > 0
  centered <- sweep(matrix[keep, , drop = FALSE], 1, feature_means[keep], "-")
  sweep(centered, 1, sqrt(feature_sd[keep]), "/")
}

qi_data <- read_qi_export(input_file)
injection_metadata <- parse_samples(colnames(qi_data$normalized))
colnames(qi_data$normalized) <- injection_metadata$injection
collapsed_data <- collapse_injections(qi_data$normalized, injection_metadata)
biological_matrix <- collapsed_data$matrix
sample_metadata <- collapsed_data$metadata

if (nrow(biological_matrix) != 867L) stop("Expected 867 LC-MS features")
if (ncol(biological_matrix) != 45L) stop("Expected 45 biological samples")

qc_columns <- injection_metadata[sample_type == "QC", injection]
qc_matrix <- qi_data$normalized[, qc_columns, drop = FALSE]
qc_detected_n <- rowSums(qc_matrix > 0, na.rm = TRUE)
qc_detection_fraction <- qc_detected_n / ncol(qc_matrix)
qc_mean_positive <- apply(qc_matrix, 1, function(x) mean(x[is.finite(x) & x > 0]))
qc_sd_positive <- apply(qc_matrix, 1, function(x) sd(x[is.finite(x) & x > 0]))
qc_cv_percent <- 100 * qc_sd_positive / qc_mean_positive

condition_detection <- sapply(levels(sample_metadata$condition), function(condition_name) {
  columns <- sample_metadata[condition == condition_name, biological_sample]
  rowMeans(biological_matrix[, columns, drop = FALSE] > 0, na.rm = TRUE)
})
maximum_condition_detection <- apply(condition_detection, 1, max, na.rm = TRUE)
analysis_filter <- qc_detection_fraction >= 0.80 & maximum_condition_detection >= 2 / 3
names(analysis_filter) <- rownames(biological_matrix)
analysis_feature_ids <- names(analysis_filter)[analysis_filter]

analysis_abundance <- biological_matrix[analysis_feature_ids, , drop = FALSE]
analysis_log2 <- half_min_log2(analysis_abundance)
analysis_pareto <- pareto_scale(analysis_log2)

pca_model <- prcomp(t(analysis_pareto), center = FALSE, scale. = FALSE)
pca_variance_percent <- 100 * pca_model$sdev^2 / sum(pca_model$sdev^2)
pca_scores <- as.data.table(pca_model$x[, 1:5, drop = FALSE], keep.rownames = "biological_sample")
pca_scores <- sample_metadata[pca_scores, on = "biological_sample"]

passage_colors <- c(
  "Early Passage" = "#3B6FB6",
  "Middle Passage" = "#E69F00",
  "Late Passage" = "#B24C63"
)
day_shapes <- c(D0 = 21, D1 = 22, D3 = 24, D6 = 25, D9 = 8)
pca_figure <- ggplot(pca_scores, aes(PC1, PC2, color = passage, shape = day)) +
  geom_point(size = 3.2, stroke = 0.9, fill = "white") +
  scale_color_manual(values = passage_colors) +
  scale_shape_manual(values = day_shapes) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pca_variance_percent[1]),
    y = sprintf("PC2 (%.1f%%)", pca_variance_percent[2]),
    color = "Passage", shape = "Day"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold")
  )

ggsave(file.path(figure_directory, "01_lipidomics_PCA_biological_samples.png"), pca_figure, width = 6.5, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(figure_directory, "01_lipidomics_PCA_biological_samples.pdf"), pca_figure, width = 6.5, height = 5.2, device = cairo_pdf)

combined_abundance <- cbind(analysis_abundance, qc_matrix[analysis_feature_ids, , drop = FALSE])
combined_log2 <- half_min_log2(combined_abundance)
combined_pareto <- pareto_scale(combined_log2)
qc_pca_model <- prcomp(t(combined_pareto), center = FALSE, scale. = FALSE)
qc_pca_variance_percent <- 100 * qc_pca_model$sdev^2 / sum(qc_pca_model$sdev^2)
qc_pca_scores <- data.table(
  sample = rownames(qc_pca_model$x),
  PC1 = qc_pca_model$x[, 1],
  PC2 = qc_pca_model$x[, 2],
  sample_type = c(rep("Biological", ncol(analysis_abundance)), rep("QC", ncol(qc_matrix))),
  group = c(as.character(sample_metadata$passage), rep("QC", ncol(qc_matrix)))
)
qc_pca_scores[, group := factor(
  group,
  levels = c("Early Passage", "Middle Passage", "Late Passage", "QC")
)]

qc_pca_figure <- ggplot(qc_pca_scores, aes(PC1, PC2, color = group, shape = sample_type)) +
  geom_point(size = 2.7, alpha = 0.85) +
  scale_color_manual(values = c(passage_colors, QC = "#222222")) +
  scale_shape_manual(values = c(Biological = 16, QC = 4)) +
  labs(
    x = sprintf("PC1 (%.1f%%)", qc_pca_variance_percent[1]),
    y = sprintf("PC2 (%.1f%%)", qc_pca_variance_percent[2]),
    color = "Group", shape = "Sample type"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold")
  )

ggsave(file.path(figure_directory, "01_lipidomics_PCA_with_pooled_QC.png"), qc_pca_figure, width = 8, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(figure_directory, "01_lipidomics_PCA_with_pooled_QC.pdf"), qc_pca_figure, width = 8, height = 5.2, device = cairo_pdf)

feature_qc_table <- data.table(
  met_id = rownames(biological_matrix),
  qc_detected_n = qc_detected_n,
  qc_injections = ncol(qc_matrix),
  qc_detection_fraction = qc_detection_fraction,
  qc_cv_percent = qc_cv_percent,
  overall_biological_detection_fraction = rowMeans(biological_matrix > 0, na.rm = TRUE),
  maximum_condition_detection_fraction = maximum_condition_detection,
  analysis_filter_pass = analysis_filter
)
feature_qc_table <- merge(
  qi_data$feature_info, feature_qc_table,
  by = "met_id", all.x = TRUE, sort = FALSE
)

pca_summary <- data.table(
  metric = c(
    "input_format", "ionization_mode",
    "input_features", "biological_samples", "pooled_QC_injections",
    "technical_injections_per_biological_sample", "analysis_filter_features",
    "QC_detection_cutoff", "biological_detection_cutoff", "QC_CV_used_as_filter",
    "biological_PCA_PC1_percent", "biological_PCA_PC2_percent",
    "pooled_QC_PCA_PC1_percent", "pooled_QC_PCA_PC2_percent"
  ),
  value = as.character(c(
    "Progenesis QI normalized abundance", "positive",
    nrow(biological_matrix), ncol(biological_matrix), ncol(qc_matrix), 3,
    length(analysis_feature_ids), 0.80, 2 / 3, FALSE,
    pca_variance_percent[1:2], qc_pca_variance_percent[1:2]
  ))
)

fwrite(
  injection_metadata[, !"input_injection"],
  file.path(table_directory, "01_injection_metadata.csv"),
  bom = TRUE
)
fwrite(sample_metadata, file.path(table_directory, "01_sample_metadata.csv"), bom = TRUE)
fwrite(feature_qc_table, file.path(table_directory, "01_feature_filter_all_features.csv"), bom = TRUE)
fwrite(feature_qc_table[analysis_filter_pass == TRUE], file.path(table_directory, "01_expression_filter_retained_features.csv"), bom = TRUE)
fwrite(pca_scores, file.path(table_directory, "01_PCA_scores_biological.csv"), bom = TRUE)
fwrite(qc_pca_scores, file.path(table_directory, "01_PCA_scores_with_pooled_QC.csv"), bom = TRUE)
fwrite(pca_summary, file.path(table_directory, "01_PCA_summary.csv"), bom = TRUE)

cat("Code 1 completed successfully.\n")
cat("Input features:", nrow(biological_matrix), "\n")
cat("Biological samples:", ncol(biological_matrix), "\n")
cat("Pooled-QC injections:", ncol(qc_matrix), "\n")
cat("Features after QC/detection filter:", length(analysis_feature_ids), "\n")
cat(sprintf("Biological PCA PC1: %.1f%%; PC2: %.1f%%\n", pca_variance_percent[1], pca_variance_percent[2]))
