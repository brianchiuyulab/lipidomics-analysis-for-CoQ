# Global heatmap of all lipidomics features selected in Code 2.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("data.table", "ComplexHeatmap", "circlize")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

# Load Code 2 objects when needed.
code2_objects <- c(
  "analysis_log2", "model_metadata", "selected_results", "analysis_root"
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

selected_feature_ids <- intersect(selected_results$met_id, rownames(analysis_log2))
heatmap_matrix <- analysis_log2[selected_feature_ids, , drop = FALSE]
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), , drop = FALSE]
heatmap_matrix[heatmap_matrix > 2] <- 2
heatmap_matrix[heatmap_matrix < -2] <- -2

heatmap_metadata <- copy(model_metadata)
heatmap_metadata <- heatmap_metadata[match(colnames(heatmap_matrix), biological_sample)]
setorder(heatmap_metadata, passage, day_numeric, biological_replicate)
heatmap_matrix <- heatmap_matrix[
  ,
  heatmap_metadata$biological_sample,
  drop = FALSE
]

heatmap_feature_table <- selected_results[
  match(rownames(heatmap_matrix), met_id),
  .(
    met_id, trajectory_direction,
    average_post_log2FC, average_post_fdr
  )
]

annotation_file <- file.path(table_directory, "03_lipid_feature_annotation_audit.csv")
if (!file.exists(annotation_file)) stop("Run Code 3 before Code 4")
annotation_audit <- fread(annotation_file)
heatmap_feature_table <- annotation_audit[heatmap_feature_table, on = "met_id"]
if (anyNA(heatmap_feature_table$annotation_evidence)) stop("Missing annotation audit rows")
heatmap_feature_table[, annotation_status := fifelse(
  startsWith(annotation_evidence, "Putative"),
  "Putatively annotated",
  "Unassigned/ambiguous"
)]
heatmap_feature_table[, figure_lipid_class := fcase(
  broad_class == "Glycerolipid candidate", "Glycerolipid",
  broad_class == "Phosphatidylcholine candidate", "Phosphatidylcholine",
  broad_class == "Sphingolipid candidate", "Sphingolipid",
  default = broad_class
)]
heatmap_feature_table[, figure_label := fifelse(
  annotation_status == "Putatively annotated" & nzchar(selected_putative_name),
  paste0(met_id, " (putative ", selected_putative_name, ")"),
  met_id
)]

feature_annotation <- heatmap_feature_table[match(rownames(heatmap_matrix), met_id)]
lipid_class_palette <- c(
  "Free fatty acid" = "#E6AB02",
  "Glycerolipid" = "#D95F02",
  "Lysophosphatidylcholine" = "#1B9E77",
  "Phosphatidylcholine" = "#E7298A",
  "Phosphatidylethanolamine" = "#7570B3",
  "Phosphatidylglycerol" = "#66A61E",
  "Sphingolipid" = "#1F9E89",
  "Sphingomyelin" = "#E76F51",
  "Ubiquinone" = "#6A3D9A",
  "Unassigned" = "#BDBDBD"
)
row_annotation <- rowAnnotation(
  `Putative lipid class` = feature_annotation$figure_lipid_class,
  col = list(
    `Putative lipid class` = lipid_class_palette
  ),
  simple_anno_size = grid::unit(3.5, "mm")
)

passage_colors <- c(
  "Early Passage" = "#3B6FB6",
  "Middle Passage" = "#E69F00",
  "Late Passage" = "#B24C63"
)
day_colors <- c(D0 = "#D9D9D9", D1 = "#7FC97F", D3 = "#BEAED4", D6 = "#FDCDAC", D9 = "#FFFF99")
column_annotation <- HeatmapAnnotation(
  Day = heatmap_metadata$day,
  Passage = heatmap_metadata$passage_group,
  col = list(Day = day_colors, Passage = passage_colors),
  annotation_name_gp = grid::gpar(fontsize = 9)
)

heatmap_object <- row_annotation + Heatmap(
  heatmap_matrix,
  name = "Row Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  top_annotation = column_annotation,
  column_split = heatmap_metadata$passage_group,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  show_row_names = TRUE,
  row_labels = feature_annotation$figure_label,
  row_names_gp = grid::gpar(fontsize = 6.5, fontface = "bold"),
  show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 6.5, fontface = "bold"),
  column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE
)

png(
  file.path(figure_directory, "04_all_selected_features_heatmap.png"),
  width = 17, height = 11.5, units = "in", res = 300, type = "cairo-png"
)
drawn_heatmap <- draw(
  heatmap_object,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

pdf(
  file.path(figure_directory, "04_all_selected_features_heatmap.pdf"),
  width = 17, height = 11.5, useDingbats = FALSE
)
draw(heatmap_object, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

row_order_index <- unlist(row_order(drawn_heatmap), use.names = FALSE)
column_order_index <- unlist(column_order(drawn_heatmap), use.names = FALSE)
heatmap_row_order <- data.table(
  heatmap_order = seq_along(row_order_index),
  met_id = rownames(heatmap_matrix)[row_order_index]
)
heatmap_column_order <- data.table(
  heatmap_order = seq_along(column_order_index),
  biological_sample = colnames(heatmap_matrix)[column_order_index]
)
heatmap_column_order <- merge(
  heatmap_column_order,
  heatmap_metadata[, .(biological_sample, passage, passage_group, day, biological_replicate)],
  by = "biological_sample",
  all.x = TRUE,
  sort = FALSE
)
setorder(heatmap_column_order, heatmap_order)

heatmap_feature_table[, heatmap_order := match(met_id, heatmap_row_order$met_id)]
setorder(heatmap_feature_table, heatmap_order)
setcolorder(
  heatmap_feature_table,
  c("heatmap_order", setdiff(names(heatmap_feature_table), "heatmap_order"))
)

fwrite(heatmap_feature_table, file.path(table_directory, "04_heatmap_feature_table.csv"), bom = TRUE)
fwrite(heatmap_row_order, file.path(table_directory, "04_heatmap_feature_order.csv"), bom = TRUE)
fwrite(heatmap_column_order, file.path(table_directory, "04_heatmap_sample_order.csv"), bom = TRUE)

analysis_summary <- data.table(
  metric = c(
    "heatmap_features", "heatmap_samples", "row_clustering",
    "column_clustering", "column_order", "row_labels", "row_annotation"
  ),
  value = as.character(c(
    nrow(heatmap_matrix), ncol(heatmap_matrix), "complete linkage",
    "not clustered", "passage, day, biological replicate",
    "MET identifiers with putative compound names where assigned",
    "putative lipid class"
  ))
)
fwrite(analysis_summary, file.path(table_directory, "04_global_heatmap_summary.csv"), bom = TRUE)

cat("Code 4 completed successfully.\n")
cat("Heatmap features:", nrow(heatmap_matrix), "\n")
cat("Heatmap samples:", ncol(heatmap_matrix), "\n")
