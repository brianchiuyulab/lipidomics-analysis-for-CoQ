# Systematic annotation audit for all lipidomics features selected in Code 2.

options(stringsAsFactors = FALSE, scipen = 999, timeout = 60)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("data.table", "curl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(curl)
})

working_directory <- getwd()
if (basename(working_directory) == "RStudio") {
  analysis_root <- dirname(working_directory)
} else if (dir.exists(file.path(working_directory, "Tables"))) {
  analysis_root <- working_directory
} else {
  stop("Open C2C12_Lipidomics_analysis.Rproj")
}

table_directory <- file.path(analysis_root, "Tables")
input_file <- file.path(table_directory, "02_selected_consistent_direction_features.csv")
if (!file.exists(input_file)) stop("Run Code 2 first")
features <- fread(input_file)

ppm_tolerance <- 5
common_positive_adducts <- c("M+H", "M+NH4", "M+Na", "M+K")

fetch_text <- function(url, attempts = 3L) {
  for (attempt in seq_len(attempts)) {
    response <- try(curl_fetch_memory(url), silent = TRUE)
    if (!inherits(response, "try-error") && response$status_code == 200L) {
      return(rawToChar(response$content))
    }
    Sys.sleep(attempt)
  }
  warning("Query failed: ", url)
  ""
}

parse_tsv <- function(text, strip_pre = FALSE) {
  if (!nzchar(text)) return(data.table())
  if (strip_pre) text <- gsub("</?pre>", "", text, ignore.case = TRUE)
  text <- trimws(text)
  if (!nzchar(text) || !grepl("\\n", text)) return(data.table())
  tryCatch(fread(text = text, sep = "\t", fill = TRUE, quote = ""), error = function(e) data.table())
}

search_neutral_mass <- function(met_id, neutral_mass) {
  tolerance_da <- neutral_mass * ppm_tolerance / 1000000
  url <- paste0(
    "https://www.lipidmaps.org/data/structure/LMSDSearch.php?",
    "Mode=ProcessTextSearch&OutputMode=File&ExactMass=",
    sprintf("%.8f", neutral_mass),
    "&ExactMassOffSet=", sprintf("%.8f", tolerance_da)
  )
  result <- parse_tsv(fetch_text(url))
  if (!nrow(result)) return(data.table())
  setnames(result, names(result), tolower(names(result)))
  result[, candidate_name := fifelse(
    !is.na(common_name) & common_name != "" & common_name != "-",
    common_name,
    systematic_name
  )]
  result[, `:=`(
    met_id = met_id,
    query_basis = "neutral exact mass",
    query_mass = neutral_mass,
    query_adduct = NA_character_,
    matched_mass = as.numeric(mass),
    delta_da = as.numeric(mass) - neutral_mass,
    delta_ppm = 1000000 * (as.numeric(mass) - neutral_mass) / neutral_mass,
    source_database = "LIPID MAPS LMSD"
  )]
  result[, .(
    met_id, query_basis, query_mass, query_adduct, matched_mass,
    delta_da, delta_ppm, candidate_name, formula,
    category, main_class, sub_class, lm_id, source_database
  )]
}

search_mz_adduct <- function(met_id, mz, adduct) {
  tolerance_da <- mz * ppm_tolerance / 1000000
  encoded_adduct <- utils::URLencode(adduct, reserved = TRUE)
  url <- paste0(
    "https://www.lipidmaps.org/rest/moverz/LMSD/",
    sprintf("%.8f", mz), "/", encoded_adduct, "/",
    sprintf("%.8f", tolerance_da), "/txt"
  )
  result <- parse_tsv(fetch_text(url), strip_pre = TRUE)
  if (!nrow(result)) return(data.table())
  setnames(result, names(result), make.names(tolower(names(result))))
  result[, `:=`(
    met_id = met_id,
    query_basis = "observed m/z plus adduct",
    query_mass = mz,
    query_adduct = adduct,
    matched_mass = as.numeric(matched.m.z),
    delta_da = as.numeric(delta),
    delta_ppm = 1000000 * as.numeric(delta) / mz,
    candidate_name = name,
    source_database = "LIPID MAPS LMSD"
  )]
  if (!"lm_id" %in% names(result)) result[, lm_id := NA_character_]
  result[, .(
    met_id, query_basis, query_mass, query_adduct, matched_mass,
    delta_da, delta_ppm, candidate_name, formula,
    category, main.class, sub.class, lm_id, source_database
  )][, setnames(.SD, c("main.class", "sub.class"), c("main_class", "sub_class"))]
}

candidate_cache_file <- file.path(table_directory, "03_all_LMSD_mass_candidates.csv")
if (file.exists(candidate_cache_file)) {
  candidates <- fread(candidate_cache_file)
} else {
  candidate_list <- vector("list", nrow(features))
  for (i in seq_len(nrow(features))) {
    feature_hits <- list()
    if (!is.na(features$neutral_mass[i])) {
      feature_hits[[length(feature_hits) + 1L]] <- search_neutral_mass(
        features$met_id[i], features$neutral_mass[i]
      )
    }
    reported_adducts <- trimws(unlist(strsplit(features$adducts[i], ",", fixed = TRUE)))
    reported_adducts <- reported_adducts[nzchar(reported_adducts)]
    adducts_to_search <- if (length(reported_adducts)) reported_adducts else common_positive_adducts
    for (adduct in adducts_to_search) {
      feature_hits[[length(feature_hits) + 1L]] <- search_mz_adduct(
        features$met_id[i], features$mz[i], adduct
      )
    }
    candidate_list[[i]] <- rbindlist(feature_hits, fill = TRUE)
    Sys.sleep(0.05)
  }
  candidates <- unique(rbindlist(candidate_list, fill = TRUE))
}
if (nrow(candidates)) {
  candidates[, absolute_delta_ppm := abs(delta_ppm)]
  candidates[, candidate_review_flag := fcase(
    grepl("-d[0-9]+($|[^:])|13C|15N|18O", candidate_name),
    "Excluded from automatic selection: isotope-labelled compound",
    default = "Eligible mass-only candidate"
  )]
  candidates[, eligible_for_unique_name := candidate_review_flag == "Eligible mass-only candidate"]
  candidates[, candidate_broad_class := fcase(
    grepl("CoQ|Coenzyme Q|Ubiquinone", candidate_name, ignore.case = TRUE), "Ubiquinone",
    grepl("Glycerophosphocholine|Phosphatidylcholine", main_class, ignore.case = TRUE), "Phosphatidylcholine candidate",
    grepl("Sphingomyelin", main_class, ignore.case = TRUE), "Sphingomyelin candidate",
    grepl("Diradylglycerol|Glycerolipid", main_class, ignore.case = TRUE), "Glycerolipid candidate",
    grepl("Ceramide|Sphingolipid", main_class, ignore.case = TRUE), "Sphingolipid candidate",
    grepl("Fatty Acid|Fatty acid", category), "Fatty-acyl candidate",
    default = "Unassigned"
  )]
  setorder(candidates, met_id, absolute_delta_ppm, candidate_name)
}

concise_qi_name <- function(compound_id, description) {
  if (is.na(compound_id) || !nzchar(compound_id)) return(NA_character_)
  if (grepl("^\\[[^]]+\\]", description)) {
    return(sub("^\\[([^]]+)\\].*$", "\\1", description))
  }
  if (grepl("^LPC", compound_id)) {
    chain <- regmatches(description, regexpr("[0-9]+:[0-9]+", description))
    if (length(chain) && nzchar(chain)) return(paste("LPC", chain))
  }
  if (grepl("^SM", compound_id)) {
    chains <- regmatches(description, regexpr("d[0-9]+:[0-9]+/[0-9]+:[0-9]+", description))
    if (length(chains) && nzchar(chains)) return(paste("SM", chains))
  }
  if (grepl("^PG", compound_id)) {
    chains <- regmatches(description, regexpr("[0-9]+:[0-9]+-[0-9]+:[0-9]+", description))
    if (length(chains) && nzchar(chains)) return(paste("PG", gsub("-", "/", chains, fixed = TRUE)))
  }
  if (grepl("^FFA", compound_id)) {
    chain <- regmatches(description, regexpr("C?[0-9]+:[0-9]+", description))
    if (length(chain) && nzchar(chain)) return(paste("FFA", sub("^C", "", chain)))
  }
  if (grepl("^PC\\(P-", description)) return(description)
  description
}

candidate_summary <- if (nrow(candidates)) candidates[, .(
  lmsd_candidate_count = .N,
  unique_candidate_names = uniqueN(candidate_name),
  eligible_candidate_count = sum(eligible_for_unique_name),
  candidate_names = paste(unique(candidate_name), collapse = " | "),
  eligible_candidate_names = paste(unique(candidate_name[eligible_for_unique_name]), collapse = " | "),
  candidate_categories = paste(unique(category), collapse = " | "),
  candidate_main_classes = paste(unique(main_class), collapse = " | "),
  candidate_sub_classes = paste(unique(sub_class), collapse = " | "),
  candidate_broad_classes = paste(unique(candidate_broad_class[eligible_for_unique_name]), collapse = " | "),
  unique_candidate_broad_class = if (
    uniqueN(candidate_broad_class[eligible_for_unique_name]) == 1L &&
      unique(candidate_broad_class[eligible_for_unique_name])[1] != "Unassigned"
  ) unique(candidate_broad_class[eligible_for_unique_name])[1] else NA_character_,
  best_absolute_delta_ppm = min(absolute_delta_ppm, na.rm = TRUE),
  has_neutral_mass_hit = any(query_basis == "neutral exact mass"),
  neutral_mass_unique_names = uniqueN(candidate_name[
    query_basis == "neutral exact mass" & eligible_for_unique_name
  ]),
  unique_neutral_candidate = if (
    uniqueN(candidate_name[query_basis == "neutral exact mass" & eligible_for_unique_name]) == 1L
  ) unique(candidate_name[query_basis == "neutral exact mass" & eligible_for_unique_name])[1] else NA_character_,
  unique_mz_candidate = if (
    uniqueN(candidate_name[query_basis == "observed m/z plus adduct" & eligible_for_unique_name]) == 1L
  ) unique(candidate_name[query_basis == "observed m/z plus adduct" & eligible_for_unique_name])[1] else NA_character_
), by = met_id] else data.table(met_id = character())

audit <- merge(features, candidate_summary, by = "met_id", all.x = TRUE, sort = FALSE)
audit[, qi_concise_name := mapply(concise_qi_name, accepted_compound_id, accepted_description)]
audit[, annotation_evidence := fcase(
  !is.na(accepted_compound_id) & accepted_compound_id != "", "Putative QI annotation",
  !is.na(unique_neutral_candidate), "Putative LMSD candidate from QI-inferred neutral mass",
  !is.na(unique_mz_candidate), "Putative LMSD candidate from m/z-adduct search",
  !is.na(lmsd_candidate_count) & lmsd_candidate_count > 0 & eligible_candidate_count == 0,
  "Excluded isotope-labelled LMSD match",
  !is.na(lmsd_candidate_count) & lmsd_candidate_count > 0, "Ambiguous LMSD candidates",
  default = "Unassigned feature"
)]
audit[, selected_putative_name := fcase(
  annotation_evidence == "Putative QI annotation", qi_concise_name,
  annotation_evidence == "Putative LMSD candidate from QI-inferred neutral mass", unique_neutral_candidate,
  annotation_evidence == "Putative LMSD candidate from m/z-adduct search", unique_mz_candidate,
  default = NA_character_
)]
audit[selected_putative_name == "Coenzyme Q9", selected_putative_name := "CoQ9"]
audit[, display_label := fcase(
  annotation_evidence %in% c(
    "Putative QI annotation",
    "Putative LMSD candidate from QI-inferred neutral mass",
    "Putative LMSD candidate from m/z-adduct search"
  ),
  paste0(met_id, " (putative ", selected_putative_name, ")"),
  default = met_id
)]

audit[, broad_class := fcase(
  grepl("^LPC", accepted_compound_id), "Lysophosphatidylcholine",
  grepl("^PC", accepted_compound_id), "Phosphatidylcholine",
  grepl("^PE", accepted_compound_id), "Phosphatidylethanolamine",
  grepl("^PG", accepted_compound_id), "Phosphatidylglycerol",
  grepl("^SM", accepted_compound_id), "Sphingomyelin",
  grepl("^FFA", accepted_compound_id), "Free fatty acid",
  grepl("CoQ|Coenzyme Q|Ubiquinone", selected_putative_name, ignore.case = TRUE), "Ubiquinone",
  !is.na(unique_candidate_broad_class), unique_candidate_broad_class,
  default = "Unassigned"
)]

audit[, identification_note := fcase(
  annotation_evidence == "Putative QI annotation",
  "QI-reported compound annotation; putative because no authentic-standard confirmation is available.",
  annotation_evidence == "Putative LMSD candidate from QI-inferred neutral mass",
  "QI-inferred neutral mass returned one LIPID MAPS LMSD name within 5 ppm; this mass/adduct-based assignment is putative, not an independently measured neutral species or confirmed identity.",
  annotation_evidence == "Putative LMSD candidate from m/z-adduct search",
  "One LIPID MAPS LMSD name was returned across prespecified positive-mode adduct searches within 5 ppm; the assignment is putative and depends on the assumed adduct.",
  annotation_evidence == "Ambiguous LMSD candidates",
  "Multiple exact-mass candidates were returned; no compound identity was selected.",
  annotation_evidence == "Excluded isotope-labelled LMSD match",
  "Only an isotope-labelled database entry matched; it was not selected because this was not a tracer-labelled experiment.",
  default = "No QI annotation or LIPID MAPS LMSD exact-mass candidate within 5 ppm."
)]

audit_output <- audit[, .(
  met_id, feature, neutral_mass, mz, rt_min,
  accepted_compound_id, accepted_description, adducts, formula,
  score, fragmentation_score,
  qi_concise_name, annotation_evidence, selected_putative_name,
  display_label, broad_class, lmsd_candidate_count,
  unique_candidate_names, eligible_candidate_count, candidate_names,
  eligible_candidate_names, candidate_categories,
  candidate_main_classes, candidate_sub_classes, candidate_broad_classes,
  has_neutral_mass_hit, neutral_mass_unique_names,
  unique_neutral_candidate, unique_mz_candidate,
  best_absolute_delta_ppm, identification_note
)]

annotation_summary <- audit_output[, .N, by = .(annotation_evidence, broad_class)]
setorder(annotation_summary, annotation_evidence, broad_class)

fwrite(candidates, file.path(table_directory, "03_all_LMSD_mass_candidates.csv"), bom = TRUE)
fwrite(audit_output, file.path(table_directory, "03_lipid_feature_annotation_audit.csv"), bom = TRUE)
fwrite(annotation_summary, file.path(table_directory, "03_lipid_annotation_summary.csv"), bom = TRUE)

cat("Code 3 completed successfully.\n")
cat("Features audited:", nrow(audit_output), "\n")
cat("LMSD candidate rows:", nrow(candidates), "\n")
print(annotation_summary)
