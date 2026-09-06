# C2C12 lipidomics analysis

Open `RStudio/C2C12_Lipidomics_analysis.Rproj` and run `00_run_all.R`, or run Codes 01–05 in order. Passage groups are reported as Early Passage, Middle Passage and Late Passage throughout the processed data, statistical contrasts and output tables.

## Analysis workflow

1. `01_load_data_and_PCA.R` loads Progenesis QI normalized abundance, applies the prespecified detection filter and performs global PCA with and without pooled QC injections.
2. `02_post_differentiation_limma.R` fits limma models to log2 normalized abundance. Final features require BH-FDR below 0.05 for the equal-weight average post-differentiation contrast, an absolute average change of at least 1.5-fold, and the same direction in all eight Middle/Late-versus-Early contrasts. The joint F-test is retained as a global screening result.
3. `03_systematic_lipid_annotation.R` retains QI Accepted annotations and searches LIPID MAPS LMSD at 5 ppm using QI-inferred neutral mass or observed m/z with specified positive-mode adducts. The search is mass-based and does not confirm compound identity.
4. `04_global_heatmap.R` draws the 60-feature global heatmap. Rows retain MET identifiers, with putative compound names shown in parentheses where assigned, and a putative lipid-class annotation. Detailed annotation evidence remains in the tables.
5. `05_selected_PCA_and_three_group_PLSDA.R` performs selected-feature PCA and three-group PLS-DA with repeated cross-validation and day-stratified label permutation.

## Directories

- `RStudio`: working analysis scripts and RStudio project.
- `Code availability`: manuscript-ready copies of the same scripts.
- `Figure`: complete PNG and PDF outputs.
- `Tables`: complete CSV results and the annotation workbook.
- `Critical`: main figures and key tables for review.

All figures use enlarged bold axis, legend and annotation text. Dense heatmap sample labels use the abbreviated forms `Early P.`, `Middle P.` and `Late P.`.

## Annotation reporting

QI Accepted Compound IDs are reported as putative QI annotations. The original export does not identify the QI search database, so HMDB is not claimed. Unique LIPID MAPS exact-mass candidates are also putative. Ambiguous and unmatched results remain MET identifiers. No authentic-standard, retention-time or MS/MS confirmation is available in this dataset.
