# C2C12 lipidomics analysis

Run 00_run_all.R from RStudio or Code availability. Both directories contain
identical R scripts. Codes 01-06 reproduce the main outputs, including the
requested original passage-model figure. Set run_exploratory <- TRUE in Code 00
to also execute Codes 07-08.

Private input: ../Lipidomics_normalized_abundance_positive.csv relative to this
analysis directory; MD5 0c1aa6b5304897d2178794ee9bc62a31. The input is not published.

## Requested original-model figure

Figure/Original_PLSDA_day_horizontal_passage_vertical.png and .svg are generated
by 06_original_passage_PLSDA.R and copied to Critical. Horizontal is negative
original component 2 (44.3% X variance); vertical is original component 1 (8.2%).
This is the original passage-supervised model with a display-only axis/sign
change. The axes retain their actual PLS component identities.

## Workflow

1. 01_load_data_and_PCA.R: normalized-abundance import, averaging three technical
   injections per biological sample, QC/detection filtering, half-minimum
   imputation, log2 transform and Pareto scaling. Global models use 485 features
   and 45 samples (3 passages x 5 Days x 3 replicates).
2. 02_post_differentiation_limma.R: passage contrasts excluding D0. The 60 final
   features require BH-FDR <0.05 for the equal-weight average post-differentiation
   contrast, absolute average change >=1.5-fold, and the same direction in all
   eight Middle/Late-versus-Early contrasts.
3. 03_systematic_lipid_annotation.R: QI accepted annotations and LIPID MAPS
   mass searches at 5 ppm, reported as putative annotations.
4. 04_global_heatmap.R: 60-feature heatmap with MET IDs, putative labels and
   lipid-class annotation.
5. 05_selected_PCA_and_day_PLSDA.R: selected-feature PCA and Day-supervised
   PLS-DA (D0/D1/D3/D6/D9). The PLS-DA uses all 485 features. Component 1 is
   48.5769% and component 2 is 10.6444% X variance; outputs 05_day_PLSDA.
6. 06_original_passage_PLSDA.R: original passage-supervised PLS-DA using the
   same 485 features and 45 samples; requested figure above and result tables
   named 06_original_passage_PLSDA_*.
7. 07_exploratory_PCA_axes.R: descriptive review of PCs 1-10 and PC1-vs-PC2 /
   PC1-vs-PC5 plots. The PC5 view was chosen after inspecting passage associations,
   not as an independent validation. Outputs: Exploratory/PCA.
8. 08_exploratory_ASCA_scores.R: balanced Day + passage + interaction
   decomposition and ASCA-guided factor-focused illustration. Outputs:
   Exploratory/ASCA. This is not global PCA or ordinary PLS-DA.

## Validation and axis interpretation

Each PLS-DA target has separately computed 100 repetitions of balanced 3-fold
cross-validation and 999 fixed-fold permutations. Imputation and Pareto scaling
are fitted in training folds; the 485-feature QC/detection filter is fixed from
the complete input. Prediction uses mixOmics max.dist. Day labels are permuted
within passage in Code 05; passage labels within Day in Code 06.

Code 06 reproduces all original saved coordinates, all CV rows and all
permutation rows. The original passage model's two-component mean balanced CV
accuracy is 0.652222; a separated Early group is not perfect three-class
prediction. Both PLS-DA models retain their natural extraction order, regardless
of display axes. Ellipses are descriptive 95% t ellipses.

The true global PCA has PC1 48.6045%, PC2 10.7194%; the separate selected-feature
PCA has PC1 53.5461%, PC2 14.2610%. PLS percentages are not PCA percentages.
See docs/MULTIVARIATE_METHODS.md for the exploratory ASCA coordinate definitions
and variance denominators.

## Directories and requirements

- RStudio: executable scripts and RStudio project.
- Code availability: identical copies of all nine R files (00-08).
- Figure: primary PNG/PDF outputs and requested original-model PNG/SVG.
- Critical: selected main figures and key tables.
- Tables: primary result tables.
- Exploratory: explicitly exploratory diagnostics and score illustrations.
- docs: method and synchronization notes.

Local Archive discussion/history directories are excluded from GitHub.
Earlier public versions remain in Git history. Raw abundance files and fitted
RDS objects containing individual abundance values are not published.

Packages: data.table, ggplot2, limma, mixOmics, curl, ComplexHeatmap, circlize.
Code 03 uses LIPID MAPS web services; future external search results may differ
from saved annotation tables. All input-dependent computations require the
private input.

Passage labels: P11 = Early, P22 = Middle, P33 = Late. QI/LIPID MAPS assignments
are putative: the export does not establish the original QI search database,
authentic-standard retention matching or MS/MS confirmation.
