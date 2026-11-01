# Multivariate methods and display definitions

## Shared inputs and separate models

All global models use the same 485 retained features and 45 biological samples,
technical-injection averaging, half-minimum imputation, log2 and Pareto scaling.
The 60 selected features apply only to the selected PCA and heatmap.

Ordinary PCA: PC1/PC2 = 48.604511%/10.719382%. PC5 = 3.987364% and shows partial
passage structure. Its displayed PC1-vs-PC5 view was chosen after looking at the
first ten components and is exploratory, not independent validation.

Day-supervised PLS-DA: 48.576855%/10.644389% X variance.
Original passage-supervised PLS-DA: 8.216946%/44.253320%.
The requested original figure shows negative original component 2 horizontally
and original component 1 vertically. No component is relabeled as PC1 or PC2.

Both targets have 100 repetitions of balanced 3-fold CV and 999 permutations,
with train-fold preprocessing and fixed complete-input detection filtering.
Day labels are permuted within passage; passage labels are permuted within Day.
The original passage model's two-component mean balanced accuracy is 0.652222,
including Middle sensitivity 0.372; visual separation is not perfect prediction.

## ASCA-guided joint score illustration

Code 08 decomposes centered X into D + P + I + E: Day main effect, passage main
effect, interaction, and within-condition residual. PCA is fitted separately to
D and P; the first loading of each effect defines its direction.

The custom combined view is explicitly:
horizontal = (X_centered - P) times the first Day-effect loading;
vertical = (X_centered - D) times the first passage-effect loading.

Each coordinate uses a different adjusted matrix. Interaction and within-condition
variation remain. This is not a common orthogonal projection onto two global PCs;
plot distances are not global PCA distances. Known factor labels are used to
learn directions, so alignment with those labels is part of the modeling.
No permutation test, predictive validation or p-value was generated for this
illustrative view. Middle/Late overlap remains. The 95% ellipses are descriptive.

Total predictor sums of squares:
Day 59.248695%, passage 4.873711%, interaction 12.209179%, within-condition 23.668415%.
The first Day-effect component captures 78.60007% of the Day effect; the first
passage-effect component captures 69.61190% of the passage effect. These have
different denominators. Their fitted main-effect contributions are 46.569515%
and 3.392683% of total X SS, respectively; these are not the total variances of
the residual-augmented plotted points.

The decomposition is independently checked against a sum-coded factorial OLS
fit for all features, and components reconstruct centered X to numerical tolerance.
Fitted matrices/RDS files with individual abundance values are not published.

Sources:
- ASCA: https://doi.org/10.1093/bioinformatics/bti476
- Decomposition and variance denominators:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4559107/
- PCA/PLS-DA ordering: https://mixomics.org/faq/
