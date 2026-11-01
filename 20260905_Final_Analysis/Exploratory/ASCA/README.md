# Day / passage factor-focused score illustration (2026-11-02)

This is an exploratory illustration requested to put Day-related variation on
the horizontal axis and passage-related variation on the vertical axis.
It is based on a balanced two-factor ASCA decomposition, not a global PCA or
a standard two-component PLS-DA fit.

Input:
- Same normalized-abundance CSV; MD5 0c1aa6b5304897d2178794ee9bc62a31.
- Same 485 QC-retained features, 45 biological samples, 15 conditions, 3 samples
  per condition. Technical-injection averaging, half-minimum imputation, log2,
  and Pareto scaling are unchanged.
- Script runs the import/preprocessing portion of production Code 01, stopping
  before production plots/tables are written.

Model:
X_centered = D + P + I + E.
D is the Day main-effect matrix, P the passage main-effect matrix, I their
interaction matrix, and E the within-condition residual matrix.
PCA is fitted separately to D and P. The first loading vector of each effect
is used. Sign conventions only orient later days rightward and Late above Early.
These choices do not guarantee that Middle will lie between Early and Late.

Joint display (custom, ASCA-guided):
horizontal_i = (X_centered - P)_i times the first Day-effect loading;
vertical_i = (X_centered - D)_i times the first Passage-effect loading.
Thus each coordinate uses a different adjusted version of the sample.
Interaction and within-condition variability remain on both coordinates.
This is NOT a common orthogonal projection of the original X matrix onto two
global principal components, and plot distances are not global PCA distances.
No sample was manually moved, removed, jittered, or averaged into a plotted point.
The circles are descriptive 95% t ellipses, not classifier boundaries.

Interpretation:
The axes are learned using the known Day and passage labels. Better alignment
with those factors is therefore part of the modeling, not independent evidence
of predictive accuracy or a confirmation of the RNA-seq PCA pattern.
The display still shows Middle/Late overlap. All 45 sample coordinates are
provided. A standard, separate ASCA effect-plus-residual diagnostic is also saved.
No permutation tests, classification validation, or p-values were generated.
The existing PCA/PLS-DA formal outputs have not been replaced.

Variance:
Total predictor SS: Day 59.248695%, Passage 4.873711%, interaction 12.209179%,
within-condition variation 23.668415%.
First Day-effect PC captures 78.60007% of the DAY EFFECT, corresponding to
46.569515% of total X SS in that fitted main-effect contribution.
First Passage-effect PC captures 69.61190% of the PASSAGE EFFECT, corresponding
to 3.392683% of total X SS in that fitted main-effect contribution.
Percentages within different effect matrices must not be treated as percentages
of the same denominator or as the variances of the residual-augmented points.

Verification:
Mean-based decomposition matches an independent sum-to-zero-coded factorial
OLS fit for all 485 features to <1e-9; components exactly reconstruct centered X;
factor sums of squares sum to total SS. Raw hash, all sample IDs, dimensions and
the balanced condition counts are checked by the executable script.

Sources:
Smilde et al. 2005, ASCA:
https://academic.oup.com/bioinformatics/article/21/13/3043/197836
Timmerman et al. 2015, effect decomposition and variance denominators:
https://pmc.ncbi.nlm.nih.gov/articles/PMC4559107/
The combined opposite-main-effect-adjusted score display is explicitly defined
above; the references support ASCA decomposition, not a claim that this custom
joint display is an ordinary PCA score plot.
