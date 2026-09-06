# C2C12 lipidomics analysis

This directory is the GitHub-ready copy of the finalized lipidomics workflow.

- Input expected locally: `Lipidomics_normalized_abundance_positive.csv`
- Analysis: `20260905_Final_Analysis/`
- Reproducible entry point: `20260905_Final_Analysis/RStudio/00_run_all.R`
- Submission scripts: `20260905_Final_Analysis/Code availability/`
- Figures: `20260905_Final_Analysis/Figure/`
- Tables: `20260905_Final_Analysis/Tables/`
- Key outputs: `20260905_Final_Analysis/Critical/`

Open `20260905_Final_Analysis/RStudio/C2C12_Lipidomics_analysis.Rproj` in RStudio and run `00_run_all.R`.

The normalized abundance input is not included in the public code release because it is an unpublished in-house dataset. Place the file beside this README before rerunning the workflow.

Feature identities are reported as putative unless supported by an authentic standard and diagnostic MS/MS evidence. The detailed provenance and confidence fields are retained in `Tables/03_lipid_annotation_tables.xlsx`.
