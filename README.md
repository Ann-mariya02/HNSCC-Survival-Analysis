# HNSCC Survival Analysis Using SEER Data

An R-based retrospective cohort analysis comparing **surgery alone** with **surgery plus postoperative radiotherapy (PORT)** among patients with head and neck squamous cell carcinoma (HNSCC). The workflow demonstrates data preparation, missing-data handling, propensity score matching, survival modelling and model diagnostics.

> **Research-code portfolio project.** This repository contains analysis code and documentation only. It does not contain the SEER dataset or patient-level records. The analysis has not been independently executed or statistically validated as part of this repository packaging.

## Project workflow

1. Inspect and prepare the SEER extract; harmonise TNM-stage variables across coding periods.
2. Construct the study cohort and define the treatment groups.
3. Assess missingness and impute T stage using `mice` with multinomial logistic regression (`polyreg`; five imputations).
4. Perform within-imputation nearest-neighbour propensity score matching, using a 0.2 caliper and exact matching on age stratum, as implemented in the script.
5. Assess covariate balance using standardised mean differences and Love plots.
6. Fit Cox proportional hazards models and pool estimates using the Rubin-rule functions in the script.
7. Assess proportional-hazards assumptions using Schoenfeld residuals and perform 60-month restricted mean survival time (RMST) analysis.

See [`documentation/methodology.md`](documentation/methodology.md) for a concise description of the implemented workflow.

## Repository structure

```text
HNSCC-Survival-Analysis/
├── README.md
├── .gitignore
├── packages.md
├── setup/
│   └── install_packages.R
├── R/
│   └── 01_hnscc_seer_analysis.R
└── documentation/
    └── methodology.md
```

The `data/` and `figures/` directories are local working directories and are intentionally not populated in this repository.

## Requirements and setup

- R (use an up-to-date R installation compatible with the listed packages).
- RStudio is recommended but not required.
- The packages listed in [`packages.md`](packages.md).

From the repository root, run the setup script once:

```r
source("setup/install_packages.R")
```

## Data access and execution

1. Obtain the SEER extract through the appropriate authorised process and follow all applicable SEER data-use terms and institutional requirements.
2. Create a local folder named `data` in the repository root.
3. Place the authorised Excel file in that folder with the filename `HNSCC_SEER_Dataset.xlsx`.
4. Open `R/01_hnscc_seer_analysis.R` and run it from the repository root, in sequence.

The script reads the file from `data/HNSCC_SEER_Dataset.xlsx`. The `data/` directory and common data-file formats are excluded by `.gitignore`. **Never commit the SEER extract, patient-level records, direct identifiers, restricted extracts, credentials or confidential RCC materials.** The script refers to a `Patient ID` field for a uniqueness check; no patient-level ID values should be committed. Do not publish derived tables or figures unless their release is permitted and disclosure risk has been reviewed.

## Reproducibility notes

- The imputation call specifies five imputations, ten iterations and a seed of `123`.
- The script prints `sessionInfo()` at the end so the R version and installed package versions can be recorded.
- Package versions are not pinned in this repository. For a more tightly controlled environment, consider adopting `renv` and committing an approved lockfile after confirming the environment used for the analysis.
- Reproducibility also depends on the authorised source-data version and local execution environment.

## Scope and limitations

This repository packaging separates package installation from the analysis script and documents the existing workflow. It does **not** change the underlying statistical analysis or establish that the modelling choices are appropriate for every use. The script was not run against the source dataset during packaging. Review cohort definitions, assumptions, diagnostics and outputs with the research supervisor before drawing substantive conclusions.

This is not clinical advice and is not a validated clinical decision-support tool.


