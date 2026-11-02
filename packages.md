# R package requirements

The analysis script loads the packages listed below. `optmatch` is also included because it is installed in the original workflow.

| Package | Role in the workflow |
|---|---|
| `MatchThem` | Propensity score matching with multiply imputed data |
| `mice` | Multiple imputation |
| `gtsummary` | Summary tables |
| `survival` | Survival objects and Cox regression |
| `survminer` | Survival-analysis visualisation |
| `ggplot2` | Plotting |
| `dplyr` | Data manipulation |
| `flextable` | Formatted tables |
| `MatchIt` | Matching methods used by the workflow |
| `cobalt` | Covariate-balance assessment and Love plots |
| `naniar` | Missing-data exploration |
| `nnet` | Multinomial modelling support used by `mice`'s `polyreg` method |
| `readxl` | Reading the SEER Excel extract |
| `survRM2` | Restricted mean survival time analysis |
| `optmatch` | Matching dependency included in the original setup |

## Installation

From the repository root in R or RStudio, run:

```r
source("setup/install_packages.R")
```

The setup script installs only packages that are not already present. Package versions are not pinned, so installations made at different times may differ.

## Recording the software environment

The main analysis script prints `sessionInfo()` at the end of execution. To capture it separately, run:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

Keep the resulting file with your analysis records. Only commit it to this repository if sharing the software environment is appropriate and permitted.
