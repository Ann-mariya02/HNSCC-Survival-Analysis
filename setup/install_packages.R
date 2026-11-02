# Install packages required by the HNSCC analysis workflow.
# Run this script once before running R/01_hnscc_seer_analysis.R.

required_packages <- c(
  "MatchThem",
  "mice",
  "gtsummary",
  "survival",
  "survminer",
  "ggplot2",
  "dplyr",
  "flextable",
  "MatchIt",
  "cobalt",
  "naniar",
  "nnet",
  "readxl",
  "survRM2",
  "optmatch"
)

installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
} else {
  message("All required packages are already installed.")
}

# Display versions after installation for reproducibility.
package_versions <- installed.packages()
print(package_versions[intersect(required_packages, rownames(package_versions)),
                       c("Package", "Version"), drop = FALSE])
