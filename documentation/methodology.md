# Methodology overview

This note summarises the workflow implemented in the supplied R script. It does not independently verify the validity of each analytic choice.

1. **Data inspection and preparation:** inspect dimensions, variable names and missingness; harmonise T, N and M stage fields across SEER/AJCC coding periods; group primary sites; rename and arrange variables.
2. **Cohort construction:** restrict to M0 records, define surgery-only and surgery-plus-PORT groups, apply survival-time and staging filters, and derive age strata and event status.
3. **Missing-data work:** examine missing T/N stage patterns and impute T stage with `mice` using `polyreg` across five imputations.
4. **Confounder assessment and matching:** compare full and reduced Cox models as coded in the script, then perform within-imputation nearest-neighbour propensity score matching with a 0.2 standardised caliper and exact matching on age stratum.
5. **Outcome analysis:** fit univariate and multivariable Cox models, including a treatment-by-age-stratum interaction model and a main-effects model; pool estimates using the custom Rubin-rule functions in the script.
6. **Diagnostics and complementary analysis:** assess proportional-hazards assumptions using Schoenfeld residuals and `cox.zph`; estimate 60-month RMST differences and pool the estimates.

Refer to the R script for the exact variable definitions, filters, formulas and implementation details.
