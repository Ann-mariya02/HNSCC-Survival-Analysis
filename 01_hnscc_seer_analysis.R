# HNSCC survival analysis using SEER data
# Before running this script, execute setup/install_packages.R once.
# Package installation is kept separate from the analysis workflow.
library(MatchThem)
library(mice)
library(gtsummary)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(flextable)
library(MatchIt)
library(cobalt)
library(naniar)
library(nnet)
library(readxl)
library(survRM2)

# Read the Excel file
df <- read_excel("data/HNSCC_SEER_Dataset.xlsx")  
# View the data
head(df)        
names(df)       
dim(df)     
# Missing values per column
colSums(is.na(df))
# Percentage missing per column
colMeans(is.na(df)) * 100

#Merging all T columns into one "T stage",N columns into "N stage"..
all_names <- names(df)
t_cols <- grep("^Derived EOD 2018 T|^Derived AJCC T|^Derived SEER Combined T|^T value", all_names, value = TRUE)
n_cols <- grep("^Derived EOD 2018 N|^Derived AJCC N|^Derived SEER Combined N|^N value", all_names, value = TRUE)
m_cols <- grep("^Derived EOD 2018 M|^Derived AJCC M|^Derived SEER Combined M|^M value", all_names, value = TRUE)
t_cols
n_cols
m_cols
cat("\n--- Derived EOD 2018 T (2018+) ---\n")
print(unique(df[[t_cols[1]]]))
cat("\n--- Derived AJCC T, 7th ed (2010-2015) ---\n")
print(unique(df[[t_cols[2]]]))
cat("\n--- Derived SEER Combined T (2016-2017) ---\n")
print(unique(df[[t_cols[3]]]))
cat("\n--- Derived AJCC T, 6th ed (2004-2015) ---\n")
print(unique(df[[t_cols[4]]]))
cat("\n--- T value - based on AJCC 3rd (1988-2003) ---\n")
print(unique(df[[t_cols[5]]]))
cat("\n--- Derived EOD 2018 N (2018+) ---\n")
print(unique(df[[n_cols[1]]]))
cat("\n--- Derived AJCC N, 7th ed (2010-2015) ---\n")
print(unique(df[[n_cols[2]]]))
cat("\n--- Derived SEER Combined N (2016-2017) ---\n")
print(unique(df[[n_cols[3]]]))
cat("\n--- Derived AJCC N, 6th ed (2004-2015) ---\n")
print(unique(df[[n_cols[4]]]))
cat("\n--- N value - based on AJCC 3rd (1988-2003) ---\n")
print(unique(df[[n_cols[5]]]))
cat("\n--- Derived EOD 2018 M (2018+) ---\n")
print(unique(df[[m_cols[1]]]))
cat("\n--- Derived AJCC M, 7th ed (2010-2015) ---\n")
print(unique(df[[m_cols[2]]]))
cat("\n--- Derived SEER Combined M (2016-2017) ---\n")
print(unique(df[[m_cols[3]]]))
cat("\n--- Derived AJCC M, 6th ed (2004-2015) ---\n")
print(unique(df[[m_cols[4]]]))
cat("\n--- M value - based on AJCC 3rd (1988-2003) ---\n")
print(unique(df[[m_cols[5]]]))

#Derived EOD 2018 T (2018+)[T1-T1,T1a,T1b,..]
df <- df %>%
  mutate(
    `Derived EOD 2018 T (2018+)` = case_when(
      `Derived EOD 2018 T (2018+)` %in% c("T1", "T1a", "T1b") ~ "T1",
      `Derived EOD 2018 T (2018+)` %in% c("T4", "T4a", "T4b") ~ "T4",
      TRUE ~ `Derived EOD 2018 T (2018+)`
    )
  )
unique(df$`Derived EOD 2018 T (2018+)`)
#Derived AJCC T, 7th ed (2010-2015)[T1-T1,T1a,T1b,..]
df <- df %>%
  mutate(
    `Derived AJCC T, 7th ed (2010-2015)` = case_when(
      `Derived AJCC T, 7th ed (2010-2015)` %in% c("T1", "T1a", "T1b", "T1NOS") ~ "T1",
      `Derived AJCC T, 7th ed (2010-2015)` %in% c("T4", "T4a", "T4b", "T4NOS") ~ "T4",
      TRUE ~ `Derived AJCC T, 7th ed (2010-2015)`
    )
  )
unique(df$`Derived AJCC T, 7th ed (2010-2015)`)
#Derived SEER Combined T (2016-2017)[T1-T1,T1a,T1b,..]
df <- df %>%
  mutate(
    `Derived SEER Combined T (2016-2017)` = case_when(
      `Derived SEER Combined T (2016-2017)` %in% c("c0") ~ "T0",
      `Derived SEER Combined T (2016-2017)` %in% c("c1", "c1A", "c1B", "p1", "p1A", "p1B") ~ "T1",
      `Derived SEER Combined T (2016-2017)` %in% c("c2", "p2") ~ "T2",
      `Derived SEER Combined T (2016-2017)` %in% c("c3", "p3") ~ "T3",
      `Derived SEER Combined T (2016-2017)` %in% c("c4A", "c4B", "p4A", "p4B") ~ "T4",
      `Derived SEER Combined T (2016-2017)` %in% c("cX", "pX") ~ "TX",
      `Derived SEER Combined T (2016-2017)` %in% c("Blank(s)", "Not applicable") ~ NA_character_,
      TRUE ~ `Derived SEER Combined T (2016-2017)`
    )
  )
unique(df$`Derived SEER Combined T (2016-2017)`)
#Derived AJCC T, 6th ed (2004-2015)[T1-T1,T1a,T1b,..]
df <- df %>%
  mutate(
    `Derived AJCC T, 6th ed (2004-2015)` = case_when(
      `Derived AJCC T, 6th ed (2004-2015)` %in% c("T1", "T1a", "T1b", "T1NOS") ~ "T1",
      `Derived AJCC T, 6th ed (2004-2015)` %in% c("T4", "T4a", "T4b", "T4NOS") ~ "T4",
      TRUE ~ `Derived AJCC T, 6th ed (2004-2015)`
    )
  )
unique(df$`Derived AJCC T, 6th ed (2004-2015)`)
#Derived EOD 2018 N (2018+)[N2,N2a,N2b,N2c into N2..]
df <- df %>%
  mutate(
    `Derived EOD 2018 N (2018+)` = case_when(
      `Derived EOD 2018 N (2018+)` %in% c("N2", "N2a", "N2b", "N2c") ~ "N2",
      `Derived EOD 2018 N (2018+)` %in% c("N3", "N3a", "N3b") ~ "N3",
      TRUE ~ `Derived EOD 2018 N (2018+)`
    )
  )
unique(df$`Derived EOD 2018 N (2018+)`)
#Derived AJCC N, 7th ed (2010-2015)[N2,N2a,N2b,N2c into N2..]
df <- df %>%
  mutate(
    `Derived AJCC N, 7th ed (2010-2015)` = case_when(
      `Derived AJCC N, 7th ed (2010-2015)` %in% c("N2", "N2a", "N2b", "N2c", "N2NOS") ~ "N2",
      TRUE ~ `Derived AJCC N, 7th ed (2010-2015)`
    )
  )
unique(df$`Derived AJCC N, 7th ed (2010-2015)`)
#Derived SEER Combined N (2016-2017)[N2,N2a,N2b,N2c into N2..]
df <- df %>%
  mutate(
    `Derived SEER Combined N (2016-2017)` = case_when(
      `Derived SEER Combined N (2016-2017)` %in% c("c0", "p0") ~ "N0",
      `Derived SEER Combined N (2016-2017)` %in% c("c1", "p1") ~ "N1",
      `Derived SEER Combined N (2016-2017)` %in% c("c2", "c2A", "c2B", "c2C", "p2", "p2A", "p2B", "p2C") ~ "N2",
      `Derived SEER Combined N (2016-2017)` %in% c("c3", "p3") ~ "N3",
      `Derived SEER Combined N (2016-2017)` %in% c("cX", "pX") ~ "NX",
      `Derived SEER Combined N (2016-2017)` %in% c("Blank(s)", "Not applicable") ~ NA_character_,
      TRUE ~ `Derived SEER Combined N (2016-2017)`
    )
  )
unique(df$`Derived SEER Combined N (2016-2017)`)
#Derived AJCC N, 6th ed (2004-2015)[N2,N2a,N2b,N2c into N2..]
df <- df %>%
  mutate(
    `Derived AJCC N, 6th ed (2004-2015)` = case_when(
      `Derived AJCC N, 6th ed (2004-2015)` %in% c("N2", "N2a", "N2b", "N2c", "N2NOS") ~ "N2",
      TRUE ~ `Derived AJCC N, 6th ed (2004-2015)`
    )
  )
unique(df$`Derived AJCC N, 6th ed (2004-2015)`)
#Derived SEER Combined M (2016-2017)[C0 into M0..]
df <- df %>%
  mutate(
    `Derived SEER Combined M (2016-2017)` = case_when(
      `Derived SEER Combined M (2016-2017)` %in% c("c0") ~ "M0",
      `Derived SEER Combined M (2016-2017)` %in% c("c1", "p1") ~ "M1",
      `Derived SEER Combined M (2016-2017)` %in% c("Blank(s)", "Not applicable") ~ NA_character_,
      TRUE ~ `Derived SEER Combined M (2016-2017)`
    )
  )
unique(df$`Derived SEER Combined M (2016-2017)`)
           
all_tnm_cols <- c(t_cols, n_cols, m_cols)
# Step 1: convert all missing-like strings to real NA
df <- df %>%
  mutate(across(all_of(all_tnm_cols), ~na_if(., "Blank(s)"))) %>%
  mutate(across(all_of(all_tnm_cols), ~na_if(., "NA"))) %>%
  mutate(across(all_of(all_tnm_cols), ~na_if(., "Not applicable"))) %>%
  mutate(across(all_of(all_tnm_cols), ~na_if(., "88")))
# Step 2: coalesce each group into a single column 
df <- df %>%
  mutate(
    `T stage` = coalesce(!!!syms(t_cols)),
    `N stage` = coalesce(!!!syms(n_cols)),
    `M stage` = coalesce(!!!syms(m_cols))
  )
# Step 3: drop the original 15 columns, keep only the 3 final ones
df <- df %>%
  select(-all_of(all_tnm_cols))
table(df$`T stage`, useNA = "always")
table(df$`N stage`, useNA = "always")
table(df$`M stage`, useNA = "always")

#Based on Exclusion criteria- No Distant Metastasis- So Include only M0
df <- df %>%
  filter(`M stage` == "M0")
table(df$`M stage`, useNA = "always")
dim(df)
table(df$`T stage`, useNA = "always")
table(df$`N stage`, useNA = "always")
table(df$`M stage`, useNA = "always")

#Primary site converting to Oral Cavity, Oropharynx,Hypopharynx,Larynx
table(df$`Primary Site - labeled`)
oral_cavity_codes <- c("C02.0","C02.1","C02.2","C02.3","C02.8","C02.9",
                       "C03.0","C03.1","C03.9",
                       "C04.0","C04.1","C04.8","C04.9",
                       "C05.0","C05.8","C05.9",
                       "C06.0","C06.1","C06.2","C06.8","C06.9")

oropharynx_codes <- c("C01.9","C02.4","C05.1","C05.2",
                      "C09.0","C09.1","C09.8","C09.9",
                      "C10.0","C10.2","C10.3","C10.8","C10.9")

hypopharynx_codes <- c("C12.9","C13.0","C13.1","C13.2","C13.8","C13.9")

larynx_codes <- c("C32.0","C32.1","C32.2","C32.3","C32.8","C32.9")

# Create the new grouping column named "Primary sites"
df <- df %>%
  mutate(
    `Primary sites` = case_when(
      sub("-.*", "", `Primary Site - labeled`) %in% oral_cavity_codes ~ "Oral cavity",
      sub("-.*", "", `Primary Site - labeled`) %in% oropharynx_codes ~ "Oropharynx",
      sub("-.*", "", `Primary Site - labeled`) %in% hypopharynx_codes ~ "Hypopharynx",
      sub("-.*", "", `Primary Site - labeled`) %in% larynx_codes ~ "Larynx",
      TRUE ~ NA_character_
    )
  )
table(df$`Primary sites`, useNA = "always")
 

# Drop the original "Primary Site - labeled" column
df <- df %>%
  select(-`Primary Site - labeled`)
table(df$`Primary sites`, useNA = "always")

# Check if Patient ID is unique (no duplicates)
length(df$`Patient ID`) == length(unique(df$`Patient ID`))

#Renaming the columns 
df <- df %>%
  rename(
    `Age category` = `Age recode with <1 year olds`,
    `Race` = `Race recode (W, B, AI, API)`,
    `Surgery` = `RX Summ--Surg Prim Site (1998+)`,
    `Radiation Sequence with surgery` = `RX Summ--Surg/Rad Seq`,
    `Radiation` = `Radiation recode`,
    `Vital status` = `Vital status recode (study cutoff used)`,
    `Year of followup` = `Year of follow-up recode`,
    `Year of death` = `Year of death recode`,
    `Age` = `Age recode with single ages and 85+`
  )

# Rearrange columns in a logical order
df <- df %>%
  select(
    `Patient ID`,
    Sex,
    Age,
    `Age category`,
    Race,
    `Year of diagnosis`,
    `Primary sites`,
    `T stage`,
    `N stage`,
    `M stage`,
    Surgery,
    `Radiation Sequence with surgery`,
    Radiation,
    `Survival months`,
    `Survival months flag`,
    `Vital status`,
    `Year of followup`,
    `Year of death`
  )
for (col in names(df)) {
  cat("\n---", col, "---\n")
  print(unique(df[[col]]))
}
#Based on Inclusion criteria filter only- Surgery+PORT and Surgery only patients
table(df$Surgery, useNA = "always")
#Filtering surgery Column to get only patients who surgery and the procedure is confirmed
df <- df %>%
  filter(!Surgery %in% 90 & !is.na(Surgery))
table(df$Surgery, useNA = "always")
dim(df)
#Converting all the numeric values in the surgery column to Yes
df <- df %>%
  mutate(Surgery = ifelse(!is.na(Surgery), "Yes", NA))
table(df$Surgery, useNA = "always")
#Filtering the Columns Radiation sequence with surgery-No radiation and Radiation after surgery
table(df$Surgery, useNA = "always")
table(df$`Radiation Sequence with surgery`, useNA = "always")
table(df$Radiation, useNA = "always")
df <- df %>%
  filter(`Radiation Sequence with surgery` %in% c(
    "No radiation and/or no surgery; unknown if surgery and/or radiation given",
    "Radiation after surgery"
  ))
table(df$`Radiation Sequence with surgery`, useNA = "always")
dim(df)
#Filtering for Surgery only and Surgery+PORT group as Treatment variable
df <- df %>%
  mutate(
    Treatment = case_when(
      Surgery == "Yes" & 
        `Radiation Sequence with surgery` == "Radiation after surgery" & 
        Radiation == "Beam radiation" ~ "Surgery + PORT",
      
      Surgery == "Yes" & 
        `Radiation Sequence with surgery` == "No radiation and/or no surgery; unknown if surgery and/or radiation given" &
        Radiation == "None/Unknown" ~ "Surgery only",
      
      TRUE ~ NA_character_
    )
  )
# Explicitly set factor levels, "Surgery only" as the reference category
df <- df %>%
  mutate(Treatment = factor(Treatment, levels = c("Surgery only", "Surgery + PORT")))
table(df$Treatment, useNA = "always")
df <- df %>% filter(!is.na(Treatment))
df <- df %>%
  select(-c(Surgery, `Radiation Sequence with surgery`, Radiation))
dim(df)
# Confirm the factor levels survived filtering/column dropping
levels(df$Treatment)

# Filtering for Survival month flag -- keep patients with complete dates (>0 days survival)
# AND patients with incomplete dates where zero follow-up is impossible

table(df$`Survival months flag`, useNA = "always")

df <- df%>%
  filter(`Survival months flag` %in% c(
    "Complete dates are available and there are more than 0 days of survival",
    "Incomplete dates are available and there cannot be zero days of follow-up"
  ))

dim(df)
table(df$`Survival months flag`, useNA = "always")  # verify only these two remain
sum(df$`Survival months` == 0, na.rm = TRUE)
round(100 * mean(df$`Survival months` == 0, na.rm = TRUE), 2)
# Remove patients with 0 survival months (not eligible to receive PORT)
n_before <- nrow(df)
df <- df %>%
  filter(`Survival months` > 0)
n_after <- nrow(df)
cat("Removed", n_before - n_after, "patients with 0 survival months\n")
dim(df)


#Checking the values in T stage and N stage - filtered only T1,T2,T3,T4 and N0,N1,N2,N3
table(df$`T stage`, useNA = "always")
table(df$`N stage`, useNA = "always")
n_before <- nrow(df)
n_before
df <- df %>%
  filter(
    !`T stage` %in% c("Tis", "T0"),
    !is.na(`T stage`),
    !is.na(`N stage`)
  )

n_after <- nrow(df)
cat("Removed", n_before - n_after, "patients (Tis and T0/ NA in T or N stage)\n")
table(df$`T stage`, useNA = "always")
table(df$`N stage`, useNA = "always")
dim(df)
# Total number of patients
n_total <- nrow(df)
# Count of TX and NX
n_TX <- sum(df$`T stage` == "TX", na.rm = TRUE)
n_NX <- sum(df$`N stage` == "NX", na.rm = TRUE)
# Percentage of TX (unknown T stage)
pct_TX <- round(100 * n_TX / n_total, 2)
# Percentage of NX (unknown N stage)
pct_NX <- round(100 * n_NX / n_total, 2)
cat("Total patients:", n_total, "\n")
cat("TX count:", n_TX, "| TX percentage:", pct_TX, "%\n")
cat("NX count:", n_NX, "| NX percentage:", pct_NX, "%\n")
# Cross-tab of T stage vs N stage (counts)
t_n_combinations <- df %>%
  count(`T stage`, `N stage`) %>%
  arrange(`T stage`, `N stage`)
print(t_n_combinations, n = Inf)   
# Counts matrix
t_n_matrix <- table(`T stage` = df$`T stage`, `N stage` = df$`N stage`)
print(t_n_matrix)
#Removing TX-N0 , TX-NX combinations
n_before <- nrow(df)

df <- df %>%
  filter(
    !(`T stage` == "T0" & `N stage` == "NX"),
    !(`T stage` == "TX" & `N stage` == "N0"),
    !(`T stage` == "TX" & `N stage` == "NX")
  )

n_after <- nrow(df)
cat("Removed", n_before - n_after, "patients (TX-N0, TX-NX combinations)\n")

# Cross-tab of T stage vs N stage (counts)
t_n_matrix <- table(`T stage` = df$`T stage`, `N stage` = df$`N stage`)
print(t_n_matrix)

# Create indicator flags before converting to NA
df <- df %>%
  mutate(
    N_stage_missing = (`T stage` %in% c("T1", "T2", "T3", "T4") & `N stage` == "NX"),
    T_stage_missing = (`N stage` %in% c("N1", "N2", "N3") & `T stage` == "TX")
  )

# Check how many rows each flag catches
table(df$N_stage_missing)
table(df$T_stage_missing)

# Cross-tab check
t_n_matrix <- table(`T stage` = df$`T stage`, `N stage` = df$`N stage`, useNA = "always")
print(t_n_matrix)

#Converting age to category as new variable- Age stratum
df <- df %>%
  mutate(Age_numeric = as.numeric(gsub("[^0-9]", "", Age)))
df %>%
  filter(Age_numeric == 85) %>%
  count(Age)  
df <- df %>%
  mutate(Age_stratum = case_when(
    Age_numeric >= 40 & Age_numeric <= 64 ~ "Middle-aged (40-64)",
    Age_numeric >= 65 & Age_numeric <= 74 ~ "Young old (65-74)",
    Age_numeric >= 75 & Age_numeric <= 84 ~ "Older old (75-84)",
    Age_numeric >= 85                     ~ "Oldest old (85+)",
    TRUE ~ NA_character_
  ))
df <- df %>%
  mutate(Age_stratum = factor(Age_stratum,
                              levels = c("Middle-aged (40-64)",
                                         "Young old (65-74)",
                                         "Older old (75-84)",
                                         "Oldest old (85+)")))
table(df$Age_stratum)
#We find missing data in T stage and N stage 
#Checking whether that missingness is due to MCAR/MAR 
#T stage 
#Convert two groups - one with complete data and one with missing data
df_test <- df %>%
  mutate(`T stage` = ifelse(`T stage` == "TX", NA, `T stage`))

group_A <- df_test %>% filter(!is.na(`T stage`))
group_B <- df_test %>% filter(is.na(`T stage`))

nrow(group_A)   # should be your total minus TX count
nrow(group_B)   # should equal your TX count
#Fitting multinomial regression
model_A <- multinom(`T stage` ~ Sex + Age_stratum + Race + `Primary sites`, data = group_A)
summary(model_A)
#Generate predicted probabilities for both groups, using the same fitted model
pred_probs_A <- predict(model_A, newdata = group_A, type = "probs")
pred_probs_B <- predict(model_A, newdata = group_B, type = "probs")
dim(pred_probs_A)   # rows = nrow(group_A), columns = number of T-stage categories
dim(pred_probs_B)   # rows = nrow(group_B), columns = same categories
head(pred_probs_A)
head(pred_probs_B)
#Performing ks test
categories <- colnames(pred_probs_A)

ks_results <- data.frame(Category = categories, D_statistic = NA, p_value = NA)

for (i in seq_along(categories)) {
  cat_name <- categories[i]
  ks_test <- ks.test(pred_probs_A[, cat_name], pred_probs_B[, cat_name])
  ks_results$D_statistic[i] <- ks_test$statistic
  ks_results$p_value[i] <- ks_test$p.value
}

print(ks_results)
ks_results$p_adjusted_BH <- p.adjust(ks_results$p_value, method = "BH")
print(ks_results)
if (all(ks_results$p_adjusted_BH > 0.05)) {
  cat("All categories: no significant distribution difference -> consistent with MCAR\n")
} else {
  cat("At least one category shows a significant difference -> NOT MCAR (consistent with MAR)\n")
  print(ks_results[ks_results$p_adjusted_BH <= 0.05, ])
}
#N stage
# Convert NX to NA for this diagnostic
df_test_N <- df %>%
  mutate(`N stage` = ifelse(`N stage` == "NX", NA, `N stage`))

group_A_N <- df_test_N %>% filter(!is.na(`N stage`))
group_B_N <- df_test_N %>% filter(is.na(`N stage`))

nrow(group_A_N)   # should be 21166+3562+3455+391 = 28,574
nrow(group_B_N)   # should equal your NX count (571)

#  Fit multinomial model on Group A only
model_A_N <- multinom(`N stage` ~ Sex + Age_stratum + Race + `Primary sites`,
                      data = group_A_N)
summary(model_A_N)

# Predict for both groups
pred_probs_A_N <- predict(model_A_N, newdata = group_A_N, type = "probs")
pred_probs_B_N <- predict(model_A_N, newdata = group_B_N, type = "probs")

# KS test per category
categories_N <- colnames(pred_probs_A_N)
ks_results_N <- data.frame(Category = categories_N, D_statistic = NA, p_value = NA)

for (i in seq_along(categories_N)) {
  cat_name <- categories_N[i]
  ks_test <- ks.test(pred_probs_A_N[, cat_name], pred_probs_B_N[, cat_name])
  ks_results_N$D_statistic[i] <- ks_test$statistic
  ks_results_N$p_value[i] <- ks_test$p.value
}

print(ks_results_N)

# BH correction
ks_results_N$p_adjusted_BH <- p.adjust(ks_results_N$p_value, method = "BH")
print(ks_results_N)

# Decision
if (all(ks_results_N$p_adjusted_BH > 0.05)) {
  cat("All categories: no significant distribution difference -> consistent with MCAR\n")
} else {
  cat("At least one category shows a significant difference -> NOT MCAR (consistent with MAR)\n")
  print(ks_results_N[ks_results_N$p_adjusted_BH <= 0.05, ])
}

#Removing missing N stage patients because  regression-based categorical test did not provide evidence against MCAR

n_before <- nrow(df)
cat("Before filtering:", n_before, "rows\n")
df <- df %>%
  filter(!(`N stage` %in% c("NX") | is.na(`N stage`)))
n_after <- nrow(df)
cat("Removed", n_before - n_after, "patients with missing N stage (NX or NA)\n")
table(df$`N stage`, useNA = "always")
dim(df)


#Statistical Analysis


#1 MICE imputation for T stage 
#Multinomial logistic regression (MI-GLM / polyreg) 

# Create the event variable from Vital status
df <- df %>%
  mutate(event = ifelse(`Vital status` == "Dead", 1, 0))
table(df$event, useNA = "always")   # confirm all 0/1, no NA

# Convert TX to genuine NA for imputation
df_impute <- df %>%
  mutate(`T stage` = ifelse(`T stage` == "TX", NA, `T stage`))

# Convert key categorical variables to factors
df_impute <- df_impute %>%
  mutate(
    `T stage` = factor(`T stage`, levels = c("T1", "T2", "T3", "T4")),
    `N stage` = factor(`N stage`, levels = c("N0", "N1", "N2", "N3")),
    Sex = factor(Sex, levels = c("Female", "Male")),                    # explicit: Female = reference
    Race = factor(Race, levels = c("White", "Black", "Asian or Pacific Islander", 
                                   "American Indian/Alaska Native")),   # explicit: White = reference (often preferred as largest group)
    `Primary sites` = factor(`Primary sites`, levels = c("Oral cavity", "Oropharynx",
                                                         "Larynx", "Hypopharynx")),  # explicit: Oral cavity = reference
    Treatment = factor(Treatment, levels = c("Surgery only", "Surgery + PORT"))
  )
levels(df_impute$Treatment)
levels(df_impute$`T stage`)
levels(df_impute$`N stage`)

# Confirm conversions worked
sapply(df_impute[c("T stage", "Sex", "Race", "Primary sites", "N stage", "Age_stratum","Treatment")], class)
table(df_impute$`T stage`, useNA = "always")

# Compute the Nelson-Aalen estimator of cumulative hazard
df_impute$nelsonaalen_H <- nelsonaalen(df_impute, `Survival months`, event)
summary(df_impute$nelsonaalen_H)   

# Set up default predictor matrix and method list 
imp_setup <- mice(df_impute, maxit = 0)
meth <- imp_setup$method
meth["T stage"] <- "polyreg"
meth[names(meth) != "T stage"] <- ""   # only impute T stage

# Build predictor matrix: Sex, Age, Race, Primary sites, N stage, outcome components

pred_matrix <- imp_setup$predictorMatrix
pred_matrix[,] <- 0
pred_matrix["T stage", c("Sex", "Age_stratum", "Race", "Primary sites",
                         "N stage","Treatment", "nelsonaalen_H", "event")] <- 1

# Confirm T_stage_missing and N_stage_missing are NOT included
pred_matrix["T stage", c("T_stage_missing", "N_stage_missing")]   # should both be 0

# Run multiple imputation with M = 5
imp_polyreg <- mice(df_impute, method = meth, predictorMatrix = pred_matrix,
                    m = 5, maxit = 10, seed = 123)
imp_polyreg$loggedEvents   # should be NULL, or much shorter than before

# Check the imputed T-stage values are no longer NA
sum(is.na(complete(imp_polyreg, 1)$`T stage`))   # should be 0 now
table(complete(imp_polyreg, 1)$`T stage`, useNA = "always")
table(complete(imp_polyreg, 3)$`T stage`, useNA = "always")
table(complete(imp_polyreg, 5)$`T stage`, useNA = "always")

# Compare distribution of imputed values (previously-missing rows only) across imputations
for (i in 1:5) {
  print(table(complete(imp_polyreg, i)$`T stage`[df_impute$T_stage_missing]))
}

#Finding the confounders using Regression method 

# ---------------------------------------------------------
# STEP 0: Setup
# ---------------------------------------------------------
exposure_var <- "Treatment"
candidates <- c("`T stage`", "`N stage`", "Sex", "Age_stratum", "Race", "`Primary sites`")
outcome_lhs <- "Surv(`Survival months`, event)"

n_imp <- imp_polyreg$m   # should be 5
all_imputations <- complete(imp_polyreg, "all")  # list of 5 completed data frames

# ---------------------------------------------------------
# STEP 1: Helper function to fit a model formula across all imputations,
# pool the coefficient (log-HR) via Rubin's rules, then convert to HR
# ---------------------------------------------------------
fit_and_pool_HR <- function(formula, imputed_list, coef_pattern) {
  
  # Fit the same model on each imputed dataset
  model_list <- lapply(imputed_list, function(d) coxph(formula, data = d))
  
  # Convert to a 'mira' object so mice::pool() can combine them
  mira_obj <- as.mira(model_list)
  pooled <- pool(mira_obj)
  
  # Extract the pooled coefficient (this is on the log-HR scale)
  pooled_summary <- summary(pooled)
  coef_row <- grep(coef_pattern, pooled_summary$term, value = FALSE)
  
  if (length(coef_row) == 0) {
    stop(paste("No coefficient matched pattern:", coef_pattern))
  }
  
  pooled_logHR <- pooled_summary$estimate[coef_row]
  pooled_se    <- pooled_summary$std.error[coef_row]
  term_name    <- pooled_summary$term[coef_row]
  
  # Convert to Hazard Ratio scale
  pooled_HR    <- exp(pooled_logHR)
  HR_lower_95  <- exp(pooled_logHR - 1.96 * pooled_se)
  HR_upper_95  <- exp(pooled_logHR + 1.96 * pooled_se)
  
  list(
    term      = term_name,
    logHR     = pooled_logHR,
    HR        = pooled_HR,
    HR_lower  = HR_lower_95,
    HR_upper  = HR_upper_95
  )
}

# ---------------------------------------------------------
# STEP 2: Fit FULL model (all candidates included) across imputations, pool it
# ---------------------------------------------------------
full_formula <- as.formula(
  paste0(outcome_lhs, " ~ ", exposure_var, " + ", paste(candidates, collapse = " + "))
)

full_pooled <- fit_and_pool_HR(full_formula, all_imputations, coef_pattern = paste0("^", exposure_var))
HR_full <- full_pooled$HR

cat("Full model pooled exposure HR (", full_pooled$term, "):", round(HR_full, 4),
    " [95% CI:", round(full_pooled$HR_lower, 4), "-", round(full_pooled$HR_upper, 4), "]\n")

# ---------------------------------------------------------
# STEP 3: Loop over candidates, drop one at a time, refit + pool, compute % change in HR
# ---------------------------------------------------------
cie_results <- data.frame(
  Variable    = character(),
  HR_Full     = numeric(),
  HR_Reduced  = numeric(),
  Pct_Change  = numeric(),
  stringsAsFactors = FALSE
)

for (v in candidates) {
  
  remaining <- setdiff(candidates, v)
  
  reduced_formula <- as.formula(
    paste0(outcome_lhs, " ~ ", exposure_var, " + ", paste(remaining, collapse = " + "))
  )
  
  reduced_pooled <- fit_and_pool_HR(reduced_formula, all_imputations, coef_pattern = paste0("^", exposure_var))
  HR_reduced <- reduced_pooled$HR
  
  # % change computed on the HR scale, NOT the raw coefficient
  pct_change <- abs(HR_full - HR_reduced) / abs(HR_full - 1) * 100
  # Note: denominator uses (HR_full - 1) because HR=1 is the "null effect" reference point,
  # analogous to how 0 is the null reference point for a log-HR/coefficient.
  # Some papers instead just use HR_full as the denominator - see note below.
  
  cie_results <- rbind(cie_results, data.frame(
    Variable   = gsub("`", "", v),
    HR_Full    = HR_full,
    HR_Reduced = HR_reduced,
    Pct_Change = pct_change
  ))
  
  cat("Dropped:", gsub("`", "", v),
      "-> HR_reduced:", round(HR_reduced, 4),
      " % change:", round(pct_change, 2), "%\n")
}

# ---------------------------------------------------------
# STEP 4: Flag confounders using 10% cutoff (adjust per Lee 2014 guidance)
# Given large sample size (~29,000), a lower cutoff is more appropriate
# per Lee (2014): larger sample size -> lower cutoff needed
# ---------------------------------------------------------
cutoff <- 10  # primary/reported cutoff

cie_results <- cie_results %>%
  mutate(Flag_Confounder = Pct_Change > cutoff) %>%
  arrange(desc(Pct_Change))

print(cie_results)

# ---------------------------------------------------------
# STEP 5: Sensitivity check across multiple cutoffs
# Using 5%, 10%, 15% given large sample size (n ~ 29,000)
# ---------------------------------------------------------
cie_sensitivity <- cie_results %>%
  mutate(
    Flag_5pct  = Pct_Change > 5,
    Flag_10pct = Pct_Change > 10,
    Flag_15pct = Pct_Change > 15
  ) %>%
  select(Variable, HR_Full, HR_Reduced, Pct_Change, Flag_5pct, Flag_10pct, Flag_15pct)

print(cie_sensitivity)

#2 Propensity Score Matching
#Run PSM for 5 imputed datasets 
matched_imp <- matchthem(
  Treatment ~ Age_stratum + `Primary sites` + `T stage` + `N stage`,
  datasets = imp_polyreg,
  approach = "within",
  method = "nearest",
  distance = "glm",
  caliper = 0.2,
  std.caliper = TRUE,        # confirms caliper is in SD-of-logit-PS units (Austin 2011)
  exact = ~ Age_stratum,     # matching performed within each age stratum
  ratio = 1                  # 1:1 matching
)
summary(matched_imp)

# Standardized Mean Difference, before vs after, across all 5 imputed-matched datasets
bal.tab(matched_imp, abs = TRUE, thresholds = c(m = 0.10), un = TRUE)

# Love plot: before vs after matching
love.plot(matched_imp, threshold = 0.10, abs = TRUE, var.order = "unadjusted")


# Standard Rubin's Rules pooling (normal approximation)
# est_mat / se_mat: M x p matrices (rows = imputations, cols = terms)
pool_rubin_simple <- function(est_mat, se_mat) {
  M <- nrow(est_mat)
  
  pooled_est  <- colMeans(est_mat)
  within_var  <- colMeans(se_mat^2)
  between_var <- apply(est_mat, 2, var)
  total_var   <- within_var + (1 + 1/M) * between_var
  pooled_se   <- sqrt(total_var)
  
  pooled_p <- 2 * pnorm(-abs(pooled_est / pooled_se))
  
  data.frame(
    term      = names(pooled_est),
    estimate  = pooled_est,
    HR        = exp(pooled_est),
    std.error = pooled_se,
    CI_lower  = exp(pooled_est - 1.96 * pooled_se),
    CI_upper  = exp(pooled_est + 1.96 * pooled_se),
    p.value   = pooled_p,
    row.names = NULL
  )
}


#3 Univariate Cox Regression Model

covariates <- c( "Age_stratum", "`Primary sites`", "`T stage`", "`N stage`")

univariate_results <- list()
M <- 5
matched_list <- complete(matched_imp, action = "all")
for (cov in covariates) {
  form <- as.formula(paste0("Surv(`Survival months`, event) ~ ", cov))
  
  log_hr_list <- list()
  se_hr_list <- list()
  
  for (i in 1:M) {
    d <- matched_list[[i]]
    d <- d[d$weights > 0, ]
    
    fit <- coxph(form, data = d, weights = weights)
    log_hr_list[[i]] <- coef(fit)
    se_hr_list[[i]] <- sqrt(diag(vcov(fit)))
  }
  
  est_mat <- do.call(rbind, log_hr_list)
  se_mat  <- do.call(rbind, se_hr_list)
  
  univariate_results[[cov]] <- pool_rubin_simple(est_mat, se_mat)
  cat("Completed:", cov, "\n")
}

# View all univariate results together
univariate_results

# Quick check: which covariates have ANY term with p < 0.10?
for (cov in names(univariate_results)) {
  min_p <- min(univariate_results[[cov]]$p.value)
  cat(cov, "- minimum p-value across levels:", round(min_p, 4), 
      ifelse(min_p < 0.10, "-> ENTERS multivariate model", "-> excluded"), "\n")
}

#4 Primary analysis - Interaction model(Cox Regression)

selected_covariates <- c("`Primary sites`", "`T stage`", "`N stage`")

final_formula <- as.formula(
  paste0("Surv(`Survival months`, event) ~ Treatment * Age_stratum + ",
         paste(selected_covariates, collapse = " + "))
)

print(final_formula)   # confirm the formula looks right before fitting

log_hr_list <- list()
se_hr_list <- list()

for (i in 1:M) {
  d <- matched_list[[i]]
  d <- d[d$weights > 0, ]
  
  fit <- coxph(final_formula, data = d, weights = weights, robust = TRUE)
  
  # Safety check: flag any dataset with estimation problems
  if (any(is.na(coef(fit)))) {
    cat("WARNING: NA coefficient(s) in imputation", i, "\n")
  }
  
  log_hr_list[[i]] <- coef(fit)
  se_hr_list[[i]] <- sqrt(diag(vcov(fit)))
}

est_mat <- do.call(rbind, log_hr_list)
se_mat  <- do.call(rbind, se_hr_list)

final_pooled <- pool_rubin_simple(est_mat, se_mat)

# Full model results table
print(final_pooled[, c("term", "HR", "CI_lower", "CI_upper", "p.value")], row.names = FALSE)

interaction_terms <- final_pooled[grepl("Treatment.*Age_stratum", final_pooled$term), ]
print(interaction_terms[, c("term", "HR", "CI_lower", "CI_upper", "p.value")], row.names = FALSE)

#5 Secondary analysis - Main effects model (no interaction term)

secondary_formula <- as.formula(
  paste0("Surv(`Survival months`, event) ~ Treatment + Age_stratum + ",
         paste(selected_covariates, collapse = " + "))
)

print(secondary_formula)   # confirm the formula looks right before fitting

log_hr_list_secondary <- list()
se_hr_list_secondary  <- list()

for (i in 1:M) {
  d <- matched_list[[i]]
  d <- d[d$weights > 0, ]
  d$Treatment <- factor(d$Treatment, levels = c("Surgery only", "Surgery + PORT"))  # safeguard
  
  fit <- coxph(secondary_formula, data = d, weights = weights, robust = TRUE)
  
  # Safety check: flag any dataset with estimation problems
  if (any(is.na(coef(fit)))) {
    cat("WARNING: NA coefficient(s) in imputation", i, "\n")
  }
  
  log_hr_list_secondary[[i]] <- coef(fit)
  se_hr_list_secondary[[i]]  <- sqrt(diag(vcov(fit)))
}

est_mat_secondary <- do.call(rbind, log_hr_list_secondary)
se_mat_secondary  <- do.call(rbind, se_hr_list_secondary)

secondary_pooled <- pool_rubin_simple(est_mat_secondary, se_mat_secondary)

# Full model results table
print(secondary_pooled[, c("term", "HR", "CI_lower", "CI_upper", "p.value")], row.names = FALSE)

# Extract just the main Treatment effect (PORT vs Surgery only) — this is now your
# single pooled HR representative across all age strata, since no interaction is modeled
treatment_effect <- secondary_pooled[grepl("^Treatment", secondary_pooled$term), ]
print(treatment_effect[, c("term", "HR", "CI_lower", "CI_upper", "p.value")], row.names = FALSE)


#6 PH check — Primary (interaction) model

schoenfeld_list <- vector("list", M)
zph_list        <- vector("list", M)
fit_list        <- vector("list", M)

for (i in 1:M) {
  
  d <- matched_list[[i]]
  d <- d[d$weights > 0, ]
  d$Treatment <- factor(d$Treatment, levels = c("Surgery only", "Surgery + PORT"))  # safeguard, added
  
  fit <- coxph(
    final_formula,
    data = d,
    weights = weights,
    robust = TRUE
  )
  
  if (any(is.na(coef(fit)))) {
    cat("WARNING: NA coefficient(s) in imputation", i, "- skipping\n")
    next
  }
  
  fit_list[[i]] <- fit
  schoenfeld_list[[i]] <- residuals(fit, type = "schoenfeld")
  zph_list[[i]] <- cox.zph(fit)
  
  cat("Imputation", i, "- extracted", nrow(schoenfeld_list[[i]]), "Schoenfeld residuals\n")
}

# Find the first NON-NULL zph result, in case dataset 1 was skipped
first_valid <- which(!sapply(zph_list, is.null))[1]

if (is.na(first_valid)) {
  stop("No valid cox.zph results were produced — check for errors above.")
}

print(zph_list[[first_valid]])
resid_terms <- colnames(zph_list[[first_valid]]$y)

# Adjust grid size to match M (not hardcoded 2x3)
n_cols <- min(3, M)
n_rows <- ceiling(M / n_cols)

for (term in resid_terms) {
  
  safe_name <- gsub("[^A-Za-z0-9]+", "_", term)
  
  png(
    filename = paste0("Schoenfeld_Primary_", safe_name, ".png"),
    width = 1200,
    height = 800,
    res = 120
  )
  
  par(mfrow = c(n_rows, n_cols))
  
  for (i in 1:M) {
    if (!is.null(zph_list[[i]])) {
      plot(
        zph_list[[i]],
        var = term,
        resid = TRUE,
        main = paste("Imputation", i),
        col = "grey70",
        lwd = 2
      )
      abline(h = 0, col = "red", lty = 2)
    }
  }
  
  par(mfrow = c(1, 1))
  dev.off()
}

# Numeric summary: global PH test across all valid imputations
global_pvals_primary <- sapply(zph_list, function(z) if (!is.null(z)) z$table["GLOBAL", "p"] else NA)
cat("Primary model - PH violated (p<0.05) in",
    sum(global_pvals_primary < 0.05, na.rm = TRUE), "out of",
    sum(!is.na(global_pvals_primary)), "valid imputed datasets\n")
print(global_pvals_primary)


# PH check — Secondary (main-effects) model

schoenfeld_list_sec <- vector("list", M)
zph_list_sec        <- vector("list", M)
fit_list_sec        <- vector("list", M)

for (i in 1:M) {
  
  d <- matched_list[[i]]
  d <- d[d$weights > 0, ]
  d$Treatment <- factor(d$Treatment, levels = c("Surgery only", "Surgery + PORT"))  # safeguard
  
  fit <- coxph(
    secondary_formula,
    data = d,
    weights = weights,
    robust = TRUE
  )
  
  if (any(is.na(coef(fit)))) {
    cat("WARNING: NA coefficient(s) in imputation", i, "- skipping\n")
    next
  }
  
  fit_list_sec[[i]] <- fit
  schoenfeld_list_sec[[i]] <- residuals(fit, type = "schoenfeld")
  zph_list_sec[[i]] <- cox.zph(fit)
  
  cat("Imputation", i, "- extracted", nrow(schoenfeld_list_sec[[i]]), "Schoenfeld residuals\n")
}

first_valid_sec <- which(!sapply(zph_list_sec, is.null))[1]

if (is.na(first_valid_sec)) {
  stop("No valid cox.zph results were produced for the secondary model — check for errors above.")
}

print(zph_list_sec[[first_valid_sec]])
resid_terms_sec <- colnames(zph_list_sec[[first_valid_sec]]$y)

for (term in resid_terms_sec) {
  
  safe_name <- gsub("[^A-Za-z0-9]+", "_", term)
  
  png(
    filename = paste0("Schoenfeld_Secondary_", safe_name, ".png"),
    width = 1200,
    height = 800,
    res = 120
  )
  
  par(mfrow = c(n_rows, n_cols))
  
  for (i in 1:M) {
    if (!is.null(zph_list_sec[[i]])) {
      plot(
        zph_list_sec[[i]],
        var = term,
        resid = TRUE,
        main = paste("Imputation", i),
        col = "grey70",
        lwd = 2
      )
      abline(h = 0, col = "red", lty = 2)
    }
  }
  
  par(mfrow = c(1, 1))
  dev.off()
}

global_pvals_secondary <- sapply(zph_list_sec, function(z) if (!is.null(z)) z$table["GLOBAL", "p"] else NA)
cat("Secondary model - PH violated (p<0.05) in",
    sum(global_pvals_secondary < 0.05, na.rm = TRUE), "out of",
    sum(!is.na(global_pvals_secondary)), "valid imputed datasets\n")
print(global_pvals_secondary)

#RMST --For 5 year( only for Treatment)

M <- 5
tau <- 60   # 5 years in months

# ===========================================================
# Overall RMST (Treatment effect, unadjusted, across all patients)
# ===========================================================

rmst_results <- vector("list", M)

for (i in 1:M) {
  d <- matched_list[[i]]
  d <- d[d$weights > 0, ]
  d$Treatment <- factor(d$Treatment, levels = c("Surgery only", "Surgery + PORT"))
  
  arm_numeric <- ifelse(d$Treatment == "Surgery + PORT", 1, 0)
  
  rmst_fit <- rmst2(
    time   = d$`Survival months`,
    status = d$event,
    arm    = arm_numeric,
    tau    = tau
  )
  
  rmst_results[[i]] <- rmst_fit
  cat("Completed RMST for imputation", i, "\n")
}

# Check the row names of the result table once, before extracting
print(rmst_results[[1]]$unadjusted.result)

# ---------------------------------------------------------
# Extract difference estimate + SE per imputation
# ---------------------------------------------------------
rmst_diff_est <- numeric(M)
rmst_diff_se  <- numeric(M)

for (i in 1:M) {
  tbl <- rmst_results[[i]]$unadjusted.result
  
  rmst_diff_est[i] <- tbl["RMST (arm=1)-(arm=0)", "Est."]
  lower <- tbl["RMST (arm=1)-(arm=0)", "lower .95"]
  upper <- tbl["RMST (arm=1)-(arm=0)", "upper .95"]
  rmst_diff_se[i]  <- (upper - lower) / (2 * 1.96)
}

rmst_by_imputation <- data.frame(
  Imputation = paste0("Imputation_", 1:M),
  RMST_diff  = round(rmst_diff_est, 3),
  SE         = round(rmst_diff_se, 3)
)

print(rmst_by_imputation, row.names = FALSE)

# ---------------------------------------------------------
# Pool RMST results across imputations (Rubin's rules)
# ---------------------------------------------------------
pool_rubin <- function(est_vec, se_vec) {
  M <- length(est_vec)
  
  Q_bar <- mean(est_vec)               # pooled estimate
  U_bar <- mean(se_vec^2)              # within-imputation variance
  B     <- var(est_vec)                # between-imputation variance
  T_var <- U_bar + (1 + 1/M) * B       # total variance (Rubin 1987)
  pooled_se <- sqrt(T_var)
  
  df <- (M - 1) * (1 + U_bar / ((1 + 1/M) * B))^2
  
  t_crit   <- qt(0.975, df)
  CI_lower <- Q_bar - t_crit * pooled_se
  CI_upper <- Q_bar + t_crit * pooled_se
  t_stat   <- Q_bar / pooled_se
  p_value  <- 2 * pt(-abs(t_stat), df)
  
  data.frame(
    estimate = Q_bar,
    se       = pooled_se,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    df       = df,
    p.value  = p_value
  )
}

pooled_result <- pool_rubin(rmst_diff_est, rmst_diff_se)
cat("\nOverall pooled RMST difference (Surgery+PORT vs Surgery only), tau =", tau, "months:\n")
print(pooled_result)



# Record the R and package environment for reproducibility.
# Save this console output with the analysis notes when reporting results.
cat("\n--- R session information ---\n")
print(sessionInfo())
