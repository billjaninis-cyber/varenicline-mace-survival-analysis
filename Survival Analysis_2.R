```r
############################################################
# Comparative Survival Analysis:
# Varenicline vs Bupropion and Risk of MACE
#
# Purpose:
# Evaluate differences in time to major adverse cardiac events
# (MACE) among patients receiving varenicline versus bupropion.
#
# Methods:
# - Data review and missingness assessment
# - Descriptive statistics
# - Crude incidence rates
# - Kaplan-Meier survival analysis
# - Cox proportional hazards regression
# - Proportional hazards diagnostics
# - Propensity score estimation
# - 1:1 nearest-neighbor propensity score matching
# - Covariate balance assessment
# - Cluster-robust Cox regression in matched cohort
#
# Author: Vasilios (Bill) Janinis
############################################################


############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

library(tableone)
library(epiDisplay)
library(survival)
library(survminer)
library(ggplot2)
library(MatchIt)
library(cobalt)


############################################################
# 2. IMPORT DATA
############################################################

# The original dataset is not included in this repository.
# Place the dataset in the project directory before running.

d1 <- read.csv(
  "analysis_assignment.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)


############################################################
# 3. INITIAL DATA REVIEW
############################################################

# Examine structure and sample observations

str(d1)

head(d1)

summary(d1)


############################################################
# 4. DATA PREPARATION
############################################################

# Convert outcome, follow-up time, and treatment variables
# to numeric format.

d1$survival_time <- as.numeric(d1$survival_time)

d1$MACE <- as.numeric(d1$MACE)

d1$varenicline <- as.numeric(d1$varenicline)


# Create treatment variable
#
# 0 = Bupropion
# 1 = Varenicline

d1$treatment <- d1$varenicline


# Create readable treatment labels

d1$treatment_label <- factor(
  d1$treatment,
  levels = c(0, 1),
  labels = c("Bupropion", "Varenicline")
)


# Convert categorical covariates to factors

factor_vars <- c(
  "sex",
  "Race",
  "T2DM",
  "Coronary_revascularization",
  "Myocardial_infaction",
  "Stable_angina",
  "Statins"
)

d1[factor_vars] <- lapply(
  d1[factor_vars],
  factor
)


############################################################
# 5. DEFINE ANALYSIS VARIABLES
############################################################

# Variables displayed in descriptive Table 1

all_vars <- c(
  "Age",
  "sex",
  "Race",
  "CCI",
  "FRAILTY",
  "T2DM",
  "Coronary_revascularization",
  "Myocardial_infaction",
  "Stable_angina",
  "Statins",
  "Vital_BMI",
  "Vital_Systolic",
  "Pack_Year_Aggregated"
)


# Variables treated as categorical in TableOne

categorical_vars <- c(
  "sex",
  "Race",
  "T2DM",
  "Coronary_revascularization",
  "Myocardial_infaction",
  "Stable_angina",
  "Statins"
)


############################################################
# 5a. MISSINGNESS ASSESSMENT
############################################################

# Define all variables required for the primary analyses.

required_vars <- c(
  all_vars,
  "treatment",
  "survival_time",
  "MACE"
)


# Identify incomplete observations.

n_incomplete <- sum(
  !complete.cases(d1[, required_vars])
)

cat(
  "\nRows with missing data on required analysis variables: ",
  n_incomplete,
  " of ",
  nrow(d1),
  " (",
  round(100 * n_incomplete / nrow(d1), 2),
  "%)\n",
  sep = ""
)


############################################################
# 5b. CREATE COMPLETE-CASE ANALYTIC COHORT
############################################################

# Use the same analytic population throughout the primary
# survival and propensity-score analyses.
#
# This prevents different models from silently using
# different sample sizes because of missing data.

analysis_data <- d1[
  complete.cases(d1[, required_vars]),
]


cat(
  "\nOriginal sample size: ",
  nrow(d1),
  "\nComplete-case analytic sample size: ",
  nrow(analysis_data),
  "\n",
  sep = ""
)


############################################################
# 5c. RACE CATEGORY / SEPARATION CHECK
############################################################

# Check for sparse Race categories within treatment groups.
#
# Very small categories can cause quasi-complete separation
# in logistic regression and unstable coefficients in Cox
# regression.

race_by_treatment <- table(
  analysis_data$Race,
  analysis_data$treatment_label
)

cat("\nRace distribution by treatment arm:\n")

print(race_by_treatment)


# Flag a Race level if fewer than 10 patients are present
# in either treatment group.

sparse_threshold <- 10

sparse_race_flag <- apply(
  race_by_treatment,
  1,
  function(x) any(x < sparse_threshold)
)


if (any(sparse_race_flag)) {
  
  cat(
    "\nWARNING: The following Race level(s) have fewer than ",
    sparse_threshold,
    " patients in at least one treatment arm:\n",
    sep = ""
  )
  
  print(
    rownames(race_by_treatment)[sparse_race_flag]
  )
  
  cat(
    "\nThese categories will be collapsed into 'Other' ",
    "for regression and propensity-score modeling.\n"
  )
}


############################################################
# 5d. CREATE COLLAPSED RACE VARIABLE FOR MODELING
############################################################

# Preserve the original Race variable for descriptive tables.
#
# Use Race_collapsed only for regression and propensity-score
# modeling to reduce instability from sparse categories.

analysis_data$Race_collapsed <- as.character(
  analysis_data$Race
)


sparse_levels <- rownames(
  race_by_treatment
)[sparse_race_flag]


analysis_data$Race_collapsed[
  analysis_data$Race_collapsed %in% sparse_levels
] <- "Other"


analysis_data$Race_collapsed <- factor(
  analysis_data$Race_collapsed
)


cat(
  "\nRace categories used in regression/PS models:\n"
)

print(
  table(
    analysis_data$Race_collapsed,
    analysis_data$treatment_label
  )
)


############################################################
# 6. BASELINE CHARACTERISTICS
############################################################

# NOTE:
# Original Race categories are displayed in Table 1.
#
# Race_collapsed is used only in regression and propensity
# score modeling.


############################################################
# 6a. OVERALL COHORT
############################################################

overall_table <- CreateTableOne(
  vars = all_vars,
  data = analysis_data,
  factorVars = categorical_vars
)


print(
  overall_table,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE
)


############################################################
# 6b. BASELINE CHARACTERISTICS BY TREATMENT
############################################################

treatment_table <- CreateTableOne(
  vars = all_vars,
  strata = "treatment_label",
  data = analysis_data,
  factorVars = categorical_vars
)


print(
  treatment_table,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  smd = TRUE
)


############################################################
# 7. CRUDE INCIDENCE RATES OF MACE
############################################################

# Calculate number of MACE events by treatment group.

events <- rowsum(
  analysis_data$MACE,
  analysis_data$treatment
)


# Convert follow-up time from days to person-years.

person_years <- rowsum(
  analysis_data$survival_time / 365.25,
  analysis_data$treatment
)


# Estimate crude incidence rates and 95% confidence intervals.

incidence_rates <- ci.poisson(
  events,
  person_years,
  alpha = 0.05
)


# Label treatment groups explicitly.

rownames(incidence_rates) <- c(
  "Bupropion",
  "Varenicline"
)


cat("\nCrude MACE incidence rates:\n")

print(incidence_rates)


############################################################
# 8. KAPLAN-MEIER SURVIVAL ANALYSIS
############################################################

# Fit Kaplan-Meier curves for MACE-free survival.

km_fit <- survfit(
  Surv(survival_time, MACE) ~ treatment_label,
  data = analysis_data
)


# Display Kaplan-Meier estimates.

summary(km_fit)


############################################################
# 8a. KAPLAN-MEIER PLOT
############################################################

km_plot <- ggsurvplot(
  km_fit,
  data = analysis_data,
  ylim = c(0.5, 1),
  xlim = c(0, 365),
  break.time.by = 90,
  censor = FALSE,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.title = "Treatment Group",
  legend.labs = c(
    "Bupropion",
    "Varenicline"
  ),
  title = "Kaplan-Meier Curves for Time to MACE",
  xlab = "Follow-up Time (Days)",
  ylab = "MACE-Free Survival Probability",
  linetype = c(
    "solid",
    "dashed"
  ),
  size = 0.7
)


print(km_plot)


############################################################
# 9. UNADJUSTED COX PROPORTIONAL HAZARDS MODEL
############################################################

# HR > 1:
# Higher hazard of MACE with varenicline relative to bupropion.
#
# HR < 1:
# Lower hazard of MACE with varenicline relative to bupropion.

unadjusted_cox <- coxph(
  Surv(survival_time, MACE) ~ treatment,
  data = analysis_data,
  ties = "breslow"
)


summary(unadjusted_cox)


############################################################
# 10. MULTIVARIABLE-ADJUSTED COX MODEL
############################################################

# Race_collapsed is used instead of Race to reduce instability
# caused by sparse categories.

adjusted_cox <- coxph(
  Surv(survival_time, MACE) ~
    treatment +
    Age +
    Race_collapsed +
    sex +
    Myocardial_infaction +
    T2DM +
    Coronary_revascularization +
    Stable_angina +
    Statins +
    CCI +
    FRAILTY +
    Vital_BMI +
    Vital_Systolic +
    Pack_Year_Aggregated,
  data = analysis_data,
  ties = "breslow"
)


summary(adjusted_cox)


############################################################
# 11. PROPORTIONAL HAZARDS ASSUMPTION
############################################################

# Test proportional hazards assumption using
# Schoenfeld residuals.

ph_test <- cox.zph(
  adjusted_cox
)


print(ph_test)


# Plot scaled Schoenfeld residuals.

plot(ph_test)


############################################################
# 12. PROPENSITY SCORE ESTIMATION
############################################################

# Estimate the probability of receiving varenicline using
# multivariable logistic regression.
#
# Only pre-treatment baseline covariates should be included
# in a propensity-score model.

ps_model <- glm(
  treatment ~
    Age +
    Race_collapsed +
    sex +
    Myocardial_infaction +
    T2DM +
    Coronary_revascularization +
    Stable_angina +
    Statins +
    CCI +
    FRAILTY +
    Vital_BMI +
    Vital_Systolic +
    Pack_Year_Aggregated,
  family = binomial(
    link = "logit"
  ),
  data = analysis_data
)


summary(ps_model)


############################################################
# 12a. SAVE PROPENSITY SCORES
############################################################

# Because analysis_data is already a complete-case cohort,
# the prediction vector has the same number of observations
# as analysis_data.

analysis_data$propensity_score <- predict(
  ps_model,
  newdata = analysis_data,
  type = "response"
)


# Examine PS distribution numerically.

summary(
  analysis_data$propensity_score
)


############################################################
# 13. PROPENSITY SCORE DISTRIBUTION BEFORE MATCHING
############################################################

ps_before_plot <- ggplot(
  analysis_data,
  aes(
    x = propensity_score,
    fill = treatment_label,
    y = after_stat(scaled)
  )
) +
  geom_density(
    alpha = 0.6
  ) +
  scale_x_continuous(
    name = "Propensity Score",
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    name = "Scaled Density"
  ) +
  labs(
    fill = "Treatment Group",
    title = "Propensity Score Distribution Before Matching"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(
      size = 12
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    )
  )


print(ps_before_plot)


############################################################
# 14. PROPENSITY SCORE MATCHING
############################################################

# Perform 1:1 nearest-neighbor propensity-score matching.
#
# IMPORTANT:
# The numeric propensity scores estimated above are supplied
# directly to MatchIt.
#
# Therefore, the same PS model is used for:
#
# 1. The pre-matching PS distribution
# 2. Matching
# 3. The post-matching PS distribution
#
# This avoids fitting the propensity score model twice.
#
# estimand = "ATT":
# Estimates the treatment effect among patients who received
# varenicline.
#
# ratio = 1:
# One bupropion patient per varenicline patient.
#
# replace = FALSE:
# Each control patient can be used only once.
#
# caliper = 0.1:
# Restricts matches to patients with sufficiently similar
# propensity scores.
#
# std.caliper = TRUE:
# The caliper is measured in standard deviation units.

match_model <- matchit(
  treatment ~
    Age +
    Race_collapsed +
    sex +
    Myocardial_infaction +
    T2DM +
    Coronary_revascularization +
    Stable_angina +
    Statins +
    CCI +
    FRAILTY +
    Vital_BMI +
    Vital_Systolic +
    Pack_Year_Aggregated,
  data = analysis_data,
  method = "nearest",
  distance = analysis_data$propensity_score,
  estimand = "ATT",
  ratio = 1,
  caliper = 0.1,
  std.caliper = TRUE,
  replace = FALSE
)


############################################################
# 14a. REVIEW MATCHING RESULTS
############################################################

summary(match_model)


# Report matched, unmatched, and discarded sample sizes.

match_nn <- summary(match_model)$nn


cat(
  "\nMatching sample sizes:\n"
)

print(match_nn)


############################################################
# 14b. EXTRACT MATCHED COHORT
############################################################

matched_data <- match.data(
  match_model
)


# Re-create readable treatment label.

matched_data$treatment_label <- factor(
  matched_data$treatment,
  levels = c(0, 1),
  labels = c(
    "Bupropion",
    "Varenicline"
  )
)


cat(
  "\nMatched treatment sample sizes:\n"
)

print(
  table(
    matched_data$treatment_label
  )
)


############################################################
# 15. COVARIATE BALANCE AFTER MATCHING
############################################################

# Descriptive table after matching.
#
# Again, the original Race categories are displayed for
# descriptive interpretation.

matched_table <- CreateTableOne(
  vars = all_vars,
  strata = "treatment_label",
  data = matched_data,
  factorVars = categorical_vars
)


print(
  matched_table,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  smd = TRUE
)


############################################################
# 16. COVARIATE BALANCE BEFORE AND AFTER MATCHING
############################################################

# Standardized mean differences are preferred over
# significance tests for evaluating balance.

balance_results <- bal.tab(
  match_model,
  un = TRUE,
  thresholds = c(
    m = 0.1
  )
)


print(balance_results)


############################################################
# 16a. LOVE PLOT
############################################################

# Absolute SMD < 0.10 is commonly interpreted as evidence
# of acceptable covariate balance.

love_plot <- love.plot(
  match_model,
  stats = "mean.diffs",
  abs = TRUE,
  threshold = 0.1,
  var.order = "unadjusted",
  binary = "std",
  sample.names = c(
    "Before Matching",
    "After Matching"
  ),
  title =
    "Covariate Balance Before and After Propensity Score Matching"
)


print(love_plot)


############################################################
# 17. PROPENSITY SCORE DISTRIBUTION AFTER MATCHING
############################################################

# MatchIt stores the propensity score used for matching
# in the distance variable in the matched dataset.

ps_after_plot <- ggplot(
  matched_data,
  aes(
    x = distance,
    fill = treatment_label,
    y = after_stat(scaled)
  )
) +
  geom_density(
    alpha = 0.6
  ) +
  scale_x_continuous(
    name = "Propensity Score",
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    name = "Scaled Density"
  ) +
  labs(
    fill = "Treatment Group",
    title =
      "Propensity Score Distribution After 1:1 Matching"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(
      size = 12
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    )
  )


print(ps_after_plot)


############################################################
# 18. COX MODEL IN PROPENSITY-SCORE-MATCHED COHORT
############################################################

# MatchIt assigns each matched pair to a subclass.
#
# Observations within a matched pair are not independent.
#
# cluster = subclass therefore provides pair-aware
# robust standard errors.

matched_data$subclass <- factor(
  matched_data$subclass
)


matched_cox <- coxph(
  Surv(survival_time, MACE) ~ treatment,
  data = matched_data,
  weights = weights,
  cluster = subclass,
  ties = "breslow"
)


summary(matched_cox)


############################################################
# 19. PROPORTIONAL HAZARDS CHECK IN MATCHED COHORT
############################################################

matched_ph_test <- cox.zph(
  matched_cox
)


print(matched_ph_test)


plot(matched_ph_test)


############################################################
# 20. EXTRACT MODEL RESULTS
############################################################

############################################################
# 20a. UNADJUSTED COX MODEL
############################################################

unadjusted_results <- data.frame(
  HR = exp(
    coef(unadjusted_cox)
  ),
  exp(
    confint(unadjusted_cox)
  )
)


colnames(unadjusted_results) <- c(
  "HR",
  "CI_lower",
  "CI_upper"
)


############################################################
# 20b. MULTIVARIABLE-ADJUSTED COX MODEL
############################################################

adjusted_results <- data.frame(
  HR = exp(
    coef(adjusted_cox)
  ),
  exp(
    confint(adjusted_cox)
  )
)


colnames(adjusted_results) <- c(
  "HR",
  "CI_lower",
  "CI_upper"
)


############################################################
# 20c. PROPENSITY-SCORE-MATCHED COX MODEL
############################################################

matched_results <- data.frame(
  HR = exp(
    coef(matched_cox)
  ),
  exp(
    confint(matched_cox)
  )
)


colnames(matched_results) <- c(
  "HR",
  "CI_lower",
  "CI_upper"
)


############################################################
# 20d. DISPLAY MODEL RESULTS
############################################################

cat(
  "\nUnadjusted Cox model:\n"
)

print(
  unadjusted_results
)


cat(
  "\nMultivariable-adjusted Cox model:\n"
)

print(
  adjusted_results
)


cat(
  "\nPropensity-score-matched Cox model:\n"
)

print(
  matched_results
)


############################################################
# 21. PRIMARY TREATMENT EFFECT SUMMARY
############################################################

# Extract treatment effect from each model.
#
# Because:
#
# treatment = 0 -> Bupropion
# treatment = 1 -> Varenicline
#
# these HRs represent:
#
# Varenicline vs Bupropion


unadjusted_treatment <- unadjusted_results[
  "treatment",
  ,
  drop = FALSE
]


adjusted_treatment <- adjusted_results[
  "treatment",
  ,
  drop = FALSE
]


matched_treatment <- matched_results[
  "treatment",
  ,
  drop = FALSE
]


treatment_effect_summary <- rbind(
  Unadjusted = unadjusted_treatment,
  Adjusted = adjusted_treatment,
  PS_Matched = matched_treatment
)


cat(
  "\nTreatment Effect Summary: Varenicline vs Bupropion\n"
)


print(
  treatment_effect_summary
)


############################################################
# 22. SAMPLE SIZE SUMMARY
############################################################

sample_size_summary <- data.frame(
  Population = c(
    "Original dataset",
    "Complete-case analytic cohort",
    "Propensity-score-matched cohort"
  ),
  N = c(
    nrow(d1),
    nrow(analysis_data),
    nrow(matched_data)
  )
)


cat(
  "\nSample Size Summary:\n"
)


print(
  sample_size_summary,
  row.names = FALSE
)


############################################################
# 23. EVENT SUMMARY
############################################################

event_summary <- aggregate(
  MACE ~ treatment_label,
  data = analysis_data,
  FUN = function(x) {
    c(
      Events = sum(x),
      Total = length(x),
      Percent = 100 * mean(x)
    )
  }
)


cat(
  "\nMACE Event Summary:\n"
)


print(event_summary)


############################################################
# 24. SESSION INFORMATION
############################################################

# Record R and package versions to support reproducibility.

sessionInfo()