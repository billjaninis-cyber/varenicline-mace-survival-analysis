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
# - Descriptive statistics
# - Crude incidence rates
# - Kaplan-Meier survival analysis
# - Cox proportional hazards regression
# - Propensity score estimation
# - 1:1 nearest-neighbor propensity score matching
# - Covariate balance assessment
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

# Convert key variables to appropriate formats

d1$survival_time <- as.numeric(d1$survival_time)

d1$MACE <- as.numeric(d1$MACE)

d1$varenicline <- as.numeric(d1$varenicline)


# Create treatment variables

# 0 = Bupropion
# 1 = Varenicline

d1$treatment <- d1$varenicline

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
# 6. BASELINE CHARACTERISTICS
############################################################

# Overall cohort

overall_table <- CreateTableOne(
  vars = all_vars,
  data = d1,
  factorVars = categorical_vars
)

print(
  overall_table,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE
)


# Baseline characteristics stratified by treatment

treatment_table <- CreateTableOne(
  vars = all_vars,
  strata = "treatment_label",
  data = d1,
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

# Calculate number of MACE events by treatment group

events <- rowsum(
  d1$MACE,
  d1$treatment
)


# Convert follow-up time from days to person-years

person_years <- rowsum(
  d1$survival_time / 365.25,
  d1$treatment
)


# Estimate crude incidence rates and 95% confidence intervals

incidence_rates <- ci.poisson(
  events,
  person_years,
  alpha = 0.05
)

incidence_rates


############################################################
# 8. KAPLAN-MEIER SURVIVAL ANALYSIS
############################################################

# Fit Kaplan-Meier survival curves by treatment

km_fit <- survfit(
  Surv(survival_time, MACE) ~ treatment_label,
  data = d1
)


# Display Kaplan-Meier estimates

summary(km_fit)


# Plot Kaplan-Meier survival curves

km_plot <- ggsurvplot(
  km_fit,
  data = d1,
  ylim = c(0.5, 1),
  xlim = c(0, 365),
  break.time.by = 90,
  censor = FALSE,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.title = "Treatment Group",
  legend.labs = c("Bupropion", "Varenicline"),
  title = "Kaplan-Meier Curves for Time to MACE",
  xlab = "Follow-up Time (Days)",
  ylab = "MACE-Free Survival Probability",
  linetype = c("solid", "dashed"),
  size = 0.7
)

print(km_plot)


############################################################
# 9. UNADJUSTED COX PROPORTIONAL HAZARDS MODEL
############################################################

unadjusted_cox <- coxph(
  Surv(survival_time, MACE) ~ treatment,
  data = d1,
  ties = "breslow"
)

summary(unadjusted_cox)


############################################################
# 10. MULTIVARIABLE-ADJUSTED COX MODEL
############################################################

adjusted_cox <- coxph(
  Surv(survival_time, MACE) ~
    treatment +
    Age +
    Race +
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
  data = d1,
  ties = "breslow"
)

summary(adjusted_cox)


############################################################
# 11. PROPORTIONAL HAZARDS ASSUMPTION
############################################################

# Evaluate proportional hazards assumption using
# Schoenfeld residuals

ph_test <- cox.zph(adjusted_cox)

print(ph_test)

plot(ph_test)


############################################################
# 12. PROPENSITY SCORE ESTIMATION
############################################################

# Estimate probability of receiving varenicline using
# multivariable logistic regression

ps_model <- glm(
  treatment ~
    Age +
    Race +
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
  family = binomial(link = "logit"),
  data = d1
)

summary(ps_model)


# Save predicted propensity scores

d1$propensity_score <- predict(
  ps_model,
  type = "response"
)


############################################################
# 13. PROPENSITY SCORE DISTRIBUTION BEFORE MATCHING
############################################################

ps_before_plot <- ggplot(
  d1,
  aes(
    x = propensity_score,
    fill = treatment_label,
    y = after_stat(scaled)
  )
) +
  geom_density(alpha = 0.6) +
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
    axis.title = element_text(size = 12),
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

# Perform 1:1 nearest-neighbor propensity score matching

match_model <- matchit(
  treatment ~
    Age +
    Race +
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
  data = d1,
  method = "nearest",
  distance = "glm",
  caliper = 0.1,
  replace = FALSE,
  na.action = na.omit
)


# Review matching results

summary(match_model)


# Extract matched cohort

matched_data <- match.data(match_model)


# Create readable treatment label

matched_data$treatment_label <- factor(
  matched_data$treatment,
  levels = c(0, 1),
  labels = c("Bupropion", "Varenicline")
)


# Display matched sample sizes

table(matched_data$treatment_label)


############################################################
# 15. COVARIATE BALANCE AFTER MATCHING
############################################################

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
# 16. LOVE PLOT: COVARIATE BALANCE BEFORE AND AFTER MATCHING
############################################################

# Calculate covariate balance before and after matching

balance_results <- bal.tab(
  match_model,
  un = TRUE,
  thresholds = c(m = 0.1)
)

# Display numerical balance results

print(balance_results)


# Create Love plot using absolute standardized mean differences

love_plot <- love.plot(
  match_model,
  stats = "mean.diffs",
  abs = TRUE,
  threshold = 0.1,
  var.order = "unadjusted",
  binary = "std",
  sample.names = c("Before Matching", "After Matching"),
  title = "Covariate Balance Before and After Propensity Score Matching"
)

print(love_plot)

############################################################
# 17. PROPENSITY SCORE DISTRIBUTION AFTER MATCHING
############################################################

ps_after_plot <- ggplot(
  matched_data,
  aes(
    x = distance,
    fill = treatment_label,
    y = after_stat(scaled)
  )
) +
  geom_density(alpha = 0.6) +
  scale_x_continuous(
    name = "Propensity Score",
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    name = "Scaled Density"
  ) +
  labs(
    fill = "Treatment Group",
    title = "Propensity Score Distribution After 1:1 Matching"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
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

matched_cox <- coxph(
  Surv(survival_time, MACE) ~ treatment,
  data = matched_data,
  weights = weights,
  ties = "breslow"
)

summary(matched_cox)


############################################################
# 19. PROPORTIONAL HAZARDS CHECK IN MATCHED COHORT
############################################################

matched_ph_test <- cox.zph(matched_cox)

print(matched_ph_test)

plot(matched_ph_test)


############################################################
# 20. MODEL RESULTS
############################################################

# Extract hazard ratios and 95% confidence intervals

unadjusted_results <- data.frame(
  HR = exp(coef(unadjusted_cox)),
  exp(confint(unadjusted_cox))
)

adjusted_results <- data.frame(
  HR = exp(coef(adjusted_cox)),
  exp(confint(adjusted_cox))
)

matched_results <- data.frame(
  HR = exp(coef(matched_cox)),
  exp(confint(matched_cox))
)


# Display model results

unadjusted_results

adjusted_results

matched_results


############################################################
# 21. SESSION INFORMATION
############################################################

# Record R and package versions to support reproducibility

sessionInfo()
```
