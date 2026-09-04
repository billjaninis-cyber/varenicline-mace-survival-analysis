# Varenicline vs. Bupropion: Survival Analysis of Major Adverse Cardiac Events

R-based observational healthcare outcomes analysis comparing **varenicline and bupropion** with respect to time to **major adverse cardiac events (MACE)**.

The project applies survival-analysis and propensity-score methods to evaluate treatment-associated differences in MACE risk while addressing measured baseline confounding.

## Project Overview

This analysis evaluates whether MACE-free survival differs between patients receiving varenicline versus bupropion.

The analytical workflow includes:

* Cohort and missing-data assessment
* Descriptive baseline characteristics
* Crude MACE incidence rates
* Kaplan-Meier survival analysis
* Log-rank testing
* Unadjusted Cox proportional hazards regression
* Multivariable-adjusted Cox regression
* Proportional hazards diagnostics
* Propensity score estimation
* 1:1 nearest-neighbor propensity score matching
* Covariate balance assessment using standardized mean differences
* Love plot visualization
* Propensity score overlap assessment
* Cluster-robust Cox regression in the matched cohort

## Study Design

The analysis treats:

* **Bupropion** as the reference treatment
* **Varenicline** as the exposure of interest
* **MACE** as the event outcome
* **Follow-up time** as the time-to-event measure

A complete-case analytic cohort is defined before modeling so that the primary survival and propensity-score analyses use a consistent population.

Sparse race categories are assessed before regression modeling. Original race categories are retained for descriptive reporting, while sparse categories may be collapsed for regression and propensity-score estimation to improve model stability.

## Survival Analysis

### Kaplan-Meier Analysis

Kaplan-Meier curves are used to estimate MACE-free survival over follow-up for the two treatment groups.

The analysis includes:

* Survival curves by treatment
* 95% confidence intervals
* Risk tables
* Log-rank comparison
* Follow-up visualization through 365 days

### Cox Proportional Hazards Models

Two conventional Cox models are estimated:

1. **Unadjusted Cox model**

   * Estimates the crude hazard ratio for varenicline versus bupropion.

2. **Multivariable-adjusted Cox model**

   * Adjusts for measured baseline characteristics including demographics, comorbidity measures, cardiovascular history, medication use, BMI, systolic blood pressure, smoking exposure, and frailty.

The proportional hazards assumption is evaluated using **Schoenfeld residuals**.

## Propensity Score Analysis

A logistic regression model estimates each patient's probability of receiving varenicline based on observed baseline covariates.

Covariates include:

* Age
* Sex
* Race
* Charlson Comorbidity Index
* Frailty
* Type 2 diabetes
* Coronary revascularization
* Myocardial infarction
* Stable angina
* Statin use
* BMI
* Systolic blood pressure
* Pack-year smoking exposure

The estimated propensity scores are used for **1:1 nearest-neighbor matching without replacement**.

The matching specification targets the **average treatment effect among the treated (ATT)** and uses a caliper to restrict matches to patients with sufficiently similar propensity scores.

## Covariate Balance Assessment

Balance before and after matching is evaluated using **standardized mean differences (SMDs)**.

Diagnostics include:

* Numerical balance statistics
* Pre- and post-matching SMD comparison
* Love plot visualization
* Propensity score density plots before matching
* Propensity score density plots after matching

An absolute SMD below approximately **0.10** is used as a practical benchmark for acceptable balance.

## Post-Matching Survival Analysis

A Cox proportional hazards model is fit in the propensity-score-matched cohort.

Because observations are organized into matched pairs, the model uses the MatchIt subclass identifier to obtain **cluster-robust standard errors**, accounting for dependence within matched pairs.

The resulting hazard ratio estimates the relative hazard of MACE for varenicline versus bupropion within the matched study population.

## R Packages

The analysis uses:

* `survival`
* `survminer`
* `MatchIt`
* `cobalt`
* `tableone`
* `epiDisplay`
* `ggplot2`

## Skills Demonstrated

This project demonstrates applied skills in:

* R programming
* Epidemiologic study design
* Real-world evidence analysis
* Healthcare outcomes research
* Time-to-event analysis
* Kaplan-Meier estimation
* Cox proportional hazards regression
* Confounding adjustment
* Propensity score modeling
* Propensity score matching
* Covariate balance diagnostics
* Standardized mean differences
* Love plots
* Survival-model diagnostics
* Robust inference for matched observational data
* Reproducible statistical programming

## Reproducibility

The analysis script records the R session and package versions using `sessionInfo()` to support reproducibility.

## Data Availability

The original dataset is not included in this repository due to academic and data-use restrictions.

The repository is intended to demonstrate the statistical analysis workflow, epidemiologic methodology, and R programming approach rather than distribute the underlying patient-level data.

## Disclaimer

This project is intended for educational and portfolio purposes. Results should not be interpreted as clinical guidance or as establishing a causal treatment effect without consideration of the assumptions and limitations inherent to observational data.
