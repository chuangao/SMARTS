# ============================================================================
# Verify Hazard Ratios After Fix
# ============================================================================
# This script verifies that the Weibull scale parameter fix produces
# the correct hazard ratios using Cox proportional hazards models

library(survival)

# Load the confounder simulation data
source("simulate_survival_confounder.R")

set.seed(123)
data <- simulate_survival_data_confounder(
  n_pairs = 10000,
  beta_treatment = log(0.5),
  beta_confounder = log(1.3),
  lambda_0 = 10,
  shape = 1.5,
  T_min = 2,
  T_max = 6,
  switch_start = 0.25,
  switch_end = 0.75,
  confounder_interval = 0.5,
  confounder_baseline_mean = 2.5,
  confounder_gap_baseline = 0.5,
  confounder_gap_peak = 1.6,
  confounder_gap_floor = 0.8,
  confounder_sd = 0.8
)




cat("=== Testing Hazard Ratio Estimation ===\n\n")

# Expected parameters
beta_treatment_true <- log(0.5)
beta_confounder_true <- log(1.3)

cat("True parameters:\n")
cat("  Treatment effect: beta =", beta_treatment_true,
    ", HR =", exp(beta_treatment_true), "\n")
cat("  Confounder effect:  beta =", beta_confounder_true,
    ", HR =", exp(beta_confounder_true), "\n\n")

# ============================================================================
# Test 1: Pre-switching period - verify confounder effect
# ============================================================================

cat("=== Test 1: Pre-switching Confounder Effect ===\n")

# Create dataset for pre-switching events
data_pre <- data.frame(
  id = data$id,
  cohort = data$cohort,
  time = ifelse(data$pre_event == 1, data$pre_event_time, data$switch_time),
  event = data$pre_event,
  confounder = data$confounder_at_switch
)

cox_pre <- coxph(Surv(time, event) ~ cohort, data = data_pre)
cox_pre

# Fit Cox model
cox_pre <- coxph(Surv(time, event) ~ confounder, data = data_pre)
cat("\nCox model: Surv(time, event) ~ confounder\n")
print(summary(cox_pre)$coefficients)

cat("\nExpected HR for confounder:", exp(beta_confounder_true), "\n")
cat("Observed HR for confounder:", exp(coef(cox_pre)["confounder"]), "\n")
cat("Difference:", abs(exp(coef(cox_pre)["confounder"]) - exp(beta_confounder_true)), "\n")

# ============================================================================
# Test 2: Post-switching period - verify treatment + confounder effect
# ============================================================================

cat("\n\n=== Test 2: Post-switching Treatment Effect (Naive) ===\n")

# Create dataset for post-switching events (using time from switching)
data_post <- data.frame(
  id = data$id,
  cohort = data$cohort,
  time = ifelse(data$post_event == 1,
                data$post_event_time - data$switch_time,
                data$T_max - data$switch_time),
  event = data$post_event,
  confounder = data$confounder_at_switch,
  is_switcher = as.numeric(data$cohort == "switcher")
)

# Naive model (no adjustment)
cox_naive <- coxph(Surv(time, event) ~ is_switcher, data = data_post)
cat("\nNaive Cox model: Surv(time, event) ~ is_switcher\n")
print(summary(cox_naive)$coefficients)

cat("\nExpected treatment HR:", exp(beta_treatment_true), "\n")
cat("Observed naive HR:", exp(coef(cox_naive)["is_switcher"]), "\n")
cat("Note: Should be BIASED (not equal to true HR) due to confounder confounding\n")

# ============================================================================
# Test 3: Post-switching with confounder adjustment
# ============================================================================

cat("\n\n=== Test 3: Post-switching Treatment Effect (Adjusted) ===\n")

# Adjusted model
cox_adj <- coxph(Surv(time, event) ~ is_switcher + confounder, data = data_post)
cat("\nAdjusted Cox model: Surv(time, event) ~ is_switcher + confounder\n")
print(summary(cox_adj)$coefficients)

cat("\nExpected treatment HR:", exp(beta_treatment_true), "\n")
cat("Observed adjusted HR:", exp(coef(cox_adj)["is_switcher"]), "\n")
cat("Difference:", abs(exp(coef(cox_adj)["is_switcher"]) - exp(beta_treatment_true)), "\n")

cat("\nExpected confounder HR:", exp(beta_confounder_true), "\n")
cat("Observed confounder HR:", exp(coef(cox_adj)["confounder"]), "\n")
cat("Difference:", abs(exp(coef(cox_adj)["confounder"]) - exp(beta_confounder_true)), "\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n\n=== Summary ===\n")
cat("If the scale parameter fix is correct:\n")
cat("1. Pre-switching confounder HR should be close to 1.3\n")
cat("2. Naive treatment HR should be biased (NOT 0.5)\n")
cat("3. Adjusted treatment HR should recover ~0.5\n")
cat("4. Adjusted confounder HR should be close to 1.3\n")
