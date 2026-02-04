# ============================================================================
# Factorial Simulation: Piecewise Constant Hazard
# WITH TIME-VARYING CONFOUNDING
# ============================================================================
# Two clinical scenarios:
# 1. Sicker switchers, treatment helps (HR=0.7)
#    Gap: +0.5 -> +1.5 (peak) -> +0.8
# 2. Healthier switchers, treatment harms (HR=1.5)
#    Gap: -0.3 -> -1.5 (minimum) -> -0.5
#
# Compares Baseline vs SMARTS across 3 hazard trends
# Methods: Naive, Adjusted, PSM, IPTW, SMR
# - IPTW estimates ATE (Average Treatment Effect)
# - SMR estimates ATT (Average Treatment Effect in Treated/Switchers)
# ============================================================================

library(survival)
library(MASS)
library(SMARTS)
library(MatchIt)
library(dplyr)
library(parallel)

source("simulate_piecewise_hazard.R")

# ============================================================================
# Parameters
# ============================================================================

n_simulations <- 30
n_pairs <- 2000
n_cores <- 4

# Define clinical scenarios
clinical_scenarios <- list(
  list(
    name = "Sicker_Helps",
    description = "Sicker switchers, treatment helps",
    true_hr = 0.7,
    gap_baseline = 0.5,     # Switchers slightly sicker at baseline
    gap_at_switch = 1.5,    # Gap peaks at switch (switchers much sicker)
    gap_end = 0.8           # Gap shrinks after treatment helps
  ),
  list(
    name = "Healthier_Harms",
    description = "Healthier switchers, treatment harms",
    true_hr = 1.5,
    gap_baseline = -0.3,    # Switchers slightly healthier at baseline
    gap_at_switch = -1.5,   # Gap at minimum (switchers much healthier)
    gap_end = -0.5          # Gap shrinks (less negative) as treatment harms
  )
)

hazard_trends <- list(
  list(name = "Constant", trend = "constant", strength = 0),
  list(name = "Increasing", trend = "increasing", strength = 0.5),
  list(name = "Decreasing", trend = "decreasing", strength = 0.25)
)

# ============================================================================
# Analysis Functions
# ============================================================================

run_baseline_analysis <- function(data) {
  # Baseline: no pseudo-switch assignment
  # Switchers: post-switch events from switch time
  # Continuers: all events from time 0

  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  first_event_cont <- ifelse(continuers$has_pre_event,
                             continuers$first_pre_event_time,
                             ifelse(continuers$has_post_event,
                                    continuers$first_post_event_time, NA))

  baseline_data <- data.frame(
    time = c(
      ifelse(switchers$has_post_event,
             switchers$first_post_event_time - switchers$switch_time,
             switchers$T_max - switchers$switch_time),
      ifelse(continuers$n_total_events > 0,
             first_event_cont,
             continuers$T_max)
    ),
    event = c(
      as.numeric(switchers$has_post_event),
      as.numeric(continuers$n_total_events > 0)
    ),
    treated = c(rep(1, nrow(switchers)), rep(0, nrow(continuers))),
    confounder = c(switchers$confounder_at_switch, continuers$conf_t0)
  )

  baseline_data <- baseline_data[baseline_data$time > 0, ]

  results <- list()

  # Naive
  tryCatch({
    cox <- coxph(Surv(time, event) ~ treated, data = baseline_data)
    results$Naive <- exp(coef(cox)["treated"])
  }, error = function(e) { results$Naive <- NA })

  # Adjusted
  tryCatch({
    cox <- coxph(Surv(time, event) ~ treated + confounder, data = baseline_data)
    results$Adj <- exp(coef(cox)["treated"])
  }, error = function(e) { results$Adj <- NA })

  # PSM
  tryCatch({
    ps_model <- glm(treated ~ confounder, data = baseline_data, family = binomial)
    baseline_data$ps <- predict(ps_model, type = "response")
    match_out <- matchit(treated ~ confounder, data = baseline_data,
                         method = "nearest", caliper = 0.2)
    matched_data <- match.data(match_out)
    cox <- coxph(Surv(time, event) ~ treated, data = matched_data)
    results$PSM <- exp(coef(cox)["treated"])
  }, error = function(e) { results$PSM <- NA })

  # IPTW (ATE)
  tryCatch({
    if (is.null(baseline_data$ps)) {
      ps_model <- glm(treated ~ confounder, data = baseline_data, family = binomial)
      baseline_data$ps <- predict(ps_model, type = "response")
    }
    p_treated <- mean(baseline_data$treated)
    baseline_data$weight <- ifelse(baseline_data$treated == 1,
                                    p_treated / baseline_data$ps,
                                    (1 - p_treated) / (1 - baseline_data$ps))
    cox <- coxph(Surv(time, event) ~ treated, data = baseline_data, weights = weight)
    results$IPTW <- exp(coef(cox)["treated"])
  }, error = function(e) { results$IPTW <- NA })

  # SMR (ATT - Average Treatment Effect in Treated/Switchers)
  # Switchers get weight=1, Continuers get weight=PS/(1-PS)
  tryCatch({
    if (is.null(baseline_data$ps)) {
      ps_model <- glm(treated ~ confounder, data = baseline_data, family = binomial)
      baseline_data$ps <- predict(ps_model, type = "response")
    }
    baseline_data$smr_weight <- ifelse(baseline_data$treated == 1,
                                        1,
                                        baseline_data$ps / (1 - baseline_data$ps))
    cox <- coxph(Surv(time, event) ~ treated, data = baseline_data, weights = smr_weight)
    results$SMR <- exp(coef(cox)["treated"])
  }, error = function(e) { results$SMR <- NA })

  return(results)
}

run_smarts_analysis <- function(data) {
  # SMARTS: assign pseudo-switch times to continuers

  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  # SMARTS assignment
  switchers_smarts <- switchers
  switchers_smarts$swi_yrs <- switchers_smarts$switch_time
  switchers_smarts$fup_yrs <- switchers_smarts$T_max

  continuers_smarts <- continuers
  continuers_smarts$swi_yrs <- NA
  continuers_smarts$fup_yrs <- continuers_smarts$T_max

  smarts_result <- tryCatch({
    random_assign(list(cont = continuers_smarts, swi = switchers_smarts),
                  nbin = 10, seed = sample(1:10000, 1),
                  swi_time = "swi_yrs", cens_time = "fup_yrs")
  }, error = function(e) NULL)

  if (is.null(smarts_result)) {
    return(list(Naive = NA, Adj = NA, PSM = NA, IPTW = NA, SMR = NA))
  }

  cont_assigned <- smarts_result$assigned$cont
  swi_assigned <- smarts_result$assigned$swi

  if (nrow(cont_assigned) == 0 || nrow(swi_assigned) == 0) {
    return(list(Naive = NA, Adj = NA, PSM = NA, IPTW = NA, SMR = NA))
  }

  # Re-derive events
  cont_assigned$new_switch_time <- cont_assigned$swi_yrs
  swi_assigned$new_switch_time <- swi_assigned$switch_time

  cont_assigned <- rederive_events_recurring(cont_assigned, "new_switch_time")
  swi_assigned <- rederive_events_recurring(swi_assigned, "new_switch_time")

  # Create analysis dataset
  smarts_data <- data.frame(
    time = c(
      ifelse(swi_assigned$new_has_post_event,
             swi_assigned$new_post_event_time_from_switch,
             swi_assigned$T_max - swi_assigned$new_switch_time),
      ifelse(cont_assigned$new_has_post_event,
             cont_assigned$new_post_event_time_from_switch,
             cont_assigned$T_max - cont_assigned$new_switch_time)
    ),
    event = c(
      as.numeric(swi_assigned$new_has_post_event),
      as.numeric(cont_assigned$new_has_post_event)
    ),
    treated = c(rep(1, nrow(swi_assigned)), rep(0, nrow(cont_assigned))),
    confounder = c(swi_assigned$confounder_at_switch,
                   cont_assigned$confounder_at_switch)
  )

  smarts_data <- smarts_data[smarts_data$time > 0 & !is.na(smarts_data$time), ]

  if (nrow(smarts_data) < 100) {
    return(list(Naive = NA, Adj = NA, PSM = NA, IPTW = NA, SMR = NA))
  }

  results <- list()

  # Naive
  tryCatch({
    cox <- coxph(Surv(time, event) ~ treated, data = smarts_data)
    results$Naive <- exp(coef(cox)["treated"])
  }, error = function(e) { results$Naive <- NA })

  # Adjusted
  tryCatch({
    cox <- coxph(Surv(time, event) ~ treated + confounder, data = smarts_data)
    results$Adj <- exp(coef(cox)["treated"])
  }, error = function(e) { results$Adj <- NA })

  # PSM
  tryCatch({
    ps_model <- glm(treated ~ confounder, data = smarts_data, family = binomial)
    smarts_data$ps <- predict(ps_model, type = "response")
    match_out <- matchit(treated ~ confounder, data = smarts_data,
                         method = "nearest", caliper = 0.2)
    matched_data <- match.data(match_out)
    cox <- coxph(Surv(time, event) ~ treated, data = matched_data)
    results$PSM <- exp(coef(cox)["treated"])
  }, error = function(e) { results$PSM <- NA })

  # IPTW (ATE)
  tryCatch({
    if (is.null(smarts_data$ps)) {
      ps_model <- glm(treated ~ confounder, data = smarts_data, family = binomial)
      smarts_data$ps <- predict(ps_model, type = "response")
    }
    p_treated <- mean(smarts_data$treated)
    smarts_data$weight <- ifelse(smarts_data$treated == 1,
                                  p_treated / smarts_data$ps,
                                  (1 - p_treated) / (1 - smarts_data$ps))
    cox <- coxph(Surv(time, event) ~ treated, data = smarts_data, weights = weight)
    results$IPTW <- exp(coef(cox)["treated"])
  }, error = function(e) { results$IPTW <- NA })

  # SMR (ATT - Average Treatment Effect in Treated/Switchers)
  tryCatch({
    if (is.null(smarts_data$ps)) {
      ps_model <- glm(treated ~ confounder, data = smarts_data, family = binomial)
      smarts_data$ps <- predict(ps_model, type = "response")
    }
    smarts_data$smr_weight <- ifelse(smarts_data$treated == 1,
                                      1,
                                      smarts_data$ps / (1 - smarts_data$ps))
    cox <- coxph(Surv(time, event) ~ treated, data = smarts_data, weights = smr_weight)
    results$SMR <- exp(coef(cox)["treated"])
  }, error = function(e) { results$SMR <- NA })

  return(results)
}

# ============================================================================
# Run Simulations
# ============================================================================

cat("================================================================\n")
cat("FACTORIAL SIMULATION: Piecewise Constant Hazard\n")
cat("Two Clinical Scenarios x 3 Hazard Trends\n")
cat("================================================================\n\n")

all_results <- data.frame()

for (scenario in clinical_scenarios) {
  cat("\n================================================================\n")
  cat("Clinical Scenario:", scenario$description, "\n")
  cat("True HR:", scenario$true_hr, "\n")
  cat("Gap trajectory:", scenario$gap_baseline, "->", scenario$gap_at_switch,
      "->", scenario$gap_end, "\n")
  cat("================================================================\n\n")

  for (h in hazard_trends) {
    cat("Hazard trend:", h$name, "\n")

    for (sim in 1:n_simulations) {
      set.seed(1000 * which(sapply(hazard_trends, function(x) x$name) == h$name) +
               100 * which(sapply(clinical_scenarios, function(x) x$name) == scenario$name) + sim)

      # Simulate data with time-varying confounding
      data <- simulate_piecewise_hazard(
        n_pairs = n_pairs,
        beta_treatment = log(scenario$true_hr),
        beta_confounder = log(2.0),
        base_hazard = 0.02,
        hazard_trend = h$trend,
        trend_strength = h$strength,
        segment_length = 0.5,
        confounder_gap_baseline = scenario$gap_baseline,
        confounder_gap_at_switch = scenario$gap_at_switch,
        confounder_gap_end = scenario$gap_end,
        rho = 0.9
      )

      # Run Baseline analysis
      baseline_results <- run_baseline_analysis(data)

      # Run SMARTS analysis
      smarts_results <- run_smarts_analysis(data)

      # Store results
      for (method in c("Naive", "Adj", "PSM", "IPTW", "SMR")) {
        all_results <- rbind(all_results, data.frame(
          clinical_scenario = scenario$name,
          true_hr = scenario$true_hr,
          hazard_trend = h$name,
          simulation = sim,
          situation = "Baseline",
          method = method,
          hr = baseline_results[[method]]
        ))

        all_results <- rbind(all_results, data.frame(
          clinical_scenario = scenario$name,
          true_hr = scenario$true_hr,
          hazard_trend = h$name,
          simulation = sim,
          situation = "SMARTS",
          method = method,
          hr = smarts_results[[method]]
        ))
      }

      if (sim %% 10 == 0) cat("  Completed", sim, "of", n_simulations, "\n")
    }
  }
}

# Save raw results
saveRDS(all_results, "factorial_results_raw.rds")
cat("\nResults saved to factorial_results_raw.rds\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n================================================================\n")
cat("SUMMARY\n")
cat("================================================================\n\n")

summary_stats <- all_results %>%
  group_by(clinical_scenario, true_hr, hazard_trend, situation, method) %>%
  summarise(
    mean_hr = round(mean(hr, na.rm = TRUE), 3),
    bias = round(mean(hr, na.rm = TRUE) - first(true_hr), 3),
    sd = round(sd(hr, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_stats, n = 50)
