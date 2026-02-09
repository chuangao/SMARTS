# ============================================================================
# Factorial Simulation Study: SMARTS vs Baseline Approach (Parallel Version)
# ============================================================================
#
# Design:
# - 8 clinical scenarios representing different reasons for treatment switching
# - Each scenario has its own true treatment effect (HR)
# - 10 data sets per scenario
# - 1000 random assignments per data set (parallelized)
#
# The 8 Clinical Scenarios (Proportional Hazards Version):
# ──────────────────────────────────────────────────────────────────────────────
# | # | Shape | Pre HR | Post HR | Tx HR | Clinical Story                      |
# |---|-------|--------|---------|-------|-------------------------------------|
# | 1 |  0.9  |  1.5   |   1.3   | 0.87  | Improving, sick switch, tx helps    |
# | 2 |  0.9  |  1.5   |   0.5   | 0.33  | Improving, sick switch, tx helps lot|
# | 3 |  0.9  |  0.67  |   1.3   | 1.94  | Improving, healthy switch, regret   |
# | 4 |  0.9  |  0.67  |   0.5   | 0.75  | Improving, healthy switch, stays ok |
# | 5 |  1.2  |  1.5   |   1.3   | 0.87  | Worsening, sick switch, tx helps    |
# | 6 |  1.2  |  1.5   |   0.5   | 0.33  | Worsening, sick switch, tx helps lot|
# | 7 |  1.2  |  0.67  |   1.3   | 1.94  | Worsening, healthy switch, regret   |
# | 8 |  1.2  |  0.67  |   0.5   | 0.75  | Worsening, healthy switch, stays ok |
# ──────────────────────────────────────────────────────────────────────────────
#
# Key Features:
# - Proportional Hazards maintained (same shape for all groups)
# - MVN-correlated confounders: Corr(t1,t2) = rho^|t1-t2|
# - gap_at_switch is the extreme point (max for sicker, min for healthier)
# - Treatment effect varies by scenario
#
# MVN Correlation Structure:
# - rho = 0: Independent confounders (original behavior)
# - rho = 0.8 (default): Strong temporal correlation
# - rho → 1: Nearly identical confounders across all times
#
# This makes SMARTS more effective because:
# - SMARTS assigns pseudo-switch times close to actual switch times
# - With MVN correlation, confounder at pseudo-switch ≈ confounder at switch
# - Baseline uses t=0 confounder, which is more different from switch time
#
# Methods: Naive, Adj, PSM, IPTW, SMR
# - SMR (Standardized Morbidity Ratio) weighting targets ATT among switchers
#
# Output: Data frame with mean, SD, bias for each scenario and data set
# ============================================================================

library(survival)
library(MatchIt)
library(dplyr)
library(SMARTS)
library(parallel)

source("simulate_survival_confounder.R")
source("validate_smarts_method.R")

# ============================================================================
# Single Random Assignment Function (for parallel execution)
# ============================================================================

run_single_assignment <- function(seed, switchers, continuers, confounder_interval = 0.5) {

  # Prepare for SMARTS
  switchers_smarts <- switchers
  switchers_smarts$swi_yrs <- switchers_smarts$switch_time
  switchers_smarts$fup_yrs <- switchers_smarts$T_max

  continuers_smarts <- continuers
  continuers_smarts$swi_yrs <- NA
  continuers_smarts$fup_yrs <- continuers_smarts$T_max

  smarts_input <- list(cont = continuers_smarts, swi = switchers_smarts)

  # Run random_assign
  smarts_result <- tryCatch({
    random_assign(smarts_input, nbin = 10, seed = seed,
                  swi_time = "swi_yrs", cens_time = "fup_yrs")
  }, error = function(e) { return(NULL) })

  if (is.null(smarts_result)) {
    return(c(naive = NA, adj = NA, psm = NA, iptw = NA, smr = NA))
  }

  cont_assigned <- smarts_result$assigned$cont
  swi_assigned <- smarts_result$assigned$swi

  if (nrow(cont_assigned) == 0 || nrow(swi_assigned) == 0) {
    return(c(naive = NA, adj = NA, psm = NA, iptw = NA, smr = NA))
  }

  cont_assigned$new_switch_time <- cont_assigned$swi_yrs
  swi_assigned$new_switch_time <- swi_assigned$switch_time

  # Rederive events
  cont_assigned <- rederive_events(cont_assigned, "new_switch_time")
  swi_assigned <- rederive_events(swi_assigned, "new_switch_time")

  # Lookup confounder
  cont_assigned <- lookup_confounder_at_time(cont_assigned, "new_switch_time", confounder_interval)
  swi_assigned <- lookup_confounder_at_time(swi_assigned, "new_switch_time", confounder_interval)

  # Combine
  cols_keep <- c("id", "cohort", "new_post_event", "new_post_event_time", "new_confounder_at_switch")

  smarts_data <- tryCatch({
    rbind(cont_assigned[, cols_keep], swi_assigned[, cols_keep])
  }, error = function(e) { return(NULL) })

  if (is.null(smarts_data)) {
    return(c(naive = NA, adj = NA, psm = NA, iptw = NA, smr = NA))
  }

  smarts_data$treated <- as.numeric(smarts_data$cohort == "switcher")
  smarts_data <- smarts_data[smarts_data$new_post_event_time > 0, ]

  if (nrow(smarts_data) < 50) {
    return(c(naive = NA, adj = NA, psm = NA, iptw = NA, smr = NA))
  }

  # Naive
  hr_naive <- tryCatch({
    cox <- coxph(Surv(new_post_event_time, new_post_event) ~ treated, data = smarts_data)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  # Adjusted
  hr_adj <- tryCatch({
    cox <- coxph(Surv(new_post_event_time, new_post_event) ~ treated + new_confounder_at_switch, data = smarts_data)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  # PSM
  hr_psm <- tryCatch({
    match_out <- matchit(treated ~ new_confounder_at_switch, data = smarts_data, method = "nearest", caliper = 0.1)
    data_matched <- match.data(match_out)
    if (nrow(data_matched) > 20 && length(unique(data_matched$subclass)) > 1) {
      cox <- coxph(Surv(new_post_event_time, new_post_event) ~ treated + strata(subclass), data = data_matched)
      exp(coef(cox)["treated"])
    } else {
      NA
    }
  }, error = function(e) { NA })

  # IPTW (ATE - Average Treatment Effect)
  hr_iptw <- tryCatch({
    ps <- glm(treated ~ new_confounder_at_switch, data = smarts_data, family = binomial)
    smarts_data$ps <- predict(ps, type = "response")
    smarts_data$weight <- ifelse(smarts_data$treated == 1, 1/smarts_data$ps, 1/(1-smarts_data$ps))
    cox <- coxph(Surv(new_post_event_time, new_post_event) ~ treated, data = smarts_data, weights = weight)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  # SMR (ATT - Average Treatment Effect in the Treated/Switchers)
  # Switchers get weight=1, Continuers get weight=PS/(1-PS)
  hr_smr <- tryCatch({
    if (is.null(smarts_data$ps)) {
      ps <- glm(treated ~ new_confounder_at_switch, data = smarts_data, family = binomial)
      smarts_data$ps <- predict(ps, type = "response")
    }
    smarts_data$smr_weight <- ifelse(smarts_data$treated == 1, 1, smarts_data$ps / (1 - smarts_data$ps))
    cox <- coxph(Surv(new_post_event_time, new_post_event) ~ treated, data = smarts_data, weights = smr_weight)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  return(c(naive = hr_naive, adj = hr_adj, psm = hr_psm, iptw = hr_iptw, smr = hr_smr))
}

# ============================================================================
# Run Single Data Set (with 1000 parallel random assignments)
# ============================================================================
#
# Parameters:
#   data         - Simulated data from simulate_survival_data_confounder()
#   true_hr      - The true hazard ratio for this scenario (for bias calculation)
#   n_assignments - Number of SMARTS random assignments to run
#   n_cores      - Number of cores for parallel execution
#
# ============================================================================

run_single_dataset <- function(data, true_hr = 0.5, n_assignments = 1000, n_cores = 12) {

  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  # =========================================
  # BASELINE APPROACH (no random assignment)
  # =========================================

  # Switchers: from switch time
  baseline_swi <- data.frame(
    id = switchers$id,
    treated = 1,
    confounder = switchers$confounder_at_switch,
    event = switchers$post_event,
    event_time = switchers$post_event_time - switchers$switch_time
  )

  # Continuers: from t=0
  baseline_cont <- data.frame(
    id = continuers$id,
    treated = 0,
    confounder = continuers$confounder_0,
    event = pmax(continuers$pre_event, continuers$post_event),
    event_time = ifelse(continuers$pre_event == 1,
                        continuers$pre_event_time,
                        continuers$post_event_time)
  )

  baseline_data <- rbind(baseline_swi, baseline_cont)
  baseline_data <- baseline_data[baseline_data$event_time > 0, ]

  # Baseline methods
  hr_b_naive <- tryCatch({
    cox <- coxph(Surv(event_time, event) ~ treated, data = baseline_data)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  hr_b_adj <- tryCatch({
    cox <- coxph(Surv(event_time, event) ~ treated + confounder, data = baseline_data)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  hr_b_psm <- tryCatch({
    match_out <- matchit(treated ~ confounder, data = baseline_data, method = "nearest", caliper = 0.1)
    data_matched <- match.data(match_out)
    if (nrow(data_matched) > 20 && length(unique(data_matched$subclass)) > 1) {
      cox <- coxph(Surv(event_time, event) ~ treated + strata(subclass), data = data_matched)
      exp(coef(cox)["treated"])
    } else {
      NA
    }
  }, error = function(e) { NA })

  hr_b_iptw <- tryCatch({
    ps <- glm(treated ~ confounder, data = baseline_data, family = binomial)
    baseline_data$ps <- predict(ps, type = "response")
    baseline_data$weight <- ifelse(baseline_data$treated == 1, 1/baseline_data$ps, 1/(1-baseline_data$ps))
    cox <- coxph(Surv(event_time, event) ~ treated, data = baseline_data, weights = weight)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  hr_b_smr <- tryCatch({
    if (is.null(baseline_data$ps)) {
      ps <- glm(treated ~ confounder, data = baseline_data, family = binomial)
      baseline_data$ps <- predict(ps, type = "response")
    }
    baseline_data$smr_weight <- ifelse(baseline_data$treated == 1, 1, baseline_data$ps / (1 - baseline_data$ps))
    cox <- coxph(Surv(event_time, event) ~ treated, data = baseline_data, weights = smr_weight)
    exp(coef(cox)["treated"])
  }, error = function(e) { NA })

  # =========================================
  # SMARTS APPROACH (1000 parallel random assignments)
  # =========================================

  seeds <- 1:n_assignments

  smarts_results <- mclapply(seeds, function(s) {
    run_single_assignment(s, switchers, continuers, 0.5)
  }, mc.cores = n_cores)

  smarts_df <- do.call(rbind, smarts_results)
  smarts_df <- as.data.frame(smarts_df)

  # Calculate mean and SD for each method
  smarts_summary <- data.frame(
    method = c("Naive", "Adj", "PSM", "IPTW", "SMR"),
    mean_hr = c(mean(smarts_df$naive, na.rm = TRUE),
                mean(smarts_df$adj, na.rm = TRUE),
                mean(smarts_df$psm, na.rm = TRUE),
                mean(smarts_df$iptw, na.rm = TRUE),
                mean(smarts_df$smr, na.rm = TRUE)),
    sd_hr = c(sd(smarts_df$naive, na.rm = TRUE),
              sd(smarts_df$adj, na.rm = TRUE),
              sd(smarts_df$psm, na.rm = TRUE),
              sd(smarts_df$iptw, na.rm = TRUE),
              sd(smarts_df$smr, na.rm = TRUE))
  )

  # Combine results
  results <- data.frame(
    situation = c(rep("Baseline", 5), rep("SMARTS", 5)),
    method = rep(c("Naive", "Adj", "PSM", "IPTW", "SMR"), 2),
    mean_hr = c(hr_b_naive, hr_b_adj, hr_b_psm, hr_b_iptw, hr_b_smr,
                smarts_summary$mean_hr),
    sd_hr = c(NA, NA, NA, NA, NA, smarts_summary$sd_hr),
    true_hr = true_hr,
    bias = c(hr_b_naive - true_hr, hr_b_adj - true_hr, hr_b_psm - true_hr,
             hr_b_iptw - true_hr, hr_b_smr - true_hr,
             smarts_summary$mean_hr - true_hr)
  )

  return(results)
}

# ============================================================================
# Run Single Scenario (10 data sets)
# ============================================================================

# ============================================================================
# Run Single Scenario (using scenario configuration)
# ============================================================================
#
# This function runs a single scenario configuration across multiple datasets.
# It uses the get_scenario_config() function from simulate_survival_confounder.R
# to get pre-defined clinical scenario parameters.
#
# The scenario_config now includes beta_treatment (true treatment effect),
# which varies by scenario to match the clinical story.
#
# ============================================================================
run_scenario <- function(
  scenario_config,           # Configuration from get_scenario_config()
  beta_confounder = log(1.3),
  rho = 0.8,                 # MVN temporal correlation parameter
  n_pairs = 5000,
  n_datasets = 10,
  n_assignments = 1000,
  n_cores = 12,
  base_seed = 123
) {

  all_results <- list()

  # Get true HR from scenario config
  true_hr <- exp(scenario_config$beta_treatment)

  for (ds in 1:n_datasets) {
    cat("    Data set", ds, "/", n_datasets, "\n")

    set.seed(base_seed + ds * 1000)

    # Simulate data using scenario configuration
    # NOTE: beta_treatment comes from scenario_config, not hardcoded
    # NOTE: rho controls MVN temporal correlation of confounders
    data <- simulate_survival_data_confounder(
      n_pairs = n_pairs,
      beta_treatment = scenario_config$beta_treatment,  # From scenario config
      beta_confounder = beta_confounder,
      lambda_0 = 10,
      shape_continuer = scenario_config$shape_continuer,
      shape_switcher_pre = scenario_config$shape_switcher_pre,
      shape_switcher_post = scenario_config$shape_switcher_post,
      T_min = 2,
      T_max = 6,
      switch_start = 0.25,
      switch_end = 0.75,
      confounder_interval = 0.5,
      confounder_baseline_mean = 2.5,
      confounder_change_magnitude = scenario_config$confounder_change_magnitude,
      confounder_gap_baseline = scenario_config$gap_baseline,
      confounder_gap_at_switch = scenario_config$gap_at_switch,
      confounder_gap_end = scenario_config$gap_end,
      confounder_sd = 0.8,
      rho = rho
    )

    # Run analysis with scenario-specific true HR
    results <- run_single_dataset(data, true_hr, n_assignments, n_cores)
    results$data_set <- ds

    all_results[[ds]] <- results
  }

  return(do.call(rbind, all_results))
}

# ============================================================================
# Run Full Factorial Simulation
# ============================================================================

# ============================================================================
# Run All 8 Clinical Scenarios
# ============================================================================
#
# This function runs all 8 pre-defined clinical scenarios.
# Each scenario now has its own true treatment effect (beta_treatment/HR).
#
# | # | True HR | Continuer | Swi Pre | Swi Post | Clinical Story              |
# |---|---------|-----------|---------|----------|------------------------------|
# | 1 |   0.5   | ↓ (0.9)   | ↑ (1.2) | ≈ (0.95) | Classic medical switch       |
# | 2 |   1.2   | ↓ (0.9)   | ↑ (1.2) | ↑ (1.15) | Failed medical switch        |
# | 3 |   1.0   | ↓ (0.9)   | ↓ (0.9) | ↓ (0.9)  | Financial switch (healthy)   |
# | 4 |   1.5   | ↓ (0.9)   | ↓ (0.85)| ↑ (1.1)  | Regret switch                |
# | 5 |   0.5   | ↑ (1.2)   | ↑ (1.4) | ↑ (1.1)  | Everyone declining, helps    |
# | 6 |   1.0   | ↑ (1.2)   | ↑ (1.3) | ↑ (1.25) | No good options              |
# | 7 |   0.5   | ↑ (1.2)   | ≈ (1.0) | ≈ (1.0)  | Healthy responders switch    |
# | 8 |   1.3   | ↑ (1.2)   | ≈ (1.0) | ↑ (1.25) | Healthy switch, regret       |
#
# ============================================================================
run_factorial_parallel <- function(
  n_pairs = 5000,
  n_datasets = 10,
  n_assignments = 1000,
  n_cores = 12,
  rho = 0.8,                # MVN temporal correlation (0=independent, 1=identical)
  scenarios_to_run = 1:8    # Which scenarios to run (default: all 8)
) {

  all_results <- list()
  total_scenarios <- length(scenarios_to_run)
  scenario_idx <- 0

  start_time <- Sys.time()

  for (scenario_num in scenarios_to_run) {
    scenario_idx <- scenario_idx + 1

    # Get pre-defined scenario configuration
    config <- get_scenario_config(scenario_num)
    true_hr <- exp(config$beta_treatment)

    cat("\n========================================\n")
    cat("Scenario", scenario_num, "(", scenario_idx, "/", total_scenarios, ")\n")
    cat("Name:", config$name, "\n")
    cat("True HR:", round(true_hr, 2),
        ifelse(true_hr < 1, "(beneficial)",
               ifelse(true_hr > 1, "(harmful)", "(neutral)")), "\n")
    cat("----------------------------------------\n")
    cat("Continuer shape:", config$shape_continuer,
        ifelse(config$shape_continuer < 1, "(improving)", "(worsening)"), "\n")
    cat("Switcher pre-switch shape:", config$shape_switcher_pre,
        ifelse(config$shape_switcher_pre < 1, "(improving)", "(worsening)"), "\n")
    cat("Switcher post-switch shape:", config$shape_switcher_post,
        ifelse(config$shape_switcher_post < 1, "(improving)", "(worsening)"), "\n")
    cat("----------------------------------------\n")
    cat("Gap trajectory:", config$gap_baseline, "→",
        config$gap_at_switch, "→", config$gap_end, "\n")
    cat("========================================\n")

    result <- run_scenario(
      scenario_config = config,
      beta_confounder = log(1.3),
      rho = rho,
      n_pairs = n_pairs,
      n_datasets = n_datasets,
      n_assignments = n_assignments,
      n_cores = n_cores
    )

    # Add scenario identifiers
    result$scenario_num <- scenario_num
    result$scenario_name <- config$name
    result$shape_continuer <- config$shape_continuer
    result$shape_switcher_pre <- config$shape_switcher_pre
    result$shape_switcher_post <- config$shape_switcher_post

    all_results[[scenario_idx]] <- result

    # Progress update
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat("  Elapsed time:", round(elapsed, 1), "minutes\n")
  }

  final_results <- do.call(rbind, all_results)
  rownames(final_results) <- NULL

  total_time <- difftime(Sys.time(), start_time, units = "mins")
  cat("\n========================================\n")
  cat("COMPLETED in", round(total_time, 1), "minutes\n")
  cat("========================================\n")

  return(final_results)
}

# ============================================================================
# Summarize Results
# ============================================================================

summarize_results <- function(results) {

  # Aggregate over data sets
  summary_df <- results %>%
    group_by(scenario_num, scenario_name, true_hr, situation, method) %>%
    summarise(
      mean_hr_avg = mean(mean_hr, na.rm = TRUE),
      mean_hr_sd = sd(mean_hr, na.rm = TRUE),
      sd_hr_avg = mean(sd_hr, na.rm = TRUE),
      bias_avg = mean(bias, na.rm = TRUE),
      bias_sd = sd(bias, na.rm = TRUE),
      .groups = "drop"
    )

  return(summary_df)
}

# ============================================================================
# Main Execution
# ============================================================================

if (interactive()) {

  cat("========================================\n")
  cat("FACTORIAL SIMULATION STUDY (Parallel)\n")
  cat("========================================\n")
  cat("Settings:\n")
  cat("  - Scenarios: 8 clinical scenarios\n")
  cat("  - Each scenario has its own true HR\n")
  cat("  - Data sets per scenario: 10\n")
  cat("  - Random assignments per data set: 1000\n")
  cat("  - Cores: 12\n")
  cat("  - Sample size: 10,000 patients\n")
  cat("========================================\n\n")

  # Print scenario summary
  print_scenario_summary()

  # Run simulation
  results <- run_factorial_parallel(
    n_pairs = 5000,
    n_datasets = 10,
    n_assignments = 1000,
    n_cores = 12
  )

  # Save raw results to output folder
  saveRDS(results, "../output/factorial_results_raw.rds")
  write.csv(results, "../output/factorial_results_raw.csv", row.names = FALSE)

  # Summarize
  summary <- summarize_results(results)

  # Save summary to output folder
  saveRDS(summary, "../output/factorial_results_summary.rds")
  write.csv(summary, "../output/factorial_results_summary.csv", row.names = FALSE)

  # Display summary
  cat("\n========================================\n")
  cat("SUMMARY RESULTS\n")
  cat("(True HR varies by scenario)\n")
  cat("========================================\n\n")

  print(summary, n = 100)

  # IPTW comparison
  cat("\n========================================\n")
  cat("IPTW COMPARISON (Baseline vs SMARTS)\n")
  cat("========================================\n\n")

  iptw_summary <- summary %>%
    filter(method == "IPTW") %>%
    select(scenario_num, scenario_name, true_hr, situation, mean_hr_avg, bias_avg)

  print(iptw_summary, n = 100)
}
