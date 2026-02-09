# ============================================================================
# Factorial Simulation Study: SMARTS vs Baseline Approach (Parallel Version)
# ============================================================================
#
# Design:
# - 8 scenarios: 2 intervals × 2 confounder × 2 hazard
# - 10 data sets per scenario
# - 1000 random assignments per data set (parallelized)
#
# Methods: Naive, Adj, PSM, IPTW, SMR
# - SMR (Standardized Morbidity Ratio) weighting targets ATT among switchers
#
# Output: Data frame with mean, SD, bias for each data set
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

run_single_dataset <- function(data, n_assignments = 1000, n_cores = 12) {

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
    bias = c(hr_b_naive - 0.5, hr_b_adj - 0.5, hr_b_psm - 0.5, hr_b_iptw - 0.5, hr_b_smr - 0.5,
             smarts_summary$mean_hr - 0.5)
  )

  return(results)
}

# ============================================================================
# Run Single Scenario (10 data sets)
# ============================================================================

run_scenario <- function(
  switch_start, switch_end,
  gap_baseline, gap_peak, gap_end,
  shape_continuer, shape_switcher_pre, shape_switcher_post,
  n_pairs = 1000,
  n_datasets = 10,
  n_assignments = 1000,
  n_cores = 12,
  base_seed = 123
) {

  all_results <- list()

  for (ds in 1:n_datasets) {
    cat("    Data set", ds, "/", n_datasets, "\n")

    set.seed(base_seed + ds * 1000)

    # Simulate data
    data <- simulate_survival_data_confounder(
      n_pairs = n_pairs,
      beta_treatment = log(0.5),
      beta_confounder = log(1.3),
      lambda_0 = 10,
      shape_continuer = shape_continuer,
      shape_switcher_pre = shape_switcher_pre,
      shape_switcher_post = shape_switcher_post,
      T_min = 2,
      T_max = 6,
      switch_start = switch_start,
      switch_end = switch_end,
      confounder_interval = 0.5,
      confounder_baseline_mean = 2.5,
      confounder_gap_baseline = gap_baseline,
      confounder_gap_peak = gap_peak,
      confounder_gap_end = gap_end,
      confounder_sd = 0.8
    )

    # Run analysis
    results <- run_single_dataset(data, n_assignments, n_cores)
    results$data_set <- ds

    all_results[[ds]] <- results
  }

  return(do.call(rbind, all_results))
}

# ============================================================================
# Run Full Factorial Simulation
# ============================================================================

run_factorial_parallel <- function(
  n_pairs = 5000,
  n_datasets = 10,
  n_assignments = 10,
  n_cores = 12
) {

  # Define factor levels (removed 0.35-0.65)
  switch_intervals <- list(
    c(0.25, 0.75),
    c(0.45, 0.55)
  )

  confounder_trajectories <- list(
    variable = c(gap_baseline = 0.5, gap_peak = 1.6, gap_end = 0.8),  # Match KM curves
    homogeneous = c(gap_baseline = 0.9, gap_peak = 1.1, gap_end = 0.9)
  )

  # Hazard configurations with differentiated shapes
  hazard_configs <- list(
    realistic = c(shape_continuer = 0.9, shape_switcher_pre = 1.4, shape_switcher_post = 0.7),
    constant = c(shape_continuer = 1.0, shape_switcher_pre = 1.0, shape_switcher_post = 1.0)
  )

  all_results <- list()
  scenario_num <- 0
  total_scenarios <- length(switch_intervals) * length(confounder_trajectories) * length(hazard_configs)

  start_time <- Sys.time()

  for (sw in switch_intervals) {
    for (conf_name in names(confounder_trajectories)) {
      for (hz_name in names(hazard_configs)) {
        scenario_num <- scenario_num + 1

        conf <- confounder_trajectories[[conf_name]]
        hz <- hazard_configs[[hz_name]]

        cat("\n========================================\n")
        cat("Scenario", scenario_num, "/", total_scenarios, "\n")
        cat("Switch:", sw[1], "-", sw[2], "\n")
        cat("Confounder:", conf_name, "\n")
        cat("Hazard:", hz_name, "\n")
        cat("  - Continuer shape:", hz["shape_continuer"], "\n")
        cat("  - Switcher pre-switch shape:", hz["shape_switcher_pre"], "\n")
        cat("  - Switcher post-switch shape:", hz["shape_switcher_post"], "\n")
        cat("========================================\n")

        result <- run_scenario(
          switch_start = sw[1],
          switch_end = sw[2],
          gap_baseline = conf["gap_baseline"],
          gap_peak = conf["gap_peak"],
          gap_end = conf["gap_end"],
          shape_continuer = hz["shape_continuer"],
          shape_switcher_pre = hz["shape_switcher_pre"],
          shape_switcher_post = hz["shape_switcher_post"],
          n_pairs = n_pairs,
          n_datasets = n_datasets,
          n_assignments = n_assignments,
          n_cores = n_cores
        )

        result$switch_interval <- paste0(sw[1], "-", sw[2])
        result$confounder_type <- conf_name
        result$hazard_type <- hz_name

        all_results[[scenario_num]] <- result

        # Progress update
        elapsed <- difftime(Sys.time(), start_time, units = "mins")
        cat("  Elapsed time:", round(elapsed, 1), "minutes\n")
      }
    }
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
    group_by(switch_interval, confounder_type, hazard_type, situation, method) %>%
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
  cat("  - Scenarios: 8\n")
  cat("  - Data sets per scenario: 10\n")
  cat("  - Random assignments per data set: 1000\n")
  cat("  - Cores: 12\n")
  cat("  - Sample size: 10,000 patients\n")
  cat("========================================\n\n")

  # Run simulation
  results <- run_factorial_parallel(
    n_pairs = 1000,
    n_datasets = 10,
    n_assignments = 10,
    n_cores = 12
  )

  # Save raw results to output folder (unique filename for this design)
  saveRDS(results, "../output/factorial_results_raw_vary_switch_interval_confounder_gap.rds")
  write.csv(results, "../output/factorial_results_raw_vary_switch_interval_confounder_gap.csv", row.names = FALSE)

  # Summarize
  summary <- summarize_results(results)

  # Save summary to output folder
  saveRDS(summary, "../output/factorial_results_summary_vary_switch_interval_confounder_gap.rds")
  write.csv(summary, "../output/factorial_results_summary_vary_switch_interval_confounder_gap.csv", row.names = FALSE)

  # Display summary
  cat("\n========================================\n")
  cat("SUMMARY RESULTS\n")
  cat("True HR = 0.500\n")
  cat("========================================\n\n")

  print(summary, n = 100)

  # IPTW comparison
  cat("\n========================================\n")
  cat("IPTW COMPARISON (Baseline vs SMARTS)\n")
  cat("========================================\n\n")

  iptw_summary <- summary %>%
    filter(method == "IPTW") %>%
    select(switch_interval, confounder_type, hazard_type, situation, mean_hr_avg, bias_avg)

  print(iptw_summary, n = 100)
}
