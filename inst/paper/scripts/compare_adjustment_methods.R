# ============================================================================
# Compare Adjustment Methods for Treatment Effect Estimation
# ============================================================================
# This function compares 4 methods for estimating treatment effects:
# 1. Naive (unadjusted)
# 2. Covariate adjustment (confounder as covariate)
# 3. Propensity Score Matching (PSM)
# 4. Inverse Probability of Treatment Weighting (IPTW)
# ============================================================================

library(survival)
library(MatchIt)

source("simulate_survival_confounder.R")

compare_adjustment_methods <- function(
  n_pairs = 1000,
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
) {

  # Simulate data
  data <- simulate_survival_data_confounder(
    n_pairs = n_pairs,
    beta_treatment = beta_treatment,
    beta_confounder = beta_confounder,
    lambda_0 = lambda_0,
    shape = shape,
    T_min = T_min,
    T_max = T_max,
    switch_start = switch_start,
    switch_end = switch_end,
    confounder_interval = confounder_interval,
    confounder_baseline_mean = confounder_baseline_mean,
    confounder_gap_baseline = confounder_gap_baseline,
    confounder_gap_peak = confounder_gap_peak,
    confounder_gap_floor = confounder_gap_floor,
    confounder_sd = confounder_sd
  )

  # Prepare post-switching dataset
  data_post <- data.frame(
    id = data$id,
    time = ifelse(data$post_event == 1,
                  data$post_event_time - data$switch_time,
                  data$T_max - data$switch_time),
    event = data$post_event,
    confounder = data$confounder_at_switch,
    treated = as.numeric(data$cohort == "switcher")
  )

  # 1. Naive (unadjusted)
  cox_naive <- coxph(Surv(time, event) ~ treated, data = data_post)
  beta_naive <- coef(cox_naive)["treated"]

  # 2. Covariate adjustment
  cox_adj <- coxph(Surv(time, event) ~ treated + confounder, data = data_post)
  beta_adj <- coef(cox_adj)["treated"]

  # 3. Propensity Score Matching (with caliper + stratified Cox)
  ps_model <- glm(treated ~ confounder, data = data_post, family = binomial)
  data_post$ps <- predict(ps_model, type = "response")

  match_out <- matchit(treated ~ confounder, data = data_post,
                       method = "nearest", caliper = 0.1)
  data_matched <- match.data(match_out)

  # Use stratified Cox to account for matched pair structure
  if (nrow(data_matched) > 0 && length(unique(data_matched$subclass)) > 1) {
    cox_psm <- coxph(Surv(time, event) ~ treated + strata(subclass),
                     data = data_matched)
    beta_psm <- coef(cox_psm)["treated"]
  } else {
    beta_psm <- NA
  }

  # 4. IPTW
  data_post$weight <- ifelse(data_post$treated == 1,
                              1 / data_post$ps,
                              1 / (1 - data_post$ps))

  cox_iptw <- coxph(Surv(time, event) ~ treated, data = data_post, weights = weight)
  beta_iptw <- coef(cox_iptw)["treated"]

  # Compile results
  results <- data.frame(
    method = c("Naive", "Covariate Adj", "PSM", "IPTW"),
    log_hr = c(beta_naive, beta_adj, beta_psm, beta_iptw),
    hr = exp(c(beta_naive, beta_adj, beta_psm, beta_iptw)),
    true_log_hr = beta_treatment,
    true_hr = exp(beta_treatment),
    bias = c(beta_naive, beta_adj, beta_psm, beta_iptw) - beta_treatment
  )

  return(results)
}

# ============================================================================
# Run Simulation Study (Repeat n times, with optional parallel computing)
# ============================================================================

run_simulation_study <- function(
  n_sim = 100,
  seed_start = 1,
  n_pairs = 1000,
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
  confounder_sd = 0.8,
  n_cores = 1,
  verbose = TRUE
) {

  # Function to run a single simulation
  run_single_sim <- function(i) {
    set.seed(seed_start + i - 1)

    sim_result <- compare_adjustment_methods(
      n_pairs = n_pairs,
      beta_treatment = beta_treatment,
      beta_confounder = beta_confounder,
      lambda_0 = lambda_0,
      shape = shape,
      T_min = T_min,
      T_max = T_max,
      switch_start = switch_start,
      switch_end = switch_end,
      confounder_interval = confounder_interval,
      confounder_baseline_mean = confounder_baseline_mean,
      confounder_gap_baseline = confounder_gap_baseline,
      confounder_gap_peak = confounder_gap_peak,
      confounder_gap_floor = confounder_gap_floor,
      confounder_sd = confounder_sd
    )

    sim_result$sim_id <- i
    sim_result$seed <- seed_start + i - 1
    return(sim_result)
  }

  # Run simulations (parallel or sequential)
  if (n_cores > 1) {
    if (verbose) {
      cat("Running", n_sim, "simulations on", n_cores, "cores...\n")
    }

    # Use parallel package (works on Mac/Linux; Windows uses 1 core)
    if (.Platform$OS.type == "windows") {
      warning("Parallel computing with mclapply not supported on Windows. Using 1 core.")
      n_cores <- 1
    }

    results_list <- parallel::mclapply(
      1:n_sim,
      run_single_sim,
      mc.cores = n_cores
    )

    if (verbose) {
      cat("Done.\n")
    }
  } else {
    # Sequential execution with progress
    results_list <- vector("list", n_sim)

    for (i in 1:n_sim) {
      if (verbose && i %% 10 == 0) {
        cat("Simulation", i, "of", n_sim, "\n")
      }
      results_list[[i]] <- run_single_sim(i)
    }
  }

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  return(results_df)
}

# ============================================================================
# Summarize Simulation Results
# ============================================================================

summarize_simulation <- function(results_df) {
  summary_df <- aggregate(
    cbind(log_hr, bias) ~ method,
    data = results_df,
    FUN = function(x) c(
      mean = mean(x),
      sd = sd(x),
      median = median(x),
      q025 = quantile(x, 0.025),
      q975 = quantile(x, 0.975)
    )
  )

  # Flatten the matrix columns
  summary_flat <- data.frame(
    method = summary_df$method,
    mean_log_hr = summary_df$log_hr[, "mean"],
    sd_log_hr = summary_df$log_hr[, "sd"],
    median_log_hr = summary_df$log_hr[, "median"],
    mean_bias = summary_df$bias[, "mean"],
    sd_bias = summary_df$bias[, "sd"],
    coverage_95 = NA
  )

  # Calculate 95% CI coverage (does true value fall within 2.5-97.5 percentile?)
  true_log_hr <- results_df$true_log_hr[1]
  for (m in unique(results_df$method)) {
    method_data <- results_df[results_df$method == m, "log_hr"]
    ci_lower <- quantile(method_data, 0.025)
    ci_upper <- quantile(method_data, 0.975)
    summary_flat$coverage_95[summary_flat$method == m] <-
      (true_log_hr >= ci_lower) & (true_log_hr <= ci_upper)
  }

  summary_flat$true_log_hr <- true_log_hr
  summary_flat$mean_hr <- exp(summary_flat$mean_log_hr)
  summary_flat$true_hr <- exp(true_log_hr)

  return(summary_flat)
}

# ============================================================================
# Example Usage
# ============================================================================

if (interactive()) {
  # Single run
  set.seed(123)
  results <- compare_adjustment_methods(n_pairs = 2000)
  print(results)

  # Simulation study (100 iterations)
  cat("\n=== Running Simulation Study ===\n")
  sim_results <- run_simulation_study(n_sim = 1000, n_pairs = 1000, n_cores = 10)

  cat("\n=== Summary ===\n")
  summary <- summarize_simulation(sim_results)
  print(summary)
}

