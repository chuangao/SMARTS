# ============================================================================
# Validate SMARTS Method: Pseudo-Switching Time Assignment
# ============================================================================
#
# This script validates the SMARTS random_assign() function by:
# 1. Simulating survival data with time-varying confounding
# 2. Using SMARTS to assign pseudo-switching times to continuers
# 3. Re-deriving events based on new pseudo-switching times
# 4. Looking up confounders at new pseudo-switching times
# 5. Balancing cohorts and running Cox PH
# 6. Comparing estimated HR to true treatment effect
# ============================================================================

library(survival)
library(MatchIt)
library(dplyr)
library(SMARTS)

source("simulate_survival_confounder.R")

# ============================================================================
# Helper Function: Re-derive events based on new pseudo-switching time
# ============================================================================

rederive_events <- function(data, new_switch_time_col = "new_switch_time") {

  data$new_switch_time <- data[[new_switch_time_col]]

  # Initialize new columns
  data$new_pre_event <- NA
  data$new_pre_event_time <- NA
  data$new_post_event <- NA
  data$new_post_event_time <- NA

  for (i in seq_len(nrow(data))) {
    new_switch <- data$new_switch_time[i]
    pre_event <- data$pre_event[i]
    pre_event_time <- data$pre_event_time[i]
    post_event <- data$post_event[i]
    post_event_time <- data$post_event_time[i]
    t_max <- data$T_max[i]

    # Case 1: new_switch_time < pre_event_time
    if (new_switch < pre_event_time) {
      # New pre-switch: no event before new switch
      data$new_pre_event[i] <- 0
      data$new_pre_event_time[i] <- new_switch  # censored at new switch

      # New post-switch: first event after new switch is the original pre_event
      data$new_post_event[i] <- max(pre_event, post_event)
      data$new_post_event_time[i] <- ifelse(
        pre_event == 1,
        pre_event_time - new_switch,
        post_event_time - new_switch
      )
      
    }
    # Case 2: pre_event_time <= new_switch_time < post_event_time
    else if (new_switch >= pre_event_time && new_switch < post_event_time) {
      # New pre-switch: same as original
      data$new_pre_event[i] <- pre_event
      data$new_pre_event_time[i] <- pre_event_time

      # New post-switch: original post_event
      data$new_post_event[i] <- post_event
      data$new_post_event_time[i] <- post_event_time - new_switch
    }
    # Case 3: new_switch_time >= post_event_time
    else {
      # New pre-switch: event occurred (the post_event is now pre-switch)
      data$new_pre_event[i] <- max(pre_event, post_event)
      data$new_pre_event_time[i] <- ifelse(
        pre_event == 1,
        pre_event_time,
        ifelse(
          post_event == 1,
          post_event_time,
          new_switch
        )
      )
      #data$new_pre_event_time[i] <- post_event_time

      # New post-switch: censored
      data$new_post_event[i] <- 0
      data$new_post_event_time[i] <- t_max - new_switch
    }
  }

  return(data)
}

# ============================================================================
# Helper Function: Look up confounder at new pseudo-switching time
# ============================================================================

lookup_confounder_at_time <- function(data, time_col = "new_switch_time",
                                       confounder_interval = 0.5) {

  # Get confounder column names
  confounder_cols <- grep("^confounder_[0-9]", colnames(data), value = TRUE)

  # Extract time points from column names
  time_points <- as.numeric(gsub("confounder_", "", confounder_cols))

  data$new_confounder_at_switch <- NA

  for (i in seq_len(nrow(data))) {
    t_switch <- data[[time_col]][i]

    # Find last time point before switching
    idx_before <- which(time_points < t_switch)

    if (length(idx_before) > 0) {
      idx_use <- max(idx_before)
      data$new_confounder_at_switch[i] <- data[i, confounder_cols[idx_use]]
    } else {
      # If no measurement before switch, use baseline (first column)
      data$new_confounder_at_switch[i] <- data[i, confounder_cols[1]]
    }
  }

  return(data)
}

# ============================================================================
# Main Validation Function
# ============================================================================


n_pairs = 1000
beta_treatment = log(0.5)
beta_confounder = log(1.3)
lambda_0 = 10
shape = 1.5
T_min = 2
T_max = 6
switch_start = 0.25
switch_end = 0.75
confounder_interval = 0.5
confounder_baseline_mean = 2.5
confounder_gap_baseline = 0.5
confounder_gap_peak = 1.6
confounder_gap_end = 0.8
confounder_sd = 0.8
nbin = 10
seed = 123


validate_smarts <- function(
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
  confounder_gap_end = 0.8,
  confounder_sd = 0.8,
  nbin = 10,
  seed = 123
) {

  set.seed(seed)

  # -------------------------------------------------------------------------
  # Step 1: Simulate data
  # -------------------------------------------------------------------------
  cat("Step 1: Simulating survival data...\n")

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
    confounder_gap_end = confounder_gap_end,
    confounder_sd = confounder_sd
  )

  cat("  Simulated", nrow(data), "patients\n")

  # -------------------------------------------------------------------------
  # Step 2: Split into switchers and continuers
  # -------------------------------------------------------------------------
  cat("\nStep 2: Splitting into switchers and continuers...\n")

  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  cat("  Switchers:", nrow(switchers), "\n")
  cat("  Continuers:", nrow(continuers), "\n")

  # Store original pseudo-switching time for continuers (from simulation)
  continuers$original_switch_time <- continuers$switch_time

  # Prepare data for SMARTS (continuers should not have switch_time for assignment)
  # We'll keep the column but SMARTS will overwrite it

  # -------------------------------------------------------------------------
  # Step 3: Use SMARTS random_assign() to assign pseudo-switching times
  # -------------------------------------------------------------------------
  cat("\nStep 3: Using SMARTS to assign pseudo-switching times...\n")

  # Prepare data list for SMARTS
  # Need to rename columns to match SMARTS expected format
  switchers_smarts <- switchers
  switchers_smarts$swi_yrs <- switchers_smarts$switch_time
  switchers_smarts$fup_yrs <- switchers_smarts$T_max

  continuers_smarts <- continuers
  continuers_smarts$swi_yrs <- NA  # Will be assigned by SMARTS
  continuers_smarts$fup_yrs <- continuers_smarts$T_max

  smarts_input <- list(
    cont = continuers_smarts,
    swi = switchers_smarts
  )

  # Run SMARTS random_assign
  smarts_result <- random_assign(
    smarts_input,
    nbin = nbin,
    seed = seed,
    swi_time = "swi_yrs",
    cens_time = "fup_yrs"
  )

  # Get assigned data
  continuers_assigned <- smarts_result$assigned$cont
  switchers_assigned <- smarts_result$assigned$swi

  cat("  Assigned continuers:", nrow(continuers_assigned), "\n")
  cat("  Unassigned continuers:", nrow(smarts_result$unassigned$cont), "\n")
  cat("  Used switchers:", nrow(switchers_assigned), "\n")

  # New pseudo-switching time for continuers
  continuers_assigned$new_switch_time <- continuers_assigned$swi_yrs

  # For switchers, new_switch_time is same as original
  switchers_assigned$new_switch_time <- switchers_assigned$switch_time

  # -------------------------------------------------------------------------
  # Step 4: Re-derive events based on new pseudo-switching times
  # -------------------------------------------------------------------------
  cat("\nStep 4: Re-deriving events based on new pseudo-switching times...\n")

  continuers_assigned <- rederive_events(continuers_assigned, "new_switch_time")
  switchers_assigned <- rederive_events(switchers_assigned, "new_switch_time")

  # -------------------------------------------------------------------------
  # Step 5: Look up confounder at new pseudo-switching time
  # -------------------------------------------------------------------------
  cat("\nStep 5: Looking up confounders at new pseudo-switching times...\n")

  continuers_assigned <- lookup_confounder_at_time(
    continuers_assigned, "new_switch_time", confounder_interval
  )
  switchers_assigned <- lookup_confounder_at_time(
    switchers_assigned, "new_switch_time", confounder_interval
  )

  # -------------------------------------------------------------------------
  # Step 6: Combine data for analysis
  # -------------------------------------------------------------------------
  cat("\nStep 6: Preparing data for analysis...\n")

  # Select relevant columns
  cols_to_keep <- c("id", "cohort", "T_max", "new_switch_time",
                    "new_pre_event", "new_pre_event_time",
                    "new_post_event", "new_post_event_time",
                    "new_confounder_at_switch")

  analysis_data <- rbind(
    continuers_assigned[, cols_to_keep],
    switchers_assigned[, cols_to_keep]
  )

  # Create treatment indicator
  analysis_data$treated <- as.numeric(analysis_data$cohort == "switcher")

  # Remove rows with invalid post-event times (should be > 0)
  analysis_data <- analysis_data[analysis_data$new_post_event_time > 0, ]

  cat("  Analysis sample size:", nrow(analysis_data), "\n")

  # -------------------------------------------------------------------------
  # Step 7: Run Cox PH models
  # -------------------------------------------------------------------------
  cat("\nStep 7: Running Cox PH models...\n")

  # 7a. Naive (unadjusted)
  cox_naive <- coxph(
    Surv(new_post_event_time, new_post_event) ~ treated,
    data = analysis_data
  )
  beta_naive <- coef(cox_naive)["treated"]

  # 7b. Adjusted for confounder only
  cox_adj_conf <- coxph(
    Surv(new_post_event_time, new_post_event) ~ treated + new_confounder_at_switch,
    data = analysis_data
  )
  beta_adj_conf <- coef(cox_adj_conf)["treated"]

  # 7c. Adjusted for confounder + pre-event
  cox_adj_full <- coxph(
    Surv(new_post_event_time, new_post_event) ~ treated + new_confounder_at_switch + new_pre_event,
    data = analysis_data
  )
  beta_adj_full <- coef(cox_adj_full)["treated"]

  # 7d. IPTW
  ps_model <- glm(
    treated ~ new_confounder_at_switch,
    data = analysis_data,
    family = binomial
  )
  analysis_data$ps <- predict(ps_model, type = "response")
  analysis_data$weight <- ifelse(
    analysis_data$treated == 1,
    1 / analysis_data$ps,
    1 / (1 - analysis_data$ps)
  )

  cox_iptw <- coxph(
    Surv(new_post_event_time, new_post_event) ~ treated,
    data = analysis_data,
    weights = weight
  )
  beta_iptw <- coef(cox_iptw)["treated"]

  # 7e. Propensity Score Matching (PSM)
  match_out <- matchit(
    treated ~ new_confounder_at_switch,
    data = analysis_data,
    method = "nearest",
    caliper = 0.1
  )
  data_matched <- match.data(match_out)

  # Use stratified Cox to account for matched pair structure
  if (nrow(data_matched) > 0 && length(unique(data_matched$subclass)) > 1) {
    cox_psm <- coxph(
      Surv(new_post_event_time, new_post_event) ~ treated + strata(subclass),
      data = data_matched
    )
    beta_psm <- coef(cox_psm)["treated"]
    n_matched <- nrow(data_matched)
  } else {
    beta_psm <- NA
    n_matched <- 0
  }

  cat("  PSM matched pairs:", n_matched / 2, "\n")

  # -------------------------------------------------------------------------
  # Step 8: Compile results
  # -------------------------------------------------------------------------
  cat("\nStep 8: Compiling results...\n")

  results <- data.frame(
    method = c("Naive", "Adj (Confounder)", "Adj (Conf + PreEvent)", "IPTW", "PSM"),
    log_hr = c(beta_naive, beta_adj_conf, beta_adj_full, beta_iptw, beta_psm),
    hr = exp(c(beta_naive, beta_adj_conf, beta_adj_full, beta_iptw, beta_psm)),
    true_log_hr = beta_treatment,
    true_hr = exp(beta_treatment),
    bias = c(beta_naive, beta_adj_conf, beta_adj_full, beta_iptw, beta_psm) - beta_treatment
  )

  return(list(
    results = results,
    analysis_data = analysis_data,
    smarts_result = smarts_result
  ))
}

# ============================================================================
# Run Validation
# ============================================================================

if (interactive()) {

  cat("=== SMARTS Method Validation ===\n\n")

  # Single run
  validation <- validate_smarts(
    n_pairs = 5000,
    beta_treatment = log(0.5),  # True HR = 0.5
    seed = 123
  )

  cat("\n=== Results ===\n")
  print(validation$results)

  cat("\n=== Interpretation ===\n")
  cat("True treatment HR:", exp(log(0.5)), "\n")
  cat("True log(HR):", log(0.5), "\n\n")

  cat("Method comparison:\n")
  for (i in seq_len(nrow(validation$results))) {
    cat(sprintf("  %s: HR = %.3f (bias = %.3f)\n",
                validation$results$method[i],
                validation$results$hr[i],
                validation$results$bias[i]))
  }

  # Check distribution of pseudo-switching times
  cat("\n=== Pseudo-Switching Time Distribution ===\n")
  cat("Switchers (actual):\n")
  print(summary(validation$analysis_data$new_switch_time[validation$analysis_data$treated == 1]))
  cat("\nContinuers (SMARTS-assigned):\n")
  print(summary(validation$analysis_data$new_switch_time[validation$analysis_data$treated == 0]))

  # QQ plot of switching times
  cat("\n=== QQ Plot ===\n")
  cat("Use qqplot() to compare switcher vs continuer pseudo-switching times\n")

  swi_times <- validation$analysis_data$new_switch_time[validation$analysis_data$treated == 1]
  cont_times <- validation$analysis_data$new_switch_time[validation$analysis_data$treated == 0]

  qqplot(swi_times, cont_times,
         main = "QQ Plot: Switching Times",
         xlab = "Switchers (actual)",
         ylab = "Continuers (SMARTS-assigned)")
  abline(0, 1, col = "red", lty = 2)
}
