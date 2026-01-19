# ============================================================================
# Simulate Survival Data with Time-Varying Confounding (Version 2)
# ============================================================================
#
# Key change from v1:
# - Continuers: CONSTANT confounder over time (stable on same treatment)
# - Switchers: DYNAMIC confounder (worsens -> peaks at switch -> improves)
#
# This is more clinically realistic:
# - Switchers deteriorate (reason to switch), then improve with new treatment
# - Continuers remain stable on their current treatment
# ============================================================================

# Helper function: Calculate switcher's confounder deviation from baseline
calculate_switcher_deviation <- function(t, t_switch, t_max,
                                          gap_baseline, gap_peak, gap_end) {
  if (t <= t_switch) {
    # Pre-switching: linear increase from baseline gap to peak
    proportion <- t / t_switch
    deviation <- gap_baseline + (gap_peak - gap_baseline) * proportion
  } else {
    # Post-switching: linear decrease from peak to end (can go to 0 or negative)
    time_after_switch <- t - t_switch
    time_remaining <- t_max - t_switch
    if (time_remaining > 0) {
      proportion <- time_after_switch / time_remaining
      deviation <- gap_peak - (gap_peak - gap_end) * proportion
    } else {
      deviation <- gap_end
    }
  }
  return(deviation)
}

simulate_survival_data_confounder <- function(
  n_pairs,                 # Number of matched pairs
  beta_treatment,          # Treatment effect (e.g., log(0.5) for HR=0.5)
  beta_confounder,         # Confounder effect on hazard (e.g., log(1.3) per point)
  lambda_0,                # Baseline hazard scale (higher = lower risk)
  shape,                   # Weibull shape parameter (>1 = accelerating hazard)
  T_min,                   # Minimum follow-up time
  T_max,                   # Maximum follow-up time
  switch_start,            # Start of switching window (e.g., 0.25)
  switch_end,              # End of switching window (e.g., 0.75)
  confounder_interval,     # Confounder measurement interval (e.g., 0.5 years)
  confounder_baseline_mean,# Baseline mean confounder (same for both cohorts at t=0)
  confounder_gap_baseline, # Switcher deviation at t=0 (e.g., 0.2)
  confounder_gap_peak,     # Switcher deviation at switch time (e.g., 1.5)
  confounder_gap_end,      # Switcher deviation at end (e.g., 0 or -0.2)
  confounder_sd            # SD for confounder random variation (e.g., 0.8)
) {

  n_total <- n_pairs * 2

  # Generate paired structure
  pair_id <- rep(1:n_pairs, each = 2)
  cohort <- rep(c("switcher", "continuer"), times = n_pairs)

  # Generate individual follow-up times (same within pair)
  T_i_pairs <- runif(n_pairs, T_min, T_max)
  T_i <- rep(T_i_pairs, each = 2)

  # Generate switching times (same within pair)
  switch_time_pairs <- runif(n_pairs, switch_start * T_i_pairs, switch_end * T_i_pairs)
  switch_time <- rep(switch_time_pairs, each = 2)

  # Time points for confounder measurements
  max_time <- max(T_i)
  time_points <- seq(0, max_time, by = confounder_interval)

  # Initialize result data frame
  result <- data.frame(
    id = 1:n_total,
    pair_id = pair_id,
    cohort = cohort,
    T_max = T_i,
    switch_time = switch_time,
    confounder_at_switch = NA,  # Will be filled
    pre_event = 0,
    pre_event_time = switch_time,
    post_event = 0,
    post_event_time = T_i
  )

  # Generate confounder scores for each patient
  confounder_matrix <- matrix(NA, nrow = n_total, ncol = length(time_points))
  colnames(confounder_matrix) <- paste0("confounder_", round(time_points, 2))

  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max <- T_i[i]

    for (j in 1:length(time_points)) {
      t <- time_points[j]

      if (is_switcher) {
        # Switchers: dynamic confounder (worsens -> peaks -> improves)
        deviation <- calculate_switcher_deviation(t, t_switch, t_max,
                                                   confounder_gap_baseline,
                                                   confounder_gap_peak,
                                                   confounder_gap_end)
        mean_confounder <- confounder_baseline_mean + deviation
      } else {
        # Continuers: constant confounder at baseline
        mean_confounder <- confounder_baseline_mean
      }

      confounder_matrix[i, j] <- rnorm(1, mean = mean_confounder, sd = confounder_sd)
    }

    # Find confounder at last time point before switching
    idx_before_switch <- which(time_points < t_switch)
    if (length(idx_before_switch) > 0) {
      idx_use <- max(idx_before_switch)
      result$confounder_at_switch[i] <- confounder_matrix[i, idx_use]
    } else {
      # If no measurement before switch, use baseline
      result$confounder_at_switch[i] <- confounder_matrix[i, 1]
    }
  }

  # Simulate events for each patient
  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max <- T_i[i]
    confounder_score <- result$confounder_at_switch[i]

    # PRE-SWITCHING PERIOD: Simulate first event
    # Hazard depends on confounder score
    scale_pre <- lambda_0 / exp(beta_confounder * confounder_score / shape)

    # Generate time to first pre-switch event
    time_to_event <- rweibull(1, shape = shape, scale = scale_pre)

    if (time_to_event < t_switch) {
      result$pre_event[i] <- 1
      result$pre_event_time[i] <- time_to_event
    }

    # POST-SWITCHING PERIOD: Simulate first event after switching
    linear_predictor <- beta_confounder * confounder_score
    if (is_switcher) {
      linear_predictor <- linear_predictor + beta_treatment
    }

    scale_post <- lambda_0 / exp(linear_predictor / shape)

    # Conditional Weibull sampling
    u <- runif(1)
    F_switch <- pweibull(t_switch, shape = shape, scale = scale_post)
    time_to_post_event <- qweibull(u * (1 - F_switch) + F_switch,
                                    shape = shape,
                                    scale = scale_post)

    if (time_to_post_event < t_max) {
      result$post_event[i] <- 1
      result$post_event_time[i] <- time_to_post_event
    }
  }

  # Combine result with confounder matrix
  result <- cbind(result, confounder_matrix)

  return(result)
}
