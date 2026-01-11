# ============================================================================
# Simulate Survival Data with Time-Varying Confounding (Paired Design)
# ============================================================================
#
# This function simulates time-to-first-event data for matched pairs:
# - Each pair: 1 switcher + 1 continuer
# - Same switching/pseudo-switching time within pair
# - Time-varying confounder scores
#
# Confounder pattern:
# - Pre-switching: Gap widens (switchers deteriorate faster)
# - Peak gap at switching time
# - Post-switching: Gap narrows to floor (switchers improve but remain worse)
#
# Output: One row per patient with confounder scores at each time interval
# ============================================================================

# Helper function: Calculate mean gap at time t relative to switching time
calculate_gap <- function(t, t_switch, t_max,
                          gap_baseline, gap_peak, gap_floor) {
  if (t <= t_switch) {
    # Pre-switching: linear increase from baseline to peak
    proportion <- t / t_switch
    gap <- gap_baseline + (gap_peak - gap_baseline) * proportion
  } else {
    # Post-switching: linear decrease from peak to floor
    time_after_switch <- t - t_switch
    time_remaining <- t_max - t_switch
    proportion <- time_after_switch / time_remaining
    gap <- gap_peak - (gap_peak - gap_floor) * proportion
  }
  return(gap)
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
  confounder_baseline_mean,# Baseline mean confounder for continuers (e.g., 2.5)
  confounder_gap_baseline, # Initial gap at t=0 (e.g., 0.5)
  confounder_gap_peak,     # Peak gap at switching (e.g., 1.6)
  confounder_gap_floor,    # Floor gap post-switching (e.g., 0.8)
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

      # Calculate mean gap at this time point
      gap <- calculate_gap(t, t_switch, t_max,
                           confounder_gap_baseline,
                           confounder_gap_peak,
                           confounder_gap_floor)

      # Generate confounder score
      if (is_switcher) {
        mean_confounder <- confounder_baseline_mean + gap
      } else {
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
    # For Weibull proportional hazards: scale = lambda_0 / exp(beta * X / shape)
    scale_pre <- lambda_0 / exp(beta_confounder * confounder_score / shape)

    # Generate time to first pre-switch event
    time_to_event <- rweibull(1, shape = shape, scale = scale_pre)

    if (time_to_event < t_switch) {
      # Event occurred before switching
      result$pre_event[i] <- 1
      result$pre_event_time[i] <- time_to_event
    }

    # POST-SWITCHING PERIOD: Simulate first event after switching
    # Hazard depends on confounder and treatment (for switchers)
    linear_predictor <- beta_confounder * confounder_score
    if (is_switcher) {
      linear_predictor <- linear_predictor + beta_treatment
    }

    # Calculate Weibull scale parameter for post-switching
    # Must divide by shape to maintain proportional hazards
    scale_post <- lambda_0 / exp(linear_predictor / shape)

    # Conditional Weibull sampling
    u <- runif(1)
    F_switch <- pweibull(t_switch, shape = shape, scale = scale_post)
    time_to_post_event <- qweibull(u * (1 - F_switch) + F_switch,
                                    shape = shape,
                                    scale = scale_post)

    if (time_to_post_event < t_max) {
      # Event occurred after switching
      result$post_event[i] <- 1
      result$post_event_time[i] <- time_to_post_event
    }
  }

  # Combine result with confounder matrix
  result <- cbind(result, confounder_matrix)

  return(result)
}

# ============================================================================
# Example Usage
# ============================================================================

set.seed(123)
data <- simulate_survival_data_confounder(
  n_pairs = 200,                      # 200 pairs = 400 patients
  beta_treatment = log(0.5),          # True treatment HR = 0.5
  beta_confounder = log(1.3),         # Confounder effect HR = 1.3 per point
  lambda_0 = 10,                      # Baseline hazard scale
  shape = 1.5,                        # Accelerating hazard
  T_min = 2,
  T_max = 6,
  switch_start = 0.25,
  switch_end = 0.75,
  confounder_interval = 0.5,          # Measure every 0.5 years
  confounder_baseline_mean = 2.5,     # Continuer baseline mean
  confounder_gap_baseline = 0.5,      # Initial gap
  confounder_gap_peak = 1.6,          # Peak gap at switching
  confounder_gap_floor = 0.8,         # Floor gap post-switching
  confounder_sd = 0.8                 # Random variation
)

# View first few rows (core variables)
cat("=== First 10 Patients (Core Variables) ===\n")
core_vars <- c("id", "pair_id", "cohort", "T_max", "switch_time",
               "confounder_at_switch", "pre_event", "post_event")
print(head(data[, core_vars], 10))

# Summary statistics
cat("\n=== Summary Statistics ===\n")
print(summary(data[, core_vars]))

# Confounder at switching by cohort
cat("\n=== Confounder at Switching by Cohort ===\n")
print(aggregate(confounder_at_switch ~ cohort, data,
                function(x) c(mean = mean(x), sd = sd(x))))

# Check pre_event distribution by cohort
cat("\n=== Pre-Event Distribution by Cohort ===\n")
print(table(data$cohort, data$pre_event))

# Check post_event distribution by cohort
cat("\n=== Post-Event Distribution by Cohort ===\n")
print(table(data$cohort, data$post_event))

# Event rates
cat("\n=== Event Rates ===\n")
cat("Pre-switch event rate:\n")
print(aggregate(pre_event ~ cohort, data, mean))
cat("\nPost-switch event rate:\n")
print(aggregate(post_event ~ cohort, data, mean))

# Verify pairing structure
cat("\n=== Verify Pairing (First 3 Pairs) ===\n")
for (p in 1:3) {
  pair_data <- data[data$pair_id == p, c("id", "cohort", "switch_time", "confounder_at_switch")]
  cat("\nPair", p, ":\n")
  print(pair_data, row.names = FALSE)
}
