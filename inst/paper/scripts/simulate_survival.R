# ============================================================================
# Simulate Survival Data with Treatment Switching and Time-Dependent Confounding
# ============================================================================
#
# This function simulates time-to-first-event data for two cohorts:
# - Switchers: Start on treatment A, switch to treatment B
# - Continuers: Stay on treatment A throughout
#
# Key features:
# - Pre-switching: Switchers have higher event risk (confounding)
# - Post-switching: Treatment effect + confounding from event history
# - pre_event: Binary indicator of first pre-switch event (measured confounder)
#
# Output: One row per patient with the following columns:
# - id: Patient ID
# - cohort: "switcher" or "continuer"
# - T_max: Total follow-up time
# - switch_time: Time of switching (or pseudo-switching for continuers)
# - pre_event: 1 if first event occurred before switching, 0 otherwise
# - pre_event_time: Time to first pre-switch event (= switch_time if no event)
# - post_event: 1 if first event occurred after switching, 0 otherwise
# - post_event_time: Time to first post-switch event (= T_max if no event)
# ============================================================================

simulate_survival_data <- function(
  n_switchers,           # Sample size for switcher cohort
  n_continuers,          # Sample size for continuer cohort
  beta_treatment,        # Treatment effect (e.g., log(0.5) for HR=0.5)
  beta_confounder,       # Confounder effect (makes switchers higher risk pre-switch)
  beta_event_flag,       # Effect of having pre-switch event on post-switch hazard
  lambda_0,              # Baseline hazard scale (higher = lower risk, e.g., 10)
  shape,                 # Weibull shape parameter (>1 = accelerating hazard, e.g., 1.5)
  T_min,                 # Minimum follow-up time
  T_max,                 # Maximum follow-up time
  switch_start,          # Start of switching window (e.g., 0.25)
  switch_end             # End of switching window (e.g., 0.75)
) {

  # Total sample size
  n_total <- n_switchers + n_continuers

  # Generate individual follow-up times
  T_i <- runif(n_total, T_min, T_max)

  # Generate switching times (proportional to individual follow-up)
  switch_time <- runif(n_total, switch_start * T_i, switch_end * T_i)

  # Assign cohorts
  cohort <- c(rep("switcher", n_switchers), rep("continuer", n_continuers))

  # Initialize result data frame
  result <- data.frame(
    id = 1:n_total,
    cohort = cohort,
    T_max = T_i,
    switch_time = switch_time,
    pre_event = 0,
    pre_event_time = switch_time,  # Default: censored at switching
    post_event = 0,
    post_event_time = T_i          # Default: censored at end of follow-up
  )

  # Simulate events for each patient
  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max <- T_i[i]

    # PRE-SWITCHING PERIOD: Simulate first event
    # Switchers have higher hazard (beta_confounder)
    # For Weibull proportional hazards: scale_i = lambda_0 / exp(beta * X_i / shape)
    # Higher hazard → smaller scale
    scale_pre <- ifelse(is_switcher,
                        lambda_0 / exp(beta_confounder / shape),  # Switchers: higher risk
                        lambda_0)                                 # Continuers: baseline risk

    # Generate time to first pre-switch event using Weibull
    time_to_event <- rweibull(1, shape = shape, scale = scale_pre)

    if (time_to_event < t_switch) {
      # Event occurred before switching
      result$pre_event[i] <- 1
      result$pre_event_time[i] <- time_to_event
      event_flag <- 1
    } else {
      # No event before switching (censored at switch_time)
      event_flag <- 0
    }

    # POST-SWITCHING PERIOD: Simulate first event after switching
    # Hazard depends on event_flag and treatment (for switchers)
    # Calculate linear predictor: beta * X
    linear_predictor <- beta_event_flag * event_flag
    if (is_switcher) {
      linear_predictor <- linear_predictor + beta_treatment
    }

    # Calculate Weibull scale parameter for post-switching period
    # Must divide by shape to maintain proportional hazards
    scale_post <- lambda_0 / exp(linear_predictor / shape)

    # Use conditional Weibull sampling to maintain continuous baseline hazard
    # Sample from P(T | T > t_switch) where T ~ Weibull(shape, scale_post)
    # This ensures hazard is based on calendar time, not time-since-switch
    u <- runif(1)
    F_switch <- pweibull(t_switch, shape = shape, scale = scale_post)

    # Conditional sampling: map U to truncated range [F_switch, 1]
    time_to_post_event <- qweibull(u * (1 - F_switch) + F_switch,
                                    shape = shape,
                                    scale = scale_post)

    if (time_to_post_event < t_max) {
      # Event occurred after switching
      result$post_event[i] <- 1
      result$post_event_time[i] <- time_to_post_event
    }
    # Otherwise remains censored at T_max (already set)
  }

  return(result)
}

# ============================================================================
# Example Usage
# ============================================================================

set.seed(123)
data <- simulate_survival_data(
  n_switchers = 200,
  n_continuers = 200,
  beta_treatment = log(0.5),      # True treatment HR = 0.5
  beta_confounder = log(2),        # Switchers vs continuers pre-switch HR = 2
  beta_event_flag = log(2),        # pre_event effect HR = 2
  lambda_0 = 10,                   # Baseline scale (larger = lower event rate)
  shape = 1.5,                     # Weibull shape (>1 = accelerating hazard)
  T_min = 2,
  T_max = 6,
  switch_start = 0.25,
  switch_end = 0.75
)

# View first few rows
cat("=== First 20 Patients ===\n")
print(head(data, 20))

# Summary statistics
cat("\n=== Summary Statistics ===\n")
print(summary(data))

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
