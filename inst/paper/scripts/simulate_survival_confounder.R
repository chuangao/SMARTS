# ============================================================================
# Simulate Survival Data with Time-Varying Confounding (Version 5 - MVN)
# ============================================================================
#
# This simulation generates realistic survival data where:
# - Hazard trajectories differ between continuers and switchers
# - Confounder trajectories are CONSISTENT with hazard trajectories
# - Both continuers and switchers can have time-varying confounders
# - Confounders at nearby time points are CORRELATED via MVN
#
# ============================================================================
# KEY PRINCIPLE: Confounder must match hazard trajectory
# ============================================================================
#
#   Shape > 1 (worsening hazard) --> Confounder INCREASING
#   Shape < 1 (improving hazard) --> Confounder DECREASING
#   Shape = 1 (constant hazard)  --> Confounder FLAT
#
# ============================================================================
# MVN CORRELATION STRUCTURE
# ============================================================================
#
# The confounder at each time point is:
#   confounder(t) = mean_trajectory(t) + random_deviation(t)
#
# Where random deviations are drawn from MVN with:
#   Corr(t1, t2) = ρ^|t1 - t2|
#
# This means:
#   - Nearby time points are highly correlated (ρ^0.5 for 0.5 apart)
#   - Distant time points are less correlated (ρ^5 for 5 apart)
#   - When ρ = 0, confounders are independent (original behavior)
#   - When ρ → 1, confounders are nearly identical across time
#
# WHY THIS MATTERS FOR SMARTS:
#   - SMARTS assigns pseudo-switching times to continuers
#   - These pseudo-switch times are near actual switch times
#   - With MVN correlation, confounder at pseudo-switch ≈ confounder at switch
#   - This makes SMARTS adjustment more effective
#
# ============================================================================
# THE 8 CLINICAL SCENARIOS (Proportional Hazards Version)
# ============================================================================
#
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
#
# Legend: Pre HR = confounding effect, Tx HR = treatment effect
#
# ============================================================================
# GAP STRUCTURE: gap_at_switch IS THE EXTREME POINT
# ============================================================================
#
# For switchers who are SICKER (positive gap):
#   - gap_at_switch is the MAXIMUM
#   - Gap trajectory: gap_baseline → gap_at_switch (max) → gap_end
#
# For switchers who are HEALTHIER (negative gap):
#   - gap_at_switch is the MINIMUM
#   - Gap trajectory: gap_baseline → gap_at_switch (min) → gap_end
#
# VISUAL EXAMPLE (Positive gap - switcher sicker):
#
#          t=0        t_switch        t_max
#           │            │              │
# Gap:   [+0.5]  →    [+1.5]   →    [+0.8]
#                     (max)
#
# VISUAL EXAMPLE (Negative gap - switcher healthier):
#
#          t=0        t_switch        t_max
#           │            │              │
# Gap:   [-0.3]  →    [-1.5]   →    [-0.8]
#                     (min)
#
# ============================================================================


# ============================================================================
# Helper Function: Calculate continuer's confounder at time t
# ============================================================================
#
# The continuer's confounder trajectory is determined by their hazard shape:
# - shape < 1: Patient improving → confounder DECREASES over time
# - shape > 1: Patient worsening → confounder INCREASES over time
# - shape = 1: Patient stable → confounder FLAT
#
# Parameters:
#   t                  - Current time point
#   t_max              - Maximum follow-up time for this patient
#   baseline_mean      - Confounder value at t=0
#   shape_continuer    - Weibull shape parameter for continuer
#   change_magnitude   - How much confounder changes over full follow-up
#
# Returns:
#   Mean confounder value for continuer at time t
#
# ============================================================================
calculate_continuer_confounder <- function(t, t_max, baseline_mean,
                                            shape_continuer, change_magnitude) {


  # Calculate the proportion of follow-up completed (0 to 1)
  time_proportion <- t / t_max

  if (shape_continuer > 1) {
    # ──────────────────────────────────────────────────────────────────────────
    # WORSENING CONTINUER (shape > 1)
    # ──────────────────────────────────────────────────────────────────────────
    # Applies to Scenarios 5, 6, 7, 8
    # Patient staying on treatment but gradually declining
    # Example: Chronic progressive disease, treatment slowing but not stopping
    # Confounder INCREASES over time
    # ──────────────────────────────────────────────────────────────────────────
    trend <- change_magnitude * time_proportion

  } else if (shape_continuer < 1) {
    # ──────────────────────────────────────────────────────────────────────────
    # IMPROVING CONTINUER (shape < 1)
    # ──────────────────────────────────────────────────────────────────────────
    # Applies to Scenarios 1, 2, 3, 4
    # Patient staying on treatment because it's working well
    # Confounder DECREASES over time
    # ──────────────────────────────────────────────────────────────────────────
    trend <- -change_magnitude * time_proportion

  } else {
    # ──────────────────────────────────────────────────────────────────────────
    # STABLE CONTINUER (shape = 1)
    # ──────────────────────────────────────────────────────────────────────────
    # Constant hazard, confounder stays flat
    # ──────────────────────────────────────────────────────────────────────────
    trend <- 0
  }

  return(baseline_mean + trend)
}


# ============================================================================
# Helper Function: Calculate switcher's gap from continuer at time t
# ============================================================================
#
# The "gap" represents how different the switcher is from the continuer:
#   - Positive gap: Switcher is SICKER than continuer
#   - Negative gap: Switcher is HEALTHIER than continuer
#   - Zero gap: Switcher and continuer have same confounder
#
# IMPORTANT: gap_at_switch is the EXTREME point of the trajectory:
#   - For positive gaps (switcher sicker): gap_at_switch is the MAXIMUM
#   - For negative gaps (switcher healthier): gap_at_switch is the MINIMUM
#
# Parameters:
#   t              - Current time point
#   t_switch       - Time of switching
#   t_max          - Maximum follow-up time
#   gap_baseline   - Gap at t=0 (switcher - continuer)
#   gap_at_switch  - Gap at t=t_switch (extreme point)
#   gap_end        - Gap at t=t_max
#
# Returns:
#   The gap value (switcher confounder - continuer confounder) at time t
#
# ============================================================================
calculate_switcher_gap <- function(t, t_switch, t_max,
                                    gap_baseline, gap_at_switch, gap_end) {

  if (t <= t_switch) {
    # ──────────────────────────────────────────────────────────────────────────
    # PRE-SWITCH PERIOD: Linear interpolation from baseline to gap_at_switch
    # ──────────────────────────────────────────────────────────────────────────
    #
    # For positive gap scenarios (switcher sicker):
    #   gap_at_switch > gap_baseline → gap INCREASES toward maximum
    #
    # For negative gap scenarios (switcher healthier):
    #   gap_at_switch < gap_baseline → gap DECREASES toward minimum
    #
    # ──────────────────────────────────────────────────────────────────────────
    if (t_switch > 0) {
      proportion <- t / t_switch
    } else {
      proportion <- 1
    }
    gap <- gap_baseline + (gap_at_switch - gap_baseline) * proportion

  } else {
    # ──────────────────────────────────────────────────────────────────────────
    # POST-SWITCH PERIOD: Linear interpolation from gap_at_switch to end
    # ──────────────────────────────────────────────────────────────────────────
    #
    # For positive gap scenarios:
    #   gap_end < gap_at_switch → gap DECREASES from maximum
    #
    # For negative gap scenarios:
    #   gap_end > gap_at_switch → gap INCREASES from minimum
    #
    # ──────────────────────────────────────────────────────────────────────────
    time_after_switch <- t - t_switch
    time_remaining <- t_max - t_switch

    if (time_remaining > 0) {
      proportion <- time_after_switch / time_remaining
      gap <- gap_at_switch + (gap_end - gap_at_switch) * proportion
    } else {
      gap <- gap_end
    }
  }

  return(gap)
}


# ============================================================================
# Main Simulation Function
# ============================================================================
#
# Generates paired survival data with:
# - Time-varying confounders for BOTH continuers and switchers
# - Differentiated hazard shapes by cohort and period
# - Consistent relationship between confounder and hazard trajectories
# - MVN-correlated random deviations (nearby time points are correlated)
#
# ============================================================================
simulate_survival_data_confounder <- function(
    n_pairs,                   # Number of matched pairs (total N = n_pairs * 2)
    beta_treatment,            # Treatment effect: log(HR), e.g., log(0.5) for HR=0.5
    beta_confounder,           # Confounder effect: log(HR per unit), e.g., log(1.3)
    lambda_0,                  # Baseline hazard scale (higher = lower baseline risk)

    # ──────────────────────────────────────────────────────────────────────────
    # HAZARD SHAPE PARAMETERS
    # ──────────────────────────────────────────────────────────────────────────
    # shape < 1: Decreasing hazard (patient improving)
    # shape = 1: Constant hazard (patient stable)
    # shape > 1: Increasing hazard (patient worsening)
    # ──────────────────────────────────────────────────────────────────────────
    shape_continuer,           # Continuer hazard shape (same for entire follow-up)
    shape_switcher_pre,        # Switcher hazard shape BEFORE switch
    shape_switcher_post,       # Switcher hazard shape AFTER switch

    # ──────────────────────────────────────────────────────────────────────────
    # FOLLOW-UP TIME PARAMETERS
    # ──────────────────────────────────────────────────────────────────────────
    T_min,                     # Minimum follow-up time (e.g., 2 years)
    T_max,                     # Maximum follow-up time (e.g., 6 years)
    switch_start,              # Switch window start as proportion (e.g., 0.25)
    switch_end,                # Switch window end as proportion (e.g., 0.75)

    # ──────────────────────────────────────────────────────────────────────────
    # CONFOUNDER PARAMETERS
    # ──────────────────────────────────────────────────────────────────────────
    confounder_interval,       # Measurement interval (e.g., 0.5 years)
    confounder_baseline_mean,  # Mean confounder at t=0 for continuers
    confounder_change_magnitude, # How much continuer confounder changes over follow-up

    # Gap parameters: switcher confounder = continuer confounder + gap
    # Positive gap = switcher sicker; Negative gap = switcher healthier
    # gap_at_switch is the EXTREME point (max for positive, min for negative)
    confounder_gap_baseline,   # Gap at t=0
    confounder_gap_at_switch,  # Gap at switch time (extreme point)
    confounder_gap_end,        # Gap at end of follow-up

    confounder_sd,             # Standard deviation of random variation

    # ──────────────────────────────────────────────────────────────────────────
    # MVN CORRELATION PARAMETER
    # ──────────────────────────────────────────────────────────────────────────
    # rho controls temporal autocorrelation: Corr(t1, t2) = rho^|t1 - t2|
    # - rho = 0: Independent confounders at each time (original behavior)
    # - rho = 0.7-0.9: Strong correlation between nearby measurements
    # - rho → 1: Nearly identical confounders across all time points
    # ──────────────────────────────────────────────────────────────────────────
    rho = 0.8                  # Temporal correlation parameter (default 0.8)
) {

  # ============================================================================
  # SETUP: Generate basic structure
  # ============================================================================

  n_total <- n_pairs * 2

  # Paired structure: each pair has one switcher and one continuer
  pair_id <- rep(1:n_pairs, each = 2)
  cohort <- rep(c("switcher", "continuer"), times = n_pairs)

  # Individual follow-up times (same within pair for comparability)
  T_i_pairs <- runif(n_pairs, T_min, T_max)
  T_i <- rep(T_i_pairs, each = 2)

  # Switching times (same within pair)
  # Switch occurs between switch_start and switch_end proportion of follow-up
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
    confounder_at_switch = NA,
    pre_event = 0,
    pre_event_time = switch_time,
    post_event = 0,
    post_event_time = T_i
  )

  # ============================================================================
  # BUILD MVN COVARIANCE MATRIX
  # ============================================================================
  #
  # The covariance matrix has structure:
  #   Cov(t_i, t_j) = sigma^2 * rho^|t_i - t_j|
  #
  # This is an AR(1)-like correlation structure where:
  #   - Variance is constant (sigma^2) at all time points
  #   - Correlation decays exponentially with time distance
  #   - Corr(t, t+1) = rho, Corr(t, t+2) = rho^2, etc.
  #
  # ============================================================================

  n_timepoints <- length(time_points)

  # Build correlation matrix: Corr[i,j] = rho^|t_i - t_j|
  time_dist_matrix <- abs(outer(time_points, time_points, "-"))
  corr_matrix <- rho ^ time_dist_matrix

  # Convert to covariance matrix: Cov = sigma^2 * Corr
  cov_matrix <- confounder_sd^2 * corr_matrix

  # Pre-generate MVN random deviations for all patients
  # Each row is one patient's random deviations across all time points
  if (rho > 0) {
    # Use MVN for correlated deviations
    random_deviations <- MASS::mvrnorm(n = n_total, mu = rep(0, n_timepoints),
                                        Sigma = cov_matrix)
  } else {
    # rho = 0: Fall back to independent normal (original behavior)
    random_deviations <- matrix(rnorm(n_total * n_timepoints, 0, confounder_sd),
                                 nrow = n_total, ncol = n_timepoints)
  }

  # ============================================================================
  # GENERATE CONFOUNDER TRAJECTORIES
  # ============================================================================
  #
  # For each patient at each time point, calculate confounder value:
  #
  # CONTINUER:
  #   confounder(t) = baseline_mean + trend(t) + random_deviation(t)
  #   where trend depends on shape_continuer
  #
  # SWITCHER:
  #   confounder(t) = continuer_confounder(t) + gap(t) + random_deviation(t)
  #   where gap interpolates: gap_baseline → gap_at_switch → gap_end
  #
  # Random deviations are drawn from MVN with Corr(t1,t2) = rho^|t1-t2|
  # This ensures nearby time points have similar random variation.
  #
  # ============================================================================

  confounder_matrix <- matrix(NA, nrow = n_total, ncol = n_timepoints)
  colnames(confounder_matrix) <- paste0("confounder_", round(time_points, 2))

  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max_i <- T_i[i]

    for (j in seq_along(time_points)) {
      t <- time_points[j]

      # ────────────────────────────────────────────────────────────────────────
      # Step 1: Calculate the mean trajectory at time t
      # ────────────────────────────────────────────────────────────────────────
      continuer_confounder_at_t <- calculate_continuer_confounder(
        t = t,
        t_max = t_max_i,
        baseline_mean = confounder_baseline_mean,
        shape_continuer = shape_continuer,
        change_magnitude = confounder_change_magnitude
      )

      if (is_switcher) {
        # ──────────────────────────────────────────────────────────────────────
        # SWITCHER: mean = continuer trajectory + gap
        # ──────────────────────────────────────────────────────────────────────
        gap <- calculate_switcher_gap(
          t = t,
          t_switch = t_switch,
          t_max = t_max_i,
          gap_baseline = confounder_gap_baseline,
          gap_at_switch = confounder_gap_at_switch,
          gap_end = confounder_gap_end
        )
        mean_confounder <- continuer_confounder_at_t + gap

      } else {
        # ──────────────────────────────────────────────────────────────────────
        # CONTINUER: mean = continuer trajectory
        # ──────────────────────────────────────────────────────────────────────
        mean_confounder <- continuer_confounder_at_t
      }

      # ────────────────────────────────────────────────────────────────────────
      # Step 2: Add MVN-correlated random deviation
      # ────────────────────────────────────────────────────────────────────────
      # The random deviation at this time point is correlated with deviations
      # at nearby time points. This means if a patient happens to have a
      # higher-than-expected confounder at t=2, they'll also tend to have
      # higher-than-expected values at t=1.5 and t=2.5.
      # ────────────────────────────────────────────────────────────────────────
      confounder_matrix[i, j] <- mean_confounder + random_deviations[i, j]
    }

    # ──────────────────────────────────────────────────────────────────────────
    # Record confounder at switch time (last measurement before switch)
    # ──────────────────────────────────────────────────────────────────────────
    idx_before_switch <- which(time_points < t_switch)
    if (length(idx_before_switch) > 0) {
      idx_use <- max(idx_before_switch)
      result$confounder_at_switch[i] <- confounder_matrix[i, idx_use]
    } else {
      result$confounder_at_switch[i] <- confounder_matrix[i, 1]
    }
  }

  # ============================================================================
  # SIMULATE EVENTS
  # ============================================================================
  #
  # For each patient, simulate:
  # 1. Pre-switch event (using shape_pre)
  # 2. Post-switch event (using shape_post, with treatment effect for switchers)
  #
  # The hazard shape determines whether risk increases or decreases over time:
  # - shape > 1: Risk accelerates (worsening patient)
  # - shape < 1: Risk decelerates (improving patient)
  # - shape = 1: Constant risk
  #
  # ============================================================================

  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max_i <- T_i[i]
    confounder_score <- result$confounder_at_switch[i]

    # ──────────────────────────────────────────────────────────────────────────
    # Determine hazard shape based on cohort and period
    # ──────────────────────────────────────────────────────────────────────────
    #
    # CONTINUERS: Same shape throughout (shape_continuer)
    #   - Scenarios 1-4: shape=0.9 (improving)
    #   - Scenarios 5-8: shape=1.2 (worsening)
    #
    # SWITCHERS PRE-SWITCH (shape_switcher_pre):
    #   - Scenarios 1,2,5,6: shape=1.4 (worsening → reason to switch)
    #   - Scenarios 3,4,7,8: shape=0.8 (improving → switch for other reasons)
    #
    # SWITCHERS POST-SWITCH (shape_switcher_post):
    #   - Scenarios 1,3,5,7: shape=0.7 (improving → switch helped or stayed well)
    #   - Scenarios 2,4,6,8: shape=1.3 (worsening → switch didn't help or regret)
    #
    # ──────────────────────────────────────────────────────────────────────────
    shape_pre <- if (is_switcher) shape_switcher_pre else shape_continuer
    shape_post <- if (is_switcher) shape_switcher_post else shape_continuer

    # ──────────────────────────────────────────────────────────────────────────
    # PRE-SWITCH EVENT SIMULATION
    # ──────────────────────────────────────────────────────────────────────────
    scale_pre <- lambda_0 / exp(beta_confounder * confounder_score / shape_pre)
    time_to_event <- rweibull(1, shape = shape_pre, scale = scale_pre)

    if (time_to_event < t_switch) {
      result$pre_event[i] <- 1
      result$pre_event_time[i] <- time_to_event
    }

    # ──────────────────────────────────────────────────────────────────────────
    # POST-SWITCH EVENT SIMULATION
    # ──────────────────────────────────────────────────────────────────────────
    # Switchers receive the treatment effect (beta_treatment)
    # Continuers do not receive treatment effect
    # ──────────────────────────────────────────────────────────────────────────
    linear_predictor <- beta_confounder * confounder_score
    if (is_switcher) {
      linear_predictor <- linear_predictor + beta_treatment
    }

    scale_post <- lambda_0 / exp(linear_predictor / shape_post)

    # Conditional Weibull sampling (event must occur after switch time)
    u <- runif(1)
    F_switch <- pweibull(t_switch, shape = shape_post, scale = scale_post)
    time_to_post_event <- qweibull(u * (1 - F_switch) + F_switch,
                                   shape = shape_post,
                                   scale = scale_post)

    if (time_to_post_event < t_max_i) {
      result$post_event[i] <- 1
      result$post_event_time[i] <- time_to_post_event
    }
  }

  # Combine result with confounder matrix
  result <- cbind(result, confounder_matrix)

  return(result)
}


# ============================================================================
# SCENARIO CONFIGURATIONS (Proportional Hazards Version)
# ============================================================================
#
# This version maintains the PROPORTIONAL HAZARDS assumption by using:
# - Same shape for all groups (continuers and switchers, pre and post)
# - Constant confounder gap (same at baseline, peak, and end)
# - Explicit target HR values with back-calculated parameters
#
# Framework:
# - 2 shapes: 0.9 (improving population) and 1.2 (worsening population)
# - 4 combinations per shape based on pre/post HR directions
# - Total: 8 scenarios
#
# Key equations:
# - Pre HR  = exp(beta_confounder × gap)           [confounding only]
# - Post HR = exp(beta_treatment + beta_confounder × gap)
#           = HR_treatment × Pre HR
#
# Back-calculation:
# - gap = log(Pre HR) / beta_confounder
# - beta_treatment = log(Post HR / Pre HR) = log(HR_treatment)
#
# ============================================================================

get_scenario_config <- function(scenario_num, beta_confounder = log(1.3)) {

  # ──────────────────────────────────────────────────────────────────────────
  # Target HR values
  # ──────────────────────────────────────────────────────────────────────────
  HR_pre_high <- 1.5    # Switcher sicker (positive gap)
  HR_pre_low <- 0.67    # Switcher healthier (negative gap)
  HR_post_high <- 1.3   # Post-switch HR > 1
  HR_post_low <- 0.5    # Post-switch HR < 1

  # Helper: calculate gap and beta_treatment from target HRs
  calc_params <- function(pre_hr, post_hr, beta_conf = beta_confounder) {
    gap <- log(pre_hr) / beta_conf
    beta_trt <- log(post_hr / pre_hr)
    hr_trt <- post_hr / pre_hr
    return(list(gap = gap, beta_treatment = beta_trt, hr_treatment = hr_trt))
  }

  scenarios <- list(

    # ════════════════════════════════════════════════════════════════════════
    # SCENARIOS 1-4: IMPROVING POPULATION (shape = 0.9)
    # ════════════════════════════════════════════════════════════════════════

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 1: Pre HR > 1, Post HR > 1
    # ────────────────────────────────────────────────────────────────────────
    # Switcher sicker, treatment helps slightly but not enough
    # Pre:  Switcher 1.5x hazard (sicker)
    # Post: Switcher 1.3x hazard (still worse off)
    # Treatment HR = 1.3/1.5 = 0.87 (slightly beneficial)
    # ────────────────────────────────────────────────────────────────────────
    s1 = local({
      p <- calc_params(HR_pre_high, HR_post_high)
      list(
        name = "Improving, sick switch, tx helps slightly",
        shape = 0.9,
        target_pre_hr = HR_pre_high,
        target_post_hr = HR_post_high,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.3
      )
    }),

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 2: Pre HR > 1, Post HR < 1
    # ────────────────────────────────────────────────────────────────────────
    # Switcher sicker, treatment helps a lot
    # Pre:  Switcher 1.5x hazard
    # Post: Switcher 0.5x hazard
    # Treatment HR = 0.5/1.5 = 0.33 (very beneficial)
    # ────────────────────────────────────────────────────────────────────────
    s2 = local({
      p <- calc_params(HR_pre_high, HR_post_low)
      list(
        name = "Improving, sick switch, tx helps a lot",
        shape = 0.9,
        target_pre_hr = HR_pre_high,
        target_post_hr = HR_post_low,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.3
      )
    }),

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 3: Pre HR < 1, Post HR > 1
    # ────────────────────────────────────────────────────────────────────────
    # Switcher healthier, treatment harms (regret)
    # Pre:  Switcher 0.67x hazard (healthier)
    # Post: Switcher 1.3x hazard (now worse)
    # Treatment HR = 1.3/0.67 = 1.94 (very harmful)
    # ────────────────────────────────────────────────────────────────────────
    s3 = local({
      p <- calc_params(HR_pre_low, HR_post_high)
      list(
        name = "Improving, healthy switch, tx harms (regret)",
        shape = 0.9,
        target_pre_hr = HR_pre_low,
        target_post_hr = HR_post_high,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.3
      )
    }),

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 4: Pre HR < 1, Post HR < 1
    # ────────────────────────────────────────────────────────────────────────
    # Switcher healthier, stays healthy
    # Pre:  Switcher 0.67x hazard
    # Post: Switcher 0.5x hazard
    # Treatment HR = 0.5/0.67 = 0.75 (slightly beneficial)
    # ────────────────────────────────────────────────────────────────────────
    s4 = local({
      p <- calc_params(HR_pre_low, HR_post_low)
      list(
        name = "Improving, healthy switch, stays healthy",
        shape = 0.9,
        target_pre_hr = HR_pre_low,
        target_post_hr = HR_post_low,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.3
      )
    }),

    # ════════════════════════════════════════════════════════════════════════
    # SCENARIOS 5-8: WORSENING POPULATION (shape = 1.2)
    # ════════════════════════════════════════════════════════════════════════

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 5: Pre HR > 1, Post HR > 1
    # ────────────────────────────────────────────────────────────────────────
    s5 = local({
      p <- calc_params(HR_pre_high, HR_post_high)
      list(
        name = "Worsening, sick switch, tx helps slightly",
        shape = 1.2,
        target_pre_hr = HR_pre_high,
        target_post_hr = HR_post_high,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.5
      )
    }),

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 6: Pre HR > 1, Post HR < 1
    # ────────────────────────────────────────────────────────────────────────
    s6 = local({
      p <- calc_params(HR_pre_high, HR_post_low)
      list(
        name = "Worsening, sick switch, tx helps a lot",
        shape = 1.2,
        target_pre_hr = HR_pre_high,
        target_post_hr = HR_post_low,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.5
      )
    }),

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 7: Pre HR < 1, Post HR > 1
    # ────────────────────────────────────────────────────────────────────────
    s7 = local({
      p <- calc_params(HR_pre_low, HR_post_high)
      list(
        name = "Worsening, healthy switch, tx harms (regret)",
        shape = 1.2,
        target_pre_hr = HR_pre_low,
        target_post_hr = HR_post_high,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.5
      )
    }),

    # ────────────────────────────────────────────────────────────────────────
    # SCENARIO 8: Pre HR < 1, Post HR < 1
    # ────────────────────────────────────────────────────────────────────────
    s8 = local({
      p <- calc_params(HR_pre_low, HR_post_low)
      list(
        name = "Worsening, healthy switch, stays healthy",
        shape = 1.2,
        target_pre_hr = HR_pre_low,
        target_post_hr = HR_post_low,
        hr_treatment = p$hr_treatment,
        beta_treatment = p$beta_treatment,
        gap = p$gap,
        confounder_change_magnitude = 0.5
      )
    })
  )

  scenario_key <- paste0("s", scenario_num)
  if (!scenario_key %in% names(scenarios)) {
    stop("Invalid scenario number. Must be 1-8.")
  }

  config <- scenarios[[scenario_key]]

  # ──────────────────────────────────────────────────────────────────────────
  # Add PH-compliant parameters: uniform shape, constant gap
  # ──────────────────────────────────────────────────────────────────────────
  config$shape_continuer <- config$shape
  config$shape_switcher_pre <- config$shape
  config$shape_switcher_post <- config$shape

  # For PH-compliant scenarios, gap is constant throughout
  # gap_at_switch is the extreme point (max for positive, min for negative)
  config$gap_baseline <- config$gap
  config$gap_at_switch <- config$gap
  config$gap_end <- config$gap

  return(config)
}


# ============================================================================
# UTILITY: Print scenario summary
# ============================================================================
print_scenario_summary <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat("       8 SCENARIOS (PH-Compliant with MVN Correlation)                          \n")
  cat("================================================================================\n")
  cat("\n")
  cat("Key features:\n")
  cat("  - Same shape for all groups (PH maintained)\n")
  cat("  - Constant confounder gap (gap_at_switch = gap_baseline = gap_end)\n")
  cat("  - Target HRs with back-calculated parameters\n")
  cat("  - MVN-correlated confounders: Corr(t1,t2) = rho^|t1-t2|\n")
  cat("\n")
  cat("Target values: Pre HR high=1.5, low=0.67 | Post HR high=1.3, low=0.5\n")
  cat("\n")
  cat("| # | Shape | Pre HR | Post HR | Tx HR | Clinical Story                      |\n")
  cat("|---|-------|--------|---------|-------|-------------------------------------|\n")
  cat("| 1 |  0.9  |  1.5   |   1.3   | 0.87  | Improving, sick switch, tx helps    |\n")
  cat("| 2 |  0.9  |  1.5   |   0.5   | 0.33  | Improving, sick switch, tx helps lot|\n")
  cat("| 3 |  0.9  |  0.67  |   1.3   | 1.94  | Improving, healthy switch, regret   |\n")
  cat("| 4 |  0.9  |  0.67  |   0.5   | 0.75  | Improving, healthy switch, stays ok |\n")
  cat("| 5 |  1.2  |  1.5   |   1.3   | 0.87  | Worsening, sick switch, tx helps    |\n")
  cat("| 6 |  1.2  |  1.5   |   0.5   | 0.33  | Worsening, sick switch, tx helps lot|\n")
  cat("| 7 |  1.2  |  0.67  |   1.3   | 1.94  | Worsening, healthy switch, regret   |\n")
  cat("| 8 |  1.2  |  0.67  |   0.5   | 0.75  | Worsening, healthy switch, stays ok |\n")
  cat("\n")
  cat("Pre HR  = exp(beta_confounder × gap)  [confounding only]\n")
  cat("Post HR = Tx HR × Pre HR              [confounding + treatment]\n")
  cat("\n")
  cat("MVN Correlation: rho controls temporal autocorrelation\n")
  cat("  - rho=0: Independent confounders (old behavior)\n")
  cat("  - rho=0.8 (default): Strong temporal correlation\n")
  cat("  - Nearby time points have similar confounder values\n")
  cat("================================================================================\n")
}
