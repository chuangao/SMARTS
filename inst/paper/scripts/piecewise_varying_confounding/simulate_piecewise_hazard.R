# ============================================================================
# Simulate Recurring Events with Piecewise Constant Hazard
# WITH TIME-VARYING CONFOUNDING
# ============================================================================
#
# Time is divided into small segments, each with its own constant hazard.
# Hazard can gradually increase or decrease across segments.
# Both switchers and continuers experience the same time-varying baseline.
# Treatment effect applies only to switchers after their switch time.
#
# TIME-VARYING CONFOUNDING:
# - Continuers have stable confounder throughout follow-up
# - Switchers have confounder gap that changes over time:
#     gap_baseline → gap_at_switch (peak) → gap_end
# - Confounder values at nearby time points are correlated (MVN)
#
# ============================================================================

library(MASS)  # For mvrnorm

# ============================================================================
# Helper: Calculate switcher's gap from continuer at time t
# ============================================================================
calculate_switcher_gap <- function(t, t_switch, t_max,
                                    gap_baseline, gap_at_switch, gap_end) {
  # gap_at_switch is the EXTREME point:
  #   - For positive gaps (switcher sicker): gap_at_switch is the MAXIMUM
  #   - For negative gaps (switcher healthier): gap_at_switch is the MINIMUM

  if (t <= t_switch) {
    # PRE-SWITCH: Linear interpolation from baseline to gap_at_switch
    if (t_switch > 0) {
      proportion <- t / t_switch
    } else {
      proportion <- 1
    }
    gap <- gap_baseline + (gap_at_switch - gap_baseline) * proportion

  } else {
    # POST-SWITCH: Linear interpolation from gap_at_switch to end
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
# Generate events with piecewise constant hazard
# Now uses time-varying confounder
# ============================================================================
generate_piecewise_events_varying <- function(t_start, t_end, segment_hazards,
                                               segment_times,
                                               confounder_values,  # vector of confounder at each segment
                                               confounder_times,   # time points for confounder
                                               beta_confounder,
                                               treatment_effect = 0,
                                               treatment_start = Inf) {

  if (t_end <= t_start) return(numeric(0))

  event_times <- c()
  t_current <- t_start

  while (t_current < t_end) {
    # Find which segment we're in
    seg_idx <- findInterval(t_current, segment_times, rightmost.closed = TRUE)
    seg_idx <- max(1, min(seg_idx, length(segment_hazards)))

    # Get segment end time
    seg_end <- if (seg_idx < length(segment_times)) segment_times[seg_idx + 1] else Inf
    seg_end <- min(seg_end, t_end)

    # Base hazard for this segment
    base_hazard <- segment_hazards[seg_idx]

    # Get confounder value for current time (use nearest time point)
    conf_idx <- findInterval(t_current, confounder_times, rightmost.closed = TRUE)
    conf_idx <- max(1, min(conf_idx, length(confounder_values)))
    current_confounder <- confounder_values[conf_idx]

    # Apply confounder effect
    linear_pred <- beta_confounder * current_confounder
    hazard <- base_hazard * exp(linear_pred)

    # Apply treatment effect if we're past treatment start
    if (t_current >= treatment_start) {
      hazard <- hazard * exp(treatment_effect)
    }

    # Generate next event time using exponential with this hazard
    if (hazard > 0) {
      inter_arrival <- rexp(1, rate = hazard)
    } else {
      inter_arrival <- Inf
    }

    t_next <- t_current + inter_arrival

    # Determine the effective boundary
    effective_end <- min(seg_end, t_end)
    if (t_current < treatment_start) {
      effective_end <- min(effective_end, treatment_start)
    }

    if (t_next < effective_end) {
      event_times <- c(event_times, t_next)
      t_current <- t_next
    } else {
      t_current <- effective_end
    }
  }
  return(event_times)
}

# ============================================================================
# Create hazard schedule (increasing or decreasing over time)
# ============================================================================
create_hazard_schedule <- function(max_time, segment_length = 0.5,
                                    base_hazard = 0.01,
                                    hazard_trend = "constant",
                                    trend_strength = 0.1) {

  segment_times <- seq(0, max_time + segment_length, by = segment_length)
  n_segments <- length(segment_times) - 1

  segment_hazards <- numeric(n_segments)

  for (i in 1:n_segments) {
    mid_time <- (segment_times[i] + segment_times[i + 1]) / 2

    if (hazard_trend == "constant") {
      segment_hazards[i] <- base_hazard
    } else if (hazard_trend == "increasing") {
      segment_hazards[i] <- base_hazard * (1 + trend_strength * mid_time)
    } else if (hazard_trend == "decreasing") {
      segment_hazards[i] <- base_hazard * max(0.1, 1 - trend_strength * mid_time)
    }
  }

  return(list(times = segment_times, hazards = segment_hazards))
}

# ============================================================================
# Main Simulation Function with Time-Varying Confounding
# ============================================================================
simulate_piecewise_hazard <- function(
    n_pairs,
    beta_treatment,            # Treatment effect (log HR)
    beta_confounder,           # Confounder effect (log HR per unit)

    # Hazard parameters
    base_hazard = 0.01,
    hazard_trend = "constant",
    trend_strength = 0.1,
    segment_length = 0.5,

    # Follow-up parameters
    T_min = 2,
    T_max = 6,
    switch_start = 0.25,
    switch_end = 0.75,

    # Confounder parameters
    confounder_interval = 0.5,      # Measurement interval
    confounder_baseline_mean = 2.5,
    confounder_sd = 0.8,

    # Gap parameters for time-varying confounding
    # gap = switcher confounder - continuer confounder
    confounder_gap_baseline = 0.5,  # Gap at t=0
    confounder_gap_at_switch = 1.5, # Gap peaks at switch (maximum for sicker switchers)
    confounder_gap_end = 0.8,       # Gap at end of follow-up

    # MVN correlation for confounder over time
    rho = 0.9                       # Temporal correlation: Corr(t1,t2) = rho^|t1-t2|
) {

  n_total <- n_pairs * 2

  # Setup paired structure
  pair_id <- rep(1:n_pairs, each = 2)
  cohort <- rep(c("switcher", "continuer"), times = n_pairs)

  # Follow-up times (same within pair)
  T_i_pairs <- runif(n_pairs, T_min, T_max)
  T_i <- rep(T_i_pairs, each = 2)

  # Switch times (same within pair)
  switch_time_pairs <- runif(n_pairs, switch_start * T_i_pairs, switch_end * T_i_pairs)
  switch_time <- rep(switch_time_pairs, each = 2)

  # Create hazard schedule
  max_follow_up <- max(T_i)
  hazard_schedule <- create_hazard_schedule(
    max_time = max_follow_up,
    segment_length = segment_length,
    base_hazard = base_hazard,
    hazard_trend = hazard_trend,
    trend_strength = trend_strength
  )

  # Time points for confounder measurements
  time_points <- seq(0, max_follow_up, by = confounder_interval)
  n_timepoints <- length(time_points)

  # ============================================================================
  # BUILD MVN COVARIANCE MATRIX for correlated confounders
  # ============================================================================
  time_dist_matrix <- abs(outer(time_points, time_points, "-"))
  corr_matrix <- rho ^ time_dist_matrix
  cov_matrix <- confounder_sd^2 * corr_matrix

  # Pre-generate MVN random deviations for all patients
  if (rho > 0) {
    random_deviations <- mvrnorm(n = n_total, mu = rep(0, n_timepoints),
                                  Sigma = cov_matrix)
  } else {
    random_deviations <- matrix(rnorm(n_total * n_timepoints, 0, confounder_sd),
                                 nrow = n_total, ncol = n_timepoints)
  }

  # ============================================================================
  # GENERATE CONFOUNDER TRAJECTORIES
  # ============================================================================
  confounder_matrix <- matrix(NA, nrow = n_total, ncol = n_timepoints)
  colnames(confounder_matrix) <- paste0("confounder_", round(time_points, 2))
  confounder_at_switch <- numeric(n_total)

  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max_i <- T_i[i]

    for (j in seq_along(time_points)) {
      t <- time_points[j]

      # Continuer confounder: stable at baseline mean
      continuer_confounder <- confounder_baseline_mean

      if (is_switcher) {
        # Switcher: continuer trajectory + time-varying gap
        gap <- calculate_switcher_gap(
          t = t,
          t_switch = t_switch,
          t_max = t_max_i,
          gap_baseline = confounder_gap_baseline,
          gap_at_switch = confounder_gap_at_switch,
          gap_end = confounder_gap_end
        )
        mean_confounder <- continuer_confounder + gap
      } else {
        mean_confounder <- continuer_confounder
      }

      # Add MVN-correlated random deviation
      confounder_matrix[i, j] <- mean_confounder + random_deviations[i, j]
    }

    # Record confounder at switch time (last measurement before switch)
    idx_before_switch <- which(time_points <= t_switch)
    if (length(idx_before_switch) > 0) {
      idx_use <- max(idx_before_switch)
      confounder_at_switch[i] <- confounder_matrix[i, idx_use]
    } else {
      confounder_at_switch[i] <- confounder_matrix[i, 1]
    }
  }

  # ============================================================================
  # GENERATE RECURRING EVENTS
  # ============================================================================
  all_events <- vector("list", n_total)

  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max_i <- T_i[i]

    # Get this person's confounder trajectory
    conf_values <- confounder_matrix[i, ]

    # Treatment effect and start time
    if (is_switcher) {
      treatment_effect <- beta_treatment
      treatment_start <- t_switch
    } else {
      treatment_effect <- 0
      treatment_start <- Inf
    }

    # Generate all events for this person
    events <- generate_piecewise_events_varying(
      t_start = 0,
      t_end = t_max_i,
      segment_hazards = hazard_schedule$hazards,
      segment_times = hazard_schedule$times,
      confounder_values = conf_values,
      confounder_times = time_points,
      beta_confounder = beta_confounder,
      treatment_effect = treatment_effect,
      treatment_start = treatment_start
    )

    # Split into pre and post switch
    pre_events <- events[events < t_switch]
    post_events <- events[events >= t_switch]

    all_events[[i]] <- list(pre = pre_events, post = post_events)
  }

  # ============================================================================
  # CREATE RESULT DATA FRAME
  # ============================================================================
  result <- data.frame(
    id = 1:n_total,
    pair_id = pair_id,
    cohort = cohort,
    T_max = T_i,
    switch_time = switch_time,
    confounder_at_switch = confounder_at_switch
  )

  result$n_pre_events <- sapply(all_events, function(x) length(x$pre))
  result$n_post_events <- sapply(all_events, function(x) length(x$post))
  result$n_total_events <- result$n_pre_events + result$n_post_events

  result$first_pre_event_time <- sapply(all_events, function(x) {
    if (length(x$pre) > 0) min(x$pre) else NA
  })

  result$first_post_event_time <- sapply(all_events, function(x) {
    if (length(x$post) > 0) min(x$post) else NA
  })

  result$has_pre_event <- result$n_pre_events > 0
  result$has_post_event <- result$n_post_events > 0

  # Store confounder at t=0 for compatibility
  result$conf_t0 <- confounder_matrix[, 1]

  # Store attributes for later use
  attr(result, "all_events") <- all_events
  attr(result, "hazard_schedule") <- hazard_schedule
  attr(result, "confounder_matrix") <- confounder_matrix
  attr(result, "confounder_times") <- time_points

  return(result)
}

# ============================================================================
# Helper: Lookup confounder at a specific time
# ============================================================================
lookup_confounder_at_time <- function(data, time_col, confounder_interval = 0.5) {
  confounder_matrix <- attr(data, "confounder_matrix")
  confounder_times <- attr(data, "confounder_times")

  if (is.null(confounder_matrix) || is.null(confounder_times)) {
    # Fall back to confounder_at_switch if no time-varying data
    return(data$confounder_at_switch)
  }

  result <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    t <- data[[time_col]][i]
    if (is.na(t)) {
      result[i] <- NA
      next
    }

    # Find the last measurement at or before this time
    idx <- which(confounder_times <= t)
    if (length(idx) > 0) {
      result[i] <- confounder_matrix[data$id[i], max(idx)]
    } else {
      result[i] <- confounder_matrix[data$id[i], 1]
    }
  }

  return(result)
}

# ============================================================================
# Helper: Re-derive events based on new switch time
# ============================================================================
rederive_events_recurring <- function(data, new_switch_col) {
  all_events <- attr(data, "all_events")
  if (is.null(all_events)) {
    stop("Data must have 'all_events' attribute")
  }

  data$new_first_pre_event_time <- NA
  data$new_first_post_event_time <- NA
  data$new_has_pre_event <- FALSE
  data$new_has_post_event <- FALSE
  data$new_post_event_time_from_switch <- NA

  for (i in 1:nrow(data)) {
    new_switch <- data[[new_switch_col]][i]
    if (is.na(new_switch)) next

    pre_events <- all_events[[data$id[i]]]$pre
    post_events <- all_events[[data$id[i]]]$post
    all_event_times <- sort(c(pre_events, post_events))

    if (length(all_event_times) == 0) next

    events_before <- all_event_times[all_event_times < new_switch]
    if (length(events_before) > 0) {
      data$new_first_pre_event_time[i] <- min(events_before)
      data$new_has_pre_event[i] <- TRUE
    }

    events_after <- all_event_times[all_event_times >= new_switch]
    if (length(events_after) > 0) {
      data$new_first_post_event_time[i] <- min(events_after)
      data$new_has_post_event[i] <- TRUE
      data$new_post_event_time_from_switch[i] <- min(events_after) - new_switch
    }
  }

  # Also lookup confounder at the new switch time
  data$confounder_at_switch <- lookup_confounder_at_time(data, new_switch_col)

  return(data)
}
