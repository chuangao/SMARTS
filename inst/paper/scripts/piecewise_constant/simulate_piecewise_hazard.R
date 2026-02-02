# ============================================================================
# Simulate Recurring Events with Piecewise Constant Hazard
# ============================================================================
#
# Time is divided into small segments, each with its own constant hazard.
# Hazard can gradually increase or decrease across segments.
# Both switchers and continuers experience the same time-varying baseline.
# Treatment effect applies only to switchers after their switch time.
#
# ============================================================================

# ============================================================================
# Generate events with piecewise constant hazard
# ============================================================================
generate_piecewise_events <- function(t_start, t_end, segment_hazards,
                                       segment_times, linear_pred,
                                       treatment_effect = 0,
                                       treatment_start = Inf) {
  # segment_hazards: vector of hazard rates for each segment
  # segment_times: vector of segment boundaries (including 0 and max time)
  # treatment_effect: log HR for treatment (applied after treatment_start)
  # treatment_start: time when treatment effect begins (for switchers)

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

    # Apply confounder effect
    hazard <- base_hazard * exp(linear_pred)

    # Apply treatment effect if we're past treatment start
    if (t_current >= treatment_start) {
      hazard <- hazard * exp(treatment_effect)
    }

    # Generate next event time using exponential with this hazard
    # hazard = 1/scale, so scale = 1/hazard
    if (hazard > 0) {
      inter_arrival <- rexp(1, rate = hazard)
    } else {
      inter_arrival <- Inf
    }

    t_next <- t_current + inter_arrival

    # Check if event occurs in this segment
    if (t_next < seg_end && t_next < t_end) {
      # Check if we cross treatment boundary
      if (t_current < treatment_start && t_next >= treatment_start) {
        # Event crosses treatment boundary - need to handle carefully
        # Simulate remaining time in pre-treatment period
        time_to_boundary <- treatment_start - t_current
        prob_event_before <- 1 - exp(-hazard * time_to_boundary)

        if (runif(1) < prob_event_before) {
          # Event before treatment starts
          t_event <- t_current + rexp(1, rate = hazard)
          t_event <- min(t_event, treatment_start - 1e-10)
          if (t_event < t_end) {
            event_times <- c(event_times, t_event)
          }
          t_current <- t_event
        } else {
          # No event before treatment, continue from boundary
          t_current <- treatment_start
        }
      } else {
        # Normal case - event in same treatment period
        event_times <- c(event_times, t_next)
        t_current <- t_next
      }
    } else {
      # No event in this segment, move to next segment
      t_current <- seg_end
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
  # hazard_trend: "constant", "increasing", "decreasing"
  # trend_strength: how much hazard changes per time unit

  segment_times <- seq(0, max_time + segment_length, by = segment_length)
  n_segments <- length(segment_times) - 1

  # Calculate hazard for each segment based on trend
  segment_hazards <- numeric(n_segments)

  for (i in 1:n_segments) {
    mid_time <- (segment_times[i] + segment_times[i + 1]) / 2

    if (hazard_trend == "constant") {
      segment_hazards[i] <- base_hazard
    } else if (hazard_trend == "increasing") {
      # Hazard increases linearly with time
      segment_hazards[i] <- base_hazard * (1 + trend_strength * mid_time)
    } else if (hazard_trend == "decreasing") {
      # Hazard decreases linearly with time (but stays positive)
      segment_hazards[i] <- base_hazard * max(0.1, 1 - trend_strength * mid_time)
    }
  }

  return(list(times = segment_times, hazards = segment_hazards))
}

# ============================================================================
# Main Simulation Function
# ============================================================================
simulate_piecewise_hazard <- function(
    n_pairs,
    beta_treatment,            # Treatment effect (log HR)
    beta_confounder,           # Confounder effect (log HR per unit)

    # Hazard parameters
    base_hazard = 0.01,        # Base hazard rate
    hazard_trend = "constant", # "constant", "increasing", "decreasing"
    trend_strength = 0.1,      # How much hazard changes over time
    segment_length = 0.5,      # Length of each piecewise segment

    # Follow-up parameters
    T_min = 2,
    T_max = 6,
    switch_start = 0.25,
    switch_end = 0.75,

    # Confounder parameters
    confounder_baseline_mean = 2.5,
    confounder_sd = 0.8,
    confounder_gap = 0.5       # Gap between switcher and continuer confounders
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

  # Generate confounder values (simplified - single baseline value per person)
  # Switchers have confounder_gap added to their baseline
  confounder <- numeric(n_total)
  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    gap <- if (is_switcher) confounder_gap else 0
    confounder[i] <- confounder_baseline_mean + gap + rnorm(1, 0, confounder_sd)
  }

  # Generate recurring events
  all_events <- vector("list", n_total)

  for (i in 1:n_total) {
    is_switcher <- cohort[i] == "switcher"
    t_switch <- switch_time[i]
    t_max_i <- T_i[i]

    # Linear predictor from confounder
    linear_pred <- beta_confounder * confounder[i]

    # Treatment effect and start time
    if (is_switcher) {
      treatment_effect <- beta_treatment
      treatment_start <- t_switch
    } else {
      treatment_effect <- 0
      treatment_start <- Inf  # Continuers never get treatment
    }

    # Generate all events for this person
    events <- generate_piecewise_events(
      t_start = 0,
      t_end = t_max_i,
      segment_hazards = hazard_schedule$hazards,
      segment_times = hazard_schedule$times,
      linear_pred = linear_pred,
      treatment_effect = treatment_effect,
      treatment_start = treatment_start
    )

    # Split into pre and post switch
    pre_events <- events[events < t_switch]
    post_events <- events[events >= t_switch]

    all_events[[i]] <- list(pre = pre_events, post = post_events)
  }

  # Create result data frame
  result <- data.frame(
    id = 1:n_total,
    pair_id = pair_id,
    cohort = cohort,
    T_max = T_i,
    switch_time = switch_time,
    confounder_at_switch = confounder  # Same as baseline in this simplified model
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

  # Store confounder as conf_t0 for compatibility
  result$conf_t0 <- confounder

  attr(result, "all_events") <- all_events
  attr(result, "hazard_schedule") <- hazard_schedule

  return(result)
}

# ============================================================================
# Helper functions for SMARTS analysis
# ============================================================================
lookup_confounder <- function(data, time_col, confounder_interval = 0.5) {
  # Simplified version - confounder doesn't change over time
  return(data$confounder_at_switch)
}

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

  return(data)
}
