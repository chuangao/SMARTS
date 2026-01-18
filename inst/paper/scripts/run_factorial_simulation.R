# ============================================================================
# Factorial Simulation Study: SMARTS vs Baseline Approach
# ============================================================================
#
# Factors varied:
# 1. Switching interval: 0.25-0.75, 0.35-0.65, 0.45-0.55
# 2. Confounder trajectory: Variable vs Homogeneous
# 3. Hazard shape: Accelerating (1.5) vs Constant (1)
#
# Two situations compared:
# 1. Baseline: Continuers from t=0, Switchers from switch time
# 2. SMARTS: Both from (pseudo-)switching time
#
# Methods: Naive, Covariate Adj, PSM, IPTW
# ============================================================================

library(survival)
library(MatchIt)
library(dplyr)
library(SMARTS)

source("simulate_survival_confounder_v2.R")
source("validate_smarts_method.R")

# ============================================================================
# Run Single Scenario
# ============================================================================

run_scenario <- function(
  switch_start, switch_end,
  gap_baseline, gap_peak, gap_end,  # Confounder trajectory
  shape,                             # Hazard shape
  n_pairs = 5000,
  n_iter = 10,
  seed = 123
) {

  set.seed(seed)

  # Simulate data
  data <- simulate_survival_data_confounder_v2(
    n_pairs = n_pairs,
    beta_treatment = log(0.5),
    beta_confounder = log(1.3),
    lambda_0 = 10,
    shape = shape,
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

  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  # =========================================
  # SITUATION 1: BASELINE APPROACH
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
  # Naive
  cox_b_naive <- coxph(Surv(event_time, event) ~ treated, data = baseline_data)
  hr_b_naive <- exp(coef(cox_b_naive)["treated"])

  # Adjusted
  cox_b_adj <- coxph(Surv(event_time, event) ~ treated + confounder, data = baseline_data)
  hr_b_adj <- exp(coef(cox_b_adj)["treated"])

  # PSM
  tryCatch({
    match_b <- matchit(treated ~ confounder, data = baseline_data, method = "nearest", caliper = 0.1)
    data_b_matched <- match.data(match_b)
    if (nrow(data_b_matched) > 10 && length(unique(data_b_matched$subclass)) > 1) {
      cox_b_psm <- coxph(Surv(event_time, event) ~ treated + strata(subclass), data = data_b_matched)
      hr_b_psm <- exp(coef(cox_b_psm)["treated"])
    } else {
      hr_b_psm <- NA
    }
  }, error = function(e) { hr_b_psm <<- NA })

  # IPTW
  ps_b <- glm(treated ~ confounder, data = baseline_data, family = binomial)
  baseline_data$ps <- predict(ps_b, type = "response")
  baseline_data$weight <- ifelse(baseline_data$treated == 1, 1/baseline_data$ps, 1/(1-baseline_data$ps))
  cox_b_iptw <- coxph(Surv(event_time, event) ~ treated, data = baseline_data, weights = weight)
  hr_b_iptw <- exp(coef(cox_b_iptw)["treated"])

  # =========================================
  # SITUATION 2: SMARTS APPROACH (n_iter times)
  # =========================================

  smarts_results <- list()

  for (iter in 1:n_iter) {
    # Prepare for SMARTS
    switchers_smarts <- switchers
    switchers_smarts$swi_yrs <- switchers_smarts$switch_time
    switchers_smarts$fup_yrs <- switchers_smarts$T_max

    continuers_smarts <- continuers
    continuers_smarts$swi_yrs <- NA
    continuers_smarts$fup_yrs <- continuers_smarts$T_max

    smarts_input <- list(cont = continuers_smarts, swi = switchers_smarts)

    smarts_result <- random_assign(smarts_input, nbin = 10, seed = iter * 100,
                                    swi_time = "swi_yrs", cens_time = "fup_yrs")

    cont_assigned <- smarts_result$assigned$cont
    swi_assigned <- smarts_result$assigned$swi

    cont_assigned$new_switch_time <- cont_assigned$swi_yrs
    swi_assigned$new_switch_time <- swi_assigned$switch_time

    cont_assigned <- rederive_events(cont_assigned, "new_switch_time")
    swi_assigned <- rederive_events(swi_assigned, "new_switch_time")

    cont_assigned <- lookup_confounder_at_time(cont_assigned, "new_switch_time", 0.5)
    swi_assigned <- lookup_confounder_at_time(swi_assigned, "new_switch_time", 0.5)

    cols_keep <- c("id", "cohort", "new_post_event", "new_post_event_time", "new_confounder_at_switch")

    smarts_data <- rbind(cont_assigned[, cols_keep], swi_assigned[, cols_keep])
    smarts_data$treated <- as.numeric(smarts_data$cohort == "switcher")
    smarts_data <- smarts_data[smarts_data$new_post_event_time > 0, ]

    # Naive
    cox_s_naive <- coxph(Surv(new_post_event_time, new_post_event) ~ treated, data = smarts_data)
    hr_s_naive <- exp(coef(cox_s_naive)["treated"])

    # Adjusted
    cox_s_adj <- coxph(Surv(new_post_event_time, new_post_event) ~ treated + new_confounder_at_switch, data = smarts_data)
    hr_s_adj <- exp(coef(cox_s_adj)["treated"])

    # PSM
    tryCatch({
      match_s <- matchit(treated ~ new_confounder_at_switch, data = smarts_data, method = "nearest", caliper = 0.1)
      data_s_matched <- match.data(match_s)
      if (nrow(data_s_matched) > 10 && length(unique(data_s_matched$subclass)) > 1) {
        cox_s_psm <- coxph(Surv(new_post_event_time, new_post_event) ~ treated + strata(subclass), data = data_s_matched)
        hr_s_psm <- exp(coef(cox_s_psm)["treated"])
      } else {
        hr_s_psm <- NA
      }
    }, error = function(e) { hr_s_psm <<- NA })

    # IPTW
    ps_s <- glm(treated ~ new_confounder_at_switch, data = smarts_data, family = binomial)
    smarts_data$ps <- predict(ps_s, type = "response")
    smarts_data$weight <- ifelse(smarts_data$treated == 1, 1/smarts_data$ps, 1/(1-smarts_data$ps))
    cox_s_iptw <- coxph(Surv(new_post_event_time, new_post_event) ~ treated, data = smarts_data, weights = weight)
    hr_s_iptw <- exp(coef(cox_s_iptw)["treated"])

    smarts_results[[iter]] <- data.frame(
      naive = hr_s_naive,
      adj = hr_s_adj,
      psm = hr_s_psm,
      iptw = hr_s_iptw
    )
  }

  smarts_df <- do.call(rbind, smarts_results)

  # Return results
  return(data.frame(
    situation = c(rep("Baseline", 4), rep("SMARTS", 4)),
    method = rep(c("Naive", "Adj", "PSM", "IPTW"), 2),
    hr = c(hr_b_naive, hr_b_adj, hr_b_psm, hr_b_iptw,
           mean(smarts_df$naive), mean(smarts_df$adj),
           mean(smarts_df$psm, na.rm = TRUE), mean(smarts_df$iptw)),
    bias = c(hr_b_naive - 0.5, hr_b_adj - 0.5, hr_b_psm - 0.5, hr_b_iptw - 0.5,
             mean(smarts_df$naive) - 0.5, mean(smarts_df$adj) - 0.5,
             mean(smarts_df$psm, na.rm = TRUE) - 0.5, mean(smarts_df$iptw) - 0.5)
  ))
}

# ============================================================================
# Run Full Factorial Simulation
# ============================================================================

run_factorial <- function() {

  # Define factor levels
  switch_intervals <- list(
    c(0.25, 0.75),
    c(0.35, 0.65),
    c(0.45, 0.55)
  )

  confounder_trajectories <- list(
    variable = c(gap_baseline = 0.2, gap_peak = 1.5, gap_end = 0),      # Variable
    homogeneous = c(gap_baseline = 0.9, gap_peak = 1.1, gap_end = 0.9)  # Homogeneous
  )

  hazard_shapes <- c(accelerating = 1.5, constant = 1)

  all_results <- list()
  scenario_num <- 0

  for (sw in switch_intervals) {
    for (conf_name in names(confounder_trajectories)) {
      for (hz_name in names(hazard_shapes)) {
        scenario_num <- scenario_num + 1

        conf <- confounder_trajectories[[conf_name]]
        hz <- hazard_shapes[[hz_name]]

        cat("Scenario", scenario_num, ": Switch=", sw[1], "-", sw[2],
            ", Confounder=", conf_name, ", Hazard=", hz_name, "\n")

        result <- run_scenario(
          switch_start = sw[1],
          switch_end = sw[2],
          gap_baseline = conf["gap_baseline"],
          gap_peak = conf["gap_peak"],
          gap_end = conf["gap_end"],
          shape = hz,
          n_pairs = 5000,
          n_iter = 10
        )

        result$switch_interval <- paste0(sw[1], "-", sw[2])
        result$confounder_type <- conf_name
        result$hazard_type <- hz_name

        all_results[[scenario_num]] <- result
      }
    }
  }

  final_results <- do.call(rbind, all_results)
  return(final_results)
}

# ============================================================================
# Run and Display Results
# ============================================================================

if (interactive()) {
  cat("Starting Factorial Simulation Study...\n\n")

  results <- run_factorial()

  cat("\n========================================\n")
  cat("FACTORIAL SIMULATION RESULTS\n")
  cat("True HR = 0.500\n")
  cat("========================================\n\n")

  # Display results by scenario
  for (sw in unique(results$switch_interval)) {
    for (conf in unique(results$confounder_type)) {
      for (hz in unique(results$hazard_type)) {
        subset <- results[results$switch_interval == sw &
                           results$confounder_type == conf &
                           results$hazard_type == hz, ]

        cat("Switch:", sw, "| Confounder:", conf, "| Hazard:", hz, "\n")
        cat("-------------------------------------------------\n")
        cat(sprintf("%-10s %-8s %8s %8s\n", "Situation", "Method", "HR", "Bias"))

        for (i in 1:nrow(subset)) {
          cat(sprintf("%-10s %-8s %8.3f %8.3f\n",
                      subset$situation[i], subset$method[i],
                      subset$hr[i], subset$bias[i]))
        }
        cat("\n")
      }
    }
  }

  # Summary: IPTW comparison
  cat("\n========================================\n")
  cat("SUMMARY: IPTW COMPARISON (Baseline vs SMARTS)\n")
  cat("========================================\n\n")

  iptw_results <- results[results$method == "IPTW", ]
  iptw_wide <- reshape(iptw_results[, c("switch_interval", "confounder_type", "hazard_type", "situation", "hr")],
                       idvar = c("switch_interval", "confounder_type", "hazard_type"),
                       timevar = "situation", direction = "wide")

  cat(sprintf("%-12s %-12s %-12s %10s %10s\n",
              "Switch", "Confounder", "Hazard", "Baseline", "SMARTS"))
  cat("----------------------------------------------------------------\n")

  for (i in 1:nrow(iptw_wide)) {
    cat(sprintf("%-12s %-12s %-12s %10.3f %10.3f\n",
                iptw_wide$switch_interval[i],
                iptw_wide$confounder_type[i],
                iptw_wide$hazard_type[i],
                iptw_wide$hr.Baseline[i],
                iptw_wide$hr.SMARTS[i]))
  }
}
