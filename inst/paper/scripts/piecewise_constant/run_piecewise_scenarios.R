# ============================================================================
# 4 Clinical Scenarios x 3 Hazard Trends (Piecewise Constant)
# ============================================================================

library(survival)
library(MASS)
library(SMARTS)
library(dplyr)
library(parallel)

source("simulate_piecewise_hazard.R")

# ============================================================================
# Analysis Functions
# ============================================================================
run_smarts_analysis <- function(data, n_assignments = 100, n_cores = 4) {
  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  run_single <- function(seed) {
    switchers_smarts <- switchers
    switchers_smarts$swi_yrs <- switchers_smarts$switch_time
    switchers_smarts$fup_yrs <- switchers_smarts$T_max

    continuers_smarts <- continuers
    continuers_smarts$swi_yrs <- NA
    continuers_smarts$fup_yrs <- continuers_smarts$T_max

    smarts_input <- list(cont = continuers_smarts, swi = switchers_smarts)

    smarts_result <- tryCatch({
      random_assign(smarts_input, nbin = 10, seed = seed,
                    swi_time = "swi_yrs", cens_time = "fup_yrs")
    }, error = function(e) NULL)

    if (is.null(smarts_result)) return(NA)

    cont_assigned <- smarts_result$assigned$cont
    swi_assigned <- smarts_result$assigned$swi
    if (nrow(cont_assigned) == 0 || nrow(swi_assigned) == 0) return(NA)

    cont_assigned$new_switch_time <- cont_assigned$swi_yrs
    swi_assigned$new_switch_time <- swi_assigned$switch_time

    cont_assigned <- rederive_events_recurring(cont_assigned, "new_switch_time")
    swi_assigned <- rederive_events_recurring(swi_assigned, "new_switch_time")

    cont_assigned$new_confounder <- lookup_confounder(cont_assigned, "new_switch_time", 0.5)
    swi_assigned$new_confounder <- lookup_confounder(swi_assigned, "new_switch_time", 0.5)

    cont_data <- data.frame(
      treated = 0,
      event = as.numeric(cont_assigned$new_has_post_event),
      time = ifelse(cont_assigned$new_has_post_event,
                    cont_assigned$new_post_event_time_from_switch,
                    cont_assigned$T_max - cont_assigned$new_switch_time),
      confounder = cont_assigned$new_confounder
    )

    swi_data <- data.frame(
      treated = 1,
      event = as.numeric(swi_assigned$new_has_post_event),
      time = ifelse(swi_assigned$new_has_post_event,
                    swi_assigned$new_post_event_time_from_switch,
                    swi_assigned$T_max - swi_assigned$new_switch_time),
      confounder = swi_assigned$new_confounder
    )

    smarts_data <- rbind(cont_data, swi_data)
    smarts_data <- smarts_data[smarts_data$time > 0 & !is.na(smarts_data$time), ]
    if (nrow(smarts_data) < 50) return(NA)

    tryCatch({
      cox <- coxph(Surv(time, event) ~ treated + confounder, data = smarts_data)
      as.numeric(exp(coef(cox)["treated"]))
    }, error = function(e) NA)
  }

  results <- mclapply(1:n_assignments, run_single, mc.cores = n_cores)
  mean(unlist(results), na.rm = TRUE)
}

run_baseline_analysis <- function(data) {
  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  baseline_swi <- data.frame(
    treated = 1,
    event = as.numeric(switchers$has_post_event),
    time = ifelse(switchers$has_post_event,
                  switchers$first_post_event_time - switchers$switch_time,
                  switchers$T_max - switchers$switch_time),
    confounder = switchers$confounder_at_switch
  )

  first_event_cont <- ifelse(continuers$has_pre_event,
                              continuers$first_pre_event_time,
                              ifelse(continuers$has_post_event,
                                     continuers$first_post_event_time, NA))
  baseline_cont <- data.frame(
    treated = 0,
    event = as.numeric(continuers$n_total_events > 0),
    time = ifelse(continuers$n_total_events > 0, first_event_cont, continuers$T_max),
    confounder = continuers$conf_t0
  )

  baseline_data <- rbind(baseline_swi, baseline_cont)
  baseline_data <- baseline_data[baseline_data$time > 0, ]

  tryCatch({
    cox <- coxph(Surv(time, event) ~ treated + confounder, data = baseline_data)
    as.numeric(exp(coef(cox)["treated"]))
  }, error = function(e) NA)
}

# ============================================================================
# Scenarios
# ============================================================================
clinical_scenarios <- list(
  list(name = "Sicker, tx harmful", gap = 0.5, true_hr = 1.3),
  list(name = "Sicker, tx helps", gap = 0.5, true_hr = 0.7),
  list(name = "Healthier, stay well", gap = -0.5, true_hr = 0.7),
  list(name = "Healthier, tx harmful", gap = -0.5, true_hr = 1.3)
)

hazard_trends <- list(
  list(name = "Constant", trend = "constant", strength = 0),
  list(name = "Increasing", trend = "increasing", strength = 0.3),
  list(name = "Decreasing", trend = "decreasing", strength = 0.15)
)

# ============================================================================
# Main
# ============================================================================
cat("================================================================\n")
cat("PIECEWISE CONSTANT HAZARD: 4 SCENARIOS x 3 TRENDS\n")
cat("================================================================\n\n")

all_results <- data.frame()

for (h in hazard_trends) {
  cat("\n================================================================\n")
  cat("HAZARD TREND:", h$name, "\n")
  cat("================================================================\n")

  for (i in 1:4) {
    s <- clinical_scenarios[[i]]

    set.seed(123 + i)

    data <- simulate_piecewise_hazard(
      n_pairs = 3000,
      beta_treatment = log(s$true_hr),
      beta_confounder = log(2.0),
      base_hazard = 0.02,
      hazard_trend = h$trend,
      trend_strength = h$strength,
      segment_length = 0.5,
      confounder_gap = s$gap
    )

    event_rate <- round(mean(data$has_post_event) * 100, 1)

    hr_baseline <- run_baseline_analysis(data)
    hr_smarts <- run_smarts_analysis(data, n_assignments = 100, n_cores = 8)

    winner <- ifelse(abs(hr_smarts - s$true_hr) < abs(hr_baseline - s$true_hr), "SMARTS", "Baseline")

    cat(s$name, "| True:", s$true_hr,
        "| Event%:", event_rate,
        "| Base:", round(hr_baseline, 2), "(", round(hr_baseline - s$true_hr, 2), ")",
        "| SMARTS:", round(hr_smarts, 2), "(", round(hr_smarts - s$true_hr, 2), ")",
        "|", winner, "\n")

    all_results <- rbind(all_results, data.frame(
      hazard_trend = h$name,
      scenario = s$name,
      true_hr = s$true_hr,
      event_rate = event_rate,
      baseline_hr = round(hr_baseline, 3),
      smarts_hr = round(hr_smarts, 3),
      baseline_bias = hr_baseline - s$true_hr,
      smarts_bias = hr_smarts - s$true_hr,
      winner = winner
    ))
  }
}

# ============================================================================
# Summary
# ============================================================================
cat("\n================================================================\n")
cat("SUMMARY BY HAZARD TREND\n")
cat("================================================================\n\n")

summary_df <- all_results %>%
  group_by(hazard_trend) %>%
  summarise(
    smarts_wins = sum(winner == "SMARTS"),
    baseline_mean_abs_bias = round(mean(abs(baseline_bias)), 3),
    smarts_mean_abs_bias = round(mean(abs(smarts_bias)), 3),
    .groups = "drop"
  )

print(summary_df)

cat("\n================================================================\n")
cat("FULL RESULTS\n")
cat("================================================================\n\n")

print(all_results %>%
  mutate(baseline_bias = round(baseline_bias, 3), smarts_bias = round(smarts_bias, 3)) %>%
  select(hazard_trend, scenario, true_hr, event_rate, baseline_hr, smarts_hr, winner))
