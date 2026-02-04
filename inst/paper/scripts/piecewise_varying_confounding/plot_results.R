# ============================================================================
# Plot Results: 4 Scenarios x 3 Hazard Trends
# Boxplot comparing Baseline vs SMARTS across multiple simulation runs
# ============================================================================

library(survival)
library(MASS)
library(SMARTS)
library(dplyr)
library(ggplot2)
library(parallel)

source("simulate_piecewise_hazard.R")

# ============================================================================
# Analysis Functions
# ============================================================================

run_smarts_analysis <- function(data, n_assignments = 50, n_cores = 4) {
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

    cont_data <- data.frame(
      treated = 0,
      event = as.numeric(cont_assigned$new_has_post_event),
      time = ifelse(cont_assigned$new_has_post_event,
                    cont_assigned$new_post_event_time_from_switch,
                    cont_assigned$T_max - cont_assigned$new_switch_time),
      confounder = cont_assigned$confounder_at_switch
    )

    swi_data <- data.frame(
      treated = 1,
      event = as.numeric(swi_assigned$new_has_post_event),
      time = ifelse(swi_assigned$new_has_post_event,
                    swi_assigned$new_post_event_time_from_switch,
                    swi_assigned$T_max - swi_assigned$new_switch_time),
      confounder = swi_assigned$confounder_at_switch
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
# Scenarios and Trends
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
# Run Multiple Simulations
# ============================================================================

n_simulations <- 20
all_results <- data.frame()

cat("Running", n_simulations, "simulations per scenario...\n")

for (h in hazard_trends) {
  cat("\nHazard trend:", h$name, "\n")

  for (s in clinical_scenarios) {
    cat("  Scenario:", s$name, "")

    for (sim in 1:n_simulations) {
      set.seed(1000 * which(sapply(hazard_trends, function(x) x$name) == h$name) +
               100 * which(sapply(clinical_scenarios, function(x) x$name) == s$name) +
               sim)

      data <- simulate_piecewise_hazard(
        n_pairs = 2000,
        beta_treatment = log(s$true_hr),
        beta_confounder = log(2.0),
        base_hazard = 0.02,
        hazard_trend = h$trend,
        trend_strength = h$strength,
        segment_length = 0.5,
        confounder_gap = s$gap
      )

      hr_baseline <- run_baseline_analysis(data)
      hr_smarts <- run_smarts_analysis(data, n_assignments = 50, n_cores = 4)

      all_results <- rbind(all_results, data.frame(
        hazard_trend = h$name,
        scenario = s$name,
        true_hr = s$true_hr,
        simulation = sim,
        method = "Baseline",
        hr = hr_baseline
      ))

      all_results <- rbind(all_results, data.frame(
        hazard_trend = h$name,
        scenario = s$name,
        true_hr = s$true_hr,
        simulation = sim,
        method = "SMARTS",
        hr = hr_smarts
      ))

      cat(".")
    }
    cat("\n")
  }
}

# Save raw results
saveRDS(all_results, "simulation_results_raw.rds")
cat("\nResults saved to simulation_results_raw.rds\n")

# ============================================================================
# Create Plot
# ============================================================================

# Calculate bias
all_results$bias <- all_results$hr - all_results$true_hr

# Set factor levels
all_results$method <- factor(all_results$method, levels = c("Baseline", "SMARTS"))
all_results$hazard_trend <- factor(all_results$hazard_trend,
                                    levels = c("Constant", "Increasing", "Decreasing"))
all_results$scenario <- factor(all_results$scenario,
                                levels = c("Sicker, tx harmful", "Sicker, tx helps",
                                           "Healthier, stay well", "Healthier, tx harmful"))

# Colors
colors <- c("Baseline" = "#E69F00", "SMARTS" = "#0072B2")

# Create boxplot
p <- ggplot(all_results, aes(x = scenario, y = hr, color = method)) +
  # Reference line at true HR (will vary by facet, so we add separately)
  geom_hline(data = data.frame(
    scenario = rep(levels(all_results$scenario), 3),
    hazard_trend = rep(levels(all_results$hazard_trend), each = 4),
    true_hr = rep(c(1.3, 0.7, 0.7, 1.3), 3)
  ), aes(yintercept = true_hr), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Boxplot
  geom_boxplot(outlier.shape = NA, fill = NA, position = position_dodge(width = 0.75)) +
  # Individual points
  geom_point(aes(shape = method),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 1.2, alpha = 0.6, stroke = 0.5, fill = NA) +
  # Facet by hazard trend
  facet_wrap(~ hazard_trend, nrow = 1) +
  # Colors and shapes
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("Baseline" = 1, "SMARTS" = 2)) +
  # Labels
  labs(x = "Scenario", y = "Estimated Hazard Ratio",
       color = "Method", shape = "Method",
       title = "SMARTS vs Baseline: Piecewise Constant Hazard Simulation",
       subtitle = "Dashed line = True HR. Each point = 1 simulation run.") +
  # Theme
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 2))),
         shape = "none")

ggsave("plot_results_boxplot.pdf", p, width = 12, height = 6)
ggsave("plot_results_boxplot.png", p, width = 12, height = 6, dpi = 300)

cat("\nPlots saved:\n")
cat("  - plot_results_boxplot.pdf\n")
cat("  - plot_results_boxplot.png\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n\n")

summary_stats <- all_results %>%
  group_by(hazard_trend, method) %>%
  summarise(
    mean_hr = round(mean(hr, na.rm = TRUE), 3),
    mean_bias = round(mean(bias, na.rm = TRUE), 3),
    mean_abs_bias = round(mean(abs(bias), na.rm = TRUE), 3),
    sd_hr = round(sd(hr, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_stats)
