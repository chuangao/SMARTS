# ============================================================================
# Plot Kaplan-Meier Curves for Piecewise Constant Hazard Simulation
# ============================================================================
# Publication-ready KM curves showing:
# - Row 1: Pre-switch period
# - Row 2: Post-switch period
# Columns: True simulation, Baseline method, SMARTS method
# ============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(MASS)
library(cowplot)
library(patchwork)  # Load after MASS so patchwork::area is not masked
library(SMARTS)
library(dplyr)

source("simulate_piecewise_hazard.R")

# ============================================================================
# Simulate Data with Increasing Hazard Trend
# ============================================================================

set.seed(123)
data <- simulate_piecewise_hazard(

  n_pairs = 3000,
  beta_treatment = log(0.7),       # True HR = 0.7
  beta_confounder = log(2.0),
  base_hazard = 0.02,
  hazard_trend = "increasing",
  trend_strength = 0.3,
  segment_length = 0.5,
  confounder_gap = 0.5             # Switchers sicker (higher confounder)
)

# Add group labels
data$group <- factor(ifelse(data$cohort == "switcher", "Switcher", "Continuer"),
                     levels = c("Continuer", "Switcher"))
data$treated <- as.numeric(data$cohort == "switcher")

switchers <- data[data$cohort == "switcher", ]
continuers <- data[data$cohort == "continuer", ]

cat("=== Simulation Summary ===\n")
cat("True HR:", 0.7, "\n")
cat("Hazard trend: Increasing\n")
cat("Switchers:", nrow(switchers), "| Post-event rate:", round(mean(switchers$has_post_event)*100, 1), "%\n")
cat("Continuers:", nrow(continuers), "| Post-event rate:", round(mean(continuers$has_post_event)*100, 1), "%\n")

# ============================================================================
# Prepare Datasets for KM Curves
# ============================================================================

# --- Dataset A: True Pre-switch (time to first pre-switch event) ---
data_A <- data.frame(
  time = ifelse(data$has_pre_event, data$first_pre_event_time, data$switch_time),
  event = as.numeric(data$has_pre_event),
  group = data$group
)

# --- Dataset B: True Post-switch (time from switch to first post-switch event) ---
data_B <- data.frame(
  time = ifelse(data$has_post_event,
                data$first_post_event_time - data$switch_time,
                data$T_max - data$switch_time),
  event = as.numeric(data$has_post_event),
  group = data$group
)

# --- Dataset C: Baseline method (misaligned - continuers from time 0) ---
# Switchers: post-switch events from switch time
# Continuers: any event from time 0
first_event_cont <- ifelse(continuers$has_pre_event,
                           continuers$first_pre_event_time,
                           ifelse(continuers$has_post_event,
                                  continuers$first_post_event_time, NA))

data_C <- data.frame(
  time = c(
    ifelse(switchers$has_post_event,
           switchers$first_post_event_time - switchers$switch_time,
           switchers$T_max - switchers$switch_time),
    ifelse(continuers$n_total_events > 0,
           first_event_cont,
           continuers$T_max)
  ),
  event = c(
    as.numeric(switchers$has_post_event),
    as.numeric(continuers$n_total_events > 0)
  ),
  group = factor(c(rep("Switcher", nrow(switchers)),
                   rep("Continuer", nrow(continuers))),
                 levels = c("Continuer", "Switcher")),
  confounder = c(switchers$confounder_at_switch,
                 continuers$conf_t0)
)

# --- SMARTS: Assign pseudo-switch times to continuers ---
switchers_smarts <- switchers
switchers_smarts$swi_yrs <- switchers_smarts$switch_time
switchers_smarts$fup_yrs <- switchers_smarts$T_max

continuers_smarts <- continuers
continuers_smarts$swi_yrs <- NA
continuers_smarts$fup_yrs <- continuers_smarts$T_max

smarts_input <- list(cont = continuers_smarts, swi = switchers_smarts)

smarts_result <- random_assign(smarts_input, nbin = 10, seed = 123,
                               swi_time = "swi_yrs", cens_time = "fup_yrs")

cont_assigned <- smarts_result$assigned$cont
swi_assigned <- smarts_result$assigned$swi

cat("\nSMARTS assignment:\n")
cat("  Assigned continuers:", nrow(cont_assigned), "\n")
cat("  Assigned switchers:", nrow(swi_assigned), "\n")

# Set new switch times
cont_assigned$new_switch_time <- cont_assigned$swi_yrs
swi_assigned$new_switch_time <- swi_assigned$switch_time

# Re-derive events based on pseudo-switch times
cont_assigned <- rederive_events_recurring(cont_assigned, "new_switch_time")
swi_assigned <- rederive_events_recurring(swi_assigned, "new_switch_time")

# --- Dataset D: SMARTS Pre-switch ---
data_D <- data.frame(
  time = c(
    ifelse(swi_assigned$new_has_pre_event,
           swi_assigned$new_first_pre_event_time,
           swi_assigned$new_switch_time),
    ifelse(cont_assigned$new_has_pre_event,
           cont_assigned$new_first_pre_event_time,
           cont_assigned$new_switch_time)
  ),
  event = c(
    as.numeric(swi_assigned$new_has_pre_event),
    as.numeric(cont_assigned$new_has_pre_event)
  ),
  group = factor(c(rep("Switcher", nrow(swi_assigned)),
                   rep("Continuer", nrow(cont_assigned))),
                 levels = c("Continuer", "Switcher"))
)

# --- Dataset E: SMARTS Post-switch ---
data_E <- data.frame(
  time = c(
    ifelse(swi_assigned$new_has_post_event,
           swi_assigned$new_post_event_time_from_switch,
           swi_assigned$T_max - swi_assigned$new_switch_time),
    ifelse(cont_assigned$new_has_post_event,
           cont_assigned$new_post_event_time_from_switch,
           cont_assigned$T_max - cont_assigned$new_switch_time)
  ),
  event = c(
    as.numeric(swi_assigned$new_has_post_event),
    as.numeric(cont_assigned$new_has_post_event)
  ),
  group = factor(c(rep("Switcher", nrow(swi_assigned)),
                   rep("Continuer", nrow(cont_assigned))),
                 levels = c("Continuer", "Switcher")),
  confounder = c(swi_assigned$confounder_at_switch,
                 cont_assigned$confounder_at_switch)
)

# Filter out invalid times
data_D <- data_D[data_D$time > 0, ]
data_E <- data_E[data_E$time > 0, ]

# ============================================================================
# Custom Theme and Colors
# ============================================================================

color_continuer <- "#E69F00"
color_switcher <- "#0072B2"

theme_km <- function() {
  theme_pubr(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = "none",
      plot.margin = margin(5, 10, 5, 5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

# ============================================================================
# Plot Function
# ============================================================================

create_km_plot <- function(data, title, panel_label,
                           xlim = c(0, 5), show_ylabel = TRUE) {

  fit <- survfit(Surv(time, event) ~ group, data = data)

  p <- ggsurvplot(
    fit,
    data = data,
    palette = c(color_continuer, color_switcher),
    linetype = c("solid", "dashed"),
    legend.labs = c("Continuer", "Switcher"),
    size = 1.0,
    alpha = 0.8,
    censor = FALSE,
    legend = "none",
    xlim = xlim,
    ylim = c(0.4, 1.05),
    break.x.by = 1,
    break.y.by = 0.2,
    xlab = "",
    ylab = ""
  )

  plot <- p$plot +
    theme_km() +
    ggtitle(title) +
    scale_y_continuous(
      limits = c(0.4, 1.05),
      breaks = c(0.4, 0.6, 0.8, 1.0),
      labels = c("0.40", "0.60", "0.80", "1.00")
    ) +
    annotate("text", x = xlim[1], y = 1.02, label = panel_label,
             hjust = 0, vjust = 0, size = 4.5, fontface = "bold")

  if (show_ylabel) {
    plot <- plot + ylab("Event-free probability")
  } else {
    plot <- plot + theme(axis.title.y = element_blank())
  }

  return(plot)
}

# ============================================================================
# Generate All Plots
# ============================================================================

cat("\nGenerating plots...\n")

# Row 1: Pre-switch plots
p_A <- create_km_plot(data_A, "True simulation\npre-switch", "A",
                      xlim = c(0, 4), show_ylabel = TRUE)

p_D <- create_km_plot(data_D, "SMARTS\npre-switch", "D",
                      xlim = c(0, 4), show_ylabel = FALSE)

# Row 2: Post-switch plots
p_B <- create_km_plot(data_B, "True simulation\npost-switch", "B",
                      xlim = c(0, 5), show_ylabel = TRUE)

p_C <- create_km_plot(data_C, "Baseline method\n(misaligned)", "C",
                      xlim = c(0, 6), show_ylabel = FALSE)

p_E <- create_km_plot(data_E, "SMARTS\npost-switch", "E",
                      xlim = c(0, 5), show_ylabel = FALSE)

# ============================================================================
# Combine Plots
# ============================================================================

# Layout: 2 rows x 3 columns
# Row 1: A (True pre), placeholder, D (SMARTS pre)
# Row 2: B (True post), C (Baseline), E (SMARTS post)

layout <- c(
  area(1, 1),  # A
  area(1, 2),  # placeholder (will use spacer)
  area(1, 3),  # D
  area(2, 1),  # B
  area(2, 2),  # C
  area(2, 3)   # E
)

# Create a spacer plot
spacer <- ggplot() + theme_void()

main_plot <- p_A + spacer + p_D + p_B + p_C + p_E +
  plot_layout(design = layout)

# Create legend
legend_plot <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) +
  geom_line(aes(color = "Continuer", linetype = "Continuer"), linewidth = 1.2, alpha = 0.8) +
  geom_line(aes(color = "Switcher", linetype = "Switcher"), linewidth = 1.2, alpha = 0.8) +
  scale_color_manual(
    name = NULL,
    values = c("Continuer" = color_continuer, "Switcher" = color_switcher),
    guide = guide_legend(override.aes = list(linetype = c("solid", "dashed")))
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Continuer" = "solid", "Switcher" = "dashed"),
    guide = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.key.width = unit(2, "cm")
  )

legend_grob <- cowplot::get_legend(legend_plot)
legend_wrap <- wrap_elements(legend_grob)

final_plot <- legend_wrap / main_plot +
  plot_layout(heights = c(0.05, 1))

# ============================================================================
# Save Plot
# ============================================================================

ggsave("km_curves_piecewise.pdf", final_plot,
       width = 12, height = 6, units = "in", dpi = 300)

ggsave("km_curves_piecewise.png", final_plot,
       width = 12, height = 6, units = "in", dpi = 300)

cat("\nPlots saved to:\n")
cat("  - km_curves_piecewise.pdf\n")
cat("  - km_curves_piecewise.png\n")

# ============================================================================
# Summary: Hazard Ratios from Each Method
# ============================================================================

cat("\n=== Hazard Ratio Estimates ===\n")
cat("True HR: 0.7\n\n")

# Baseline method
baseline_data <- data_C
baseline_data <- baseline_data[baseline_data$time > 0, ]
baseline_data$treated <- as.numeric(baseline_data$group == "Switcher")

cox_baseline_naive <- coxph(Surv(time, event) ~ treated, data = baseline_data)
cox_baseline_adj <- coxph(Surv(time, event) ~ treated + confounder, data = baseline_data)

cat("Baseline (misaligned):\n")
cat("  Naive:    ", round(exp(coef(cox_baseline_naive)["treated"]), 3), "\n")
cat("  Adjusted: ", round(exp(coef(cox_baseline_adj)["treated"]), 3), "\n")

# SMARTS method
smarts_data <- data_E
smarts_data <- smarts_data[smarts_data$time > 0 & !is.na(smarts_data$time), ]
smarts_data$treated <- as.numeric(smarts_data$group == "Switcher")

cox_smarts_naive <- coxph(Surv(time, event) ~ treated, data = smarts_data)
cox_smarts_adj <- coxph(Surv(time, event) ~ treated + confounder, data = smarts_data)

cat("\nSMARTS (aligned):\n")
cat("  Naive:    ", round(exp(coef(cox_smarts_naive)["treated"]), 3), "\n")
cat("  Adjusted: ", round(exp(coef(cox_smarts_adj)["treated"]), 3), "\n")
