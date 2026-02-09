# ============================================================================
# Plot Kaplan-Meier Curves for Piecewise Constant Hazard Simulation
# ============================================================================
# Scenario: Sicker switchers, treatment helps (gap=0.5, HR=0.7)
# Using INCREASING hazard trend to demonstrate SMARTS advantage
#
# Layout: A,B (True), C (No pseudo-switch), D,E (pseudo-switch),
#         F,G (PSM), H,I (IPTW)
# ============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(MatchIt)
library(MASS)
library(cowplot)
library(patchwork)
library(SMARTS)
library(dplyr)

source("simulate_piecewise_hazard.R")

# ============================================================================
# Simulate Data
# ============================================================================

set.seed(123)
data <- simulate_piecewise_hazard(
  n_pairs = 5000,
  beta_treatment = log(0.7),
  beta_confounder = log(2.0),
  base_hazard = 0.01,             # Lower base hazard
  hazard_trend = "increasing",
  trend_strength = 0.5,           # Stronger trend for visible variation
  segment_length = 0.5,
  confounder_gap = 0.5  # Sicker switchers
)

# Add treatment indicator and confounder column
data$treated <- as.numeric(data$cohort == "switcher")
data$confounder <- data$confounder_at_switch
data$group <- factor(ifelse(data$treated == 1, "Switcher", "Continuer"),
                     levels = c("Continuer", "Switcher"))

switchers <- data[data$cohort == "switcher", ]
continuers <- data[data$cohort == "continuer", ]

cat("=== Simulation Summary ===\n")
cat("Scenario: Sicker switchers, treatment helps\n")
cat("True HR:", 0.7, "| Confounder HR:", 2.0, "| Confounder gap:", 0.5, "\n")
cat("Hazard trend: Increasing\n")
cat("Switchers:", nrow(switchers), "| Continuers:", nrow(continuers), "\n")

# ============================================================================
# Prepare Datasets
# ============================================================================

# --- Dataset A: Pre-switching (True simulation) ---
data_A <- data.frame(
  time = ifelse(data$has_pre_event, data$first_pre_event_time, data$switch_time),
  event = as.numeric(data$has_pre_event),
  group = data$group
)

# --- Dataset B: Post-switching (True simulation) ---
data_B <- data.frame(
  time = ifelse(data$has_post_event,
                data$first_post_event_time - data$switch_time,
                data$T_max - data$switch_time),
  event = as.numeric(data$has_post_event),
  group = data$group
)

# --- Dataset C: No pseudo-switch (misaligned Baseline) ---
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
                 levels = c("Continuer", "Switcher"))
)

# ============================================================================
# SMARTS: Assign pseudo-switch times
# ============================================================================

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

# --- Datasets D, E: Assign pseudo-switch ---
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
                 levels = c("Continuer", "Switcher"))
)

# ============================================================================
# Propensity Score Matching
# ============================================================================

# Combine assigned data for PSM
smarts_combined <- rbind(
  data.frame(
    id = swi_assigned$id,
    treated = 1,
    confounder = swi_assigned$confounder_at_switch,
    new_has_pre_event = swi_assigned$new_has_pre_event,
    new_first_pre_event_time = swi_assigned$new_first_pre_event_time,
    new_switch_time = swi_assigned$new_switch_time,
    new_has_post_event = swi_assigned$new_has_post_event,
    new_post_event_time_from_switch = swi_assigned$new_post_event_time_from_switch,
    T_max = swi_assigned$T_max
  ),
  data.frame(
    id = cont_assigned$id,
    treated = 0,
    confounder = cont_assigned$confounder_at_switch,
    new_has_pre_event = cont_assigned$new_has_pre_event,
    new_first_pre_event_time = cont_assigned$new_first_pre_event_time,
    new_switch_time = cont_assigned$new_switch_time,
    new_has_post_event = cont_assigned$new_has_post_event,
    new_post_event_time_from_switch = cont_assigned$new_post_event_time_from_switch,
    T_max = cont_assigned$T_max
  )
)

ps_model <- glm(treated ~ confounder, data = smarts_combined, family = binomial)
smarts_combined$ps <- predict(ps_model, type = "response")

match_out <- matchit(treated ~ confounder, data = smarts_combined,
                     method = "nearest", caliper = 0.1)
data_matched <- match.data(match_out)
data_matched$group <- factor(ifelse(data_matched$treated == 1, "Switcher", "Continuer"),
                             levels = c("Continuer", "Switcher"))

cat("  PSM matched:", nrow(data_matched), "\n")

# --- Dataset F: Pre-switching (PSM) ---
data_F <- data.frame(
  time = ifelse(data_matched$new_has_pre_event,
                data_matched$new_first_pre_event_time,
                data_matched$new_switch_time),
  event = as.numeric(data_matched$new_has_pre_event),
  group = data_matched$group
)

# --- Dataset G: Post-switching (PSM) ---
data_G <- data.frame(
  time = ifelse(data_matched$new_has_post_event,
                data_matched$new_post_event_time_from_switch,
                data_matched$T_max - data_matched$new_switch_time),
  event = as.numeric(data_matched$new_has_post_event),
  group = data_matched$group
)

# ============================================================================
# IPTW Weights
# ============================================================================

p_treated <- mean(smarts_combined$treated)
smarts_combined$weight <- ifelse(smarts_combined$treated == 1,
                                  p_treated / smarts_combined$ps,
                                  (1 - p_treated) / (1 - smarts_combined$ps))

smarts_combined$group <- factor(ifelse(smarts_combined$treated == 1, "Switcher", "Continuer"),
                                levels = c("Continuer", "Switcher"))

# --- Dataset H: Pre-switching (IPTW) ---
data_H <- data.frame(
  time = ifelse(smarts_combined$new_has_pre_event,
                smarts_combined$new_first_pre_event_time,
                smarts_combined$new_switch_time),
  event = as.numeric(smarts_combined$new_has_pre_event),
  group = smarts_combined$group,
  weight = smarts_combined$weight
)

# --- Dataset I: Post-switching (IPTW) ---
data_I <- data.frame(
  time = ifelse(smarts_combined$new_has_post_event,
                smarts_combined$new_post_event_time_from_switch,
                smarts_combined$T_max - smarts_combined$new_switch_time),
  event = as.numeric(smarts_combined$new_has_post_event),
  group = smarts_combined$group,
  weight = smarts_combined$weight
)

# ============================================================================
# Custom Theme and Colors
# ============================================================================

color_continuer <- "#E69F00"
color_switcher <- "#0072B2"
colors <- c(color_continuer, color_switcher)

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

create_km_plot <- function(data, title, panel_label, weighted = FALSE,
                           xlim = c(0, 4), show_ylabel = TRUE) {

  data <- data[data$time > 0 & !is.na(data$time), ]

  if (weighted && "weight" %in% names(data)) {
    fit <- survfit(Surv(time, event) ~ group, data = data, weights = weight)
  } else {
    fit <- survfit(Surv(time, event) ~ group, data = data)
  }

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
    ylim = c(0.25, 1.05),
    break.x.by = 1,
    break.y.by = 0.25,
    xlab = "",
    ylab = ""
  )

  plot <- p$plot +
    theme_km() +
    ggtitle(title) +
    scale_y_continuous(
      limits = c(0.25, 1.05),
      breaks = c(0.25, 0.50, 0.75, 1.00),
      labels = c("0.25", "0.50", "0.75", "1.00")
    ) +
    annotate("text", x = xlim[1], y = 1.02, label = panel_label,
             hjust = 0, vjust = 0, size = 4.5, fontface = "bold")

  if (show_ylabel) {
    plot <- plot + ylab("Survival probability")
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

p_C <- create_km_plot(data_C, "No pseudo-switch", "C",
                      xlim = c(0, 6), show_ylabel = FALSE)

p_D <- create_km_plot(data_D, "Assign pseudo-switch\npre-switch", "D",
                      xlim = c(0, 5), show_ylabel = FALSE)

p_F <- create_km_plot(data_F, "Assign pseudo-switch\npre-switch + PSM", "F",
                      xlim = c(0, 4), show_ylabel = FALSE)

p_H <- create_km_plot(data_H, "Assign pseudo-switch\npre-switch + IPTW", "H",
                      weighted = TRUE, xlim = c(0, 4), show_ylabel = FALSE)

# Row 2: Post-switch plots
p_B <- create_km_plot(data_B, "True simulation\npost-switch", "B",
                      xlim = c(0, 4), show_ylabel = TRUE)

p_E <- create_km_plot(data_E, "Assign pseudo-switch\npost-switch", "E",
                      xlim = c(0, 5), show_ylabel = FALSE)

p_G <- create_km_plot(data_G, "Assign pseudo-switch\npost-switch + PSM", "G",
                      xlim = c(0, 4), show_ylabel = FALSE)

p_I <- create_km_plot(data_I, "Assign pseudo-switch\npost-switch + IPTW", "I",
                      weighted = TRUE, xlim = c(0, 4), show_ylabel = FALSE)

# ============================================================================
# Combine Plots using patchwork
# ============================================================================

# Layout: C (No pseudo-switch) spans both rows in column 2
layout <- c(
  area(1, 1),       # A: row 1, col 1
  area(2, 1),       # B: row 2, col 1
  area(1, 2, 2, 2), # C: rows 1-2, col 2 (spans both rows)
  area(1, 3),       # D: row 1, col 3
  area(2, 3),       # E: row 2, col 3
  area(1, 4),       # F: row 1, col 4
  area(2, 4),       # G: row 2, col 4
  area(1, 5),       # H: row 1, col 5
  area(2, 5)        # I: row 2, col 5
)

main_plot <- p_A + p_B + p_C + p_D + p_E + p_F + p_G + p_H + p_I +
  plot_layout(design = layout)

# Create a standalone legend
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
       width = 14, height = 5.5, units = "in", dpi = 300)

ggsave("km_curves_piecewise.png", final_plot,
       width = 14, height = 5.5, units = "in", dpi = 300)

cat("\nPlots saved to:\n")
cat("  - km_curves_piecewise.pdf\n")
cat("  - km_curves_piecewise.png (300 dpi)\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")
cat("True HR = 0.7, Confounder HR = 2.0, Hazard trend = Increasing\n")
cat("\nSample sizes:\n")
cat("  All subjects:", nrow(data), "\n")
cat("  SMARTS assigned:", nrow(smarts_combined), "\n")
cat("  PSM matched:", nrow(data_matched), "\n")
