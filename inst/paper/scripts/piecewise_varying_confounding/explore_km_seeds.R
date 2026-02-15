# ============================================================================
# Quick exploration: try multiple seeds and base hazards for KM illustration
# Pick the best-looking one, then use those parameters in plot_km_curves.R
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
# Colors and theme (same as main script)
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

create_km_plot <- function(data, title, panel_label, weighted = FALSE,
                           xlim = c(0, 4), show_ylabel = TRUE,
                           ylim = c(0.25, 1.05)) {

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
    ylim = ylim,
    break.x.by = 1,
    break.y.by = 0.25,
    xlab = "",
    ylab = ""
  )

  y_breaks <- seq(ylim[1], 1.00, by = 0.25)
  y_labels <- sprintf("%.2f", y_breaks)

  plot <- p$plot +
    theme_km() +
    ggtitle(title) +
    scale_y_continuous(
      limits = ylim,
      breaks = y_breaks,
      labels = y_labels
    ) +
    annotate("text", x = xlim[1], y = ylim[2] - 0.03, label = panel_label,
             hjust = 0, vjust = 0, size = 4.5, fontface = "bold")

  if (show_ylabel) {
    plot <- plot + ylab("Survival probability")
  } else {
    plot <- plot + theme(axis.title.y = element_blank())
  }

  return(plot)
}

# ============================================================================
# Quick function: simulate + generate full 9-panel figure for one seed
# ============================================================================

try_seed <- function(seed, true_hr, gap_baseline, gap_at_switch, gap_end,
                     base_hazard, scenario_label) {

  set.seed(seed)
  data <- simulate_piecewise_hazard(
    n_pairs = 5000,
    beta_treatment = log(true_hr),
    beta_confounder = log(2.0),
    base_hazard = base_hazard,
    hazard_trend = "increasing",
    trend_strength = 0.5,
    segment_length = 0.5,
    confounder_gap_baseline = gap_baseline,
    confounder_gap_at_switch = gap_at_switch,
    confounder_gap_end = gap_end
  )

  data$treated <- as.numeric(data$cohort == "switcher")
  data$confounder <- data$confounder_at_switch
  data$group <- factor(ifelse(data$treated == 1, "Switcher", "Continuer"),
                       levels = c("Continuer", "Switcher"))

  switchers <- data[data$cohort == "switcher", ]
  continuers <- data[data$cohort == "continuer", ]

  # A: Pre-switch (True)
  data_A <- data.frame(
    time = ifelse(data$has_pre_event, data$first_pre_event_time, data$switch_time),
    event = as.numeric(data$has_pre_event),
    group = data$group
  )

  # B: Post-switch (True)
  data_B <- data.frame(
    time = ifelse(data$has_post_event,
                  data$first_post_event_time - data$switch_time,
                  data$T_max - data$switch_time),
    event = as.numeric(data$has_post_event),
    group = data$group
  )

  # C: No pseudo-switch
  first_event_cont <- ifelse(continuers$has_pre_event,
                             continuers$first_pre_event_time,
                             ifelse(continuers$has_post_event,
                                    continuers$first_post_event_time, NA))
  data_C <- data.frame(
    time = c(
      ifelse(switchers$has_post_event,
             switchers$first_post_event_time - switchers$switch_time,
             switchers$T_max - switchers$switch_time),
      ifelse(continuers$n_total_events > 0, first_event_cont, continuers$T_max)
    ),
    event = c(as.numeric(switchers$has_post_event),
              as.numeric(continuers$n_total_events > 0)),
    group = factor(c(rep("Switcher", nrow(switchers)),
                     rep("Continuer", nrow(continuers))),
                   levels = c("Continuer", "Switcher"))
  )

  # SMARTS
  switchers_smarts <- switchers
  switchers_smarts$swi_yrs <- switchers_smarts$switch_time
  switchers_smarts$fup_yrs <- switchers_smarts$T_max
  continuers_smarts <- continuers
  continuers_smarts$swi_yrs <- NA
  continuers_smarts$fup_yrs <- continuers_smarts$T_max

  smarts_result <- random_assign(
    list(cont = continuers_smarts, swi = switchers_smarts),
    nbin = 10, seed = seed, swi_time = "swi_yrs", cens_time = "fup_yrs"
  )

  cont_assigned <- smarts_result$assigned$cont
  swi_assigned <- smarts_result$assigned$swi
  cont_assigned$new_switch_time <- cont_assigned$swi_yrs
  swi_assigned$new_switch_time <- swi_assigned$switch_time
  cont_assigned <- rederive_events_recurring(cont_assigned, "new_switch_time")
  swi_assigned <- rederive_events_recurring(swi_assigned, "new_switch_time")

  # D,E: SMARTS only
  data_D <- data.frame(
    time = c(ifelse(swi_assigned$new_has_pre_event, swi_assigned$new_first_pre_event_time, swi_assigned$new_switch_time),
             ifelse(cont_assigned$new_has_pre_event, cont_assigned$new_first_pre_event_time, cont_assigned$new_switch_time)),
    event = c(as.numeric(swi_assigned$new_has_pre_event), as.numeric(cont_assigned$new_has_pre_event)),
    group = factor(c(rep("Switcher", nrow(swi_assigned)), rep("Continuer", nrow(cont_assigned))),
                   levels = c("Continuer", "Switcher"))
  )
  data_E <- data.frame(
    time = c(ifelse(swi_assigned$new_has_post_event, swi_assigned$new_post_event_time_from_switch, swi_assigned$T_max - swi_assigned$new_switch_time),
             ifelse(cont_assigned$new_has_post_event, cont_assigned$new_post_event_time_from_switch, cont_assigned$T_max - cont_assigned$new_switch_time)),
    event = c(as.numeric(swi_assigned$new_has_post_event), as.numeric(cont_assigned$new_has_post_event)),
    group = factor(c(rep("Switcher", nrow(swi_assigned)), rep("Continuer", nrow(cont_assigned))),
                   levels = c("Continuer", "Switcher"))
  )

  # PSM
  smarts_combined <- rbind(
    data.frame(id = swi_assigned$id, treated = 1, confounder = swi_assigned$confounder_at_switch,
               new_has_pre_event = swi_assigned$new_has_pre_event, new_first_pre_event_time = swi_assigned$new_first_pre_event_time,
               new_switch_time = swi_assigned$new_switch_time, new_has_post_event = swi_assigned$new_has_post_event,
               new_post_event_time_from_switch = swi_assigned$new_post_event_time_from_switch, T_max = swi_assigned$T_max),
    data.frame(id = cont_assigned$id, treated = 0, confounder = cont_assigned$confounder_at_switch,
               new_has_pre_event = cont_assigned$new_has_pre_event, new_first_pre_event_time = cont_assigned$new_first_pre_event_time,
               new_switch_time = cont_assigned$new_switch_time, new_has_post_event = cont_assigned$new_has_post_event,
               new_post_event_time_from_switch = cont_assigned$new_post_event_time_from_switch, T_max = cont_assigned$T_max)
  )

  ps_model <- glm(treated ~ confounder, data = smarts_combined, family = binomial)
  smarts_combined$ps <- predict(ps_model, type = "response")
  match_out <- matchit(treated ~ confounder, data = smarts_combined, method = "nearest", caliper = 0.1)
  data_matched <- match.data(match_out)
  data_matched$group <- factor(ifelse(data_matched$treated == 1, "Switcher", "Continuer"),
                               levels = c("Continuer", "Switcher"))

  data_F <- data.frame(time = ifelse(data_matched$new_has_pre_event, data_matched$new_first_pre_event_time, data_matched$new_switch_time),
                       event = as.numeric(data_matched$new_has_pre_event), group = data_matched$group)
  data_G <- data.frame(time = ifelse(data_matched$new_has_post_event, data_matched$new_post_event_time_from_switch, data_matched$T_max - data_matched$new_switch_time),
                       event = as.numeric(data_matched$new_has_post_event), group = data_matched$group)

  # IPTW
  p_treated <- mean(smarts_combined$treated)
  smarts_combined$weight <- ifelse(smarts_combined$treated == 1, p_treated / smarts_combined$ps, (1 - p_treated) / (1 - smarts_combined$ps))
  smarts_combined$group <- factor(ifelse(smarts_combined$treated == 1, "Switcher", "Continuer"),
                                  levels = c("Continuer", "Switcher"))

  data_H <- data.frame(time = ifelse(smarts_combined$new_has_pre_event, smarts_combined$new_first_pre_event_time, smarts_combined$new_switch_time),
                       event = as.numeric(smarts_combined$new_has_pre_event), group = smarts_combined$group, weight = smarts_combined$weight)
  data_I <- data.frame(time = ifelse(smarts_combined$new_has_post_event, smarts_combined$new_post_event_time_from_switch, smarts_combined$T_max - smarts_combined$new_switch_time),
                       event = as.numeric(smarts_combined$new_has_post_event), group = smarts_combined$group, weight = smarts_combined$weight)

  # Generate plots
  p_A <- create_km_plot(data_A, "True simulation\npre-switch", "A", xlim = c(0, 4), show_ylabel = TRUE)
  p_B <- create_km_plot(data_B, "True simulation\npost-switch", "B", xlim = c(0, 4), show_ylabel = TRUE)
  p_C <- create_km_plot(data_C, "No pseudo-switch", "C", xlim = c(0, 6), show_ylabel = FALSE)
  p_D <- create_km_plot(data_D, "Assign pseudo-switch\npre-switch", "D", xlim = c(0, 5), show_ylabel = FALSE)
  p_E <- create_km_plot(data_E, "Assign pseudo-switch\npost-switch", "E", xlim = c(0, 5), show_ylabel = FALSE)
  p_F <- create_km_plot(data_F, "Assign pseudo-switch\npre-switch + PSM", "F", xlim = c(0, 4), show_ylabel = FALSE)
  p_G <- create_km_plot(data_G, "Assign pseudo-switch\npost-switch + PSM", "G", xlim = c(0, 4), show_ylabel = FALSE)
  p_H <- create_km_plot(data_H, "Assign pseudo-switch\npre-switch + IPTW", "H", weighted = TRUE, xlim = c(0, 4), show_ylabel = FALSE)
  p_I <- create_km_plot(data_I, "Assign pseudo-switch\npost-switch + IPTW", "I", weighted = TRUE, xlim = c(0, 4), show_ylabel = FALSE)

  layout <- c(
    area(1, 1), area(2, 1), area(1, 2, 2, 2),
    area(1, 3), area(2, 3), area(1, 4), area(2, 4), area(1, 5), area(2, 5)
  )

  main_plot <- p_A + p_B + p_C + p_D + p_E + p_F + p_G + p_H + p_I +
    plot_layout(design = layout)

  legend_plot <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) +
    geom_line(aes(color = "Continuer", linetype = "Continuer"), linewidth = 1.2, alpha = 0.8) +
    geom_line(aes(color = "Switcher", linetype = "Switcher"), linewidth = 1.2, alpha = 0.8) +
    scale_color_manual(name = NULL, values = c("Continuer" = color_continuer, "Switcher" = color_switcher),
                       guide = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
    scale_linetype_manual(name = NULL, values = c("Continuer" = "solid", "Switcher" = "dashed"), guide = "none") +
    theme_void() + theme(legend.position = "top", legend.direction = "horizontal",
                         legend.text = element_text(size = 12), legend.key.width = unit(2, "cm"))

  legend_grob <- cowplot::get_legend(legend_plot)
  legend_wrap <- wrap_elements(legend_grob)
  title_wrap <- wrap_elements(
    grid::textGrob(paste0(scenario_label, "  [seed=", seed, ", base_hazard=", base_hazard, "]"),
                   gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )

  final_plot <- title_wrap / legend_wrap / main_plot +
    plot_layout(heights = c(0.03, 0.04, 1))

  filename <- paste0("explore_km_", gsub("[^a-zA-Z0-9]", "_", scenario_label),
                     "_seed", seed, "_hz", base_hazard, ".png")
  ggsave(filename, final_plot, width = 14, height = 5.5, units = "in", dpi = 200)
  cat("  Saved:", filename, "\n")
}

# ============================================================================
# Try combinations
# ============================================================================

seeds <- c(42, 99, 123, 200, 314, 777)

cat("\n=== Scenario 1: Sicker switchers, treatment helps (HR=0.7) ===\n")
cat("Trying base_hazard = 0.005 and 0.007\n\n")
for (hz in c(0.005, 0.007)) {
  for (s in seeds) {
    try_seed(seed = s, true_hr = 0.7,
             gap_baseline = 0.5, gap_at_switch = 1.5, gap_end = 0.8,
             base_hazard = hz, scenario_label = "Sicker helps")
  }
}

cat("\n=== Scenario 2: Healthier switchers, treatment harms (HR=1.5) ===\n")
cat("Trying base_hazard = 0.005 and 0.007\n\n")
for (hz in c(0.005, 0.007)) {
  for (s in seeds) {
    try_seed(seed = s, true_hr = 1.5,
             gap_baseline = -0.3, gap_at_switch = -1.5, gap_end = -0.5,
             base_hazard = hz, scenario_label = "Healthier harms")
  }
}

cat("\nDone! Review the PNGs and pick the best seed + hazard for each scenario.\n")
