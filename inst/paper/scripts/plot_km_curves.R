# ============================================================================
# Plot Kaplan-Meier Curves for Treatment Switching Study
# ============================================================================
# Publication-ready KM curves using ggsurvplot
# Layout: A,B (True), C (No pseudo-switch), D,E (pseudo-switch),
#         F,G (PSM), H,I (IPTW)
# ============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(MatchIt)
library(patchwork)
library(cowplot)

source("simulate_survival_confounder.R")

# Working directory set by caller
# ============================================================================
# Simulate Data
# ============================================================================

set.seed(123)
data <- simulate_survival_data_confounder(
  n_pairs = 5000,
  beta_treatment = log(0.5),
  beta_confounder = log(1.3),
  lambda_0 = 10,
  shape = 1.5,
  T_min = 2,
  T_max = 6,
  switch_start = 0.25,
  switch_end = 0.75,
  confounder_interval = 0.5,
  confounder_baseline_mean = 2.5,
  confounder_gap_baseline = 0.5,
  confounder_gap_peak = 1.6,
  confounder_gap_end = 0.8,
  confounder_sd = 0.8
)

# Add treatment indicator and confounder column
data$treated <- as.numeric(data$cohort == "switcher")
data$confounder <- data$confounder_at_switch
data$group <- factor(ifelse(data$treated == 1, "Switcher", "Continuer"),
                     levels = c("Continuer", "Switcher"))

# ============================================================================
# Prepare Datasets
# ============================================================================

# --- Dataset A: Pre-switching (True simulation) ---
data_A <- data.frame(
  time = data$pre_event_time,
  event = data$pre_event,
  group = data$group
)

# --- Dataset B: Post-switching (True simulation) ---
data_B <- data.frame(
  time = data$post_event_time - data$switch_time,
  event = data$post_event,
  group = data$group
)

# --- Dataset C: No pseudo-switch (misaligned) ---
data_C <- data.frame(
  time = ifelse(data$treated == 1,
                ifelse(data$post_event == 1,
                       data$post_event_time - data$switch_time,
                       data$T_max - data$switch_time),
                ifelse(data$pre_event == 1,
                       data$pre_event_time,
                       ifelse(data$post_event == 1,
                              data$post_event_time,
                              data$T_max))),
  event = ifelse(data$treated == 1,
                 data$post_event,
                 pmax(data$pre_event, data$post_event)),
  group = data$group
)

# --- Datasets D, E: Assign pseudo-switch (same as A, B for this simulation) ---
data_D <- data_A
data_E <- data_B

# ============================================================================
# Propensity Score Matching
# ============================================================================

ps_model <- glm(treated ~ confounder, data = data, family = binomial)
data$ps <- predict(ps_model, type = "response")

match_out <- matchit(treated ~ confounder, data = data,
                     method = "nearest", caliper = 0.1)
data_matched <- match.data(match_out)
data_matched$group <- factor(ifelse(data_matched$cohort == "switcher",
                                    "Switcher", "Continuer"),
                             levels = c("Continuer", "Switcher"))

# --- Dataset F: Pre-switching (PSM) ---
data_F <- data.frame(
  time = data_matched$pre_event_time,
  event = data_matched$pre_event,
  group = data_matched$group
)

# --- Dataset G: Post-switching (PSM) ---
data_G <- data.frame(
  time = data_matched$post_event_time - data_matched$switch_time,
  event = data_matched$post_event,
  group = data_matched$group
)

# ============================================================================
# IPTW Weights
# ============================================================================

p_treated <- mean(data$treated)
data$weight <- ifelse(data$treated == 1,
                      p_treated / data$ps,
                      (1 - p_treated) / (1 - data$ps))

# --- Dataset H: Pre-switching (IPTW) ---
data_H <- data.frame(
  time = data$pre_event_time,
  event = data$pre_event,
  group = data$group,
  weight = data$weight
)

# --- Dataset I: Post-switching (IPTW) ---
data_I <- data.frame(
  time = data$post_event_time - data$switch_time,
  event = data$post_event,
  group = data$group,
  weight = data$weight
)

# ============================================================================
# Custom Theme and Colors (matching reference image)
# ============================================================================

# Colors: Orange for Continuer, Blue for Switcher (matching boxplot colors)
color_continuer <- "#E69F00"
color_switcher <- "#0072B2"
colors <- c(color_continuer, color_switcher)

# Publication-ready theme with larger fonts, no top/right borders, with grid
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

  if (weighted && "weight" %in% names(data)) {
    fit <- survfit(Surv(time, event) ~ group, data = data, weights = weight)
  } else {
    fit <- survfit(Surv(time, event) ~ group, data = data)
  }

  # Use ggsurvplot with explicit color mapping and line types
  p <- ggsurvplot(
    fit,
    data = data,
    palette = c(color_continuer, color_switcher),
    linetype = c("solid", "dashed"),  # Continuer = solid, Switcher = dashed
    legend.labs = c("Continuer", "Switcher"),
    size = 1.0,
    alpha = 0.8,  # Add transparency to match boxplot appearance
    censor = FALSE,
    legend = "none",
    xlim = xlim,
    ylim = c(0.25, 1.05),
    break.x.by = 1,
    break.y.by = 0.25,
    xlab = "",
    ylab = ""
  )

  # Extract and customize plot
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

cat("Generating plots...\n")

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
# Column 1: A (top), B (bottom)
# Column 2: C (spans full height)
# Column 3: D (top), E (bottom)
# Column 4: F (top), G (bottom)
# Column 5: H (top), I (bottom)

# Use explicit layout with area()
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

# Combine plots in layout order: A, B, C, D, E, F, G, H, I
main_plot <- p_A + p_B + p_C + p_D + p_E + p_F + p_G + p_H + p_I +
  plot_layout(design = layout)

# Create a standalone legend as a ggplot
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

# Extract just the legend using cowplot
legend_grob <- cowplot::get_legend(legend_plot)

# Wrap legend as a patchwork element
legend_wrap <- wrap_elements(legend_grob)

# Combine: legend on top, main plot below
final_plot <- legend_wrap / main_plot +
  plot_layout(heights = c(0.05, 1))

# ============================================================================
# Save Plot
# ============================================================================

ggsave("../output/km_curves_publication.pdf", final_plot,
       width = 14, height = 5.5, units = "in", dpi = 300)

ggsave("../output/km_curves_publication.png", final_plot,
       width = 14, height = 5.5, units = "in", dpi = 300)

cat("\nPlots saved to:\n")
cat("  - ../output/km_curves_publication.pdf\n")
cat("  - ../output/km_curves_publication.png (300 dpi)\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")
cat("True HR = 0.5, Confounder HR = 1.3\n")
cat("\nSample sizes:\n")
cat("  All subjects:", nrow(data), "\n")
cat("  PSM matched: ", nrow(data_matched), "\n")
