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
library(cowplot)

source("simulate_survival_confounder.R")

setwd("/Users/cg253/iCloud/com\~apple\~CloudDocs/Projects/Athens_method")
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
  confounder_gap_floor = 0.8,
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

# Colors: Yellow for Continuer, Blue for Switcher
color_continuer <- "#E6A817"
color_switcher <- "#1E88E5"
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

  # Use ggsurvplot with explicit color mapping
  p <- ggsurvplot(
    fit,
    data = data,
    palette = c(color_continuer, color_switcher),
    legend.labs = c("Continuer", "Switcher"),
    size = 0.7,
    censor = FALSE,
    legend = "none",
    xlim = xlim,
    ylim = c(0, 1.05),
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
      limits = c(0, 1.05),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = c("0", "0.25", "0.50", "0.75", "1.00")
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

# Empty placeholder for C column row 2
p_empty <- ggplot() + theme_void()

# ============================================================================
# Create Legend
# ============================================================================

legend_df <- data.frame(
  x = c(1, 2, 1, 2),
  y = c(1, 1, 2, 2),
  group = factor(rep(c("Continuer", "Switcher"), each = 2),
                 levels = c("Continuer", "Switcher"))
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = group, group = group)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Continuer" = color_continuer,
                                "Switcher" = color_switcher),
                     name = "") +
  theme_void() +
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        legend.key.width = unit(2, "cm"))

legend <- get_legend(legend_plot)

# ============================================================================
# Combine Plots
# ============================================================================

# Row 1
row1 <- plot_grid(p_A, p_C, p_D, p_F, p_H,
                  nrow = 1, rel_widths = c(1, 1, 1, 1, 1),
                  align = "h", axis = "tb")

# Row 2
row2 <- plot_grid(p_B, p_empty, p_E, p_G, p_I,
                  nrow = 1, rel_widths = c(1, 1, 1, 1, 1),
                  align = "h", axis = "tb")

# Common x-axis label
xlab <- ggdraw() +
  draw_label("Time", size = 14, hjust = 0.5)

# Combine all
final_plot <- plot_grid(
  legend,
  row1,
  row2,
  xlab,
  ncol = 1,
  rel_heights = c(0.06, 1, 1, 0.06)
)

# ============================================================================
# Save Plot
# ============================================================================

ggsave("../output/km_curves_publication.pdf", final_plot,
       width = 14, height = 5.5, units = "in", dpi = 300)

ggsave("../output/km_curves_publication.png", final_plot,
       width = 14, height = 5.5, units = "in", dpi = 300)

cat("\nPlots saved to:\n")
cat("  - ../output/km_curves_publication.pdf\n")
cat("  - ../output/km_curves_publication.png\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")
cat("True HR = 0.5, Confounder HR = 1.3\n")
cat("\nSample sizes:\n")
cat("  All subjects:", nrow(data), "\n")
cat("  PSM matched: ", nrow(data_matched), "\n")
