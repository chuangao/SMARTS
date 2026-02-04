# ============================================================================
# Plot Factorial Simulation Results - Piecewise Constant Hazard
# Two Clinical Scenarios x 3 Hazard Trends
# ============================================================================

library(ggplot2)
library(dplyr)

# Load raw results
raw <- readRDS("factorial_results_raw.rds")

# Set factor levels for proper ordering
raw$method <- factor(raw$method, levels = c("Naive", "Adj", "PSM", "IPTW", "SMR"))
raw$situation <- factor(raw$situation, levels = c("Baseline", "SMARTS"))
raw$hazard_trend <- factor(raw$hazard_trend,
                           levels = c("Constant", "Increasing", "Decreasing"))

# Create nice labels for clinical scenarios
raw$clinical_scenario_label <- ifelse(
  raw$clinical_scenario == "Sicker_Helps",
  "Sicker switchers, treatment helps (HR=0.7)",
  "Healthier switchers, treatment harms (HR=1.5)"
)
raw$clinical_scenario_label <- factor(
  raw$clinical_scenario_label,
  levels = c("Sicker switchers, treatment helps (HR=0.7)",
             "Healthier switchers, treatment harms (HR=1.5)")
)

# ============================================================================
# Faceted Boxplot - Both Scenarios
# ============================================================================

p <- ggplot(raw, aes(x = method, y = hr, color = situation)) +
  # Reference line at true HR (varies by scenario)
  geom_hline(aes(yintercept = true_hr), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Boxplot with no fill
  geom_boxplot(outlier.shape = NA, fill = NA, position = position_dodge(width = 0.75)) +
  # Individual points
  geom_point(aes(shape = situation),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 1.0, alpha = 0.5, stroke = 0.4, fill = NA) +
  # Facet by scenario (rows) and hazard trend (columns)
  facet_grid(clinical_scenario_label ~ hazard_trend, scales = "free_y") +
  # Colors
  scale_color_manual(values = c("Baseline" = "#E69F00", "SMARTS" = "#0072B2")) +
  scale_shape_manual(values = c("Baseline" = 1, "SMARTS" = 2)) +
  # Labels
  labs(x = "Method", y = "Hazard Ratio",
       color = "Approach", shape = "Approach",
       title = "Estimated HR by Clinical Scenario, Hazard Trend, and Method",
       subtitle = "Dashed line = True HR. Each point = 1 simulation.") +
  # Theme
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"),
        strip.text.y = element_text(size = 8),
        plot.title = element_text(size = 12, face = "bold")) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 2))),
         shape = "none")

ggsave("plot_factorial_boxplot.png", p, width = 12, height = 8, dpi = 300)
ggsave("plot_factorial_boxplot.pdf", p, width = 12, height = 8)

cat("Saved: plot_factorial_boxplot.png (300 dpi)\n")
cat("Saved: plot_factorial_boxplot.pdf\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n\n")

summary_stats <- raw %>%
  group_by(clinical_scenario, true_hr, hazard_trend, situation, method) %>%
  summarise(
    mean_hr = round(mean(hr, na.rm = TRUE), 3),
    bias = round(mean(hr, na.rm = TRUE) - first(true_hr), 3),
    abs_bias = round(abs(mean(hr, na.rm = TRUE) - first(true_hr)), 3),
    sd = round(sd(hr, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_stats, n = 50)

# Summary by approach
cat("\n=== Mean Absolute Bias by Approach ===\n\n")

approach_summary <- raw %>%
  group_by(clinical_scenario, hazard_trend, situation) %>%
  summarise(
    mean_abs_bias = round(mean(abs(hr - true_hr), na.rm = TRUE), 3),
    .groups = "drop"
  )

print(approach_summary)

# ============================================================================
# Separate plots for each scenario (optional)
# ============================================================================

for (scenario in unique(raw$clinical_scenario)) {
  scenario_data <- raw[raw$clinical_scenario == scenario, ]
  true_hr_val <- unique(scenario_data$true_hr)

  p_scenario <- ggplot(scenario_data, aes(x = method, y = hr, color = situation)) +
    geom_hline(yintercept = true_hr_val, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, fill = NA, position = position_dodge(width = 0.75)) +
    geom_point(aes(shape = situation),
               position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
               size = 1.2, alpha = 0.6, stroke = 0.5, fill = NA) +
    facet_wrap(~ hazard_trend, nrow = 1) +
    scale_color_manual(values = c("Baseline" = "#E69F00", "SMARTS" = "#0072B2")) +
    scale_shape_manual(values = c("Baseline" = 1, "SMARTS" = 2)) +
    labs(x = "Method", y = "Hazard Ratio",
         color = "Approach", shape = "Approach",
         title = paste("Estimated HR:", unique(scenario_data$clinical_scenario_label)),
         subtitle = paste("Dashed line = True HR (", true_hr_val, "). Each point = 1 simulation.", sep = "")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold"),
          plot.title = element_text(size = 12, face = "bold")) +
    guides(color = guide_legend(override.aes = list(shape = c(1, 2))),
           shape = "none")

  filename <- paste0("plot_factorial_boxplot_", scenario, ".png")
  ggsave(filename, p_scenario, width = 10, height = 6, dpi = 300)
  cat("Saved:", filename, "\n")
}
