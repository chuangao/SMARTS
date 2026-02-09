# ============================================================================
# Plot Factorial Simulation Results - Piecewise Constant Hazard
# Boxplot with Individual Points
# ============================================================================

library(ggplot2)
library(dplyr)

# Load raw results
raw <- readRDS("factorial_results_raw.rds")

# Set factor levels for proper ordering
raw$method <- factor(raw$method, levels = c("Naive", "Adj", "PSM", "IPTW"))
raw$situation <- factor(raw$situation, levels = c("Baseline", "SMARTS"))
raw$hazard_trend <- factor(raw$hazard_trend,
                           levels = c("Constant", "Increasing", "Decreasing"))

# True HR
true_hr <- 0.7

# ============================================================================
# Faceted Boxplot with Individual Points
# ============================================================================

p <- ggplot(raw, aes(x = method, y = hr, color = situation)) +
  # Reference line at true HR
  geom_hline(yintercept = true_hr, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Boxplot with no fill
  geom_boxplot(outlier.shape = NA, fill = NA, position = position_dodge(width = 0.75)) +
  # Individual points
  geom_point(aes(shape = situation),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 1.2, alpha = 0.6, stroke = 0.5, fill = NA) +
  # Facet by hazard trend
  facet_wrap(~ hazard_trend, nrow = 1) +
  # Colors
  scale_color_manual(values = c("Baseline" = "#E69F00", "SMARTS" = "#0072B2")) +
  scale_shape_manual(values = c("Baseline" = 1, "SMARTS" = 2)) +
  # Labels
  labs(x = "Method", y = "Hazard Ratio",
       color = "Approach", shape = "Approach",
       title = "Estimated HR by Hazard Trend and Method",
       subtitle = "Scenario: Sicker switchers, treatment helps. Dashed line = True HR (0.7). Each point = 1 simulation.") +
  # Theme
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 2))),
         shape = "none")

ggsave("plot_factorial_boxplot.png", p, width = 10, height = 6, dpi = 300)
ggsave("plot_factorial_boxplot.pdf", p, width = 10, height = 6)

cat("Saved: plot_factorial_boxplot.png (300 dpi)\n")
cat("Saved: plot_factorial_boxplot.pdf\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n\n")

summary_stats <- raw %>%
  group_by(hazard_trend, situation, method) %>%
  summarise(
    mean_hr = round(mean(hr, na.rm = TRUE), 3),
    bias = round(mean(hr, na.rm = TRUE) - true_hr, 3),
    abs_bias = round(abs(mean(hr, na.rm = TRUE) - true_hr), 3),
    sd = round(sd(hr, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_stats, n = 30)

# Summary by approach
cat("\n=== Mean Absolute Bias by Approach ===\n\n")

approach_summary <- raw %>%
  group_by(hazard_trend, situation) %>%
  summarise(
    mean_abs_bias = round(mean(abs(hr - true_hr), na.rm = TRUE), 3),
    .groups = "drop"
  )

print(approach_summary)
