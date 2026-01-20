# ============================================================================
# Plot Factorial Simulation Results - Switch Interval × Confounder Gap Design
# Boxplot with Individual Points (matching main factorial plot style)
# ============================================================================

library(ggplot2)
library(dplyr)

# Load raw results from output folder
raw <- readRDS("../output/factorial_results_raw_vary_switch_interval_confounder_gap.rds")

# Set factor levels for proper ordering
raw$method <- factor(raw$method, levels = c("Naive", "Adj", "PSM", "IPTW", "SMR"))
raw$situation <- factor(raw$situation, levels = c("Baseline", "SMARTS"))
raw$hazard_type <- factor(raw$hazard_type,
                          levels = c("accelerating", "constant"),
                          labels = c("Accelerating Hazard", "Constant Hazard"))
raw$confounder_type <- factor(raw$confounder_type,
                              levels = c("variable", "homogeneous"),
                              labels = c("Variable Confounder", "Homogeneous Confounder"))

# ============================================================================
# Faceted Boxplot with Individual Points (no fill, border only)
# ============================================================================

p <- ggplot(raw, aes(x = method, y = mean_hr, color = situation)) +
  # Reference line at true HR
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Boxplot with no fill
  geom_boxplot(outlier.shape = NA, fill = NA, position = position_dodge(width = 0.75)) +
  # Individual points (open shapes: circle for Baseline, triangle for SMARTS)
  geom_point(aes(shape = situation),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 1.2, alpha = 0.6, stroke = 0.5, fill = NA) +
  # Facet by scenario
  facet_grid(hazard_type ~ switch_interval + confounder_type) +
  # Colors (matching KM curves and main factorial plot)
  scale_color_manual(values = c("Baseline" = "#E69F00", "SMARTS" = "#0072B2")) +
  scale_shape_manual(values = c("Baseline" = 1, "SMARTS" = 2)) +  # 1 = open circle, 2 = open triangle
  # Labels
  labs(x = "Method", y = "Hazard Ratio",
       color = "Approach", shape = "Approach",
       title = "Estimated HR by Scenario and Method",
       subtitle = "Dashed line = True HR (0.5). Each point = 1 dataset.") +
  # Theme
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 9)) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 2))),
         shape = "none")

ggsave("../output/plot_factorial_boxplot_vary_switch_interval_confounder_gap.png", p, width = 14, height = 8, dpi = 300)
ggsave("../output/plot_factorial_boxplot_vary_switch_interval_confounder_gap.pdf", p, width = 14, height = 8)
cat("Saved: ../output/plot_factorial_boxplot_vary_switch_interval_confounder_gap.png (300 dpi)\n")
cat("Saved: ../output/plot_factorial_boxplot_vary_switch_interval_confounder_gap.pdf\n")
