# ============================================================================
# Visualize Piecewise Constant Hazard Schedules
# ============================================================================

source("simulate_piecewise_hazard.R")

# Create hazard schedules for each trend type
max_time <- 6

constant <- create_hazard_schedule(max_time, segment_length = 0.5,
                                    base_hazard = 0.02,
                                    hazard_trend = "constant", trend_strength = 0)

increasing <- create_hazard_schedule(max_time, segment_length = 0.5,
                                      base_hazard = 0.02,
                                      hazard_trend = "increasing", trend_strength = 0.3)

decreasing <- create_hazard_schedule(max_time, segment_length = 0.5,
                                      base_hazard = 0.02,
                                      hazard_trend = "decreasing", trend_strength = 0.15)

# Print hazard values at key time points
cat("================================================================\n")
cat("HAZARD SCHEDULES\n")
cat("================================================================\n\n")

cat("Time points:", head(constant$times, 13), "\n\n")

cat("Constant hazards:   ", round(constant$hazards, 4), "\n")
cat("Increasing hazards: ", round(increasing$hazards, 4), "\n")
cat("Decreasing hazards: ", round(decreasing$hazards, 4), "\n")

# Calculate ratio of early to late hazard
cat("\n================================================================\n")
cat("HAZARD RATIO (time=5 vs time=1)\n")
cat("================================================================\n\n")

# Time 1 is roughly segment 2, time 5 is roughly segment 10
cat("Constant:   ", round(constant$hazards[10] / constant$hazards[2], 2), "\n")
cat("Increasing: ", round(increasing$hazards[10] / increasing$hazards[2], 2), "\n")
cat("Decreasing: ", round(decreasing$hazards[10] / decreasing$hazards[2], 2), "\n")

# Plot
pdf("hazard_schedules.pdf", width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# Helper to plot step function
plot_hazard <- function(schedule, title, col = "steelblue") {
  times <- schedule$times
  hazards <- schedule$hazards
  n <- length(hazards)

  plot(NA, xlim = c(0, max_time), ylim = c(0, max(hazards) * 1.2),
       xlab = "Time", ylab = "Hazard", main = title)

  for (i in 1:n) {
    segments(times[i], hazards[i], times[i+1], hazards[i], col = col, lwd = 2)
    if (i < n) {
      segments(times[i+1], hazards[i], times[i+1], hazards[i+1], col = col, lwd = 1, lty = 2)
    }
  }

  # Add average switch time
  abline(v = 2, col = "red", lty = 2)
  text(2.1, max(hazards) * 1.1, "Avg switch\ntime", adj = 0, cex = 0.8, col = "red")
}

plot_hazard(constant, "Constant Hazard")
plot_hazard(increasing, "Increasing Hazard")
plot_hazard(decreasing, "Decreasing Hazard")

dev.off()

cat("\nPlot saved to hazard_schedules.pdf\n")

# ============================================================================
# Demonstrate the problem with Baseline method
# ============================================================================
cat("\n================================================================\n")
cat("WHY BASELINE METHOD IS BIASED\n")
cat("================================================================\n\n")

cat("With INCREASING hazard:\n")
cat("- Continuers' events counted from time 0: low hazard period\n")
cat("- Switchers' events counted from switch time (~t=2): higher hazard period\n")
cat("- Baseline method attributes more switcher events to treatment, not time\n")
cat("- Result: Overestimates treatment HR\n\n")

cat("With DECREASING hazard:\n")
cat("- Continuers' events counted from time 0: high hazard period\n")
cat("- Switchers' events counted from switch time (~t=2): lower hazard period\n")
cat("- Baseline method attributes fewer switcher events to treatment effect\n")
cat("- Result: Underestimates treatment HR\n\n")

cat("SMARTS method aligns pseudo-switch times, so both groups are\n")
cat("compared in similar calendar time windows, reducing this bias.\n")
