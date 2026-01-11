# ============================================================================
# Verify Confounder Gap Pattern Over Time
# ============================================================================

# Source the simulation to get the helper function
source("simulate_survival_confounder.R")

# Parameters matching the simulation
confounder_baseline_mean <- 2.5
confounder_gap_baseline <- 0.5
confounder_gap_peak <- 1.6
confounder_gap_floor <- 0.8
t_switch <- 2.0
t_max <- 6.0

# Calculate expected gap at various time points
times <- seq(0, 6, by = 0.5)
gaps <- sapply(times, function(t) {
  calculate_gap(t, t_switch, t_max,
                confounder_gap_baseline,
                confounder_gap_peak,
                confounder_gap_floor)
})

cat("=== Expected Confounder Gap Pattern ===\n\n")
cat("Switching time:", t_switch, "\n")
cat("Max follow-up:", t_max, "\n\n")

cat(sprintf("%8s %12s %15s %15s %12s\n",
            "Time", "Gap", "Switcher_mean", "Continuer_mean", "Period"))
cat(sprintf("%8s %12s %15s %15s %12s\n",
            "----", "---", "-------------", "---------------", "------"))

for (i in 1:length(times)) {
  t <- times[i]
  gap <- gaps[i]
  switcher_mean <- confounder_baseline_mean + gap
  continuer_mean <- confounder_baseline_mean
  period <- ifelse(t < t_switch, "Pre",
                   ifelse(t == t_switch, "Switch", "Post"))

  cat(sprintf("%8.1f %12.2f %15.2f %15.2f %12s\n",
              t, gap, switcher_mean, continuer_mean, period))
}

cat("\n=== Gap Pattern Description ===\n")
cat("Pre-switching (0.0 to 2.0):\n")
cat("  - Gap increases from", confounder_gap_baseline, "to", confounder_gap_peak, "\n")
cat("  - Switchers deteriorate faster (confounding by indication)\n\n")

cat("Post-switching (2.0 to 6.0):\n")
cat("  - Gap decreases from", confounder_gap_peak, "to", confounder_gap_floor, "\n")
cat("  - Switchers improve on new treatment but remain worse\n")
cat("  - Floor prevents switchers from becoming better than continuers\n\n")

# Now verify with actual simulated data
cat("\n=== Verification with Simulated Data ===\n\n")

# Extract confounder columns
confounder_cols <- grep("^confounder_[0-9]", names(data), value = TRUE)

# Calculate mean confounder by cohort at each time point
for (col in confounder_cols[1:min(13, length(confounder_cols))]) {  # First 13 time points
  time <- as.numeric(gsub("confounder_", "", col))

  switcher_mean <- mean(data[data$cohort == "switcher", col], na.rm = TRUE)
  continuer_mean <- mean(data[data$cohort == "continuer", col], na.rm = TRUE)
  observed_gap <- switcher_mean - continuer_mean

  cat(sprintf("Time %.1f: Switcher=%.2f, Continuer=%.2f, Gap=%.2f\n",
              time, switcher_mean, continuer_mean, observed_gap))
}

cat("\nNote: Observed gaps vary due to:\n")
cat("  1. Random variation in confounder scores (SD = 0.8)\n")
cat("  2. Variable switching times across pairs\n")
cat("  3. Different follow-up times (not all patients observed at all times)\n")
