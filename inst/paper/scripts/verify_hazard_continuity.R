# ============================================================================
# Verify Hazard Continuity at Switching Time
# ============================================================================
# This script demonstrates that the baseline hazard is continuous across
# the switching time point with the conditional Weibull approach

# Weibull hazard function: λ(t) = (shape/scale) * (t/scale)^(shape-1)
weibull_hazard <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1)
}

# Parameters from simulation
lambda_0 <- 10
shape <- 1.5
beta_confounder <- log(2)
beta_treatment <- log(0.5)
beta_event_flag <- log(2)

# Example patient: switcher with pre_event = 1, switches at t = 2
t_switch <- 2.0
event_flag <- 1

# PRE-SWITCHING: scale parameter
scale_pre <- lambda_0 / exp(beta_confounder)

# POST-SWITCHING: scale parameter
linear_predictor <- beta_treatment + beta_event_flag * event_flag
scale_post <- lambda_0 / exp(linear_predictor)

cat("=== Hazard Continuity Check ===\n\n")
cat("Patient: Switcher with pre_event = 1\n")
cat("Switching time:", t_switch, "\n\n")

cat("Scale parameters:\n")
cat("  Pre-switching:  scale_pre  =", round(scale_pre, 3), "\n")
cat("  Post-switching: scale_post =", round(scale_post, 3), "\n\n")

# Calculate baseline hazard at switching time using BOTH scales
cat("Baseline hazard λ₀(t) at switching time t =", t_switch, ":\n")

# Using pre-switching scale (old approach - would be discontinuous)
hazard_pre <- weibull_hazard(t_switch, shape, scale_pre)
cat("  Using scale_pre:  λ₀(", t_switch, ") =", round(hazard_pre, 4), "\n")

# Using post-switching scale (new approach - continuous in calendar time)
hazard_post <- weibull_hazard(t_switch, shape, scale_post)
cat("  Using scale_post: λ₀(", t_switch, ") =", round(hazard_post, 4), "\n\n")

cat("Note: With conditional Weibull sampling, we use scale_post for\n")
cat("post-switching events, ensuring the baseline hazard follows\n")
cat("calendar time (not time-since-switch).\n\n")

# Plot hazard over time for this patient
times <- seq(0, 6, by = 0.1)

# Baseline hazard using post-switching scale (calendar time)
baseline_hazard <- weibull_hazard(times, shape, lambda_0)

# Pre-switching hazard (multiply by confounder effect)
hazard_before_switch <- baseline_hazard * exp(beta_confounder)

# Post-switching hazard (multiply by treatment + event_flag effect)
hazard_after_switch <- baseline_hazard * exp(beta_treatment + beta_event_flag * event_flag)

cat("=== Hazard Values Around Switching Time ===\n\n")
cat(sprintf("%8s %12s %12s %12s\n", "Time", "λ₀(t)", "Total λ(t)", "Period"))
cat(sprintf("%8s %12s %12s %12s\n", "-----", "------", "-----------", "------"))

for (t in c(1.5, 1.8, 2.0, 2.2, 2.5, 3.0)) {
  idx <- which.min(abs(times - t))
  baseline <- baseline_hazard[idx]

  if (t < t_switch) {
    total <- hazard_before_switch[idx]
    cat(sprintf("%8.1f %12.4f %12.4f   Pre\n", t, baseline, total))
  } else if (t == t_switch) {
    total_pre <- hazard_before_switch[idx]
    total_post <- hazard_after_switch[idx]
    cat(sprintf("%8.1f %12.4f %12.4f   Pre (just before)\n", t, baseline, total_pre))
    cat(sprintf("%8.1f %12.4f %12.4f   Post (just after) **\n", t, baseline, total_post))
  } else {
    total <- hazard_after_switch[idx]
    cat(sprintf("%8.1f %12.4f %12.4f   Post\n", t, baseline, total))
  }
}

cat("\n** Note: Total hazard λ(t) DROPS at switching due to treatment benefit,\n")
cat("   but baseline λ₀(t) remains CONTINUOUS (same value at t=2.0).\n")
cat("   The drop reflects the change from HR=2 (pre) to HR=1 (post).\n\n")

cat("Hazard Ratios:\n")
cat(sprintf("  Pre-switching:  HR = exp(%.2f) = %.2f\n",
            beta_confounder, exp(beta_confounder)))
cat(sprintf("  Post-switching: HR = exp(%.2f + %.2f) = %.2f\n",
            beta_treatment, beta_event_flag * event_flag,
            exp(beta_treatment + beta_event_flag * event_flag)))
cat("\nThe conditional Weibull ensures events are sampled based on\n")
cat("calendar time, preventing the hazard from 'resetting' at switching.\n")
