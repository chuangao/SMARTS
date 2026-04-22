# SMARTS (Statistical Methods for Assessment of Real-world Treatment Switching)

SMARTS is a tool for assigning pseudo-switching time to patients who initiate a treatment and stay on the same treatment (continuers), in comparison to patients who initiate a treatment and switch to another treatment (switchers). The assigned pseudo-switching time is used in analysis to reduce bias caused by time-varying confounding and timing bias.

## Installation

### Install using devtools

```r
library(devtools)
install_github("chuangao/SMARTS")
```

### Install from source

```bash
git clone https://github.com/chuangao/SMARTS
R CMD INSTALL SMARTS
```

## Quick Start with Mock Data

For a quick demonstration using the built-in mock data:

```r
library(SMARTS)
library(survival)

# Load mock data
head(mock_data)

# Convert to list format
data <- list(
  cont = mock_data[mock_data$swi == 0, ],
  swi = mock_data[mock_data$swi == 1, ]
)

# Assign pseudo-switching times
result <- random_assign(data, nbin = 10, seed = 123,
                        swi_time = "swi_yrs", cens_time = "fup_yrs")

# Check the distribution match
qqplot(result$assigned$swi$swi_yrs,
       result$assigned$cont$swi_yrs,
       main = "Switching Time Distribution",
       xlab = "Switchers", ylab = "Continuers (assigned)")
abline(0, 1, col = "red")
```

## Example: Assigning Pseudo-Switching Times with Simulated Data

This example demonstrates the full SMARTS workflow using simulated data with piecewise constant hazard and time-varying confounding. We simulate a dataset, use `random_assign()` to assign pseudo-switching times to continuers, verify the assignment quality, and run a propensity score matched (PSM) analysis. We then show how to repeat the assignment multiple times and aggregate results for more stable estimates.

```r
library(SMARTS)
library(survival)
library(MatchIt)
library(dplyr)

# Source the simulation and helper functions
source("inst/paper/scripts/piecewise_varying_confounding/simulate_piecewise_hazard.R")

# =============================================================================
# Step 1: Simulate survival data
# =============================================================================
# Simulates paired switcher-continuer data with:
#   - Piecewise constant baseline hazard (can be constant, increasing, or decreasing)
#   - Time-varying confounding: the confounder gap between switchers and continuers
#     evolves over time (gap_baseline -> gap_at_switch -> gap_end)
#   - Temporally correlated confounder trajectories (MVN with AR(1) structure)

set.seed(42)

data <- simulate_piecewise_hazard(
  n_pairs = 2000,
  beta_treatment = log(0.7),         # True treatment HR = 0.7
  beta_confounder = log(2.0),        # Confounder effect on hazard
  base_hazard = 0.02,
  hazard_trend = "increasing",       # Baseline hazard increases over time
  trend_strength = 0.25,
  segment_length = 0.5,
  confounder_gap_baseline = 0.5,     # Switchers slightly sicker at baseline
  confounder_gap_at_switch = 1.5,    # Gap peaks at switch time
  confounder_gap_end = 0.8,          # Gap shrinks after treatment helps
  rho = 0.9                          # Temporal correlation of confounder
)

switchers <- data[data$cohort == "switcher", ]
continuers <- data[data$cohort == "continuer", ]

cat("Switchers:", nrow(switchers), "\n")
cat("Continuers:", nrow(continuers), "\n")

# =============================================================================
# Step 2: Assign pseudo-switching times using SMARTS
# =============================================================================
# random_assign() requires:
#   - A list with $cont (continuers) and $swi (switchers)
#   - swi_time: column with actual switch times (switchers) or NA (continuers)
#   - cens_time: column with end-of-follow-up times

switchers$swi_yrs <- switchers$switch_time
switchers$fup_yrs <- switchers$T_max

continuers$swi_yrs <- NA
continuers$fup_yrs <- continuers$T_max

smarts_result <- random_assign(
  list(cont = continuers, swi = switchers),
  nbin = 10,
  seed = 123,
  swi_time = "swi_yrs",
  cens_time = "fup_yrs"
)

cat("Assigned continuers:", nrow(smarts_result$assigned$cont), "\n")
cat("Unassigned continuers:", nrow(smarts_result$unassigned$cont), "\n")

# =============================================================================
# Step 3: Verify the assignment with a QQ plot
# =============================================================================
# The assigned pseudo-switching times should closely match the actual
# switching time distribution of switchers.

qqplot(smarts_result$assigned$swi$swi_yrs,
       smarts_result$assigned$cont$swi_yrs,
       main = "QQ Plot: Switching Times",
       xlab = "Switchers (actual)",
       ylab = "Continuers (SMARTS-assigned)")
abline(0, 1, col = "red", lty = 2)

# =============================================================================
# Step 4: Re-derive events and build analysis dataset
# =============================================================================
# After assigning pseudo-switching times, re-derive which events fall
# before vs. after the (pseudo-)switch, and look up the confounder value
# at the assigned switch time.

cont_assigned <- smarts_result$assigned$cont
swi_assigned <- smarts_result$assigned$swi

cont_assigned$new_switch_time <- cont_assigned$swi_yrs
swi_assigned$new_switch_time <- swi_assigned$switch_time

cont_assigned <- rederive_events_recurring(cont_assigned, "new_switch_time")
swi_assigned <- rederive_events_recurring(swi_assigned, "new_switch_time")

# Build analysis data: post-switch survival time from switch
analysis_data <- data.frame(
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
  treated = c(rep(1, nrow(swi_assigned)), rep(0, nrow(cont_assigned))),
  confounder = c(swi_assigned$confounder_at_switch,
                 cont_assigned$confounder_at_switch)
)

analysis_data <- analysis_data[analysis_data$time > 0 & !is.na(analysis_data$time), ]

# =============================================================================
# Step 5: Propensity score matched analysis
# =============================================================================

ps_model <- glm(treated ~ confounder, data = analysis_data, family = binomial)
analysis_data$ps <- predict(ps_model, type = "response")

match_out <- matchit(treated ~ confounder, data = analysis_data,
                     method = "nearest", caliper = 0.2)
matched_data <- match.data(match_out)

cox_psm <- coxph(Surv(time, event) ~ treated, data = matched_data)
cat("PSM HR:", round(exp(coef(cox_psm)["treated"]), 3), "\n")
cat("True HR: 0.7\n")
```

### Repeated Assignment with Aggregation

Because `random_assign()` involves random sampling, running it multiple times with different seeds and averaging the results gives more stable HR estimates.

```r
n_reps <- 10
hr_results <- numeric(n_reps)

for (r in 1:n_reps) {
  # Re-assign with a different seed each time
  smarts_r <- random_assign(
    list(cont = continuers, swi = switchers),
    nbin = 10,
    seed = r * 100,
    swi_time = "swi_yrs",
    cens_time = "fup_yrs"
  )

  cont_r <- smarts_r$assigned$cont
  swi_r <- smarts_r$assigned$swi
  if (nrow(cont_r) == 0 || nrow(swi_r) == 0) next

  # Re-derive events
  cont_r$new_switch_time <- cont_r$swi_yrs
  swi_r$new_switch_time <- swi_r$switch_time
  cont_r <- rederive_events_recurring(cont_r, "new_switch_time")
  swi_r <- rederive_events_recurring(swi_r, "new_switch_time")

  # Build analysis data
  ad <- data.frame(
    time = c(
      ifelse(swi_r$new_has_post_event,
             swi_r$new_post_event_time_from_switch,
             swi_r$T_max - swi_r$new_switch_time),
      ifelse(cont_r$new_has_post_event,
             cont_r$new_post_event_time_from_switch,
             cont_r$T_max - cont_r$new_switch_time)
    ),
    event = c(as.numeric(swi_r$new_has_post_event),
              as.numeric(cont_r$new_has_post_event)),
    treated = c(rep(1, nrow(swi_r)), rep(0, nrow(cont_r))),
    confounder = c(swi_r$confounder_at_switch, cont_r$confounder_at_switch)
  )
  ad <- ad[ad$time > 0 & !is.na(ad$time), ]

  # PSM
  tryCatch({
    ps <- glm(treated ~ confounder, data = ad, family = binomial)
    ad$ps <- predict(ps, type = "response")
    m <- matchit(treated ~ confounder, data = ad, method = "nearest", caliper = 0.2)
    md <- match.data(m)
    cox <- coxph(Surv(time, event) ~ treated, data = md)
    hr_results[r] <- exp(coef(cox)["treated"])
  }, error = function(e) {
    hr_results[r] <<- NA
  })
}

cat("Individual HRs:", round(hr_results, 3), "\n")
cat("Mean HR:", round(mean(hr_results, na.rm = TRUE), 3), "\n")
cat("SD:", round(sd(hr_results, na.rm = TRUE), 3), "\n")
cat("True HR: 0.7\n")
```

## Paper Reproduction

Scripts for reproducing results from the SMARTS paper are located in `inst/paper/scripts/`.

### Piecewise Constant Hazard with Time-Varying Confounding (Primary)

The `inst/paper/scripts/piecewise_varying_confounding/` folder contains the main simulation study. It uses piecewise constant baseline hazard with time-varying confounding, where the confounder gap between switchers and continuers evolves over time (baseline → peak at switch → post-switch).

| Script | Description |
|--------|-------------|
| `simulate_piecewise_hazard.R` | Simulates recurring events with piecewise hazard and time-varying confounding |
| `run_factorial_simulation.R` | Factorial design: 2 clinical scenarios × 3 hazard trends, comparing Baseline vs SMARTS |
| `plot_factorial_results.R` | Boxplot visualization of simulation results |
| `plot_km_curves.R` | Kaplan-Meier curves for the paper |
| `generate_tables.R` | Generates summary tables |
| `explore_km_seeds.R` | Explores seeds for representative KM curves |

### Additional Simulation Scripts

| Folder | Description |
|--------|-------------|
| `piecewise_constant/` | Simulations with time-varying baseline hazard but constant confounding |
| Root scripts (`simulate_survival_confounder.R`, etc.) | Earlier Weibull-based simulations |

### Output Files

All simulation results and plots are saved to `inst/paper/output/`.

## Key Concepts

1. **Timing Bias**: When comparing switchers to continuers, the time at which patients switch introduces bias because switchers are compared from their switch time while continuers are compared from baseline.

2. **Pseudo-Switching Time**: SMARTS assigns a "pseudo-switching time" to continuers that matches the distribution of actual switching times in switchers, aligning the comparison.

3. **Time-Varying Confounding**: Patients who switch often have different confounder trajectories compared to continuers. The gap between switcher and continuer confounders typically widens leading up to the switch (e.g., switchers getting sicker), then narrows afterward. SMARTS helps address this by:
   - Aligning time origins via pseudo-switching time assignment
   - Enabling proper confounder adjustment at the (pseudo-)switching time
   - Supporting repeated random assignment for stable inference

## License

GPL-3
