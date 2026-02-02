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

## Paper Reproduction

Scripts for generating results presented in the SMARTS method paper are located in:

```
inst/paper/scripts/
```

### Key Scripts

| Script | Description |
|--------|-------------|
| `simulate_survival_confounder.R` | Simulates survival data with time-varying confounding |
| `validate_smarts_method.R` | Validates SMARTS method with helper functions |
| `run_factorial_simulation_parallel.R` | Main factorial simulation (hazard × confounding design) |
| `plot_factorial_results.R` | Generates boxplot visualization of simulation results |
| `plot_km_curves.R` | Generates Kaplan-Meier curves for the paper |

### Piecewise Constant Hazard Simulation

The `inst/paper/scripts/piecewise_constant/` folder contains simulations with **time-varying baseline hazard**. This tests SMARTS under realistic conditions where event rates change over calendar time.

| Script | Description |
|--------|-------------|
| `simulate_piecewise_hazard.R` | Simulates recurring events with piecewise constant hazard |
| `run_piecewise_scenarios.R` | Runs 4 clinical scenarios × 3 hazard trends |
| `plot_km_curves.R` | Generates KM curves comparing methods |
| `plot_results.R` | Generates boxplot comparing Baseline vs SMARTS |

**Key findings**: When baseline hazard changes over time (increasing or decreasing), the Baseline method is biased because it compares events from different time periods. SMARTS aligns comparison periods and substantially reduces this bias.

### Output Files

All simulation results and plots are saved to `inst/paper/output/`.

## Example: Simulating Data and Assigning Pseudo-Switching Times

This example demonstrates the complete workflow: simulating survival data with confounding, and using SMARTS to assign pseudo-switching times to continuers.

```r
library(SMARTS)
library(survival)
library(dplyr)

# Source the simulation function (from inst/paper/scripts/)
source("inst/paper/scripts/simulate_survival_confounder.R")
source("inst/paper/scripts/validate_smarts_method.R")

# =============================================================================
# Step 1: Simulate survival data with time-varying confounding
# =============================================================================

set.seed(123)

data <- simulate_survival_data_confounder(
  n_pairs = 1000,                    # 1000 switcher-continuer pairs
  beta_treatment = log(0.5),         # True treatment HR = 0.5
  beta_confounder = log(1.3),        # Confounder HR = 1.3
  lambda_0 = 10,                     # Baseline hazard scale
  shape = 1.5,                       # Weibull shape (>1 = accelerating hazard)
  T_min = 2,                         # Min follow-up time
  T_max = 6,                         # Max follow-up time
  switch_start = 0.25,               # Earliest switch at 25% of follow-up
  switch_end = 0.75,                 # Latest switch at 75% of follow-up
  confounder_interval = 0.5,         # Confounder measured every 0.5 years
  confounder_baseline_mean = 2.5,    # Mean confounder at baseline
  confounder_gap_baseline = 0.5,     # Switcher-continuer gap at baseline
  confounder_gap_peak = 1.6,         # Gap peaks at switch time
  confounder_gap_end = 0.8,          # Gap after treatment effect
  confounder_sd = 0.8                # Confounder standard deviation
)

# View the data structure
head(data)

# Split into switchers and continuers
switchers <- data[data$cohort == "switcher", ]
continuers <- data[data$cohort == "continuer", ]

cat("Switchers:", nrow(switchers), "\n")
cat("Continuers:", nrow(continuers), "\n")

# =============================================================================
# Step 2: Prepare data for SMARTS random_assign()
# =============================================================================

# SMARTS requires specific column names
switchers_smarts <- switchers
switchers_smarts$swi_yrs <- switchers_smarts$switch_time
switchers_smarts$fup_yrs <- switchers_smarts$T_max

continuers_smarts <- continuers
continuers_smarts$swi_yrs <- NA  # Will be assigned by SMARTS
continuers_smarts$fup_yrs <- continuers_smarts$T_max

# Create input list
smarts_input <- list(
  cont = continuers_smarts,
  swi = switchers_smarts
)

# =============================================================================
# Step 3: Assign pseudo-switching times using SMARTS
# =============================================================================

smarts_result <- random_assign(
  smarts_input,
  nbin = 10,              # Number of bins for matching
  seed = 123,             # Random seed for reproducibility
  swi_time = "swi_yrs",   # Column name for switch time
  cens_time = "fup_yrs"   # Column name for censoring time
)

# Check assignment results
cat("Assigned continuers:", nrow(smarts_result$assigned$cont), "\n")
cat("Unassigned continuers:", nrow(smarts_result$unassigned$cont), "\n")

# =============================================================================
# Step 4: Verify pseudo-switching time distribution matches switchers
# =============================================================================

# QQ plot to compare distributions
swi_times <- smarts_result$assigned$swi$swi_yrs
cont_times <- smarts_result$assigned$cont$swi_yrs

qqplot(swi_times, cont_times,
       main = "QQ Plot: Switching Times",
       xlab = "Switchers (actual)",
       ylab = "Continuers (SMARTS-assigned)")
abline(0, 1, col = "red", lty = 2)

# =============================================================================
# Step 5: Re-derive events based on pseudo-switching times
# =============================================================================

# Get assigned data
cont_assigned <- smarts_result$assigned$cont
swi_assigned <- smarts_result$assigned$swi

# Set new switch time
cont_assigned$new_switch_time <- cont_assigned$swi_yrs
swi_assigned$new_switch_time <- swi_assigned$switch_time

# Re-derive events (function from validate_smarts_method.R)
cont_assigned <- rederive_events(cont_assigned, "new_switch_time")
swi_assigned <- rederive_events(swi_assigned, "new_switch_time")

# Look up confounder at new switch time
cont_assigned <- lookup_confounder_at_time(cont_assigned, "new_switch_time", 0.5)
swi_assigned <- lookup_confounder_at_time(swi_assigned, "new_switch_time", 0.5)

# =============================================================================
# Step 6: Combine and analyze
# =============================================================================

# Select columns for analysis
cols_keep <- c("id", "cohort", "new_post_event", "new_post_event_time",
               "new_confounder_at_switch")

analysis_data <- rbind(
  cont_assigned[, cols_keep],
  swi_assigned[, cols_keep]
)

# Create treatment indicator
analysis_data$treated <- as.numeric(analysis_data$cohort == "switcher")

# Remove invalid observations
analysis_data <- analysis_data[analysis_data$new_post_event_time > 0, ]

# =============================================================================
# Step 7: Fit Cox models
# =============================================================================

# Naive (unadjusted)
cox_naive <- coxph(
  Surv(new_post_event_time, new_post_event) ~ treated,
  data = analysis_data
)
cat("\nNaive HR:", exp(coef(cox_naive)["treated"]), "\n")

# Adjusted for confounder
cox_adj <- coxph(
  Surv(new_post_event_time, new_post_event) ~ treated + new_confounder_at_switch,
  data = analysis_data
)
cat("Adjusted HR:", exp(coef(cox_adj)["treated"]), "\n")

# IPTW
ps_model <- glm(treated ~ new_confounder_at_switch,
                data = analysis_data, family = binomial)
analysis_data$ps <- predict(ps_model, type = "response")
analysis_data$weight <- ifelse(analysis_data$treated == 1,
                                1 / analysis_data$ps,
                                1 / (1 - analysis_data$ps))

cox_iptw <- coxph(
  Surv(new_post_event_time, new_post_event) ~ treated,
  data = analysis_data,
  weights = weight
)
cat("IPTW HR:", exp(coef(cox_iptw)["treated"]), "\n")
cat("\nTrue HR: 0.5\n")
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

## Key Concepts

1. **Timing Bias**: When comparing switchers to continuers, the time at which patients switch introduces bias because switchers are compared from their switch time while continuers are compared from baseline.

2. **Pseudo-Switching Time**: SMARTS assigns a "pseudo-switching time" to continuers that matches the distribution of actual switching times in switchers, aligning the comparison.

3. **Time-Varying Confounding**: Patients who switch often have different characteristics at the time of switching compared to those who continue. SMARTS helps address this by:
   - Aligning time origins
   - Enabling proper confounder adjustment at the (pseudo-)switching time

## License

GPL-3
