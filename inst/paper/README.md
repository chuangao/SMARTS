# Reproducing the Paper Results

This directory contains scripts to reproduce the results from:

> [Author et al. (2024). Statistical Methods for Assessment of Real-world Treatment Switching (SMARTS). Journal. DOI: xxx]

## Requirements

- R >= 4.0
- SMARTS package (this package)

### R Package Dependencies

```r
install.packages(c("survival", "survminer", "ggplot2", "ggpubr",
                   "MatchIt", "cowplot", "parallel"))
```

## Installation

```r
# Install SMARTS from GitHub
devtools::install_github("chuangao/SMARTS")

# Or install from local source (if you cloned the repo)
devtools::install("path/to/SMARTS")
```

## Scripts Overview

| Script | Description | Output |
|--------|-------------|--------|
| `simulate_survival_confounder.R` | Simulate survival data with time-varying confounding | Data objects |
| `simulate_survival.R` | Basic survival simulation with treatment switching | Data objects |
| `compare_adjustment_methods.R` | Compare Naive, Covariate Adj, PSM, IPTW methods | Tables, summaries |
| `plot_km_curves.R` | Generate Kaplan-Meier curves (Figure X) | `output/km_curves_publication.pdf` |
| `verify_hazard_ratios.R` | Verify hazard ratio estimation accuracy | Console output |
| `verify_hazard_continuity.R` | Verify hazard continuity at switching time | Console output |
| `verify_hasbled_pattern.R` | Verify confounder pattern over time | Console output |

## Reproduce All Results

```r
# Set working directory to this folder
setwd(system.file("paper/scripts", package = "SMARTS"))

# Run scripts in order
source("simulate_survival_confounder.R")
source("compare_adjustment_methods.R")
source("plot_km_curves.R")
```

## Output

Generated figures and tables are saved to `output/`:

- `km_curves_publication.pdf` - Kaplan-Meier curves (Figure X in paper)
- `km_curves_publication.png` - PNG version for preview

## Reproducing Specific Figures

### Figure X: Kaplan-Meier Curves

```r
source("simulate_survival_confounder.R")
source("plot_km_curves.R")
# Output: output/km_curves_publication.pdf
```

### Table X: Method Comparison

```r
source("simulate_survival_confounder.R")
source("compare_adjustment_methods.R")
# Results printed to console
```

## Notes

- All scripts use `set.seed()` for reproducibility
- Simulation study with 1000 iterations takes approximately 30 minutes
- Use `n_cores` parameter in `run_simulation_study()` for parallel computing

## Contact

For questions about reproduction, please open an issue at:
https://github.com/chuangao/SMARTS/issues
