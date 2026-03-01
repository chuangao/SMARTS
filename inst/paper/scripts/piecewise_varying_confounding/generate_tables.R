# ============================================================================
# Generate Summary Tables for Manuscript (CSV/Excel format)
# ============================================================================

library(dplyr)
library(tidyr)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Load raw results
raw <- readRDS("factorial_results_raw.rds")

# Set factor levels
raw$method <- factor(raw$method, levels = c("Naive", "Adj", "PSM", "IPTW", "SMR"))
raw$situation <- factor(raw$situation, levels = c("Baseline", "SMARTS"),
                       labels = c("Conventional", "SMARTS"))
raw$hazard_trend <- factor(raw$hazard_trend, levels = c("Constant", "Increasing", "Decreasing"))

# ============================================================================
# Compute summary statistics
# ============================================================================

summary_stats <- raw %>%
  group_by(clinical_scenario, true_hr, hazard_trend, situation, method) %>%
  summarise(
    mean_hr = mean(hr, na.rm = TRUE),
    sd_hr = sd(hr, na.rm = TRUE),
    median_hr = median(hr, na.rm = TRUE),
    bias = mean(hr, na.rm = TRUE) - first(true_hr),
    abs_bias = abs(mean(hr, na.rm = TRUE) - first(true_hr)),
    pct_bias = (mean(hr, na.rm = TRUE) - first(true_hr)) / first(true_hr) * 100,
    rmse = sqrt(mean((hr - first(true_hr))^2, na.rm = TRUE)),
    n_sims = sum(!is.na(hr)),
    .groups = "drop"
  )

# ============================================================================
# TABLE 1: Main results - Mean HR (SD) format (matches manuscript)
# One table per scenario, methods as columns
# ============================================================================

for (scenario in c("Sicker_Helps", "Healthier_Harms")) {
  scenario_data <- summary_stats %>%
    filter(clinical_scenario == scenario) %>%
    mutate(hr_sd = sprintf("%.2f (%.2f)", mean_hr, sd_hr)) %>%
    select(hazard_trend, situation, method, hr_sd) %>%
    pivot_wider(names_from = method, values_from = hr_sd)

  true_hr_val <- ifelse(scenario == "Sicker_Helps", 0.70, 1.50)
  scenario_label <- ifelse(scenario == "Sicker_Helps",
                           "Sicker Switchers - Treatment Helps",
                           "Healthier Switchers - Treatment Harms")

  filename <- paste0("results/tables/table_hr_sd_", scenario, ".csv")
  write.csv(scenario_data, filename, row.names = FALSE)
  cat("Saved:", filename, "(True HR =", true_hr_val, ")\n")
}

# ============================================================================
# TABLE 2: Bias table - methods as columns
# ============================================================================

for (scenario in c("Sicker_Helps", "Healthier_Harms")) {
  scenario_data <- summary_stats %>%
    filter(clinical_scenario == scenario) %>%
    mutate(bias_fmt = sprintf("%+.3f", bias)) %>%
    select(hazard_trend, situation, method, bias_fmt) %>%
    pivot_wider(names_from = method, values_from = bias_fmt)

  filename <- paste0("results/tables/table_bias_", scenario, ".csv")
  write.csv(scenario_data, filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
}

# ============================================================================
# TABLE 3: Combined table - Mean HR (SD) and Bias side by side
# One table per scenario
# ============================================================================

for (scenario in c("Sicker_Helps", "Healthier_Harms")) {
  scenario_data <- summary_stats %>%
    filter(clinical_scenario == scenario) %>%
    mutate(hr_sd_bias = sprintf("%.2f (%.2f) [%+.3f]", mean_hr, sd_hr, bias)) %>%
    select(hazard_trend, situation, method, hr_sd_bias) %>%
    pivot_wider(names_from = method, values_from = hr_sd_bias)

  filename <- paste0("results/tables/table_combined_", scenario, ".csv")
  write.csv(scenario_data, filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
}

# ============================================================================
# TABLE 4: Full detailed statistics (all numbers)
# ============================================================================

detailed <- summary_stats %>%
  mutate(
    mean_hr = round(mean_hr, 3),
    sd_hr = round(sd_hr, 3),
    median_hr = round(median_hr, 3),
    bias = round(bias, 3),
    abs_bias = round(abs_bias, 3),
    pct_bias = round(pct_bias, 1),
    rmse = round(rmse, 3)
  ) %>%
  arrange(clinical_scenario, hazard_trend, situation, method)

write.csv(detailed, "results/tables/table_detailed_all.csv", row.names = FALSE)
cat("Saved: results/tables/table_detailed_all.csv\n")

# ============================================================================
# TABLE 5: SMARTS improvement table
# Shows how much SMARTS reduces bias compared to Baseline
# ============================================================================

improvement <- summary_stats %>%
  select(clinical_scenario, true_hr, hazard_trend, situation, method, bias, abs_bias) %>%
  pivot_wider(
    names_from = situation,
    values_from = c(bias, abs_bias),
    names_sep = "_"
  ) %>%
  mutate(
    bias_reduction = abs_bias_Conventional - abs_bias_SMARTS,
    pct_improvement = round((abs_bias_Conventional - abs_bias_SMARTS) / abs_bias_Conventional * 100, 1)
  ) %>%
  mutate(across(where(is.numeric) & !matches("pct_improvement|true_hr"), ~ round(., 3))) %>%
  arrange(clinical_scenario, hazard_trend, method)

write.csv(improvement, "results/tables/table_smarts_improvement.csv", row.names = FALSE)
cat("Saved: results/tables/table_smarts_improvement.csv\n")

# ============================================================================
# TABLE 6: Compact manuscript table (both scenarios combined)
# Format: Hazard Trend | Approach | Naive | Adj | PSM | IPTW | SMR
# ============================================================================

compact_all <- summary_stats %>%
  mutate(hr_sd = sprintf("%.2f (%.2f)", mean_hr, sd_hr)) %>%
  select(clinical_scenario, true_hr, hazard_trend, situation, method, hr_sd) %>%
  pivot_wider(names_from = method, values_from = hr_sd) %>%
  arrange(clinical_scenario, hazard_trend, situation)

# Add scenario header rows
compact_final <- data.frame()
for (scenario in c("Sicker_Helps", "Healthier_Harms")) {
  true_hr_val <- ifelse(scenario == "Sicker_Helps", 0.70, 1.50)
  label <- ifelse(scenario == "Sicker_Helps",
                  "Scenario 1: Sicker Switchers, Treatment Helps (True HR = 0.70)",
                  "Scenario 2: Healthier Switchers, Treatment Harms (True HR = 1.50)")

  # Header row
  header <- data.frame(
    clinical_scenario = label,
    true_hr = "", hazard_trend = "", situation = "",
    Naive = "", Adj = "", PSM = "", IPTW = "", SMR = "",
    stringsAsFactors = FALSE
  )

  scenario_rows <- compact_all %>%
    filter(clinical_scenario == scenario) %>%
    mutate(across(everything(), as.character))

  compact_final <- bind_rows(compact_final, header, scenario_rows)
}

write.csv(compact_final, "results/tables/table_manuscript_compact.csv", row.names = FALSE)
cat("Saved: results/tables/table_manuscript_compact.csv\n")

# ============================================================================
# Print preview of main tables
# ============================================================================

cat("\n============================================================\n")
cat("PREVIEW: Scenario 1 - Sicker Switchers (True HR = 0.70)\n")
cat("============================================================\n")
preview1 <- summary_stats %>%
  filter(clinical_scenario == "Sicker_Helps") %>%
  mutate(hr_sd = sprintf("%.2f (%.2f)", mean_hr, sd_hr)) %>%
  select(hazard_trend, situation, method, hr_sd) %>%
  pivot_wider(names_from = method, values_from = hr_sd)
print(as.data.frame(preview1))

cat("\n============================================================\n")
cat("PREVIEW: Scenario 2 - Healthier Switchers (True HR = 1.50)\n")
cat("============================================================\n")
preview2 <- summary_stats %>%
  filter(clinical_scenario == "Healthier_Harms") %>%
  mutate(hr_sd = sprintf("%.2f (%.2f)", mean_hr, sd_hr)) %>%
  select(hazard_trend, situation, method, hr_sd) %>%
  pivot_wider(names_from = method, values_from = hr_sd)
print(as.data.frame(preview2))

cat("\n============================================================\n")
cat("PREVIEW: SMARTS Improvement (Bias Reduction)\n")
cat("============================================================\n")
improve_preview <- improvement %>%
  select(clinical_scenario, hazard_trend, method,
         bias_Conventional, bias_SMARTS, bias_reduction, pct_improvement)
print(as.data.frame(improve_preview))

cat("\nAll tables saved as CSV files. Open in Excel to copy into Word.\n")
