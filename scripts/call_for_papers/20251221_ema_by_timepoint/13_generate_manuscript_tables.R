# ==============================================================================
# 13_generate_manuscript_tables.R
# Generate publication-ready tables for manuscript
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(knitr)
  library(kableExtra)
  library(gt)
  library(flextable)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("GENERATING MANUSCRIPT TABLES\n")
cat(rep("=", 70), "\n\n")

dir.create("tables", showWarnings = FALSE)

# ==============================================================================
# TABLE 1: Sample Characteristics
# ==============================================================================

cat("Creating Table 1: Sample Characteristics...\n")

df <- readRDS("results/df_analysis.rds")

# Person-level data
person_data <- df %>%
  distinct(ID, .keep_all = TRUE)

table1 <- tibble(
  Characteristic = c(
    "N participants",
    "N observations (voice recordings)",
    "Age, M (SD)",
    "Sex, N (%) female",
    "",
    "**PID-5 Dimensions (EMA), M (SD)**",
    "  Negative Affectivity",
    "  Detachment",
    "  Antagonism",
    "  Disinhibition",
    "  Psychoticism",
    "",
    "**PID-5 Dimensions (Baseline), M (SD)**",
    "  Negative Affectivity",
    "  Detachment",
    "  Antagonism",
    "  Disinhibition",
    "  Psychoticism",
    "",
    "**EMA Compliance**",
    "  Total EMA assessments, M (SD)",
    "  Completion rate, M (SD) %"
  ),

  Value = c(
    sprintf("%d", n_distinct(df$ID)),
    sprintf("%d", nrow(df)),
    sprintf(
      "%.1f (%.1f)",
      mean(person_data$age, na.rm = TRUE),
      sd(person_data$age, na.rm = TRUE)
    ),
    sprintf(
      "%d (%.1f%%)",
      sum(person_data$sex == "Femmina", na.rm = TRUE),
      100 * mean(person_data$sex == "Femmina", na.rm = TRUE)
    ),
    "",
    "",
    sprintf(
      "%.2f (%.2f)",
      mean(person_data$pid5_negative_affectivity, na.rm = TRUE),
      sd(person_data$pid5_negative_affectivity, na.rm = TRUE)
    ),
    sprintf(
      "%.2f (%.2f)",
      mean(person_data$pid5_detachment, na.rm = TRUE),
      sd(person_data$pid5_detachment, na.rm = TRUE)
    ),
    sprintf(
      "%.2f (%.2f)",
      mean(person_data$pid5_antagonism, na.rm = TRUE),
      sd(person_data$pid5_antagonism, na.rm = TRUE)
    ),
    sprintf(
      "%.2f (%.2f)",
      mean(person_data$pid5_disinhibition, na.rm = TRUE),
      sd(person_data$pid5_disinhibition, na.rm = TRUE)
    ),
    sprintf(
      "%.2f (%.2f)",
      mean(person_data$pid5_psychoticism, na.rm = TRUE),
      sd(person_data$pid5_psychoticism, na.rm = TRUE)
    ),
    "",
    "",
    sprintf(
      "%.1f (%.1f)",
      mean(person_data$domain_negative_affect_baseline, na.rm = TRUE),
      sd(person_data$domain_negative_affect_baseline, na.rm = TRUE)
    ),
    sprintf(
      "%.1f (%.1f)",
      mean(person_data$domain_detachment_baseline, na.rm = TRUE),
      sd(person_data$domain_detachment_baseline, na.rm = TRUE)
    ),
    sprintf(
      "%.1f (%.1f)",
      mean(person_data$domain_antagonism_baseline, na.rm = TRUE),
      sd(person_data$domain_antagonism_baseline, na.rm = TRUE)
    ),
    sprintf(
      "%.1f (%.1f)",
      mean(person_data$domain_disinhibition_baseline, na.rm = TRUE),
      sd(person_data$domain_disinhibition_baseline, na.rm = TRUE)
    ),
    sprintf(
      "%.1f (%.1f)",
      mean(person_data$domain_psychoticism_baseline, na.rm = TRUE),
      sd(person_data$domain_psychoticism_baseline, na.rm = TRUE)
    ),
    "",
    "",
    sprintf(
      "%.1f (%.1f)",
      mean(person_data$n_ema_obs, na.rm = TRUE),
      sd(person_data$n_ema_obs, na.rm = TRUE)
    ),
    sprintf(
      "%.1f (%.1f)",
      100 * mean(person_data$n_ema_obs / 31, na.rm = TRUE),
      100 * sd(person_data$n_ema_obs / 31, na.rm = TRUE)
    )
  )
)

# Save as CSV
write_csv(table1, "tables/table1_sample_characteristics.csv")

# Create formatted Word table
ft1 <- flextable(table1) %>%
  set_caption("Table 1. Sample Characteristics") %>%
  width(j = 1, width = 3.5) %>%
  width(j = 2, width = 1.5) %>%
  align(j = 2, align = "center") %>%
  bold(i = ~ grepl("\\*\\*", Characteristic)) %>%
  italic(i = ~ grepl("^  ", Characteristic))

save_as_docx(ft1, path = "tables/table1_sample_characteristics.docx")

cat("✓ Table 1 saved\n")

# ==============================================================================
# TABLE 2: Stress Main Effects
# ==============================================================================

cat("Creating Table 2: Stress Main Effects...\n")

stress_effects <- read_csv(
  "results/stress_main_effects.csv",
  show_col_types = FALSE
)

table2 <- stress_effects %>%
  mutate(
    # Format stress effect
    stress_formatted = sprintf(
      "%.2f [%.2f, %.2f]",
      stress_beta,
      stress_ci_lower,
      stress_ci_upper
    ),
    stress_formatted = if_else(
      stress_credible,
      paste0(stress_formatted, "*"),
      stress_formatted
    ),

    # Format Cohen's d
    stress_d_formatted = sprintf("%.2f", stress_d),

    # Format recovery effect
    recovery_formatted = sprintf(
      "%.2f [%.2f, %.2f]",
      recovery_beta,
      recovery_ci_lower,
      recovery_ci_upper
    ),
    recovery_formatted = if_else(
      recovery_credible,
      paste0(recovery_formatted, "*"),
      recovery_formatted
    ),

    # Format Cohen's d
    recovery_d_formatted = sprintf("%.2f", recovery_d)
  ) %>%
  select(
    `Acoustic Feature` = outcome,
    Unit = unit,
    `Stress Effect β [95% CI]` = stress_formatted,
    `d` = stress_d_formatted,
    `Recovery Effect β [95% CI]` = recovery_formatted,
    `d ` = recovery_d_formatted
  )

write_csv(table2, "tables/table2_stress_main_effects.csv")

ft2 <- flextable(table2) %>%
  set_caption(
    "Table 2. Main Effects of Academic Stress on Acoustic Features"
  ) %>%
  add_footer_lines(
    "Note. * = 95% credible interval excludes zero. Stress = baseline → pre-exam. Recovery = pre-exam → post-exam. d = Cohen's d standardized effect size."
  ) %>%
  width(j = 1, width = 1.5) %>%
  width(j = 2:6, width = 1.0) %>%
  align(j = 2:6, align = "center") %>%
  bold(part = "header")

save_as_docx(ft2, path = "tables/table2_stress_main_effects.docx")

cat("✓ Table 2 saved\n")

# ==============================================================================
# TABLE 3: Personality Moderation Effects (Main Results)
# ==============================================================================

cat("Creating Table 3: Personality Moderation Effects...\n")

# Load results from script 03
if (file.exists("results/moderation_results_summary.csv")) {
  moderation_results <- read_csv(
    "results/moderation_results_summary.csv",
    show_col_types = FALSE
  )
} else {
  # Extract from models manually
  cat("  Extracting from models...\n")

  # This is a placeholder - you'll need to adapt based on your actual results structure
  # I'm assuming you have the results from script 03_moderation_analysis_FINAL.R
  moderation_results <- tibble(
    parameter = character(),
    outcome = character(),
    vowel = character(),
    estimate = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    significant = logical()
  )
}

# Filter to significant interactions only
table3 <- moderation_results %>%
  filter(grepl(":", parameter)) %>% # Only interactions
  filter(significant == TRUE) %>%
  mutate(
    # Clean parameter names
    parameter_clean = parameter %>%
      str_remove("pid5_") %>%
      str_remove("_c") %>%
      str_replace("c1_stress", "Stress") %>%
      str_replace("c2_recovery", "Recovery") %>%
      str_replace("negative_affectivity", "Neg. Affect.") %>%
      str_replace("detachment", "Detach.") %>%
      str_replace("antagonism", "Antag.") %>%
      str_replace("disinhibition", "Disinhib.") %>%
      str_replace("psychoticism", "Psychot.") %>%
      str_replace(":", " × "),

    # Format outcome
    outcome_clean = paste0(
      str_replace(outcome, "_", " "),
      " /",
      vowel,
      "/"
    ),

    # Format estimate and CI
    estimate_formatted = sprintf(
      "%.2f [%.2f, %.2f]",
      estimate,
      ci_lower,
      ci_upper
    )
  ) %>%
  select(
    Interaction = parameter_clean,
    `Acoustic Feature` = outcome_clean,
    `β [95% CI]` = estimate_formatted
  ) %>%
  arrange(Interaction, `Acoustic Feature`)

write_csv(table3, "tables/table3_personality_moderation.csv")

ft3 <- flextable(table3) %>%
  set_caption(
    "Table 3. Personality × Stress Interactions on Acoustic Features"
  ) %>%
  add_footer_lines(
    "Note. Only credible interactions shown (95% CI excludes zero). β = unstandardized regression coefficient."
  ) %>%
  width(j = 1, width = 2.0) %>%
  width(j = 2, width = 1.5) %>%
  width(j = 3, width = 1.5) %>%
  align(j = 3, align = "center") %>%
  bold(part = "header")

save_as_docx(ft3, path = "tables/table3_personality_moderation.docx")

cat("✓ Table 3 saved\n")

# ==============================================================================
# TABLE 4: Model Comparison EMA vs Baseline
# ==============================================================================

cat("Creating Table 4: EMA vs Baseline Comparison...\n")

comparison <- readRDS("results/ema_vs_baseline_comparison.rds")

table4 <- comparison$r2_comparison %>%
  mutate(
    # Format R²
    R2_EMA_formatted = sprintf("%.3f", R2_EMA),
    R2_Baseline_formatted = sprintf("%.3f", R2_Baseline),
    R2_Diff_formatted = sprintf("%.3f", R2_Diff),

    # Add interpretation
    Interpretation = case_when(
      abs(R2_Diff) < 0.02 ~ "Equivalent",
      R2_Diff > 0.02 ~ "EMA superior",
      R2_Diff < -0.02 ~ "Baseline superior"
    )
  ) %>%
  select(
    `Outcome` = Outcome,
    `EMA R²` = R2_EMA_formatted,
    `Baseline R²` = R2_Baseline_formatted,
    `Difference` = R2_Diff_formatted,
    Interpretation
  )

write_csv(table4, "tables/table4_ema_vs_baseline.csv")

ft4 <- flextable(table4) %>%
  set_caption(
    "Table 4. Predictive Accuracy: EMA vs Baseline PID-5 Assessment"
  ) %>%
  add_footer_lines(
    "Note. EMA = 15 items, ~30 assessments. Baseline = 220 items, single assessment. Difference = EMA - Baseline."
  ) %>%
  width(j = 1:5, width = 1.0) %>%
  align(j = 2:5, align = "center") %>%
  bold(part = "header")

save_as_docx(ft4, path = "tables/table4_ema_vs_baseline.docx")

cat("✓ Table 4 saved\n")

# ==============================================================================
# TABLE S1: Within-Person Variability Effects (Supplementary)
# ==============================================================================

cat("Creating Table S1: Variability Effects...\n")

if (file.exists("results/variability/version_a_results.csv")) {
  variability_results <- read_csv(
    "results/variability/version_a_results.csv",
    show_col_types = FALSE
  )

  tableS1 <- variability_results %>%
    filter(sd_stress_credible | sd_recovery_credible) %>%
    mutate(
      # Format SD × Stress
      sd_stress_formatted = sprintf(
        "%.2f [%.2f, %.2f]%s",
        sd_stress_est,
        sd_stress_ci_lower,
        sd_stress_ci_upper,
        if_else(sd_stress_credible, "*", "")
      ),

      # Format SD × Recovery
      sd_recovery_formatted = sprintf(
        "%.2f [%.2f, %.2f]%s",
        sd_recovery_est,
        sd_recovery_ci_lower,
        sd_recovery_ci_upper,
        if_else(sd_recovery_credible, "*", "")
      )
    ) %>%
    select(
      Outcome = outcome,
      Domain = domain,
      `SD × Stress` = sd_stress_formatted,
      `SD × Recovery` = sd_recovery_formatted
    )

  write_csv(tableS1, "tables/tableS1_variability_effects.csv")

  ftS1 <- flextable(tableS1) %>%
    set_caption("Table S1. Within-Person Variability × Stress Interactions") %>%
    add_footer_lines(
      "Note. * = credible effect. SD = within-subject standard deviation in PID-5 across ~30 EMA assessments."
    ) %>%
    width(j = 1:4, width = 1.25) %>%
    align(j = 3:4, align = "center") %>%
    bold(part = "header")

  save_as_docx(ftS1, path = "tables/tableS1_variability_effects.docx")

  cat("✓ Table S1 saved\n")
} else {
  cat("⚠ Variability results not found, skipping Table S1\n")
}

# ==============================================================================
# TABLE S2: Full Model Parameters (Supplementary)
# ==============================================================================

cat("Creating Table S2: Full Model Parameters...\n")

# Example for one model - expand as needed
if (file.exists("models/m_f0_mean_a.rds")) {
  m <- readRDS("models/m_f0_mean_a.rds")
  fe <- fixef(m)

  tableS2 <- as_tibble(fe, rownames = "Parameter") %>%
    mutate(
      Parameter = str_remove(Parameter, "b_"),
      `β [95% CI]` = sprintf(
        "%.3f [%.3f, %.3f]",
        Estimate,
        Q2.5,
        Q97.5
      ),
      Credible = if_else(sign(Q2.5) == sign(Q97.5), "Yes", "No")
    ) %>%
    select(
      Parameter,
      `β [95% CI]`,
      Credible
    )

  write_csv(tableS2, "tables/tableS2_full_parameters_f0mean_a.csv")

  ftS2 <- flextable(tableS2) %>%
    set_caption("Table S2. Full Model Parameters: F0 Mean /a/ (Example)") %>%
    width(j = 1, width = 2.5) %>%
    width(j = 2, width = 1.5) %>%
    width(j = 3, width = 1.0) %>%
    align(j = 2:3, align = "center") %>%
    bold(part = "header")

  save_as_docx(ftS2, path = "tables/tableS2_full_parameters_example.docx")

  cat("✓ Table S2 saved\n")
}

cat("\n", rep("=", 70), "\n")
cat("TABLE GENERATION COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Tables saved in tables/ directory:\n")
cat("  - table1_sample_characteristics (CSV + DOCX)\n")
cat("  - table2_stress_main_effects (CSV + DOCX)\n")
cat("  - table3_personality_moderation (CSV + DOCX)\n")
cat("  - table4_ema_vs_baseline (CSV + DOCX)\n")
cat("  - tableS1_variability_effects (CSV + DOCX)\n")
cat("  - tableS2_full_parameters_example (CSV + DOCX)\n\n")
