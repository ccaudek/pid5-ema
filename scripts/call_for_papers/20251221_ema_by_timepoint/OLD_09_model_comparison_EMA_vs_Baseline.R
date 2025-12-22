# ==============================================================================
# 09_model_comparison_EMA_vs_Baseline.R
# Do EMA measures add predictive value beyond baseline PID-5?
# ==============================================================================
# OBIETTIVO:
#   Confrontare modelli che usano:
#   - EMA aggregato (15 item, media su 2.5 mesi)
#   - PID-5 full baseline (220 item, single assessment)
#
# MOTIVAZIONE:
#   L'abstract afferma: "EMA measures added predictive value beyond baseline"
#   Questa analisi supporta empiricamente questo claim.
#
# APPROCCIO:
#   1. Fit modelli di moderazione con EMA (già fatto in script 03)
#   2. Fit modelli identici con PID-5 full baseline
#   3. Confronta R², LOO-IC, e effetti specifici
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayestestR)
  library(loo)
})

options(brms.backend = "cmdstanr")

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODEL COMPARISON: EMA vs Baseline PID-5\n")
cat("Do EMA measures add predictive value?\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

if (file.exists("results/df_analysis.rds")) {
  df <- readRDS("results/df_analysis.rds")
  cat("✓ Loaded analysis dataset\n")
} else {
  stop("Esegui prima 02_voice_personality_analysis_FINAL.R")
}

# Verify we have both EMA and baseline PID-5
cat("\nChecking available PID-5 measures:\n")

ema_vars <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

baseline_vars <- c(
  "domain_negative_affect_baseline",
  "domain_detachment_baseline",
  "domain_antagonism_baseline",
  "domain_disinhibition_baseline",
  "domain_psychoticism_baseline"
)

cat("EMA variables present:", all(ema_vars %in% names(df)), "\n")
cat("Baseline variables present:", all(baseline_vars %in% names(df)), "\n")

if (!all(baseline_vars %in% names(df))) {
  stop("Baseline PID-5 variables not found. Check data preparation.")
}

# ==============================================================================
# 2. PREPARE BASELINE PID-5 PREDICTORS
# ==============================================================================

cat("\n=== PREPARING BASELINE PID-5 PREDICTORS ===\n")

df <- df %>%
  mutate(
    # Center baseline PID-5 (for fair comparison with EMA)
    domain_negative_affect_baseline_c = scale(domain_negative_affect_baseline)[, 1],
    domain_detachment_baseline_c = scale(domain_detachment_baseline)[, 1],
    domain_antagonism_baseline_c = scale(domain_antagonism_baseline)[, 1],
    domain_disinhibition_baseline_c = scale(domain_disinhibition_baseline)[, 1],
    domain_psychoticism_baseline_c = scale(domain_psychoticism_baseline)[, 1]
  )

# Add contrasts if not present
if (!"c1_stress" %in% names(df)) {
  df <- df %>%
    mutate(
      c1_stress = case_when(
        timepoint == "baseline" ~ -0.5,
        timepoint == "pre" ~ 0.5,
        timepoint == "post" ~ 0
      ),
      c2_recovery = case_when(
        timepoint == "baseline" ~ 0,
        timepoint == "pre" ~ -0.5,
        timepoint == "post" ~ 0.5
      )
    )
}

cat("Baseline PID-5 centered and ready\n")

# ==============================================================================
# 3. DESCRIPTIVE COMPARISON
# ==============================================================================

cat("\n=== DESCRIPTIVE COMPARISON: EMA vs Baseline PID-5 ===\n\n")

# Get person-level data
person_data <- df %>%
  distinct(ID, .keep_all = TRUE)

# Compare means and SDs
comparison_desc <- tibble(
  Domain = c("Negative Affectivity", "Detachment", "Antagonism", 
             "Disinhibition", "Psychoticism")
) %>%
  mutate(
    EMA_M = c(
      mean(person_data$pid5_negative_affectivity, na.rm = TRUE),
      mean(person_data$pid5_detachment, na.rm = TRUE),
      mean(person_data$pid5_antagonism, na.rm = TRUE),
      mean(person_data$pid5_disinhibition, na.rm = TRUE),
      mean(person_data$pid5_psychoticism, na.rm = TRUE)
    ),
    EMA_SD = c(
      sd(person_data$pid5_negative_affectivity, na.rm = TRUE),
      sd(person_data$pid5_detachment, na.rm = TRUE),
      sd(person_data$pid5_antagonism, na.rm = TRUE),
      sd(person_data$pid5_disinhibition, na.rm = TRUE),
      sd(person_data$pid5_psychoticism, na.rm = TRUE)
    ),
    Baseline_M = c(
      mean(person_data$domain_negative_affect_baseline, na.rm = TRUE),
      mean(person_data$domain_detachment_baseline, na.rm = TRUE),
      mean(person_data$domain_antagonism_baseline, na.rm = TRUE),
      mean(person_data$domain_disinhibition_baseline, na.rm = TRUE),
      mean(person_data$domain_psychoticism_baseline, na.rm = TRUE)
    ),
    Baseline_SD = c(
      sd(person_data$domain_negative_affect_baseline, na.rm = TRUE),
      sd(person_data$domain_detachment_baseline, na.rm = TRUE),
      sd(person_data$domain_antagonism_baseline, na.rm = TRUE),
      sd(person_data$domain_disinhibition_baseline, na.rm = TRUE),
      sd(person_data$domain_psychoticism_baseline, na.rm = TRUE)
    )
  )

print(comparison_desc)

# ==============================================================================
# 4. FIT MODELS WITH BASELINE PID-5
# ==============================================================================

cat("\n=== FITTING MODELS WITH BASELINE PID-5 ===\n")
cat("Using same model structure as EMA models for direct comparison\n\n")

# Model formula with baseline predictors
baseline_traits <- "domain_negative_affect_baseline_c + domain_detachment_baseline_c + domain_antagonism_baseline_c + domain_disinhibition_baseline_c + domain_psychoticism_baseline_c"

build_formula_baseline <- function(outcome) {
  as.formula(
    paste0(
      outcome, " ~ c1_stress * (", baseline_traits, ") + ",
      "c2_recovery * (", baseline_traits, ") + ",
      "(1 + c1_stress + c2_recovery || ID)"
    )
  )
}

# Priors (same as EMA models)
prior_f0mean <- c(
  prior(student_t(3, 220, 30), class = Intercept),
  prior(normal(0, 10), class = b),
  prior(exponential(0.3), class = sigma),
  prior(exponential(0.3), class = sd)
)

prior_f2mean <- c(
  prior(student_t(3, 1500, 150), class = Intercept),
  prior(normal(0, 50), class = b),
  prior(exponential(0.01), class = sigma),
  prior(exponential(0.01), class = sd),
  prior(gamma(2, 0.1), class = nu)
)

# Fit baseline models for key outcomes
cat("Fitting baseline model: F0 Mean /a/...\n")
m_baseline_f0mean_a <- brm(
  formula = build_formula_baseline("f0_mean_a"),
  data = df,
  family = gaussian(),
  prior = prior_f0mean,
  iter = 5000, warmup = 2500, chains = 4, cores = 4, seed = 123,
  control = list(adapt_delta = 0.995, max_treedepth = 18),
  file = "models/m_baseline_f0mean_a"
)

cat("Fitting baseline model: F2 Mean /a/...\n")
m_baseline_f2mean_a <- brm(
  formula = build_formula_baseline("f2_mean_a"),
  data = df,
  family = student(),
  prior = prior_f2mean,
  iter = 5000, warmup = 2500, chains = 4, cores = 4, seed = 123,
  control = list(adapt_delta = 0.995, max_treedepth = 18),
  file = "models/m_baseline_f2mean_a"
)

# ==============================================================================
# 5. LOAD EMA MODELS FOR COMPARISON
# ==============================================================================

cat("\n=== LOADING EMA MODELS ===\n")

if (file.exists("models/m_f0_mean_a.rds")) {
  m_ema_f0mean_a <- readRDS("models/m_f0_mean_a.rds")
  cat("✓ Loaded EMA model: F0 Mean /a/\n")
} else {
  stop("EMA model not found. Run 03_moderation_analysis_FINAL.R first.")
}

if (file.exists("models/m_f2_mean_a.rds")) {
  m_ema_f2mean_a <- readRDS("models/m_f2_mean_a.rds")
  cat("✓ Loaded EMA model: F2 Mean /a/\n")
} else {
  stop("EMA model not found. Run 03_moderation_analysis_FINAL.R first.")
}

# ==============================================================================
# 6. MODEL COMPARISON: R²
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODEL COMPARISON: Bayesian R²\n")
cat(rep("=", 70), "\n\n")

# F0 Mean
cat("--- F0 Mean /a/ ---\n")
r2_ema_f0 <- bayes_R2(m_ema_f0mean_a)
r2_baseline_f0 <- bayes_R2(m_baseline_f0mean_a)

cat("EMA model R²:      ", sprintf("%.3f [%.3f, %.3f]", 
    mean(r2_ema_f0), quantile(r2_ema_f0, 0.025), quantile(r2_ema_f0, 0.975)), "\n")
cat("Baseline model R²: ", sprintf("%.3f [%.3f, %.3f]",
    mean(r2_baseline_f0), quantile(r2_baseline_f0, 0.025), quantile(r2_baseline_f0, 0.975)), "\n")
cat("Difference:        ", sprintf("%.3f", mean(r2_ema_f0) - mean(r2_baseline_f0)), "\n\n")

# F2 Mean
cat("--- F2 Mean /a/ ---\n")
r2_ema_f2 <- bayes_R2(m_ema_f2mean_a)
r2_baseline_f2 <- bayes_R2(m_baseline_f2mean_a)

cat("EMA model R²:      ", sprintf("%.3f [%.3f, %.3f]",
    mean(r2_ema_f2), quantile(r2_ema_f2, 0.025), quantile(r2_ema_f2, 0.975)), "\n")
cat("Baseline model R²: ", sprintf("%.3f [%.3f, %.3f]",
    mean(r2_baseline_f2), quantile(r2_baseline_f2, 0.025), quantile(r2_baseline_f2, 0.975)), "\n")
cat("Difference:        ", sprintf("%.3f", mean(r2_ema_f2) - mean(r2_baseline_f2)), "\n\n")

# ==============================================================================
# 7. MODEL COMPARISON: LOO-IC
# ==============================================================================

cat(rep("=", 70), "\n", sep = "")
cat("MODEL COMPARISON: LOO Information Criterion\n")
cat(rep("=", 70), "\n\n")

cat("Computing LOO for F0 Mean models...\n")
loo_ema_f0 <- loo(m_ema_f0mean_a)
loo_baseline_f0 <- loo(m_baseline_f0mean_a)
comp_f0 <- loo_compare(loo_ema_f0, loo_baseline_f0)

cat("\n--- F0 Mean /a/ LOO Comparison ---\n")
print(comp_f0)

cat("\nComputing LOO for F2 Mean models...\n")
loo_ema_f2 <- loo(m_ema_f2mean_a)
loo_baseline_f2 <- loo(m_baseline_f2mean_a)
comp_f2 <- loo_compare(loo_ema_f2, loo_baseline_f2)

cat("\n--- F2 Mean /a/ LOO Comparison ---\n")
print(comp_f2)

# ==============================================================================
# 8. EFFECT SIZE COMPARISON
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("EFFECT SIZE COMPARISON: Key Interactions\n")
cat(rep("=", 70), "\n\n")

# Extract key interaction effects
cat("--- Negative Affectivity × Stress → F0 Mean ---\n")

# EMA model
fe_ema <- fixef(m_ema_f0mean_a)
if ("c1_stress:pid5_negative_affectivity_c" %in% rownames(fe_ema)) {
  beta_ema <- fe_ema["c1_stress:pid5_negative_affectivity_c", ]
  cat("EMA:      β =", sprintf("%.3f", beta_ema["Estimate"]),
      ", 95% CI [", sprintf("%.3f", beta_ema["Q2.5"]), ",",
      sprintf("%.3f", beta_ema["Q97.5"]), "]\n")
}

# Baseline model
fe_baseline <- fixef(m_baseline_f0mean_a)
if ("c1_stress:domain_negative_affect_baseline_c" %in% rownames(fe_baseline)) {
  beta_baseline <- fe_baseline["c1_stress:domain_negative_affect_baseline_c", ]
  cat("Baseline: β =", sprintf("%.3f", beta_baseline["Estimate"]),
      ", 95% CI [", sprintf("%.3f", beta_baseline["Q2.5"]), ",",
      sprintf("%.3f", beta_baseline["Q97.5"]), "]\n")
}

# ==============================================================================
# 9. SUMMARY & CONCLUSIONS
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY: Do EMA measures add value?\n")
cat(rep("=", 70), "\n\n")

# Create summary table
comparison_summary <- tibble(
  Outcome = c("F0 Mean", "F2 Mean"),
  R2_EMA = c(mean(r2_ema_f0), mean(r2_ema_f2)),
  R2_Baseline = c(mean(r2_baseline_f0), mean(r2_baseline_f2)),
  R2_Diff = c(mean(r2_ema_f0) - mean(r2_baseline_f0),
              mean(r2_ema_f2) - mean(r2_baseline_f2)),
  LOO_Diff = c(
    comp_f0[2, "elpd_diff"],
    comp_f2[2, "elpd_diff"]
  )
)

print(comparison_summary)

cat("\nInterpretation:\n")

# R² interpretation
if (all(comparison_summary$R2_Diff > 0)) {
  cat("✓ EMA models show consistently higher R² than baseline models\n")
  cat("  → EMA measures capture additional variance in acoustic outcomes\n")
} else if (all(comparison_summary$R2_Diff < 0)) {
  cat("⚠ Baseline models show higher R² than EMA models\n")
  cat("  → Full PID-5 (220 items) outperforms brief EMA (15 items)\n")
} else {
  cat("• Mixed results: EMA advantage varies by outcome\n")
}

# LOO interpretation
if (all(comparison_summary$LOO_Diff < -2)) {
  cat("\n✓ STRONG evidence for EMA models (LOO difference > 2)\n")
  cat("  → EMA measures improve predictive accuracy\n")
} else if (all(comparison_summary$LOO_Diff > 2)) {
  cat("\n⚠ STRONG evidence for Baseline models (LOO difference > 2)\n")
  cat("  → Full PID-5 provides better predictions\n")
} else {
  cat("\n• Modest differences in predictive accuracy\n")
  cat("  → Both EMA and baseline provide similar information\n")
}

cat("\nCONCLUSION:\n")
if (mean(comparison_summary$R2_Diff) > 0.01 | mean(comparison_summary$LOO_Diff) < -2) {
  cat("EMA measures ADD PREDICTIVE VALUE beyond baseline assessment.\n")
  cat("The intensive longitudinal assessment captures dynamic trait expression.\n")
} else if (mean(comparison_summary$R2_Diff) < -0.01 | mean(comparison_summary$LOO_Diff) > 2) {
  cat("Baseline PID-5 (220 items) outperforms brief EMA (15 items).\n")
  cat("However, EMA still provides ecological validity and temporal resolution.\n")
} else {
  cat("EMA and baseline measures provide EQUIVALENT information.\n")
  cat("EMA offers the advantage of temporal resolution without loss of validity.\n")
}

# ==============================================================================
# 10. SAVE RESULTS
# ==============================================================================

comparison_results <- list(
  descriptive_comparison = comparison_desc,
  r2_comparison = comparison_summary,
  loo_f0 = comp_f0,
  loo_f2 = comp_f2
)

saveRDS(comparison_results, "results/ema_vs_baseline_comparison.rds")

cat("\n✓ Results saved: results/ema_vs_baseline_comparison.rds\n")

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODEL COMPARISON COMPLETE\n")
cat(rep("=", 70), "\n\n")
