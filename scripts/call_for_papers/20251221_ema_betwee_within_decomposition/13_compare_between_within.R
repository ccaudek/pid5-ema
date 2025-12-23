# ==============================================================================
# 13_compare_between_within.R
# Compare trait-only vs trait+state models
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(loo)
  library(bayesplot)
  library(patchwork)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODEL COMPARISON: Trait-Only vs Trait+State\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD MODELS
# ==============================================================================

cat("Loading models...\n\n")

# Trait-only model (from original analysis)
if (file.exists("models/f0_mean_a_moderation.rds")) {
  fit_trait_only <- readRDS("models/f0_mean_a_moderation.rds")
  cat("✓ Trait-only model loaded\n")
} else {
  cat("⚠ Trait-only model not found at models/f0_mean_a_moderation.rds\n")
  cat("  Attempting alternative location...\n")

  # Try alternative naming
  trait_files <- list.files(
    "models",
    pattern = "f0.*mean.*a",
    full.names = TRUE
  )
  if (length(trait_files) > 0) {
    fit_trait_only <- readRDS(trait_files[1])
    cat("✓ Loaded:", trait_files[1], "\n")
  } else {
    stop(
      "Cannot find trait-only model. Run 03_moderation_analysis_FINAL.R first."
    )
  }
}

# Trait+State model (from between-within analysis)
if (file.exists("results/between_within/fit_f0_mean_a_bw.rds")) {
  fit_trait_state <- readRDS("results/between_within/fit_f0_mean_a_bw.rds")
  cat("✓ Trait+State model loaded\n\n")
} else {
  stop("Run 12_fit_between_within_models.R first")
}

# Check sample sizes
n_obs_trait <- nobs(fit_trait_only)
n_obs_trait_state <- nobs(fit_trait_state)

cat("Sample sizes:\n")
cat("  Trait-only:", n_obs_trait, "observations\n")
cat("  Trait+State:", n_obs_trait_state, "observations\n")

if (n_obs_trait != n_obs_trait_state) {
  cat("\n⚠ WARNING: Models have different sample sizes!\n")
  cat("  This may be due to different missing data handling.\n")
  cat("  LOO comparison will not be possible.\n")
  cat("  We'll proceed with R² comparison only.\n\n")
  can_compare_loo <- FALSE
} else {
  can_compare_loo <- TRUE
  cat("  ✓ Sample sizes match\n\n")
}

# ==============================================================================
# 2. BAYESIAN R² COMPARISON
# ==============================================================================

cat(rep("=", 70), "\n")
cat("BAYESIAN R² COMPARISON\n")
cat(rep("=", 70), "\n\n")

# Compute Bayesian R²
r2_trait <- bayes_R2(fit_trait_only)
r2_trait_state <- bayes_R2(fit_trait_state)

cat("Trait-only model:\n")
cat(sprintf(
  "  R² = %.3f [%.3f, %.3f]\n",
  mean(r2_trait),
  quantile(r2_trait, 0.025),
  quantile(r2_trait, 0.975)
))

cat("\nTrait+State model:\n")
cat(sprintf(
  "  R² = %.3f [%.3f, %.3f]\n",
  mean(r2_trait_state),
  quantile(r2_trait_state, 0.025),
  quantile(r2_trait_state, 0.975)
))

# Difference
r2_diff <- r2_trait_state - r2_trait
cat("\nDifference (Trait+State - Trait-only):\n")
cat(sprintf(
  "  ΔR² = %.3f [%.3f, %.3f]\n",
  mean(r2_diff),
  quantile(r2_diff, 0.025),
  quantile(r2_diff, 0.975)
))

# Probability that Trait+State is better
prob_better <- mean(r2_diff > 0)
cat(sprintf("\nP(Trait+State better) = %.3f\n", prob_better))

if (prob_better > 0.95) {
  cat("→ STRONG evidence for Trait+State model\n\n")
  r2_conclusion <- "strong_trait_state"
} else if (prob_better > 0.90) {
  cat("→ MODERATE evidence for Trait+State model\n\n")
  r2_conclusion <- "moderate_trait_state"
} else if (prob_better < 0.10) {
  cat("→ STRONG evidence for Trait-only model\n\n")
  r2_conclusion <- "strong_trait_only"
} else if (prob_better < 0.20) {
  cat("→ MODERATE evidence for Trait-only model\n\n")
  r2_conclusion <- "moderate_trait_only"
} else {
  cat("→ INCONCLUSIVE: Models perform similarly\n\n")
  r2_conclusion <- "inconclusive"
}

# ==============================================================================
# 3. LOO-IC COMPARISON
# ==============================================================================

if (can_compare_loo) {
  cat(rep("=", 70), "\n")
  cat("LOO-IC COMPARISON (Lower is better)\n")
  cat(rep("=", 70), "\n\n")

  cat("Computing LOO for trait-only model...\n")
  loo_trait <- loo(fit_trait_only, cores = 4, moment_match = TRUE)

  cat("Computing LOO for trait+state model...\n")
  loo_trait_state <- loo(fit_trait_state, cores = 4, moment_match = TRUE)

  cat("\nTrait-only model:\n")
  print(loo_trait)

  cat("\nTrait+State model:\n")
  print(loo_trait_state)

  # Formal comparison
  cat("\n", rep("-", 70), "\n")
  cat("FORMAL COMPARISON\n")
  cat(rep("-", 70), "\n\n")

  loo_comp <- loo_compare(loo_trait, loo_trait_state)
  print(loo_comp)

  cat("\nInterpretation:\n")
  elpd_diff <- loo_comp[2, "elpd_diff"]
  se_diff <- loo_comp[2, "se_diff"]

  cat(sprintf("ELPD difference: %.1f (SE = %.1f)\n", elpd_diff, se_diff))

  if (abs(elpd_diff) < se_diff) {
    cat("Models are statistically equivalent (|diff| < SE)\n")
    cat("→ Use simpler model (Trait-only)\n\n")
    loo_conclusion <- "equivalent_prefer_simpler"
  } else if (elpd_diff < -2 * se_diff) {
    cat("Second model is clearly better (diff > 2*SE)\n")
    if (rownames(loo_comp)[2] == "fit_trait_state") {
      cat("→ Trait+State model preferred\n\n")
      loo_conclusion <- "trait_state_better"
    } else {
      cat("→ Trait-only model preferred\n\n")
      loo_conclusion <- "trait_only_better"
    }
  } else {
    cat("Weak evidence for difference\n")
    cat("→ Consider practical significance\n\n")
    loo_conclusion <- "weak_evidence"
  }

  # Extract p_loo
  p_loo_trait <- loo_trait$estimates["p_loo", "Estimate"]
  p_loo_trait_state <- loo_trait_state$estimates["p_loo", "Estimate"]

  loo_available <- TRUE
} else {
  cat(rep("=", 70), "\n")
  cat("LOO-IC COMPARISON: SKIPPED\n")
  cat(rep("=", 70), "\n\n")
  cat("Cannot compare models with different sample sizes.\n\n")

  loo_conclusion <- "not_comparable"
  elpd_diff <- NA
  se_diff <- NA
  p_loo_trait <- NA
  p_loo_trait_state <- NA
  loo_available <- FALSE
}

# ==============================================================================
# 4. PARAMETER EFFICIENCY
# ==============================================================================

cat(rep("=", 70), "\n")
cat("PARAMETER EFFICIENCY\n")
cat(rep("=", 70), "\n\n")

# Count parameters
n_params_trait <- length(variables(fit_trait_only))
n_params_trait_state <- length(variables(fit_trait_state))

cat("Number of parameters:\n")
cat("  Trait-only:", n_params_trait, "\n")
cat("  Trait+State:", n_params_trait_state, "\n")
cat("  Additional parameters:", n_params_trait_state - n_params_trait, "\n\n")

if (loo_available) {
  cat("Effective number of parameters (p_loo):\n")
  cat("  Trait-only:", round(p_loo_trait, 1), "\n")
  cat("  Trait+State:", round(p_loo_trait_state, 1), "\n")
  cat("  Difference:", round(p_loo_trait_state - p_loo_trait, 1), "\n\n")

  cat("Interpretation:\n")
  cat("If p_loo >> nominal parameters: Model may be overfitting\n")
  cat(
    "If Trait+State has much higher p_loo: Additional complexity not justified\n\n"
  )
}

# ==============================================================================
# 5. WITHIN-PERSON EFFECTS
# ==============================================================================

cat(rep("=", 70), "\n")
cat("WITHIN-PERSON EFFECTS\n")
cat(rep("=", 70), "\n\n")

# Count credible within-person effects
if (file.exists("results/between_within/f0_mean_a_comparison.rds")) {
  comparison <- readRDS("results/between_within/f0_mean_a_comparison.rds")

  n_within_main <- sum(comparison$within_main_credible, na.rm = TRUE)
  n_within_stress <- sum(comparison$within_stress_credible, na.rm = TRUE)
  n_within_recovery <- sum(comparison$within_recovery_credible, na.rm = TRUE)
  n_within_effects <- n_within_main + n_within_stress + n_within_recovery

  n_between_main <- sum(comparison$between_main_credible, na.rm = TRUE)
  n_between_stress <- sum(comparison$between_stress_credible, na.rm = TRUE)
  n_between_recovery <- sum(comparison$between_recovery_credible, na.rm = TRUE)
  n_between_effects <- n_between_main + n_between_stress + n_between_recovery

  cat("Credible effects (95% CI excludes 0):\n")
  cat("  Between-person:\n")
  cat("    Main effects:", n_between_main, "\n")
  cat("    Stress interactions:", n_between_stress, "\n")
  cat("    Recovery interactions:", n_between_recovery, "\n")
  cat("    TOTAL:", n_between_effects, "\n\n")

  cat("  Within-person:\n")
  cat("    Main effects:", n_within_main, "\n")
  cat("    Stress interactions:", n_within_stress, "\n")
  cat("    Recovery interactions:", n_within_recovery, "\n")
  cat("    TOTAL:", n_within_effects, "\n\n")
} else {
  n_within_effects <- 0
  n_between_effects <- 0
  cat(
    "⚠ Comparison file not found. Run 12_fit_between_within_models.R first.\n\n"
  )
}

# ==============================================================================
# 6. SUMMARY AND RECOMMENDATION
# ==============================================================================

cat(rep("=", 70), "\n")
cat("SUMMARY AND RECOMMENDATION\n")
cat(rep("=", 70), "\n\n")

# Collect evidence
r2_favors_state <- mean(r2_diff) > 0.01 # At least 1% improvement

if (loo_available) {
  loo_favors_state <- (loo_conclusion == "trait_state_better")
} else {
  loo_favors_state <- FALSE
}

cat("Evidence summary:\n")
cat(
  "  Bayesian R²:",
  ifelse(r2_favors_state, "✓ Favors Trait+State", "✗ No advantage"),
  "\n"
)

if (loo_available) {
  cat(
    "  LOO-IC:",
    ifelse(
      loo_favors_state,
      "✓ Favors Trait+State",
      "✗ Favors Trait-only or equivalent"
    ),
    "\n"
  )
} else {
  cat("  LOO-IC: (not comparable)\n")
}

cat("  Within-person effects:", n_within_effects, "credible effects\n")
cat("  Between-person effects:", n_between_effects, "credible effects\n\n")

# Decision
if (
  (loo_available && loo_favors_state && n_within_effects >= 2) ||
    (!loo_available && r2_favors_state && n_within_effects >= 2)
) {
  cat("RECOMMENDATION: Include Trait+State in main manuscript\n\n")
  cat("Rationale:\n")
  cat("- Clear improvement in model fit\n")
  cat("- Multiple credible within-person effects\n")
  cat("- Theoretical value: Demonstrates state-trait dissociation\n\n")
  decision <- "main_manuscript"
} else if (n_within_effects >= 1 && mean(r2_diff) > 0.005) {
  cat("RECOMMENDATION: Include in supplementary materials\n\n")
  cat("Rationale:\n")
  cat("- Some evidence for within-person effects\n")
  cat("- Model fit improvement is marginal\n")
  cat("- Interesting for specialists but not essential for main narrative\n\n")
  decision <- "supplementary"
} else {
  cat("RECOMMENDATION: Supplementary materials only, brief mention\n\n")
  cat("Rationale:\n")
  cat("- No substantial improvement over trait-only model\n")
  cat("- Trait-only model is simpler and sufficient\n")
  cat("- Between-person analysis (current manuscript) tells main story\n\n")
  decision <- "supplementary_brief"
}

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

cat("Saving comparison results...\n")

dir.create("results/between_within", showWarnings = FALSE, recursive = TRUE)

comparison_summary <- tibble(
  metric = c(
    "R2_trait",
    "R2_trait_state",
    "R2_diff",
    "P_trait_state_better",
    "LOO_trait",
    "LOO_trait_state",
    "ELPD_diff",
    "SE_diff",
    "p_loo_trait",
    "p_loo_trait_state",
    "n_within_effects",
    "n_between_effects",
    "decision"
  ),
  value = c(
    mean(r2_trait),
    mean(r2_trait_state),
    mean(r2_diff),
    prob_better,
    ifelse(loo_available, loo_trait$estimates["looic", "Estimate"], NA),
    ifelse(loo_available, loo_trait_state$estimates["looic", "Estimate"], NA),
    ifelse(loo_available, elpd_diff, NA),
    ifelse(loo_available, se_diff, NA),
    p_loo_trait,
    p_loo_trait_state,
    n_within_effects,
    n_between_effects,
    NA
  ),
  ci_lower = c(
    quantile(r2_trait, 0.025),
    quantile(r2_trait_state, 0.025),
    quantile(r2_diff, 0.025),
    NA,
    rep(NA, 9)
  ),
  ci_upper = c(
    quantile(r2_trait, 0.975),
    quantile(r2_trait_state, 0.975),
    quantile(r2_diff, 0.975),
    NA,
    rep(NA, 9)
  ),
  interpretation = c(
    "Trait-only variance explained",
    "Trait+State variance explained",
    "Improvement from adding state",
    "Probability Trait+State better",
    "Trait-only predictive accuracy",
    "Trait+State predictive accuracy",
    "Difference in expected log predictive density",
    "Standard error of difference",
    "Effective parameters trait-only",
    "Effective parameters trait+state",
    "Number of credible within-person effects",
    "Number of credible between-person effects",
    decision
  )
)

saveRDS(
  comparison_summary,
  "results/between_within/model_comparison_summary.rds"
)
rio::export(
  comparison_summary,
  "results/between_within/model_comparison_summary.csv"
)

# Detailed report
sink("results/between_within/model_comparison_report.txt")
cat("MODEL COMPARISON REPORT: Trait-Only vs Trait+State\n")
cat(rep("=", 70), "\n\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

cat("BAYESIAN R²\n")
cat(rep("-", 70), "\n")
cat(sprintf(
  "Trait-only: %.3f [%.3f, %.3f]\n",
  mean(r2_trait),
  quantile(r2_trait, 0.025),
  quantile(r2_trait, 0.975)
))
cat(sprintf(
  "Trait+State: %.3f [%.3f, %.3f]\n",
  mean(r2_trait_state),
  quantile(r2_trait_state, 0.025),
  quantile(r2_trait_state, 0.975)
))
cat(sprintf(
  "Difference: %.3f [%.3f, %.3f]\n",
  mean(r2_diff),
  quantile(r2_diff, 0.025),
  quantile(r2_diff, 0.975)
))
cat(sprintf("P(Trait+State better): %.3f\n", prob_better))
cat(sprintf("Conclusion: %s\n\n", r2_conclusion))

if (loo_available) {
  cat("LOO-IC COMPARISON\n")
  cat(rep("-", 70), "\n")
  print(loo_comp)
  cat(sprintf("\nConclusion: %s\n\n", loo_conclusion))
} else {
  cat("LOO-IC COMPARISON\n")
  cat(rep("-", 70), "\n")
  cat("Not available (different sample sizes)\n\n")
}

cat("WITHIN-PERSON EFFECTS\n")
cat(rep("-", 70), "\n")
cat(sprintf("Credible within-person effects: %d\n", n_within_effects))
cat(sprintf("Credible between-person effects: %d\n\n", n_between_effects))

cat("RECOMMENDATION\n")
cat(rep("-", 70), "\n")
cat(decision, "\n")
sink()

cat("✓ Summary saved: results/between_within/model_comparison_summary.csv\n")
cat("✓ Report saved: results/between_within/model_comparison_report.txt\n\n")

cat(rep("=", 70), "\n")
cat("COMPARISON COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Next steps:\n")
if (decision == "main_manuscript") {
  cat("1. Run 14_visualize_between_within.R for publication figures\n")
  cat("2. Extend to other outcomes (F0 /i/ /u/, F2, etc.)\n")
  cat("3. Update manuscript Methods and Results sections\n\n")
} else {
  cat("1. Run 14_visualize_between_within.R for supplementary figures\n")
  cat("2. Decide whether to extend to other outcomes\n")
  cat("3. Prepare brief supplementary note\n\n")
}
