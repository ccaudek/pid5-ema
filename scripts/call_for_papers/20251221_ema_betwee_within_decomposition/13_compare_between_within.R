# ==============================================================================
# 13_compare_between_within.R
# Formal model comparison: Trait-only vs Trait+State
# ==============================================================================
# PURPOSE:
#   Determine whether adding within-person state components improves
#   model fit beyond trait-only models
#
# METHODS:
#   - Bayesian R² comparison
#   - LOO-IC (Leave-One-Out Information Criterion)
#   - WAIC (Widely Applicable Information Criterion)
#   - Bayes Factor (if feasible)
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
  trait_files <- list.files("models", pattern = "f0.*mean.*a", full.names = TRUE)
  if (length(trait_files) > 0) {
    fit_trait_only <- readRDS(trait_files[1])
    cat("✓ Loaded:", trait_files[1], "\n")
  } else {
    stop("Cannot find trait-only model. Run 03_moderation_analysis_FINAL.R first.")
  }
}

# Trait+State model (from between-within analysis)
if (file.exists("results/between_within/fit_f0_mean_a_bw.rds")) {
  fit_trait_state <- readRDS("results/between_within/fit_f0_mean_a_bw.rds")
  cat("✓ Trait+State model loaded\n\n")
} else {
  stop("Run 12_fit_between_within_models.R first")
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
cat(sprintf("  R² = %.3f [%.3f, %.3f]\n", 
            mean(r2_trait), 
            quantile(r2_trait, 0.025),
            quantile(r2_trait, 0.975)))

cat("\nTrait+State model:\n")
cat(sprintf("  R² = %.3f [%.3f, %.3f]\n",
            mean(r2_trait_state),
            quantile(r2_trait_state, 0.025),
            quantile(r2_trait_state, 0.975)))

# Difference
r2_diff <- r2_trait_state - r2_trait
cat("\nDifference (Trait+State - Trait-only):\n")
cat(sprintf("  ΔR² = %.3f [%.3f, %.3f]\n",
            mean(r2_diff),
            quantile(r2_diff, 0.025),
            quantile(r2_diff, 0.975)))

# Probability that Trait+State is better
prob_better <- mean(r2_diff > 0)
cat(sprintf("\nP(Trait+State better) = %.3f\n", prob_better))

if (prob_better > 0.95) {
  cat("→ STRONG evidence for Trait+State model\n\n")
} else if (prob_better > 0.90) {
  cat("→ MODERATE evidence for Trait+State model\n\n")
} else if (prob_better < 0.10) {
  cat("→ STRONG evidence for Trait-only model\n\n")
} else if (prob_better < 0.20) {
  cat("→ MODERATE evidence for Trait-only model\n\n")
} else {
  cat("→ INCONCLUSIVE: Models perform similarly\n\n")
}

# ==============================================================================
# 3. LOO-IC COMPARISON
# ==============================================================================

cat(rep("=", 70), "\n")
cat("LOO-IC COMPARISON (Lower is better)\n")
cat(rep("=", 70), "\n\n")

cat("Computing LOO for trait-only model...\n")
loo_trait <- loo(fit_trait_only, cores = 4)

cat("Computing LOO for trait+state model...\n")
loo_trait_state <- loo(fit_trait_state, cores = 4)

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

if (abs(elpd_diff) < se_diff) {
  cat("Models are statistically equivalent (|diff| < SE)\n")
  cat("→ Use simpler model (Trait-only)\n\n")
} else if (elpd_diff < -2*se_diff) {
  cat("Second model is clearly better (diff > 2*SE)\n")
  if (rownames(loo_comp)[2] == "fit_trait_state") {
    cat("→ Trait+State model preferred\n\n")
  } else {
    cat("→ Trait-only model preferred\n\n")
  }
} else {
  cat("Weak evidence for difference\n")
  cat("→ Consider practical significance\n\n")
}

# ==============================================================================
# 4. WAIC COMPARISON
# ==============================================================================

cat(rep("=", 70), "\n")
cat("WAIC COMPARISON (Lower is better)\n")
cat(rep("=", 70), "\n\n")

cat("Computing WAIC for trait-only model...\n")
waic_trait <- waic(fit_trait_only)

cat("Computing WAIC for trait+state model...\n")
waic_trait_state <- waic(fit_trait_state)

cat("\nTrait-only model:\n")
print(waic_trait)

cat("\nTrait+State model:\n")
print(waic_trait_state)

waic_comp <- loo_compare(waic_trait, waic_trait_state)

cat("\nComparison:\n")
print(waic_comp)

# ==============================================================================
# 5. PARAMETER EFFICIENCY
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PARAMETER EFFICIENCY\n")
cat(rep("=", 70), "\n\n")

# Count parameters
n_params_trait <- length(variables(fit_trait_only))
n_params_trait_state <- length(variables(fit_trait_state))

cat("Number of parameters:\n")
cat("  Trait-only:", n_params_trait, "\n")
cat("  Trait+State:", n_params_trait_state, "\n")
cat("  Additional parameters:", n_params_trait_state - n_params_trait, "\n\n")

# Extract effective parameters from LOO
p_loo_trait <- loo_trait$estimates["p_loo", "Estimate"]
p_loo_trait_state <- loo_trait_state$estimates["p_loo", "Estimate"]

cat("Effective number of parameters (p_loo):\n")
cat("  Trait-only:", round(p_loo_trait, 1), "\n")
cat("  Trait+State:", round(p_loo_trait_state, 1), "\n")
cat("  Difference:", round(p_loo_trait_state - p_loo_trait, 1), "\n\n")

cat("Interpretation:\n")
cat("If p_loo >> nominal parameters: Model may be overfitting\n")
cat("If Trait+State has much higher p_loo: Additional complexity not justified\n\n")

# ==============================================================================
# 6. SUMMARY AND RECOMMENDATION
# ==============================================================================

cat(rep("=", 70), "\n")
cat("SUMMARY AND RECOMMENDATION\n")
cat(rep("=", 70), "\n\n")

# Collect evidence
r2_favors_state <- mean(r2_diff) > 0.01  # At least 1% improvement
loo_favors_state <- elpd_diff < -2*se_diff && rownames(loo_comp)[1] == "fit_trait_state"

# Count credible within-person effects
if (file.exists("results/between_within/f0_mean_a_comparison.rds")) {
  comparison <- readRDS("results/between_within/f0_mean_a_comparison.rds")
  n_within_effects <- sum(
    comparison$within_main_credible,
    comparison$within_stress_credible,
    comparison$within_recovery_credible
  )
} else {
  n_within_effects <- 0
}

cat("Evidence summary:\n")
cat("  Bayesian R²:", ifelse(r2_favors_state, "✓ Favors Trait+State", "✗ No advantage"), "\n")
cat("  LOO-IC:", ifelse(loo_favors_state, "✓ Favors Trait+State", "✗ No clear winner"), "\n")
cat("  Within-person effects:", n_within_effects, "credible effects\n\n")

# Decision
if (loo_favors_state && n_within_effects >= 2) {
  cat("RECOMMENDATION: Include Trait+State in main manuscript\n\n")
  cat("Rationale:\n")
  cat("- Clear improvement in predictive accuracy (LOO)\n")
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

comparison_summary <- tibble(
  metric = c("R2_trait", "R2_trait_state", "R2_diff", "P_trait_state_better",
             "LOO_trait", "LOO_trait_state", "ELPD_diff", "SE_diff",
             "p_loo_trait", "p_loo_trait_state",
             "n_within_effects", "decision"),
  value = c(
    mean(r2_trait), mean(r2_trait_state), mean(r2_diff), prob_better,
    loo_trait$estimates["looic", "Estimate"],
    loo_trait_state$estimates["looic", "Estimate"],
    elpd_diff, se_diff,
    p_loo_trait, p_loo_trait_state,
    n_within_effects, NA
  ),
  ci_lower = c(
    quantile(r2_trait, 0.025), quantile(r2_trait_state, 0.025),
    quantile(r2_diff, 0.025), NA,
    rep(NA, 8)
  ),
  ci_upper = c(
    quantile(r2_trait, 0.975), quantile(r2_trait_state, 0.975),
    quantile(r2_diff, 0.975), NA,
    rep(NA, 8)
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
    decision
  )
)

saveRDS(comparison_summary, "results/between_within/model_comparison_summary.rds")
rio::export(comparison_summary, "results/between_within/model_comparison_summary.csv")

# Detailed report
sink("results/between_within/model_comparison_report.txt")
cat("MODEL COMPARISON REPORT: Trait-Only vs Trait+State\n")
cat(rep("=", 70), "\n\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

cat("BAYESIAN R²\n")
cat(rep("-", 70), "\n")
cat(sprintf("Trait-only: %.3f [%.3f, %.3f]\n", mean(r2_trait), 
            quantile(r2_trait, 0.025), quantile(r2_trait, 0.975)))
cat(sprintf("Trait+State: %.3f [%.3f, %.3f]\n", mean(r2_trait_state),
            quantile(r2_trait_state, 0.025), quantile(r2_trait_state, 0.975)))
cat(sprintf("Difference: %.3f [%.3f, %.3f]\n", mean(r2_diff),
            quantile(r2_diff, 0.025), quantile(r2_diff, 0.975)))
cat(sprintf("P(Trait+State better): %.3f\n\n", prob_better))

cat("LOO-IC COMPARISON\n")
cat(rep("-", 70), "\n")
print(loo_comp)
cat("\n")

cat("WAIC COMPARISON\n")
cat(rep("-", 70), "\n")
print(waic_comp)
cat("\n")

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
