# ==============================================================================
# 12_stress_main_effects_summary.R
# Extract and visualize MAIN EFFECTS of stress on voice
# ==============================================================================
# PURPOSE:
#   Quantify how much stress ITSELF changes acoustic features
#   (before considering personality moderation)
#
# OUTPUT:
#   - Table of stress effect sizes for each acoustic feature
#   - Visualization showing magnitude and direction of stress effects
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(bayestestR)
  library(ggplot2)
  library(patchwork)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("STRESS MAIN EFFECTS: Effect Size Summary\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD MODELS
# ==============================================================================

cat("Loading fitted models...\n")

# Define model files
model_files <- list(
  list(
    file = "models/m_f0_mean_a.rds",
    name = "F0 Mean /a/",
    var = "f0_mean_a",
    unit = "Hz"
  ),
  list(
    file = "models/m_f0_mean_i.rds",
    name = "F0 Mean /i/",
    var = "f0_mean_i",
    unit = "Hz"
  ),
  list(
    file = "models/m_f0_mean_u.rds",
    name = "F0 Mean /u/",
    var = "f0_mean_u",
    unit = "Hz"
  ),

  list(
    file = "models/m_f0_std_a.rds",
    name = "F0 SD /a/",
    var = "f0_std_a",
    unit = "Hz"
  ),
  list(
    file = "models/m_f0_std_i.rds",
    name = "F0 SD /i/",
    var = "f0_std_i",
    unit = "Hz"
  ),
  list(
    file = "models/m_f0_std_u.rds",
    name = "F0 SD /u/",
    var = "f0_std_u",
    unit = "Hz"
  ),

  list(
    file = "models/m_f2_mean_a.rds",
    name = "F2 Mean /a/",
    var = "f2_mean_a",
    unit = "Hz"
  ),
  list(
    file = "models/m_f2_std_a.rds",
    name = "F2 SD /a/",
    var = "f2_std_a",
    unit = "Hz"
  ),

  list(
    file = "models/m_jitter_a.rds",
    name = "Jitter /a/",
    var = "jitter_a",
    unit = "%"
  ),
  list(
    file = "models/m_nne_a.rds",
    name = "NNE /a/",
    var = "nne_a",
    unit = "dB"
  )
)

# Load data for effect size standardization
df <- readRDS("results/df_analysis.rds")

# Extract stress effects
stress_effects <- map_df(model_files, function(model_info) {
  if (!file.exists(model_info$file)) {
    cat("⚠ Model not found:", model_info$file, "\n")
    return(NULL)
  }

  cat("Processing:", model_info$name, "\n")

  model <- readRDS(model_info$file)
  fe <- fixef(model)

  # Get outcome variable SD for standardization
  outcome_sd <- sd(df[[model_info$var]], na.rm = TRUE)

  # Extract c1_stress (baseline → pre-exam)
  if ("c1_stress" %in% rownames(fe)) {
    stress_est <- fe["c1_stress", "Estimate"]
    stress_lower <- fe["c1_stress", "Q2.5"]
    stress_upper <- fe["c1_stress", "Q97.5"]
    stress_se <- fe["c1_stress", "Est.Error"]

    # Standardized effect size (Cohen's d)
    stress_d <- stress_est / outcome_sd
  } else {
    stress_est <- stress_lower <- stress_upper <- stress_se <- stress_d <- NA
  }

  # Extract c2_recovery (pre-exam → post-exam)
  if ("c2_recovery" %in% rownames(fe)) {
    recovery_est <- fe["c2_recovery", "Estimate"]
    recovery_lower <- fe["c2_recovery", "Q2.5"]
    recovery_upper <- fe["c2_recovery", "Q97.5"]
    recovery_se <- fe["c2_recovery", "Est.Error"]

    # Standardized effect size
    recovery_d <- recovery_est / outcome_sd
  } else {
    recovery_est <- recovery_lower <- recovery_upper <- recovery_se <- recovery_d <- NA
  }

  tibble(
    outcome = model_info$name,
    outcome_var = model_info$var,
    unit = model_info$unit,

    # Stress effect (baseline → pre)
    stress_beta = stress_est,
    stress_ci_lower = stress_lower,
    stress_ci_upper = stress_upper,
    stress_se = stress_se,
    stress_d = stress_d,
    stress_credible = sign(stress_lower) == sign(stress_upper),

    # Recovery effect (pre → post)
    recovery_beta = recovery_est,
    recovery_ci_lower = recovery_lower,
    recovery_ci_upper = recovery_upper,
    recovery_se = recovery_se,
    recovery_d = recovery_d,
    recovery_credible = sign(recovery_lower) == sign(recovery_upper)
  )
})

cat("\n✓ Extraction complete\n\n")

# ==============================================================================
# 2. SUMMARY TABLE
# ==============================================================================

cat(rep("=", 70), "\n")
cat("STRESS MAIN EFFECTS: Summary Table\n")
cat(rep("=", 70), "\n\n")

# Format for display
summary_table <- stress_effects %>%
  mutate(
    # Format stress effect
    stress_effect = sprintf(
      "%.2f [%.2f, %.2f]%s (d = %.2f)",
      stress_beta,
      stress_ci_lower,
      stress_ci_upper,
      ifelse(stress_credible, "*", ""),
      stress_d
    ),

    # Format recovery effect
    recovery_effect = sprintf(
      "%.2f [%.2f, %.2f]%s (d = %.2f)",
      recovery_beta,
      recovery_ci_lower,
      recovery_ci_upper,
      ifelse(recovery_credible, "*", ""),
      recovery_d
    )
  ) %>%
  select(
    Outcome = outcome,
    Unit = unit,
    `Stress Effect (baseline → pre)` = stress_effect,
    `Recovery Effect (pre → post)` = recovery_effect
  )

print(summary_table, n = Inf)

cat("\n* = 95% CI excludes zero (credible effect)\n")
cat("d = Cohen's d (standardized effect size)\n\n")

# Interpretation guide
cat("Effect size interpretation (Cohen's d):\n")
cat("  |d| < 0.2  : Negligible\n")
cat("  |d| = 0.2  : Small\n")
cat("  |d| = 0.5  : Medium\n")
cat("  |d| = 0.8  : Large\n\n")

# ==============================================================================
# 3. COUNT CREDIBLE EFFECTS
# ==============================================================================

cat(rep("=", 70), "\n")
cat("CREDIBLE EFFECTS COUNT\n")
cat(rep("=", 70), "\n\n")

n_stress_credible <- sum(stress_effects$stress_credible, na.rm = TRUE)
n_recovery_credible <- sum(stress_effects$recovery_credible, na.rm = TRUE)

cat("Stress effects (baseline → pre-exam):", n_stress_credible, "credible\n")
cat("Recovery effects (pre → post-exam):", n_recovery_credible, "credible\n\n")

if (n_stress_credible > 0) {
  cat("Credible stress effects:\n")
  credible_stress <- stress_effects %>%
    filter(stress_credible) %>%
    select(outcome, stress_beta, stress_d) %>%
    arrange(desc(abs(stress_d)))
  print(credible_stress, n = Inf)
  cat("\n")
}

if (n_recovery_credible > 0) {
  cat("Credible recovery effects:\n")
  credible_recovery <- stress_effects %>%
    filter(recovery_credible) %>%
    select(outcome, recovery_beta, recovery_d) %>%
    arrange(desc(abs(recovery_d)))
  print(credible_recovery, n = Inf)
  cat("\n")
}

# ==============================================================================
# 4. VISUALIZATION
# ==============================================================================

cat("Creating visualization...\n")

# Prepare data for plotting
plot_data <- stress_effects %>%
  select(
    outcome,
    stress_d,
    stress_ci_lower,
    stress_ci_upper,
    recovery_d,
    recovery_ci_lower,
    recovery_ci_upper
  ) %>%
  pivot_longer(
    cols = -outcome,
    names_to = c("contrast", ".value"),
    names_pattern = "(.+)_(d|ci_lower|ci_upper)"
  ) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c("stress", "recovery"),
      labels = c("Stress\n(baseline → pre)", "Recovery\n(pre → post)")
    ),
    # Standardize CI names for plotting
    ci_lower_d = ci_lower,
    ci_upper_d = ci_upper
  )

# Forest plot
p_forest <- ggplot(plot_data, aes(x = d, y = outcome, color = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray70") +
  geom_linerange(
    aes(xmin = ci_lower_d, xmax = ci_upper_d),
    position = position_dodge(width = 0.5),
    linewidth = 1
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "Stress\n(baseline → pre)" = "#D32F2F",
      "Recovery\n(pre → post)" = "#1976D2"
    ),
    name = "Contrast"
  ) +
  labs(
    title = "Stress Main Effects: Standardized Effect Sizes",
    subtitle = "How much do acoustic features change under academic stress?",
    x = "Cohen's d (standardized effect size)",
    y = "Acoustic Feature"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/stress_main_effects_forest.png",
  p_forest,
  width = 10,
  height = 8,
  dpi = 300
)

cat("✓ Saved: figures/stress_main_effects_forest.png\n\n")

# Effect size magnitude plot
plot_data_magnitude <- plot_data %>%
  mutate(
    effect_size = abs(d),
    magnitude = case_when(
      effect_size < 0.2 ~ "Negligible",
      effect_size < 0.5 ~ "Small",
      effect_size < 0.8 ~ "Medium",
      TRUE ~ "Large"
    ),
    magnitude = factor(
      magnitude,
      levels = c("Negligible", "Small", "Medium", "Large")
    )
  )

p_magnitude <- ggplot(
  plot_data_magnitude,
  aes(x = contrast, y = effect_size, fill = magnitude)
) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_viridis_d(option = "plasma", direction = -1, name = "Magnitude") +
  labs(
    title = "Effect Size Magnitude Distribution",
    x = NULL,
    y = "Effect Size (|Cohen's d|)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  facet_wrap(~outcome, ncol = 2)

ggsave(
  "figures/stress_effect_size_magnitude.png",
  p_magnitude,
  width = 12,
  height = 10,
  dpi = 300
)

cat("✓ Saved: figures/stress_effect_size_magnitude.png\n\n")

# ==============================================================================
# 5. SAVE RESULTS
# ==============================================================================

cat("Saving results...\n")

dir.create("results", showWarnings = FALSE)

saveRDS(stress_effects, "results/stress_main_effects.rds")
write_csv(stress_effects, "results/stress_main_effects.csv")
write_csv(summary_table, "results/stress_main_effects_formatted.csv")

cat("✓ Saved: results/stress_main_effects.rds\n")
cat("✓ Saved: results/stress_main_effects.csv\n")
cat("✓ Saved: results/stress_main_effects_formatted.csv\n\n")

# ==============================================================================
# 6. MANUSCRIPT TEXT GENERATOR
# ==============================================================================

cat(rep("=", 70), "\n")
cat("SUGGESTED MANUSCRIPT TEXT\n")
cat(rep("=", 70), "\n\n")

cat("=== For Results Section ===\n\n")

cat(
  "Academic stress showed credible main effects on",
  n_stress_credible,
  "acoustic parameters.\n"
)

if (n_stress_credible > 0) {
  # Get strongest effects
  strongest <- stress_effects %>%
    filter(stress_credible) %>%
    arrange(desc(abs(stress_d))) %>%
    slice(1:min(3, n()))

  for (i in 1:nrow(strongest)) {
    cat(sprintf(
      "%s %s from baseline to pre-exam (β = %.2f %s, 95%% CI [%.2f, %.2f], d = %.2f).\n",
      strongest$outcome[i],
      ifelse(strongest$stress_beta[i] > 0, "increased", "decreased"),
      abs(strongest$stress_beta[i]),
      strongest$unit[i],
      strongest$stress_ci_lower[i],
      strongest$stress_ci_upper[i],
      strongest$stress_d[i]
    ))
  }
}

if (n_recovery_credible > 0) {
  cat(
    "\nRecovery effects were observed for",
    n_recovery_credible,
    "parameters.\n"
  )

  strongest_rec <- stress_effects %>%
    filter(recovery_credible) %>%
    arrange(desc(abs(recovery_d))) %>%
    slice(1:min(2, n()))

  for (i in 1:nrow(strongest_rec)) {
    cat(sprintf(
      "%s showed %s from pre- to post-exam (β = %.2f %s, d = %.2f).\n",
      strongest_rec$outcome[i],
      ifelse(strongest_rec$recovery_beta[i] > 0, "recovery", "further change"),
      abs(strongest_rec$recovery_beta[i]),
      strongest_rec$unit[i],
      strongest_rec$recovery_d[i]
    ))
  }
}

cat("\n=== For Methods Section ===\n\n")
cat("We examined stress main effects using orthogonal contrasts:\n")
cat("(1) c1_stress: baseline vs pre-exam (stress induction)\n")
cat("(2) c2_recovery: pre-exam vs post-exam (stress recovery)\n")
cat("Effect sizes were computed as Cohen's d using outcome variable SD.\n")

cat("\n", rep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n")
