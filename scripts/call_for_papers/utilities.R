# ==============================================================================
# UTILITIES: Tables, Diagnostics, and Sensitivity Analyses
# Supporting analyses for manuscript
# ==============================================================================

library(tidyverse)
library(brms)
library(bayesplot)
library(bayestestR)
library(knitr)
library(kableExtra)
library(ggdist)
library(patchwork)

# Load all analyses
load("results/moderation_analysis_complete.RData")

# ==============================================================================
# 1. TABLE 1: SAMPLE CHARACTERISTICS
# ==============================================================================

cat("\n=== Creating Table 1: Sample Characteristics ===\n")

# Calculate descriptive statistics for PID-5 at baseline
table1_pid5 <- df_analysis %>%
  filter(timepoint == "baseline") %>%
  select(ID, ends_with("_bl")) %>%
  select(-ends_with("_bl_c")) %>% # Remove centered versions
  distinct() %>%
  pivot_longer(cols = -ID, names_to = "trait", values_to = "score") %>%
  group_by(trait) %>%
  summarise(
    M = mean(score, na.rm = TRUE),
    SD = sd(score, na.rm = TRUE),
    Min = min(score, na.rm = TRUE),
    Max = max(score, na.rm = TRUE),
    N = sum(!is.na(score)),
    .groups = "drop"
  ) %>%
  mutate(
    trait = str_remove(trait, "pid5_") %>%
      str_remove("_bl") %>%
      str_replace("na", "Negative Affectivity") %>%
      str_replace("det", "Detachment") %>%
      str_replace("ant", "Antagonism") %>%
      str_replace("dis", "Disinhibition") %>%
      str_replace("psy", "Psychoticism")
  )

# Calculate descriptives for acoustic variables
table1_acoustics <- df_analysis %>%
  group_by(timepoint) %>%
  summarise(
    across(
      c(f0_mean, f0_std, jitter, nne),
      list(M = ~ mean(., na.rm = TRUE), SD = ~ sd(., na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    n = n(),
    .groups = "drop"
  )

# Print and save
cat("\nPID-5 Traits (Baseline):\n")
print(table1_pid5)

cat("\nAcoustic Features by Timepoint:\n")
print(table1_acoustics)

# Create formatted table for manuscript
table1_formatted <- table1_pid5 %>%
  mutate(
    `M (SD)` = sprintf("%.2f (%.2f)", M, SD),
    Range = sprintf("%.2f - %.2f", Min, Max)
  ) %>%
  select(Trait = trait, `M (SD)`, Range, N)

# Save as CSV
write_csv(table1_formatted, "results/table1_sample_characteristics.csv")
write_csv(table1_acoustics, "results/table1_acoustics_descriptives.csv")

cat("\n✓ Table 1 saved to results/table1_*.csv\n")

# ==============================================================================
# 2. TABLE 2: MAIN EFFECTS SUMMARY
# ==============================================================================

cat("\n=== Creating Table 2: Main Effects ===\n")

create_main_effects_table <- function(models, outcome_names) {
  results <- map2_dfr(models, outcome_names, function(model, outcome) {
    # Get fixed effects
    fixed <- fixef(model)

    # Extract timepoint effects
    tibble(
      Outcome = outcome,
      `PRE vs BASE β` = fixed["timepointpre", "Estimate"],
      `PRE vs BASE 95% CI` = sprintf(
        "[%.2f, %.2f]",
        fixed["timepointpre", "Q2.5"],
        fixed["timepointpre", "Q97.5"]
      ),
      `POST vs BASE β` = fixed["timepointpost", "Estimate"],
      `POST vs BASE 95% CI` = sprintf(
        "[%.2f, %.2f]",
        fixed["timepointpost", "Q2.5"],
        fixed["timepointpost", "Q97.5"]
      )
    )
  })

  return(results)
}

table2 <- create_main_effects_table(
  list(m1_f0mean, m2_f0std, m3_jitter, m4_nne),
  c("F0 Mean", "F0 SD", "Jitter", "NNE")
)

print(table2)
write_csv(table2, "results/table2_main_effects.csv")

cat("\n✓ Table 2 saved to results/table2_main_effects.csv\n")

# ==============================================================================
# 3. TABLE 3: MODERATION EFFECTS SUMMARY (Enhanced)
# ==============================================================================

cat("\n=== Creating Table 3: Moderation Effects (Enhanced) ===\n")

create_moderation_table_enhanced <- function(models, outcome_names) {
  results <- map2_dfr(models, outcome_names, function(model, outcome) {
    fixed <- fixef(model)

    # Define interactions
    interactions <- c(
      "c1_stress:pid5_na_bl_c" = "Stress × NA",
      "c1_stress:pid5_dis_bl_c" = "Stress × DIS",
      "c1_stress:pid5_psy_bl_c" = "Stress × PSY",
      "c2_recovery:pid5_na_bl_c" = "Recovery × NA",
      "c2_recovery:pid5_dis_bl_c" = "Recovery × DIS",
      "c2_recovery:pid5_psy_bl_c" = "Recovery × PSY"
    )

    # Extract for each interaction
    map_dfr(names(interactions), function(param) {
      if (param %in% rownames(fixed)) {
        # Get posterior for additional metrics
        post <- as.data.frame(model)[, paste0("b_", param)]
        pd <- bayestestR::pd(post)
        rope <- bayestestR::rope(post, range = c(-0.1, 0.1), ci = 0.95)

        tibble(
          Outcome = outcome,
          Interaction = interactions[param],
          β = fixed[param, "Estimate"],
          SE = fixed[param, "Est.Error"],
          `95% CI Lower` = fixed[param, "Q2.5"],
          `95% CI Upper` = fixed[param, "Q97.5"],
          `PD (%)` = pd * 100,
          `% in ROPE` = rope$ROPE_Percentage
        )
      }
    })
  })

  return(results)
}

table3 <- create_moderation_table_enhanced(
  list(m1_f0mean_mod, m2_f0std_mod, m3_jitter_mod, m4_nne_mod),
  c("F0 Mean", "F0 SD", "Jitter", "NNE")
)

# Format for publication
table3_formatted <- table3 %>%
  mutate(
    `β [95% CI]` = sprintf(
      "%.3f [%.3f, %.3f]",
      β,
      `95% CI Lower`,
      `95% CI Upper`
    ),
    `PD` = sprintf("%.1f%%", `PD (%)`),
    Evidence = case_when(
      `PD (%)` > 97.5 & `95% CI Lower` > 0 ~ "Strong positive",
      `PD (%)` > 97.5 & `95% CI Upper` < 0 ~ "Strong negative",
      `PD (%)` > 95 ~ "Moderate",
      TRUE ~ "Weak/Inconclusive"
    )
  ) %>%
  select(Outcome, Interaction, `β [95% CI]`, PD, Evidence)

print(table3_formatted)
write_csv(table3_formatted, "results/table3_moderation_effects.csv")

cat("\n✓ Table 3 saved to results/table3_moderation_effects.csv\n")

# ==============================================================================
# 4. DIAGNOSTIC PLOTS
# ==============================================================================

cat("\n=== Creating Diagnostic Plots ===\n")

# Function to create comprehensive diagnostics
create_diagnostics_plot <- function(model, model_name) {
  # 1. Trace plots
  p_trace <- mcmc_trace(model, pars = vars(starts_with("b_"))) +
    labs(title = paste(model_name, "- Trace Plots"))

  # 2. R-hat values
  rhat_vals <- rhat(model)
  p_rhat <- mcmc_rhat(rhat_vals) +
    labs(title = paste(model_name, "- R-hat Diagnostics"))

  # 3. Effective sample size
  neff_vals <- neff_ratio(model)
  p_neff <- mcmc_neff(neff_vals) +
    labs(title = paste(model_name, "- Effective Sample Size"))

  # 4. Posterior predictive checks
  p_pp <- pp_check(model, ndraws = 100) +
    labs(title = paste(model_name, "- Posterior Predictive Check"))

  # Combine
  p_combined <- (p_trace | p_rhat) /
    (p_neff | p_pp) +
    plot_annotation(
      title = paste("Diagnostic Plots:", model_name),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  return(p_combined)
}

# Create diagnostics for key models
diag_f0mean <- create_diagnostics_plot(m1_f0mean_mod, "F0 Mean Moderation")
ggsave(
  "figures/diagnostics_f0mean.png",
  diag_f0mean,
  width = 14,
  height = 10,
  dpi = 300
)

diag_f0std <- create_diagnostics_plot(m2_f0std_mod, "F0 SD Moderation")
ggsave(
  "figures/diagnostics_f0std.png",
  diag_f0std,
  width = 14,
  height = 10,
  dpi = 300
)

cat("\n✓ Diagnostic plots saved to figures/diagnostics_*.png\n")

# Print diagnostic summary
cat("\n=== Diagnostic Summary ===\n")

models_list <- list(
  "F0 Mean Main" = m1_f0mean,
  "F0 SD Main" = m2_f0std,
  "Jitter Main" = m3_jitter,
  "NNE Main" = m4_nne,
  "F0 Mean Mod" = m1_f0mean_mod,
  "F0 SD Mod" = m2_f0std_mod,
  "Jitter Mod" = m3_jitter_mod,
  "NNE Mod" = m4_nne_mod
)

diagnostic_summary <- map_dfr(names(models_list), function(name) {
  model <- models_list[[name]]

  tibble(
    Model = name,
    `Max R-hat` = max(rhat(model), na.rm = TRUE),
    `Min ESS ratio` = min(neff_ratio(model), na.rm = TRUE),
    `All R-hat < 1.01` = all(rhat(model) < 1.01, na.rm = TRUE),
    `All ESS > 0.1` = all(neff_ratio(model) > 0.1, na.rm = TRUE)
  )
})

print(diagnostic_summary)
write_csv(diagnostic_summary, "results/diagnostic_summary.csv")

# ==============================================================================
# 5. SENSITIVITY ANALYSES
# ==============================================================================

cat("\n=== Running Sensitivity Analyses ===\n")

# Sensitivity 1: Alternative priors (more skeptical for interactions)
cat("\nSensitivity 1: Skeptical priors for interactions...\n")

priors_skeptical <- c(
  prior(student_t(3, 0, 10), class = "Intercept"),
  prior(normal(0, 5), class = "b", coef = "c1_stress"),
  prior(normal(0, 5), class = "b", coef = "c2_recovery"),
  prior(normal(0, 3), class = "b", coef = "pid5_na_bl_c"),
  prior(normal(0, 3), class = "b", coef = "pid5_dis_bl_c"),
  prior(normal(0, 3), class = "b", coef = "pid5_psy_bl_c"),
  # More skeptical priors for interactions (narrower)
  prior(normal(0, 1), class = "b", coef = "c1_stress:pid5_na_bl_c"),
  prior(normal(0, 1), class = "b", coef = "c2_recovery:pid5_na_bl_c"),
  prior(normal(0, 1), class = "b", coef = "c1_stress:pid5_dis_bl_c"),
  prior(normal(0, 1), class = "b", coef = "c2_recovery:pid5_dis_bl_c"),
  prior(normal(0, 1), class = "b", coef = "c1_stress:pid5_psy_bl_c"),
  prior(normal(0, 1), class = "b", coef = "c2_recovery:pid5_psy_bl_c"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# Fit one model as example (F0 mean)
m1_f0mean_skeptical <- brm(
  formula_f0mean_mod,
  data = df_analysis,
  family = gaussian(),
  prior = priors_skeptical,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 456,
  control = list(adapt_delta = 0.95),
  file = "models/m1_f0mean_skeptical"
)

cat("\nComparing original vs skeptical priors:\n")
comparison <- tibble(
  Parameter = c("Stress × NA", "Stress × DIS", "Stress × PSY"),
  Original = c(
    fixef(m1_f0mean_mod)["c1_stress:pid5_na_bl_c", "Estimate"],
    fixef(m1_f0mean_mod)["c1_stress:pid5_dis_bl_c", "Estimate"],
    fixef(m1_f0mean_mod)["c1_stress:pid5_psy_bl_c", "Estimate"]
  ),
  Skeptical = c(
    fixef(m1_f0mean_skeptical)["c1_stress:pid5_na_bl_c", "Estimate"],
    fixef(m1_f0mean_skeptical)["c1_stress:pid5_dis_bl_c", "Estimate"],
    fixef(m1_f0mean_skeptical)["c1_stress:pid5_psy_bl_c", "Estimate"]
  ),
  Difference = abs(Original - Skeptical)
)

print(comparison)
write_csv(comparison, "results/sensitivity_priors.csv")

cat(
  "\n✓ Sensitivity to priors: estimates",
  ifelse(max(comparison$Difference) < 0.1, "ROBUST", "SENSITIVE"),
  "\n"
)

# Sensitivity 2: Exclude potential outliers
cat("\nSensitivity 2: Outlier influence...\n")

# Identify potential outliers (>3 SD from mean)
outliers <- df_analysis %>%
  group_by(timepoint) %>%
  mutate(
    across(
      c(f0_mean, f0_std, jitter, nne),
      list(z = ~ abs(scale(.))),
      .names = "{.col}_z"
    )
  ) %>%
  filter(if_any(ends_with("_z"), ~ . > 3))

cat("Number of potential outlier observations:", nrow(outliers), "\n")
cat(
  "Percentage:",
  sprintf("%.2f%%", nrow(outliers) / nrow(df_analysis) * 100),
  "\n"
)

if (nrow(outliers) > 0) {
  # Refit without outliers
  df_no_outliers <- df_analysis %>%
    anti_join(outliers, by = c("ID", "timepoint"))

  m1_f0mean_no_outliers <- update(
    m1_f0mean_mod,
    newdata = df_no_outliers,
    file = "models/m1_f0mean_no_outliers"
  )

  # Compare
  cat("\nKey interaction estimates (with vs without outliers):\n")
  comparison_outliers <- tibble(
    Parameter = c("Stress × NA", "Stress × DIS", "Stress × PSY"),
    `Full Sample` = c(
      fixef(m1_f0mean_mod)["c1_stress:pid5_na_bl_c", "Estimate"],
      fixef(m1_f0mean_mod)["c1_stress:pid5_dis_bl_c", "Estimate"],
      fixef(m1_f0mean_mod)["c1_stress:pid5_psy_bl_c", "Estimate"]
    ),
    `No Outliers` = c(
      fixef(m1_f0mean_no_outliers)["c1_stress:pid5_na_bl_c", "Estimate"],
      fixef(m1_f0mean_no_outliers)["c1_stress:pid5_dis_bl_c", "Estimate"],
      fixef(m1_f0mean_no_outliers)["c1_stress:pid5_psy_bl_c", "Estimate"]
    )
  )

  print(comparison_outliers)
  write_csv(comparison_outliers, "results/sensitivity_outliers.csv")
}

# ==============================================================================
# 6. SUPPLEMENTARY FIGURES
# ==============================================================================

cat("\n=== Creating Supplementary Figures ===\n")

# Supplementary Figure: Individual trajectories by PID-5 tertiles
create_individual_trajectories <- function(
  data,
  outcome_var,
  pid5_var,
  outcome_label,
  pid5_label
) {
  # Create tertiles of PID-5
  data_plot <- data %>%
    group_by(ID) %>%
    mutate(
      pid5_tertile = ntile({{ pid5_var }}, 3)
    ) %>%
    ungroup() %>%
    mutate(
      pid5_group = factor(
        pid5_tertile,
        levels = 1:3,
        labels = c("Low", "Medium", "High")
      )
    )

  # Plot individual trajectories
  ggplot(
    data_plot,
    aes(x = timepoint, y = {{ outcome_var }}, group = ID, color = pid5_group)
  ) +
    geom_line(alpha = 0.3) +
    stat_summary(
      aes(group = pid5_group),
      fun = mean,
      geom = "line",
      size = 1.5
    ) +
    stat_summary(
      aes(group = pid5_group),
      fun = mean,
      geom = "point",
      size = 3
    ) +
    facet_wrap(~pid5_group, ncol = 3) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = paste(outcome_label, "by", pid5_label, "Levels"),
      subtitle = "Individual trajectories (thin lines) and group means (thick lines)",
      x = "Timepoint",
      y = outcome_label,
      color = pid5_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )
}

p_supp1 <- create_individual_trajectories(
  df_analysis,
  f0_mean,
  pid5_na_bl_c,
  "F0 Mean",
  "Negative Affectivity"
)

ggsave(
  "figures/supplementary_trajectories_NA.png",
  p_supp1,
  width = 14,
  height = 5,
  dpi = 300
)

p_supp2 <- create_individual_trajectories(
  df_analysis,
  f0_std,
  pid5_dis_bl_c,
  "F0 SD",
  "Disinhibition"
)

ggsave(
  "figures/supplementary_trajectories_DIS.png",
  p_supp2,
  width = 14,
  height = 5,
  dpi = 300
)

cat("\n✓ Supplementary figures saved\n")

# ==============================================================================
# 7. EFFECT SIZE CALCULATIONS
# ==============================================================================

cat("\n=== Calculating Standardized Effect Sizes ===\n")

# Calculate within-person effect sizes (Cohen's d)
calculate_within_d <- function(data, outcome_var, time1, time2) {
  data_wide <- data %>%
    filter(timepoint %in% c(time1, time2)) %>%
    select(ID, timepoint, {{ outcome_var }}) %>%
    pivot_wider(names_from = timepoint, values_from = {{ outcome_var }})

  # Calculate difference scores
  diff <- data_wide[[time2]] - data_wide[[time1]]

  # Cohen's d for within-person design
  d <- mean(diff, na.rm = TRUE) / sd(diff, na.rm = TRUE)

  return(d)
}

effect_sizes <- tibble(
  Outcome = rep(c("F0 Mean", "F0 SD", "Jitter", "NNE"), each = 2),
  Contrast = rep(c("PRE vs BASE", "POST vs PRE"), 4),
  `Cohen's d` = c(
    calculate_within_d(df_analysis, f0_mean, "baseline", "pre"),
    calculate_within_d(df_analysis, f0_mean, "pre", "post"),
    calculate_within_d(df_analysis, f0_std, "baseline", "pre"),
    calculate_within_d(df_analysis, f0_std, "pre", "post"),
    calculate_within_d(df_analysis, jitter, "baseline", "pre"),
    calculate_within_d(df_analysis, jitter, "pre", "post"),
    calculate_within_d(df_analysis, nne, "baseline", "pre"),
    calculate_within_d(df_analysis, nne, "pre", "post")
  )
) %>%
  mutate(
    Magnitude = case_when(
      abs(`Cohen's d`) < 0.2 ~ "Negligible",
      abs(`Cohen's d`) < 0.5 ~ "Small",
      abs(`Cohen's d`) < 0.8 ~ "Medium",
      TRUE ~ "Large"
    )
  )

print(effect_sizes)
write_csv(effect_sizes, "results/effect_sizes.csv")

# ==============================================================================
# SAVE ALL
# ==============================================================================

save.image("results/utilities_complete.RData")

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("UTILITIES COMPLETE\n")
cat(rep("=", 80), "\n", sep = "")
cat("\n✓ All supplementary analyses complete\n")
cat("✓ Tables formatted and saved\n")
cat("✓ Diagnostics created\n")
cat("✓ Sensitivity analyses conducted\n")
cat("✓ Effect sizes calculated\n\n")
