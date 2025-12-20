# ==============================================================================
# 2. SET UP PRIORS FOR MODERATION MODELS
# ==============================================================================

# Weakly informative priors for moderation with ALL 5 PID-5 domains
priors_moderation <- c(
  prior(student_t(3, 0, 10), class = "Intercept"),
  prior(normal(0, 5), class = "b"), # Applies to all fixed effects
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# ==============================================================================
# 3. FIT MODERATION MODELS
# ==============================================================================

cat("\n=== FITTING MODERATION MODELS (ALL 5 PID-5 DOMAINS) ===\n")
cat(
  "Domains: Negative Affectivity, Detachment, Antagonism, Disinhibition, Psychoticism\n"
)
cat("This will take ~15-20 minutes per model...\n\n")

# ------------------------------------------------------------------------------
# Model 1: F0 Mean (vowel /a/) ~ (Stress Contrasts) × PID-5 (all 5)
# ------------------------------------------------------------------------------

cat("Fitting: F0 Mean (/a/) moderation model...\n")

formula_f0mean_mod_a <- bf(
  f0_mean_a ~
    c1_stress *
      (pid5_negative_affectivity_bl_c +
        pid5_detachment_bl_c +
        pid5_antagonism_bl_c +
        pid5_disinhibition_bl_c +
        pid5_psychoticism_bl_c) +
      c2_recovery *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
      (1 + c1_stress + c2_recovery | ID)
)

m1_f0mean_mod_a <- brm(
  formula_f0mean_mod_a,
  data = df_analysis,
  family = gaussian(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
  # file = "models/m1_f0mean_a_moderation_5domains"
)
pp_check(m1_f0mean_mod_a)
summary(m1_f0mean_mod_a)

# ------------------------------------------------------------------------------
# Model 2: F0 SD (vowel /a/) ~ (Stress Contrasts) × PID-5 (all 5)
# ------------------------------------------------------------------------------

cat("Fitting: F0 SD (/a/) moderation model...\n")

formula_f0std_mod_a <- bf(
  f0_std_a ~
    c1_stress *
      (pid5_negative_affectivity_bl_c +
        pid5_detachment_bl_c +
        pid5_antagonism_bl_c +
        pid5_disinhibition_bl_c +
        pid5_psychoticism_bl_c) +
      c2_recovery *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
      (1 + c1_stress + c2_recovery | ID)
)

m2_f0std_mod_a <- brm(
  formula_f0std_mod_a,
  data = df_analysis,
  family = lognormal(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
  # file = "models/m2_f0std_a_moderation_5domains"
)
pp_check(m2_f0std_mod_a)
summary(m2_f0std_mod_a)


# ------------------------------------------------------------------------------
# Model 3: Jitter (vowel /a/) ~ (Stress Contrasts) × PID-5 (all 5)
# ------------------------------------------------------------------------------

cat("Fitting: Jitter (/a/) moderation model...\n")

formula_jitter_mod_a <- bf(
  jitter_a ~
    c1_stress *
      (pid5_negative_affectivity_bl_c +
        pid5_detachment_bl_c +
        pid5_antagonism_bl_c +
        pid5_disinhibition_bl_c +
        pid5_psychoticism_bl_c) +
      c2_recovery *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
      (1 + c1_stress + c2_recovery | ID)
)

m3_jitter_mod_a <- brm(
  formula_jitter_mod_a,
  data = df_analysis,
  family = lognormal(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
  # file = "models/m3_jitter_a_moderation_5domains"
)
pp_check(m3_jitter_mod_a)
summary(m3_jitter_mod_a)


# ------------------------------------------------------------------------------
# Model 4: NNE (vowel /a/) ~ (Stress Contrasts) × PID-5 (all 5)
# ------------------------------------------------------------------------------

cat("Fitting: NNE (/a/) moderation model...\n")

formula_nne_mod_a <- bf(
  nne_a ~
    c1_stress *
      (pid5_negative_affectivity_bl_c +
        pid5_detachment_bl_c +
        pid5_antagonism_bl_c +
        pid5_disinhibition_bl_c +
        pid5_psychoticism_bl_c) +
      c2_recovery *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
      (1 + c1_stress + c2_recovery | ID)
)

m4_nne_mod_a <- brm(
  formula_nne_mod_a,
  data = df_analysis,
  family = gaussian(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
  # file = "models/m4_nne_a_moderation_5domains"
)
pp_check(m4_nne_mod_a)
summary(m4_nne_mod_a)

# ------------------------------------------------------------------------------
# Model 5: F2 Mean (vowel /a/) ~ (Stress Contrasts) × PID-5 (all 5)
# ------------------------------------------------------------------------------

cat("Fitting: F2 Mean (/a/) moderation model...\n")

formula_f2mean_mod_a <- bf(
  f2_mean_a ~
    c1_stress *
      (pid5_negative_affectivity_bl_c +
        pid5_detachment_bl_c +
        pid5_antagonism_bl_c +
        pid5_disinhibition_bl_c +
        pid5_psychoticism_bl_c) +
      c2_recovery *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
      (1 + c1_stress + c2_recovery | ID)
)

m5_f2mean_mod_a <- brm(
  formula_f2mean_mod_a,
  data = df_analysis,
  family = asym_laplace(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
  # file = "models/m1_f0mean_a_moderation_5domains"
)
pp_check(m5_f2mean_mod_a)
summary(m5_f2mean_mod_a)

# ------------------------------------------------------------------------------
# Model 6: F2 SD (vowel /a/) ~ (Stress Contrasts) × PID-5 (all 5)
# ------------------------------------------------------------------------------

cat("Fitting: F2 SD (/a/) moderation model...\n")

formula_f2std_mod_a <- bf(
  f2_std_a ~
    c1_stress *
      (pid5_negative_affectivity_bl_c +
        pid5_detachment_bl_c +
        pid5_antagonism_bl_c +
        pid5_disinhibition_bl_c +
        pid5_psychoticism_bl_c) +
      c2_recovery *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
      (1 + c1_stress + c2_recovery | ID)
)

m6_f2std_mod_a <- brm(
  formula_f2std_mod_a,
  data = df_analysis,
  family = lognormal(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
  # file = "models/m2_f0std_a_moderation_5domains"
)
pp_check(m2_f0std_mod_a)
summary(m2_f0std_mod_a)

cat("\n✓ All moderation models fitted successfully\n")

# ==============================================================================
# 4. EXTRACT AND SUMMARIZE MODERATION RESULTS
# ==============================================================================

cat("\n=== MODERATION RESULTS (ALL 5 PID-5 DOMAINS) ===\n")

# Function to extract interaction effects
summarize_moderation <- function(model, outcome_name) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("OUTCOME:", outcome_name, "\n")
  cat(rep("=", 70), "\n", sep = "")

  # Full summary
  cat("\nFixed Effects Summary:\n")
  print(summary(model)$fixed, digits = 3)

  # Focus on interaction effects - STRESS REACTIVITY
  cat("\n\n=== STRESS REACTIVITY INTERACTIONS (PRE vs BASELINE) ===\n")

  interactions_stress <- c(
    "c1_stress:pid5_negative_affectivity_bl_c" = "Negative Affectivity",
    "c1_stress:pid5_detachment_bl_c" = "Detachment",
    "c1_stress:pid5_antagonism_bl_c" = "Antagonism",
    "c1_stress:pid5_disinhibition_bl_c" = "Disinhibition",
    "c1_stress:pid5_psychoticism_bl_c" = "Psychoticism"
  )

  for (i in seq_along(interactions_stress)) {
    param <- names(interactions_stress)[i]
    trait_name <- interactions_stress[i]

    h <- hypothesis(model, paste0(param, " = 0"))
    cat("\n", trait_name, ":\n", sep = "")
    print(h$hypothesis, digits = 3)

    # Probability of direction
    post_sample <- as.data.frame(model)[, paste0("b_", param)]
    pd <- bayestestR::pd(post_sample)
    cat("Probability of Direction:", sprintf("%.1f%%", pd * 100), "\n")

    # Practical significance (ROPE: ±0.1)
    rope_result <- bayestestR::rope(
      post_sample,
      range = c(-0.1, 0.1),
      ci = 0.95
    )
    cat(
      "% in ROPE (±0.1):",
      sprintf("%.1f%%", rope_result$ROPE_Percentage),
      "\n"
    )
  }

  # Focus on interaction effects - RECOVERY
  cat("\n\n=== RECOVERY INTERACTIONS (POST vs PRE) ===\n")

  interactions_recovery <- c(
    "c2_recovery:pid5_negative_affectivity_bl_c" = "Negative Affectivity",
    "c2_recovery:pid5_detachment_bl_c" = "Detachment",
    "c2_recovery:pid5_antagonism_bl_c" = "Antagonism",
    "c2_recovery:pid5_disinhibition_bl_c" = "Disinhibition",
    "c2_recovery:pid5_psychoticism_bl_c" = "Psychoticism"
  )

  for (i in seq_along(interactions_recovery)) {
    param <- names(interactions_recovery)[i]
    trait_name <- interactions_recovery[i]

    h <- hypothesis(model, paste0(param, " = 0"))
    cat("\n", trait_name, ":\n", sep = "")
    print(h$hypothesis, digits = 3)

    # Probability of direction
    post_sample <- as.data.frame(model)[, paste0("b_", param)]
    pd <- bayestestR::pd(post_sample)
    cat("Probability of Direction:", sprintf("%.1f%%", pd * 100), "\n")

    # Practical significance
    rope_result <- bayestestR::rope(
      post_sample,
      range = c(-0.1, 0.1),
      ci = 0.95
    )
    cat(
      "% in ROPE (±0.1):",
      sprintf("%.1f%%", rope_result$ROPE_Percentage),
      "\n"
    )
  }

  # Model diagnostics
  cat("\n\n=== MODEL DIAGNOSTICS ===\n")
  cat("All R-hat < 1.01:", all(rhat(model) < 1.01, na.rm = TRUE), "\n")
  cat("All ESS ratio > 0.1:", all(neff_ratio(model) > 0.1, na.rm = TRUE), "\n")

  return(invisible(model))
}

# Summarize each model
summarize_moderation(m1_f0mean_mod, "F0 Mean /a/ (Hz)")
summarize_moderation(m2_f0std_mod, "F0 SD /a/ (Hz)")
summarize_moderation(m3_jitter_mod, "Jitter /a/")
summarize_moderation(m4_nne_mod, "NNE /a/")

# ==============================================================================
# 5. VISUALIZE INTERACTIONS: CONDITIONAL EFFECTS
# ==============================================================================

cat("\n=== CREATING CONDITIONAL EFFECTS PLOTS ===\n")

# Function to plot simple slopes for key interactions
plot_conditional_slopes <- function(
  model,
  outcome_name,
  moderator_var,
  moderator_label
) {
  # Get conditional effects
  cond_effects <- conditional_effects(
    model,
    effects = paste0("c1_stress:", moderator_var),
    int_conditions = list(across = moderator_var, values = c(-1, 0, 1))
  )

  # Extract data for plotting
  plot_data <- cond_effects[[1]] %>%
    mutate(
      moderator_level = factor(
        cond__.across__,
        levels = c(-1, 0, 1),
        labels = c("Low (-1 SD)", "Mean", "High (+1 SD)")
      ),
      timepoint_label = case_when(
        c1_stress == -0.5 ~ "Baseline",
        c1_stress == 0.5 ~ "Pre-Exam",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(timepoint_label))

  # Create plot
  p <- ggplot(
    plot_data,
    aes(
      x = timepoint_label,
      y = estimate__,
      color = moderator_level,
      group = moderator_level
    )
  ) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_ribbon(
      aes(ymin = lower__, ymax = upper__, fill = moderator_level),
      alpha = 0.2,
      color = NA
    ) +
    scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73")) +
    scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73")) +
    labs(
      title = paste(outcome_name, "×", moderator_label),
      subtitle = "Stress reactivity moderated by personality pathology",
      x = "Timepoint",
      y = outcome_name,
      color = moderator_label,
      fill = moderator_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  return(p)
}

# Create plots for ALL 5 PID-5 domains (focusing on F0 mean and F0 SD)
cat("Creating interaction plots for all 5 PID-5 domains...\n")

# F0 Mean - all 5 domains
p_f0mean_na <- plot_conditional_slopes(
  m1_f0mean_mod,
  "F0 Mean",
  "pid5_negative_affectivity_bl_c",
  "Negative Affectivity"
)
p_f0mean_det <- plot_conditional_slopes(
  m1_f0mean_mod,
  "F0 Mean",
  "pid5_detachment_bl_c",
  "Detachment"
)
p_f0mean_ant <- plot_conditional_slopes(
  m1_f0mean_mod,
  "F0 Mean",
  "pid5_antagonism_bl_c",
  "Antagonism"
)
p_f0mean_dis <- plot_conditional_slopes(
  m1_f0mean_mod,
  "F0 Mean",
  "pid5_disinhibition_bl_c",
  "Disinhibition"
)
p_f0mean_psy <- plot_conditional_slopes(
  m1_f0mean_mod,
  "F0 Mean",
  "pid5_psychoticism_bl_c",
  "Psychoticism"
)

# F0 SD - all 5 domains
p_f0std_na <- plot_conditional_slopes(
  m2_f0std_mod,
  "F0 SD",
  "pid5_negative_affectivity_bl_c",
  "Negative Affectivity"
)
p_f0std_det <- plot_conditional_slopes(
  m2_f0std_mod,
  "F0 SD",
  "pid5_detachment_bl_c",
  "Detachment"
)
p_f0std_ant <- plot_conditional_slopes(
  m2_f0std_mod,
  "F0 SD",
  "pid5_antagonism_bl_c",
  "Antagonism"
)
p_f0std_dis <- plot_conditional_slopes(
  m2_f0std_mod,
  "F0 SD",
  "pid5_disinhibition_bl_c",
  "Disinhibition"
)
p_f0std_psy <- plot_conditional_slopes(
  m2_f0std_mod,
  "F0 SD",
  "pid5_psychoticism_bl_c",
  "Psychoticism"
)

# Combine plots - F0 Mean (all 5 domains)
p_f0mean_all <- (p_f0mean_na + p_f0mean_det + p_f0mean_ant) /
  (p_f0mean_dis + p_f0mean_psy + plot_spacer()) +
  plot_annotation(
    title = "F0 Mean: All 5 PID-5 Domains × Stress Reactivity",
    subtitle = "Baseline → Pre-Exam (95% CI)",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(
  "figures/moderation_f0mean_all5domains.png",
  p_f0mean_all,
  width = 18,
  height = 12,
  dpi = 300
)

# Combine plots - F0 SD (all 5 domains)
p_f0std_all <- (p_f0std_na + p_f0std_det + p_f0std_ant) /
  (p_f0std_dis + p_f0std_psy + plot_spacer()) +
  plot_annotation(
    title = "F0 SD: All 5 PID-5 Domains × Stress Reactivity",
    subtitle = "Baseline → Pre-Exam (95% CI)",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(
  "figures/moderation_f0std_all5domains.png",
  p_f0std_all,
  width = 18,
  height = 12,
  dpi = 300
)

cat("\n✓ Interaction plots saved\n")

# ==============================================================================
# 6. POSTERIOR DISTRIBUTIONS OF INTERACTIONS
# ==============================================================================

cat("\n=== Creating posterior distribution plots ===\n")

# Function to visualize posterior distributions
plot_interaction_posteriors <- function(models_list, model_names) {
  all_draws <- map2_dfr(models_list, model_names, function(model, name) {
    draws <- as_draws_df(model)

    draws %>%
      select(
        starts_with("b_c1_stress:pid5"),
        starts_with("b_c2_recovery:pid5")
      ) %>%
      pivot_longer(
        everything(),
        names_to = "parameter",
        values_to = "value"
      ) %>%
      mutate(
        outcome = name,
        parameter = str_remove(parameter, "b_"),
        interaction_type = ifelse(
          str_detect(parameter, "c1_stress"),
          "Stress Reactivity",
          "Recovery"
        ),
        pid5_trait = case_when(
          str_detect(parameter, "negative_affectivity") ~
            "Negative Affectivity",
          str_detect(parameter, "detachment") ~ "Detachment",
          str_detect(parameter, "antagonism") ~ "Antagonism",
          str_detect(parameter, "disinhibition") ~ "Disinhibition",
          str_detect(parameter, "psychoticism") ~ "Psychoticism"
        )
      )
  })

  # Plot
  ggplot(all_draws, aes(x = value, y = pid5_trait, fill = outcome)) +
    stat_halfeye(alpha = 0.7, .width = c(0.95, 0.66)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    facet_grid(interaction_type ~ outcome, scales = "free_x") +
    labs(
      title = "Posterior Distributions: All 5 PID-5 Domains",
      subtitle = "How personality pathology moderates acoustic stress responses",
      x = "Interaction coefficient (standardized)",
      y = "PID-5 Domain",
      fill = "Acoustic\nOutcome"
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

models_list <- list(m1_f0mean_mod, m2_f0std_mod, m3_jitter_mod, m4_nne_mod)
model_names <- c("F0 Mean", "F0 SD", "Jitter", "NNE")

p_posteriors <- plot_interaction_posteriors(models_list, model_names)

ggsave(
  "figures/interaction_posteriors_all5domains.png",
  p_posteriors,
  width = 18,
  height = 14,
  dpi = 300
)

cat("\n✓ Posterior distributions saved\n")

# ==============================================================================
# 7. SUMMARY TABLE FOR MANUSCRIPT
# ==============================================================================

cat("\n=== Creating summary table ===\n")

# Create publication-ready summary table
create_results_table <- function(models_list, outcome_names) {
  results <- map2_dfr(models_list, outcome_names, function(model, outcome) {
    fixed_effects <- fixef(model)

    # Extract key interactions
    interactions <- fixed_effects %>%
      as.data.frame() %>%
      rownames_to_column("parameter") %>%
      filter(str_detect(parameter, ":")) %>%
      mutate(
        outcome = outcome,
        interaction_type = ifelse(
          str_detect(parameter, "c1_stress"),
          "Reactivity",
          "Recovery"
        ),
        pid5_trait = case_when(
          str_detect(parameter, "negative_affectivity") ~ "NA",
          str_detect(parameter, "detachment") ~ "DET",
          str_detect(parameter, "antagonism") ~ "ANT",
          str_detect(parameter, "disinhibition") ~ "DIS",
          str_detect(parameter, "psychoticism") ~ "PSY"
        )
      ) %>%
      select(
        outcome,
        interaction_type,
        pid5_trait,
        Estimate,
        `Est.Error`,
        `Q2.5`,
        `Q97.5`
      )
  })

  return(results)
}

results_table <- create_results_table(models_list, model_names)

cat("\n=== SUMMARY TABLE: Moderation Effects (All 5 Domains) ===\n")
print(results_table, n = Inf)

# Save to CSV
write_csv(results_table, "results/moderation_summary_table_5domains.csv")
cat(
  "\n✓ Summary table saved to results/moderation_summary_table_5domains.csv\n"
)

# ==============================================================================
# 8. SAVE ALL RESULTS
# ==============================================================================

save.image("results/moderation_analysis_5domains_complete.RData")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS COMPLETE (ALL 5 PID-5 DOMAINS)\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n✓ All models fitted and results saved\n")
cat("✓ Figures created in figures/ directory\n")
cat(
  "✓ Workspace saved to results/moderation_analysis_5domains_complete.RData\n\n"
)
cat("Key findings ready for manuscript:\n")
cat("1. Moderation by ALL 5 PID-5 domains\n")
cat("2. Stress reactivity (PRE vs BASELINE) interactions\n")
cat("3. Recovery (POST vs PRE) interactions\n")
cat("4. Posterior distributions and credible intervals\n")
cat("5. Publication-ready summary table\n\n")
