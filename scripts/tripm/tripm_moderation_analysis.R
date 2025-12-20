# ==============================================================================
# Bayesian Moderation Analysis: TriPM Traits × Stress Reactivity
# Testing context-dependent expression of psychopathic traits
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayestestR)
  library(bayesplot)
  library(tidybayes)
  library(patchwork)
  library(ggdist)
  library(marginaleffects)
  library(tidyr)
  library(stringr)
})

options(brms.backend = "cmdstanr")
options(max.print = 5000)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS: TriPM × Stress Reactivity\n")
cat("Boldness, Meanness, Disinhibition\n")
cat(rep("=", 70), "\n", sep = "")

# ==============================================================================
# 1. PREPARE DATA FOR MODERATION MODELS - TriPM VERSION
# ==============================================================================

# Extract baseline TriPM scores (for between-person predictors)
tripm_baseline <- df_clean %>%
  dplyr::filter(timepoint == "baseline") %>%
  dplyr::select(ID, tripm_boldness, tripm_meanness, pid5_disinhibition) %>%
  rename_with(
    ~ paste0(., "_bl"),
    .cols = c(tripm_boldness, tripm_meanness, pid5_disinhibition)
  )

# Create long format with baseline TriPM as between-person predictors
df_long_tripm <- df_clean %>%
  dplyr::select(-tripm_boldness, -tripm_meanness) %>% # Remove timepoint-varying TriPM
  left_join(tripm_baseline, by = "ID") %>%
  # Center baseline TriPM scores for interpretability
  mutate(
    across(
      ends_with("_bl"),
      ~ scale(., center = TRUE, scale = TRUE)[, 1],
      .names = "{.col}_c"
    )
  )

# Remove missing cases
df_analysis_tripm <- df_long_tripm %>%
  dplyr::filter(complete.cases(.))

# Check TriPM distributions at baseline
cat("\n=== TriPM Baseline Trait Distributions ===\n")
df_analysis_tripm %>%
  dplyr::filter(timepoint == "baseline") %>%
  dplyr::select(ends_with("_bl_c")) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    M = mean(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Min = min(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# Create contrast-coded timepoint variable for interaction interpretation
# C1: PRE vs BASELINE (stress reactivity)
# C2: POST vs PRE (stress recovery)

df_analysis_tripm <- df_analysis_tripm %>%
  mutate(
    # Contrast 1: PRE vs BASELINE (test stress reactivity)
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre" ~ 0.5,
      timepoint == "post" ~ 0
    ),
    # Contrast 2: POST vs PRE (test recovery)
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre" ~ -0.5,
      timepoint == "post" ~ 0.5
    )
  )

# ==============================================================================
# FIX ZERO VALUES FOR LOGNORMAL OUTCOMES (AFTER DATA CREATION!)
# ==============================================================================

cat("\n=== Checking for zero values in lognormal outcomes ===\n")

# Variables that use lognormal family
lognormal_vars <- c(
  "f0_std_a",
  "f0_std_i",
  "f0_std_u",
  "jitter_a",
  "jitter_i",
  "jitter_u",
  "f2_std_a",
  "f2_std_i",
  "f2_std_u"
)

# Check zeros in each variable
zero_counts <- sapply(lognormal_vars, function(var) {
  sum(df_analysis_tripm[[var]] == 0, na.rm = TRUE)
})

print(data.frame(
  Variable = lognormal_vars,
  N_zeros = zero_counts,
  Percent = round(100 * zero_counts / nrow(df_analysis_tripm), 2)
))

# Replace zeros with minimum positive value (only if < 5% are zeros)
for (var in lognormal_vars) {
  n_zeros <- sum(df_analysis_tripm[[var]] == 0, na.rm = TRUE)

  if (n_zeros > 0) {
    pct_zeros <- 100 * n_zeros / nrow(df_analysis_tripm)

    if (pct_zeros < 5) {
      # Few zeros: replace with minimum positive value
      min_positive <- min(
        df_analysis_tripm[[var]][df_analysis_tripm[[var]] > 0],
        na.rm = TRUE
      )
      df_analysis_tripm[[var]] <- ifelse(
        df_analysis_tripm[[var]] == 0,
        min_positive,
        df_analysis_tripm[[var]]
      )
      cat(
        "  ",
        var,
        ": Replaced",
        n_zeros,
        "zeros with",
        round(min_positive, 4),
        "\n"
      )
    } else {
      # Many zeros: set as missing (more conservative)
      df_analysis_tripm[[var]] <- ifelse(
        df_analysis_tripm[[var]] == 0,
        NA_real_,
        df_analysis_tripm[[var]]
      )
      cat("  ", var, ": Set", n_zeros, "zeros to NA (>5% of data)\n")
    }
  }
}

cat("\n")

# ==============================================================================
# 2. DEFINE PRIORS AND MODEL SPECIFICATIONS
# ==============================================================================

# Vowel-specific empirical parameters (from your data)
vowel_means <- list(
  a = list(
    f0_mean = 189.9,
    f0_std = 6.230,
    jitter = 0.8366,
    nne = -24.42,
    f2_mean = 1207.4,
    f2_std = 59.37
  ),
  i = list(
    f0_mean = 194.4,
    f0_std = 8.888,
    jitter = 1.058,
    nne = -27.89,
    f2_mean = 2189.0,
    f2_std = 141.04
  ),
  u = list(
    f0_mean = 195.1,
    f0_std = 9.724,
    jitter = 1.286,
    nne = -28.429,
    f2_mean = 1123.2,
    f2_std = 161.57
  )
)

# Prior function
make_priors_vowel <- function(outcome, vowel) {
  vm <- vowel_means[[vowel]]
  fmt <- function(x, digits = 6) formatC(x, format = "f", digits = digits)

  switch(
    outcome,

    "f0_mean" = c(
      prior_string(
        paste0("student_t(3, ", fmt(vm$f0_mean), ", 30)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 10)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "f0_std" = c(
      prior_string(
        paste0("student_t(3, ", fmt(log(vm$f0_std)), ", 0.5)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 0.3)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "jitter" = c(
      prior_string(
        paste0("student_t(3, ", fmt(log(vm$jitter)), ", 0.5)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 0.3)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "nne" = c(
      prior_string(
        paste0("student_t(3, ", fmt(vm$nne), ", 5)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 2)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "f2_mean" = c(
      prior_string(
        paste0("student_t(3, ", fmt(vm$f2_mean), ", 150)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 50)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "f2_std" = c(
      prior_string(
        paste0("student_t(3, ", fmt(log(vm$f2_std)), ", 0.6)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 0.4)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    stop("Outcome non riconosciuto")
  )
}

# Outcome-family mapping
outcome_families <- list(
  f0_mean = gaussian(),
  f0_std = lognormal(),
  jitter = lognormal(),
  nne = gaussian(),
  f2_mean = asym_laplace(),
  f2_std = lognormal()
)

# ==============================================================================
# 3. FIT MODELS - TriPM VERSION
# ==============================================================================

# Define outcomes and vowels
outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
vowels <- c("a", "i", "u")

# Sampling parameters
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
seed <- 123
control <- list(adapt_delta = 0.95, max_treedepth = 12)

# Create models directory
dir.create("models_tripm", showWarnings = FALSE)

# Control flag
run_models <- TRUE

# Loop to fit models
fitted_models_tripm <- list()

for (v in vowels) {
  message("\n--- Vowel: /", v, "/ -------------------------------")

  for (out in outcomes) {
    # Column name in dataset
    colname <- paste0(out, "_", v)
    fam <- outcome_families[[out]]

    # Safety check
    if (!colname %in% names(df_analysis_tripm)) {
      message("  Skipping ", colname, " (column not found in dataset).")
      next
    }

    # Formula with TriPM dimensions (3 predictors instead of 5)
    fmla <- bf(
      as.formula(
        paste0(
          colname,
          " ~ c1_stress * (tripm_boldness_bl_c + tripm_meanness_bl_c + pid5_disinhibition_bl_c) + ",
          "c2_recovery * (tripm_boldness_bl_c + tripm_meanness_bl_c + pid5_disinhibition_bl_c) + ",
          "(1 + c1_stress + c2_recovery | ID)"
        )
      )
    )

    priors_here <- make_priors_vowel(out, v)
    model_name <- paste0("m_tripm_", out, "_", v)

    # Extract family name safely
    fam_name <- NULL
    if (!is.null(fam$family)) {
      fam_name <- fam$family
    } else if (!is.null(attr(fam, "family"))) {
      fam_name <- attr(fam, "family")
    } else {
      fam_name <- paste0(class(fam)[1])
    }

    message("  Model: ", model_name, "  (family: ", fam_name, ")")

    if (run_models) {
      # Fit model
      fit <- brm(
        formula = fmla,
        data = df_analysis_tripm,
        family = fam,
        prior = priors_here,
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores,
        seed = seed,
        control = control,
        file = file.path("models_tripm", model_name)
      )

      fitted_models_tripm[[model_name]] <- fit

      # Quick diagnostic check
      print(summary(fit, waic = FALSE, odds = FALSE))
      pp_check(fit)
    } else {
      message("    (run_models = FALSE) - priors that would be used:")
      print(priors_here)
    }
  }
}

# ==============================================================================
# 4. EXTRACT AND SUMMARIZE RESULTS
# ==============================================================================

# Extract summary from all models
results_table_tripm <- lapply(names(fitted_models_tripm), function(mn) {
  fit <- fitted_models_tripm[[mn]]
  fe <- fixef(fit, summary = TRUE)

  # Parse model name: m_tripm_f0_mean_a -> outcome = f0, measure = mean, vowel = a
  parts <- str_split(mn, "_", simplify = TRUE)
  outcome <- parts[3]
  if (length(parts) == 5) {
    measure <- parts[4]
    vowel <- parts[5]
    outcome_full <- paste0(outcome, "_", measure)
  } else {
    vowel <- parts[4]
    outcome_full <- outcome
  }

  # Identify term type
  terms <- rownames(fe)
  type <- ifelse(str_detect(terms, ":"), "Interaction", "Main")

  tibble(
    Outcome = outcome_full,
    Vowel = vowel,
    Term = terms,
    Type = type,
    Estimate = fe[, "Estimate"],
    SE = fe[, "Est.Error"],
    CI_low = fe[, "Q2.5"],
    CI_high = fe[, "Q97.5"],
    Credible = ifelse(fe[, "Q2.5"] > 0 | fe[, "Q97.5"] < 0, TRUE, FALSE)
  )
}) %>%
  bind_rows() %>%
  arrange(Outcome, Vowel, Type, Term)

# Display results table
cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("RESULTS TABLE: TriPM × Stress Reactivity Interactions\n")
cat(rep("=", 70), "\n", sep = "")
print(data.frame(results_table_tripm), row.names = FALSE)

# ==============================================================================
# 5. VISUALIZATIONS
# ==============================================================================

# Forest plot
cat("\n=== Creating Forest Plot ===\n")

p_forest <- ggplot(
  results_table_tripm,
  aes(
    x = Estimate,
    y = Term,
    xmin = CI_low,
    xmax = CI_high,
    color = Type,
    shape = Credible
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(Outcome ~ Vowel, scales = "free_y") +
  scale_color_manual(
    values = c("Main" = "blue", "Interaction" = "red"),
    name = "Effect Type"
  ) +
  scale_shape_manual(
    values = c(FALSE = 1, TRUE = 16),
    name = "Credible Effect",
    labels = c("No", "Yes")
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Posterior Estimate (95% CI)",
    y = "Parameter",
    title = "Bayesian Forest Plot: TriPM × Stress Reactivity",
    subtitle = "Moderation effects across acoustic outcomes and vowels"
  ) +
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_forest)

# Save plot
ggsave(
  "tripm_moderation_forest_plot.png",
  p_forest,
  width = 14,
  height = 12,
  dpi = 300
)

# ==============================================================================
# 6. EXTRACT CREDIBLE INTERACTIONS
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("CREDIBLE INTERACTIONS (95% CI excludes zero)\n")
cat(rep("=", 70), "\n", sep = "")

credible_interactions <- results_table_tripm %>%
  filter(Type == "Interaction" & Credible == TRUE) %>%
  arrange(Outcome, Vowel, Term)

if (nrow(credible_interactions) > 0) {
  print(data.frame(credible_interactions), row.names = FALSE)

  cat("\n=== Summary of Credible Interactions ===\n")
  cat("Total credible interactions:", nrow(credible_interactions), "\n")

  credible_interactions %>%
    count(Outcome, name = "N_credible") %>%
    arrange(desc(N_credible)) %>%
    print()
} else {
  cat("No credible interactions found (all 95% CIs include zero).\n")
}

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

# Save results table
write_csv(results_table_tripm, "tripm_moderation_results.csv")
write_csv(credible_interactions, "tripm_credible_interactions.csv")

# Save workspace
# save(
#   fitted_models_tripm,
#   results_table_tripm,
#   credible_interactions,
#   df_analysis_tripm,
#   file = "tripm_moderation_analysis.RData"
# )

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat("Results saved to:\n")
cat("  - tripm_moderation_results.csv\n")
cat("  - tripm_credible_interactions.csv\n")
cat("  - tripm_moderation_forest_plot.png\n")
cat("  - tripm_moderation_analysis.RData\n")
cat(rep("=", 70), "\n", sep = "")
