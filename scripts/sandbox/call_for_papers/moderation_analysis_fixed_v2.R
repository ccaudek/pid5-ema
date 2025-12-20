# ==============================================================================
# Bayesian Moderation Analysis: PID-5 Traits × Stress Reactivity
# Testing context-dependent expression of personality pathology
# CORRECTED VERSION - Automatic prior calculation from data
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
cat("MODERATION ANALYSIS: PID-5 × Stress Reactivity\n")
cat("All 5 PID-5 Domains\n")
cat(rep("=", 70), "\n", sep = "")

# ==============================================================================
# 1. VERIFY DATA EXISTS
# ==============================================================================

if (!exists("df_analysis")) {
  if (file.exists("results/df_analysis.rds")) {
    df_analysis <- readRDS("results/df_analysis.rds")
    cat("✓ Loaded df_analysis from results/df_analysis.rds\n")
  } else if (file.exists("results/main_effects_workspace.RData")) {
    load("results/main_effects_workspace.RData")
    cat("✓ Loaded df_analysis from workspace\n")
  } else {
    stop("df_analysis not found. Run voice_personality_analysis.R first.")
  }
}

# ==============================================================================
# 2. PREPARE DATA FOR MODERATION MODELS
# ==============================================================================

# Check PID-5 distributions at baseline
cat("\n=== PID-5 Baseline Trait Distributions (All 5 Domains) ===\n")
df_analysis %>%
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

# ==============================================================================
# 3. HANDLE ZEROS IN VARIABLES TO BE LOG-TRANSFORMED
# ==============================================================================

vars_to_check <- c(
  "f2_std_a",
  "f2_std_i",
  "f2_std_u",
  "f0_std_a",
  "f0_std_i",
  "f0_std_u",
  "jitter_a",
  "jitter_i",
  "jitter_u"
)

cat("\n=== Checking for zeros/negatives before log transform ===\n")
for (v in vars_to_check) {
  if (v %in% names(df_analysis)) {
    n_zero <- sum(df_analysis[[v]] == 0, na.rm = TRUE)
    n_neg <- sum(df_analysis[[v]] < 0, na.rm = TRUE)
    n_inf <- sum(is.infinite(df_analysis[[v]]), na.rm = TRUE)
    if (n_zero > 0 | n_neg > 0 | n_inf > 0) {
      cat(
        v,
        ": zeros=",
        n_zero,
        ", negatives=",
        n_neg,
        ", infinite=",
        n_inf,
        "\n"
      )
      min_pos <- min(
        df_analysis[[v]][df_analysis[[v]] > 0 & is.finite(df_analysis[[v]])],
        na.rm = TRUE
      )
      df_analysis[[v]] <- ifelse(
        df_analysis[[v]] <= 0 | !is.finite(df_analysis[[v]]),
        min_pos / 2,
        df_analysis[[v]]
      )
      cat("  -> Replaced with", min_pos / 2, "\n")
    }
  }
}

# ==============================================================================
# 4. CREATE CONTRAST CODES
# ==============================================================================

df_analysis <- df_analysis %>%
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
# 5. LOG-TRANSFORM SKEWED VARIABLES
# ==============================================================================

# Check if already transformed (logged values would be small/negative)
mean_jitter <- mean(df_analysis$jitter_a, na.rm = TRUE)
if (mean_jitter > 0.5) {
  cat("\n=== Applying log transformation ===\n")
  df_analysis <- df_analysis %>%
    mutate(across(
      matches("^(f0_std|f2_std|jitter)_[aiu]$"),
      ~ log(.x)
    ))
  cat("✓ Log-transformed f0_std, f2_std, and jitter\n")
} else {
  cat("\n(Variables appear already log-transformed)\n")
}

# Final check for infinite values
n_inf <- sum(sapply(df_analysis, function(x) sum(is.infinite(x))))
if (n_inf > 0) {
  cat("\n⚠ Warning:", n_inf, "infinite values found. Replacing with NA...\n")
  df_analysis <- df_analysis %>%
    mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), NA, .)))
}

# ==============================================================================
# 6. COMPUTE EMPIRICAL STATISTICS FOR PRIORS
# ==============================================================================

cat("\n=== Computing empirical statistics for priors ===\n")

# Define outcomes and vowels
outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
vowels <- c("a", "i", "u")

# Which outcomes are log-transformed
log_transformed <- c("f0_std", "f2_std", "jitter")

# Compute statistics for each outcome-vowel combination
empirical_stats <- expand_grid(outcome = outcomes, vowel = vowels) %>%
  rowwise() %>%
  mutate(
    colname = paste0(outcome, "_", vowel),
    # Get the column values
    values = list(df_analysis[[colname]]),
    # Compute mean and SD on observed scale
    obs_mean = mean(values, na.rm = TRUE),
    obs_sd = sd(values, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-values)

cat("✓ Computed empirical statistics:\n")
print(empirical_stats)

# ==============================================================================
# 7. DEFINE PRIOR-MAKING FUNCTION
# ==============================================================================

# Prior scale parameters (can be adjusted based on domain knowledge)
prior_scales <- list(
  f0_mean = list(intercept_scale = 30), # Hz scale

  f0_std = list(intercept_scale = 0.5), # log scale
  jitter = list(intercept_scale = 0.5), # log scale
  nne = list(intercept_scale = 5), # dB scale
  f2_mean = list(intercept_scale = 150), # Hz scale
  f2_std = list(intercept_scale = 0.6) # log scale
)

make_priors_vowel <- function(outcome, vowel, stats_df = empirical_stats) {
  # Get empirical statistics for this outcome-vowel
  row <- stats_df %>%
    filter(outcome == !!outcome, vowel == !!vowel)

  if (nrow(row) != 1) {
    stop("Could not find statistics for outcome=", outcome, ", vowel=", vowel)
  }

  # Get the observed mean (already on correct scale after log-transform)
  center <- row$obs_mean
  scale <- prior_scales[[outcome]]$intercept_scale

  # Format for prior string
  fmt <- function(x, digits = 4) {
    formatC(x, format = "f", digits = digits)
  }

  # Intercept prior centered on empirical mean
  intercept_prior <- prior_string(
    paste0("student_t(3, ", fmt(center), ", ", fmt(scale), ")"),
    class = "Intercept"
  )

  # Fixed effects: weakly informative
  # Main effects and interactions get the same prior
  # (brms doesn't support pattern matching for coef names)
  fixed_prior <- prior_string("normal(0, 1)", class = "b")

  # Random effects and residual
  random_priors <- c(
    prior_string("exponential(2)", class = "sd"),
    prior_string("exponential(2)", class = "sigma")
  )

  # Combine all priors
  c(intercept_prior, fixed_prior, random_priors)
}

# Test the function
cat("\n=== Testing prior function ===\n")
cat("Priors for jitter_a:\n")
print(make_priors_vowel("jitter", "a"))
cat("\nPriors for f0_mean_i:\n")
print(make_priors_vowel("f0_mean", "i"))

# ==============================================================================
# 8. MODEL CONFIGURATION
# ==============================================================================

# Map outcome -> family
outcome_families <- list(
  f0_mean = gaussian(),
  f0_std = gaussian(),
  jitter = gaussian(),
  nne = gaussian(),
  f2_mean = asym_laplace(),
  f2_std = gaussian()
)

# PID-5 domain names
pid5_domains <- c(
  "pid5_negative_affectivity_bl_c",
  "pid5_detachment_bl_c",
  "pid5_antagonism_bl_c",
  "pid5_disinhibition_bl_c",
  "pid5_psychoticism_bl_c"
)

# Sampling options
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
seed <- 123
control <- list(adapt_delta = 0.99, max_treedepth = 15)

# Create models directory
dir.create("models", showWarnings = FALSE)

# Flag: if TRUE run models; if FALSE just show plan
run_models <- TRUE

# ==============================================================================
# 9. FIT MODELS
# ==============================================================================

fitted_models <- list()

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("FITTING MODERATION MODELS\n")
cat(rep("=", 70), "\n", sep = "")

for (v in vowels) {
  message("\n--- Vowel: /", v, "/ -------------------------------")
  for (out in outcomes) {
    colname <- paste0(out, "_", v)
    fam <- outcome_families[[out]]

    # Check column exists
    if (!colname %in% names(df_analysis)) {
      message("  Skipping ", colname, " (column not found)")
      next
    }

    # Build formula with all 5 PID-5 domains
    fmla <- bf(
      as.formula(
        paste0(
          colname,
          " ~ c1_stress * (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c + ",
          "pid5_antagonism_bl_c + pid5_disinhibition_bl_c + pid5_psychoticism_bl_c) + ",
          "c2_recovery * (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c + ",
          "pid5_antagonism_bl_c + pid5_disinhibition_bl_c + pid5_psychoticism_bl_c) + ",
          "(1 + c1_stress + c2_recovery | ID)"
        )
      )
    )

    # Generate priors automatically
    priors_here <- make_priors_vowel(out, v)

    model_name <- paste0("m_", out, "_", v)

    # Get family name for display
    fam_name <- if (!is.null(fam$family)) {
      fam$family
    } else if (!is.null(attr(fam, "family"))) {
      attr(fam, "family")
    } else {
      class(fam)[1]
    }

    message("  Model: ", model_name, "  (family: ", fam_name, ")")

    if (run_models) {
      tryCatch(
        {
          fit <- brm(
            formula = fmla,
            data = df_analysis,
            family = fam,
            prior = priors_here,
            iter = iter,
            warmup = warmup,
            chains = chains,
            cores = cores,
            seed = seed,
            control = control,
            file = file.path("models", model_name)
          )
          fitted_models[[model_name]] <- fit

          # Quick diagnostic check
          print(summary(fit, waic = FALSE))
        },
        error = function(e) {
          message("  ERROR: ", e$message)
        }
      )
    } else {
      message("    (run_models = FALSE) - priors:")
      print(priors_here)
    }
  }
}

cat("\n✓ Fitted", length(fitted_models), "models\n")

# ==============================================================================
# 10. EXTRACT AND SUMMARIZE RESULTS
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("EXTRACTING RESULTS\n")
cat(rep("=", 70), "\n", sep = "")

results_table <- lapply(names(fitted_models), function(mn) {
  fit <- fitted_models[[mn]]
  fe <- fixef(fit, summary = TRUE)

  # Parse model name: m_f0_mean_a -> outcome = f0_mean, vowel = a
  parts <- str_match(mn, "m_(.+)_([aiu])$")
  outcome <- parts[2]
  vowel <- parts[3]

  tibble(
    outcome = outcome,
    vowel = vowel,
    parameter = rownames(fe),
    estimate = fe[, "Estimate"],
    est_error = fe[, "Est.Error"],
    ci_lower = fe[, "Q2.5"],
    ci_upper = fe[, "Q97.5"],
    significant = (ci_lower > 0) | (ci_upper < 0)
  ) %>%
    mutate(
      type = case_when(
        str_detect(parameter, ":") ~ "Interaction",
        TRUE ~ "Main"
      )
    )
}) %>%
  bind_rows()

# Display summary
cat("\n=== RESULTS SUMMARY ===\n")
cat("Total effects:", nrow(results_table), "\n")
cat("Significant (95% CI excludes 0):", sum(results_table$significant), "\n")

# Show significant effects
cat("\n=== SIGNIFICANT EFFECTS ===\n\n")
results_table %>%
  filter(significant) %>%
  arrange(outcome, vowel, type, parameter) %>%
  print(n = 100)

# Show significant interactions specifically
cat("\n=== SIGNIFICANT INTERACTIONS ===\n\n")
results_table %>%
  filter(significant, type == "Interaction") %>%
  arrange(outcome, vowel, parameter) %>%
  print(n = 50)

# ==============================================================================
# 11. SAVE RESULTS
# ==============================================================================

dir.create("results", showWarnings = FALSE)

# Save results table
write_csv(results_table, "results/moderation_effects_all.csv")

# Save empirical statistics used for priors
write_csv(empirical_stats, "results/empirical_stats_for_priors.csv")

# Save workspace with fitted models
save(
  fitted_models,
  df_analysis,
  results_table,
  empirical_stats,
  file = "results/moderation_analysis_complete.RData"
)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nFiles saved:\n")
cat("  - results/moderation_effects_all.csv\n")
cat("  - results/empirical_stats_for_priors.csv\n")
cat("  - results/moderation_analysis_complete.RData\n")
cat("\nTo reload: load('results/moderation_analysis_complete.RData')\n\n")
