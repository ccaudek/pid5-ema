# ==============================================================================
# Bayesian Moderation Analysis: PID-5 Traits × Stress Reactivity
# Testing context-dependent expression of personality pathology
# CORRECTED VERSION - Fixed prior specification
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

# Fix zeros in f2_std variables (and any other variables that need log transform)
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
      # Replace with half the minimum positive value
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

# Create contrast-coded timepoint variable for interaction interpretation
# C1: PRE vs BASELINE (stress reactivity)
# C2: POST vs PRE (stress recovery)

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
  # Not yet transformed
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
# 6. DEFINE PRIORS (CORRECTED)
# ==============================================================================

base_stats <- df_analysis %>%
  select(ID, matches("^(f0_mean|f0_std|jitter|nne|f2_mean|f2_std)_[aiu]$")) %>%
  pivot_longer(-ID) %>%
  separate(name, into = c("outcome", "vowel"), sep = "_(?=[aiu]$)") %>%
  group_by(outcome, vowel) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

log_stats <- df_analysis %>%
  select(ID, matches("^(f0_std|f2_std|jitter)_[aiu]$")) %>%
  pivot_longer(-ID) %>%
  separate(name, into = c("outcome", "vowel"), sep = "_(?=[aiu]$)") %>%
  group_by(outcome, vowel) %>%
  summarise(
    log_mean = mean(log(value[value > 0]), na.rm = TRUE),
    log_sd = sd(log(value[value > 0]), na.rm = TRUE),
    .groups = "drop"
  )

vowel_means <- base_stats %>%
  left_join(log_stats, by = c("outcome", "vowel"))

make_priors_vowel <- function(outcome, vowel) {
  row <- vowel_means %>%
    filter(outcome == !!outcome, vowel == !!vowel)

  stopifnot(nrow(row) == 1)

  m <- row$mean
  lm <- row$log_mean

  intercept_prior <- switch(
    outcome,

    # NON log-transformati → usa mean su scala originale
    "f0_mean" = prior_string(
      paste0("student_t(3, ", m, ", 30)"),
      class = "Intercept"
    ),

    "nne" = prior_string(
      paste0("student_t(3, ", m, ", 5)"),
      class = "Intercept"
    ),

    "f2_mean" = prior_string(
      paste0("student_t(3, ", m, ", 150)"),
      class = "Intercept"
    ),

    # log-transformati → usa log_mean
    "f0_std" = prior_string(
      paste0("student_t(3, ", lm, ", 0.5)"),
      class = "Intercept"
    ),

    "f2_std" = prior_string(
      paste0("student_t(3, ", lm, ", 0.6)"),
      class = "Intercept"
    ),

    "jitter" = prior_string(
      paste0("student_t(3, ", lm, ", 0.5)"),
      class = "Intercept"
    ),

    stop("Outcome non riconosciuto.")
  )

  priors <- c(
    intercept_prior,
    prior(normal(0, 1), class = "b"),
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sigma")
  )

  return(priors)
}

make_priors_vowel("jitter", "a")
make_priors_vowel("f0_mean", "i")


# --- Empirical parameters from summary statistics
# vowel_means <- list(
#   a = list(
#     f0_mean = 189.9,
#     f0_std = 6.230,
#     jitter = 0.8366,
#     nne = -24.42,
#     f2_mean = 1207.4,
#     f2_std = 59.37
#   ),
#   i = list(
#     f0_mean = 194.4,
#     f0_std = 8.888,
#     jitter = 1.058,
#     nne = -27.89,
#     f2_mean = 2189.0,
#     f2_std = 141.04
#   ),
#   u = list(
#     f0_mean = 195.1,
#     f0_std = 9.724,
#     jitter = 1.286,
#     nne = -28.429,
#     f2_mean = 1123.2,
#     f2_std = 161.57
#   )
# )

# PID-5 domain names (for generating interaction coefficient names)
pid5_domains <- c(
  "pid5_negative_affectivity_bl_c",
  "pid5_detachment_bl_c",
  "pid5_antagonism_bl_c",
  "pid5_disinhibition_bl_c",
  "pid5_psychoticism_bl_c"
)

make_priors_vowel <- function(outcome, vowel) {
  vm <- vowel_means[[vowel]]

  fmt <- function(x, digits = 6) {
    formatC(x, format = "f", digits = digits)
  }

  # -------------------------------
  # Intercept priors (outcome-scale)
  # -------------------------------
  intercept_prior <- switch(
    outcome,

    "f0_mean" = prior_string(
      paste0("student_t(3, ", fmt(vm$f0_mean), ", 30)"),
      class = "Intercept"
    ),

    "f0_std" = prior_string(
      paste0("student_t(3, ", fmt(log(vm$f0_std)), ", 0.5)"),
      class = "Intercept"
    ),

    "jitter" = prior_string(
      paste0("student_t(3, ", fmt(log(vm$jitter)), ", 0.5)"),
      class = "Intercept"
    ),

    "nne" = prior_string(
      paste0("student_t(3, ", fmt(vm$nne), ", 5)"),
      class = "Intercept"
    ),

    "f2_mean" = prior_string(
      paste0("student_t(3, ", fmt(vm$f2_mean), ", 150)"),
      class = "Intercept"
    ),

    "f2_std" = prior_string(
      paste0("student_t(3, ", fmt(log(vm$f2_std)), ", 0.6)"),
      class = "Intercept"
    ),

    stop("Outcome non riconosciuto")
  )

  # -----------------------------------
  # Fixed effects priors
  # -----------------------------------
  # CORRECTED: Use a single prior for all b coefficients
  # The coef = "pattern:" syntax does NOT work in brms
  #
  # Option 1: Same prior for all fixed effects (simple, recommended)
  # Option 2: List each interaction coefficient explicitly (tedious)
  #
  # Using Option 1 with a moderately informative prior

  priors_fixed <- prior_string("normal(0, 1)", class = "b")

  # If you want tighter priors on interactions specifically, you would need to
  # list each coefficient explicitly like this (uncomment if desired):
  #
  # interaction_priors <- c(
  #   # c1_stress interactions
  #   lapply(pid5_domains, function(d) {
  #     prior_string("normal(0, 0.5)", class = "b",
  #                  coef = paste0("c1_stress:", d))
  #   }),
  #   # c2_recovery interactions
  #   lapply(pid5_domains, function(d) {
  #     prior_string("normal(0, 0.5)", class = "b",
  #                  coef = paste0(d, ":c2_recovery"))
  #   })
  # ) %>% unlist()

  # -------------------------
  # Random effects + residual
  # -------------------------
  priors_random <- c(
    prior_string("exponential(2)", class = "sd"),
    prior_string("exponential(2)", class = "sigma")
  )

  # -------------------------
  # Final prior set
  # -------------------------
  c(
    intercept_prior,
    priors_fixed,
    priors_random
  )
}


# --- Map outcome -> family
outcome_families <- list(
  f0_mean = gaussian(),
  f0_std = gaussian(),
  jitter = gaussian(),
  nne = gaussian(),
  f2_mean = asym_laplace(),
  f2_std = gaussian()
)

# --- Outcomes and vowels
outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
vowels <- c("a", "i", "u")

# --- Sampling options
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
seed <- 123
control <- list(adapt_delta = 0.99, max_treedepth = 15)

# --- Create models directory
dir.create("models", showWarnings = FALSE)

# --- Flag: if TRUE run models; if FALSE just show plan
run_models <- TRUE

# ==============================================================================
# 7. FIT MODELS
# ==============================================================================

fitted_models <- list()

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("FITTING MODERATION MODELS\n")
cat(rep("=", 70), "\n", sep = "")

for (v in vowels) {
  message("\n--- Vowel: /", v, "/ -------------------------------")
  for (out in outcomes) {
    # Column name (e.g., f0_mean_a)
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

    priors_here <- make_priors_vowel(out, v)

    model_name <- paste0("m_", out, "_", v)

    # Get family name for display
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
# 8. EXTRACT AND SUMMARIZE RESULTS
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("EXTRACTING RESULTS\n")
cat(rep("=", 70), "\n", sep = "")

results_table <- lapply(names(fitted_models), function(mn) {
  fit <- fitted_models[[mn]]
  fe <- fixef(fit, summary = TRUE)

  # Parse model name: m_f0_mean_a -> outcome = f0_mean, vowel = a
  # Handle both m_f0_mean_a and m_f0_std_a patterns
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
# 9. SAVE RESULTS
# ==============================================================================

dir.create("results", showWarnings = FALSE)

# Save results table
write_csv(results_table, "results/moderation_effects_all.csv")

# Save workspace with fitted models
save(
  fitted_models,
  df_analysis,
  results_table,
  file = "results/moderation_analysis_complete.RData"
)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nFiles saved:\n")
cat("  - results/moderation_effects_all.csv\n")
cat("  - results/moderation_analysis_complete.RData\n")
cat("\nTo reload: load('results/moderation_analysis_complete.RData')\n\n")
