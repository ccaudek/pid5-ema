# ==============================================================================
# 06_moderation_with_variability_VERSION_A.R
# Test Mean + Variability as co-moderators of stress reactivity
# ==============================================================================
# PURPOSE:
#   Examine whether BOTH mean level (trait) AND variability (lability)
#   moderate vocal stress responses
#
# MODEL:
#   Voice = Mean + SD + Stress + Recovery +
#           Mean×Stress + SD×Stress +
#           Mean×Recovery + SD×Recovery +
#           Random effects
#
# INTERPRETATION:
#   - Mean×Stress: Does average personality level moderate stress reactivity?
#   - SD×Stress: Does emotional lability/instability moderate stress reactivity?
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayesplot)
  library(tidybayes)
  library(patchwork)
})

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS: MEAN + VARIABILITY (VERSION A)\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

if (!file.exists("results/variability/df_with_variability.rds")) {
  stop("Run 05_calculate_within_subject_variability.R first")
}

df <- readRDS("results/variability/df_with_variability.rds")

cat("Data loaded: N =", nrow(df), "observations\n")
cat("  N subjects:", n_distinct(df$ID), "\n\n")

# Filter to reliable variability estimates
df_complete <- df %>%
  filter(
    reliable_variability == TRUE,
    !is.na(pid5_negative_affectivity_mean_c),
    !is.na(pid5_negative_affectivity_sd_c)
  )

cat("After filtering for reliable variability:\n")
cat("  N observations:", nrow(df_complete), "\n")
cat("  N subjects:", n_distinct(df_complete$ID), "\n\n")

# ==============================================================================
# ADD CONTRAST CODES
# ==============================================================================

cat("Adding contrast codes for timepoint...\n")

df_complete <- df_complete %>%
  mutate(
    # Ensure timepoint is factor with correct levels
    timepoint = factor(
      timepoint,
      levels = c("baseline", "pre", "post")
    ),

    # Contrast 1: Stress (baseline vs pre-exam)
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre" ~ 0.5,
      timepoint == "post" ~ 0
    ),

    # Contrast 2: Recovery (pre-exam vs post-exam)
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre" ~ -0.5,
      timepoint == "post" ~ 0.5
    )
  )

cat("✓ Contrast codes added\n")
cat("  c1_stress: baseline (-0.5) vs pre (+0.5), post (0)\n")
cat("  c2_recovery: pre (-0.5) vs post (+0.5), baseline (0)\n\n")

# Verify contrast codes
cat("Verifying contrast codes:\n")
contrast_check <- df_complete %>%
  group_by(timepoint) %>%
  summarise(
    mean_c1 = mean(c1_stress),
    mean_c2 = mean(c2_recovery),
    n = n(),
    .groups = "drop"
  )
print(contrast_check)
cat("\n")

# ==============================================================================
# 2. MODEL SPECIFICATION
# ==============================================================================

cat("Model specification:\n")
cat("  Predictors: Mean (centered) + SD (centered)\n")
cat("  Moderators: Stress + Recovery\n")
cat("  Interactions: Mean×Stress, SD×Stress, Mean×Recovery, SD×Recovery\n\n")

# PID-5 domains
pid5_domains <- c(
  "negative_affectivity",
  "detachment",
  "antagonism",
  "disinhibition",
  "psychoticism"
)

# Function to build formula
build_formula_variability <- function(
  outcome,
  domain,
  family_spec = gaussian()
) {
  # Variable names
  mean_var <- paste0("pid5_", domain, "_mean_c")
  sd_var <- paste0("pid5_", domain, "_sd_c")

  # Formula with both mean and SD
  formula_str <- paste0(
    outcome,
    " ~ 1 + ",
    mean_var,
    " + ",
    sd_var,
    " + ",
    "c1_stress + c2_recovery + ",
    mean_var,
    ":c1_stress + ",
    sd_var,
    ":c1_stress + ",
    mean_var,
    ":c2_recovery + ",
    sd_var,
    ":c2_recovery + ",
    "(1 + c1_stress + c2_recovery || ID)"
  )

  bf(formula_str, family = family_spec)
}

# ==============================================================================
# 3. PRIORS
# ==============================================================================

# Weakly informative priors
get_priors_variability <- function(outcome_var, is_lognormal = FALSE) {
  y_mu <- mean(df_complete[[outcome_var]], na.rm = TRUE)
  y_sd <- sd(df_complete[[outcome_var]], na.rm = TRUE)

  if (is_lognormal) {
    b_prior <- "normal(0, 1)"
    intercept_location <- log(max(y_mu, 1e-06))
    intercept_scale <- 1
  } else {
    b_prior <- "normal(0, 10)"
    intercept_location <- y_mu
    intercept_scale <- 2.5 * y_sd
  }

  priors <- c(
    prior_string(b_prior, class = "b"),
    prior_string(
      paste0("student_t(3, ", intercept_location, ", ", intercept_scale, ")"),
      class = "Intercept"
    ),
    prior_string("exponential(1)", class = "sigma")
  )

  return(priors)
}

# ==============================================================================
# 4. FIT MODELS FOR EACH DOMAIN AND OUTCOME
# ==============================================================================

cat(rep("=", 70), "\n")
cat("FITTING MODELS\n")
cat(rep("=", 70), "\n\n")

dir.create(
  "models/variability_version_a",
  showWarnings = FALSE,
  recursive = TRUE
)

# Define outcomes
outcomes <- list(
  # F0 (pitch)
  list(
    var = "f0_mean_a",
    name = "F0 mean /a/",
    family = gaussian(),
    log = FALSE
  ),
  list(var = "f0_std_a", name = "F0 SD /a/", family = lognormal(), log = TRUE),

  # F2 (formant)
  list(
    var = "f2_mean_a",
    name = "F2 mean /a/",
    family = gaussian(),
    log = FALSE
  ),
  list(var = "f2_std_a", name = "F2 SD /a/", family = lognormal(), log = TRUE),

  # Voice quality
  list(var = "jitter_a", name = "Jitter /a/", family = lognormal(), log = TRUE),
  list(var = "nne_a", name = "NNE /a/", family = gaussian(), log = FALSE)
)

# Storage for results
all_results <- list()

# Fit models
for (outcome_info in outcomes) {
  outcome_var <- outcome_info$var
  outcome_name <- outcome_info$name
  outcome_family <- outcome_info$family
  is_log <- outcome_info$log

  cat("\n", rep("=", 70), "\n", sep = "")
  cat("OUTCOME:", outcome_name, "\n")
  cat(rep("=", 70), "\n\n")

  domain_results <- list()

  for (domain in pid5_domains) {
    domain_name <- str_to_title(str_replace_all(domain, "_", " "))
    cat("\n--- Domain:", domain_name, "---\n")

    # Build formula
    formula_obj <- build_formula_variability(
      outcome_var,
      domain,
      outcome_family
    )

    cat("Formula:\n")
    print(formula_obj)
    cat("\n")

    # Get priors
    priors <- get_priors_variability(outcome_var, is_log)

    # Filename for caching
    model_file <- paste0(
      "models/variability_version_a/",
      outcome_var,
      "_",
      domain
    )

    # Fit model
    cat("Fitting model...\n")

    fit <- brm(
      formula = formula_obj,
      data = df_complete,
      prior = priors,
      chains = 4,
      iter = 5000,
      warmup = 2500,
      cores = 4,
      seed = 123,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = model_file,
      backend = "cmdstanr"
    )

    cat("✓ Model fitted\n")

    # Extract key parameters
    posterior <- as_draws_df(fit)

    # Variable names
    mean_var <- paste0("pid5_", domain, "_mean_c")
    sd_var <- paste0("pid5_", domain, "_sd_c")

    # Parameter names in model
    mean_main <- paste0("b_", mean_var)
    sd_main <- paste0("b_", sd_var)
    mean_stress <- paste0("b_", mean_var, ":c1_stress")
    sd_stress <- paste0("b_", sd_var, ":c1_stress")
    mean_recovery <- paste0("b_", mean_var, ":c2_recovery")
    sd_recovery <- paste0("b_", sd_var, ":c2_recovery")

    # Extract estimates
    extract_estimate <- function(param_name) {
      if (param_name %in% names(posterior)) {
        tibble(
          mean = mean(posterior[[param_name]], na.rm = TRUE),
          sd = sd(posterior[[param_name]], na.rm = TRUE),
          ci_lower = quantile(posterior[[param_name]], 0.025, na.rm = TRUE),
          ci_upper = quantile(posterior[[param_name]], 0.975, na.rm = TRUE),
          prob_direction = max(
            mean(posterior[[param_name]] > 0, na.rm = TRUE),
            mean(posterior[[param_name]] < 0, na.rm = TRUE)
          ),
          credible = sign(ci_lower) == sign(ci_upper)
        )
      } else {
        tibble(
          mean = NA,
          sd = NA,
          ci_lower = NA,
          ci_upper = NA,
          prob_direction = NA,
          credible = FALSE
        )
      }
    }

    result <- tibble(
      outcome = outcome_name,
      domain = domain_name,

      # Mean effects
      mean_main_est = extract_estimate(mean_main)$mean,
      mean_main_ci_lower = extract_estimate(mean_main)$ci_lower,
      mean_main_ci_upper = extract_estimate(mean_main)$ci_upper,
      mean_main_credible = extract_estimate(mean_main)$credible,

      mean_stress_est = extract_estimate(mean_stress)$mean,
      mean_stress_ci_lower = extract_estimate(mean_stress)$ci_lower,
      mean_stress_ci_upper = extract_estimate(mean_stress)$ci_upper,
      mean_stress_credible = extract_estimate(mean_stress)$credible,

      mean_recovery_est = extract_estimate(mean_recovery)$mean,
      mean_recovery_ci_lower = extract_estimate(mean_recovery)$ci_lower,
      mean_recovery_ci_upper = extract_estimate(mean_recovery)$ci_upper,
      mean_recovery_credible = extract_estimate(mean_recovery)$credible,

      # SD (variability) effects
      sd_main_est = extract_estimate(sd_main)$mean,
      sd_main_ci_lower = extract_estimate(sd_main)$ci_lower,
      sd_main_ci_upper = extract_estimate(sd_main)$ci_upper,
      sd_main_credible = extract_estimate(sd_main)$credible,

      sd_stress_est = extract_estimate(sd_stress)$mean,
      sd_stress_ci_lower = extract_estimate(sd_stress)$ci_lower,
      sd_stress_ci_upper = extract_estimate(sd_stress)$ci_upper,
      sd_stress_credible = extract_estimate(sd_stress)$credible,

      sd_recovery_est = extract_estimate(sd_recovery)$mean,
      sd_recovery_ci_lower = extract_estimate(sd_recovery)$ci_lower,
      sd_recovery_ci_upper = extract_estimate(sd_recovery)$ci_upper,
      sd_recovery_credible = extract_estimate(sd_recovery)$credible
    )

    domain_results[[domain]] <- result
  }

  # Combine results for this outcome
  all_results[[outcome_var]] <- bind_rows(domain_results)
}

# Combine all results
results_combined <- bind_rows(all_results)

# ==============================================================================
# 5. SAVE RESULTS
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("SAVING RESULTS\n")
cat(rep("=", 70), "\n\n")

dir.create("results/variability", showWarnings = FALSE, recursive = TRUE)

saveRDS(results_combined, "results/variability/version_a_results.rds")
rio::export(results_combined, "results/variability/version_a_results.csv")

cat("✓ Saved: version_a_results.rds\n")
cat("✓ Saved: version_a_results.csv\n\n")

# ==============================================================================
# 6. SUMMARY: CREDIBLE VARIABILITY EFFECTS
# ==============================================================================

cat(rep("=", 70), "\n")
cat("SUMMARY: CREDIBLE VARIABILITY EFFECTS\n")
cat(rep("=", 70), "\n\n")

# Count credible effects
n_mean_main <- sum(results_combined$mean_main_credible, na.rm = TRUE)
n_mean_stress <- sum(results_combined$mean_stress_credible, na.rm = TRUE)
n_mean_recovery <- sum(results_combined$mean_recovery_credible, na.rm = TRUE)

n_sd_main <- sum(results_combined$sd_main_credible, na.rm = TRUE)
n_sd_stress <- sum(results_combined$sd_stress_credible, na.rm = TRUE)
n_sd_recovery <- sum(results_combined$sd_recovery_credible, na.rm = TRUE)

cat("MEAN (Trait) Effects:\n")
cat("  Main effects:", n_mean_main, "\n")
cat("  Stress interactions:", n_mean_stress, "\n")
cat("  Recovery interactions:", n_mean_recovery, "\n")
cat("  TOTAL:", n_mean_main + n_mean_stress + n_mean_recovery, "\n\n")

cat("SD (Variability/Lability) Effects:\n")
cat("  Main effects:", n_sd_main, "\n")
cat("  Stress interactions:", n_sd_stress, "\n")
cat("  Recovery interactions:", n_sd_recovery, "\n")
cat("  TOTAL:", n_sd_main + n_sd_stress + n_sd_recovery, "\n\n")

# Show credible SD effects
sd_credible <- results_combined %>%
  filter(sd_stress_credible | sd_recovery_credible) %>%
  select(outcome, domain, starts_with("sd_"))

if (nrow(sd_credible) > 0) {
  cat("CREDIBLE VARIABILITY × STRESS/RECOVERY INTERACTIONS:\n")
  cat(
    "(These show that emotional lability/instability moderates stress reactivity)\n\n"
  )
  print(sd_credible, n = Inf)
  cat("\n")
} else {
  cat("No credible variability × stress/recovery interactions found.\n\n")
}

# ==============================================================================
# 7. COMPARISON: MEAN VS VARIABILITY EFFECTS
# ==============================================================================

cat(rep("=", 70), "\n")
cat("KEY QUESTION: Does Variability Add Beyond Mean?\n")
cat(rep("=", 70), "\n\n")

# For each outcome-domain combination, compare mean vs SD effects
comparison <- results_combined %>%
  mutate(
    # Which has stronger stress interaction?
    stronger_stress = case_when(
      !mean_stress_credible & !sd_stress_credible ~ "Neither",
      mean_stress_credible & !sd_stress_credible ~ "Mean only",
      !mean_stress_credible & sd_stress_credible ~ "SD only",
      abs(mean_stress_est) > abs(sd_stress_est) ~ "Mean stronger",
      TRUE ~ "SD stronger"
    ),

    # Which has stronger recovery interaction?
    stronger_recovery = case_when(
      !mean_recovery_credible & !sd_recovery_credible ~ "Neither",
      mean_recovery_credible & !sd_recovery_credible ~ "Mean only",
      !mean_recovery_credible & sd_recovery_credible ~ "SD only",
      abs(mean_recovery_est) > abs(sd_recovery_est) ~ "Mean stronger",
      TRUE ~ "SD stronger"
    )
  ) %>%
  select(outcome, domain, stronger_stress, stronger_recovery)

cat("Which predictor has stronger moderation effects?\n\n")
print(comparison, n = Inf)
cat("\n")

# Summary counts
cat("Summary across all outcomes:\n")
cat("\nStress interactions:\n")
print(table(comparison$stronger_stress))
cat("\nRecovery interactions:\n")
print(table(comparison$stronger_recovery))
cat("\n")

# ==============================================================================
# 8. INTERPRETATION GUIDE
# ==============================================================================

cat(rep("=", 70), "\n")
cat("INTERPRETATION GUIDE\n")
cat(rep("=", 70), "\n\n")

cat("MEAN effects (trait level):\n")
cat("  - Main: Does average personality predict voice?\n")
cat("  - × Stress: Do people high in trait X react differently to stress?\n")
cat("  - × Recovery: Do people high in trait X recover differently?\n\n")

cat("SD effects (variability/lability):\n")
cat("  - Main: Does emotional instability predict voice?\n")
cat("  - × Stress: Do emotionally labile people react differently to stress?\n")
cat("  - × Recovery: Do emotionally labile people recover differently?\n\n")

cat("Clinical relevance:\n")
cat("  If SD × Stress is credible → Emotional instability (not just level)\n")
cat("  predicts stress vulnerability. This suggests:\n")
cat("  - Dysregulation matters beyond mean level\n")
cat("  - Intra-individual variability is a risk factor\n")
cat("  - EMA design captures clinically relevant dynamics\n\n")

cat(rep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Next steps:\n")
cat("1. Review version_a_results.csv for detailed estimates\n")
cat("2. Run 07_visualize_variability_effects.R to plot key findings\n")
cat("3. Compare with Version B (SD-only) if needed\n")
cat("4. Consider which findings to report in manuscript\n\n")
