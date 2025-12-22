# ==============================================================================
# 12_fit_between_within_models.R
# Fit Bayesian models with between-person and within-person decomposition
# ==============================================================================
# PURPOSE:
#   Test whether within-person state fluctuations in PID-5 predict voice
#   beyond stable between-person traits
#
# MODEL:
#   Voice[t] = β0 +
#              β_between × Trait[constant] +
#              β_within × State[varying] +
#              β_stress × Stress +
#              β_recovery × Recovery +
#              β_between:stress × Trait × Stress +
#              β_within:stress × State × Stress +
#              β_between:recovery × Trait × Recovery +
#              β_within:recovery × State × Recovery +
#              Random effects
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayesplot)
  library(tidybayes)
})

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

cat("\n", rep("=", 70), "\n", sep = "")
cat("BETWEEN-WITHIN MODERATION ANALYSIS\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

if (!file.exists("results/between_within/df_between_within.rds")) {
  stop("Run 11_prepare_between_within_data.R first")
}

df <- readRDS("results/between_within/df_between_within.rds")

cat(
  "Data loaded: N =",
  nrow(df),
  "observations,",
  n_distinct(df$ID),
  "subjects\n\n"
)

# Complete cases only
df_complete <- df %>%
  filter(
    !is.na(pid5_negative_affectivity_between_c),
    !is.na(pid5_negative_affectivity_within_c)
  )

cat("Complete cases: N =", nrow(df_complete), "\n\n")

# ==============================================================================
# 2. MODEL SPECIFICATION
# ==============================================================================

# We'll fit models for F0 mean (vowel /a/) as main example
# Then extend to other outcomes if results are promising

# Predictors
pid5_domains <- c(
  "negative_affectivity",
  "detachment",
  "antagonism",
  "disinhibition",
  "psychoticism"
)

# Build formula
# Fixed effects: between + within + stress + recovery + all 2-way interactions
build_formula <- function(outcome, family_spec = gaussian()) {
  # Between-person main effects
  between_terms <- paste0("pid5_", pid5_domains, "_between_c")
  between_formula <- paste(between_terms, collapse = " + ")

  # Within-person main effects
  within_terms <- paste0("pid5_", pid5_domains, "_within_c")
  within_formula <- paste(within_terms, collapse = " + ")

  # Between × Stress interactions
  between_stress <- paste0("pid5_", pid5_domains, "_between_c:c1_stress")
  between_stress_formula <- paste(between_stress, collapse = " + ")

  # Within × Stress interactions
  within_stress <- paste0("pid5_", pid5_domains, "_within_c:c1_stress")
  within_stress_formula <- paste(within_stress, collapse = " + ")

  # Between × Recovery interactions
  between_recovery <- paste0("pid5_", pid5_domains, "_between_c:c2_recovery")
  between_recovery_formula <- paste(between_recovery, collapse = " + ")

  # Within × Recovery interactions
  within_recovery <- paste0("pid5_", pid5_domains, "_within_c:c2_recovery")
  within_recovery_formula <- paste(within_recovery, collapse = " + ")

  # Complete formula
  formula_str <- paste0(
    outcome,
    " ~ 1 + c1_stress + c2_recovery + ",
    between_formula,
    " + ",
    within_formula,
    " + ",
    between_stress_formula,
    " + ",
    within_stress_formula,
    " + ",
    between_recovery_formula,
    " + ",
    within_recovery_formula,
    " + (1 + c1_stress + c2_recovery || ID)"
  )

  bf(formula_str, family = family_spec)
}

# ==============================================================================
# 3. PRIORS
# ==============================================================================

# Weakly informative priors
# More conservative for within-person effects (less data per person)

get_priors_bw <- function(
  outcome,
  is_lognormal = FALSE,
  has_cor = FALSE,
  data = df_analysis
) {
  # prendi la y giusta (se lognormal userai outcome_pos fuori, qui teniamo outcome)
  y <- data[[outcome]]
  y <- as.numeric(y)
  y_ok <- y[is.finite(y)]
  y_mu <- mean(y_ok, na.rm = TRUE)
  y_sd <- sd(y_ok, na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
  if (!is.finite(y_mu)) y_mu <- 0

  # prior "puliti"
  pri <- c(
    prior(normal(0, if (is_lognormal) 1 else 10), class = "b"),
    prior(
      student_t(
        3,
        if (is_lognormal) log(max(y_mu, 1e-6)) else y_mu,
        if (is_lognormal) 1 else 2.5 * y_sd
      ),
      class = "Intercept"
    ),
    prior(exponential(1), class = "sigma")
  )

  # LKJ solo se stimi correlazioni dei random effects (|), NON con ||
  if (has_cor) {
    pri <- c(pri, prior(lkj(2), class = "cor"))
  }

  pri
}

# ==============================================================================
# 4. FIT MODELS
# ==============================================================================

cat("Fitting between-within models...\n\n")

dir.create("models/between_within", showWarnings = FALSE, recursive = TRUE)

# Start with F0 mean /a/ as test case
outcome <- "f0_mean_a"
formula_f0 <- build_formula(outcome, gaussian())

cat("=== MODEL: F0 mean /a/ ===\n")
cat("Formula:\n")
print(formula_f0)
cat("\n")

priors_f0 <- get_priors_bw(outcome, is_lognormal = FALSE)

cat("Priors:\n")
print(priors_f0)
cat("\n")

# Fit
cat("Estimating model (this may take 10-15 minutes)...\n")

fit_f0_bw <- brm(
  formula = formula_f0,
  data = df_complete,
  prior = priors_f0,
  chains = 4,
  iter = 5000,
  warmup = 2500,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "models/between_within/f0_mean_a_bw",
  backend = "cmdstanr"
)

cat("\n✓ Model fitted\n\n")

# ==============================================================================
# 5. MODEL SUMMARY
# ==============================================================================

cat(rep("=", 70), "\n")
cat("MODEL SUMMARY: F0 mean /a/\n")
cat(rep("=", 70), "\n\n")

print(summary(fit_f0_bw))

# ==============================================================================
# 6. EXTRACT KEY PARAMETERS
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("KEY COMPARISONS: Between vs Within Effects\n")
cat(rep("=", 70), "\n\n")

# Extract posterior samples
posterior <- as_draws_df(fit_f0_bw)

# Function to extract and compare between vs within effects
compare_between_within <- function(posterior, domain) {
  # Main effects
  between_main <- paste0("b_pid5_", domain, "_between_c")
  within_main <- paste0("b_pid5_", domain, "_within_c")

  # Stress interactions
  between_stress <- paste0("b_pid5_", domain, "_between_c:c1_stress")
  within_stress <- paste0("b_pid5_", domain, "_within_c:c1_stress")

  # Recovery interactions
  between_recovery <- paste0("b_pid5_", domain, "_between_c:c2_recovery")
  within_recovery <- paste0("b_pid5_", domain, "_within_c:c2_recovery")

  # Extract
  results <- tibble(
    domain = domain,

    # Between effects
    between_main_mean = mean(posterior[[between_main]], na.rm = TRUE),
    between_main_ci_lower = quantile(
      posterior[[between_main]],
      0.025,
      na.rm = TRUE
    ),
    between_main_ci_upper = quantile(
      posterior[[between_main]],
      0.975,
      na.rm = TRUE
    ),

    # Within effects
    within_main_mean = mean(posterior[[within_main]], na.rm = TRUE),
    within_main_ci_lower = quantile(
      posterior[[within_main]],
      0.025,
      na.rm = TRUE
    ),
    within_main_ci_upper = quantile(
      posterior[[within_main]],
      0.975,
      na.rm = TRUE
    ),

    # Stress interactions
    between_stress_mean = mean(posterior[[between_stress]], na.rm = TRUE),
    between_stress_ci_lower = quantile(
      posterior[[between_stress]],
      0.025,
      na.rm = TRUE
    ),
    between_stress_ci_upper = quantile(
      posterior[[between_stress]],
      0.975,
      na.rm = TRUE
    ),

    within_stress_mean = mean(posterior[[within_stress]], na.rm = TRUE),
    within_stress_ci_lower = quantile(
      posterior[[within_stress]],
      0.025,
      na.rm = TRUE
    ),
    within_stress_ci_upper = quantile(
      posterior[[within_stress]],
      0.975,
      na.rm = TRUE
    ),

    # Recovery interactions
    between_recovery_mean = mean(posterior[[between_recovery]], na.rm = TRUE),
    between_recovery_ci_lower = quantile(
      posterior[[between_recovery]],
      0.025,
      na.rm = TRUE
    ),
    between_recovery_ci_upper = quantile(
      posterior[[between_recovery]],
      0.975,
      na.rm = TRUE
    ),

    within_recovery_mean = mean(posterior[[within_recovery]], na.rm = TRUE),
    within_recovery_ci_lower = quantile(
      posterior[[within_recovery]],
      0.025,
      na.rm = TRUE
    ),
    within_recovery_ci_upper = quantile(
      posterior[[within_recovery]],
      0.975,
      na.rm = TRUE
    )
  ) %>%
    mutate(
      # Credible effects (95% CI excludes 0)
      between_main_credible = sign(between_main_ci_lower) ==
        sign(between_main_ci_upper),
      within_main_credible = sign(within_main_ci_lower) ==
        sign(within_main_ci_upper),
      between_stress_credible = sign(between_stress_ci_lower) ==
        sign(between_stress_ci_upper),
      within_stress_credible = sign(within_stress_ci_lower) ==
        sign(within_stress_ci_upper),
      between_recovery_credible = sign(between_recovery_ci_lower) ==
        sign(between_recovery_ci_upper),
      within_recovery_credible = sign(within_recovery_ci_lower) ==
        sign(within_recovery_ci_upper)
    )

  return(results)
}

# Compare for all domains
comparison_results <- map_df(
  pid5_domains,
  ~ compare_between_within(posterior, .x)
)

cat("Between-person (trait) vs Within-person (state) effects:\n\n")
print(
  comparison_results %>%
    select(domain, ends_with("_mean"), ends_with("_credible")),
  n = Inf
)

# Save
saveRDS(comparison_results, "results/between_within/f0_mean_a_comparison.rds")
rio::export(
  comparison_results,
  "results/between_within/f0_mean_a_comparison.csv"
)

# ==============================================================================
# 7. IDENTIFY NEW FINDINGS
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("NEW FINDINGS: Within-person effects beyond trait\n")
cat(rep("=", 70), "\n\n")

# Main effects
within_main_credible <- comparison_results %>%
  filter(within_main_credible == TRUE) %>%
  select(domain, within_main_mean, within_main_ci_lower, within_main_ci_upper)

if (nrow(within_main_credible) > 0) {
  cat("WITHIN-PERSON MAIN EFFECTS (state predicts voice):\n")
  print(within_main_credible, n = Inf)
  cat("\n")
} else {
  cat("No credible within-person main effects found.\n\n")
}

# Stress interactions
within_stress_credible <- comparison_results %>%
  filter(within_stress_credible == TRUE) %>%
  select(
    domain,
    within_stress_mean,
    within_stress_ci_lower,
    within_stress_ci_upper
  )

if (nrow(within_stress_credible) > 0) {
  cat("WITHIN-PERSON × STRESS INTERACTIONS (state moderates stress):\n")
  print(within_stress_credible, n = Inf)
  cat("\n")
} else {
  cat("No credible within-person × stress interactions found.\n\n")
}

# Recovery interactions
within_recovery_credible <- comparison_results %>%
  filter(within_recovery_credible == TRUE) %>%
  select(
    domain,
    within_recovery_mean,
    within_recovery_ci_lower,
    within_recovery_ci_upper
  )

if (nrow(within_recovery_credible) > 0) {
  cat("WITHIN-PERSON × RECOVERY INTERACTIONS (state moderates recovery):\n")
  print(within_recovery_credible, n = Inf)
  cat("\n")
} else {
  cat("No credible within-person × recovery interactions found.\n\n")
}

# ==============================================================================
# 8. DECISION CRITERIA
# ==============================================================================

cat(rep("=", 70), "\n")
cat("DECISION CRITERIA: Include in main manuscript?\n")
cat(rep("=", 70), "\n\n")

n_within_effects <- sum(
  comparison_results$within_main_credible,
  comparison_results$within_stress_credible,
  comparison_results$within_recovery_credible
)

n_between_effects <- sum(
  comparison_results$between_main_credible,
  comparison_results$between_stress_credible,
  comparison_results$between_recovery_credible
)

cat("Summary:\n")
cat("  Between-person effects:", n_between_effects, "\n")
cat("  Within-person effects:", n_within_effects, "\n\n")

if (n_within_effects >= 2) {
  cat("RECOMMENDATION: Include in main manuscript\n")
  cat("Rationale: Multiple credible within-person effects suggest state\n")
  cat("fluctuations contribute meaningful variance beyond stable traits.\n\n")
} else if (n_within_effects == 1) {
  cat("RECOMMENDATION: Include in supplementary materials\n")
  cat(
    "Rationale: Single within-person effect - interesting but not robust.\n\n"
  )
} else {
  cat("RECOMMENDATION: Supplementary materials only\n")
  cat(
    "Rationale: No credible within-person effects - trait-only model sufficient.\n\n"
  )
}

# ==============================================================================
# 9. SAVE OUTPUT
# ==============================================================================

cat("Saving results...\n")

# Save model
saveRDS(fit_f0_bw, "results/between_within/fit_f0_mean_a_bw.rds")

# Summary text
sink("results/between_within/f0_mean_a_summary.txt")
cat("BETWEEN-WITHIN ANALYSIS: F0 mean /a/\n")
cat(rep("=", 70), "\n\n")
cat("Model fitted:", as.character(Sys.time()), "\n\n")
print(summary(fit_f0_bw))
cat("\n\nCOMPARISON: Between vs Within\n")
cat(rep("=", 70), "\n\n")
print(comparison_results)
sink()

cat("✓ Model saved: results/between_within/fit_f0_mean_a_bw.rds\n")
cat("✓ Comparison saved: results/between_within/f0_mean_a_comparison.csv\n")
cat("✓ Summary saved: results/between_within/f0_mean_a_summary.txt\n\n")

cat(rep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Next steps:\n")
cat("1. Review comparison results to decide manuscript vs supplementary\n")
cat("2. If promising, fit additional outcomes (other vowels, F2, etc.)\n")
cat("3. Run 13_compare_between_within.R for formal model comparison\n\n")
