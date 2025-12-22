# ==============================================================================
# 15_extend_to_all_outcomes.R
# Extend between-within analysis to all acoustic outcomes
# ==============================================================================
# PURPOSE:
#   If F0 results are promising, extend to all 18 outcomes:
#   - 6 features (f0_mean, f0_std, f2_mean, f2_std, jitter, nne)
#   - 3 vowels (/a/, /i/, /u/)
#   - Total: 18 models
#
# RUN THIS ONLY IF:
#   - F0 mean /a/ showed credible within-person effects, OR
#   - Model comparison favored trait+state model
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(furrr)
  library(progressr)
})

options(brms.backend = "cmdstanr")
plan(multisession, workers = 4)  # Parallel processing

cat("\n", rep("=", 70), "\n", sep = "")
cat("EXTENDING BETWEEN-WITHIN ANALYSIS TO ALL OUTCOMES\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 0. CHECK IF EXTENSION IS WARRANTED
# ==============================================================================

if (!file.exists("results/between_within/model_comparison_summary.rds")) {
  stop("Run 13_compare_between_within.R first to determine if extension is warranted")
}

model_comp <- readRDS("results/between_within/model_comparison_summary.rds")
decision <- model_comp %>% filter(metric == "decision") %>% pull(interpretation)

cat("Decision from F0 mean /a/ analysis:", decision, "\n\n")

if (decision == "supplementary_brief") {
  cat("⚠ WARNING: F0 analysis recommended supplementary-only status\n")
  cat("Extension to all outcomes may not be necessary.\n\n")
  cat("Do you want to continue? (y/n): ")
  response <- readline()
  if (tolower(response) != "y") {
    stop("Extension cancelled by user")
  }
}

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

df <- readRDS("results/between_within/df_between_within.rds")

# Complete cases
df_complete <- df %>%
  filter(!is.na(pid5_negative_affectivity_between_c))

cat("Data loaded: N =", nrow(df_complete), "observations\n\n")

# ==============================================================================
# 2. DEFINE OUTCOMES AND FAMILIES
# ==============================================================================

# All combinations of features and vowels
features <- c("f0_mean", "f0_std", "f2_mean", "f2_std", "jitter", "nne")
vowels <- c("a", "i", "u")

outcomes_grid <- expand.grid(
  feature = features,
  vowel = vowels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    outcome = paste0(feature, "_", vowel),
    family = case_when(
      feature %in% c("f0_std", "f2_std", "jitter") ~ "lognormal",
      feature == "f2_mean" ~ "student",  # Robust for occasional outliers
      TRUE ~ "gaussian"
    )
  ) %>%
  # Skip if outcome doesn't exist in data
  filter(outcome %in% names(df_complete))

cat("Outcomes to model:\n")
print(outcomes_grid, n = Inf)
cat("\nTotal:", nrow(outcomes_grid), "models to fit\n\n")

# ==============================================================================
# 3. MODEL SPECIFICATION FUNCTION
# ==============================================================================

build_bw_formula <- function(outcome, family_name) {
  
  pid5_domains <- c("negative_affectivity", "detachment", "antagonism",
                   "disinhibition", "psychoticism")
  
  # Between-person terms
  between_terms <- paste0("pid5_", pid5_domains, "_between_c")
  between_str <- paste(between_terms, collapse = " + ")
  
  # Within-person terms
  within_terms <- paste0("pid5_", pid5_domains, "_within_c")
  within_str <- paste(within_terms, collapse = " + ")
  
  # Interactions with stress
  between_stress <- paste0(between_terms, ":c1_stress")
  within_stress <- paste0(within_terms, ":c1_stress")
  stress_str <- paste(c(between_stress, within_stress), collapse = " + ")
  
  # Interactions with recovery
  between_recovery <- paste0(between_terms, ":c2_recovery")
  within_recovery <- paste0(within_terms, ":c2_recovery")
  recovery_str <- paste(c(between_recovery, within_recovery), collapse = " + ")
  
  # Complete formula
  formula_str <- paste0(
    outcome, " ~ 1 + c1_stress + c2_recovery + ",
    between_str, " + ", within_str, " + ",
    stress_str, " + ", recovery_str,
    " + (1 + c1_stress + c2_recovery || ID)"
  )
  
  # Family
  fam <- switch(family_name,
                "gaussian" = gaussian(),
                "lognormal" = lognormal(),
                "student" = student())
  
  bf(formula_str, family = fam)
}

# ==============================================================================
# 4. PRIORS FUNCTION
# ==============================================================================

get_priors <- function(outcome, family_name, data) {
  
  emp_mean <- if (family_name == "lognormal") {
    exp(mean(log(data[[outcome]][data[[outcome]] > 0]), na.rm = TRUE))
  } else {
    mean(data[[outcome]], na.rm = TRUE)
  }
  
  emp_sd <- if (family_name == "lognormal") {
    sd(log(data[[outcome]][data[[outcome]] > 0]), na.rm = TRUE)
  } else {
    sd(data[[outcome]], na.rm = TRUE)
  }
  
  if (family_name == "lognormal") {
    prior(student_t(3, log(emp_mean), emp_sd), class = Intercept) +
    prior(normal(0, 0.5), class = b) +
    prior(exponential(1), class = sigma) +
    prior(exponential(1), class = sd)
  } else if (family_name == "student") {
    prior(student_t(3, emp_mean, emp_sd), class = Intercept) +
    prior(normal(0, 10), class = b) +
    prior(gamma(2, 0.1), class = nu) +
    prior(exponential(0.1), class = sigma) +
    prior(exponential(1), class = sd)
  } else {
    prior(student_t(3, emp_mean, emp_sd), class = Intercept) +
    prior(normal(0, 10), class = b) +
    prior(exponential(0.1), class = sigma) +
    prior(exponential(1), class = sd)
  }
}

# ==============================================================================
# 5. FIT ALL MODELS IN PARALLEL
# ==============================================================================

cat("Fitting", nrow(outcomes_grid), "models in parallel...\n")
cat("This will take approximately", 
    round(nrow(outcomes_grid) * 10 / 4), "minutes with 4 cores\n\n")

dir.create("models/between_within/all_outcomes", showWarnings = FALSE, recursive = TRUE)

# Function to fit one model
fit_one_model <- function(i, outcomes_df, data) {
  
  outcome <- outcomes_df$outcome[i]
  family_name <- outcomes_df$family[i]
  
  # Build formula and priors
  form <- build_bw_formula(outcome, family_name)
  priors <- get_priors(outcome, family_name, data)
  
  # File path
  file_path <- paste0("models/between_within/all_outcomes/", outcome, "_bw")
  
  # Fit
  fit <- brm(
    formula = form,
    data = data,
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 1,  # 1 per model since we're parallelizing across models
    seed = 123 + i,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    file = file_path,
    backend = "cmdstanr",
    silent = 2,
    refresh = 0
  )
  
  return(list(
    outcome = outcome,
    file = file_path,
    success = TRUE
  ))
}

# Fit with progress bar
with_progress({
  p <- progressor(steps = nrow(outcomes_grid))
  
  results <- future_map(1:nrow(outcomes_grid), function(i) {
    p(sprintf("Fitting %s", outcomes_grid$outcome[i]))
    fit_one_model(i, outcomes_grid, df_complete)
  }, .options = furrr_options(seed = TRUE))
})

cat("\n✓ All models fitted\n\n")

# ==============================================================================
# 6. EXTRACT RESULTS FROM ALL MODELS
# ==============================================================================

cat("Extracting results from all models...\n")

extract_comparison <- function(outcome_name, family_name) {
  
  # Load model
  file_path <- paste0("models/between_within/all_outcomes/", outcome_name, "_bw.rds")
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  fit <- readRDS(file_path)
  posterior <- as_draws_df(fit)
  
  pid5_domains <- c("negative_affectivity", "detachment", "antagonism",
                   "disinhibition", "psychoticism")
  
  # Extract for each domain
  results <- map_df(pid5_domains, function(domain) {
    
    # Parameter names
    between_main <- paste0("b_pid5_", domain, "_between_c")
    within_main <- paste0("b_pid5_", domain, "_within_c")
    between_stress <- paste0("b_pid5_", domain, "_between_c:c1_stress")
    within_stress <- paste0("b_pid5_", domain, "_within_c:c1_stress")
    between_recovery <- paste0("b_pid5_", domain, "_between_c:c2_recovery")
    within_recovery <- paste0("b_pid5_", domain, "_within_c:c2_recovery")
    
    # Extract
    tibble(
      outcome = outcome_name,
      domain = domain,
      
      # Within main
      within_main_mean = mean(posterior[[within_main]]),
      within_main_ci_lower = quantile(posterior[[within_main]], 0.025),
      within_main_ci_upper = quantile(posterior[[within_main]], 0.975),
      within_main_credible = sign(within_main_ci_lower) == sign(within_main_ci_upper),
      
      # Within stress
      within_stress_mean = mean(posterior[[within_stress]]),
      within_stress_ci_lower = quantile(posterior[[within_stress]], 0.025),
      within_stress_ci_upper = quantile(posterior[[within_stress]], 0.975),
      within_stress_credible = sign(within_stress_ci_lower) == sign(within_stress_ci_upper),
      
      # Within recovery
      within_recovery_mean = mean(posterior[[within_recovery]]),
      within_recovery_ci_lower = quantile(posterior[[within_recovery]], 0.025),
      within_recovery_ci_upper = quantile(posterior[[within_recovery]], 0.975),
      within_recovery_credible = sign(within_recovery_ci_lower) == sign(within_recovery_ci_upper),
      
      # Between main (for comparison)
      between_main_mean = mean(posterior[[between_main]]),
      between_main_credible = sign(quantile(posterior[[between_main]], 0.025)) == 
                             sign(quantile(posterior[[between_main]], 0.975))
    )
  })
  
  return(results)
}

# Extract all
all_comparisons <- map2_df(
  outcomes_grid$outcome,
  outcomes_grid$family,
  extract_comparison
)

# Save
saveRDS(all_comparisons, "results/between_within/all_outcomes_comparison.rds")
rio::export(all_comparisons, "results/between_within/all_outcomes_comparison.csv")

cat("✓ Results extracted and saved\n\n")

# ==============================================================================
# 7. SUMMARY OF WITHIN-PERSON EFFECTS ACROSS OUTCOMES
# ==============================================================================

cat(rep("=", 70), "\n")
cat("SUMMARY: Within-Person Effects Across All Outcomes\n")
cat(rep("=", 70), "\n\n")

# Count credible effects by outcome
within_summary <- all_comparisons %>%
  group_by(outcome) %>%
  summarise(
    n_within_main = sum(within_main_credible),
    n_within_stress = sum(within_stress_credible),
    n_within_recovery = sum(within_recovery_credible),
    n_within_total = n_within_main + n_within_stress + n_within_recovery,
    n_between_main = sum(between_main_credible),
    .groups = "drop"
  ) %>%
  arrange(desc(n_within_total))

cat("Within-person effects by outcome:\n\n")
print(within_summary, n = Inf)

# Count by domain
domain_summary <- all_comparisons %>%
  group_by(domain) %>%
  summarise(
    n_within_main = sum(within_main_credible),
    n_within_stress = sum(within_stress_credible),
    n_within_recovery = sum(within_recovery_credible),
    n_within_total = n_within_main + n_within_stress + n_within_recovery,
    .groups = "drop"
  ) %>%
  arrange(desc(n_within_total))

cat("\n\nWithin-person effects by domain:\n\n")
print(domain_summary, n = Inf)

# ==============================================================================
# 8. IDENTIFY MOST INTERESTING OUTCOMES
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("MOST INTERESTING OUTCOMES FOR MANUSCRIPT\n")
cat(rep("=", 70), "\n\n")

# Outcomes with most within-person effects
top_outcomes <- within_summary %>%
  filter(n_within_total >= 2) %>%
  arrange(desc(n_within_total))

if (nrow(top_outcomes) > 0) {
  cat("Outcomes with ≥2 credible within-person effects:\n\n")
  print(top_outcomes, n = Inf)
  cat("\nRECOMMENDATION: Include these in main manuscript\n\n")
} else {
  cat("No outcomes with ≥2 credible within-person effects\n")
  cat("RECOMMENDATION: Supplementary materials only\n\n")
}

# Save summary
saveRDS(within_summary, "results/between_within/within_effects_summary.rds")
rio::export(within_summary, "results/between_within/within_effects_summary.csv")

cat("✓ Summary saved\n\n")

# ==============================================================================
# 9. CREATE MASTER SUMMARY TABLE
# ==============================================================================

cat("Creating master summary table...\n")

# All credible within-person effects
all_credible_within <- all_comparisons %>%
  filter(within_main_credible | within_stress_credible | within_recovery_credible) %>%
  select(outcome, domain, 
         starts_with("within_main"),
         starts_with("within_stress"),
         starts_with("within_recovery")) %>%
  arrange(outcome, domain)

rio::export(all_credible_within, 
            "results/between_within/all_credible_within_effects.csv")

cat("✓ Master table saved: all_credible_within_effects.csv\n\n")

# ==============================================================================
# 10. FINAL RECOMMENDATION
# ==============================================================================

cat(rep("=", 70), "\n")
cat("FINAL RECOMMENDATION\n")
cat(rep("=", 70), "\n\n")

total_within <- sum(within_summary$n_within_total)
total_between <- sum(within_summary$n_between_main)

cat("Total effects across all outcomes:\n")
cat("  Between-person:", total_between, "\n")
cat("  Within-person:", total_within, "\n\n")

if (total_within >= 10) {
  cat("STRONG RECOMMENDATION: Include between-within decomposition in main manuscript\n\n")
  cat("Rationale:\n")
  cat("- Substantial number of within-person effects (", total_within, ")\n")
  cat("- State fluctuations clearly contribute beyond stable traits\n")
  cat("- Important theoretical implication for personality-stress dynamics\n\n")
  
} else if (total_within >= 5) {
  cat("MODERATE RECOMMENDATION: Include in main manuscript or extensive supplement\n\n")
  cat("Rationale:\n")
  cat("- Moderate number of within-person effects (", total_within, ")\n")
  cat("- Some evidence for state contributions\n")
  cat("- Decision depends on journal space and theoretical focus\n\n")
  
} else {
  cat("RECOMMENDATION: Supplementary materials with brief mention\n\n")
  cat("Rationale:\n")
  cat("- Limited within-person effects (", total_within, ")\n")
  cat("- Trait-only model appears sufficient for main story\n")
  cat("- State component adds complexity without major insights\n\n")
}

cat(rep("=", 70), "\n")
cat("EXTENSION ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("All results saved in: results/between_within/\n")
cat("All models saved in: models/between_within/all_outcomes/\n\n")
