# ==============================================================================
# 08_temporal_covariation_analysis.R
# Temporal Covariation: Do EMA trait levels co-vary with acoustic changes?
# ==============================================================================
# OBIETTIVO:
#   Testare se i livelli TIME-VARYING dei tratti EMA predicono
#   i cambiamenti acustici (coupling dinamico trait-contesto)
#
# MOTIVAZIONE:
#   L'abstract afferma: "EMA measures...co-varied temporally with acoustic changes"
#   Questa analisi supporta empiricamente questo claim.
#
# APPROCCIO:
#   1. Usare valori EMA aggregati PER TIMEPOINT (non overall mean)
#   2. Testare se quando un tratto aumenta/diminuisce, anche la voce cambia
#   3. Confrontare con modelli che usano solo trait between-person
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayestestR)
  library(tidybayes)
  library(ggplot2)
})

options(brms.backend = "cmdstanr")

cat("\n", rep("=", 70), "\n", sep = "")
cat("TEMPORAL COVARIATION ANALYSIS\n")
cat("Do EMA trait fluctuations co-vary with acoustic changes?\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

# Load cleaned EMA data
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
d <- rio::import(data_path)

# Load voice data
if (file.exists("results/df_analysis.rds")) {
  voice_data <- readRDS("results/df_analysis.rds")
  cat("✓ Loaded voice data from results/df_analysis.rds\n")
} else {
  stop("Esegui prima 02_voice_personality_analysis_FINAL.R")
}

# EMA PID-5 variables
ema_pid5_vars <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

# ==============================================================================
# 2. PREPARE TIME-VARYING PID-5 MEASURES
# ==============================================================================

cat("\n=== CREATING TIME-VARYING PID-5 MEASURES ===\n")

# Recode exam_period to match voice timepoint
d_ema <- d %>%
  dplyr::rename(ID = user_id) %>%
  mutate(
    timepoint = forcats::fct_recode(
      exam_period,
      "baseline" = "baseline",
      "pre" = "pre_exam",
      "post" = "post_exam"
    )
  ) %>%
  filter(!is.na(timepoint))

# Aggregate EMA per ID × timepoint
pid5_timevarying <- d_ema %>%
  group_by(ID, timepoint) %>%
  summarise(
    across(all_of(ema_pid5_vars), ~ mean(.x, na.rm = TRUE)),
    n_ema = n(),
    .groups = "drop"
  )

cat("Time-varying PID-5 data:\n")
cat("  N observations:", nrow(pid5_timevarying), "\n")
cat("  N subjects:", n_distinct(pid5_timevarying$ID), "\n")
cat("  Mean EMA assessments per timepoint:", 
    round(mean(pid5_timevarying$n_ema, na.rm = TRUE), 1), "\n\n")

# Check coverage
coverage <- pid5_timevarying %>%
  group_by(timepoint) %>%
  summarise(n_subjects = n_distinct(ID), .groups = "drop")

cat("Coverage by timepoint:\n")
print(coverage)

# ==============================================================================
# 3. COMPUTE WITHIN-PERSON CENTERING
# ==============================================================================

cat("\n=== DECOMPOSING BETWEEN vs WITHIN-PERSON VARIANCE ===\n")

# Person means (between-person component)
pid5_between <- pid5_timevarying %>%
  group_by(ID) %>%
  summarise(
    across(all_of(ema_pid5_vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_between"),
    .groups = "drop"
  )

# Join and compute within-person deviations
pid5_decomposed <- pid5_timevarying %>%
  left_join(pid5_between, by = "ID") %>%
  mutate(
    # Within-person deviations (time-varying component)
    across(
      all_of(ema_pid5_vars),
      ~ .x - get(paste0(cur_column(), "_between")),
      .names = "{.col}_within"
    )
  )

cat("Within-person variance (SD of deviations):\n")
within_vars <- pid5_decomposed %>%
  summarise(
    across(
      ends_with("_within"),
      ~ sd(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "sd")

print(within_vars)

# ==============================================================================
# 4. JOIN WITH VOICE DATA
# ==============================================================================

cat("\n=== JOINING WITH VOICE DATA ===\n")

# Select voice outcomes (vowel /a/ for simplicity)
voice_simple <- voice_data %>%
  dplyr::select(
    ID, timepoint,
    f0_mean_a, f0_std_a, jitter_a, nne_a, f2_mean_a, f2_std_a
  )

# Join time-varying PID-5
df_temporal <- voice_simple %>%
  left_join(pid5_decomposed, by = c("ID", "timepoint"))

# Check missingness
n_complete <- sum(complete.cases(df_temporal))
cat("Complete cases:", n_complete, "/", nrow(df_temporal), "\n\n")

# ==============================================================================
# 5. BAYESIAN MODELS: TIME-VARYING EFFECTS
# ==============================================================================

cat("\n=== FITTING TIME-VARYING MODELS ===\n")
cat("Testing: Do within-person fluctuations in traits predict acoustic changes?\n\n")

# Model formula: Outcome ~ Between + Within + (1|ID)
# Between = stable individual differences
# Within = time-varying deviations (temporal covariation)

# Example: Negative Affectivity predicting F0 mean
cat("--- Model: F0 Mean ~ NA (between + within) ---\n")

m_temporal_na_f0 <- brm(
  f0_mean_a ~ 
    pid5_negative_affectivity_between +  # Between-person effect
    pid5_negative_affectivity_within +   # Within-person (time-varying) effect
    (1 | ID),
  data = df_temporal,
  family = gaussian(),
  prior = c(
    prior(student_t(3, 220, 30), class = Intercept),
    prior(normal(0, 10), class = b),
    prior(exponential(0.3), class = sigma),
    prior(exponential(0.3), class = sd)
  ),
  iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 123,
  control = list(adapt_delta = 0.95),
  file = "models/m_temporal_na_f0"
)

cat("\nResults:\n")
print(summary(m_temporal_na_f0)$fixed)

cat("\nInterpretation:\n")
fe <- fixef(m_temporal_na_f0)

if ("pid5_negative_affectivity_between" %in% rownames(fe)) {
  beta_between <- fe["pid5_negative_affectivity_between", "Estimate"]
  ci_between <- fe["pid5_negative_affectivity_between", c("Q2.5", "Q97.5")]
  
  cat("Between-person effect: β =", round(beta_between, 3), 
      ", 95% CI [", round(ci_between[1], 3), ",", round(ci_between[2], 3), "]\n")
  cat("  → People with higher NA have", 
      ifelse(beta_between > 0, "higher", "lower"), "baseline F0\n")
}

if ("pid5_negative_affectivity_within" %in% rownames(fe)) {
  beta_within <- fe["pid5_negative_affectivity_within", "Estimate"]
  ci_within <- fe["pid5_negative_affectivity_within", c("Q2.5", "Q97.5")]
  
  cat("Within-person effect: β =", round(beta_within, 3),
      ", 95% CI [", round(ci_within[1], 3), ",", round(ci_within[2], 3), "]\n")
  
  if (ci_within[1] > 0) {
    cat("  ✓ TEMPORAL COVARIATION: When NA increases, F0 increases\n")
  } else if (ci_within[2] < 0) {
    cat("  ✓ TEMPORAL COVARIATION: When NA increases, F0 decreases\n")
  } else {
    cat("  • No clear within-person coupling (CI includes zero)\n")
  }
}

# ==============================================================================
# 6. SYSTEMATIC TEST ACROSS OUTCOMES
# ==============================================================================

cat("\n\n=== SYSTEMATIC TESTING: All Traits × F0 Mean ===\n")

# Function to fit and summarize temporal model
fit_temporal_model <- function(outcome, trait_name, data) {
  between_var <- paste0(trait_name, "_between")
  within_var <- paste0(trait_name, "_within")
  
  fmla <- as.formula(paste0(
    outcome, " ~ ", between_var, " + ", within_var, " + (1 | ID)"
  ))
  
  # Use appropriate prior based on outcome
  if (outcome == "f0_mean_a") {
    prior_here <- c(
      prior(student_t(3, 220, 30), class = Intercept),
      prior(normal(0, 10), class = b),
      prior(exponential(0.3), class = sigma),
      prior(exponential(0.3), class = sd)
    )
    fam <- gaussian()
  } else if (outcome %in% c("f0_std_a", "jitter_a", "f2_std_a")) {
    prior_here <- c(
      prior(student_t(3, 0, 1), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(1), class = sigma),
      prior(exponential(0.5), class = sd)
    )
    fam <- lognormal()
  } else if (outcome == "f2_mean_a") {
    prior_here <- c(
      prior(student_t(3, 1500, 150), class = Intercept),
      prior(normal(0, 50), class = b),
      prior(exponential(0.01), class = sigma),
      prior(exponential(0.01), class = sd),
      prior(gamma(2, 0.1), class = nu)
    )
    fam <- student()
  } else {
    prior_here <- c(
      prior(student_t(3, -20, 5), class = Intercept),
      prior(normal(0, 5), class = b),
      prior(exponential(0.5), class = sigma),
      prior(exponential(0.3), class = sd)
    )
    fam <- gaussian()
  }
  
  model_name <- paste0("m_temporal_", gsub("pid5_", "", trait_name), "_", outcome)
  
  m <- brm(
    formula = fmla,
    data = data,
    family = fam,
    prior = prior_here,
    iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 123,
    control = list(adapt_delta = 0.95),
    file = file.path("models", model_name),
    silent = 2, refresh = 0
  )
  
  fe <- fixef(m)
  
  tibble(
    outcome = outcome,
    trait = trait_name,
    beta_between = fe[between_var, "Estimate"],
    ci_between_lower = fe[between_var, "Q2.5"],
    ci_between_upper = fe[between_var, "Q97.5"],
    beta_within = fe[within_var, "Estimate"],
    ci_within_lower = fe[within_var, "Q2.5"],
    ci_within_upper = fe[within_var, "Q97.5"],
    sig_between = (fe[between_var, "Q2.5"] > 0) | (fe[between_var, "Q97.5"] < 0),
    sig_within = (fe[within_var, "Q2.5"] > 0) | (fe[within_var, "Q97.5"] < 0)
  )
}

# Fit models for all traits
traits_to_test <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

temporal_results <- map_df(traits_to_test, function(trait) {
  cat("Testing:", trait, "...\n")
  fit_temporal_model("f0_mean_a", trait, df_temporal)
})

cat("\n=== TEMPORAL COVARIATION RESULTS (F0 Mean) ===\n\n")
print(temporal_results %>% select(trait, beta_within, ci_within_lower, ci_within_upper, sig_within))

# ==============================================================================
# 7. VISUALIZATION
# ==============================================================================

cat("\n=== CREATING VISUALIZATION ===\n")

# Plot within-person effects
p_temporal <- ggplot(
  temporal_results,
  aes(
    x = beta_within,
    y = trait,
    xmin = ci_within_lower,
    xmax = ci_within_upper,
    color = sig_within
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
    labels = c("FALSE" = "Not significant", "TRUE" = "Significant"),
    name = "95% CI excludes zero"
  ) +
  labs(
    title = "Temporal Covariation: Within-Person Effects on F0",
    subtitle = "Do EMA trait fluctuations predict acoustic changes?",
    x = "Within-Person Effect (β)",
    y = "PID-5 Domain"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("figures/temporal_covariation_f0.png", p_temporal, 
       width = 10, height = 6, dpi = 300)

cat("✓ Plot saved: figures/temporal_covariation_f0.png\n")

# ==============================================================================
# 8. SUMMARY & INTERPRETATION
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("TEMPORAL COVARIATION SUMMARY\n")
cat(rep("=", 70), "\n\n")

sig_within <- temporal_results %>% filter(sig_within)

if (nrow(sig_within) > 0) {
  cat("✓ EVIDENCE FOR TEMPORAL COVARIATION:\n")
  cat("The following traits show significant within-person coupling with F0:\n\n")
  
  for (i in 1:nrow(sig_within)) {
    cat(" -", sig_within$trait[i], ":\n")
    cat("   β =", round(sig_within$beta_within[i], 3),
        ", 95% CI [", round(sig_within$ci_within_lower[i], 3), ",",
        round(sig_within$ci_within_upper[i], 3), "]\n")
    cat("   Interpretation: When this trait increases, F0",
        ifelse(sig_within$beta_within[i] > 0, "increases", "decreases"), "\n\n")
  }
  
  cat("CONCLUSION: EMA measures co-vary temporally with acoustic changes.\n")
  cat("This supports the claim in the abstract.\n")
} else {
  cat("⚠ LIMITED EVIDENCE for temporal covariation.\n")
  cat("Between-person effects are present, but within-person coupling is weak.\n")
  cat("This suggests traits are relatively stable across the exam period.\n")
}

# ==============================================================================
# 9. SAVE RESULTS
# ==============================================================================

saveRDS(temporal_results, "results/temporal_covariation_results.rds")
saveRDS(df_temporal, "results/df_temporal.rds")

cat("\n✓ Results saved:\n")
cat("  - results/temporal_covariation_results.rds\n")
cat("  - results/df_temporal.rds\n\n")

cat(rep("=", 70), "\n", sep = "")
cat("TEMPORAL COVARIATION ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n")
