# ==============================================================================
# 02_voice_personality_analysis_FINAL.R
# Bayesian Analysis: Voice Acoustics, Stress, and Personality Pathology
# VERSIONE FINALE con validazione EMA ↔ PID-5 full
# ==============================================================================
# AGGIUNTE IN QUESTA VERSIONE:
#   1. Validazione convergente: EMA (15 item) vs PID-5 full (220 item)
#   2. Statistiche descrittive più complete
#   3. Sample description (female only)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(brms)
  library(cmdstanr)
  library(bayestestR)
  library(bayesplot)
  library(tidybayes)
  library(patchwork)
  library(ggdist)
  library(lubridate)
  library(stringr)
  library(missRanger)
  library(rio)
})

options(brms.backend = "cmdstanr")

dir.create("models", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

PID5_VERSION <- "ema_all"  # Media di tutti i valori EMA per soggetto

cat("\n=== CONFIGURATION ===\n")
cat("PID-5 version:", PID5_VERSION, "\n")
cat("Sample: Female university students only\n\n")

# ==============================================================================
# 1. DATA PREPARATION - Voice Data
# ==============================================================================

baseline <- read_excel(
  here::here("data", "raw", "acustic_features", "datiacustici", "AUDIO.xlsx"),
  sheet = "BASELINE"
)
pre <- read_excel(
  here::here("data", "raw", "acustic_features", "datiacustici", "AUDIO.xlsx"),
  sheet = "PRE"
)
post <- read_excel(
  here::here("data", "raw", "acustic_features", "datiacustici", "AUDIO.xlsx"),
  sheet = "POST"
)

baseline$timepoint <- "baseline"
pre$timepoint <- "pre"
post$timepoint <- "post"

df_wide <- bind_rows(baseline, pre, post)
names(df_wide) <- str_trim(names(df_wide))

# Correct ID names
df_wide <- df_wide %>%
  mutate(
    ID = case_when(
      ID == "am_bo_1988_08_24_166" ~ "an_bo_1988_08_24_166",
      ID == "as_li_2005_04_26_447" ~ "as_si_2005_04_26_447",
      ID == "cl_bo_1987_10_16_628" ~ "ca_bo_1987_10_16_628",
      ID == "hi_na_2005_03_08_339" ~ "gi_na_2005_03_08_339",
      ID == "ma_si_2003_10_31_940" ~ "si_ma_2003_10_31_940",
      TRUE ~ ID
    )
  )

cat("Voice data: N subjects =", n_distinct(df_wide$ID), "\n")

df_clean <- df_wide %>%
  dplyr::select(
    ID, Data, timepoint,
    f0_mean_a = `F0 mean Hz /a/`,
    f0_std_a = `F0 std Hz /a/`,
    jitter_a = `Jitter /a/`,
    nne_a = `NNE /a/`,
    f2_mean_a = `F2 mean Hz /a/`,
    f2_std_a = `F2 Std Hz /a/`,
    f0_mean_i = `F0 mean Hz /i/`,
    f0_std_i = `F0 std Hz /i/`,
    jitter_i = `Jitter /i/`,
    nne_i = `NNE /i/`,
    f2_mean_i = `F2 mean Hz /i/`,
    f2_std_i = `F2 Std Hz /i/`,
    f0_mean_u = `F0 mean Hz /u/`,
    f0_std_u = `F0 std Hz /u/`,
    jitter_u = `Jitter /u/`,
    nne_u = `NNE /u/`,
    f2_mean_u = `F2 mean Hz /u/`,
    f2_std_u = `F2 Std Hz /u/`
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("baseline", "pre", "post")),
    date = as.Date(Data)
  ) %>%
  dplyr::select(-Data)

# ==============================================================================
# 2. LOAD EMA DATA
# ==============================================================================

data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
d <- rio::import(data_path)

# Variabili EMA PID-5 (15 item, 3 per dominio)
ema_pid5_vars <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

# Variabili PID-5 completo baseline
baseline_pid5_vars <- c(
  "domain_negative_affect_baseline",
  "domain_detachment_baseline",
  "domain_antagonism_baseline",
  "domain_disinhibition_baseline",
  "domain_psychoticism_baseline"
)

d1 <- d %>%
  dplyr::select(
    user_id, exam_period, day, sex,
    all_of(ema_pid5_vars),
    all_of(baseline_pid5_vars)
  ) %>%
  dplyr::rename(ID = user_id)

# Filtra solo soggetti con dati vocali
allowed_id <- unique(df_clean$ID)
d1 <- d1 %>% filter(ID %in% allowed_id)

cat("EMA data: N subjects =", n_distinct(d1$ID), "\n")

# ==============================================================================
# 3. AGGREGAZIONE PID-5 EMA (media across all timepoints)
# ==============================================================================

cat("\n=== AGGREGATING EMA PID-5 ===\n")

pid5_person <- d1 %>%
  group_by(ID) %>%
  summarise(
    # EMA aggregato (media su 2.5 mesi)
    across(all_of(ema_pid5_vars), ~ mean(.x, na.rm = TRUE)),
    sex = first(na.omit(sex)),
    
    # Mantieni anche baseline PID-5 per validazione
    across(all_of(baseline_pid5_vars), ~ first(na.omit(.x))),
    
    n_ema_obs = n(),
    .groups = "drop"
  )

cat("N subjects with EMA:", nrow(pid5_person), "\n")
cat("Mean EMA observations per subject:", round(mean(pid5_person$n_ema_obs), 1), "\n")

# ==============================================================================
# 4. VALIDAZIONE: EMA vs PID-5 FULL
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("CONVERGENT VALIDATION: EMA (15 items) vs PID-5 Full (220 items)\n")
cat(rep("=", 70), "\n\n", sep = "")

# Calcola correlazioni
validation_cors <- tibble(
  domain = c("Negative Affectivity", "Detachment", "Antagonism", 
             "Disinhibition", "Psychoticism"),
  ema_var = ema_pid5_vars,
  full_var = baseline_pid5_vars
) %>%
  rowwise() %>%
  mutate(
    r = cor(
      pid5_person[[ema_var]], 
      pid5_person[[full_var]], 
      use = "complete.obs"
    ),
    n = sum(!is.na(pid5_person[[ema_var]]) & !is.na(pid5_person[[full_var]]))
  ) %>%
  ungroup()

cat("Correlations between EMA and Full PID-5:\n\n")
print(validation_cors %>% select(domain, r, n), n = Inf)

# Plot
p_validation <- ggplot(validation_cors, aes(x = domain, y = r)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = sprintf("r = %.2f", r)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Convergent Validity: EMA vs Full PID-5",
    subtitle = "Correlations between 15-item EMA and 220-item baseline assessment",
    x = "Domain",
    y = "Pearson r"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/ema_validation.png", p_validation, width = 8, height = 6, dpi = 300)

cat("\n✓ Validation plot saved: figures/ema_validation.png\n")

# Salva risultati
saveRDS(validation_cors, "results/ema_validation_correlations.rds")

# Interpretazione
cat("\nInterpretation:\n")
if (all(validation_cors$r > 0.7)) {
  cat("✓ STRONG convergence: All correlations > 0.70\n")
  cat("  EMA measures are valid proxies for full PID-5 domains.\n")
} else if (all(validation_cors$r > 0.5)) {
  cat("✓ MODERATE-STRONG convergence: All correlations > 0.50\n")
  cat("  EMA measures capture substantial variance in full domains.\n")
} else {
  cat("⚠ MIXED convergence: Some correlations < 0.50\n")
  cat("  Consider domains individually in interpretation.\n")
}

weak_domains <- validation_cors %>% 
  filter(r < 0.5) %>% 
  pull(domain)

if (length(weak_domains) > 0) {
  cat("\nDomains with weaker convergence (r < 0.50):\n")
  cat(paste(" -", weak_domains, collapse = "\n"), "\n")
}

# ==============================================================================
# 5. JOIN VOICE + PERSONALITY
# ==============================================================================

voice_df <- df_clean
joined <- voice_df %>%
  left_join(pid5_person, by = "ID")

# Verifica missing
n_missing_pid5 <- sum(is.na(joined$pid5_negative_affectivity))
if (n_missing_pid5 > 0) {
  cat("\nWarning:", n_missing_pid5, "observations missing PID-5 data\n")
}

# ==============================================================================
# 6. IMPUTATION (if needed)
# ==============================================================================

vars_to_impute <- c(
  "f0_mean_a", "f0_std_a", "jitter_a", "nne_a", "f2_mean_a", "f2_std_a",
  "f0_mean_i", "f0_std_i", "jitter_i", "nne_i", "f2_mean_i", "f2_std_i",
  "f0_mean_u", "f0_std_u", "jitter_u", "nne_u", "f2_mean_u", "f2_std_u",
  ema_pid5_vars
)

df_to_impute <- joined %>%
  dplyr::select(ID, timepoint, all_of(vars_to_impute))

n_missing_before <- sum(is.na(df_to_impute[, vars_to_impute]))

if (n_missing_before > 0) {
  cat("\nImputing", n_missing_before, "missing values...\n")
  
  set.seed(123)
  df_imputed <- missRanger(
    df_to_impute,
    pmm.k = 3,
    num.trees = 100,
    verbose = 0
  )
  
  joined_final <- joined %>%
    dplyr::select(-any_of(vars_to_impute)) %>%
    left_join(df_imputed, by = c("ID", "timepoint"))
} else {
  joined_final <- joined
}

# ==============================================================================
# 7. CREATE ANALYSIS DATASET
# ==============================================================================

df_analysis <- joined_final %>%
  mutate(
    # Center PID-5 variables (between-person)
    pid5_negative_affectivity_c = scale(pid5_negative_affectivity)[, 1],
    pid5_detachment_c = scale(pid5_detachment)[, 1],
    pid5_antagonism_c = scale(pid5_antagonism)[, 1],
    pid5_disinhibition_c = scale(pid5_disinhibition)[, 1],
    pid5_psychoticism_c = scale(pid5_psychoticism)[, 1]
  ) %>%
  filter(!is.na(ID), !is.na(timepoint))

cat("\n=== FINAL ANALYSIS DATASET ===\n")
cat("N observations:", nrow(df_analysis), "\n")
cat("N subjects:", n_distinct(df_analysis$ID), "\n")
cat("Observations per subject:", round(nrow(df_analysis) / n_distinct(df_analysis$ID), 1), "\n")

# ==============================================================================
# 8. DESCRIPTIVE STATISTICS: PID-5
# ==============================================================================

cat("\n=== DESCRIPTIVE STATISTICS: PID-5 DOMAINS ===\n\n")

pid5_descriptives <- df_analysis %>%
  distinct(ID, .keep_all = TRUE) %>%
  summarise(
    across(
      c(
        pid5_negative_affectivity,
        pid5_detachment,
        pid5_antagonism,
        pid5_disinhibition,
        pid5_psychoticism
      ),
      list(M = ~ mean(.x, na.rm = TRUE), SD = ~ sd(.x, na.rm = TRUE))
    )
  ) %>%
  pivot_longer(everything(), names_to = "stat", values_to = "value")

print(pid5_descriptives)

# ==============================================================================
# 9. DESCRIPTIVE STATISTICS: Acoustic by Timepoint
# ==============================================================================

cat("\n=== DESCRIPTIVE STATISTICS: ACOUSTIC MEASURES BY TIMEPOINT ===\n")

descriptives <- df_analysis %>%
  group_by(timepoint) %>%
  summarise(
    across(
      c(f0_mean_a, f0_std_a, jitter_a, nne_a, f2_mean_a, f2_std_a),
      list(M = ~ mean(., na.rm = TRUE), SD = ~ sd(., na.rm = TRUE))
    ),
    n = n()
  )
print(descriptives, width = Inf)

# ==============================================================================
# 10. FIX LOGNORMAL ISSUES (valori <= 0)
# ==============================================================================

lognormal_vars <- c("f0_std", "f2_std", "jitter")
vowels <- c("a", "i", "u")

for (var in lognormal_vars) {
  for (v in vowels) {
    col <- paste0(var, "_", v)
    if (!col %in% names(df_analysis)) next

    bad <- which(df_analysis[[col]] <= 0)
    if (length(bad) > 0) {
      min_pos <- min(df_analysis[[col]][df_analysis[[col]] > 0], na.rm = TRUE)
      df_analysis[[col]][bad] <- min_pos
      cat("Fixed", length(bad), "values <=0 in", col, "\n")
    }
  }
}

# ==============================================================================
# 11. BAYESIAN MODELS - MAIN EFFECTS (vowel /a/)
# ==============================================================================

cat("\n=== FITTING MAIN EFFECTS MODELS ===\n")

# Priors
prior_f0mean <- c(
  prior(student_t(3, 220, 30), class = Intercept),
  prior(normal(0, 10), class = b),
  prior(exponential(0.3), class = sigma),
  prior(exponential(0.3), class = sd)
)

prior_lognormal <- c(
  prior(student_t(3, 0, 1), class = Intercept),
  prior(normal(0, 1), class = b),
  prior(exponential(1), class = sigma),
  prior(exponential(0.5), class = sd)
)

prior_nne <- c(
  prior(student_t(3, -20, 5), class = Intercept),
  prior(normal(0, 5), class = b),
  prior(exponential(0.5), class = sigma),
  prior(exponential(0.3), class = sd)
)

prior_f2mean <- c(
  prior(student_t(3, 1500, 150), class = Intercept),
  prior(normal(0, 50), class = b),
  prior(exponential(0.01), class = sigma),
  prior(exponential(0.01), class = sd),
  prior(gamma(2, 0.1), class = nu)
)

# Model fitting settings
iter <- 6000
warmup <- 3000
chains <- 4
cores <- 4
seed <- 123

# M1: F0 Mean
cat("\nFitting m1_f0mean_a...\n")
m1_f0mean_a <- brm(
  f0_mean_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = gaussian(),
  prior = prior_f0mean,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
  file = "models/m1_f0mean_a"
)

# M2: F0 SD
cat("Fitting m2_f0std_a...\n")
m2_f0std_a <- brm(
  f0_std_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = lognormal(),
  prior = prior_lognormal,
  iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
  file = "models/m2_f0std_a"
)

# M3: Jitter
cat("Fitting m3_jitter_a...\n")
m3_jitter_a <- brm(
  jitter_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = lognormal(),
  prior = prior_lognormal,
  control = list(adapt_delta = 0.98, max_treedepth = 15),
  iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
  file = "models/m3_jitter_a"
)

# M4: NNE
cat("Fitting m4_nne_a...\n")
m4_nne_a <- brm(
  nne_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = gaussian(),
  prior = prior_nne,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
  file = "models/m4_nne_a"
)

# M5: F2 Mean
cat("Fitting m5_f2mean_a...\n")
m5_f2mean_a <- brm(
  f2_mean_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = student(),
  prior = prior_f2mean,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
  file = "models/m5_f2mean_a"
)

# M6: F2 SD
cat("Fitting m6_f2std_a...\n")
m6_f2std_a <- brm(
  f2_std_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = lognormal(),
  prior = prior_lognormal,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
  file = "models/m6_f2std_a"
)

# ==============================================================================
# 12. MAIN EFFECTS RESULTS
# ==============================================================================

cat("\n=== MAIN EFFECTS RESULTS ===\n")

summarize_main_effects <- function(model, outcome_name) {
  cat("\n", rep("-", 60), "\n", sep = "")
  cat(outcome_name, "\n")
  cat(rep("-", 60), "\n", sep = "")

  sum_mod <- summary(model)
  print(round(sum_mod$fixed, 3))

  cat("\nPlanned Contrasts:\n")
  cat("[H1] PRE vs BASELINE:\n")
  h1 <- hypothesis(model, "timepointpre = 0")
  cat("  β =", round(h1$hypothesis$Estimate, 3),
      ", 95% CI [", round(h1$hypothesis$CI.Lower, 3), ",",
      round(h1$hypothesis$CI.Upper, 3), "]\n")

  cat("[H2] POST vs PRE:\n")
  h2 <- hypothesis(model, "timepointpost - timepointpre = 0")
  cat("  β =", round(h2$hypothesis$Estimate, 3),
      ", 95% CI [", round(h2$hypothesis$CI.Lower, 3), ",",
      round(h2$hypothesis$CI.Upper, 3), "]\n")

  cat("[H3] POST vs BASELINE:\n")
  h3 <- hypothesis(model, "timepointpost = 0")
  cat("  β =", round(h3$hypothesis$Estimate, 3),
      ", 95% CI [", round(h3$hypothesis$CI.Lower, 3), ",",
      round(h3$hypothesis$CI.Upper, 3), "]\n")
}

summarize_main_effects(m1_f0mean_a, "F0 Mean (Hz)")
summarize_main_effects(m2_f0std_a, "F0 SD [lognormal]")
summarize_main_effects(m3_jitter_a, "Jitter [lognormal]")
summarize_main_effects(m4_nne_a, "NNE (dB)")
summarize_main_effects(m5_f2mean_a, "F2 Mean (Hz)")
summarize_main_effects(m6_f2std_a, "F2 SD [lognormal]")

# ==============================================================================
# 13. SAVE
# ==============================================================================

saveRDS(df_analysis, "results/df_analysis.rds")
saveRDS(validation_cors, "results/ema_validation_correlations.rds")
save.image("results/main_effects_workspace.RData")

cat("\n", rep("=", 70), "\n", sep = "")
cat("MAIN EFFECTS ANALYSIS COMPLETE\n")
cat("PID-5 version used:", PID5_VERSION, "\n")
cat("Sample: Female students only\n")
cat("EMA validation completed: see results/ema_validation_correlations.rds\n")
cat(rep("=", 70), "\n\n")
cat("Proceed to 03_moderation_analysis_FINAL.R...\n\n")
