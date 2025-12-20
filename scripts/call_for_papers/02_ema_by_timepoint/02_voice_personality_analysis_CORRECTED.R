# ==============================================================================
# Bayesian Analysis: Voice Acoustics, Stress, and Personality Pathology
# CORRECTED VERSION - Proper EMA aggregation
# ==============================================================================
# CORREZIONI PRINCIPALI:
#   1. PID-5 EMA: aggregati su TUTTI i timepoint EMA (non solo baseline)
#   2. Opzione per usare domain_*_baseline come alternativa/aggiunta
#   3. Versione "time-varying" opzionale per analisi di sensibilità
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
# CONFIGURAZIONE: Scegli quale versione dei PID-5 usare
# ==============================================================================

# Opzioni:
# "ema_all"      = Media di tutti i valori EMA per soggetto (RACCOMANDATO per abstract)
# "ema_baseline" = Solo valori EMA del periodo baseline (versione attuale - ERRATA)
# "pid5_full"    = Punteggi PID-5 completo del baseline (domain_*_baseline)
# "ema_timevarying" = Valori EMA che variano per timepoint (più rumore, meno stabilità)

PID5_VERSION <- "ema_all" # ← MODIFICA QUI SE NECESSARIO

cat("\n=== CONFIGURAZIONE ===\n")
cat("Versione PID-5 selezionata:", PID5_VERSION, "\n")

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

cat("\nVoice data: N soggetti unici =", n_distinct(df_wide$ID), "\n")

df_clean <- df_wide %>%
  dplyr::select(
    ID,
    Data,
    timepoint,
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
    user_id,
    exam_period,
    day,
    sex,
    all_of(ema_pid5_vars),
    all_of(baseline_pid5_vars)
  ) %>%
  dplyr::rename(ID = user_id)

# Filtra solo soggetti con dati vocali
allowed_id <- unique(df_clean$ID)
d1 <- d1 %>% filter(ID %in% allowed_id)

# ==============================================================================
# 3. AGGREGAZIONE PID-5 - VERSIONI DIVERSE
# ==============================================================================

cat("\n=== Preparazione PID-5 con versione:", PID5_VERSION, "===\n")

if (PID5_VERSION == "ema_all") {
  # =========================================================================
  # VERSIONE RACCOMANDATA: Media di TUTTI i valori EMA per soggetto
  # Questo crea una misura trait-like stabile basata su tutti i 2.5 mesi di EMA
  # =========================================================================

  pid5_person <- d1 %>%
    group_by(ID) %>%
    summarise(
      across(all_of(ema_pid5_vars), ~ mean(.x, na.rm = TRUE)),
      sex = first(na.omit(sex)),
      # Mantieni anche i baseline PID-5 per analisi supplementari
      across(all_of(baseline_pid5_vars), ~ first(na.omit(.x))),
      n_ema_obs = n(),
      .groups = "drop"
    )

  cat("N soggetti con EMA:", nrow(pid5_person), "\n")
  cat("Media osservazioni EMA per soggetto:", mean(pid5_person$n_ema_obs), "\n")

  # Join: ogni timepoint vocale ha gli STESSI valori PID-5 (between-person)
  voice_df <- df_clean
  joined <- voice_df %>%
    left_join(pid5_person, by = "ID")
} else if (PID5_VERSION == "ema_baseline") {
  # =========================================================================
  # VERSIONE ORIGINALE (PROBLEMATICA): Solo EMA del periodo baseline
  # =========================================================================

  d1$timepoint <- forcats::fct_recode(
    d1$exam_period,
    "baseline" = "baseline",
    "pre" = "pre_exam",
    "post" = "post_exam"
  )

  pid5_agg <- d1 %>%
    group_by(ID, timepoint) %>%
    summarise(
      across(all_of(ema_pid5_vars), ~ mean(.x, na.rm = TRUE)),
      sex = first(na.omit(sex)),
      .groups = "drop"
    )

  # Estrai solo baseline
  pid5_baseline <- pid5_agg %>%
    filter(timepoint == "baseline") %>%
    dplyr::select(-timepoint)

  voice_df <- df_clean
  joined <- voice_df %>%
    left_join(pid5_baseline, by = "ID")
} else if (PID5_VERSION == "pid5_full") {
  # =========================================================================
  # VERSIONE ALTERNATIVA: PID-5 completo (220 item) del baseline
  # =========================================================================

  pid5_person <- d1 %>%
    group_by(ID) %>%
    summarise(
      # Rinomina per coerenza con lo script
      pid5_negative_affectivity = first(na.omit(
        domain_negative_affect_baseline
      )),
      pid5_detachment = first(na.omit(domain_detachment_baseline)),
      pid5_antagonism = first(na.omit(domain_antagonism_baseline)),
      pid5_disinhibition = first(na.omit(domain_disinhibition_baseline)),
      pid5_psychoticism = first(na.omit(domain_psychoticism_baseline)),
      sex = first(na.omit(sex)),
      .groups = "drop"
    )

  voice_df <- df_clean
  joined <- voice_df %>%
    left_join(pid5_person, by = "ID")
} else if (PID5_VERSION == "ema_timevarying") {
  # =========================================================================
  # VERSIONE TIME-VARYING: Valori EMA diversi per ogni timepoint vocale
  # ATTENZIONE: Più rumore, meno stabilità - usa con cautela
  # =========================================================================

  d1$timepoint <- forcats::fct_recode(
    d1$exam_period,
    "baseline" = "baseline",
    "pre" = "pre_exam",
    "post" = "post_exam"
  )

  pid5_timevarying <- d1 %>%
    group_by(ID, timepoint) %>%
    summarise(
      across(all_of(ema_pid5_vars), ~ mean(.x, na.rm = TRUE)),
      sex = first(na.omit(sex)),
      .groups = "drop"
    )

  voice_df <- df_clean
  joined <- voice_df %>%
    left_join(pid5_timevarying, by = c("ID", "timepoint"))
} else {
  stop(
    "PID5_VERSION non riconosciuta. Usa: 'ema_all', 'ema_baseline', 'pid5_full', o 'ema_timevarying'"
  )
}

# ==============================================================================
# 4. IMPUTATION E PREPARAZIONE FINALE
# ==============================================================================

cat("\nMissing prima dell'imputazione:\n")
print(colSums(is.na(joined)))

# Fix: n_ema_obs è integer ma missRanger imputa come double
# Convertiamo a numeric prima dell'imputazione
if ("n_ema_obs" %in% names(joined)) {
  joined$n_ema_obs <- as.numeric(joined$n_ema_obs)
}

# Anche domain_*_baseline potrebbero avere lo stesso problema
# Assicuriamoci che tutte le variabili numeriche siano double
numeric_cols <- names(joined)[sapply(joined, is.numeric)]
for (col in numeric_cols) {
  if (is.integer(joined[[col]])) {
    joined[[col]] <- as.numeric(joined[[col]])
  }
}

set.seed(123)
imp <- missRanger(joined, num.trees = 200, verbose = 1)
df_clean <- imp

# Standardizza PID-5 (between-person, centrati)
df_long <- df_clean %>%
  mutate(
    pid5_negative_affectivity_c = scale(pid5_negative_affectivity)[, 1],
    pid5_detachment_c = scale(pid5_detachment)[, 1],
    pid5_antagonism_c = scale(pid5_antagonism)[, 1],
    pid5_disinhibition_c = scale(pid5_disinhibition)[, 1],
    pid5_psychoticism_c = scale(pid5_psychoticism)[, 1]
  )

df_analysis <- df_long %>%
  dplyr::filter(complete.cases(.))

# ==============================================================================
# 5. SAMPLE CHARACTERISTICS
# ==============================================================================

cat("\n=== SAMPLE CHARACTERISTICS ===\n")
cat("Total observations:", nrow(df_analysis), "\n")
cat("Unique participants:", n_distinct(df_analysis$ID), "\n")
df_analysis %>% count(timepoint) %>% print()

cat("\nPID-5 EMA descriptives:\n")
df_analysis %>%
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
  pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
  print()

# ==============================================================================
# 6. DESCRIPTIVE STATISTICS: Acoustic by Timepoint
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
# 7. FIX LOGNORMAL ISSUES (valori <= 0)
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
# 8. BAYESIAN MODELS - MAIN EFFECTS (vowel /a/)
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

# M1: F0 Mean (gaussian)
cat("\nFitting m1_f0mean_a...\n")
m1_f0mean_a <- brm(
  f0_mean_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = gaussian(),
  prior = prior_f0mean,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  seed = seed,
  file = "models/m1_f0mean_a"
)
pp_check(m1_f0mean_a)
summary(m1_f0mean_a)

# M2: F0 SD (lognormal)
cat("Fitting m2_f0std_a...\n")
m2_f0std_a <- brm(
  f0_std_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = lognormal(),
  prior = prior_lognormal,
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  seed = seed,
  file = "models/m2_f0std_a"
)
pp_check(m2_f0std_a)
summary(m2_f0std_a)

# M3: Jitter (lognormal)
cat("Fitting m3_jitter_a...\n")
m3_jitter_a <- brm(
  jitter_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = lognormal(),
  prior = prior_lognormal,
  control = list(adapt_delta = 0.98, max_treedepth = 15),
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  seed = seed,
  file = "models/m3_jitter_a"
)
pp_check(m3_jitter_a)
summary(m3_jitter_a)

# M4: NNE (gaussian)
cat("Fitting m4_nne_a...\n")
m4_nne_a <- brm(
  nne_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = gaussian(),
  prior = prior_nne,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  seed = seed,
  file = "models/m4_nne_a"
)
pp_check(m4_nne_a)
summary(m4_nne_a)

# M5: F2 Mean (student-t)
cat("Fitting m5_f2mean_a...\n")
m5_f2mean_a <- brm(
  f2_mean_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = student(),
  prior = prior_f2mean,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  seed = seed,
  file = "models/m5_f2mean_a"
)
pp_check(m5_f2mean_a)
summary(m5_f2mean_a)

# M6: F2 SD (lognormal)
cat("Fitting m6_f2std_a...\n")
m6_f2std_a <- brm(
  f2_std_a ~ timepoint + (1 + timepoint || ID),
  data = df_analysis,
  family = lognormal(),
  prior = prior_lognormal,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  seed = seed,
  file = "models/m6_f2std_a"
)
pp_check(m6_f2std_a)
summary(m6_f2std_a)

# ==============================================================================
# 9. MAIN EFFECTS RESULTS
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
  cat(
    "  β =",
    round(h1$hypothesis$Estimate, 3),
    ", 95% CI [",
    round(h1$hypothesis$CI.Lower, 3),
    ",",
    round(h1$hypothesis$CI.Upper, 3),
    "]\n"
  )

  cat("[H2] POST vs PRE:\n")
  h2 <- hypothesis(model, "timepointpost - timepointpre = 0")
  cat(
    "  β =",
    round(h2$hypothesis$Estimate, 3),
    ", 95% CI [",
    round(h2$hypothesis$CI.Lower, 3),
    ",",
    round(h2$hypothesis$CI.Upper, 3),
    "]\n"
  )

  cat("[H3] POST vs BASELINE:\n")
  h3 <- hypothesis(model, "timepointpost = 0")
  cat(
    "  β =",
    round(h3$hypothesis$Estimate, 3),
    ", 95% CI [",
    round(h3$hypothesis$CI.Lower, 3),
    ",",
    round(h3$hypothesis$CI.Upper, 3),
    "]\n"
  )
}

summarize_main_effects(m1_f0mean_a, "F0 Mean (Hz)")
summarize_main_effects(m2_f0std_a, "F0 SD [lognormal, log-scale]")
summarize_main_effects(m3_jitter_a, "Jitter [lognormal, log-scale]")
summarize_main_effects(m4_nne_a, "NNE (dB)")
summarize_main_effects(m5_f2mean_a, "F2 Mean (Hz) [student-t]")
summarize_main_effects(m6_f2std_a, "F2 SD [lognormal, log-scale]")

# ==============================================================================
# 10. SAVE
# ==============================================================================

saveRDS(df_analysis, "results/df_analysis.rds")
save.image("results/main_effects_workspace.RData")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MAIN EFFECTS ANALYSIS COMPLETE\n")
cat("PID-5 version used:", PID5_VERSION, "\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nProceed to 03_moderation_analysis_CORRECTED.R...\n\n")
