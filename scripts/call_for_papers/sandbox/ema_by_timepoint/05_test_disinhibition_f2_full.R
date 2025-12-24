# ==============================================================================
# TEST: Disinhibition × F2 con gestione outlier
# ==============================================================================
# F2 mean presenta outlier (valori > 2000 Hz probabilmente errori di estrazione)
# Questo script offre diverse strategie per gestirli
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(rio)
  library(here)
  library(bayestestR)
  library(missRanger)
})

options(brms.backend = "cmdstanr")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("TEST: Disinhibition × F2 con gestione outlier\n")
cat(rep("=", 70), "\n\n", sep = "")

# ==============================================================================
# CONFIGURAZIONE: Scegli metodo per gestire outlier
# ==============================================================================

# Opzioni:
# "winsorize"  = Sostituisce outlier con percentile (RACCOMANDATO - no rumore aggiunto)
# "impute"     = Mette NA e imputa con missRanger (proposta originale)
# "mad"        = Usa MAD per identificare outlier, poi winsorizza
# "exclude"    = Esclude completamente le osservazioni con outlier

OUTLIER_METHOD <- "winsorize" # ← MODIFICA QUI SE NECESSARIO

# Soglie
F2_UPPER_THRESHOLD <- 1800 # Valori > 1800 considerati outlier
F2_LOWER_THRESHOLD <- 600 # Valori < 600 potrebbero essere errori (opzionale)
WINSORIZE_PERCENTILE <- 0.975 # Per winsorizzazione

cat("Metodo outlier selezionato:", OUTLIER_METHOD, "\n\n")

# ==============================================================================
# 1. CARICA DATI
# ==============================================================================

if (file.exists("results/df_analysis.rds")) {
  df_orig <- readRDS("results/df_analysis.rds")
  cat("✓ Loaded df_analysis\n")
} else {
  stop("Esegui prima 02_voice_personality_analysis_CORRECTED.R")
}

# ==============================================================================
# 2. DIAGNOSTICA OUTLIER F2
# ==============================================================================

cat("\n=== DIAGNOSTICA OUTLIER F2 ===\n")

f2_vars <- c(
  "f2_mean_a",
  "f2_mean_i",
  "f2_mean_u",
  "f2_std_a",
  "f2_std_i",
  "f2_std_u"
)

for (v in f2_vars) {
  if (!v %in% names(df_orig)) next

  x <- df_orig[[v]]
  x <- x[!is.na(x)]

  # Statistiche
  cat("\n", v, ":\n", sep = "")
  cat("  N =", length(x), "\n")
  cat("  Mean =", round(mean(x), 1), ", SD =", round(sd(x), 1), "\n")
  cat("  Median =", round(median(x), 1), ", MAD =", round(mad(x), 1), "\n")
  cat("  Range: [", round(min(x), 1), ", ", round(max(x), 1), "]\n", sep = "")
  cat(
    "  Q1 =",
    round(quantile(x, 0.25), 1),
    ", Q3 =",
    round(quantile(x, 0.75), 1),
    "\n"
  )
  cat(
    "  P2.5 =",
    round(quantile(x, 0.025), 1),
    ", P97.5 =",
    round(quantile(x, 0.975), 1),
    "\n"
  )

  # Conta outlier con soglia fissa
  n_high <- sum(x > F2_UPPER_THRESHOLD)
  n_low <- sum(x < F2_LOWER_THRESHOLD)
  cat(
    "  Outlier (>",
    F2_UPPER_THRESHOLD,
    "):",
    n_high,
    "(",
    round(100 * n_high / length(x), 1),
    "%)\n"
  )
  cat(
    "  Outlier (<",
    F2_LOWER_THRESHOLD,
    "):",
    n_low,
    "(",
    round(100 * n_low / length(x), 1),
    "%)\n"
  )

  # MAD-based outlier detection
  med <- median(x)
  mad_val <- mad(x)
  n_mad_outlier <- sum(abs(x - med) > 3 * mad_val)
  cat("  MAD outlier (|x - median| > 3*MAD):", n_mad_outlier, "\n")
}

# ==============================================================================
# 3. FUNZIONI PER GESTIONE OUTLIER
# ==============================================================================

# Winsorizzazione
winsorize <- function(x, probs = c(0.025, 0.975)) {
  if (all(is.na(x))) return(x)
  q <- quantile(x, probs, na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

# Winsorizzazione con soglia fissa
winsorize_fixed <- function(x, lower = NULL, upper = NULL) {
  if (!is.null(lower)) x[x < lower] <- lower
  if (!is.null(upper)) x[x > upper] <- upper
  x
}

# MAD-based winsorizzazione
winsorize_mad <- function(x, k = 3) {
  if (all(is.na(x))) return(x)
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  lower <- med - k * mad_val
  upper <- med + k * mad_val
  x[x < lower] <- lower
  x[x > upper] <- upper
  x
}

# ==============================================================================
# 4. APPLICA METODO SCELTO
# ==============================================================================

df_clean <- df_orig

cat("\n=== Applicazione metodo:", OUTLIER_METHOD, "===\n")

if (OUTLIER_METHOD == "winsorize") {
  # Winsorizzazione con percentili
  for (v in f2_vars) {
    if (!v %in% names(df_clean)) next
    before <- sum(df_clean[[v]] > F2_UPPER_THRESHOLD, na.rm = TRUE)
    df_clean[[v]] <- winsorize(
      df_clean[[v]],
      probs = c(0.025, WINSORIZE_PERCENTILE)
    )
    after <- sum(df_clean[[v]] > F2_UPPER_THRESHOLD, na.rm = TRUE)
    cat(v, ": winsorizzati", before, "valori alti\n")
  }
} else if (OUTLIER_METHOD == "impute") {
  # NA + imputazione
  for (v in f2_vars) {
    if (!v %in% names(df_clean)) next
    n_outlier <- sum(
      df_clean[[v]] > F2_UPPER_THRESHOLD |
        df_clean[[v]] < F2_LOWER_THRESHOLD,
      na.rm = TRUE
    )
    df_clean[[v]][df_clean[[v]] > F2_UPPER_THRESHOLD] <- NA
    df_clean[[v]][df_clean[[v]] < F2_LOWER_THRESHOLD] <- NA
    cat(v, ": impostati", n_outlier, "valori a NA\n")
  }

  # Converti tutte le colonne numeriche a double per missRanger
  numeric_cols <- names(df_clean)[sapply(df_clean, is.numeric)]
  for (col in numeric_cols) {
    df_clean[[col]] <- as.numeric(df_clean[[col]])
  }

  cat("\nImputazione con missRanger...\n")
  set.seed(123)
  df_clean <- missRanger(df_clean, num.trees = 200, verbose = 1)
} else if (OUTLIER_METHOD == "mad") {
  # MAD-based winsorizzazione
  for (v in f2_vars) {
    if (!v %in% names(df_clean)) next
    x_orig <- df_clean[[v]]
    df_clean[[v]] <- winsorize_mad(df_clean[[v]], k = 3)
    n_changed <- sum(x_orig != df_clean[[v]], na.rm = TRUE)
    cat(v, ": winsorizzati", n_changed, "valori (MAD-based)\n")
  }
} else if (OUTLIER_METHOD == "exclude") {
  # Esclusione completa
  for (v in f2_vars) {
    if (!v %in% names(df_clean)) next
    n_outlier <- sum(
      df_clean[[v]] > F2_UPPER_THRESHOLD |
        df_clean[[v]] < F2_LOWER_THRESHOLD,
      na.rm = TRUE
    )
    cat(v, ":", n_outlier, "osservazioni da escludere\n")
  }

  # Crea flag per esclusione
  df_clean <- df_clean %>%
    mutate(
      f2_outlier = (f2_mean_a > F2_UPPER_THRESHOLD |
        f2_mean_a < F2_LOWER_THRESHOLD) |
        (f2_mean_i > F2_UPPER_THRESHOLD | f2_mean_i < F2_LOWER_THRESHOLD) |
        (f2_mean_u > F2_UPPER_THRESHOLD | f2_mean_u < F2_LOWER_THRESHOLD)
    ) %>%
    filter(!f2_outlier | is.na(f2_outlier)) %>%
    select(-f2_outlier)

  cat("Osservazioni rimanenti:", nrow(df_clean), "\n")
}

# ==============================================================================
# 5. VERIFICA POST-PULIZIA
# ==============================================================================

cat("\n=== Statistiche post-pulizia ===\n")

for (v in c("f2_mean_a", "f2_mean_i", "f2_mean_u")) {
  if (!v %in% names(df_clean)) next
  x <- df_clean[[v]]
  x <- x[!is.na(x)]
  cat(
    v,
    ": Mean =",
    round(mean(x), 1),
    ", SD =",
    round(sd(x), 1),
    ", Range = [",
    round(min(x), 1),
    ",",
    round(max(x), 1),
    "]\n"
  )
}

# ==============================================================================
# 6. PREPARA PREDITTORI
# ==============================================================================

# Aggiungi variabili PID-5 full baseline se non presenti
baseline_vars <- c(
  "domain_negative_affect_baseline",
  "domain_detachment_baseline",
  "domain_antagonism_baseline",
  "domain_disinhibition_baseline",
  "domain_psychoticism_baseline"
)

if (!all(baseline_vars %in% names(df_clean))) {
  cat("\nAggiunta variabili PID-5 baseline...\n")
  data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
  d <- rio::import(data_path)

  baseline_scores <- d %>%
    dplyr::select(user_id, all_of(baseline_vars)) %>%
    dplyr::rename(ID = user_id) %>%
    group_by(ID) %>%
    summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop")

  df_clean <- df_clean %>%
    left_join(baseline_scores, by = "ID")
}

# Centra e standardizza
df_test <- df_clean %>%
  mutate(
    # PID-5 Full
    pid5_na_full_c = scale(domain_negative_affect_baseline)[, 1],
    pid5_det_full_c = scale(domain_detachment_baseline)[, 1],
    pid5_ant_full_c = scale(domain_antagonism_baseline)[, 1],
    pid5_dis_full_c = scale(domain_disinhibition_baseline)[, 1],
    pid5_psy_full_c = scale(domain_psychoticism_baseline)[, 1],

    # PID-5 EMA (se presenti)
    pid5_na_ema_c = if ("pid5_negative_affectivity_c" %in% names(.))
      pid5_negative_affectivity_c else scale(pid5_negative_affectivity)[, 1],
    pid5_det_ema_c = if ("pid5_detachment_c" %in% names(.))
      pid5_detachment_c else scale(pid5_detachment)[, 1],
    pid5_ant_ema_c = if ("pid5_antagonism_c" %in% names(.))
      pid5_antagonism_c else scale(pid5_antagonism)[, 1],
    pid5_dis_ema_c = if ("pid5_disinhibition_c" %in% names(.))
      pid5_disinhibition_c else scale(pid5_disinhibition)[, 1],
    pid5_psy_ema_c = if ("pid5_psychoticism_c" %in% names(.))
      pid5_psychoticism_c else scale(pid5_psychoticism)[, 1],

    # Contrast coding
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre" ~ 0.5,
      timepoint == "post" ~ 0
    ),
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre" ~ -0.5,
      timepoint == "post" ~ 0.5
    )
  ) %>%
  filter(!is.na(pid5_dis_full_c))

cat("\nN osservazioni per analisi:", nrow(df_test), "\n")
cat("N soggetti:", n_distinct(df_test$ID), "\n")

# ==============================================================================
# 7. MODELLI F2 MEAN - CONFRONTO EMA vs FULL
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODELLI F2 MEAN - Confronto EMA vs Full PID-5\n")
cat(rep("=", 70), "\n", sep = "")

# Priors per F2
prior_f2 <- c(
  prior(student_t(3, 1200, 150), class = Intercept),
  prior(normal(0, 50), class = b),
  prior(exponential(0.01), class = sigma),
  prior(exponential(0.01), class = sd),
  prior(gamma(2, 0.1), class = nu)
)

# Settings
iter <- 5000
warmup <- 2500
chains <- 4
cores <- 4
seed <- 123
control <- list(adapt_delta = 0.995, max_treedepth = 18)

results_comparison <- list()

for (v in c("a", "i", "u")) {
  col <- paste0("f2_mean_", v)
  if (!col %in% names(df_test)) next

  cat("\n--- Vocale /", v, "/ ---\n", sep = "")

  # Modello con PID-5 Full
  cat("Fitting modello PID-5 Full...\n")

  traits_full <- "pid5_na_full_c + pid5_det_full_c + pid5_ant_full_c + pid5_dis_full_c + pid5_psy_full_c"

  formula_full <- bf(as.formula(paste0(
    col,
    " ~ c1_stress * (",
    traits_full,
    ") + c2_recovery * (",
    traits_full,
    ") + (1 + c1_stress + c2_recovery || ID)"
  )))

  m_full <- brm(
    formula = formula_full,
    data = df_test,
    family = student(),
    prior = prior_f2,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    control = control,
    file = paste0("models/m_f2_full_clean_", v)
  )

  # Modello con PID-5 EMA
  cat("Fitting modello PID-5 EMA...\n")

  traits_ema <- "pid5_na_ema_c + pid5_det_ema_c + pid5_ant_ema_c + pid5_dis_ema_c + pid5_psy_ema_c"

  formula_ema <- bf(as.formula(paste0(
    col,
    " ~ c1_stress * (",
    traits_ema,
    ") + c2_recovery * (",
    traits_ema,
    ") + (1 + c1_stress + c2_recovery || ID)"
  )))

  m_ema <- brm(
    formula = formula_ema,
    data = df_test,
    family = student(),
    prior = prior_f2,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    control = control,
    file = paste0("models/m_f2_ema_clean_", v)
  )

  # Estrai risultati Disinhibition
  fe_full <- fixef(m_full)
  fe_ema <- fixef(m_ema)

  dis_full <- fe_full[grep("dis_full", rownames(fe_full)), , drop = FALSE]
  dis_ema <- fe_ema[grep("dis_ema", rownames(fe_ema)), , drop = FALSE]

  results_comparison[[paste0("full_", v)]] <- dis_full
  results_comparison[[paste0("ema_", v)]] <- dis_ema

  cat("\nDISINHIBITION - PID-5 Full:\n")
  print(round(dis_full, 3))

  cat("\nDISINHIBITION - PID-5 EMA:\n")
  print(round(dis_ema, 3))
}

# ==============================================================================
# 8. RIEPILOGO FINALE
# ==============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("RIEPILOGO: Disinhibition × F2 (dati puliti)\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nAbstract riporta: β = 63 Hz [26, 99] per Disinhibition × F2\n\n")

cat("--- Risultati per vocale /a/ ---\n")
if (!is.null(results_comparison$full_a)) {
  cat("\nPID-5 Full:\n")
  print(round(results_comparison$full_a, 2))
}
if (!is.null(results_comparison$ema_a)) {
  cat("\nPID-5 EMA:\n")
  print(round(results_comparison$ema_a, 2))
}

cat("\n--- Risultati per vocale /i/ ---\n")
if (!is.null(results_comparison$full_i)) {
  cat("\nPID-5 Full:\n")
  print(round(results_comparison$full_i, 2))
}

cat("\n--- Risultati per vocale /u/ ---\n")
if (!is.null(results_comparison$full_u)) {
  cat("\nPID-5 Full:\n")
  print(round(results_comparison$full_u, 2))
}

# Estrai effetto specifico Disinhibition × Stress per confronto
cat("\n\n--- CONFRONTO DIRETTO: c1_stress × Disinhibition ---\n")
cat("(Questo è probabilmente l'effetto riportato nell'abstract)\n\n")

for (v in c("a", "i", "u")) {
  full_key <- paste0("full_", v)
  ema_key <- paste0("ema_", v)

  if (!is.null(results_comparison[[full_key]])) {
    fe <- results_comparison[[full_key]]
    stress_row <- grep("c1_stress", rownames(fe))
    if (length(stress_row) > 0) {
      cat(
        "Vocale /",
        v,
        "/ - PID-5 Full: β = ",
        round(fe[stress_row, "Estimate"], 1),
        " [",
        round(fe[stress_row, "Q2.5"], 1),
        ", ",
        round(fe[stress_row, "Q97.5"], 1),
        "]\n",
        sep = ""
      )
    }
  }

  if (!is.null(results_comparison[[ema_key]])) {
    fe <- results_comparison[[ema_key]]
    stress_row <- grep("c1_stress", rownames(fe))
    if (length(stress_row) > 0) {
      cat(
        "Vocale /",
        v,
        "/ - PID-5 EMA:  β = ",
        round(fe[stress_row, "Estimate"], 1),
        " [",
        round(fe[stress_row, "Q2.5"], 1),
        ", ",
        round(fe[stress_row, "Q97.5"], 1),
        "]\n",
        sep = ""
      )
    }
  }
  cat("\n")
}

cat("\n")
cat("INTERPRETAZIONE:\n")
cat("- Abstract: β = 63 Hz [26, 99]\n")
cat("- Se i risultati sono vicini → effetto replicato\n")
cat("- Se molto diversi → effetto specifico per campione originale\n")
cat(
  "- La pulizia outlier (metodo:",
  OUTLIER_METHOD,
  ") potrebbe influenzare i risultati\n"
)

# Salva
saveRDS(results_comparison, "results/disinhibition_f2_comparison.rds")
saveRDS(df_test, "results/df_test_f2_clean.rds")

cat("\n✓ Analisi completata. Risultati salvati.\n\n")
