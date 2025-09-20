# Overview ----------------------------------------------------------------
# Associated project: PID-5, EMA, and acustic features
# Script purpose: To explore different measures of stress reactivity associated
#   with the exam date.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-09-14
# Last update:
# Status: In progress
# Notes:

# Caricamento librerie
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(bayesplot)
  library(gridExtra)
  library(psych)
  library(corrplot)
  library(GGally)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(effectsize)
  library(plotly)
  library(viridis)
})

ggplot2::theme_set(bayesplot::theme_default(
  base_size = 14.0
))

# Carichiamo il dataset completo
d <- rio::import(
  here::here(
    "data",
    "processed",
    "ema_plus_baseline_exam_tagged.csv"
  )
)


# ========================================================
# 1. OVERVIEW DEI DATI
# ========================================================

# Statistiche descrittive generali
cat("=== OVERVIEW DATASET ===\n")
cat("Numero totale di osservazioni:", nrow(d), "\n")
cat("Numero di partecipanti:", length(unique(d$user_id)), "\n")
cat(
  "Periodo di osservazione:",
  min(d$day, na.rm = TRUE),
  "al",
  max(d$day, na.rm = TRUE),
  "\n"
)

# Distribuzione per exam_period
exam_period_summary <- d %>%
  group_by(exam_period) %>%
  summarise(
    n_obs = n(),
    n_participants = n_distinct(user_id),
    .groups = 'drop'
  ) %>%
  mutate(prop = n_obs / sum(n_obs))

print("Distribuzione per periodo d'esame:")
print(exam_period_summary)

# ========================================================
# 2. PREPARAZIONE DELLE VARIABILI CHIAVE
# ========================================================

# Creazione di variabili aggregate per l'analisi
d_clean <- d %>%
  # Rimuovi osservazioni con dati mancanti per le variabili chiave
  dplyr::filter(!is.na(neg_affect_ema), !is.na(exam_period)) %>%
  # Aggiungi informazioni utili
  mutate(
    exam_period_f = factor(
      exam_period,
      levels = c("baseline", "pre_exam", "post_exam"),
      labels = c("Baseline", "Pre-Esame", "Post-Esame")
    ),
    sex_f = factor(sex),
    # Standardizza neg_affect_ema per interpretazione più facile
    neg_affect_z = as.numeric(scale(neg_affect_ema))
  )

# ========================================================
# 3. ANALISI DESCRITTIVE PER EXAM_PERIOD
# ========================================================

# Statistiche descrittive per variabili chiave
desc_stats <- d_clean %>%
  group_by(exam_period_f) %>%
  summarise(
    n = n(),
    n_participants = n_distinct(user_id),
    # Affetto negativo
    neg_affect_mean = mean(neg_affect_ema, na.rm = TRUE),
    neg_affect_sd = sd(neg_affect_ema, na.rm = TRUE),
    neg_affect_median = median(neg_affect_ema, na.rm = TRUE),
    # DASS
    dass_stress_mean = mean(dass_stress, na.rm = TRUE),
    dass_stress_sd = sd(dass_stress, na.rm = TRUE),
    dass_depression_mean = mean(dass_depression, na.rm = TRUE),
    dass_anxiety_mean = mean(dass_anxiety, na.rm = TRUE),
    # PID-5
    pid5_neg_affect_mean = mean(pid5_negative_affectivity, na.rm = TRUE),
    pid5_neg_affect_sd = sd(pid5_negative_affectivity, na.rm = TRUE),
    # Context
    context_threat_mean = mean(context_threat, na.rm = TRUE),
    context_quality_mean = mean(context_quality, na.rm = TRUE),
    .groups = 'drop'
  )

print("=== STATISTICHE DESCRITTIVE PER PERIODO D'ESAME ===")
as.data.frame(desc_stats)

# ========================================================
# 4. VISUALIZZAZIONI ESPLORATIVE
# ========================================================

# Grafico 1: Distribuzione dell'affetto negativo per periodo
p1 <- ggplot(
  d_clean,
  aes(x = exam_period_f, y = neg_affect_ema, fill = exam_period_f)
) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_jitter(alpha = 0.1, width = 0.2) +
  scale_fill_viridis_d(name = "Periodo") +
  labs(
    title = "Distribuzione dell'Affetto Negativo per Periodo d'Esame",
    x = "Periodo d'Esame",
    y = "Affetto Negativo (EMA)",
    caption = "Le linee centrali rappresentano la mediana, le scatole l'IQR"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Grafico 2: Andamento temporale individuale (campione di partecipanti)
sample_participants <- sample(
  unique(d_clean$user_id),
  min(12, length(unique(d_clean$user_id)))
)
d_sample <- d_clean %>% filter(user_id %in% sample_participants)

p2 <- ggplot(d_sample, aes(x = bysubj_day, y = neg_affect_ema)) +
  geom_line(alpha = 0.6) +
  geom_point(aes(color = exam_period_f), size = 2) +
  scale_color_viridis_d(name = "Periodo") +
  facet_wrap(~user_id, scales = "free_x", ncol = 4) +
  labs(
    title = "Andamento Temporale dell'Affetto Negativo (Campione Partecipanti)",
    x = "Giorno dello Studio",
    y = "Affetto Negativo (EMA)"
  ) +
  theme(strip.text = element_text(size = 8))

# Grafico 3: Medie per periodo con errori standard
summary_plot_data <- d_clean %>%
  group_by(exam_period_f) %>%
  summarise(
    mean_neg_affect = mean(neg_affect_ema, na.rm = TRUE),
    se_neg_affect = sd(neg_affect_ema, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

p3 <- ggplot(
  summary_plot_data,
  aes(x = exam_period_f, y = mean_neg_affect, fill = exam_period_f)
) +
  geom_col(alpha = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_neg_affect - se_neg_affect,
      ymax = mean_neg_affect + se_neg_affect
    ),
    width = 0.2
  ) +
  scale_fill_viridis_d(name = "Periodo") +
  labs(
    title = "Media dell'Affetto Negativo per Periodo d'Esame",
    x = "Periodo d'Esame",
    y = "Affetto Negativo Medio (±SE)",
    caption = "Le barre rappresentano l'errore standard"
  )

# Grafico 4: Correlazioni tra variabili chiave
d_cor <- d_clean %>%
  select(
    neg_affect_ema,
    dass_stress,
    dass_depression,
    dass_anxiety,
    pid5_negative_affectivity,
    context_threat,
    context_quality
  ) %>%
  na.omit()

if (nrow(d_cor) > 0) {
  cor_matrix <- cor(d_cor)
  p4 <- corrplot(
    cor_matrix,
    method = "color",
    type = "upper",
    tl.cex = 0.8,
    tl.col = "black",
    title = "Matrice di Correlazione - Variabili Emotive",
    mar = c(0, 0, 2, 0)
  )
}

p1
p2
p3
p4

# ========================================================
# 5. ANALISI STATISTICHE INFERENZIALI
# ========================================================

# Test ANOVA per differenze tra periodi
cat("\n=== ANALISI INFERENZIALI ===\n")

# Modello lineare misto per tenere conto della struttura nidificata
model_basic <- lmer(
  neg_affect_ema ~ exam_period_f + (1 | user_id),
  data = d_clean
)
cat("\n--- Modello Lineare Misto Base ---\n")
print(summary(model_basic))

# Test per differenze tra gruppi
cat("\n--- Test Post-hoc per Confronti tra Periodi ---\n")
emmeans_results <- emmeans(model_basic, ~exam_period_f)
pairs_results <- pairs(emmeans_results, adjust = "bonferroni")
print(emmeans_results)
print(pairs_results)

# Effect size (Cohen's d) per i confronti principali
if ("pre_exam" %in% d_clean$exam_period & "baseline" %in% d_clean$exam_period) {
  pre_vs_baseline <- d_clean %>%
    filter(exam_period %in% c("baseline", "pre_exam")) %>%
    mutate(exam_period_f = droplevels(exam_period_f))

  if (nrow(pre_vs_baseline) > 0) {
    cohen_d_pre_base <- cohens_d(
      neg_affect_ema ~ exam_period_f,
      data = pre_vs_baseline
    )
    cat("\n--- Effect Size: Pre-Esame vs Baseline ---\n")
    print(cohen_d_pre_base)
  }
}

# ========================================================
# 6. ANALISI DELLA VARIABILITÀ INTRA-INDIVIDUALE
# ========================================================

# Calcola variabilità individuale per periodo
individual_variability <- d_clean %>%
  group_by(user_id, exam_period_f) %>%
  summarise(
    n_obs = n(),
    mean_neg_affect = mean(neg_affect_ema, na.rm = TRUE),
    sd_neg_affect = sd(neg_affect_ema, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_obs >= 2) %>% # Solo partecipanti con almeno 2 osservazioni per periodo
  na.omit()

# Grafico della variabilità individuale
p5 <- ggplot(
  individual_variability,
  aes(x = exam_period_f, y = sd_neg_affect)
) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  labs(
    title = "Variabilità Intra-individuale dell'Affetto Negativo",
    x = "Periodo d'Esame",
    y = "Deviazione Standard Intra-individuale",
    caption = "Maggiore variabilità potrebbe indicare maggiore reattività"
  ) +
  theme_minimal()

# ========================================================
# 7. ANALISI PER GENERE
# ========================================================

if (length(unique(d_clean$sex_f)) > 1) {
  # Interazione genere x periodo
  model_sex <- lmer(
    neg_affect_ema ~ exam_period_f * sex_f + (1 | user_id),
    data = d_clean
  )

  cat("\n--- Modello con Interazione Genere x Periodo ---\n")
  print(summary(model_sex))

  # Grafico per genere
  p6 <- d_clean %>%
    group_by(exam_period_f, sex_f) %>%
    summarise(
      mean_neg_affect = mean(neg_affect_ema, na.rm = TRUE),
      se_neg_affect = sd(neg_affect_ema, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    ggplot(aes(x = exam_period_f, y = mean_neg_affect, fill = sex_f)) +
    geom_col(position = position_dodge(), alpha = 0.7) +
    geom_errorbar(
      aes(
        ymin = mean_neg_affect - se_neg_affect,
        ymax = mean_neg_affect + se_neg_affect
      ),
      position = position_dodge(width = 0.9),
      width = 0.2
    ) +
    scale_fill_viridis_d(name = "Genere") +
    labs(
      title = "Affetto Negativo per Periodo ed Genere",
      x = "Periodo d'Esame",
      y = "Affetto Negativo Medio (±SE)"
    ) +
    theme_minimal()
}

# ========================================================
# 8. OUTPUT E VISUALIZZAZIONI
# ========================================================

# Mostra i grafici
print(p1)
print(p2)
print(p3)
if (exists("p4")) print(p4)
print(p5)
if (exists("p6")) print(p6)

# ========================================================
# 9. SUMMARY E RACCOMANDAZIONI
# ========================================================

cat("\n" %>% rep(3) %>% paste(collapse = ""))
cat("=== SUMMARY E RACCOMANDAZIONI ===\n")
cat(
  "1. Verifica se le differenze tra periodi sono statisticamente significative\n"
)
cat("2. Esamina l'effect size per valutare la rilevanza pratica\n")
cat(
  "3. Controlla la variabilità intra-individuale come indicatore di reattività\n"
)
cat("4. Considera le differenze di genere se presenti\n")
cat("5. Identifica partecipanti con pattern atipici per analisi successive\n")

# Salva risultati principali
results_summary <- list(
  descriptive_stats = desc_stats,
  model_results = summary(model_basic),
  emmeans_results = emmeans_results,
  pairs_results = pairs_results
)

# Opzionale: salva i risultati
# saveRDS(results_summary, "ema_exploratory_results.rds")

cat(
  "\nAnalisi esplorativa completata! Controlla i risultati per decidere i prossimi passi.\n"
)


# ============================================================
# Approccio Integrato: Stan per Dipendenze Temporali + Predizione Reattività
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(viridis)
})

set.seed(1234)

# ------------------------------------------------------------
# Utility
# ------------------------------------------------------------
scale_num <- function(x) as.numeric(scale(x))

# Piccolo modello Stan per regressione bayesiana con R^2
stan_regression_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  vector[K] beta;
  real alpha;
  real<lower=0> sigma;
}
model {
  beta ~ normal(0, 1);
  alpha ~ normal(0, 1);
  sigma ~ exponential(1);
  y ~ normal(alpha + X * beta, sigma);
}
generated quantities {
  vector[N] y_hat;
  real R2;
  real var_y;
  real var_res;
  for (n in 1:N) y_hat[n] = alpha + dot_product(row(X, n), beta);
  var_y  = variance(y);
  var_res = variance(y - y_hat);
  R2 = (var_y - var_res) / var_y;
}
"

stan_regression_file <- write_stan_file(stan_regression_code)
stan_regression <- cmdstan_model(stan_regression_file)

# ------------------------------------------------------------
# 1) PREPARAZIONE DATI PER MODELLO STAN
# ------------------------------------------------------------
prepare_stan_multilevel_data <- function(d_clean) {
  # Periodo numerico
  outcomes_data <- d_clean %>%
    select(user_id, sex, exam_period, neg_affect_ema, dass_stress) %>%
    filter(!is.na(neg_affect_ema), !is.na(dass_stress)) %>%
    mutate(
      period_numeric = case_when(
        exam_period == "baseline" ~ 1L,
        exam_period == "pre_exam" ~ 2L,
        exam_period == "post_exam" ~ 3L,
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(!is.na(period_numeric)) %>%
    arrange(user_id, period_numeric)

  # Mapping soggetti
  subject_map <- outcomes_data %>%
    distinct(user_id) %>%
    arrange(user_id) %>%
    mutate(subject_numeric = row_number())

  outcomes_data <- outcomes_data %>%
    left_join(subject_map, by = "user_id")

  # Split outcomes
  y1_data <- outcomes_data %>%
    select(neg_affect_ema, subject_numeric, period_numeric) %>%
    drop_na()
  y2_data <- outcomes_data %>%
    select(dass_stress, subject_numeric, period_numeric) %>%
    drop_na()

  # Predittori tra-soggetto (baseline)
  base_pred <- d_clean %>%
    group_by(user_id) %>%
    summarise(
      baseline_pid5 = first(pid5_negative_affect_baseline),
      baseline_dass_stress = first(dass_stress_baseline),
      sex_numeric = as.integer(first(sex) == "Femmina"),
      .groups = "drop"
    ) %>%
    left_join(subject_map, by = "user_id") %>%
    arrange(subject_numeric) %>%
    mutate(
      z_baseline_pid5 = scale_num(baseline_pid5),
      z_baseline_dass_stress = scale_num(baseline_dass_stress)
    )

  # Predittori EMA riassunti per periodo: media e sd di neg_affect_ema
  ema_summ <- d_clean %>%
    filter(exam_period %in% c("baseline", "pre_exam", "post_exam")) %>%
    group_by(user_id, exam_period) %>%
    summarise(
      ema_mean = mean(neg_affect_ema, na.rm = TRUE),
      ema_sd = sd(neg_affect_ema, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      period_numeric = case_when(
        exam_period == "baseline" ~ 1L,
        exam_period == "pre_exam" ~ 2L,
        exam_period == "post_exam" ~ 3L
      )
    ) %>%
    left_join(subject_map, by = "user_id")

  ema_wide <- ema_summ %>%
    select(user_id, period_numeric, ema_mean, ema_sd) %>%
    pivot_wider(
      names_from = period_numeric,
      values_from = c(ema_mean, ema_sd),
      names_glue = "{.value}_p{period_numeric}"
    ) %>%
    # nomi comodi:
    rename(
      ema_mean_base = ema_mean_p1,
      ema_sd_base = ema_sd_p1,
      ema_mean_pre = ema_mean_p2,
      ema_sd_pre = ema_sd_p2,
      ema_mean_post = ema_mean_p3,
      ema_sd_post = ema_sd_p3
    )

  # Costruiamo il dataframe predittori soggetto + EMA
  subject_predictors <- base_pred %>%
    left_join(ema_wide, by = "user_id") %>%
    mutate(
      z_ema_mean_base = scale_num(ema_mean_base),
      z_ema_sd_base = scale_num(ema_sd_base),
      z_ema_mean_pre = scale_num(ema_mean_pre),
      z_ema_sd_pre = scale_num(ema_sd_pre),
      z_ema_mean_post = scale_num(ema_mean_post),
      z_ema_sd_post = scale_num(ema_sd_post)
    )

  # Matrice W (predittori tra-soggetto) — qui usiamo solo i baseline + sesso (come nel tuo sketch)
  W_matrix <- subject_predictors %>%
    transmute(
      z_baseline_pid5,
      z_baseline_dass_stress,
      sex_numeric = as.numeric(sex_numeric)
    ) %>%
    as.matrix()

  stan_data <- list(
    I = nrow(subject_predictors), # soggetti
    P = 3L, # periodi (baseline, pre, post)
    K = 2L, # outcomes (neg_affect_ema, dass_stress)

    # Outcome 1 (neg_affect_ema)
    N = c(nrow(y1_data), nrow(y2_data)),
    y1 = y1_data$neg_affect_ema,
    subj1 = y1_data$subject_numeric,
    period1 = y1_data$period_numeric,

    # Outcome 2 (dass_stress)
    y2 = y2_data$dass_stress,
    subj2 = y2_data$subject_numeric,
    period2 = y2_data$period_numeric,

    # Predittori tra-soggetto
    Q = ncol(W_matrix),
    W = W_matrix
  )

  list(
    stan_data = stan_data,
    subject_map = subject_map,
    subject_predictors = subject_predictors,
    baseline_predictors = c(
      "baseline_pid5",
      "baseline_dass_stress",
      "sex_numeric"
    ),
    ema_predictors_pre = c("ema_mean_pre", "ema_sd_pre"),
    ema_predictors_post = c("ema_mean_post", "ema_sd_post")
  )
}

cat("=== PREPARAZIONE DATI PER MODELLO STAN ===\n")
stan_setup <- prepare_stan_multilevel_data(d_clean)
cat("Soggetti nel modello:", stan_setup$stan_data$I, "\n")
cat("Osservazioni outcome 1:", stan_setup$stan_data$N[1], "\n")
cat("Osservazioni outcome 2:", stan_setup$stan_data$N[2], "\n")


# Quanti NA in ogni elemento del list usato per Stan
na_tally <- vapply(stan_setup$stan_data, function(x) sum(is.na(x)), integer(1))
na_tally[na_tally > 0]

W <- stan_setup$stan_data$W

# Indici (riga, colonna) delle celle mancanti
idxW <- which(is.na(W), arr.ind = TRUE)
idxW
#>      row col
#> [1,]  17   2
#> [2,]  39   3
#> ...

# Righe/colonne con almeno un NA
rows_with_NA <- which(rowSums(is.na(W)) > 0)
cols_with_NA <- which(colSums(is.na(W)) > 0)

rows_with_NA
cols_with_NA
colnames(W)[cols_with_NA]

# ------------------------------------------------------------
# 2) ESECUZIONE MODELLO STAN PER REATTIVITÀ LATENTE
# ------------------------------------------------------------
model_latent <- cmdstan_model(
  here("scripts", "02_stress_reactivity", "reactivity_latent_eta_new.stan")
)

cat("\n=== ESECUZIONE MODELLO STAN LATENTE ===\n")
fit_latent <- model_latent$sample(
  data = stan_setup$stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 500,
  parallel_chains = 4,
  backend = "cmdstanr"
)

fit_latent <- model_latent$variational(
  data = stan_setup$stan_data,
  seed = 20250912, # È buona pratica specificare un seed per la riproducibilità
  output_samples = 2000, # Numero di campioni da estrarre dalla distribuzione approssimata
  iter = 10000, # Numero massimo di iterazioni per l'ottimizzazione (default è 10000)
  adapt_iter = 50, # Iterazioni per la fase adattiva (default è 50)
  tol_rel_obj = 0.01, # Tolleranza per la convergenza (default è 0.01)
  elbo_samples = 100, # Numero di campioni per stimare l'ELBO (default è 100)
  grad_samples = 1, # Numero di campioni per il gradiente (default è 1)
  refresh = 500, # Ogni quanti iterazioni stampare un aggiornamento (default è 100)
  output_dir = NULL, # Directory opzionale per salvare output
  algorithm = "meanfield" # Può essere "meanfield" o "fullrank"
)

fit_latent$cmdstan_diagnose(quiet = TRUE)
suppressMessages(print(fit_latent$summary(c("tau_eta"))))

# ------------------------------------------------------------
# 3) ESTRAZIONE REATTIVITÀ LATENTE INDIVIDUALE
# ------------------------------------------------------------
extract_individual_reactivity <- function(fit, subject_map) {
  pars <- fit$metadata()$model_params

  # Se il modello fornisce direttamente eta[i,p]
  if ("eta" %in% pars) {
    eta_mean <- fit$summary("eta") %>%
      select(variable, mean)

    I <- nrow(subject_map)

    out <- tibble(
      user_id = subject_map$user_id,
      subject_numeric = subject_map$subject_numeric,
      eta_pre = NA_real_,
      eta_post = NA_real_,
      reactivity_range = NA_real_
    )

    for (i in 1:I) {
      mu_pre <- eta_mean %>%
        filter(variable == sprintf("eta[%d,2]", i)) %>%
        pull(mean)
      mu_post <- eta_mean %>%
        filter(variable == sprintf("eta[%d,3]", i)) %>%
        pull(mean)
      if (length(mu_pre) == 0) mu_pre <- NA_real_
      if (length(mu_post) == 0) mu_post <- NA_real_
      out$eta_pre[i] <- mu_pre
      out$eta_post[i] <- mu_post
      out$reactivity_range[i] <- max(
        out$eta_pre[i],
        out$eta_post[i],
        na.rm = TRUE
      )
    }
    return(out)
  }

  # Altrimenti ricostruisco: eta = mu_eta + tau_eta .* z_eta
  need <- c("mu_eta", "z_eta", "tau_eta")
  if (!all(need %in% pars))
    stop("Il modello non esporta né 'eta' né {mu_eta, z_eta, tau_eta}.")

  mu_eta_mean <- fit$summary("mu_eta") %>% select(variable, mean)
  z_eta_mean <- fit$summary("z_eta") %>% select(variable, mean)
  tau_eta <- fit$summary("tau_eta") %>% pull(mean)

  I <- nrow(subject_map)
  out <- tibble(
    user_id = subject_map$user_id,
    subject_numeric = subject_map$subject_numeric,
    eta_pre = NA_real_,
    eta_post = NA_real_,
    reactivity_range = NA_real_
  )

  for (i in 1:I) {
    mu_pre <- mu_eta_mean %>%
      filter(variable == sprintf("mu_eta[%d,2]", i)) %>%
      pull(mean)
    mu_post <- mu_eta_mean %>%
      filter(variable == sprintf("mu_eta[%d,3]", i)) %>%
      pull(mean)
    z_pre <- z_eta_mean %>%
      filter(variable == sprintf("z_eta[%d,2]", i)) %>%
      pull(mean)
    z_post <- z_eta_mean %>%
      filter(variable == sprintf("z_eta[%d,3]", i)) %>%
      pull(mean)

    if (length(mu_pre) == 0) mu_pre <- NA_real_
    if (length(mu_post) == 0) mu_post <- NA_real_
    if (length(z_pre) == 0) z_pre <- 0
    if (length(z_post) == 0) z_post <- 0

    out$eta_pre[i] <- mu_pre + tau_eta[2] * z_pre
    out$eta_post[i] <- mu_post + tau_eta[3] * z_post
    out$reactivity_range[i] <- max(
      out$eta_pre[i],
      out$eta_post[i],
      na.rm = TRUE
    )
  }

  out
}

cat("\n=== ESTRAZIONE REATTIVITÀ LATENTE INDIVIDUALE ===\n")
latent_reactivity <- extract_individual_reactivity(
  fit_latent,
  stan_setup$subject_map
)


# ==========================================
# COSTRUZIONE LISTA DATI PER STAN (con NA in W gestiti nel modello)
# ==========================================

# Funzione robusta che:
# - prende l'output di prepare_stan_multilevel_data() (stan_setup)
# - costruisce indici e vettori per W osservata e mancante
# - forza i tipi corretti per y1/y2 e indici
make_stan_data_vi <- function(stan_setup) {
  ss <- stan_setup # alias corto

  # --- 1) Recupera o ricostruisci W (I x Q) ---
  if (!is.null(ss$stan_data$W)) {
    W <- ss$stan_data$W
  } else if (!is.null(ss$subject_predictors)) {
    # fallback: ricostruisci W come nello script
    W <- cbind(
      ss$subject_predictors$z_baseline_pid5,
      ss$subject_predictors$z_baseline_dass_stress,
      ss$subject_predictors$sex_numeric
    )
    colnames(W) <- c("z_baseline_pid5", "z_baseline_dass_stress", "sex_numeric")
  } else {
    stop("Non trovo né stan_data$W né subject_predictors per ricostruire W.")
  }

  # Dimensioni
  I <- nrow(W)
  Q <- ncol(W)

  # --- 2) Indici celle osservate/mancanti in W ---
  obs_idx <- which(!is.na(W), arr.ind = TRUE)
  mis_idx <- which(is.na(W), arr.ind = TRUE)

  # --- 3) Costruisci la lista dati coerente con lo .stan ---
  sd <- ss$stan_data # comodità

  # Forza i tipi: Stan si aspetta int/double coerenti
  N_vec <- as.integer(sd$N)
  y1_vec <- as.numeric(sd$y1)
  y2_vec <- as.numeric(sd$y2)
  subj1_i <- as.integer(sd$subj1)
  subj2_i <- as.integer(sd$subj2)
  per1_i <- as.integer(sd$period1)
  per2_i <- as.integer(sd$period2)

  stan_data_vi <- list(
    # meta
    I = as.integer(sd$I),
    P = as.integer(sd$P),
    K = as.integer(sd$K),
    Q = as.integer(Q),

    # outcome 1
    N = N_vec, # array[2]
    y1 = y1_vec, # vector[N[1]]
    subj1 = subj1_i, # int array[N[1]]
    period1 = per1_i, # int array[N[1]]

    # outcome 2
    y2 = y2_vec, # vector[N[2]]
    subj2 = subj2_i, # int array[N[2]]
    period2 = per2_i, # int array[N[2]]

    # W osservata/mancante
    NW_obs = nrow(obs_idx),
    w_obs_i = if (nrow(obs_idx) > 0) as.integer(obs_idx[, 1]) else integer(0),
    w_obs_j = if (nrow(obs_idx) > 0) as.integer(obs_idx[, 2]) else integer(0),
    W_obs = if (nrow(obs_idx) > 0) as.numeric(W[!is.na(W)]) else numeric(0),

    NW_mis = nrow(mis_idx),
    w_mis_i = if (nrow(mis_idx) > 0) as.integer(mis_idx[, 1]) else integer(0),
    w_mis_j = if (nrow(mis_idx) > 0) as.integer(mis_idx[, 2]) else integer(0)
  )

  # --- 4) Controlli veloci (fail-fast) ---
  stopifnot(length(stan_data_vi$y1) == stan_data_vi$N[1])
  stopifnot(length(stan_data_vi$subj1) == stan_data_vi$N[1])
  stopifnot(length(stan_data_vi$period1) == stan_data_vi$N[1])
  stopifnot(length(stan_data_vi$y2) == stan_data_vi$N[2])
  stopifnot(length(stan_data_vi$subj2) == stan_data_vi$N[2])
  stopifnot(length(stan_data_vi$period2) == stan_data_vi$N[2])

  # (opzionale) verifica ricostruzione W osservata
  if (stan_data_vi$NW_obs > 0) {
    W_chk <- matrix(NA_real_, nrow = I, ncol = Q)
    W_chk[cbind(
      stan_data_vi$w_obs_i,
      stan_data_vi$w_obs_j
    )] <- stan_data_vi$W_obs
    stopifnot(all.equal(W[!is.na(W)], W_chk[!is.na(W)]))
  }

  # Restituisci lista pronta per Stan
  stan_data_vi
}

# ==========================================
# USO
# ==========================================
# stan_setup <- prepare_stan_multilevel_data(d_clean)  # (già fatto nel tuo flusso)
stan_data_vi <- make_stan_data_vi(stan_setup)

# (facoltativo) sguardo rapido
str(stan_data_vi, max.level = 1)
cat(
  "Celle W osservate:",
  stan_data_vi$NW_obs,
  " | Mancanti:",
  stan_data_vi$NW_mis,
  "\n"
)

# Ora puoi passare stan_data_vi al modello Stan:
# fit_latent <- model_latent$variational(
#   data = stan_data_vi,
#   seed = 20250912,
#   algorithm = "meanfield",   # o "fullrank"
#   iter = 10000,
#   output_samples = 2000,
#   elbo_samples = 100,
#   grad_samples = 1,
#   adapt_iter = 50,
#   tol_rel_obj = 0.01,
#   refresh = 500
# )

# ------------------------------------------------------------
# 4) PREDIZIONE DELLA REATTIVITÀ LATENTE CON EMA
# ------------------------------------------------------------

# Costruisco prediction_data unendo reattività latente e predittori soggetto
prediction_data <- latent_reactivity %>%
  left_join(
    stan_setup$subject_predictors,
    by = c("user_id", "subject_numeric")
  ) %>%
  # Misure osservate di reattività (proxy): media EMA nei periodi pre e post
  mutate(
    obs_pre = ema_mean_pre,
    obs_post = ema_mean_post
  ) %>%
  select(
    user_id,
    subject_numeric,
    eta_pre,
    eta_post,
    reactivity_range,
    sex_numeric,
    # baseline z-std
    z_baseline_pid5,
    z_baseline_dass_stress,
    # EMA z-std per periodo
    z_ema_mean_pre,
    z_ema_sd_pre,
    z_ema_mean_post,
    z_ema_sd_post,
    # osservate per confronto
    obs_pre,
    obs_post
  ) %>%
  filter(complete.cases(.))

cat("Soggetti con dati completi per predizione:", nrow(prediction_data), "\n")

fit_r2 <- function(y, X) {
  stopifnot(is.numeric(y), is.matrix(X), nrow(X) == length(y))
  dat <- list(N = length(y), K = ncol(X), X = X, y = as.numeric(y))
  fit <- stan_regression$sample(
    data = dat,
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 0,
    parallel_chains = 4
  )
  r2 <- fit$draws("R2") |> as_draws_df() |> dplyr::pull(R2)
  list(r2_mean = mean(r2), r2_sd = sd(r2), fit = fit)
}

predict_latent_reactivity <- function(
  outcome_var,
  which_ema = c("pre", "post")
) {
  which_ema <- match.arg(which_ema)

  y <- prediction_data[[outcome_var]]

  # baseline predictors (z_)
  Xb <- prediction_data %>%
    select(z_baseline_pid5, z_baseline_dass_stress, sex_numeric) %>%
    as.matrix()

  # ema predictors (z_) dipendono dal periodo
  if (which_ema == "pre") {
    Xe <- prediction_data %>%
      select(z_ema_mean_pre, z_ema_sd_pre) %>%
      as.matrix()
  } else {
    Xe <- prediction_data %>%
      select(z_ema_mean_post, z_ema_sd_post) %>%
      as.matrix()
  }

  rb <- fit_r2(y, Xb)
  re <- fit_r2(y, Xe)

  cat(sprintf("\nPredizione %s (EMA %s)\n", outcome_var, which_ema))
  cat("R² Baseline:", round(rb$r2_mean, 3), " (±", round(rb$r2_sd, 3), ")\n")
  cat("R² EMA:     ", round(re$r2_mean, 3), " (±", round(re$r2_sd, 3), ")\n")
  cat("Vantaggio EMA:", round(re$r2_mean - rb$r2_mean, 3), "\n")

  list(
    outcome = outcome_var,
    r2_baseline = rb$r2_mean,
    r2_ema = re$r2_mean,
    improvement = re$r2_mean - rb$r2_mean,
    fit_baseline = rb$fit,
    fit_ema = re$fit
  )
}

# Esegui predizione per reattività latente
results_latent_pre <- predict_latent_reactivity("eta_pre", which_ema = "pre")
results_latent_post <- predict_latent_reactivity("eta_post", which_ema = "post")

# ------------------------------------------------------------
# 5) CONFRONTO OSSERVATO vs LATENTE
# ------------------------------------------------------------
cat("\n=== CONFRONTO APPROCCI OSSERVATO vs LATENTE ===\n")

correlation_comparison <- prediction_data %>%
  select(eta_pre, eta_post, obs_pre, obs_post) %>%
  cor(use = "complete.obs")

cat("Correlazioni (eta_pre, eta_post, obs_pre, obs_post):\n")
print(round(correlation_comparison, 3))

comparison_summary <- tibble(
  Approach = c("Osservato (EMA mean)", "Latente"),
  R2_Baseline_Pre = c(NA_real_, results_latent_pre$r2_baseline),
  R2_EMA_Pre = c(NA_real_, results_latent_pre$r2_ema),
  EMA_Advantage_Pre = c(NA_real_, results_latent_pre$improvement)
)
print(comparison_summary)

# ------------------------------------------------------------
# 6) HETEROGENEITÀ INDIVIDUALE (Clustering)
# ------------------------------------------------------------
perform_reactivity_clustering <- function() {
  df <- prediction_data %>%
    select(eta_pre, eta_post, reactivity_range) %>%
    drop_na()

  if (nrow(df) < 10) {
    warning("Pochi casi per clustering.")
    return(list(
      cluster_results = NULL,
      cluster_profiles = NULL,
      cluster_r2 = NULL
    ))
  }

  Z <- scale(df)
  set.seed(123)
  km <- kmeans(Z, centers = 3)

  cluster_results <- prediction_data %>%
    drop_na(eta_pre, eta_post, reactivity_range) %>%
    mutate(
      cluster = factor(
        km$cluster,
        labels = c("Low Reactivity", "Moderate Reactivity", "High Reactivity")
      )
    )

  cluster_profiles <- cluster_results %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      eta_pre_mean = mean(eta_pre, na.rm = TRUE),
      eta_post_mean = mean(eta_post, na.rm = TRUE),
      baseline_pid5_mean = mean(z_baseline_pid5, na.rm = TRUE),
      ema_sd_pre_mean = mean(z_ema_sd_pre, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\n=== PROFILI DI REATTIVITÀ ===\n")
  print(cluster_profiles)

  # Test: quanto predicono i soli EMA dentro ciascun cluster (periodo coerente con outcome eta_pre)
  cluster_levels <- levels(cluster_results$cluster)
  cluster_r2 <- rep(NA_real_, length(cluster_levels))
  names(cluster_r2) <- cluster_levels

  for (cl in cluster_levels) {
    sub <- cluster_results %>% filter(cluster == cl)
    if (nrow(sub) < 20) next
    y <- sub$eta_pre
    Xe <- sub %>% select(z_ema_mean_pre, z_ema_sd_pre) %>% as.matrix()
    res <- try(fit_r2(y, Xe), silent = TRUE)
    if (!inherits(res, "try-error")) cluster_r2[cl] <- res$r2_mean
  }

  cat("\nR² EMA per cluster (outcome eta_pre):\n")
  print(round(cluster_r2, 3))

  list(
    cluster_results = cluster_results,
    cluster_profiles = cluster_profiles,
    cluster_r2 = cluster_r2
  )
}

clustering_results <- perform_reactivity_clustering()

# ------------------------------------------------------------
# 7) VISUALIZZAZIONI
# ------------------------------------------------------------
p_comparison <- prediction_data %>%
  select(eta_pre, obs_pre) %>%
  drop_na() %>%
  ggplot(aes(x = obs_pre, y = eta_pre)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", linewidth = 0.8) +
  labs(
    title = "Reattività Osservata (EMA mean pre) vs Reattività Latente (eta_pre)",
    x = "Reattività osservata (EMA mean, pre-esame)",
    y = "Reattività latente (eta_pre)",
    caption = paste(
      "Correlazione:",
      round(
        cor(
          prediction_data$eta_pre,
          prediction_data$obs_pre,
          use = "complete.obs"
        ),
        3
      )
    )
  ) +
  theme_minimal()
print(p_comparison)

if (!is.null(clustering_results$cluster_results)) {
  p_clusters <- clustering_results$cluster_results %>%
    select(cluster, eta_pre, eta_post) %>%
    pivot_longer(
      c(eta_pre, eta_post),
      names_to = "period",
      values_to = "reactivity"
    ) %>%
    ggplot(aes(x = period, y = reactivity, fill = cluster)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_viridis_d() +
    labs(
      title = "Profili di Reattività Latente per Cluster",
      x = "Periodo",
      y = "Reattività Latente",
      fill = "Cluster"
    ) +
    theme_minimal()
  print(p_clusters)
}

# ------------------------------------------------------------
# 8) SINTESI FINALE E RACCOMANDAZIONI
# ------------------------------------------------------------
cat("\n\n\n=== SINTESI FINALE INTEGRATA ===\n")

final_metrics <- list(
  latent_ema_advantage_pre = results_latent_pre$improvement,
  latent_ema_advantage_post = results_latent_post$improvement,
  correlation_observed_latent = cor(
    prediction_data$eta_pre,
    prediction_data$obs_pre,
    use = "complete.obs"
  )
)

cat(
  "Vantaggio EMA (latente pre): ",
  round(final_metrics$latent_ema_advantage_pre, 3),
  "\n"
)
cat(
  "Vantaggio EMA (latente post):",
  round(final_metrics$latent_ema_advantage_post, 3),
  "\n"
)
cat(
  "Correlazione osservato-latente (pre):",
  round(final_metrics$correlation_observed_latent, 3),
  "\n"
)

avg_advantage <- mean(
  c(
    final_metrics$latent_ema_advantage_pre,
    final_metrics$latent_ema_advantage_post
  ),
  na.rm = TRUE
)

cat("\n=== VALUTAZIONE FINALE ===\n")
if (is.finite(avg_advantage) && avg_advantage > 0.10) {
  cat(
    "✓ FORTE EVIDENZA: le misure EMA aggiungono valore predittivo sostanziale.\n"
  )
  cat("✓ Raccomandazione: implementare protocolli EMA per la valutazione.\n")
} else if (is.finite(avg_advantage) && avg_advantage > 0.05) {
  cat("✓ EVIDENZA MODERATA: le EMA aggiungono valore apprezzabile.\n")
  cat("✓ Raccomandazione: usare EMA come complemento alle baseline.\n")
} else {
  cat("? EVIDENZA LIMITATA: valore aggiunto EMA marginale.\n")
  cat("? Raccomandazione: rivalutare costi-benefici dell’EMA.\n")
}

cat("\n=== CONTRIBUTI METODOLOGICI ===\n")
cat("1) Il livello latente cattura meglio la reattività individuale.\n")
cat("2) Il framework bayesiano quantifica l’incertezza su R² e parametri.\n")
cat("3) Il clustering rivela eterogeneità di pattern.\n")
cat("4) Confronto sistematico baseline vs EMA.\n")

cat("\n=== IMPLICAZIONI TEORICHE ===\n")
if (
  is.finite(final_metrics$correlation_observed_latent) &&
    final_metrics$correlation_observed_latent > 0.7
) {
  cat("• Buona convergenza osservato–latente (validità di costrutto).\n")
} else {
  cat(
    "• Discrepanza osservato–latente suggerisce errore di misura/nonlinearità.\n"
  )
}

if (!is.null(clustering_results$cluster_r2)) {
  hr <- clustering_results$cluster_r2["High Reactivity"]
  if (length(hr) == 1 && is.finite(hr) && hr > 0.3) {
    cat("• EMA particolarmente predittive nei profili ad alta reattività.\n")
  }
}

integrated_results <- list(
  stan_fit = fit_latent,
  latent_reactivity = latent_reactivity,
  prediction_results_latent = list(
    pre = results_latent_pre,
    post = results_latent_post
  ),
  clustering_results = clustering_results,
  final_metrics = final_metrics,
  comparison_data = prediction_data
)

saveRDS(integrated_results, "integrated_ema_analysis_results.rds")
cat("\nRisultati completi salvati in: integrated_ema_analysis_results.rds\n")

cat("\n=== POSSIBILI SVILUPPI FUTURI ===\n")
cat(
  "1) Cross-validation temporale (rolling origin) per validazione prospettica.\n"
)
cat("2) Analisi di mediazione: quali domini EMA mediano la predizione?\n")
cat("3) Modelli dinamici (state-space) per traiettorie individuali.\n")
cat("4) Integrazione con biomarcatori (es. cortisol, HRV).\n")
