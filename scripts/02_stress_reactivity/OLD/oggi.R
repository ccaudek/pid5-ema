# ============================================================
# Reattività latente con VI in Stan (imputazione missing in W)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(viridis)
  library(conflicted)
})

# Risoluzione conflitti (utile se caricherai altri pacchetti dopo)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

set.seed(20250912)


# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
scale_num <- function(x) as.numeric(scale(x))

# ------------------------------------------------------------
# 1) PREPARAZIONE DATI PER IL MODELLO STAN
#    (costruisce outcomes, mapping soggetti, W con eventuali NA)
# ------------------------------------------------------------
prepare_stan_multilevel_data <- function(d_clean) {
  # 1a) Periodo numerico e selezione outcomes
  outcomes_data <- d_clean %>%
    select(user_id, sex, exam_period, neg_affect_ema, dass_stress) %>%
    filter(!is.na(exam_period)) %>%
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

  # 1b) Mappa soggetti
  subject_map <- outcomes_data %>%
    distinct(user_id) %>%
    arrange(user_id) %>%
    mutate(subject_numeric = row_number())

  outcomes_data <- outcomes_data %>% left_join(subject_map, by = "user_id")

  # 1c) Split per outcome (tieni solo righe con outcome presente)
  y1_data <- outcomes_data %>%
    select(neg_affect_ema, subject_numeric, period_numeric) %>%
    drop_na()

  y2_data <- outcomes_data %>%
    select(dass_stress, subject_numeric, period_numeric) %>%
    drop_na()

  # 1d) Predittori tra-soggetto (baseline)
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

  # (opzionale) riassunti EMA per uso analitico successivo (non entrano in Stan)
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
    rename(
      ema_mean_base = ema_mean_p1,
      ema_sd_base = ema_sd_p1,
      ema_mean_pre = ema_mean_p2,
      ema_sd_pre = ema_sd_p2,
      ema_mean_post = ema_mean_p3,
      ema_sd_post = ema_sd_p3
    )

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

  # 1e) Matrice W (può contenere NA: verranno imputati in Stan)
  W_matrix <- subject_predictors %>%
    transmute(
      z_baseline_pid5,
      z_baseline_dass_stress,
      sex_numeric = as.numeric(sex_numeric) # evita NA qui se possibile
    ) %>%
    as.matrix()
  colnames(W_matrix) <- c(
    "z_baseline_pid5",
    "z_baseline_dass_stress",
    "sex_numeric"
  )

  # 1f) Lista base per Stan (senza ancora gli indici di W: li costruiremo a parte)
  stan_data <- list(
    I = nrow(subject_predictors),
    P = 3L,
    K = 2L,
    N = c(nrow(y1_data), nrow(y2_data)),
    y1 = as.numeric(y1_data$neg_affect_ema),
    subj1 = as.integer(y1_data$subject_numeric),
    period1 = as.integer(y1_data$period_numeric),
    y2 = as.numeric(y2_data$dass_stress),
    subj2 = as.integer(y2_data$subject_numeric),
    period2 = as.integer(y2_data$period_numeric),
    Q = ncol(W_matrix),
    W = W_matrix # utile per costruire gli indici; NON verrà passato allo Stan finale
  )

  list(
    stan_data = stan_data,
    subject_map = subject_map,
    subject_predictors = subject_predictors
  )
}

cat("=== PREPARAZIONE DATI ===\n")
stan_setup <- prepare_stan_multilevel_data(d_clean)
cat(
  "Soggetti:",
  stan_setup$stan_data$I,
  " | N y1:",
  stan_setup$stan_data$N[1],
  " | N y2:",
  stan_setup$stan_data$N[2],
  "\n"
)

# ------------------------------------------------------------
# 2) COSTRUZIONE LISTA DATI PER IL MODELLO STAN (imputazione W interna)
# ------------------------------------------------------------
make_stan_data_vi <- function(stan_setup) {
  sd <- stan_setup$stan_data
  W <- sd$W
  I <- nrow(W)
  Q <- ncol(W)

  obs_idx <- which(!is.na(W), arr.ind = TRUE)
  mis_idx <- which(is.na(W), arr.ind = TRUE)

  list(
    I = as.integer(sd$I),
    P = as.integer(sd$P),
    K = as.integer(sd$K),
    Q = as.integer(Q),

    N = as.integer(sd$N),
    y1 = as.numeric(sd$y1),
    subj1 = as.integer(sd$subj1),
    period1 = as.integer(sd$period1),
    y2 = as.numeric(sd$y2),
    subj2 = as.integer(sd$subj2),
    period2 = as.integer(sd$period2),

    NW_obs = nrow(obs_idx),
    w_obs_i = if (nrow(obs_idx) > 0) as.integer(obs_idx[, 1]) else integer(0),
    w_obs_j = if (nrow(obs_idx) > 0) as.integer(obs_idx[, 2]) else integer(0),
    W_obs = if (nrow(obs_idx) > 0) as.numeric(W[!is.na(W)]) else numeric(0),

    NW_mis = nrow(mis_idx),
    w_mis_i = if (nrow(mis_idx) > 0) as.integer(mis_idx[, 1]) else integer(0),
    w_mis_j = if (nrow(mis_idx) > 0) as.integer(mis_idx[, 2]) else integer(0)
  )
}

stan_data_vi <- make_stan_data_vi(stan_setup)
cat(
  "Celle W: osservate =",
  stan_data_vi$NW_obs,
  " | mancanti =",
  stan_data_vi$NW_mis,
  "\n"
)

# ------------------------------------------------------------
# 3) COMPILAZIONE E VARIATIONAL INFERENCE
# ------------------------------------------------------------
# Imposta il path del modello Stan (aggiorna se necessario)
stan_file <- here(
  "scripts",
  "02_stress_reactivity",
  "reactivity_latent_eta_new.stan"
)
stopifnot(file.exists(stan_file))
model_latent <- cmdstan_model(stan_file)

cat("\n=== VARIATIONAL INFERENCE (meanfield) ===\n")
fit_latent <- model_latent$variational(
  data = stan_data_vi,
  seed = 20250912,
  algorithm = "meanfield", # "fullrank" per maggiore fedeltà (più lento)
  iter = 10000,
  output_samples = 2000,
  elbo_samples = 100,
  grad_samples = 1,
  adapt_iter = 50,
  tol_rel_obj = 0.01,
  refresh = 500
)

# Riepilogo rapido di parametri chiave
print(fit_latent$summary(c("mu0", "sigma", "tau_eta", "lambda")))

# ------------------------------------------------------------
# 4) ESTRAZIONE REATTIVITÀ LATENTE (eta) E MERGE CON PREDITTORI
# ------------------------------------------------------------
extract_individual_reactivity <- function(fit, subject_map, P_expect = 3L) {
  # Stan esporta eta[i,p] in transformed parameters
  smry <- fit$summary("eta") %>% select(variable, mean)
  I <- nrow(subject_map)

  # helper per estrarre mean(eta[i,p])
  get_eta <- function(i, p) {
    v <- smry %>% filter(variable == sprintf("eta[%d,%d]", i, p)) %>% pull(mean)
    if (length(v) == 0) NA_real_ else v
  }

  tibble(
    user_id = subject_map$user_id,
    subject_numeric = subject_map$subject_numeric,
    eta_base = vapply(seq_len(I), \(i) get_eta(i, 1), numeric(1)),
    eta_pre = if (P_expect >= 2)
      vapply(seq_len(I), \(i) get_eta(i, 2), numeric(1)) else NA_real_,
    eta_post = if (P_expect >= 3)
      vapply(seq_len(I), \(i) get_eta(i, 3), numeric(1)) else NA_real_,
    reactivity_range = pmax(eta_pre, eta_post, na.rm = TRUE)
  )
}

latent_reactivity <- extract_individual_reactivity(
  fit_latent,
  stan_setup$subject_map,
  P_expect = stan_data_vi$P
)

# ------------------------------------------------------------
# 5) ANALISI ESSENZIALI
#    - Confronto osservato (EMA mean) vs latente
#    - Correlazioni, R^2 PPC dal modello
#    - Visualizzazioni
# ------------------------------------------------------------
prediction_data <- latent_reactivity %>%
  left_join(
    stan_setup$subject_predictors,
    by = c("user_id", "subject_numeric")
  ) %>%
  mutate(
    obs_pre = ema_mean_pre,
    obs_post = ema_mean_post
  )

cat("\n=== Correlazioni osservato-latente (complete.obs) ===\n")
corr_tbl <- prediction_data %>%
  select(eta_pre, eta_post, obs_pre, obs_post) %>%
  cor(use = "complete.obs")
print(round(corr_tbl, 3))

# R^2 dal generated quantities del modello (posteriori della VI)
r2_df <- fit_latent$draws(c("R2_y1", "R2_y2"), format = "df") %>%
  as_tibble()
cat("\nR^2 posteriori (medie):\n")
r2_means <- colMeans(r2_df, na.rm = TRUE)
print(round(r2_means, 3))

# ------------------------------------------------------------
# 6) PLOT
# ------------------------------------------------------------
p_scatter <- prediction_data %>%
  drop_na(eta_pre, obs_pre) %>%
  ggplot(aes(x = obs_pre, y = eta_pre)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", linewidth = 0.8) +
  labs(
    title = "Reattività osservata (EMA mean pre) vs reattività latente (eta_pre)",
    x = "EMA (media periodo pre)",
    y = "Latente (eta_pre)",
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
print(p_scatter)

# (opzionale) boxplot per periodo
if (all(c("eta_pre", "eta_post") %in% names(prediction_data))) {
  p_box <- prediction_data %>%
    select(eta_pre, eta_post) %>%
    pivot_longer(everything(), names_to = "period", values_to = "eta") %>%
    ggplot(aes(x = period, y = eta, fill = period)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_viridis_d(guide = "none") +
    labs(
      title = "Distribuzione della reattività latente per periodo",
      x = NULL,
      y = "eta"
    ) +
    theme_minimal()
  print(p_box)
}

# ------------------------------------------------------------
# 7) SALVATAGGIO RISULTATI
# ------------------------------------------------------------
integrated_results <- list(
  stan_fit = fit_latent,
  latent_reactivity = latent_reactivity,
  r2 = r2_means,
  correlation_table = corr_tbl,
  subject_predictors = stan_setup$subject_predictors
)
saveRDS(integrated_results, file = "integrated_ema_analysis_results.rds")
cat("\nSalvato: integrated_ema_analysis_results.rds\n")


# ============================================================
# Valore predittivo incrementale EMA (k-fold CV, VI in Stan)
# ============================================================

library(cmdstanr)
library(posterior)
set.seed(20250912)

# -- piccolo modello Stan di regressione (gaussiana) per VI
stan_reg_vi_code <- "
data {
  int<lower=1> N;
  int<lower=0> K;
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
"
stan_reg_vi_file <- write_stan_file(stan_reg_vi_code)
stan_reg_vi <- cmdstan_model(stan_reg_vi_file)

# -- helper: fit VI e predici su test (usa media posteriore)
vi_fit_predict <- function(
  X_train,
  y_train,
  X_test,
  iter = 4000,
  elbo_samples = 100,
  output_samples = 1000
) {
  dat <- list(N = nrow(X_train), K = ncol(X_train), X = X_train, y = y_train)
  fit <- stan_reg_vi$variational(
    data = dat,
    algorithm = "meanfield",
    iter = iter,
    elbo_samples = elbo_samples,
    output_samples = output_samples,
    grad_samples = 1,
    adapt_iter = 50,
    tol_rel_obj = 0.01,
    refresh = 0
  )
  dr <- as_draws_df(fit$draws(c("alpha", "beta"), format = "df"))
  # medie posteriori
  a_hat <- mean(dr$alpha)
  beta_cols <- grep("^beta\\[", names(dr), value = TRUE)
  b_hat <- colMeans(dr[, beta_cols, drop = FALSE])
  y_hat <- as.numeric(a_hat + X_test %*% b_hat)
  list(y_hat = y_hat)
}

# -- R^2 out-of-sample su un fold
r2_oos <- function(y_test, y_hat) {
  sse <- sum((y_test - y_hat)^2)
  sst <- sum((y_test - mean(y_test))^2) # baseline naive = media del test
  1 - sse / sst
}

# -- k-fold CV a livello di soggetto
kfold_vi_r2 <- function(df, y_var, x_vars, Kfold = 10) {
  df <- df %>% select(all_of(c("subject_numeric", y_var, x_vars))) %>% drop_na()
  stopifnot(nrow(df) >= Kfold)
  # assegna fold stratificando grossolanamente per outcome (opzionale)
  folds <- sample(rep(1:Kfold, length.out = nrow(df)))
  r2s <- numeric(Kfold)

  for (k in 1:Kfold) {
    tr <- df[folds != k, , drop = FALSE]
    te <- df[folds == k, , drop = FALSE]
    Xtr <- as.matrix(tr[, x_vars, drop = FALSE])
    Xte <- as.matrix(te[, x_vars, drop = FALSE])
    ytr <- tr[[y_var]]
    yte <- te[[y_var]]
    # se nessun predittore (K=0) evita VI e usa solo media (R^2=0 per definizione)
    if (ncol(Xtr) == 0) {
      yhat <- rep(mean(ytr), length(yte))
    } else {
      yhat <- vi_fit_predict(Xtr, ytr, Xte)$y_hat
    }
    r2s[k] <- r2_oos(yte, yhat)
  }
  c(mean = mean(r2s), sd = sd(r2s))
}

# -- specifica dei set di predittori
x_base <- c("z_baseline_pid5", "z_baseline_dass_stress", "sex_numeric")
x_ema_pre <- c("z_ema_mean_pre", "z_ema_sd_pre")
x_ema_post <- c("z_ema_mean_post", "z_ema_sd_post")

# -- esegui per eta_pre
res_pre <- tibble::tibble(
  Model = c("Baseline", "EMA-only (pre)", "Combined (pre)"),
  R2_mean = NA_real_,
  R2_sd = NA_real_
)

m1 <- kfold_vi_r2(prediction_data, "eta_pre", x_base, Kfold = 10)
m2 <- kfold_vi_r2(prediction_data, "eta_pre", x_ema_pre, Kfold = 10)
m3 <- kfold_vi_r2(prediction_data, "eta_pre", c(x_base, x_ema_pre), Kfold = 10)
res_pre$R2_mean <- c(m1["mean"], m2["mean"], m3["mean"])
res_pre$R2_sd <- c(m1["sd"], m2["sd"], m3["sd"])
res_pre <- res_pre %>%
  mutate(Delta_vs_Baseline = R2_mean - R2_mean[Model == "Baseline"])

cat("\n== R^2 out-of-sample su eta_pre (k-fold VI) ==\n")
print(res_pre)

# -- esegui per eta_post
res_post <- tibble::tibble(
  Model = c("Baseline", "EMA-only (post)", "Combined (post)"),
  R2_mean = NA_real_,
  R2_sd = NA_real_
)

m1p <- kfold_vi_r2(prediction_data, "eta_post", x_base, Kfold = 10)
m2p <- kfold_vi_r2(prediction_data, "eta_post", x_ema_post, Kfold = 10)
m3p <- kfold_vi_r2(
  prediction_data,
  "eta_post",
  c(x_base, x_ema_post),
  Kfold = 10
)
res_post$R2_mean <- c(m1p["mean"], m2p["mean"], m3p["mean"])
res_post$R2_sd <- c(m1p["sd"], m2p["sd"], m3p["sd"])
res_post <- res_post %>%
  mutate(Delta_vs_Baseline = R2_mean - R2_mean[Model == "Baseline"])

cat("\n== R^2 out-of-sample su eta_post (k-fold VI) ==\n")
print(res_post)

# -- piccolo riassunto unico
summary_oos <- dplyr::bind_rows(
  res_pre %>% mutate(Outcome = "eta_pre"),
  res_post %>% mutate(Outcome = "eta_post")
) %>%
  relocate(Outcome)
cat("\n== Riepilogo complessivo R^2 OOS ==\n")
print(summary_oos)

#######################################################
