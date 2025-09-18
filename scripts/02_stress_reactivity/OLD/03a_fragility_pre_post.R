suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tidyselect)
  library(lubridate)
  library(loo)
  library(conflicted)
})


# ============================================================
# R: prepara i dati per fragility_ema_vs_baseline.stan
# - NA EMA viene stimata *nel modello* dai 4 item ordinali
# - A Stan passiamo solo le 4 EMA osservate: det/ant/dis/psy
# - Le 5 baseline (una riga per soggetto) entrano come X_base
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(lubridate)
  library(rio)
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(Matrix)
  library(conflicted)
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("sd", "stats")
conflict_prefer("var", "stats")
conflict_prefer("mad", "stats")
conflict_prefer("match", "base")
conflict_prefer("lag", "dplyr")

set.seed(20250917)

# z-score "sicuro"
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m)
  (x - m) / s
}

# ---------- Carica il dataset "clean" ----------
d <- rio::import(
  here::here("data", "processed", "ema_plus_scales_cleaned.csv")
) |>
  dplyr::rename(
    # rinomina per coerenza con lo Stan (baseline a 5 dimensioni)
    pid5_negative_affect_baseline = domain_negative_affect_baseline,
    pid5_detachment_baseline = domain_detachment_baseline,
    pid5_antagonism_baseline = domain_antagonism_baseline,
    pid5_disinhibition_baseline = domain_disinhibition_baseline,
    pid5_psychoticism_baseline = domain_psychoticism_baseline
  )

# ---------- Funzione di preparazione ----------
prepare_fragility_data <- function(d, min_obs_per_period = 1) {
  d <- d %>% mutate(day = as.Date(day))

  # 1) Item ordinali (happy/sad/satisfied/angry) -> 1..7
  item_data <- d %>%
    filter(!is.na(exam_period)) %>%
    select(
      user_id,
      day,
      hour,
      exam_period,
      happy,
      sad,
      satisfied,
      angry,
      # EMA osservate (4 dimensioni: det/ant/dis/psy)
      pid5_detachment,
      pid5_antagonism,
      pid5_disinhibition,
      pid5_psychoticism,
      # baseline 5D (una per soggetto)
      pid5_negative_affect_baseline,
      pid5_detachment_baseline,
      pid5_antagonism_baseline,
      pid5_disinhibition_baseline,
      pid5_psychoticism_baseline
    ) %>%
    mutate(
      # mappa 0..100 -> 1..7
      happy_ord = pmax(1L, pmin(7L, as.integer(round(happy / 100 * 6) + 1))),
      sad_ord = pmax(1L, pmin(7L, as.integer(round(sad / 100 * 6) + 1))),
      satisfied_ord = pmax(
        1L,
        pmin(7L, as.integer(round(satisfied / 100 * 6) + 1))
      ),
      angry_ord = pmax(1L, pmin(7L, as.integer(round(angry / 100 * 6) + 1))),
      period_numeric = case_when(
        exam_period == "baseline" ~ 1L,
        exam_period == "pre_exam" ~ 2L,
        exam_period == "post_exam" ~ 3L,
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(!is.na(period_numeric)) %>%
    arrange(user_id, day, hour)

  # 2) indicizza soggetti, tieni solo chi ha almeno un pre e un post
  subj0 <- item_data %>%
    distinct(user_id) %>%
    arrange(user_id) %>%
    mutate(subject_numeric0 = row_number())

  obs0 <- item_data %>%
    left_join(subj0, by = "user_id") %>%
    arrange(subject_numeric0, day, hour)

  keep_ids <- obs0 %>%
    group_by(user_id) %>%
    summarise(
      has_pre = any(period_numeric == 2L),
      has_post = any(period_numeric == 3L),
      .groups = "drop"
    ) %>%
    filter(has_pre & has_post) %>%
    pull(user_id)

  obs1 <- obs0 %>% filter(user_id %in% keep_ids)

  # 3) re-index finale soggetti
  subj_map <- obs1 %>%
    distinct(user_id) %>%
    arrange(user_id) %>%
    mutate(subject_numeric = row_number())

  obs2 <- obs1 %>%
    select(-subject_numeric0) %>%
    left_join(subj_map, by = "user_id") %>%
    arrange(subject_numeric, day, hour)

  # 4) filtra osservazioni con almeno 1 item presente
  obs_items <- obs2 %>%
    filter(if_any(
      c(happy_ord, sad_ord, satisfied_ord, angry_ord),
      ~ !is.na(.)
    )) %>%
    arrange(subject_numeric, day, hour)

  # 5) vincolo minimo di osservazioni per PRE e POST
  keep_subj <- obs_items %>%
    count(subject_numeric, period_numeric, name = "n_obs") %>%
    pivot_wider(
      names_from = period_numeric,
      values_from = n_obs,
      names_prefix = "p",
      values_fill = 0
    ) %>%
    mutate(keep = (p2 >= min_obs_per_period) & (p3 >= min_obs_per_period)) %>%
    filter(keep) %>%
    pull(subject_numeric)

  obs_data <- obs_items %>%
    filter(subject_numeric %in% keep_subj) %>%
    arrange(subject_numeric, day, hour) %>%
    mutate(obs_id = row_number())

  # 6) long degli item (solo non-NA)
  items_long <- obs_data %>%
    select(obs_id, happy_ord, sad_ord, satisfied_ord, angry_ord) %>%
    pivot_longer(
      ends_with("_ord"),
      names_to = "item_name",
      values_to = "response"
    ) %>%
    filter(!is.na(response)) %>%
    mutate(
      item_id = case_when(
        item_name == "happy_ord" ~ 1L,
        item_name == "sad_ord" ~ 2L,
        item_name == "satisfied_ord" ~ 3L,
        item_name == "angry_ord" ~ 4L
      ),
      response = as.integer(response)
    ) %>%
    arrange(obs_id, item_id)

  # 7) Baseline 5D per soggetto, z-score + imputazione a 0 se NA
  base_dims <- obs_data %>%
    group_by(subject_numeric, user_id) %>%
    summarise(
      naff_b = first(pid5_negative_affect_baseline),
      det_b = first(pid5_detachment_baseline),
      ant_b = first(pid5_antagonism_baseline),
      dis_b = first(pid5_disinhibition_baseline),
      psy_b = first(pid5_psychoticism_baseline),
      .groups = "drop"
    ) %>%
    arrange(subject_numeric)

  X_base_df <- base_dims %>%
    mutate(across(
      c(naff_b, det_b, ant_b, dis_b, psy_b),
      z_,
      .names = "z_{.col}"
    )) %>%
    mutate(across(starts_with("z_"), ~ coalesce(.x, 0))) %>%
    select(starts_with("z_"))

  X_base <- as.matrix(X_base_df)
  colnames(X_base) <- c("z_naff_b", "z_det_b", "z_ant_b", "z_dis_b", "z_psy_b")

  # 8) EMA osservate (SOLO 4 dimensioni: det/ant/dis/psy) per osservazione
  ema_mat <- obs_data %>%
    transmute(
      obs_id,
      subj = subject_numeric,
      det_e = z_(pid5_detachment),
      ant_e = z_(pid5_antagonism),
      dis_e = z_(pid5_disinhibition),
      psy_e = z_(pid5_psychoticism)
    )

  ema_long <- ema_mat %>%
    pivot_longer(
      cols = c(det_e, ant_e, dis_e, psy_e),
      names_to = "dim",
      values_to = "value"
    ) %>%
    filter(is.finite(value)) %>%
    mutate(
      # Mappa alle dimensioni Stan: 1=NA (latente dagli item), 2..5 = det/ant/dis/psy
      dim_id = dplyr::recode(
        dim,
        det_e = 2L,
        ant_e = 3L,
        dis_e = 4L,
        psy_e = 5L
      ),
      dim_id = as.integer(dim_id)
    ) %>%
    arrange(obs_id, dim_id)

  # 9) controlli
  I <- n_distinct(obs_data$subject_numeric)
  N_obs <- nrow(obs_data)
  stopifnot(identical(sort(unique(obs_data$subject_numeric)), seq_len(I)))
  stopifnot(max(obs_data$subject_numeric) == I)
  stopifnot(nrow(X_base) == I, ncol(X_base) == 5)

  # 10) dati per Stan
  stan_data <- list(
    I = I,
    N_obs = N_obs,
    K = 4L, # happy, sad, satisfied, angry
    P = 3L, # baseline, pre, post
    D = 5L, # 5 dimensioni EMA a livello soggetto (1=NA latente, 2..5 osservate)
    N_items = nrow(items_long),
    y_item = as.integer(items_long$response),
    item_id = as.integer(items_long$item_id),
    obs_id = as.integer(items_long$obs_id),
    subject = as.integer(obs_data$subject_numeric),
    period = as.integer(obs_data$period_numeric),
    M_ema = nrow(ema_long),
    ema_val = as.numeric(ema_long$value),
    ema_dim = as.integer(ema_long$dim_id), # solo {2,3,4,5}
    ema_obs = as.integer(ema_long$obs_id),
    X_base = X_base,
    use_ema = 1.0
  )

  list(
    stan_data = stan_data,
    obs_data = obs_data,
    subject_map = subj_map,
    items_long = items_long
  )
}

# ---------- Prepara e fit ----------
frag_data <- prepare_fragility_data(d, min_obs_per_period = 1)

with(frag_data$stan_data, {
  cat("I =", I, "| N_obs =", N_obs, "| N_items =", N_items, "\n")
})

stan_file <- here::here(
  "scripts",
  "02_stress_reactivity",
  "fragility_ema_vs_baseline.stan"
)
model_frag <- cmdstan_model(stan_file)

# Fit A: baseline-only
data_A <- frag_data$stan_data
data_A$use_ema <- 0.0
fit_A <- model_frag$variational(
  data = data_A,
  seed = 20250912,
  algorithm = "meanfield",
  elbo_samples = 100,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 200
)

# Fit B: baseline + EMA (inclusa NA latente dagli item)
data_B <- frag_data$stan_data
data_B$use_ema <- 1.0
fit_B <- model_frag$variational(
  data = data_B,
  seed = 20250912,
  algorithm = "meanfield",
  elbo_samples = 100,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 200
)


# 5) Estrai Bayesian R2 della regressione della fragilità
R2_A <- as_draws_df(fit_A$draws("R2_frag"))$R2_frag
R2_B <- as_draws_df(fit_B$draws("R2_frag"))$R2_frag
cat("Bayesian R2 – baseline-only:     ", round(mean(R2_A), 3), "\n")
cat("Bayesian R2 – baseline + EMA:    ", round(mean(R2_B), 3), "\n")
cat(
  "ΔR2 (B − A):                     ",
  round(mean(R2_B - R2_A), 3),
  "  (95% CrI ≈ [",
  round(quantile(R2_B - R2_A, c(.025, .975)), 3),
  "] )\n",
  sep = ""
)

# =======================
# LOO: test cruciale EMA
# =======================

# 1) Estrarre log-lik per osservazione (beep)
ll_A <- fit_A$draws("log_lik_obs", format = "matrix") # draws x N_obs
ll_B <- fit_B$draws("log_lik_obs", format = "matrix")

# controllo dimensioni coerenti
stopifnot(ncol(ll_A) == ncol(ll_B))
N_obs <- ncol(ll_A)

# 2) LOO con moment matching (riduce impatto k alti)
loo_A <- loo::loo(ll_A, moment_match = TRUE)
loo_B <- loo::loo(ll_B, moment_match = TRUE)

cat("\n=== Pareto-k (beep level) ===\n")
print(loo::pareto_k_table(loo_A))
print(loo::pareto_k_table(loo_B))

# 3) Confronto modelli (beep level)
cmp_obs <- loo::loo_compare(list(B = loo_B, A = loo_A))
cat("\n=== LOO comparison (beep level) ===\n")
print(cmp_obs)

# Δelpd = elpd(B) - elpd(A)
delta_elpd_obs <- as.numeric(cmp_obs["A", "elpd_diff"]) * -1
per_obs <- delta_elpd_obs / N_obs
pct_gain <- (exp(per_obs) - 1) * 100

cat(sprintf(
  "\nΔelpd (B−A): %.1f  | per-beep: %.4f  | mean %% gain per beep: %.2f%%\n",
  delta_elpd_obs,
  per_obs,
  pct_gain
))

# 4) (Opzionale) LOO aggregato per soggetto: somma log-lik di tutti i beep dello stesso soggetto
#    È utile come analisi di robustezza (meno punti = minore varianza del PSIS)
subj_idx <- frag_data$stan_data$subject # length = N_obs, valori 1..I
I <- max(subj_idx)

# Matrice di raggruppamento JxN (J=I soggetti). Sparse per efficienza.
G <- sparseMatrix(i = subj_idx, j = seq_len(N_obs), x = 1, dims = c(I, N_obs)) # I x N_obs

# draws x N_obs  %*%  t(G)  ->  draws x I (somma per soggetto)
ll_A_subj <- ll_A %*% t(as.matrix(G))
ll_B_subj <- ll_B %*% t(as.matrix(G))

loo_A_subj <- loo::loo(ll_A_subj, moment_match = TRUE)
loo_B_subj <- loo::loo(ll_B_subj, moment_match = TRUE)

cat("\n=== Pareto-k (subject level) ===\n")
print(loo::pareto_k_table(loo_A_subj))
print(loo::pareto_k_table(loo_B_subj))

cmp_subj <- loo::loo_compare(list(B = loo_B_subj, A = loo_A_subj))
cat("\n=== LOO comparison (subject level) ===\n")
print(cmp_subj)

delta_elpd_subj <- as.numeric(cmp_subj["A", "elpd_diff"]) * -1
per_subj <- delta_elpd_subj / I
pct_gain_subj <- (exp(per_subj) - 1) * 100

cat(sprintf(
  "\nΔelpd (B−A): %.1f  | per-subject: %.4f  | mean %% gain per subject: %.2f%%\n",
  delta_elpd_subj,
  per_subj,
  pct_gain_subj
))

# 5) Riassunto finale — “decisione” per la domanda di ricerca
winner_obs <- rownames(cmp_obs)[1]
winner_subj <- rownames(cmp_subj)[1]
cat("\n=== Verdict ===\n")
cat(sprintf(
  "Beep level:   %s ha migliore ELPD (Δelpd = %.1f).\n",
  winner_obs,
  delta_elpd_obs
))
cat(sprintf(
  "Subject level:%s ha migliore ELPD (Δelpd = %.1f).\n",
  winner_subj,
  delta_elpd_subj
))

# Nota interpretativa (breve):
# - Se B è in testa e Δelpd è > ~4-10 (regola empirica) con SE non enorme, c’è
#   evidenza che l’aggiunta delle EMA migliora la predizione out-of-sample.
# - I Pareto-k > 0.7 in molti punti suggeriscono più “instabilità”; qui abbiamo
#   usato moment matching. In caso di troppi k>1, valuta k-fold CV.

#######
ll_subj_A <- fit_A$draws("log_lik_subj", format = "matrix")
ll_subj_B <- fit_B$draws("log_lik_subj", format = "matrix")

loo_subj_A <- loo::loo(ll_subj_A, moment_match = TRUE)
loo_subj_B <- loo::loo(ll_subj_B, moment_match = TRUE)

loo::pareto_k_table(loo_subj_A)
loo::pareto_k_table(loo_subj_B)

cmp_subj <- loo::loo_compare(list(A = loo_subj_A, B = loo_subj_B))
print(cmp_subj)


# # Estrai log-lik per osservazione (beep)
# ll_obs_A <- fit_A$draws("log_lik_obs", format = "matrix")
# ll_obs_B <- fit_B$draws("log_lik_obs", format = "matrix")
#
# loo_obs_A <- loo::loo(ll_obs_A, moment_match = TRUE)
# loo_obs_B <- loo::loo(ll_obs_B, moment_match = TRUE)
#
# print(loo::pareto_k_table(loo_obs_A))
# print(loo::pareto_k_table(loo_obs_B))
#
# cmp_obs <- loo::loo_compare(list(A = loo_obs_A, B = loo_obs_B))
# print(cmp_obs)

# 7) Coefficienti interpretabili
summ_A <- fit_A$summary(
  variables = c("a_frag", "b_base", "sigma_diff", "gamma_pre", "gamma_post")
)
summ_B <- fit_B$summary(
  variables = c(
    "a_frag",
    "b_base",
    "c_ema",
    "sigma_diff",
    "gamma_pre",
    "gamma_post"
  )
)
print(summ_A, n = Inf)
print(summ_B, n = Inf)

# In particolare, in B guarda c_ema[1..5]: quanto ciascuna dimensione EMA
# aggiunge predittività sulla fragilità al netto delle baseline b_base[1..5].

diff_draws <- as_draws_df(fit_A$draws("diff_frag"))
diff_mean <- diff_draws %>%
  summarise(across(everything(), mean)) %>%
  t() %>%
  as.vector()
c(mean = mean(diff_mean), sd = sd(diff_mean), iqr = IQR(diff_mean))

print(fit_B$summary(variables = "c_ema"), n = Inf)


loo::pareto_k_table(loo_items_A)
loo::pareto_k_table(loo_items_B)
loo::loo(ll_items_A, moment_match = TRUE)
loo::loo(ll_items_B, moment_match = TRUE)


delta_elpd <- as.numeric(cmp_items["A", "elpd_diff"]) * -1 # B - A
N_items <- ncol(ll_items_A)
per_item <- delta_elpd / N_items
rel_gain <- exp(per_item) - 1
c(delta_elpd = delta_elpd, per_item = per_item, pct_gain = 100 * rel_gain)


obs_data <- frag_data$obs_data # has subject_numeric, period_numeric, and items

# Count per subject × period (2=pre, 3=post)
counts_prepost <- obs_data %>%
  count(subject_numeric, period_numeric, name = "n_obs") %>%
  filter(period_numeric %in% c(2L, 3L)) %>%
  pivot_wider(
    names_from = period_numeric,
    values_from = n_obs,
    names_prefix = "p",
    values_fill = 0
  ) %>%
  rename(pre_exam = p2, post_exam = p3)

min_obs_by_subj <- counts_prepost %>%
  mutate(min_obs_per_period = pmin(pre_exam, post_exam))

# Summary + histogram
summary(min_obs_by_subj$min_obs_per_period)
table(min_obs_by_subj$min_obs_per_period)

threshold <- 2
ggplot(min_obs_by_subj, aes(x = min_obs_per_period)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
  geom_vline(xintercept = threshold, linetype = "dashed") +
  scale_x_continuous(
    breaks = 0:max(min_obs_by_subj$min_obs_per_period, na.rm = TRUE)
  ) +
  labs(
    title = "Min # of valid EMA observations per subject (pre vs post)",
    x = "min_obs_per_period = min(#pre_exam, #post_exam)",
    y = "Number of subjects"
  ) +
  theme_minimal()
