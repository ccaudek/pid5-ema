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

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("sd", "stats")
conflict_prefer("var", "stats")
conflict_prefer("mad", "stats")
conflict_prefer("match", "base")
conflict_prefer("lag", "dplyr")

set.seed(20250917)

z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m)
  (x - m) / s
}

# ============================================================
# Data preparation for: fragility_ema_vs_baseline.stan
# - include ONLY subjects with both PRE and POST
# - build ordinal items, EMA-long, baseline matrix
# ============================================================

# Carichiamo il dataset completo
d <- rio::import(
  here::here(
    "data",
    "processed",
    "ema_plus_baseline_exam_tagged_20250916.csv"
  )
)

# z-score robusto con fallback se sd==0 o NA
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m)
  (x - m) / s
}

prepare_fragility_data <- function(d, min_obs_per_period = 1) {
  d <- d %>% mutate(day = as.Date(day))

  # 1) Item-level + ordinalizzazioni
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
      # 5 dimensioni PID-5 EMA (ripetute)
      pid5_negative_affectivity,
      pid5_detachment,
      pid5_antagonism,
      pid5_disinhibition,
      pid5_psychoticism,
      # 5 dimensioni PID-5 baseline (una volta)
      pid5_negative_affect_baseline,
      pid5_detachment_baseline,
      pid5_antagonism_baseline,
      pid5_disinhibition_baseline,
      pid5_psychoticism_baseline
    ) %>%
    mutate(
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

  # 2) prima mappa (provvisoria)
  subject_map0 <- item_data %>%
    distinct(user_id) %>%
    arrange(user_id) %>%
    mutate(subject_numeric0 = row_number())

  obs0 <- item_data %>%
    left_join(subject_map0, by = "user_id") %>%
    arrange(subject_numeric0, day, hour)

  # 3) tieni solo chi ha almeno un pre e un post
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

  # 4) re-index dopo il filtro (contiguo 1..J)
  subject_map1 <- obs1 %>%
    distinct(user_id) %>%
    arrange(user_id) %>%
    mutate(subject_numeric = row_number())

  obs2 <- obs1 %>%
    select(-subject_numeric0) %>%
    left_join(subject_map1, by = "user_id") %>%
    arrange(subject_numeric, day, hour)

  # 5) osservazioni valide
  obs_data0 <- obs2 %>%
    filter(if_any(
      c(happy_ord, sad_ord, satisfied_ord, angry_ord),
      ~ !is.na(.)
    )) %>%
    arrange(subject_numeric, day, hour)

  # 5b) imponi minimo di osservazioni per PRE e POST
  keep_subj <- obs_data0 %>%
    group_by(subject_numeric, period_numeric) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = period_numeric,
      values_from = n_obs,
      names_prefix = "p",
      values_fill = 0
    ) %>%
    mutate(keep = (p2 >= min_obs_per_period) & (p3 >= min_obs_per_period)) %>%
    filter(keep) %>%
    pull(subject_numeric)

  obs_data1 <- obs_data0 %>% filter(subject_numeric %in% keep_subj)

  # *** RE-INDEX FINALE DOPO IL MINIMO PER PERIODO ***
  final_subject_map <- obs_data1 %>%
    distinct(user_id) %>%
    arrange(user_id) %>%
    mutate(subject_numeric_final = row_number())

  obs_data <- obs_data1 %>%
    left_join(final_subject_map, by = "user_id") %>%
    select(-subject_numeric) %>%
    rename(subject_numeric = subject_numeric_final) %>%
    arrange(subject_numeric, day, hour) %>%
    mutate(obs_id = row_number())

  # 6) long per item
  items_long <- obs_data %>%
    select(obs_id, happy_ord, sad_ord, satisfied_ord, angry_ord) %>%
    tidyr::pivot_longer(
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

  # 7) baseline 5D per soggetto (z-score + imputazione 0)
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

  # 8) EMA 5D ripetuti: z-score per osservazione
  ema_mat <- obs_data %>%
    transmute(
      obs_id,
      subj = subject_numeric,
      naff_e = pid5_negative_affectivity,
      det_e = pid5_detachment,
      ant_e = pid5_antagonism,
      dis_e = pid5_disinhibition,
      psy_e = pid5_psychoticism
    ) %>%
    mutate(across(c(naff_e, det_e, ant_e, dis_e, psy_e), z_))

  ema_long <- ema_mat %>%
    tidyr::pivot_longer(
      cols = c(naff_e, det_e, ant_e, dis_e, psy_e),
      names_to = "dim",
      values_to = "value"
    ) %>%
    filter(is.finite(value)) %>%
    mutate(
      dim_id = match(dim, c("naff_e", "det_e", "ant_e", "dis_e", "psy_e"))
    ) %>%
    arrange(obs_id, dim_id)

  M_ema <- nrow(ema_long)
  ema_val <- ema_long$value
  ema_dim <- ema_long$dim_id
  ema_obs <- ema_long$obs_id

  # 9) sanity checks (ora devono passare)
  I <- dplyr::n_distinct(obs_data$subject_numeric)
  N_obs <- nrow(obs_data)

  stopifnot(identical(sort(unique(obs_data$subject_numeric)), seq_len(I)))
  stopifnot(max(obs_data$subject_numeric) == I)
  stopifnot(all(ema_obs >= 1 & ema_obs <= N_obs))
  stopifnot(nrow(X_base) == I, ncol(X_base) == 5)

  # 10) lista Stan
  stan_data <- list(
    I = I,
    N_obs = N_obs,
    K = 4L,
    P = 3L,
    D = 5L,
    N_items = nrow(items_long),
    y_item = as.integer(items_long$response),
    item_id = as.integer(items_long$item_id),
    obs_id = as.integer(items_long$obs_id),
    subject = as.integer(obs_data$subject_numeric),
    period = as.integer(obs_data$period_numeric),
    M_ema = as.integer(M_ema),
    ema_val = as.numeric(ema_val),
    ema_dim = as.integer(ema_dim),
    ema_obs = as.integer(ema_obs),
    X_base = X_base,
    use_ema = 1.0
  )

  list(
    stan_data = stan_data,
    obs_data = obs_data,
    subject_map = final_subject_map %>%
      rename(subject_numeric = subject_numeric_final),
    base_dims = base_dims,
    items_long = items_long
  )
}

# 1) Prepara i dati
frag_data <- prepare_fragility_data(d, min_obs_per_period = 1)

with(frag_data$stan_data, {
  cat(
    "I =",
    I,
    " | max(subject) =",
    max(subject),
    " | unique subjects =",
    length(unique(subject)),
    "\n"
  )
  cat("N_obs =", N_obs, " | max(obs_id) =", max(obs_id), "\n")
})

# 2) Compila il modello
stan_file <- here::here(
  "scripts",
  "02_stress_reactivity",
  "fragility_ema_vs_baseline.stan"
)
model_frag <- cmdstan_model(stan_file)

# 3) Fit A: baseline-only (use_ema = 0)
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

# 4) Fit B: baseline + EMA (use_ema = 1)
data_B <- frag_data$stan_data
data_B$use_ema <- 1.0

fit_B <- model_frag$variational(
  data = data_B,
  seed = 20250912,
  algorithm = "meanfield", # puoi rifare con fullrank quando vuoi
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

# 6) LOO su log-lik della fragilità (per soggetto)
# Extract item-wise log-lik matrices: draws x N_items
ll_items_A <- fit_A$draws("log_lik_item", format = "matrix")
ll_items_B <- fit_B$draws("log_lik_item", format = "matrix")

loo_items_A <- loo::loo(ll_items_A, moment_match = TRUE)
loo_items_B <- loo::loo(ll_items_B, moment_match = TRUE)

# Name models explicitly to avoid confusion
cmp_items <- loo::loo_compare(list(A = loo_items_A, B = loo_items_B))
print(cmp_items)
attr(cmp_items, "ics") # to double-check which row is which


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
