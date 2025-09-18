# ============================================================
# Fragilità psicologica — soggetto-level (PRE-or-POST vs BASELINE)
# A: baseline-only (PID-5 baseline)
# B: baseline + EMA (residui delle medie EMA + dinamiche EMA al baseline)
# Likelihood: Student-t; confronto: PSIS-LOO (k<0.7 atteso)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(tibble)
  library(ggplot2)
  library(conflicted)
  library(rio)
  library(here)
})

# ------------------------------------------------------------
# Preferenze conflitti
# ------------------------------------------------------------
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("sd", "stats")
conflict_prefer("var", "stats")
conflict_prefer("mad", "stats")
conflict_prefer("match", "base")

# ------------------------------------------------------------
# Parametri e path
# ------------------------------------------------------------
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
stan_file_path <- here::here(
  "scripts",
  "02_stress_reactivity",
  "fragility_subject_rhs.stan"
)

MIN_ITEMS_PER_BEEP <- 2L # minimo item presenti per calcolare NA_beep
TRIM_PROP <- 0.10 # trimmed mean (10%)
MIN_OBS_BASE <- 3L # minimo beep baseline
MIN_OBS_EXAM <- 1L # minimo beep PRE o POST
USE_GENDER <- 0L # opzionale (0=off)
USE_EMA_MEAN_RES <- 1L # usa i residui delle medie EMA (incremental validity)
USE_EMA_DYN <- 1L # usa le dinamiche EMA (SD, MSSD, AR1)
SEED <- 20250918

# ------------------------------------------------------------
# Utility
# ------------------------------------------------------------
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(x - m)
  (x - m) / s
}
find_first <- function(df, candidates) {
  nm <- names(df)
  hit <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hit) == 0) NULL else hit[1]
}
trim_mean <- function(x, trim = 0.1) mean(x, trim = trim, na.rm = TRUE)

# ------------------------------------------------------------
# 1) Carica dati + normalizza nomi PID-5 baseline
# ------------------------------------------------------------
d <- rio::import(data_path) %>%
  dplyr::rename(
    pid5_negative_affect_baseline = any_of(c(
      "domain_negative_affect_baseline",
      "pid5_negative_affect_baseline"
    )),
    pid5_detachment_baseline = any_of(c(
      "domain_detachment_baseline",
      "pid5_detachment_baseline"
    )),
    pid5_antagonism_baseline = any_of(c(
      "domain_antagonism_baseline",
      "pid5_antagonism_baseline"
    )),
    pid5_disinhibition_baseline = any_of(c(
      "domain_disinhibition_baseline",
      "pid5_disinhibition_baseline"
    )),
    pid5_psychoticism_baseline = any_of(c(
      "domain_psychoticism_baseline",
      "pid5_psychoticism_baseline"
    ))
  )

# ------------------------------------------------------------
# 2) ID soggetto e periodo -> subj, per
# ------------------------------------------------------------
subj_col <- d %>%
  select(any_of(c("user_id", "id", "subject_id", "participant_id"))) %>%
  names() %>%
  .[1]
if (is.na(subj_col))
  stop(
    "Nessuna colonna soggetto trovata (user_id/id/subject_id/participant_id)."
  )

d$`__subj__` <- as.integer(factor(d[[subj_col]]))
per_col <- d %>%
  select(any_of(c("exam_period", "period", "phase"))) %>%
  names() %>%
  .[1]
if (!is.na(per_col)) {
  pr <- d[[per_col]]
  per <- if (is.numeric(pr)) as.integer(pr) else {
    lv <- tolower(as.character(pr))
    dplyr::case_when(
      lv %in% c("baseline", "base", "t0", "bl") ~ 1L,
      lv %in% c("pre", "pre_exam", "pre-exam", "preexam") ~ 2L,
      lv %in% c("post", "post_exam", "post-exam", "postexam") ~ 3L,
      TRUE ~ 1L
    )
  }
} else per <- rep(1L, nrow(d))
d$`__per__` <- as.integer(per)

d <- d %>% mutate(subj = `__subj__`, per = `__per__`)

# ------------------------------------------------------------
# 3) Genere -> female (0/1) per soggetto
# ------------------------------------------------------------
sex_col <- d %>%
  select(any_of(c("female", "sex", "gender", "sesso"))) %>%
  names() %>%
  .[1]
if (!is.na(sex_col)) {
  sx <- d[[sex_col]]
  if (is.numeric(sx)) female_by_row <- as.integer(sx != 0) else {
    lv <- tolower(as.character(sx))
    female_by_row <- as.integer(
      lv %in% c("f", "female", "donna", "femmina", "woman", "1", "true")
    )
  }
} else {
  warning("Colonna sesso assente: imposto female=1 per tutti (placeholder).")
  female_by_row <- rep(1L, nrow(d))
}
fem_by_subj <- tapply(female_by_row, d$subj, function(v) mean(v, na.rm = TRUE))
fem_by_subj[is.na(fem_by_subj)] <- 1
female_vec_all <- as.integer(fem_by_subj >= 0.5)

# ------------------------------------------------------------
# 4) Item grezzi
# ------------------------------------------------------------
get_item <- function(name_candidates) {
  nm <- find_first(d, name_candidates)
  if (is.null(nm)) return(NA_real_)
  x <- d[[nm]]
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    x[x == ""] <- NA_character_
    x <- readr::parse_number(
      x,
      locale = readr::locale(decimal_mark = ".", grouping_mark = ",")
    )
  }
  as.numeric(x)
}
happy <- get_item(c("happy", "happiness", "felice", "contento"))
sad <- get_item(c("sad", "sadness", "triste"))
satisfied <- get_item(c(
  "satisfied",
  "satisfaction",
  "soddisfatto",
  "soddisfazione"
))
angry <- get_item(c("angry", "anger", "arrabbiato", "rabbia"))

# ------------------------------------------------------------
# 5) NA_beep robusto ai missing (z-score sui disponibili; media dei segni disponibili)
# ------------------------------------------------------------
d_items_all <- d %>%
  mutate(happy = happy, sad = sad, satisfied = satisfied, angry = angry) %>%
  mutate(
    z_happy = z_(happy),
    z_sad = z_(sad),
    z_sat = z_(satisfied),
    z_angry = z_(angry)
  )

n_present <- function(h, s, t, a)
  rowSums(cbind(is.finite(h), is.finite(s), is.finite(t), is.finite(a)))

d_items <- d_items_all %>%
  mutate(
    n_items = n_present(z_happy, z_sad, z_sat, z_angry),
    NA_beep = rowMeans(cbind(-z_happy, +z_sad, -z_sat, +z_angry), na.rm = TRUE)
  ) %>%
  filter(is.finite(NA_beep), n_items >= MIN_ITEMS_PER_BEEP) %>%
  select(subj, per, NA_beep)

message(sprintf(
  "Beep utili: %d | soggetti unici: %d",
  nrow(d_items),
  dplyr::n_distinct(d_items$subj)
))
message("Distribuzione per periodi (1=baseline, 2=pre, 3=post):")
print(table(d_items$per))

# ------------------------------------------------------------
# 6) Outcome soggetto-level: y = NA_exam - NA_baseline (exam = PRE se disponibile, altrimenti POST)
# ------------------------------------------------------------
base_agg <- d_items %>%
  filter(per == 1L) %>%
  group_by(subj) %>%
  summarise(
    n_base = n(),
    NA_base = mean(NA_beep, trim = TRIM_PROP, na.rm = TRUE),
    .groups = "drop"
  )

pre_agg <- d_items %>%
  filter(per == 2L) %>%
  group_by(subj) %>%
  summarise(
    n_pre = n(),
    NA_pre = mean(NA_beep, trim = TRIM_PROP, na.rm = TRUE),
    .groups = "drop"
  )

post_agg <- d_items %>%
  filter(per == 3L) %>%
  group_by(subj) %>%
  summarise(
    n_post = n(),
    NA_post = mean(NA_beep, trim = TRIM_PROP, na.rm = TRUE),
    .groups = "drop"
  )

exam_df <- full_join(pre_agg, post_agg, by = "subj") %>%
  mutate(
    has_pre = !is.na(NA_pre) & n_pre >= MIN_OBS_EXAM,
    has_post = !is.na(NA_post) & n_post >= MIN_OBS_EXAM,
    NA_exam = dplyr::case_when(
      has_pre ~ NA_pre,
      !has_pre & has_post ~ NA_post,
      TRUE ~ NA_real_
    ),
    n_exam = dplyr::case_when(
      has_pre ~ n_pre,
      !has_pre & has_post ~ n_post,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(subj, NA_exam, n_exam)

agg_subj <- base_agg %>%
  inner_join(exam_df, by = "subj") %>%
  filter(
    n_base >= MIN_OBS_BASE,
    is.finite(NA_base),
    is.finite(NA_exam),
    n_exam >= MIN_OBS_EXAM
  ) %>%
  mutate(y = NA_exam - NA_base) %>%
  arrange(subj)

if (nrow(agg_subj) == 0)
  stop("Nessun soggetto con baseline e PRE/POST sufficienti.")

ids <- agg_subj$subj
I <- nrow(agg_subj)
y <- agg_subj$y

message(sprintf(
  "Soggetti inclusi: %d  |  beep baseline mediano: %d  |  beep exam mediano: %d",
  I,
  stats::median(agg_subj$n_base),
  stats::median(agg_subj$n_exam)
))

# ------------------------------------------------------------
# 7) X_base (I x 5) — z-score per colonna + imputazione a 0 dei NA
# ------------------------------------------------------------
base_cols <- c(
  "pid5_negative_affect_baseline",
  "pid5_detachment_baseline",
  "pid5_antagonism_baseline",
  "pid5_disinhibition_baseline",
  "pid5_psychoticism_baseline"
)

if (!all(base_cols %in% names(d))) {
  stop(
    "Mancano colonne baseline attese: ",
    paste(setdiff(base_cols, names(d)), collapse = ", ")
  )
}

base_by_subj <- d %>%
  select(subj, all_of(base_cols)) %>%
  group_by(subj) %>%
  summarise(
    across(
      all_of(base_cols),
      ~ {
        tmp <- .
        idx <- which(!is.na(tmp))[1]
        if (is.na(idx)) NA_real_ else as.numeric(tmp[idx])
      }
    ),
    .groups = "drop"
  ) %>%
  filter(subj %in% ids) %>%
  arrange(subj)

stopifnot(nrow(base_by_subj) == I)

Xb <- as.matrix(base_by_subj[, base_cols, drop = FALSE])
Xb <- apply(Xb, 2, z_)
n_na_before <- sum(!is.finite(Xb))
if (n_na_before > 0) {
  Xb[!is.finite(Xb)] <- 0
  message(sprintf(
    "Imputati a 0 (post z-score) %d valori baseline mancanti.",
    n_na_before
  ))
}
colnames(Xb) <- c("z_naff_b", "z_det_b", "z_ant_b", "z_dis_b", "z_psy_b")

# ------------------------------------------------------------
# 8) EMA baseline: (a) medie residualizzate su X_base  (b) dinamiche (SD, MSSD, AR1)
# ------------------------------------------------------------
# 8a) Medie EMA al baseline
ema_cols <- list(
  det = find_first(
    d,
    c("pid5_detachment", "ema_detachment", "detachment", "z_det", "det")
  ),
  ant = find_first(
    d,
    c("pid5_antagonism", "ema_antagonism", "antagonism", "z_ant", "ant")
  ),
  dis = find_first(
    d,
    c(
      "pid5_disinhibition",
      "ema_disinhibition",
      "disinhibition",
      "z_dis",
      "dis"
    )
  ),
  psy = find_first(
    d,
    c("pid5_psychoticism", "ema_psychoticism", "psychoticism", "z_psy", "psy")
  )
)

EMA_base <- d %>%
  filter(per == 1L, subj %in% ids) %>%
  transmute(
    subj,
    det = if (!is.null(ema_cols$det)) as.numeric(.data[[ema_cols$det]]) else
      NA_real_,
    ant = if (!is.null(ema_cols$ant)) as.numeric(.data[[ema_cols$ant]]) else
      NA_real_,
    dis = if (!is.null(ema_cols$dis)) as.numeric(.data[[ema_cols$dis]]) else
      NA_real_,
    psy = if (!is.null(ema_cols$psy)) as.numeric(.data[[ema_cols$psy]]) else
      NA_real_
  ) %>%
  group_by(subj) %>%
  summarise(
    det = mean(det, trim = TRIM_PROP, na.rm = TRUE),
    ant = mean(ant, trim = TRIM_PROP, na.rm = TRUE),
    dis = mean(dis, trim = TRIM_PROP, na.rm = TRUE),
    psy = mean(psy, trim = TRIM_PROP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(subj)

EMA_base <- right_join(tibble(subj = ids), EMA_base, by = "subj") %>%
  arrange(subj)
Xema_mean <- EMA_base %>% select(det, ant, dis, psy) %>% as.matrix()
Xema_mean <- apply(Xema_mean, 2, z_)
Xema_mean[!is.finite(Xema_mean)] <- 0
colnames(Xema_mean) <- c("z_det_e", "z_ant_e", "z_dis_e", "z_psy_e")

# 8a-bis) Residualizza le medie EMA su X_base (incremental validity pura)
residualize <- function(y, X) {
  X_ <- cbind(1, X)
  beta <- solve(t(X_) %*% X_, t(X_) %*% y)
  as.numeric(y - X_ %*% beta)
}
Xema_res <- Xema_mean
for (j in seq_len(ncol(Xema_res)))
  Xema_res[, j] <- residualize(Xema_mean[, j], Xb)
# z-score dei residui (facoltativo ma utile)
Xema_res <- apply(Xema_res, 2, z_)
Xema_res[!is.finite(Xema_res)] <- 0
colnames(Xema_res) <- paste0(colnames(Xema_mean), "_res")

# 8b) Dinamiche EMA al baseline: SD, MSSD, AR(1) per ciascun dominio
feats_per_subj <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(sd = NA, mssd = NA, ar1 = NA))
  sdv <- sd(x)
  mssd <- mean(diff(x)^2)
  ar1 <- if (length(x) >= 3)
    suppressWarnings(acf(x, lag.max = 1, plot = FALSE)$acf[2]) else NA_real_
  c(sd = sdv, mssd = mssd, ar1 = ar1)
}

EMA_base_long <- d %>%
  filter(per == 1L, subj %in% ids) %>%
  transmute(
    subj,
    det = if (!is.null(ema_cols$det)) as.numeric(.data[[ema_cols$det]]) else
      NA_real_,
    ant = if (!is.null(ema_cols$ant)) as.numeric(.data[[ema_cols$ant]]) else
      NA_real_,
    dis = if (!is.null(ema_cols$dis)) as.numeric(.data[[ema_cols$dis]]) else
      NA_real_,
    psy = if (!is.null(ema_cols$psy)) as.numeric(.data[[ema_cols$psy]]) else
      NA_real_
  )

make_dyn <- function(var) {
  EMA_base_long %>%
    group_by(subj) %>%
    summarise(t(feats_per_subj(.data[[var]])), .groups = "drop") %>%
    rename_with(~ paste0(var, "_", .x), -subj)
}
dyn_det <- make_dyn("det")
dyn_ant <- make_dyn("ant")
dyn_dis <- make_dyn("dis")
dyn_psy <- make_dyn("psy")

EMA_dyn <- reduce(
  list(dyn_det, dyn_ant, dyn_dis, dyn_psy),
  full_join,
  by = "subj"
) %>%
  right_join(tibble(subj = ids), by = "subj") %>%
  arrange(subj) %>%
  select(-subj) %>%
  as.matrix()

EMA_dyn <- apply(EMA_dyn, 2, function(col) {
  m <- mean(col, na.rm = TRUE)
  s <- sd(col, na.rm = TRUE)
  if (!is.finite(s) || s == 0) col - m else (col - m) / s
})
EMA_dyn[!is.finite(EMA_dyn)] <- 0

# 8c) Scegli le feature EMA finali
X_ema_list <- list()
if (USE_EMA_MEAN_RES == 1L) X_ema_list <- c(X_ema_list, list(Xema_res))
if (USE_EMA_DYN == 1L) X_ema_list <- c(X_ema_list, list(EMA_dyn))

if (length(X_ema_list) == 0) {
  Xema_final <- matrix(0, nrow = I, ncol = 1) # placeholder 1 colonna
  colnames(Xema_final) <- "EMA_placeholder"
} else {
  Xema_final <- do.call(cbind, X_ema_list)
}

# ------------------------------------------------------------
# 9) female allineato
# ------------------------------------------------------------
subj_levels <- sort(unique(d$subj))
female_map <- tibble(
  subj = subj_levels,
  female = as.integer(fem_by_subj >= 0.5)
)
female_vec <- female_map$female[match(ids, female_map$subj)]
female_vec[is.na(female_vec)] <- 1L

# ------------------------------------------------------------
# 10) Standardizza y (migliora geometria posterior)
# ------------------------------------------------------------
y_mean <- mean(y)
y_sd <- sd(y)
y_z <- as.numeric((y - y_mean) / y_sd)
y <- y_z

# ------------------------------------------------------------
# 11) Stan data (usa nome diverso da 'sd' per non confliggere con stats::sd)
# ------------------------------------------------------------
stan_data <- list(
  I = I,
  D_base = 5L,
  D_ema = ncol(Xema_final),
  y = as.vector(y),
  X_base = unname(Xb),
  X_ema = unname(Xema_final),
  female = as.integer(female_vec),
  use_gender = as.integer(USE_GENDER),
  use_ema = 0L # 0 per A, 1 per B
)

# ------------------------------------------------------------
# 12) Fit Stan: Modello A (use_ema=0) e B (use_ema=1)
#      - usa adapt_delta alto per evitare divergenze residue
# ------------------------------------------------------------
stopifnot(file.exists(stan_file_path))
mod <- cmdstan_model(stan_file_path)

fit_A <- mod$sample(
  data = within(stan_data, {
    use_ema <- 0L
  }),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 1500,
  adapt_delta = 0.99,
  max_treedepth = 12,
  seed = SEED
)

fit_B <- mod$sample(
  data = within(stan_data, {
    use_ema <- 1L
  }),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 1500,
  adapt_delta = 0.99,
  max_treedepth = 12,
  seed = SEED + 1
)

# ------------------------------------------------------------
# 13) PSIS-LOO (soggetto-level) e tabella finale
# ------------------------------------------------------------
ll_A <- fit_A$draws("log_lik_subj", format = "matrix")
ll_B <- fit_B$draws("log_lik_subj", format = "matrix")

loo_A <- loo::loo(ll_A, r_eff = NA)
loo_B <- loo::loo(ll_B, r_eff = NA)
cmp <- loo::loo_compare(list(B = loo_B, A = loo_A))
print(cmp)

delta_elpd <- as.numeric(cmp["A", "elpd_diff"]) * -1
se_diff <- as.numeric(cmp["A", "se_diff"])
decisivo <- ifelse(abs(delta_elpd) > 2 * se_diff, "SÌ", "NO")

tbl <- tibble::tibble(
  Model = c("A: baseline-only", "B: baseline + EMA", "ΔELPD (B−A)"),
  elpd = c(
    loo_A$estimates["elpd_loo", "Estimate"],
    loo_B$estimates["elpd_loo", "Estimate"],
    delta_elpd
  ),
  se = c(
    loo_A$estimates["elpd_loo", "SE"],
    loo_B$estimates["elpd_loo", "SE"],
    se_diff
  )
) %>%
  mutate(looic = -2 * elpd)

print(tbl, n = Inf)
cat(sprintf(
  "\nΔELPD (B−A): %.2f ± %.2f  | decisivo se |Δ|>2SE → %s\n",
  delta_elpd,
  se_diff,
  decisivo
))

# ------------------------------------------------------------
# 14) (Opzionale) grafico paired dell’errore assoluto per soggetto
# ------------------------------------------------------------
muA <- {
  a <- posterior::summarise_draws(fit_A$draws("a"))$mean
  bb <- apply(
    as_draws_matrix(fit_A$draws("b_base", format = "matrix")),
    2,
    mean
  )
  mu <- as.vector(a + stan_data$X_base %*% bb)
  if (stan_data$use_gender == 1) {
    bf <- posterior::summarise_draws(fit_A$draws("b_female"))$mean
    mu <- mu + stan_data$female * bf
  }
  mu
}
muB <- {
  a <- posterior::summarise_draws(fit_B$draws("a"))$mean
  bb <- apply(
    as_draws_matrix(fit_B$draws("b_base", format = "matrix")),
    2,
    mean
  )
  cc <- apply(as_draws_matrix(fit_B$draws("c_ema", format = "matrix")), 2, mean)
  mu <- as.vector(a + stan_data$X_base %*% bb + stan_data$X_ema %*% cc)
  if (stan_data$use_gender == 1) {
    bf <- posterior::summarise_draws(fit_B$draws("b_female"))$mean
    mu <- mu + stan_data$female * bf
  }
  mu
}

df_plot <- tibble::tibble(
  subj = ids,
  abs_err_A = abs(y - muA),
  abs_err_B = abs(y - muB)
) %>%
  tidyr::pivot_longer(
    c(abs_err_A, abs_err_B),
    names_to = "model",
    values_to = "abs_err"
  ) %>%
  mutate(model = ifelse(model == "abs_err_A", "A (baseline)", "B (+EMA)"))

gg <- ggplot(df_plot, aes(x = model, y = abs_err, group = subj)) +
  geom_line(alpha = 0.25) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  labs(
    x = NULL,
    y = "Absolute prediction error (subject-level)",
    title = "Per-subject absolute error: Model A vs Model B"
  ) +
  theme_minimal(base_size = 12)

print(gg)

# ------------------------------------------------------------
# Fine script
# ------------------------------------------------------------
