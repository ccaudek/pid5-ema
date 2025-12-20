suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(brms)
  library(loo)
  library(rio)
  library(here)
})

set.seed(1234)

# --- utilità ---------------------------------------------------------------
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) x - m else (x - m) / s
}
winsor <- function(x, p = c(.01, .99)) {
  qs <- quantile(x, probs = p, na.rm = TRUE)
  pmin(pmax(x, qs[[1]]), qs[[2]])
}
as_item_1_7 <- function(x) {
  x <- suppressWarnings(readr::parse_number(as.character(x)))
  if (all(is.na(x))) return(as.integer(x))
  if (all(x[is.finite(x)] %in% 1:7)) return(as.integer(x))
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  y <- 1 + 6 * (x - xmin) / (xmax - xmin)
  as.integer(round(pmax(1, pmin(7, y))))
}
first_non_na <- function(v) {
  i <- which(!is.na(v))[1]
  if (length(i)) v[i] else NA
}

# --- 1) dati ---------------------------------------------------------------
d <- rio::import(here::here("data", "processed", "ema_plus_scales_cleaned.csv"))

# id e periodo
id_col <- names(d)[
  names(d) %in% c("user_id", "id", "subject_id", "participant_id")
][1]
stopifnot(length(id_col) == 1)
per_col <- names(d)[names(d) %in% c("exam_period", "period", "phase")][1]
stopifnot(length(per_col) == 1)

d$subj <- as.integer(factor(d[[id_col]]))
per_raw <- tolower(as.character(d[[per_col]]))
d$period <- dplyr::case_when(
  per_raw %in% c("baseline", "base", "t0", "bl") ~ "baseline",
  per_raw %in% c("pre", "pre_exam", "pre-exam", "preexam") ~ "pre",
  per_raw %in% c("post", "post_exam", "post-exam", "postexam") ~ "post",
  TRUE ~ "baseline"
)

# item NA
happy <- names(d)[
  tolower(names(d)) %in% c("happy", "happiness", "felice", "contento")
][1]
sad <- names(d)[tolower(names(d)) %in% c("sad", "sadness", "triste")][1]
satis <- names(d)[
  tolower(names(d)) %in%
    c("satisfied", "satisfaction", "soddisfatto", "soddisfazione")
][1]
angry <- names(d)[
  tolower(names(d)) %in% c("angry", "anger", "arrabbiato", "rabbia")
][1]
stopifnot(all(!is.na(c(happy, sad, satis, angry))))

d$happy <- as_item_1_7(d[[happy]])
d$sad <- as_item_1_7(d[[sad]])
d$satis <- as_item_1_7(d[[satis]])
d$angry <- as_item_1_7(d[[angry]])

keep <- with(
  d,
  is.finite(happy) & is.finite(sad) & is.finite(satis) & is.finite(angry)
)
d <- d[keep, , drop = FALSE]

NA_score <- (7 - d$happy) + (7 - d$satis) + d$sad + d$angry
d$y <- z_(NA_score)

# sesso a livello soggetto
sex_col <- names(d)[
  tolower(names(d)) %in% c("female", "sex", "gender", "sesso")
][1]
if (length(sex_col) == 1) {
  sx <- tolower(as.character(d[[sex_col]]))
  d$female_row <- as.integer(
    sx %in% c("f", "female", "donna", "femmina", "woman", "1", "true")
  )
} else {
  d$female_row <- 1L
}
female_by_subj <- d %>%
  group_by(subj) %>%
  summarise(
    female = as.integer(mean(female_row, na.rm = TRUE) >= .5),
    .groups = "drop"
  )

# normalizza nomi baseline PID-5
d <- d %>%
  dplyr::rename(
    pid5_negative_affect_baseline = dplyr::any_of(c(
      "domain_negative_affect_baseline",
      "pid5_negative_affect_baseline"
    )),
    pid5_detachment_baseline = dplyr::any_of(c(
      "domain_detachment_baseline",
      "pid5_detachment_baseline"
    )),
    pid5_antagonism_baseline = dplyr::any_of(c(
      "domain_antagonism_baseline",
      "pid5_antagonism_baseline"
    )),
    pid5_disinhibition_baseline = dplyr::any_of(c(
      "domain_disinhibition_baseline",
      "pid5_disinhibition_baseline"
    )),
    pid5_psychoticism_baseline = dplyr::any_of(c(
      "domain_psychoticism_baseline",
      "pid5_psychoticism_baseline"
    ))
  )

# baseline PID-5 a livello soggetto
base_by_subj <- d %>%
  group_by(subj) %>%
  summarise(
    z_naff_b = z_(first_non_na(pid5_negative_affect_baseline)),
    z_det_b = z_(first_non_na(pid5_detachment_baseline)),
    z_ant_b = z_(first_non_na(pid5_antagonism_baseline)),
    z_dis_b = z_(first_non_na(pid5_disinhibition_baseline)),
    z_psy_b = z_(first_non_na(pid5_psychoticism_baseline)),
    .groups = "drop"
  )

# EMA candidati (variabili time-varying)
ema_names <- c(
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)
ema_names <- ema_names[ema_names %in% names(d)]

# filtro: almeno 2 PRE e 2 POST
subj_counts <- d %>%
  count(subj, period) %>%
  tidyr::pivot_wider(names_from = period, values_from = n, values_fill = 0)
ok_subj <- subj_counts %>%
  mutate(ok = pre >= 2 & post >= 2) %>%
  filter(ok) %>%
  pull(subj)
d2 <- d %>% filter(subj %in% ok_subj)

# outcome y_post (z tra soggetti) e covariata y_pre (z tra soggetti)
y_post <- d2 %>%
  filter(period == "post") %>%
  group_by(subj) %>%
  summarise(y_post = mean(y, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_post = z_(y_post))

y_pre <- d2 %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(y_pre = mean(y, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_pre = z_(y_pre))

# EMA PRE: media per soggetto, winsor + z
ema_pre_by_subj <- NULL
if (length(ema_names)) {
  ema_pre_by_subj <- d2 %>%
    filter(period == "pre") %>%
    group_by(subj) %>%
    summarise(
      across(
        all_of(ema_names),
        ~ mean(suppressWarnings(as.numeric(.x)), na.rm = TRUE)
      ),
      .groups = "drop"
    )
  # rinomina w_*
  names(ema_pre_by_subj)[names(ema_pre_by_subj) %in% ema_names] <-
    paste0("w_", names(ema_pre_by_subj)[names(ema_pre_by_subj) %in% ema_names])
  for (nm in setdiff(names(ema_pre_by_subj), "subj")) {
    ema_pre_by_subj[[nm]] <- z_(winsor(ema_pre_by_subj[[nm]]))
  }
}

# dataset soggetto-livello finale
df_s <- y_post %>%
  left_join(y_pre, by = "subj") %>%
  left_join(base_by_subj, by = "subj") %>%
  left_join(female_by_subj, by = "subj") %>%
  {
    if (!is.null(ema_pre_by_subj))
      left_join(., ema_pre_by_subj, by = "subj") else .
  } %>%
  drop_na(y_post, y_pre, z_naff_b, z_det_b, z_ant_b, z_dis_b, z_psy_b, female)

# --- 2) brms: modelli ridotto e completo -----------------------------------
# Nota: salviamo i parametri per LOO; family Gaussian
priors <- c(
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"), # predittori già z
  set_prior("student_t(3, 0, 1)", class = "sigma")
)

form_reduced <- bf(
  y_post ~ y_pre + z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female
)

ema_cols <- grep("^w_", names(df_s), value = TRUE)
form_full <- if (length(ema_cols)) {
  bf(as.formula(paste(
    "y_post ~ y_pre + z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female +",
    paste(ema_cols, collapse = " + ")
  )))
} else form_reduced

fit_red <- brm(
  formula = form_reduced,
  data = df_s,
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95),
  seed = 1234,
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

fit_full <- brm(
  formula = form_full,
  data = df_s,
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95),
  seed = 1235,
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

# --- 3) PSIS-LOO (con moment matching e reloo se serve) ---------------------
fit_red <- add_criterion(fit_red, "loo", moment_match = TRUE)
fit_full <- add_criterion(fit_full, "loo", moment_match = TRUE)

# Se rimangono k alti, si può forzare il refit locale:
if (any(fit_red$criteria$loo$diagnostics$pareto_k > 0.7)) {
  fit_red <- add_criterion(fit_red, "loo", reloo = TRUE)
}
if (any(fit_full$criteria$loo$diagnostics$pareto_k > 0.7)) {
  fit_full <- add_criterion(fit_full, "loo", reloo = TRUE)
}

cmp <- loo_compare(fit_full, fit_red)
print(cmp)

# Report compatto
elpd_diff <- as.numeric(cmp[1, "elpd_diff"])
se_diff <- as.numeric(cmp[1, "se_diff"])
cat(sprintf(
  "\nΔELPD (Full - Reduced) = %.1f (SE = %.1f)\n",
  elpd_diff,
  se_diff
))


# Modelli senza y_pre (entrambi)
form_reduced_nopre <- bf(
  y_post ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female
)
ema_cols <- grep("^w_", names(df_s), value = TRUE)
form_full_nopre <- bf(as.formula(
  paste(
    "y_post ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female",
    if (length(ema_cols)) paste("+", paste(ema_cols, collapse = " + ")) else ""
  )
))

fit_red2 <- brm(
  form_reduced_nopre,
  data = df_s,
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  control = list(adapt_delta = .95),
  seed = 11,
  save_pars = save_pars(all = TRUE)
)
fit_full2 <- brm(
  form_full_nopre,
  data = df_s,
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  control = list(adapt_delta = .95),
  seed = 12,
  save_pars = save_pars(all = TRUE)
)
fit_red2 <- add_criterion(fit_red2, "loo", moment_match = TRUE)
fit_full2 <- add_criterion(fit_full2, "loo", moment_match = TRUE)
loo_compare(fit_full2, fit_red2)


####### V3 ------------

# --- PRE/POST robusti (soggetto-livello) ---
rob_mean <- function(x, trim = .20) mean(x, trim = trim, na.rm = TRUE)

y_pre_rob <- d2 %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(y_pre = rob_mean(y), .groups = "drop")

y_post_rob <- d2 %>%
  filter(period == "post") %>%
  group_by(subj) %>%
  summarise(y_post = rob_mean(y), .groups = "drop")

df_delta <- y_pre_rob %>%
  inner_join(y_post_rob, by = "subj") %>%
  mutate(
    y_pre_z = as.numeric(scale(y_pre)),
    y_post_z = as.numeric(scale(y_post)),
    dy = y_post_z - y_pre_z # outcome = cambiamento
  ) %>%
  left_join(base_by_subj, by = "subj") %>%
  left_join(female_by_subj, by = "subj")


# EMA PRE: solo detachment, riassunto robusto, winsor + z
ema_det_pre <- d2 %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(w_det = rob_mean(as.numeric(pid5_detachment)), .groups = "drop") %>%
  mutate(w_det = z_(winsor(w_det)))

df_delta <- df_delta %>%
  left_join(ema_det_pre, by = "subj") %>%
  drop_na(dy, z_naff_b, z_det_b, z_ant_b, z_dis_b, z_psy_b, female, w_det)


# Priors comuni (valgono per entrambi i modelli)
priors_common <- c(
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"), # predittori z
  set_prior("student_t(3, 0, 1)", class = "sigma"),
  set_prior("gamma(2, 0.1)", class = "nu") # df della t
)

# Prior addizionale SOLO per il coefficiente dell'EMA w_det nel modello full
prior_wdet <- set_prior("normal(0, 0.35)", class = "b", coef = "w_det")

# Modelli
form_red_d <- bf(dy ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female)
form_full_d <- bf(
  dy ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female + w_det
)

# Fit (backend cmdstanr ok)
fit_red_d <- brm(
  formula = form_red_d,
  data = df_delta,
  family = student(),
  prior = priors_common, # <-- niente prior_wdet qui
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 41,
  save_pars = save_pars(all = TRUE)
)

fit_full_d <- brm(
  formula = form_full_d,
  data = df_delta,
  family = student(),
  prior = c(priors_common, prior_wdet), # <-- aggiungo solo qui
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 42,
  save_pars = save_pars(all = TRUE)
)

# LOO
fit_red_d <- add_criterion(fit_red_d, "loo", moment_match = TRUE)
fit_full_d <- add_criterion(fit_full_d, "loo", moment_match = TRUE)

loo_compare(fit_full_d, fit_red_d)


##### V4 --------------------------

d_pre <- d2 %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(
    w_det_sd = sd(pid5_detachment, na.rm = TRUE),
    w_det_slope = if (n() >= 3)
      coef(lm(pid5_detachment ~ as.numeric(day)))[2] else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    w_det_sd = z_(winsor(w_det_sd, c(.05, .95))),
    w_det_slope = z_(winsor(w_det_slope, c(.05, .95)))
  )
df_delta2 <- df_delta %>% left_join(d_pre, by = "subj")
