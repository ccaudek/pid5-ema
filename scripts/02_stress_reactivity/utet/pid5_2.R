## ------------------------------------------------------------------
## DEMO DIDATTICA — Pipeline completa (filtri, feature, brms, LOO)
## ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(brms)
  library(loo)
  library(rio)
  library(here)
})

## -------------------- 0) Utility robuste --------------------------
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) x - m else (x - m) / s
}
winsor <- function(x, p = c(.05, .95)) {
  qs <- stats::quantile(x, probs = p, na.rm = TRUE)
  pmin(pmax(x, qs[[1]]), qs[[2]])
}
rob_mean <- function(x, trim = .20) mean(x, trim = trim, na.rm = TRUE)
first_non_na <- function(v) {
  i <- which(!is.na(v))[1]
  if (length(i)) v[i] else NA
}
as_item_1_7 <- function(x) {
  x <- suppressWarnings(readr::parse_number(as.character(x)))
  if (all(is.na(x))) return(as.integer(x))
  if (all(x[is.finite(x)] %in% 1:7)) return(as.integer(x))
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  # mappa range generico (es. 0–100) su 1..7
  y <- 1 + 6 * (x - xmin) / (xmax - xmin)
  as.integer(round(pmax(1, pmin(7, y))))
}

## -------------------- 1) Parametri di filtro ----------------------
MIN_PRE <- 2 # beeps min nel PRE (per drift robusto)
MIN_POST <- 1 # beeps min nel POST (per Δy)
EARLY_FRAC <- 0.50 # frazione iniziale per "early" (es. 50% della PRE)
LATE_FRAC <- 1 / 3 # frazione finale per "late"  (es. ultimo terzo della PRE)

## -------------------- 2) Caricamento & normalizzazioni ------------
d <- rio::import(here::here("data", "processed", "ema_plus_scales_cleaned.csv"))

# ID soggetto
id_col <- names(d)[
  names(d) %in% c("user_id", "id", "subject_id", "participant_id")
][1]
stopifnot(length(id_col) == 1)
d$subj <- as.integer(factor(d[[id_col]]))

length(unique(d$subj))

# Periodo (baseline / pre / post)
per_col <- names(d)[names(d) %in% c("exam_period", "period", "phase")][1]
stopifnot(length(per_col) == 1)
per_raw <- tolower(as.character(d[[per_col]]))
d$period <- dplyr::case_when(
  per_raw %in% c("baseline", "base", "t0", "bl") ~ "baseline",
  per_raw %in% c("pre", "pre_exam", "pre-exam", "preexam") ~ "pre",
  per_raw %in% c("post", "post_exam", "post-exam", "postexam") ~ "post",
  TRUE ~ "baseline"
)

# Elementi affettivi → 1..7 e outcome y (Negative Affect; z a livello beep)
happy_nm <- names(d)[
  tolower(names(d)) %in% c("happy", "happiness", "felice", "contento")
][1]
sad_nm <- names(d)[tolower(names(d)) %in% c("sad", "sadness", "triste")][1]
sat_nm <- names(d)[
  tolower(names(d)) %in%
    c("satisfied", "satisfaction", "soddisfatto", "soddisfazione")
][1]
ang_nm <- names(d)[
  tolower(names(d)) %in% c("angry", "anger", "arrabbiato", "rabbia")
][1]
stopifnot(all(!is.na(c(happy_nm, sad_nm, sat_nm, ang_nm))))

d$happy <- as_item_1_7(d[[happy_nm]])
d$sad <- as_item_1_7(d[[sad_nm]])
d$satis <- as_item_1_7(d[[sat_nm]])
d$angry <- as_item_1_7(d[[ang_nm]])

d <- d %>%
  filter(is.finite(happy), is.finite(sad), is.finite(satis), is.finite(angry))

NA_score <- (7 - d$happy) + (7 - d$satis) + d$sad + d$angry
d$y <- z_(NA_score)

# Sesso a livello soggetto (0/1, female)
sex_col <- names(d)[
  tolower(names(d)) %in% c("female", "sex", "gender", "sesso")
][1]
if (length(sex_col) == 1) {
  sx <- tolower(as.character(d[[sex_col]]))
  d$female_row <- as.integer(
    sx %in% c("f", "female", "donna", "femmina", "woman", "1", "true")
  )
} else {
  d$female_row <- NA_integer_
}

female_by_subj <- d %>%
  group_by(subj) %>%
  summarise(
    female = as.integer(mean(female_row, na.rm = TRUE) >= .5),
    .groups = "drop"
  )
female_by_subj$female[!is.finite(female_by_subj$female)] <- 1L # fallback

# Baseline PID-5: normalizza i nomi e crea z-score tra soggetti
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

base_by_subj_raw <- d %>%
  group_by(subj) %>%
  summarise(
    naff_b = first_non_na(pid5_negative_affect_baseline),
    det_b = first_non_na(pid5_detachment_baseline),
    ant_b = first_non_na(pid5_antagonism_baseline),
    dis_b = first_non_na(pid5_disinhibition_baseline),
    psy_b = first_non_na(pid5_psychoticism_baseline),
    .groups = "drop"
  )

base_by_subj <- base_by_subj_raw %>%
  mutate(
    z_naff_b = z_(naff_b),
    z_det_b = z_(det_b),
    z_ant_b = z_(ant_b),
    z_dis_b = z_(dis_b),
    z_psy_b = z_(psy_b)
  ) %>%
  select(subj, z_naff_b, z_det_b, z_ant_b, z_dis_b, z_psy_b)

## -------------------- 3) Filtri qualità “raw → analisi” --------------
# vincolo minimo per robustezza del drift e del Δy
qc <- d %>%
  count(subj, period) %>%
  tidyr::pivot_wider(names_from = period, values_from = n, values_fill = 0)

ok_subj <- qc %>%
  mutate(ok = pre >= MIN_PRE & post >= MIN_POST) %>%
  filter(ok) %>%
  pull(subj)

d2 <- d %>% filter(subj %in% ok_subj)

## -------------------- 4) Outcome a livello soggetto: Δy --------------
y_pre <- d2 %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(y_pre = rob_mean(y), .groups = "drop")
y_post <- d2 %>%
  filter(period == "post") %>%
  group_by(subj) %>%
  summarise(y_post = rob_mean(y), .groups = "drop")

df_delta <- y_pre %>%
  inner_join(y_post, by = "subj") %>%
  mutate(
    y_pre_z = as.numeric(scale(y_pre)),
    y_post_z = as.numeric(scale(y_post)),
    dy = y_post_z - y_pre_z
  )

## -------------------- 5) Feature EMA: “drift PRE” (late − early) -----
# ordina nel PRE e taglia in early/late in base a EARLY_FRAC / LATE_FRAC
# NB: se non hai un time stamp denso, l'ordinamento day+hour è sufficiente
d_pre <- d2 %>%
  filter(period == "pre") %>%
  arrange(subj, day, hour) %>%
  group_by(subj) %>%
  mutate(r = row_number(), n = n()) %>%
  ungroup() %>%
  group_by(subj) %>%
  mutate(
    cut_early = max(1L, floor(n[1] * EARLY_FRAC)),
    cut_late = min(n[1], n[1] - ceiling(n[1] * LATE_FRAC) + 1L),
    part = case_when(
      r <= cut_early ~ "early",
      r >= cut_late ~ "late",
      TRUE ~ "mid"
    )
  ) %>%
  ungroup()

# EMA usata per il drift (qui detachment; cambia se preferisci)
ema_var <- "pid5_detachment"
stopifnot(ema_var %in% names(d_pre))

det_drift <- d_pre %>%
  group_by(subj) %>%
  summarise(
    w_det_early = rob_mean(as.numeric(.data[[ema_var]][part == "early"])),
    w_det_late = rob_mean(as.numeric(.data[[ema_var]][part == "late"])),
    w_det_drift_raw = w_det_late - w_det_early,
    .groups = "drop"
  ) %>%
  mutate(
    w_det_drift = z_(winsor(w_det_drift_raw, c(.05, .95)))
  )

# (Opzionale) modulatore di contesto: minaccia media PRE
thr_pre <- NULL
if ("context_threat" %in% names(d_pre)) {
  thr_pre <- d_pre %>%
    group_by(subj) %>%
    summarise(
      threat_pre_mean = z_(winsor(
        rob_mean(as.numeric(context_threat)),
        c(.05, .95)
      )),
      .groups = "drop"
    )
}

## -------------------- 6) Dataset soggetto-livello finale -------------
df_s <- df_delta %>%
  left_join(base_by_subj, by = "subj") %>%
  left_join(female_by_subj, by = "subj") %>%
  left_join(det_drift, by = "subj")

if (!is.null(thr_pre)) df_s <- df_s %>% left_join(thr_pre, by = "subj")

df_s <- df_s %>%
  mutate(
    female = as.integer(female),
    w_det_drift = z_(w_det_drift)
  ) %>%
  tidyr::drop_na(
    dy,
    z_naff_b,
    z_det_b,
    z_ant_b,
    z_dis_b,
    z_psy_b,
    female,
    w_det_drift
  )

cat("Soggetti finali:", dplyr::n_distinct(df_s$subj), "\n")

## -------------------- 7) BRMS: ridotto vs completo -------------------
priors_common <- c(
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"), # baseline + female
  set_prior("student_t(3, 0, 1)", class = "sigma"),
  set_prior("gamma(2, 0.1)", class = "nu")
)
prior_drift <- set_prior("normal(0, 0.35)", class = "b", coef = "w_det_drift")

form_red <- bf(dy ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female)
form_full1 <- bf(
  dy ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female + w_det_drift
)

# (Opzionale) interazione con threat PRE
use_interaction <- "threat_pre_mean" %in% names(df_s)
if (use_interaction) {
  prior_inter <- set_prior(
    "normal(0, 0.35)",
    class = "b",
    coef = "w_det_drift:threat_pre_mean"
  )
  form_full2 <- bf(
    dy ~
      z_naff_b +
        z_det_b +
        z_ant_b +
        z_dis_b +
        z_psy_b +
        female +
        w_det_drift * threat_pre_mean
  )
}

fit_red <- brm(
  form_red,
  data = df_s,
  family = student(),
  prior = priors_common,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 6101,
  save_pars = save_pars(all = TRUE)
)

fit_full1 <- brm(
  form_full1,
  data = df_s,
  family = student(),
  prior = c(priors_common, prior_drift),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 6102,
  save_pars = save_pars(all = TRUE)
)

if (use_interaction) {
  fit_full2 <- brm(
    form_full2,
    data = df_s,
    family = student(),
    prior = c(priors_common, prior_drift, prior_inter),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    control = list(adapt_delta = .97),
    backend = "cmdstanr",
    seed = 6103,
    save_pars = save_pars(all = TRUE)
  )
}

## -------------------- 8) Valutazione predittiva ---------------------
fit_red <- add_criterion(fit_red, "loo", moment_match = TRUE)
fit_full1 <- add_criterion(fit_full1, "loo", moment_match = TRUE)
print(loo::loo_compare(fit_full1, fit_red))

if (use_interaction) {
  fit_full2 <- add_criterion(fit_full2, "loo", moment_match = TRUE)
  print(loo::loo_compare(fit_full2, fit_red))
  print(loo::loo_compare(fit_full2, fit_full1))
}

## (Facoltativo) K-fold a livello soggetto per SE più stabili:
# kf_red   <- kfold(fit_red,   K=5, chains=2, cores=4)
# kf_full1 <- kfold(fit_full1, K=5, chains=2, cores=4)
# loo::loo_compare(kf_full1, kf_red)

## -------------------- 9) Diagnostiche rapide ------------------------
print(fixef(fit_red))
print(fixef(fit_full1))
pp_check(fit_red)
pp_check(fit_full1)
