suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(rio)
  library(here)
  library(brms)
  library(loo)
})

# ---------- utility ----------
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) x - m else (x - m) / s
}
rob_mean <- function(x, trim = .20) mean(x, trim = trim, na.rm = TRUE)
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

# ---------- parametri ----------
MIN_BEEPS_PER_DAY <- 1
USE_WEIGHTS <- TRUE # pesi = n_beeps nel giorno POST
USE_COURSE_FE <- TRUE # effetti fissi del corso

# ---------- 1) carica dati ----------
d <- rio::import(here::here("data", "processed", "ema_plus_scales_cleaned.csv"))

# id, periodo, corso
id_col <- names(d)[
  names(d) %in% c("user_id", "id", "subject_id", "participant_id")
][1]
stopifnot(length(id_col) == 1)
d$subj <- as.integer(factor(d[[id_col]]))

per_col <- names(d)[names(d) %in% c("exam_period", "period", "phase")][1]
stopifnot(length(per_col) == 1)
per_raw <- tolower(as.character(d[[per_col]]))
d$period <- dplyr::case_when(
  per_raw %in% c("baseline", "base", "t0", "bl") ~ "baseline",
  per_raw %in% c("pre", "pre_exam", "pre-exam", "preexam") ~ "pre",
  per_raw %in% c("post", "post_exam", "post-exam", "postexam") ~ "post",
  TRUE ~ "baseline"
)

if (!"course" %in% names(d)) d$course <- "one_course"
d$course <- as.factor(d$course)

# ---------- 2) Negative Affect per beep ----------
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

# ---------- 3) female (0/1) ----------
sex_col <- names(d)[
  tolower(names(d)) %in% c("female", "sex", "gender", "sesso")
][1]
if (length(sex_col) == 1) {
  sx <- tolower(as.character(d[[sex_col]]))
  d$female_row <- as.integer(
    sx %in% c("f", "female", "donna", "femmina", "woman", "1", "true")
  )
} else d$female_row <- NA_integer_

female_by_subj <- d %>%
  group_by(subj) %>%
  summarise(
    female = as.integer(mean(female_row, na.rm = TRUE) >= .5),
    .groups = "drop"
  )
female_by_subj$female[!is.finite(female_by_subj$female)] <- 1L

# ---------- 4) baseline PID-5 (z tra soggetti) ----------
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

# ---------- 5) PRE ultimo giorno, POST primo giorno; aggregazioni ----------
stopifnot("day" %in% names(d))
pre_days <- d %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(day_pre_last = max(day, na.rm = TRUE), .groups = "drop")
post_days <- d %>%
  filter(period == "post") %>%
  group_by(subj) %>%
  summarise(day_post_first = min(day, na.rm = TRUE), .groups = "drop")

pre_day <- d %>%
  inner_join(pre_days, by = "subj") %>%
  filter(period == "pre", day == day_pre_last)
post_day <- d %>%
  inner_join(post_days, by = "subj") %>%
  filter(period == "post", day == day_post_first)

# beeps minimi per giorno
pre_ok <- pre_day %>%
  count(subj) %>%
  filter(n >= MIN_BEEPS_PER_DAY) %>%
  pull(subj)
post_ok <- post_day %>%
  count(subj) %>%
  filter(n >= MIN_BEEPS_PER_DAY) %>%
  pull(subj)
ok <- intersect(pre_ok, post_ok)

pre_agg <- pre_day %>%
  filter(subj %in% ok) %>%
  group_by(subj) %>%
  summarise(y_pre_day = rob_mean(y), n_pre_beeps = n(), .groups = "drop")
post_agg <- post_day %>%
  filter(subj %in% ok) %>%
  group_by(subj) %>%
  summarise(y_post_day = rob_mean(y), n_post_beeps = n(), .groups = "drop")

# ---------- 6) dataset soggetto-livello per ANCOVA ----------
df_s <- pre_agg %>%
  inner_join(post_agg, by = "subj") %>%
  left_join(base_by_subj, by = "subj") %>%
  left_join(female_by_subj, by = "subj") %>%
  left_join(d %>% select(subj, course) %>% distinct(), by = "subj") %>%
  mutate(female = as.integer(female))

# pesi = numero di beep del giorno POST (più informativo sulla variabile dipendente)
if (USE_WEIGHTS) df_s$w <- df_s$n_post_beeps else df_s$w <- 1

# rimuovi NA nelle colonne chiave
df_s <- df_s %>%
  tidyr::drop_na(
    y_pre_day,
    y_post_day,
    female,
    z_naff_b,
    z_det_b,
    z_ant_b,
    z_dis_b,
    z_psy_b,
    w,
    course
  )

cat("Soggetti validi (ANCOVA):", dplyr::n_distinct(df_s$subj), "\n")


# ----- PC1 dei 5 domini baseline (tra soggetti) -----
Xpid <- df_s %>% select(z_naff_b, z_det_b, z_ant_b, z_dis_b, z_psy_b)
pc <- prcomp(Xpid, center = TRUE, scale. = TRUE)
df_s$pid5_pc1 <- as.numeric(scale(pc$x[, 1])) # general maladaptivity (z)

# ----- Modelli: ridotto vs full-PC1 -----
priors_pc1 <- c(
  set_prior("student_t(3,0,2.5)", class = "Intercept"),
  set_prior("normal(0,0.5)", class = "b"),
  set_prior("student_t(3,0,1)", class = "sigma"),
  set_prior("gamma(2,0.1)", class = "nu")
)

form_red_pc1 <- bf(y_post_day ~ y_pre_day)
form_full_pc1 <- bf(y_post_day ~ y_pre_day + female)

fit_red_pc1 <- brm(
  form_red_pc1,
  data = df_s,
  family = student(),
  prior = priors_pc1,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 9101,
  save_pars = save_pars(all = TRUE)
)
fit_full_pc1 <- brm(
  form_full_pc1,
  data = df_s,
  family = student(),
  prior = c(
    priors_pc1,
    set_prior("normal(0,0.35)", class = "b", coef = "pid5_pc1")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 9102,
  save_pars = save_pars(all = TRUE)
)

fit_red_pc1 <- add_criterion(fit_red_pc1, "loo", moment_match = TRUE)
fit_full_pc1 <- add_criterion(fit_full_pc1, "loo", moment_match = TRUE)
print(loo::loo_compare(fit_full_pc1, fit_red_pc1))


# ---------- 7) brms: ridotto vs full (Student-t, pesi, FE corso opz.) ----------
priors <- c(
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("student_t(3, 0, 1)", class = "sigma"),
  set_prior("gamma(2, 0.1)", class = "nu")
)

# formule
if (USE_COURSE_FE) {
  form_red <- bf(y_post_day ~ y_pre_day + female + course) # ridotto
  form_full <- bf(
    y_post_day ~
      y_pre_day +
        female +
        course + # full
        z_naff_b +
        z_det_b +
        z_ant_b +
        z_dis_b +
        z_psy_b
  )
} else {
  form_red <- bf(y_post_day ~ y_pre_day + female)
  form_full <- bf(
    y_post_day ~
      y_pre_day + female + z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b
  )
}

w <- df_s$w # rende 'w' visibile a brms

fit_red <- brm(
  form_red,
  data = df_s,
  family = student(),
  prior = priors,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 8101,
  save_pars = save_pars(all = TRUE)
)

fit_full <- brm(
  form_full,
  data = df_s,
  family = student(),
  prior = priors,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = .95),
  backend = "cmdstanr",
  seed = 8102,
  save_pars = save_pars(all = TRUE)
)


# confronto PSIS-LOO + breve report
fit_red <- add_criterion(fit_red, "loo", moment_match = TRUE)
fit_full <- add_criterion(fit_full, "loo", moment_match = TRUE)
print(loo::loo_compare(fit_full, fit_red))
cmp <- loo::loo_compare(fit_full, fit_red)
cat(sprintf(
  "\nΔELPD (Full − Reduced) = %.1f (SE = %.1f)\n",
  as.numeric(cmp[1, "elpd_diff"]),
  as.numeric(cmp[1, "se_diff"])
))

# diagnostiche essenziali
print(fixef(fit_red))
print(fixef(fit_full))
pp_check(fit_red)
pp_check(fit_full)
