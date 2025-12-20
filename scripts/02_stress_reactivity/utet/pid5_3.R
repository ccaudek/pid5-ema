suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(rio)
  library(here)
  library(brms)
  library(loo)
})

# ---------- utility robuste ----------
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

# ---------- parametri minimi ----------
MIN_BEEPS_PER_DAY <- 1 # beeps minimi nel giorno PRE e nel giorno POST

# ---------- 1) carica dati grezzi ----------
d <- rio::import(here::here("data", "processed", "ema_plus_scales_cleaned.csv"))

# id soggetto
id_col <- names(d)[
  names(d) %in% c("user_id", "id", "subject_id", "participant_id")
][1]
stopifnot(length(id_col) == 1)
d$subj <- as.integer(factor(d[[id_col]]))

# periodo (baseline/pre/post)
per_col <- names(d)[names(d) %in% c("exam_period", "period", "phase")][1]
stopifnot(length(per_col) == 1)
per_raw <- tolower(as.character(d[[per_col]]))
d$period <- dplyr::case_when(
  per_raw %in% c("baseline", "base", "t0", "bl") ~ "baseline",
  per_raw %in% c("pre", "pre_exam", "pre-exam", "preexam") ~ "pre",
  per_raw %in% c("post", "post_exam", "post-exam", "postexam") ~ "post",
  TRUE ~ "baseline"
)

# ---------- 2) negative affect y (scala 1..7 → composite → z) ----------
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

# ---------- 3) sesso (female 0/1 a livello soggetto) ----------
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
female_by_subj$female[!is.finite(female_by_subj$female)] <- 1L # fallback

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

# ---------- 5) selezione “giorno prima” (ultimo PRE) e “giorno dopo” (primo POST) ----------
# Assumo che 'day' sia un IDate (come nel tuo glimpse). Se è <date> va bene uguale.
stopifnot("day" %in% names(d))

# per ogni soggetto: trova ultimo giorno PRE e primo giorno POST
pre_days <- d %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(day_pre_last = max(day, na.rm = TRUE), .groups = "drop")
post_days <- d %>%
  filter(period == "post") %>%
  group_by(subj) %>%
  summarise(day_post_first = min(day, na.rm = TRUE), .groups = "drop")

# beeps per giorno selezionato
pre_day_beeps <- d %>%
  inner_join(pre_days, by = "subj") %>%
  filter(period == "pre", day == day_pre_last)
post_day_beeps <- d %>%
  inner_join(post_days, by = "subj") %>%
  filter(period == "post", day == day_post_first)

# richiedi >=2 beeps in ciascun giorno
pre_ok <- pre_day_beeps %>%
  count(subj) %>%
  filter(n >= MIN_BEEPS_PER_DAY) %>%
  pull(subj)
post_ok <- post_day_beeps %>%
  count(subj) %>%
  filter(n >= MIN_BEEPS_PER_DAY) %>%
  pull(subj)
ok_subj <- intersect(pre_ok, post_ok)

pre_day_y <- pre_day_beeps %>%
  filter(subj %in% ok_subj) %>%
  group_by(subj) %>%
  summarise(y_pre_day = rob_mean(y), .groups = "drop")
post_day_y <- post_day_beeps %>%
  filter(subj %in% ok_subj) %>%
  group_by(subj) %>%
  summarise(y_post_day = rob_mean(y), .groups = "drop")

# outcome: reattività = (giorno prima) − (giorno dopo)
df_delta <- pre_day_y %>%
  inner_join(post_day_y, by = "subj") %>%
  mutate(delta = y_pre_day - y_post_day)

# ---------- 6) dataset finale soggetto-livello ----------
df_s <- df_delta %>%
  left_join(base_by_subj, by = "subj") %>%
  left_join(female_by_subj, by = "subj") %>%
  mutate(female = as.integer(female)) %>%
  tidyr::drop_na(delta, z_naff_b, z_det_b, z_ant_b, z_dis_b, z_psy_b, female)

cat(
  "Soggetti con PRE-last & POST-first validi:",
  dplyr::n_distinct(df_s$subj),
  "\n"
)

# ---------- 7) brms: reduced vs full ----------
priors <- c(
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("student_t(3, 0, 1)", class = "sigma"),
  set_prior("gamma(2, 0.1)", class = "nu")
)

form_red <- bf(delta ~ female) # solo genere
form_full <- bf(
  delta ~ female + z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b
) # genere + 5 tratti

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
  seed = 7201,
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
  seed = 7202,
  save_pars = save_pars(all = TRUE)
)

fit_red <- add_criterion(fit_red, "loo", moment_match = TRUE)
fit_full <- add_criterion(fit_full, "loo", moment_match = TRUE)

print(loo::loo_compare(fit_full, fit_red))

# report breve
cmp <- loo::loo_compare(fit_full, fit_red)
cat(sprintf(
  "\nΔELPD (Full - Reduced) = %.1f (SE = %.1f)\n",
  as.numeric(cmp[1, "elpd_diff"]),
  as.numeric(cmp[1, "se_diff"])
))

# diagnostiche (opzionali per slide)
print(fixef(fit_red))
print(fixef(fit_full))
pp_check(fit_red)
pp_check(fit_full)
