suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(loo)
  library(stringr)
  library(rio)
})

# === 0) Utilità ==============================================================
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
  # fallback 0..100 → 1..7
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  y <- 1 + 6 * (x - xmin) / (xmax - xmin)
  as.integer(round(pmax(1, pmin(7, y))))
}
first_non_na <- function(v) {
  i <- which(!is.na(v))[1]
  if (length(i)) v[i] else NA
}

# === 1) Carica dati e colonne chiave (adatta path) ===========================

# Usa direttamente il tuo CSV "ema_plus_scales_cleaned.csv"
d <- rio::import(here::here("data", "processed", "ema_plus_scales_cleaned.csv"))

# Normalizza i nomi dei domini PID-5 baseline
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

# Identificativi e periodo
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

# Items per NA score
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

# Sesso (0/1)
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

# Baseline PID-5
base_map <- c(
  naff = "pid5_negative_affect_baseline",
  det = "pid5_detachment_baseline",
  ant = "pid5_antagonism_baseline",
  dis = "pid5_disinhibition_baseline",
  psy = "pid5_psychoticism_baseline"
)
stopifnot(all(base_map %in% names(d)))

base_by_subj <- d %>%
  group_by(subj) %>%
  summarise(
    z_naff_b = z_(first_non_na(.data[[base_map["naff"]]])),
    z_det_b = z_(first_non_na(.data[[base_map["det"]]])),
    z_ant_b = z_(first_non_na(.data[[base_map["ant"]]])),
    z_dis_b = z_(first_non_na(.data[[base_map["dis"]]])),
    z_psy_b = z_(first_non_na(.data[[base_map["psy"]]])),
    .groups = "drop"
  )

# EMA PID-5 (entro-soggetto → tra-soggetto: media PRE)
ema_candidates <- list(
  det = names(d)[
    tolower(names(d)) %in%
      c("ema_detachment", "pid5_detachment", "detachment", "z_det", "det")
  ][1],
  ant = names(d)[
    tolower(names(d)) %in%
      c("ema_antagonism", "pid5_antagonism", "antagonism", "z_ant", "ant")
  ][1],
  dis = names(d)[
    tolower(names(d)) %in%
      c(
        "ema_disinhibition",
        "pid5_disinhibition",
        "disinhibition",
        "z_dis",
        "dis"
      )
  ][1],
  psy = names(d)[
    tolower(names(d)) %in%
      c("ema_psychoticism", "pid5_psychoticism", "psychoticism", "z_psy", "psy")
  ][1]
)
ema_names <- unlist(ema_candidates[!vapply(ema_candidates, is.na, TRUE)])

# Filtri minimi PRE/POST per soggetto (almeno 2 osservazioni)
subj_counts <- d %>%
  count(subj, period) %>%
  tidyr::pivot_wider(names_from = period, values_from = n, values_fill = 0)
ok_subj <- subj_counts %>%
  mutate(ok = pre >= 2 & post >= 2) %>%
  filter(ok) %>%
  pull(subj)

d2 <- d %>% filter(subj %in% ok_subj)

# Outcome per soggetto: media POST (z tra soggetti)
y_by_subj <- d2 %>%
  filter(period == "post") %>%
  group_by(subj) %>%
  summarise(y_post = mean(y, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_post = z_(y_post))

# (Opzionale) covariata y_pre per controllo del livello
y_pre_by_subj <- d2 %>%
  filter(period == "pre") %>%
  group_by(subj) %>%
  summarise(y_pre = mean(y, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_pre = z_(y_pre))

# EMA PRE → medie per dominio + winsor + z tra soggetti
ema_pre_by_subj <- NULL
if (length(ema_names) > 0) {
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
  # rinomina in w_*
  names(ema_pre_by_subj)[names(ema_pre_by_subj) %in% ema_names] <- paste0(
    "w_",
    names(ema_pre_by_subj)[names(ema_pre_by_subj) %in% ema_names]
  )
  # winsor + z
  for (nm in setdiff(names(ema_pre_by_subj), "subj")) {
    ema_pre_by_subj[[nm]] <- z_(winsor(ema_pre_by_subj[[nm]]))
  }
}

# Merge soggetto-livello
df_s <- y_by_subj %>%
  left_join(y_pre_by_subj, by = "subj") %>%
  left_join(base_by_subj, by = "subj") %>%
  left_join(female_by_subj, by = "subj") %>%
  {
    if (!is.null(ema_pre_by_subj))
      left_join(., ema_pre_by_subj, by = "subj") else .
  } %>%
  drop_na(y_post, z_naff_b, z_det_b, z_ant_b, z_dis_b, z_psy_b, female)

# === 2) Modelli semplici =====================================================

# Ridotto: y_post ~ baseline + female (+ y_pre opzionale)
form_reduced <- as.formula(
  "y_post ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female + y_pre"
)
# Completo: aggiunge le EMA PRE
ema_cols <- grep("^w_", names(df_s), value = TRUE)
form_full <- as.formula(
  paste(
    "y_post ~ z_naff_b + z_det_b + z_ant_b + z_dis_b + z_psy_b + female + y_pre",
    if (length(ema_cols)) paste("+", paste(ema_cols, collapse = " + ")) else ""
  )
)

m_red <- lm(form_reduced, data = df_s)
m_full <- lm(form_full, data = df_s)

# === 3) PSIS-LOO su lm() ====================================================

log_lik_lm <- function(fit) {
  y <- fit$model[[1]]
  mu <- fitted(fit)
  sig <- sqrt(sum(residuals(fit)^2) / length(y))
  matrix(dnorm(y, mean = mu, sd = sig, log = TRUE), ncol = 1)
}

ll_red <- log_lik_lm(m_red)
ll_full <- log_lik_lm(m_full)

loo_red <- loo::loo(ll_red)
loo_full <- loo::loo(ll_full)

print(loo::loo_compare(list(Full = loo_full, Reduced = loo_red)))
# Tipicamente: elpd_diff > 0 a favore del Full

# Output sintetico
cat("\n--- Summary -------------------------------------------------\n")
cat(sprintf("n subjects = %d\n", nrow(df_s)))
cat("Reduced:", deparse(form_reduced), "\n")
cat("Full   :", deparse(form_full), "\n\n")
print(loo::pareto_k_table(loo_red))
print(loo::pareto_k_table(loo_full))
