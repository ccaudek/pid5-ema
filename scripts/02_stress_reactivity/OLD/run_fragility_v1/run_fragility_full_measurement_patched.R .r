suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(tibble)
  library(conflicted)
  library(rio)
  library(here)
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("sd", "stats")

# ------------------------------------------------------------
# Path e opzioni
# ------------------------------------------------------------
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
stan_file_path <- here::here(
  "scripts",
  "02_stress_reactivity",
  "run_fragility_v1",
  "fragility_full_measurement_patched.stan"
)

SEED <- 20250918

# ------------------------------------------------------------
# Utility
# ------------------------------------------------------------
as_item_to_1_7 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    x[x == ""] <- NA_character_
    x <- readr::parse_number(
      x,
      locale = readr::locale(decimal_mark = ".", grouping_mark = ",")
    )
  }
  x <- as.numeric(x)
  if (all(is.na(x))) return(as.integer(x))
  if (all(x[is.finite(x)] %in% 1:7)) return(as.integer(x))
  xmin <- suppressWarnings(min(x, na.rm = TRUE))
  xmax <- suppressWarnings(max(x, na.rm = TRUE))
  if (is.finite(xmin) && is.finite(xmax) && xmin >= 0 && xmax <= 100) {
    brk <- seq(0, 100, length.out = 8)
    return(as.integer(findInterval(
      x,
      brk,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )))
  }
  y <- 1 + 6 * (x - xmin) / (xmax - xmin)
  y <- pmax(1, pmin(7, y))
  as.integer(round(y))
}
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) x - m else (x - m) / s
}
find_first <- function(df, candidates) {
  nm <- names(df)
  hit <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hit) == 0) NULL else hit[1]
}

# ------------------------------------------------------------
# 1) Caricamento dati + normalizzazione nomi baseline
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
# 2) Indici soggetto e periodo
# ------------------------------------------------------------
subj_col <- d %>%
  select(any_of(c("user_id", "id", "subject_id", "participant_id"))) %>%
  names() %>%
  .[1]
stopifnot(!is.na(subj_col))
d$.__subj__ <- as.integer(factor(d[[subj_col]]))

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
d$.__per__ <- as.integer(per)

# ------------------------------------------------------------
# 3) Genere → female (0/1) a livello soggetto
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
fem_mean <- tapply(
  female_by_row,
  d$.__subj__,
  function(v) mean(v, na.rm = TRUE)
)
fem_mean[is.na(fem_mean)] <- 1
I_all <- max(d$.__subj__)
female_vec <- integer(I_all)
female_vec[as.integer(names(fem_mean))] <- as.integer(fem_mean >= 0.5)

# ------------------------------------------------------------
# 4) Item ordinali (1..7) — long
# ------------------------------------------------------------
happy_col <- find_first(d, c("happy", "happiness", "felice", "contento"))
sad_col <- find_first(d, c("sad", "sadness", "triste"))
satisfied_col <- find_first(
  d,
  c("satisfied", "satisfaction", "soddisfatto", "soddisfazione")
)
angry_col <- find_first(d, c("angry", "anger", "arrabbiato", "rabbia"))
stopifnot(
  !is.null(happy_col),
  !is.null(sad_col),
  !is.null(satisfied_col),
  !is.null(angry_col)
)

d$.__happy__ <- as_item_to_1_7(d[[happy_col]])
d$.__sad__ <- as_item_to_1_7(d[[sad_col]])
d$.__satisfied__ <- as_item_to_1_7(d[[satisfied_col]])
d$.__angry__ <- as_item_to_1_7(d[[angry_col]])

keep_items <- with(
  d,
  is.finite(.__happy__) &
    is.finite(.__sad__) &
    is.finite(.__satisfied__) &
    is.finite(.__angry__)
)
d2 <- d[keep_items, , drop = FALSE]
stopifnot(nrow(d2) > 0)

I <- max(d2$.__subj__)
P <- 3L
K <- 4L
D <- 5L

N_obs <- nrow(d2)
obs_id <- seq_len(N_obs)
subject <- as.integer(d2$.__subj__)
period <- as.integer(d2$.__per__)

y_wide <- cbind(d2$.__happy__, d2$.__sad__, d2$.__satisfied__, d2$.__angry__)
N_items <- 4L * N_obs
y_item <- as.integer(c(t(y_wide)))
item_id <- as.integer(rep(1:4, times = N_obs))
obs_id_long <- as.integer(rep(obs_id, each = 4L))

# ------------------------------------------------------------
# 5) EMA per-beep: SOLO baseline (period==1), z-score
# ------------------------------------------------------------
ema_cols <- list(
  det = find_first(
    d2,
    c("pid5_detachment", "ema_detachment", "detachment", "z_det", "det")
  ),
  ant = find_first(
    d2,
    c("pid5_antagonism", "ema_antagonism", "antagonism", "z_ant", "ant")
  ),
  dis = find_first(
    d2,
    c(
      "pid5_disinhibition",
      "ema_disinhibition",
      "disinhibition",
      "z_dis",
      "dis"
    )
  ),
  psy = find_first(
    d2,
    c("pid5_psychoticism", "ema_psychoticism", "psychoticism", "z_psy", "psy")
  )
)

M_ema <- 0L
ema_val <- numeric(0)
ema_dim <- integer(0)
ema_obs <- integer(0)
if (!all(sapply(ema_cols, is.null))) {
  idx_baseline <- which(period == 1L) # << SOLO baseline
  EMA_df_base <- tibble::tibble(
    det = if (!is.null(ema_cols$det)) as.numeric(d2[[ema_cols$det]]) else
      NA_real_,
    ant = if (!is.null(ema_cols$ant)) as.numeric(d2[[ema_cols$ant]]) else
      NA_real_,
    dis = if (!is.null(ema_cols$dis)) as.numeric(d2[[ema_cols$dis]]) else
      NA_real_,
    psy = if (!is.null(ema_cols$psy)) as.numeric(d2[[ema_cols$psy]]) else
      NA_real_
  )[idx_baseline, , drop = FALSE] %>%
    mutate(det = z_(det), ant = z_(ant), dis = z_(dis), psy = z_(psy))

  add_dim <- function(vals, d_index) {
    ok <- which(is.finite(vals))
    if (length(ok) > 0) {
      ema_val <<- c(ema_val, vals[ok])
      ema_dim <<- c(ema_dim, rep(d_index, length(ok))) # 2..5
      ema_obs <<- c(ema_obs, idx_baseline[ok]) # mappa agli indici N_obs originali
    }
  }
  add_dim(EMA_df_base$det, 2L)
  add_dim(EMA_df_base$ant, 3L)
  add_dim(EMA_df_base$dis, 4L)
  add_dim(EMA_df_base$psy, 5L)
  M_ema <- length(ema_val)
}

# ------------------------------------------------------------
# 6) X_base (I x 5) a livello soggetto — z-score; NA -> 0
# ------------------------------------------------------------
base_cols <- c(
  "pid5_negative_affect_baseline",
  "pid5_detachment_baseline",
  "pid5_antagonism_baseline",
  "pid5_disinhibition_baseline",
  "pid5_psychoticism_baseline"
)
stopifnot(all(base_cols %in% names(d)))

base_by_subj <- d %>%
  dplyr::select(.__subj__, dplyr::all_of(base_cols)) %>%
  dplyr::group_by(.__subj__) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(base_cols),
      ~ {
        tmp <- .
        idx <- which(!is.na(tmp))[1]
        if (is.na(idx)) NA_real_ else as.numeric(tmp[idx])
      }
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.__subj__)

Xb <- as.matrix(base_by_subj[, base_cols, drop = FALSE]) %>% apply(2, z_)
if (nrow(Xb) < I) {
  Xb_full <- matrix(NA_real_, nrow = I, ncol = 5)
  Xb_full[base_by_subj$.__subj__, ] <- Xb
  Xb <- apply(Xb_full, 2, z_)
}
Xb[!is.finite(Xb)] <- 0
colnames(Xb) <- c("z_naff_b", "z_det_b", "z_ant_b", "z_dis_b", "z_psy_b")

female <- as.integer(female_vec[seq_len(I)])

# ------------------------------------------------------------
# 7) Lista dati Stan (interazioni spente)
# ------------------------------------------------------------
stan_data <- list(
  I = I,
  N_obs = N_obs,
  K = K,
  P = P,
  D = D,
  N_items = N_items,
  y_item = y_item,
  item_id = item_id,
  obs_id = obs_id_long,
  subject = subject,
  period = period,
  M_ema = M_ema,
  ema_val = as.vector(ema_val),
  ema_dim = as.integer(ema_dim),
  ema_obs = as.integer(ema_obs),
  X_base = unname(Xb),
  female = female,
  N_interact = 0L, # << niente interazioni
  use_interact = 0L, # << switch off
  use_ema = 0.0
)

# ------------------------------------------------------------
# 8) Compila e lancia A/B
# ------------------------------------------------------------
stopifnot(file.exists(stan_file_path))
mod <- cmdstan_model(stan_file_path)

# fit_A <- mod$sample(
#   data = within(stan_data, {
#     use_ema <- 0.0
#   }),
#   seed = SEED,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1500,
#   iter_sampling = 1500,
#   adapt_delta = 0.99,
#   max_treedepth = 12
# )
#
# fit_B <- mod$sample(
#   data = within(stan_data, {
#     use_ema <- 1.0
#   }),
#   seed = SEED + 1,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1500,
#   iter_sampling = 1500,
#   adapt_delta = 0.99,
#   max_treedepth = 12
# )

fit_A <- mod$variational(
  data = within(stan_data, {
    use_ema <- 0.0
  }),
  seed = 20250917,
  algorithm = "meanfield",
  elbo_samples = 100,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 0
)

fit_B <- mod$variational(
  data = within(stan_data, {
    use_ema <- 1L
  }),
  seed = 20250917,
  algorithm = "meanfield",
  elbo_samples = 100,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 0
)

# ------------------------------------------------------------
# 9) LOO: per osservazione (items/beep) e per soggetto (fragilità)
# ------------------------------------------------------------
# Per osservazione
ll_A_obs <- fit_A$draws("log_lik_obs", format = "matrix")
ll_B_obs <- fit_B$draws("log_lik_obs", format = "matrix")
loo_A_o <- loo::loo(ll_A_obs, moment_match = TRUE)
loo_B_o <- loo::loo(ll_B_obs, moment_match = TRUE)
cmp_o <- loo::loo_compare(list(B = loo_B_o, A = loo_A_o))
print(cmp_o)

# Per soggetto (target)
ll_A_s <- fit_A$draws("log_lik_subj", format = "matrix")
ll_B_s <- fit_B$draws("log_lik_subj", format = "matrix")
loo_A_s <- loo::loo(ll_A_s, r_eff = NA)
loo_B_s <- loo::loo(ll_B_s, r_eff = NA)
cmp_s <- loo::loo_compare(list(B = loo_B_s, A = loo_A_s))
print(cmp_s)

# Tabelline
delta_elpd_o <- as.numeric(cmp_o["A", "elpd_diff"]) * -1
se_diff_o <- as.numeric(cmp_o["A", "se_diff"])
delta_elpd_s <- as.numeric(cmp_s["A", "elpd_diff"]) * -1
se_diff_s <- as.numeric(cmp_s["A", "se_diff"])

tbl <- tibble::tibble(
  Level = c(
    "Observation (items)",
    "Subject (fragility)",
    "ΔELPD (B−A) obs",
    "ΔELPD (B−A) subj"
  ),
  elpd = c(
    loo_A_o$estimates["elpd_loo", "Estimate"],
    loo_A_s$estimates["elpd_loo", "Estimate"],
    delta_elpd_o,
    delta_elpd_s
  ),
  se = c(
    loo_A_o$estimates["elpd_loo", "SE"],
    loo_A_s$estimates["elpd_loo", "SE"],
    se_diff_o,
    se_diff_s
  )
) %>%
  mutate(looic = -2 * elpd)
print(tbl, n = Inf)

cat(sprintf(
  "\nΔELPD_obs (B−A): %.2f ± %.2f | decisivo se |Δ|>2SE → %s",
  delta_elpd_o,
  se_diff_o,
  ifelse(abs(delta_elpd_o) > 2 * se_diff_o, "SÌ", "NO")
))
cat(sprintf(
  "\nΔELPD_subj (B−A): %.2f ± %.2f | decisivo se |Δ|>2SE → %s\n",
  delta_elpd_s,
  se_diff_s,
  ifelse(abs(delta_elpd_s) > 2 * se_diff_s, "SÌ", "NO")
))
