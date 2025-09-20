# ------------------------------------------------------------
# Librerie
# ------------------------------------------------------------
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
  library(stringr)
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
  "run_fragility_v2",
  "subject_loo_reducesum_normal.stan"
)
SEED <- 20250919

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
fem_by_subj <- tapply(
  female_by_row,
  d$.__subj__,
  function(v) mean(v, na.rm = TRUE)
)
fem_by_subj[is.na(fem_by_subj)] <- 1
I_all <- max(d$.__subj__)
female_vec <- integer(I_all)
female_vec[as.integer(names(fem_by_subj))] <- as.integer(fem_by_subj >= 0.5)

# ------------------------------------------------------------
# 4) Item 1..7 + outcome y = Negative Affect (z)
#    NA = (7 - happy) + (7 - satisfied) + sad + angry
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

keep <- with(
  d,
  is.finite(.__happy__) &
    is.finite(.__sad__) &
    is.finite(.__satisfied__) &
    is.finite(.__angry__)
)
d2 <- d[keep, , drop = FALSE]
stopifnot(nrow(d2) > 0)

# Outcome (continua, approx. gaussiana)
NA_score <- (7 - d2$.__happy__) +
  (7 - d2$.__satisfied__) +
  d2$.__sad__ +
  d2$.__angry__
y <- z_(NA_score)

# ------------------------------------------------------------
# 5) Predittori BETWEEN: PID-5 baseline (z a livello soggetto), female
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
  select(.__subj__, all_of(base_cols)) %>%
  group_by(.__subj__) %>%
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
  arrange(.__subj__)

Xb <- as.matrix(base_by_subj[, base_cols, drop = FALSE]) %>% apply(2, z_)
if (nrow(Xb) < max(d2$.__subj__)) {
  Xb_full <- matrix(NA_real_, nrow = max(d2$.__subj__), ncol = 5)
  Xb_full[base_by_subj$.__subj__, ] <- Xb
  Xb <- apply(Xb_full, 2, z_)
}
Xb[!is.finite(Xb)] <- 0
colnames(Xb) <- c("z_naff_b", "z_det_b", "z_ant_b", "z_dis_b", "z_psy_b")

female <- as.integer(female_vec[seq_len(max(d2$.__subj__))])

# Espandi sui beep
subj_idx <- as.integer(d2$.__subj__)
Xb_exp <- Xb[subj_idx, , drop = FALSE]
female_exp <- female[subj_idx]

# ------------------------------------------------------------
# 6) Predittori WITHIN opzionali (EMA “dinamici” z-score) — se disponibili
#    (qui uso colonne tipo ema_detachment, ema_antagonism, ecc. se esistono)
# ------------------------------------------------------------
ema_within <- list(
  det = find_first(
    d2,
    c("ema_detachment", "pid5_detachment", "detachment", "z_det", "det")
  ),
  ant = find_first(
    d2,
    c("ema_antagonism", "pid5_antagonism", "antagonism", "z_ant", "ant")
  ),
  dis = find_first(
    d2,
    c(
      "ema_disinhibition",
      "pid5_disinhibition",
      "disinhibition",
      "z_dis",
      "dis"
    )
  ),
  psy = find_first(
    d2,
    c("ema_psychoticism", "pid5_psychoticism", "psychoticism", "z_psy", "psy")
  )
)

W_list <- list()
if (!is.null(ema_within$det))
  W_list$w_det <- z_(as.numeric(d2[[ema_within$det]]))
if (!is.null(ema_within$ant))
  W_list$w_ant <- z_(as.numeric(d2[[ema_within$ant]]))
if (!is.null(ema_within$dis))
  W_list$w_dis <- z_(as.numeric(d2[[ema_within$dis]]))
if (!is.null(ema_within$psy))
  W_list$w_psy <- z_(as.numeric(d2[[ema_within$psy]]))

W <- if (length(W_list)) as.matrix(as_tibble(W_list)) else
  matrix(numeric(0), nrow = nrow(d2), ncol = 0)

# ------------------------------------------------------------
# 7) Fattori di periodo a livello osservazione (dummy, senza intercetta)
# ------------------------------------------------------------
period <- as.integer(d2$.__per__)
per_levels <- sort(unique(period))
# Dummy per: baseline=ref (nessuna colonna), pre e post come dummies se esistono
P_cols <- NULL
if (2L %in% per_levels) P_cols <- cbind(P_cols, as.integer(period == 2L))
if (3L %in% per_levels) P_cols <- cbind(P_cols, as.integer(period == 3L))
if (is.null(P_cols)) {
  P <- matrix(numeric(0), nrow = nrow(d2), ncol = 0)
  colnames(P) <- character(0)
} else {
  P <- as.matrix(P_cols)
  colnames(P) <- c(
    if (2L %in% per_levels) "is_pre" else NULL,
    if (3L %in% per_levels) "is_post" else NULL
  )
}

# ------------------------------------------------------------
# 8) DATASET COERENTE df_mod (una riga per osservazione)
# ------------------------------------------------------------
is_pre <- as.integer(period == 2L)
is_post <- as.integer(period == 3L)

df_mod <- tibble::tibble(
  y = as.numeric(y),
  subj = as.integer(subj_idx),
  z_naff_b = Xb_exp[, "z_naff_b"],
  z_det_b = Xb_exp[, "z_det_b"],
  z_ant_b = Xb_exp[, "z_ant_b"],
  z_dis_b = Xb_exp[, "z_dis_b"],
  z_psy_b = Xb_exp[, "z_psy_b"],
  female = as.integer(female_exp),
  is_pre = is_pre,
  is_post = is_post
)

# EMA within opzionali
Wdf <- NULL
if (exists("W") && is.matrix(W) && nrow(W) == nrow(df_mod) && ncol(W) > 0) {
  Wdf <- as_tibble(W)
  names(Wdf) <- make.names(names(Wdf))
  df_mod <- bind_cols(df_mod, Wdf)
}

# Centra/scala predittori (NA -> 0)
pred_cols <- setdiff(names(df_mod), c("y", "subj"))
df_mod[pred_cols] <- lapply(df_mod[pred_cols], function(x) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  x <- if (is.finite(s) && s > 0) (x - m) / s else (x - m)
  x[!is.finite(x)] <- 0
  x
})

# ------------------------------------------------------------
# 9) Funzione generica per costruire i data list (A o B)
#    - Filtra righe in base alla completezza di X
#    - Rimappa subj -> 1..S densi dopo il filtro
#    - Ricrea obs_by_subj/start_subj/n_per_subj
# ------------------------------------------------------------
build_stan_data <- function(df, x_cols, grainsize_hint = NULL) {
  stopifnot(all(x_cols %in% names(df)))
  X <- as.matrix(df[, x_cols, drop = FALSE])

  # Filtra righe: y finita + tutte le colonne X finite
  ok <- is.finite(df$y) & apply(is.finite(X), 1, all)
  df <- df[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]

  # Rimappa soggetti a 1..S DENSI DOPO il filtro
  subj_dense <- as.integer(factor(df$subj, levels = sort(unique(df$subj))))
  N <- nrow(df)
  S <- max(subj_dense)

  # Indici per soggetto
  n_per_subj <- tabulate(subj_dense, nbins = S)
  start_subj <- integer(S)
  obs_by_subj <- integer(N)
  pos <- 1L
  for (s in 1:S) {
    ns <- n_per_subj[s]
    if (ns > 0) {
      ids <- which(subj_dense == s)
      start_subj[s] <- pos
      obs_by_subj[pos:(pos + ns - 1L)] <- ids
      pos <- pos + ns
    } else {
      start_subj[s] <- 1L
    }
  }

  # grainsize robusto
  if (is.null(grainsize_hint)) {
    grainsize <- max(
      250L,
      as.integer(N / max(8L, parallel::detectCores(logical = TRUE) / 2))
    )
  } else {
    grainsize <- as.integer(grainsize_hint)
  }

  # Check difensivi
  stopifnot(
    nrow(X) == N,
    length(subj_dense) == N,
    identical(sort(unique(subj_dense)), seq_len(S))
  )

  list(
    N = N,
    S = S,
    K = ncol(X),
    y = as.numeric(df$y),
    X = unname(X),
    subj = as.integer(subj_dense),
    grainsize = grainsize,
    obs_by_subj = as.integer(obs_by_subj),
    start_subj = as.integer(start_subj),
    n_per_subj = as.integer(n_per_subj)
  )
}

# ------------------------------------------------------------
# 10) Colonne X per A e per B + creazione data list
# ------------------------------------------------------------
x_cols_A <- c(
  "z_naff_b",
  "z_det_b",
  "z_ant_b",
  "z_dis_b",
  "z_psy_b",
  "female",
  "is_pre",
  "is_post"
)

if (!is.null(Wdf)) {
  x_cols_B <- c(x_cols_A, names(Wdf))
} else {
  x_cols_B <- x_cols_A
}

data_A <- build_stan_data(df_mod, x_cols_A)
data_B <- build_stan_data(df_mod, x_cols_B)

cat(sprintf(
  "A: N=%d, S=%d, K=%d | B: N=%d, S=%d, K=%d\n",
  data_A$N,
  data_A$S,
  data_A$K,
  data_B$N,
  data_B$S,
  data_B$K
))
stopifnot(
  length(data_A$y) == data_A$N,
  nrow(data_A$X) == data_A$N,
  length(data_A$subj) == data_A$N,
  max(data_A$subj) == data_A$S,
  length(data_B$y) == data_B$N,
  nrow(data_B$X) == data_B$N,
  length(data_B$subj) == data_B$N,
  max(data_B$subj) == data_B$S
)

str(list(
  N_A = data_A$N,
  len_y_A = length(data_A$y),
  nrow_X_A = nrow(data_A$X),
  len_subj_A = length(data_A$subj),
  range_subj_A = range(data_A$subj),
  N_B = data_B$N,
  len_y_B = length(data_B$y),
  nrow_X_B = nrow(data_B$X),
  len_subj_B = length(data_B$subj),
  range_subj_B = range(data_B$subj)
))


# ------------------------------------------------------------
# 11) Compilazione modello
# ------------------------------------------------------------
stopifnot(file.exists(stan_file_path))
cat(readLines(stan_file_path, n = 40), sep = "\n")
mod <- cmdstan_model(
  stan_file_path,
  cpp_options = list(stan_threads = TRUE),
  force_recompile = TRUE
)

# ------------------------------------------------------------
# 12) Fit A/B (scegli VI *oppure* MCMC breve)
# ------------------------------------------------------------

stan_file_path_no_rs <- here::here(
  "scripts",
  "02_stress_reactivity",
  "run_fragility_v2",
  "subject_loo_normal_no_reduce_sum.stan"
)
mod_nr <- cmdstan_model(stan_file_path_no_rs, force_recompile = TRUE)

# VI (se vuoi restare su VI)
fit_A <- mod_nr$variational(
  data = data_A,
  seed = SEED,
  algorithm = "meanfield",
  elbo_samples = 1000,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 0
)
fit_B <- mod_nr$variational(
  data = data_B,
  seed = SEED + 1,
  algorithm = "meanfield",
  elbo_samples = 1000,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 0
)

# Oppure MCMC breve, consigliato:
# fit_A <- mod_nr$sample(data=data_A, seed=SEED, chains=4, parallel_chains=4,
#                        iter_warmup=500, iter_sampling=500, adapt_delta=0.92)
# fit_B <- mod_nr$sample(data=data_B, seed=SEED+1, chains=4, parallel_chains=4,
#                        iter_warmup=500, iter_sampling=500, adapt_delta=0.92)

# LOO per soggetto
ll_A_s <- fit_A$draws("log_lik_subj", format = "matrix")
ll_B_s <- fit_B$draws("log_lik_subj", format = "matrix")
loo_A_s <- loo::loo(ll_A_s, r_eff = NA, moment_match = TRUE)
loo_B_s <- loo::loo(ll_B_s, r_eff = NA, moment_match = TRUE)
loo::loo_compare(list(B = loo_B_s, A = loo_A_s))

# ------------------------------------------------------------
# 13.1) Convergenza e dimensione campionaria effettiva
# ------------------------------------------------------------

posterior::summarise_draws(
  fit_B$draws(c("beta", "sigma", "sigma_b")),
  "mean",
  "sd",
  "rhat",
  "ess_bulk",
  "ess_tail"
) %>%
  print(n = Inf)


# ------------------------------------------------------------
# 13.2) PPC generali (densità, QQ, ECDF)
# ------------------------------------------------------------

library(bayesplot)

# --- prepara oggetti utili ---
X <- data_B$X
subj <- data_B$subj
y <- data_B$y
N <- length(y)
K <- ncol(X)
S <- max(subj)

# Estrai un sottoinsieme di draw
nd <- 300
dr <- posterior::as_draws_matrix(fit_B$draws(c("beta", "b_subj", "sigma")))
i_samp <- sample(seq_len(nrow(dr)), min(nd, nrow(dr)))

beta_draws <- dr[i_samp, paste0("beta[", 1:K, "]"), drop = FALSE] # nd x K
bsubj_draws <- dr[i_samp, paste0("b_subj[", 1:S, "]"), drop = FALSE] # nd x S
sigma_draws <- as.numeric(dr[i_samp, "sigma"])

# Costruisci yrep (nd x N)
yrep <- matrix(NA_real_, nrow = length(i_samp), ncol = N)
for (d in seq_along(i_samp)) {
  beta_d <- as.numeric(beta_draws[d, , drop = TRUE]) # vettore K
  mu_fix <- as.numeric(X %*% beta_d) # vettore N
  b_i <- as.numeric(bsubj_draws[d, subj, drop = TRUE]) # vettore N
  mu_d <- mu_fix + b_i
  yrep[d, ] <- rnorm(N, mean = mu_d, sd = sigma_draws[d])
}

# --- PPC densità complessiva ---
bayesplot::ppc_dens_overlay(y, yrep[1:50, ])

# --- PPC ECDF overlay (calibrazione globale) ---
bayesplot::ppc_ecdf_overlay(y, yrep[1:50, ])

# Assumo che tu abbia già creato: X, subj, y, yrep (nd x N)

# --- PPC-QQ: yrep vs y (replica "ppc_qq") ---
ppc_qq_obs <- function(y, yrep, ndraws = 50) {
  stopifnot(is.numeric(y), is.matrix(yrep))
  nd <- min(ndraws, nrow(yrep))
  N <- length(y)
  # quantili osservati
  q_obs <- sort(y) # quantili empirici con p = (1..N)/(N+1) impliciti
  # dataframe lungo con quantili delle repliche
  dfs <- lapply(seq_len(nd), function(d) {
    data.frame(
      draw = paste0("draw_", d),
      q_obs = q_obs,
      q_rep = sort(yrep[d, ])
    )
  })
  df <- do.call(rbind, dfs)
  ggplot2::ggplot(df, ggplot2::aes(x = q_obs, y = q_rep, group = draw)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_line(alpha = 0.25) +
    ggplot2::labs(
      x = "Quantili osservati",
      y = "Quantili predetti",
      title = "PPC QQ: yrep vs y (linee = draw posteriori)"
    ) +
    ggplot2::theme_minimal()
}

# Esegui
ppc_qq_obs(y, yrep, ndraws = 50)


# ------------------------------------------------------------
# 13.3) PPC “per soggetto” (statistiche raggruppate)
# ------------------------------------------------------------

# Media per soggetto dei dati reali
ybar_by_s <- tapply(y, subj, mean)

# Media per soggetto delle repliche
yrep_by_s <- apply(yrep, 1, function(yr) tapply(yr, subj, mean)) # S x nd
yrep_by_s <- t(yrep_by_s) # nd x S

# Overlay distribuzioni delle medie per soggetto (campiona ~30 soggetti per visibilità)
set.seed(1)
s_plot <- sort(sample(seq_len(S), size = min(30, S)))
par(mfrow = c(5, 6), mar = c(2, 2, 2, 1))
for (s in s_plot) {
  bayesplot::ppc_dens_overlay(
    y = ybar_by_s[s],
    yrep = matrix(yrep_by_s[, s], nrow = nd, ncol = 1)
  )
  mtext(paste("Soggetto", s), line = 0.2, cex = 0.8)
}
par(mfrow = c(1, 1))
