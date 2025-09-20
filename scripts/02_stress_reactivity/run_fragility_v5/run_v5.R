# ============================================================
# R script — Modello t eteroschedastico + LOO per soggetto
# ============================================================

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
  library(bayesplot)
  library(ggplot2)
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("sd", "stats")
conflict_prefer("rhat", "posterior")

# ------------------------------------------------------------
# Path e opzioni
# ------------------------------------------------------------
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
stan_file_thet <- here::here(
  "scripts",
  "02_stress_reactivity",
  "run_fragility_v5",
  "subject_loo_normal_rslope_marginal.stan"
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
# 7) Fattori di periodo a livello osservazione (dummy 0/1)
# ------------------------------------------------------------
period <- as.integer(d2$.__per__)
is_pre <- as.integer(period == 2L)
is_post <- as.integer(period == 3L)
# versioni BINARIE che NON vanno scalate (servono a Stan per la scala eterosched.)
is_pre_bin <- is_pre
is_post_bin <- is_post

# ------------------------------------------------------------
# 8) DATASET COERENTE df_mod (una riga per osservazione)
# ------------------------------------------------------------
df_mod <- tibble::tibble(
  y = as.numeric(y),
  subj = as.integer(subj_idx),
  z_naff_b = Xb_exp[, "z_naff_b"],
  z_det_b = Xb_exp[, "z_det_b"],
  z_ant_b = Xb_exp[, "z_ant_b"],
  z_dis_b = Xb_exp[, "z_dis_b"],
  z_psy_b = Xb_exp[, "z_psy_b"],
  female = as.integer(female_exp),
  is_pre = is_pre, # in X (può essere scalato)
  is_post = is_post
)
# Aggiungi le colonne binarie NON scalate per Stan (eteroschedasticità)
df_mod$is_pre_bin <- is_pre_bin
df_mod$is_post_bin <- is_post_bin

# EMA within opzionali → append
if (is.matrix(W) && nrow(W) == nrow(df_mod) && ncol(W) > 0) {
  Wdf <- as_tibble(W)
  names(Wdf) <- make.names(names(Wdf))
  df_mod <- bind_cols(df_mod, Wdf)
} else {
  Wdf <- NULL
}

# Centra/scala predittori: escludi le BIN
pred_cols <- setdiff(names(df_mod), c("y", "subj", "is_pre_bin", "is_post_bin"))
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
#    - RITORNA anche rows_kept per riallineare altre variabili
# ------------------------------------------------------------
build_stan_data <- function(df, x_cols, grainsize_hint = NULL) {
  stopifnot(all(x_cols %in% names(df)))
  X0 <- as.matrix(df[, x_cols, drop = FALSE])
  ok <- is.finite(df$y) & apply(is.finite(X0), 1, all)
  df <- df[ok, , drop = FALSE]
  X <- X0[ok, , drop = FALSE]
  rows_kept <- which(ok)

  subj_dense <- as.integer(factor(df$subj, levels = sort(unique(df$subj))))
  N <- nrow(df)
  S <- max(subj_dense)

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
    } else start_subj[s] <- 1L
  }

  grainsize <- if (is.null(grainsize_hint)) {
    max(
      250L,
      as.integer(N / max(8L, parallel::detectCores(logical = TRUE) / 2))
    )
  } else as.integer(grainsize_hint)

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
    n_per_subj = as.integer(n_per_subj),
    rows_kept = as.integer(rows_kept) # << per riallineare dummies BIN
  )
}

# ------------------------------------------------------------
# 10) Colonne X per A e per B + creazione data list
# ------------------------------------------------------------
# 10) Colonne X per A e B (fissi)
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
x_cols_B <- if (!is.null(Wdf)) c(x_cols_A, names(Wdf)) else x_cols_A

# Costruisci data list base (senza Z)
data_A <- build_stan_data(df_mod, x_cols_A)
data_B <- build_stan_data(df_mod, x_cols_B)

# --- COSTRUISCI Z (random effects design) allineata a ciascun data_? ----
# Scegli i random slopes: inizia con intercept + periodo
make_Z <- function(df, rows_kept, add_ema_slopes = FALSE) {
  Zlist <- list(`(Intercept)` = rep(1, nrow(df)))
  Zlist$is_pre <- df$is_pre
  Zlist$is_post <- df$is_post
  if (add_ema_slopes && !is.null(Wdf)) {
    for (nm in names(Wdf)) Zlist[[nm]] <- df[[nm]]
  }
  Z <- as.matrix(as_tibble(Zlist))
  Z[rows_kept, , drop = FALSE]
}
Z_A <- make_Z(df_mod, data_A$rows_kept, add_ema_slopes = TRUE) # TRUE = anche EMA random (facoltativo)
Z_B <- make_Z(df_mod, data_B$rows_kept, add_ema_slopes = TRUE)

data_A$R <- ncol(Z_A)
data_A$Z <- Z_A
data_B$R <- ncol(Z_B)
data_B$Z <- Z_B

# ------------------------------------------------------------
# 11) Compilazione modello (forza ricompilazione)
# ------------------------------------------------------------
stan_file_rs <- here::here(
  "scripts",
  "02_stress_reactivity",
  "run_fragility_v5",
  "subject_loo_normal_rslope_marginal.stan"
)
mod_rs <- cmdstanr::cmdstan_model(stan_file_rs, force_recompile = TRUE)


# ------------------------------------------------------------
# 12) Fit A/B (MCMC breve)
# ------------------------------------------------------------

fit_A <- mod_rs$sample(
  data = within(data_A, {
    rm(Q, gh_nodes, gh_weights, is_pre, is_post)
  }),
  seed = SEED,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95
)
fit_B <- mod_rs$sample(
  data = within(data_B, {
    rm(Q, gh_nodes, gh_weights, is_pre, is_post)
  }),
  seed = SEED + 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95
)

# ------------------------------------------------------------
# 13) LOO per soggetto (MARGINALIZZATO)
# ------------------------------------------------------------

ll_A_s <- fit_A$draws("log_lik_subj", format = "matrix")
ll_B_s <- fit_B$draws("log_lik_subj", format = "matrix")
loo_A_s <- loo::loo(ll_A_s, r_eff = NA, moment_match = TRUE)
loo_B_s <- loo::loo(ll_B_s, r_eff = NA, moment_match = TRUE)
print(loo::loo_compare(list(B = loo_B_s, A = loo_A_s)))
print(loo::pareto_k_table(loo_A_s))
print(loo::pareto_k_table(loo_B_s))


# ------------------------------------------------------------
# 14) Diagnostiche di convergenza
# ------------------------------------------------------------
posterior::summarise_draws(
  fit_B$draws(c("beta", "sigma", "sigma_b", "nu", "gamma_pre", "gamma_post")),
  "mean",
  "sd",
  "rhat",
  "ess_bulk",
  "ess_tail"
) %>%
  print(n = Inf)

# ------------------------------------------------------------
# 15) PPC generali (densità, QQ, ECDF) — modello t eteroschedastico
# ------------------------------------------------------------
X <- data_B$X
Z <- data_B$Z
subj <- data_B$subj
y <- data_B$y
N <- length(y)
K <- ncol(X)
R <- ncol(Z)
S <- max(subj)

dr <- posterior::as_draws_matrix(fit_B$draws(c("beta", "b_subj", "sigma")))
nd <- min(300, nrow(dr))
i_samp <- sample(seq_len(nrow(dr)), nd)

beta_draws <- dr[i_samp, paste0("beta[", 1:K, "]"), drop = FALSE] # nd x K
sigma_draws <- as.numeric(dr[i_samp, "sigma"])
# Estrai b_subj come matrice nd x (S*R), poi riassembla per indexing
b_mat <- dr[i_samp, grep("^b_subj\\[", colnames(dr)), drop = FALSE] # nd x (S*R)

# Funzione per ottenere (nd x N) di contributo random: (Z * b_subj[subj,]) per ogni draw
re_contrib_all <- function(b_mat, Z, subj, S, R) {
  nd <- nrow(b_mat)
  N <- nrow(Z)
  out <- matrix(NA_real_, nd, N)
  for (d in 1:nd) {
    # ricostruisci B_s (S x R) per il draw d
    Bsr <- matrix(NA_real_, nrow = S, ncol = R)
    # colonne in ordine b_subj[s,r] by r fastest
    # colnames tipiche: b_subj[1,1], b_subj[1,2], ..., b_subj[S,R]
    # usiamo indice posizionale:
    Bsr[] <- as.numeric(b_mat[d, ])
    # contributo per ogni osservazione: row(Z,n) %*% Bsr[subj[n], ]
    for (n in 1:N) {
      out[d, n] <- sum(Z[n, ] * Bsr[subj[n], ])
    }
  }
  out
}

mu_fix_all <- beta_draws %*% t(X) # nd x N
re_all <- re_contrib_all(b_mat, Z, subj, S, R) # nd x N
mu_all <- mu_fix_all + re_all

yrep <- matrix(NA_real_, nd, N)
for (d in 1:nd) yrep[d, ] <- rnorm(N, mu_all[d, ], sigma_draws[d])

# PPC densità globale / QQ / ECDF
bayesplot::ppc_dens_overlay(y, yrep[1:50, ])
bayesplot::ppc_ecdf_overlay(y, yrep[1:50, ])
ppc_qq_obs <- function(y, yrep, ndraws = 50) {
  nd <- min(ndraws, nrow(yrep))
  q_obs <- sort(y)
  df <- do.call(
    rbind,
    lapply(
      seq_len(nd),
      function(d) data.frame(draw = d, q_obs = q_obs, q_rep = sort(yrep[d, ]))
    )
  )
  ggplot(df, aes(q_obs, q_rep, group = draw)) +
    geom_abline(lty = 2) +
    geom_line(alpha = .25) +
    labs(x = "Quantili osservati", y = "Quantili predetti", title = "PPC QQ") +
    theme_minimal()
}
ppc_qq_obs(y, yrep, 50)

# PPC sulla SD entro-soggetto (come il tuo check)
sd_by_s_obs <- tapply(y, subj, sd)
sd_by_s_reps <- apply(yrep, 1, function(yr) tapply(yr, subj, sd)) # S x nd
bayesplot::ppc_dens_overlay(sd_by_s_obs, t(sd_by_s_reps)[1:50, ])

# ricava media/sd usate nello z-score della somma (già calcolate quando hai creato y)
mu_sum <- mean(NA_score)
sd_sum <- sd(NA_score)

# porta y e yrep alla scala sum-score e coarsena a interi 4..28
y_sum <- pmin(28, pmax(4, round(y * sd_sum + mu_sum)))
yrep_sum <- pmin(28, pmax(4, round(yrep * sd_sum + mu_sum))) # yrep è nd x N

bayesplot::ppc_hist(y_sum, yrep_sum[1:50, ])
bayesplot::ppc_bars(y_sum, yrep_sum[1:50, ])

# medie per soggetto
ybar_obs <- tapply(y, subj, mean)
ybar_reps <- apply(yrep, 1, function(yr) tapply(yr, subj, mean)) # S x nd
bayesplot::ppc_dens_overlay(ybar_obs, t(ybar_reps)[1:50, ])

# sd entro-soggetto separata per periodo
is_pre_vec <- as.integer(data_B$is_pre)
is_post_vec <- as.integer(data_B$is_post)

sd_pre_obs <- tapply(y[is_pre_vec == 1], subj[is_pre_vec == 1], sd)
sd_post_obs <- tapply(y[is_post_vec == 1], subj[is_post_vec == 1], sd)

sd_pre_rep <- apply(
  yrep[, is_pre_vec == 1, drop = FALSE],
  1,
  function(yr) tapply(yr, subj[is_pre_vec == 1], sd)
)
sd_post_rep <- apply(
  yrep[, is_post_vec == 1, drop = FALSE],
  1,
  function(yr) tapply(yr, subj[is_post_vec == 1], sd)
)

bayesplot::ppc_dens_overlay(sd_pre_obs, t(sd_pre_rep)[1:50, ])
bayesplot::ppc_dens_overlay(sd_post_obs, t(sd_post_rep)[1:50, ])


posterior::summarise_draws(
  fit_B$draws(c("beta", "sigma", "tau", "L_Omega")),
  "mean",
  "sd",
  "rhat",
  "ess_bulk",
  "ess_tail"
) %>%
  print(n = Inf)
