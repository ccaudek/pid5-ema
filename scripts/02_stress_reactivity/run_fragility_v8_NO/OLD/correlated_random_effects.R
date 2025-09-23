# Overview ----------------------------------------------------------------
# Associated project: PID-5 and EMA
# Script purpose: Subject-level LOO with random slopes & subject-level
#   heteroskedasticity. Model: Gaussian outcomes, random intercepts + random
#   slopes, per-subject scale (sigma_s). LOO at subject level computed via
#   Gauss–Hermite marginalization (log_lik_subj_gh).
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-09-20
# Last update:
# Status: In progress
# Notes:

# ------------------------------------------------------------
# Libraries
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
# Paths & options
# ------------------------------------------------------------
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
# Stan file: random slopes + per-subject scale (sigma_s), with GH LOO per subject
stan_file_rs_sigmaS <- here::here(
  "scripts",
  "02_stress_reactivity",
  "run_fragility_v7",
  "subject_loo_normal_rslope_sigmaS.stan"
)
SEED <- 20250919
set.seed(SEED) # ensure reproducible subsampling in PPC

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------
# Convert vectors to 1..7 ordinal coding; robust to messy input
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

# z-score with NA safety; if sd=0 returns centered
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) x - m else (x - m) / s
}

# find first matching column ignoring case
find_first <- function(df, candidates) {
  nm <- names(df)
  hit <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hit) == 0) NULL else hit[1]
}

# ------------------------------------------------------------
# 1) Load data + normalize baseline variable names
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
# 2) Subject and period indices
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
# 3) Sex → female (0/1) at subject level
#    (robustly detect, then majority vote per subject)
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
  warning("Sex column not found: set female=1 for all (placeholder).")
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
# 4) Items 1..7 + outcome y = Negative Affect (z)
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
# 5) BETWEEN predictors: PID-5 baseline (subject-level z), female
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

# Expand to beeps
subj_idx <- as.integer(d2$.__subj__)
Xb_exp <- Xb[subj_idx, , drop = FALSE]
female_exp <- female[subj_idx]

# ------------------------------------------------------------
# 6) Optional WITHIN predictors (dynamic EMA z-scores), if present
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
# 7) Period factors at the observation level (0/1 dummies)
#    Also keep binary copies (not scaled) for PPC by subject×period
# ------------------------------------------------------------
period <- as.integer(d2$.__per__)
is_pre <- as.integer(period == 2L)
is_post <- as.integer(period == 3L)
is_pre_bin <- is_pre
is_post_bin <- is_post

# ------------------------------------------------------------
# 8) Cohesive long data frame (one row per observation)
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
  is_pre = is_pre,
  is_post = is_post
)
# Add non-scaled binaries for PPC grouping
df_mod$is_pre_bin <- is_pre_bin
df_mod$is_post_bin <- is_post_bin

# Append EMA within (if present)
if (is.matrix(W) && nrow(W) == nrow(df_mod) && ncol(W) > 0) {
  Wdf <- as_tibble(W)
  names(Wdf) <- make.names(names(Wdf))
  df_mod <- bind_cols(df_mod, Wdf)
} else {
  Wdf <- NULL
}

# Center/scale predictors; exclude binary copies
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
# 9) Build data lists for Stan (A or B)
#    - Filter rows with complete X
#    - Remap subj -> dense 1..S after filtering
#    - Recreate obs_by_subj / start_subj / n_per_subj
#    - Return rows_kept for downstream alignment (Z, PPC grouping)
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
    rows_kept = as.integer(rows_kept)
  )
}

# ------------------------------------------------------------
# 10) Choose X columns for models A and B; create data lists
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
x_cols_B <- if (!is.null(Wdf)) c(x_cols_A, names(Wdf)) else x_cols_A

data_A <- build_stan_data(df_mod, x_cols_A)
data_B <- build_stan_data(df_mod, x_cols_B)

# --- Build Z (random-effects design) aligned with rows_kept for each dataset ---
# Random slopes: start with Intercept + period; optionally EMA slopes too
make_Z <- function(df, rows_kept, add_ema_slopes = TRUE) {
  Zlist <- list(`(Intercept)` = rep(1, nrow(df)))
  Zlist$is_pre <- df$is_pre
  Zlist$is_post <- df$is_post
  if (add_ema_slopes && !is.null(Wdf)) {
    for (nm in names(Wdf)) Zlist[[nm]] <- df[[nm]]
  }
  Z <- as.matrix(as_tibble(Zlist))
  Z[rows_kept, , drop = FALSE]
}
Z_A <- make_Z(df_mod, data_A$rows_kept, add_ema_slopes = TRUE)
Z_B <- make_Z(df_mod, data_B$rows_kept, add_ema_slopes = TRUE)

data_A$R <- ncol(Z_A)
data_A$Z <- Z_A
data_B$R <- ncol(Z_B)
data_B$Z <- Z_B

# Gauss–Hermite nodes/weights for subject-level marginal likelihood
Q <- 20L
if (requireNamespace("fastGHQuad", quietly = TRUE)) {
  gh <- fastGHQuad::gaussHermiteData(Q)
  gh_nodes <- as.numeric(gh$x)
  gh_weights <- as.numeric(gh$w)
} else {
  gh <- statmod::gauss.quad(Q, kind = "hermite")
  gh_nodes <- as.numeric(gh$nodes)
  gh_weights <- as.numeric(gh$weights)
}
data_A$Q <- Q
data_A$gh_nodes <- gh_nodes
data_A$gh_weights <- gh_weights
data_B$Q <- Q
data_B$gh_nodes <- gh_nodes
data_B$gh_weights <- gh_weights

cat(sprintf(
  "A: N=%d, S=%d, K=%d, R=%d | B: N=%d, S=%d, K=%d, R=%d\n",
  data_A$N,
  data_A$S,
  data_A$K,
  data_A$R,
  data_B$N,
  data_B$S,
  data_B$K,
  data_B$R
))

# ------------------------------------------------------------
# 11) Compile Stan model
# ------------------------------------------------------------
stopifnot(file.exists(stan_file_rs_sigmaS))
mod_rs_sig <- cmdstan_model(stan_file_rs_sigmaS, force_recompile = FALSE)

# ------------------------------------------------------------
# 12) Fit A/B (short MCMC)
# ------------------------------------------------------------
fit_A <- mod_rs_sig$sample(
  data = data_A,
  seed = SEED + 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95
)
fit_B <- mod_rs_sig$sample(
  data = data_B,
  seed = SEED + 2,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95
)

# Saving fitted model objects

# Load posterior draws into the fitted model object and omit other output.
fit_A$draws()
# Save the object to a file.
qs::qsave(x = fit_A, file = "fit_A.qs")
# Read the object.
# fit_A <- qs::qread("fit_A.qs")

fit_B$draws()
qs::qsave(x = fit_B, file = "fit_B.qs")
# Read the object.
# fit_B <- qs::qread("fit_B.qs")

# ------------------------------------------------------------
# 13) Subject-level LOO (marginalized)
# ------------------------------------------------------------
# Compute chain_id-aware r_eff to reduce bias
# A
ll_A_s <- fit_A$draws("log_lik_subj_gh", format = "matrix") # (draws_total x S)
n_chains_A <- fit_A$num_chains()
n_draws_A <- nrow(ll_A_s)
stopifnot(n_draws_A %% n_chains_A == 0)
n_iter_A <- n_draws_A / n_chains_A
chain_id_A <- rep(seq_len(n_chains_A), each = n_iter_A)
r_eff_A <- loo::relative_eff(exp(ll_A_s), chain_id = chain_id_A)

# B
ll_B_s <- fit_B$draws("log_lik_subj_gh", format = "matrix")
n_chains_B <- fit_B$num_chains()
n_draws_B <- nrow(ll_B_s)
stopifnot(n_draws_B %% n_chains_B == 0)
n_iter_B <- n_draws_B / n_chains_B
chain_id_B <- rep(seq_len(n_chains_B), each = n_iter_B)
r_eff_B <- loo::relative_eff(exp(ll_B_s), chain_id = chain_id_B)

# PSIS-LOO with moment matching
loo_A_s <- loo::loo(ll_A_s, r_eff = r_eff_A, moment_match = TRUE)
loo_B_s <- loo::loo(ll_B_s, r_eff = r_eff_B, moment_match = TRUE)

print(loo::loo_compare(list(B = loo_B_s, A = loo_A_s)))
#  elpd_diff se_diff
# B    0.0       0.0
# A -132.9       9.8
loo::pareto_k_table(loo_A_s)
# All Pareto k estimates are good (k < 0.7).
loo::pareto_k_table(loo_B_s)
# All Pareto k estimates are good (k < 0.7).

# Persist LOO objects for reproducibility
# saveRDS(loo_A_s, file = "loo_A_subject.rds")
# saveRDS(loo_B_s, file = "loo_B_subject.rds")

# ------------------------------------------------------------
# 14) Convergence diagnostics (concise)
# ------------------------------------------------------------
posterior::summarise_draws(
  fit_B$draws(c("beta", "sigma", "sd_delta")),
  "mean",
  "sd",
  "rhat",
  "ess_bulk",
  "ess_tail"
) %>%
  print(n = Inf)
# variable    mean      sd  rhat ess_bulk ess_tail
# <chr>      <dbl>   <dbl> <dbl>    <dbl>    <dbl>
# 1 beta[1]   0.181  0.0706  1.00     1613.    2154.
# 2 beta[2]  -0.125  0.0550  1.00     1744.    2460.
# 3 beta[3]   0.0637 0.0542  1.00     1410.    2298.
# 4 beta[4]  -0.0193 0.0503  1.00     1207.    2096.
# 5 beta[5]  -0.115  0.0633  1.00     1336.    2521.
# 6 beta[6]   0.0292 0.0330  1.00     1410.    2272.
# 7 beta[7]   0.0713 0.00905 1.00     4205.    3181.
# 8 beta[8]  -0.0333 0.00940 1.00     5373.    3267.
# 9 beta[9]   0.428  0.0236  1.00     2810.    3234.
# 10 beta[10] -0.127  0.0203  1.00     3064.    2651.
# 11 beta[11]  0.256  0.0180  1.00     3468.    3008.
# 12 beta[12] -0.0197 0.0183  1.000    4190.    3652.
# 13 sigma     0.646  0.0150  1.00     1084.    1961.
# 14 sd_delta  0.312  0.0181  1.00     1550.    2588.

# ------------------------------------------------------------
# 15) Posterior Predictive Checks (global + targeted)
#     - Global: density, ECDF, QQ
#     - Targeted: within-subject SD; subject×period SD; subject means
# ------------------------------------------------------------

# Helper: reconstruct random-effects contribution for each observation
# b_mat:      nd x (S*R), columns named "b_subj[s,r]"
# Z:          N x R  (random-effects design)
# subj:       length-N vector with indices 1..S
# Reorders columns by r then s if names are present, to match matrix fill order
re_contrib_all <- function(b_mat, Z, subj, S, R, b_colnames = NULL) {
  stopifnot(is.matrix(b_mat), is.matrix(Z), length(subj) == nrow(Z))
  if (!is.null(b_colnames)) {
    rx <- regexec("^b_subj\\[(\\d+),(\\d+)\\]$", b_colnames)
    mt <- regmatches(b_colnames, rx)
    ok <- vapply(mt, function(x) length(x) == 3, TRUE)
    if (all(ok)) {
      s_idx <- as.integer(vapply(mt, `[`, "", 2))
      r_idx <- as.integer(vapply(mt, `[`, "", 3))
      ord <- order(r_idx, s_idx) # order by r, then s
      b_mat <- b_mat[, ord, drop = FALSE]
    }
  }
  nd <- nrow(b_mat)
  N <- nrow(Z)
  out <- matrix(NA_real_, nd, N)
  for (d in seq_len(nd)) {
    Bsr <- matrix(as.numeric(b_mat[d, ]), nrow = S, ncol = R, byrow = FALSE)
    B_obs <- Bsr[subj, , drop = FALSE] # N x R (row picks subject for each obs)
    out[d, ] <- rowSums(Z * B_obs) # elementwise row dot product
  }
  out
}

# Extract objects from current environment
X <- data_B$X
Z <- data_B$Z
subj <- data_B$subj
y <- data_B$y
N <- length(y)
K <- ncol(X)
R <- ncol(Z)
S <- max(subj)

# Draws needed for PPC
dr <- posterior::as_draws_matrix(fit_B$draws(c(
  "beta",
  "b_subj",
  "sigma",
  "sd_delta",
  "delta_raw"
)))
nd <- min(300, nrow(dr)) # number of posterior draws for PPC
set.seed(SEED + 42)
i_samp <- sample(seq_len(nrow(dr)), nd)

beta_draws <- dr[i_samp, paste0("beta[", 1:K, "]"), drop = FALSE]
sigma_draws <- as.numeric(dr[i_samp, "sigma"])
sdD_draws <- as.numeric(dr[i_samp, "sd_delta"])
delta_mat <- dr[i_samp, grep("^delta_raw\\[", colnames(dr)), drop = FALSE] # nd x S
b_cols <- grep("^b_subj\\[", colnames(dr), value = TRUE)
b_mat <- dr[i_samp, b_cols, drop = FALSE] # nd x (S*R)

# Fixed part and random-effects contribution
mu_fix_all <- beta_draws %*% t(X) # nd x N
re_all <- re_contrib_all(b_mat, Z, subj, S, R, b_colnames = b_cols) # nd x N
mu_all <- mu_fix_all + re_all

# Observation-specific sigma: sigma_s per subject, broadcast to observations
sigma_obs <- matrix(NA_real_, nd, N)
for (d in 1:nd) {
  sigma_s <- sigma_draws[d] * exp(sdD_draws[d] * as.numeric(delta_mat[d, 1:S]))
  sigma_obs[d, ] <- sigma_s[subj]
}

# Replicated datasets
yrep <- matrix(NA_real_, nd, N)
for (d in 1:nd) yrep[d, ] <- rnorm(N, mu_all[d, ], sigma_obs[d, ])

# --- PPC: Global overlays ---
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
  ggplot2::ggplot(df, ggplot2::aes(q_obs, q_rep, group = draw)) +
    ggplot2::geom_abline(lty = 2) +
    ggplot2::geom_line(alpha = .25) +
    ggplot2::labs(
      x = "Observed quantiles",
      y = "Replicated quantiles",
      title = "PPC QQ"
    ) +
    ggplot2::theme_minimal()
}
ppc_qq_obs(y, yrep, 50)

# --- PPC: Within-subject SD ---
sd_by_s_obs <- tapply(y, subj, sd)
sd_by_s_reps <- apply(yrep, 1, function(yr) tapply(yr, subj, sd)) # S x nd
bayesplot::ppc_dens_overlay(sd_by_s_obs, t(sd_by_s_reps)[1:50, ])

# --- PPC: Subject means ---
ybar_obs <- tapply(y, subj, mean)
ybar_reps <- apply(yrep, 1, function(yr) tapply(yr, subj, mean))
bayesplot::ppc_dens_overlay(ybar_obs, t(ybar_reps)[1:50, ])

# --- PPC: SD by subject × period (filter groups with n >= 2) ---
per <- df_mod$is_pre_bin + 2L * df_mod$is_post_bin
per <- per[data_B$rows_kept] # align to data_B
grp <- interaction(subj, per, drop = TRUE)
n_by_grp <- tapply(y, grp, function(v) sum(is.finite(v)))
keep_lvls <- names(n_by_grp)[n_by_grp >= 2L]
f_keep <- factor(grp, levels = keep_lvls)
idx_keep <- !is.na(f_keep)

sd_by_sp_obs <- tapply(y[idx_keep], f_keep[idx_keep], sd)

G <- length(keep_lvls)
sd_by_sp_reps <- t(vapply(
  X = seq_len(nrow(yrep)),
  FUN = function(d) tapply(yrep[d, idx_keep], f_keep[idx_keep], sd),
  FUN.VALUE = numeric(G)
))
bayesplot::ppc_dens_overlay(sd_by_sp_obs, sd_by_sp_reps[1:50, ])

# ------------------------------------------------------------
# 16) Session info (for reproducibility)
# ------------------------------------------------------------
print(sessionInfo())
