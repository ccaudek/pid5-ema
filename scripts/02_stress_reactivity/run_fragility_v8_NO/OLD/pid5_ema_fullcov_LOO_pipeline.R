# Overview ----------------------------------------------------------------
# Associated project: PID-5 and EMA
# Script purpose: Subject-level LOO with random slopes & subject-level
#   heteroskedasticity. THIS VERSION targets the *correlated random-effects*
#   Stan model (with LKJ prior) where the subject-wise marginal log-likelihood
#   is computed exactly in generated quantities (no GH quadrature).
#
# Key changes from previous script:
#   - Data list now matches the new Stan model: {N, S, K, R, X, Z, sid, y}.
#   - Removed Gauss–Hermite nodes/weights (no longer used).
#   - LOO now reads 'subj_loglik' (vector[S]) instead of 'log_lik_subj_gh'.
#   - PPC draws now use: beta, b (R x S transformed parameter), sigma_s (S-vector).
#   - Random-effects contribution rebuilt from 'b' with the Z design and subject ids.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-09-22 (adapted to full-covariance RE model)
# Status: In progress
# Notes: Assumes your Stan file implements the correlated RE structure and
#        computes 'subj_loglik' and exposes transformed parameters 'b' and 'sigma_s'.

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

# Stan file: FULL-COV random slopes + per-subject scale (sigma_s);
# *subject-wise* exact marginal LOO provided as 'subj_loglik' (length S).
stan_file_fullcov <- here::here(
  "scripts",
  "02_stress_reactivity",
  "run_fragility_v7",
  "subject_loo_normal_rslope_sigmaS.stan"  # <- same filename, but updated contents
)

SEED <- 20250922
set.seed(SEED)

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
# 3) Sex → female (0/1) at subject level (majority vote)
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
# 9) Build data lists for Stan
#    - Filter rows with complete X
#    - Remap subj -> dense 1..S after filtering
# ------------------------------------------------------------
build_stan_data <- function(df, x_cols) {
  stopifnot(all(x_cols %in% names(df)))
  X0 <- as.matrix(df[, x_cols, drop = FALSE])
  ok <- is.finite(df$y) & apply(is.finite(X0), 1, all)
  df <- df[ok, , drop = FALSE]
  X <- X0[ok, , drop = FALSE]
  rows_kept <- which(ok)

  sid <- as.integer(factor(df$subj, levels = sort(unique(df$subj))))
  N <- nrow(df)
  S <- max(sid)

  stopifnot(
    nrow(X) == N,
    length(sid) == N,
    identical(sort(unique(sid)), seq_len(S))
  )

  list(
    N = N,
    S = S,
    K = ncol(X),
    y = as.numeric(df$y),
    X = unname(X),
    sid = as.integer(sid),  # <-- matches Stan 'sid'
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
# Random slopes: Intercept + period; optionally EMA slopes too
make_Z <- function(df, rows_kept, add_ema_slopes = TRUE) {
  Zlist <- list(`(Intercept)` = rep(1, nrow(df)))
  Zlist$is_pre <- df$is_pre
  Zlist$is_post <- df$is_post
  if (add_ema_slopes && !is.null(Wdf)) {
    for (nm in names(Wdf)) Zlist[[nm]] <- df[[nm]]
  }
  Z_full <- as.matrix(as_tibble(Zlist))
  Z_full[rows_kept, , drop = FALSE]
}

Z_A <- make_Z(df_mod, data_A$rows_kept, add_ema_slopes = TRUE)
Z_B <- make_Z(df_mod, data_B$rows_kept, add_ema_slopes = TRUE)

data_A$R <- ncol(Z_A); data_A$Z <- Z_A
data_B$R <- ncol(Z_B); data_B$Z <- Z_B

cat(sprintf(
  "A: N=%d, S=%d, K=%d, R=%d | B: N=%d, S=%d, K=%d, R=%d\n",
  data_A$N, data_A$S, data_A$K, data_A$R,
  data_B$N, data_B$S, data_B$K, data_B$R
))

# ------------------------------------------------------------
# 11) Compile Stan model
# ------------------------------------------------------------
stopifnot(file.exists(stan_file_fullcov))
mod_fullcov <- cmdstan_model(stan_file_fullcov, force_recompile = FALSE)

# ------------------------------------------------------------
# 12) Fit A/B (short MCMC)
# ------------------------------------------------------------
fit_A <- mod_fullcov$sample(
  data = within(data_A, rm(rows_kept)),  # Stan doesn't need rows_kept
  seed = SEED + 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95
)

fit_B <- mod_fullcov$sample(
  data = within(data_B, rm(rows_kept)),
  seed = SEED + 2,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95
)

# Optionally save lightweight fitted objects (QS or RDS)
# qs::qsave(x = fit_A, file = "fit_A_fullcov.qs")
# qs::qsave(x = fit_B, file = "fit_B_fullcov.qs")

# ------------------------------------------------------------
# 13) Subject-level LOO (exact marginal, from subj_loglik)
# ------------------------------------------------------------
# A
ll_A_s <- fit_A$draws("subj_loglik", format = "matrix") # (draws_total x S)
n_chains_A <- fit_A$num_chains()
n_draws_A <- nrow(ll_A_s)
stopifnot(n_draws_A %% n_chains_A == 0)
n_iter_A <- n_draws_A / n_chains_A
chain_id_A <- rep(seq_len(n_chains_A), each = n_iter_A)
r_eff_A <- loo::relative_eff(exp(ll_A_s), chain_id = chain_id_A)

# B
ll_B_s <- fit_B$draws("subj_loglik", format = "matrix")
n_chains_B <- fit_B$num_chains()
n_draws_B <- nrow(ll_B_s)
stopifnot(n_draws_B %% n_chains_B == 0)
n_iter_B <- n_draws_B / n_chains_B
chain_id_B <- rep(seq_len(n_chains_B), each = n_iter_B)
r_eff_B <- loo::relative_eff(exp(ll_B_s), chain_id = chain_id_B)

# PSIS-LOO with moment matching (more robust)
loo_A_s <- loo::loo(ll_A_s, r_eff = r_eff_A, moment_match = TRUE)
loo_B_s <- loo::loo(ll_B_s, r_eff = r_eff_B, moment_match = TRUE)

print(loo::loo_compare(list(B = loo_B_s, A = loo_A_s)))
loo::pareto_k_table(loo_A_s)
loo::pareto_k_table(loo_B_s)

# ------------------------------------------------------------
# 14) Convergence diagnostics (concise)
# ------------------------------------------------------------
posterior::summarise_draws(
  fit_B$draws(c("beta", "tau", "alpha_sigma", "sigma_delta")),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail"
) %>% print(n = Inf)

# ------------------------------------------------------------
# 15) Posterior Predictive Checks (global + targeted)
#     Uses: beta (K), b (R x S transformed), sigma_s (S)
# ------------------------------------------------------------
# Helper: compute random-effects contribution given b (R x S), Z (N x R), sid (N)
re_contrib_one <- function(b_RxS, Z, sid) {
  # b_RxS: numeric matrix R x S
  B_obs <- b_RxS[, sid, drop = FALSE]         # R x N
  colSums(Z * t(B_obs))                        # length N
}

# Extract objects from current environment for model B
X <- data_B$X
Z <- data_B$Z
sid <- data_B$sid
y  <- data_B$y
N <- length(y)
K <- ncol(X)
R <- ncol(Z)
S <- max(sid)

# Posterior draws for PPC
dr_beta  <- posterior::as_draws_matrix(fit_B$draws(paste0("beta[", 1:K, "]")))
# 'b' is R x S; in draws it appears as b[ r , s ]
dr_b     <- posterior::as_draws_matrix(fit_B$draws("b"))
dr_sigS  <- posterior::as_draws_matrix(fit_B$draws("sigma_s"))

nd <- min(300, nrow(dr_beta), nrow(dr_b), nrow(dr_sigS))
set.seed(SEED + 100)
i_samp <- sample(seq_len(nd), nd)

beta_draws <- dr_beta[i_samp, , drop = FALSE]     # nd x K
# Rebuild b as a 3D array: nd x R x S
# Column names like "b[1,1]", "b[2,1]", ..., "b[R,S]"
b_coln <- colnames(dr_b)
get_idx <- function(nm) {
  # extract indices inside brackets
  m <- regmatches(nm, regexec("^b\\[(\\d+),(\\d+)\\]$", nm))[[1]]
  c(as.integer(m[2]), as.integer(m[3]))
}
idx_mat <- t(vapply(b_coln, get_idx, FUN.VALUE = integer(2)))
ord <- order(idx_mat[,1], idx_mat[,2])  # order by r then s
dr_b <- dr_b[, ord, drop = FALSE]
idx_mat <- idx_mat[ord, , drop = FALSE]

# Build array
b_arr <- array(NA_real_, dim = c(nd, R, S))
for (j in seq_len(ncol(dr_b))) {
  r <- idx_mat[j, 1]; s <- idx_mat[j, 2]
  b_arr[, r, s] <- dr_b[i_samp, j]
}
# sigma_s: nd x S
sigS_mat <- as.matrix(dr_sigS[i_samp, paste0("sigma_s[", 1:S, "]"), drop = FALSE])

# Fixed part and random-effects contribution
mu_fix_all <- beta_draws %*% t(X)  # nd x N
mu_all <- matrix(NA_real_, nd, N)
for (d in 1:nd) {
  mu_all[d, ] <- mu_fix_all[d, ] + re_contrib_one(b_arr[d, , ], Z, sid)
}

# Replicated datasets
yrep <- matrix(NA_real_, nd, N)
for (d in 1:nd) {
  sigma_obs <- sigS_mat[d, sid]  # length N
  yrep[d, ] <- rnorm(N, mu_all[d, ], sigma_obs)
}

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
sd_by_s_obs <- tapply(y, sid, sd)
sd_by_s_reps <- apply(yrep, 1, function(yr) tapply(yr, sid, sd)) # S x nd
bayesplot::ppc_dens_overlay(sd_by_s_obs, t(sd_by_s_reps)[1:50, ])

# --- PPC: Subject means ---
ybar_obs  <- tapply(y, sid, mean)
ybar_reps <- apply(yrep, 1, function(yr) tapply(yr, sid, mean))
bayesplot::ppc_dens_overlay(ybar_obs, t(ybar_reps)[1:50, ])

# --- PPC: SD by subject × period (filter groups with n >= 2) ---
per <- df_mod$is_pre_bin + 2L * df_mod$is_post_bin
per <- per[data_B$rows_kept] # align to data_B
grp <- interaction(sid, per, drop = TRUE)
n_by_grp <- tapply(y, grp, function(v) sum(is.finite(v)))
keep_lvls <- names(n_by_grp)[n_by_grp >= 2L]
f_keep <- factor(grp, levels = keep_lvls)
idx_keep <- !is.na(f_keep)

sd_by_sp_obs <- tapply(y[idx_keep], f_keep[idx_keep], sd)

G <- length(keep_lvls)
sd_by_sp_reps <- t(vapply(
  X = seq_len(nrow(yrep)),
  FUN = function(d) tapply(yrep[d, idx_keep], f_keep[idx_keep], sd),
  FUN_VALUE = numeric(G)
))
bayesplot::ppc_dens_overlay(sd_by_sp_obs, sd_by_sp_reps[1:50, ])

# ------------------------------------------------------------
# 16) Session info (for reproducibility)
# ------------------------------------------------------------
print(sessionInfo())
