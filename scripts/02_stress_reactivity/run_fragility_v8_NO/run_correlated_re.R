# PID-5 EMA — Full-cov RE + Student-t residuals
# - Within/Between centering for EMA covariates
# - Automatic pruning of random slopes by simple pre-fit variance screens
# - Subject-grouped PSIS-LOO using summed observation log-lik (subj_loglik_sum)

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
  library(rlang)
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("rhat", "posterior")

SEED <- 20250922
set.seed(SEED)

# ---------- Utilities ----------
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

# ---------- Paths ----------
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
stan_file <- here::here(
  "scripts",
  "02_stress_reactivity",
  "run_fragility_v8",
  "subject_loo_student_t_fullcov_rslope_sigmaS.stan"
)

# ---------- 1) Load & normalize ----------
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

# ---------- 2) Subject, period ----------
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

# ---------- 3) Sex (subject level) ----------
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
  warning("Sex column not found: set female=1 (placeholder).")
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

# ---------- 4) Outcome: Negative Affect (z) ----------
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

# ---------- 5) Between predictors (subject-level z) ----------
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

# Expand to observations
subj_idx <- as.integer(d2$.__subj__)
Xb_exp <- Xb[subj_idx, , drop = FALSE]
female_exp <- female[subj_idx]

# ---------- 6) EMA WITHIN/BETWEEN centering ----------
ema_candidates <- list(
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

EMA_df <- NULL
if (any(!sapply(ema_candidates, is.null))) {
  # 1) pull EMA columns + subject id
  ema_names <- na.omit(unlist(ema_candidates))
  tmp <- d2 %>%
    select(.__subj__, !!!rlang::syms(ema_names)) %>%
    rename(sid_subj = .__subj__)

  # 2) standardize raw EMA columns first (stabilizes means)
  cols_to_scale <- setdiff(names(tmp), "sid_subj")
  tmp_scaled <- tmp %>%
    mutate(across(all_of(cols_to_scale), ~ z_(as.numeric(.))))

  # 3) within/between split (explicit column list)
  EMA_df <- tmp_scaled %>%
    group_by(sid_subj) %>%
    mutate(
      across(
        all_of(cols_to_scale),
        list(
          mean = ~ mean(., na.rm = TRUE),
          within = ~ . - mean(., na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    ungroup()

  # 4) Z-score the *between* terms that actually exist
  between_cols <- paste0(cols_to_scale, "_mean")
  existing_between <- intersect(between_cols, names(EMA_df))
  EMA_df[existing_between] <- lapply(EMA_df[existing_between], z_)

  # (Optional) quick sanity check
  # print(setdiff(between_cols, names(EMA_df)))  # should be character(0)
}


# ---------- 7) Period dummies ----------
period <- as.integer(d2$.__per__)
is_pre <- as.integer(period == 2L)
is_post <- as.integer(period == 3L)
is_pre_bin <- is_pre
is_post_bin <- is_post

# ---------- 8) Assemble modeling frame ----------
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

# Append EMA within/between
if (!is.null(EMA_df)) {
  within_cols <- grep("_within$", names(EMA_df), value = TRUE)
  mean_cols <- grep("_mean$", names(EMA_df), value = TRUE)
  df_mod <- bind_cols(df_mod, EMA_df[, c(within_cols, mean_cols), drop = FALSE])
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

# ---------- 9) Build Stan data with automatic random-slope pruning ----------
# Heuristic: keep RE for predictors with enough within-subject variation.
# Rule: median(var_by_subj) >= 0.02 AND at least 70% subjects have variance > 0.
variance_screen <- function(x, sid) {
  v_by_s <- tapply(x, sid, stats::var, na.rm = TRUE)
  med_v <- median(v_by_s[is.finite(v_by_s)], na.rm = TRUE)
  prop_nz <- mean((v_by_s %>% replace_na(0)) > 0)
  keep <- is.finite(med_v) && med_v >= 0.02 && prop_nz >= 0.70
  list(keep = keep, med_v = med_v, prop_nz = prop_nz)
}

build_designs <- function(df) {
  # Fixed effects: base + period + EMA_between + EMA_within (all as fixed)
  X_cols_base <- c(
    "z_naff_b",
    "z_det_b",
    "z_ant_b",
    "z_dis_b",
    "z_psy_b",
    "female",
    "is_pre",
    "is_post"
  )
  ema_between <- grep("_mean$", names(df), value = TRUE)
  ema_within <- grep("_within$", names(df), value = TRUE)

  X_cols <- c(X_cols_base, ema_between, ema_within)
  X <- as.matrix(df[, X_cols, drop = FALSE])

  # Random effects: Intercept + period + EMA_within (subset by variance screen)
  sid <- as.integer(df$subj)

  cand_re <- c("is_pre", "is_post", ema_within)
  kept <- c("(Intercept)")

  variance_screen <- function(x, sid) {
    v_by_s <- tapply(x, sid, stats::var, na.rm = TRUE)
    med_v <- median(v_by_s[is.finite(v_by_s)], na.rm = TRUE)
    prop_nz <- mean((v_by_s %>% tidyr::replace_na(0)) > 0)
    keep <- is.finite(med_v) && med_v >= 0.02 && prop_nz >= 0.70
    list(keep = keep, med_v = med_v, prop_nz = prop_nz)
  }

  # decide which candidates to keep
  for (nm in cand_re) {
    sc <- variance_screen(df[[nm]], sid)
    if (isTRUE(sc$keep)) kept <- c(kept, nm)
  }

  # Build Z explicitly, adding EMA-within columns from df if they were kept
  Zdf <- data.frame(`(Intercept)` = rep(1, nrow(df)))
  if ("is_pre" %in% kept) Zdf$is_pre <- df$is_pre
  if ("is_post" %in% kept) Zdf$is_post <- df$is_post
  # add any kept EMA-within columns
  kept_within <- intersect(ema_within, kept)
  for (nm in kept_within) Zdf[[nm]] <- df[[nm]]

  Z <- as.matrix(Zdf)

  list(X = X, Z = Z, X_cols = X_cols, Z_cols = names(Zdf), sid = sid)
}


design_A <- build_designs(df_mod)
design_B <- design_A # (A/B identical here; adjust if you want an alternate spec)

data_list <- function(design, y) {
  list(
    N = nrow(design$X),
    S = max(design$sid),
    K = ncol(design$X),
    R = ncol(design$Z),
    X = unname(design$X),
    Z = unname(design$Z),
    sid = as.integer(design$sid),
    y = as.numeric(y)
  )
}
data_A <- data_list(design_A, df_mod$y)
data_B <- data_list(design_B, df_mod$y)

cat(sprintf(
  "Design kept RE columns: %s\n",
  paste(design_B$Z_cols, collapse = ", ")
))

# ---------- 10) Compile & sample ----------
stopifnot(file.exists(stan_file))
mod <- cmdstan_model(stan_file, force_recompile = TRUE)

fit_A <- mod$sample(
  data = data_A,
  seed = SEED + 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95
)
fit_B <- mod$sample(
  data = data_B,
  seed = SEED + 2,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95
)

# ---------- 11) Subject-grouped PSIS-LOO ----------
# We use subj_loglik_sum produced by Stan (sum of per-observation t log-liks by subject).
ll_A_s <- fit_A$draws("subj_loglik_sum", format = "matrix")
ll_B_s <- fit_B$draws("subj_loglik_sum", format = "matrix")

# r_eff by chain
n_chains_A <- fit_A$num_chains()
n_iter_A <- nrow(ll_A_s) / n_chains_A
chain_id_A <- rep(seq_len(n_chains_A), each = n_iter_A)
r_eff_A <- loo::relative_eff(exp(ll_A_s), chain_id = chain_id_A)

n_chains_B <- fit_B$num_chains()
n_iter_B <- nrow(ll_B_s) / n_chains_B
chain_id_B <- rep(seq_len(n_chains_B), each = n_iter_B)
r_eff_B <- loo::relative_eff(exp(ll_B_s), chain_id = chain_id_B)

loo_A_s <- loo::loo(ll_A_s, r_eff = r_eff_A, moment_match = TRUE)
loo_B_s <- loo::loo(ll_B_s, r_eff = r_eff_B, moment_match = TRUE)

print(loo::loo_compare(list(B = loo_B_s, A = loo_A_s)))
loo::pareto_k_table(loo_A_s)
loo::pareto_k_table(loo_B_s)

# ---------- 12) Diagnostics ----------
posterior::summarise_draws(
  fit_B$draws(c("beta", "tau", "alpha_sigma", "sigma_delta", "nu")),
  "mean",
  "sd",
  "rhat",
  "ess_bulk",
  "ess_tail"
) %>%
  print(n = Inf)

# ---------- 13) Posterior Predictive Checks ----------
# Rebuild mu quickly for subset of draws
X <- data_B$X
Z <- data_B$Z
sid <- data_B$sid
y <- data_B$y
N <- length(y)
K <- ncol(X)
R <- ncol(Z)
S <- max(sid)

dr_beta <- posterior::as_draws_matrix(fit_B$draws(paste0("beta[", 1:K, "]")))
dr_b <- posterior::as_draws_matrix(fit_B$draws("b"))
dr_sigS <- posterior::as_draws_matrix(fit_B$draws("sigma_s"))
dr_nu <- as.numeric(posterior::as_draws_matrix(fit_B$draws("nu"))[, 1])

nd <- min(300, nrow(dr_beta), nrow(dr_b), nrow(dr_sigS))
set.seed(SEED + 100)
i_samp <- sample(seq_len(nd), nd)
beta_draws <- dr_beta[i_samp, , drop = FALSE]

# Rebuild b array: nd x R x S
b_coln <- colnames(dr_b)
get_idx <- function(nm) {
  m <- regmatches(nm, regexec("^b\\[(\\d+),(\\d+)\\]$", nm))[[1]]
  c(as.integer(m[2]), as.integer(m[3]))
}
idx_mat <- t(vapply(b_coln, get_idx, FUN.VALUE = integer(2)))
ord <- order(idx_mat[, 1], idx_mat[, 2])
dr_b <- dr_b[, ord, drop = FALSE]
idx_mat <- idx_mat[ord, , drop = FALSE]
b_arr <- array(NA_real_, dim = c(nd, R, S))
for (j in seq_len(ncol(dr_b))) {
  r <- idx_mat[j, 1]
  s <- idx_mat[j, 2]
  b_arr[, r, s] <- dr_b[i_samp, j]
}
sigS_mat <- as.matrix(dr_sigS[
  i_samp,
  paste0("sigma_s[", 1:S, "]"),
  drop = FALSE
])
nu_vec <- dr_nu[i_samp]

mu_fix_all <- beta_draws %*% t(X) # nd x N
mu_all <- matrix(NA_real_, nd, N)
for (d in 1:nd) {
  # random-effects contribution for draw d: b_arr[d, , sid] is R x N
  mu_all[d, ] <- mu_fix_all[d, ] + colSums(Z * t(b_arr[d, , sid, drop = FALSE]))
}

# Generate yrep from Student-t
yrep <- matrix(NA_real_, nd, N)
for (d in 1:nd) {
  scale_vec <- sigS_mat[d, sid]
  yrep[d, ] <- stats::rt(N, df = nu_vec[d]) * scale_vec + mu_all[d, ]
}

# PPC overlays
bayesplot::ppc_dens_overlay(y, yrep[1:50, ])
bayesplot::ppc_ecdf_overlay(y, yrep[1:50, ])

sd_by_s_obs <- tapply(y, sid, sd)
sd_by_s_reps <- apply(yrep, 1, function(yr) tapply(yr, sid, sd))
bayesplot::ppc_dens_overlay(sd_by_s_obs, t(sd_by_s_reps)[1:50, ])

ybar_obs <- tapply(y, sid, mean)
ybar_reps <- apply(yrep, 1, function(yr) tapply(yr, sid, mean))
bayesplot::ppc_dens_overlay(ybar_obs, t(ybar_reps)[1:50, ])

# Period × subject SD
per <- df_mod$is_pre_bin + 2L * df_mod$is_post_bin
grp <- interaction(sid, per, drop = TRUE)
n_by_grp <- tapply(y, grp, function(v) sum(is.finite(v)))
keep_lvls <- names(n_by_grp)[n_by_grp >= 2L]
f_keep <- factor(grp, levels = keep_lvls)
idx_keep <- !is.na(f_keep)
sd_by_sp_obs <- tapply(y[idx_keep], f_keep[idx_keep], sd)
G <- length(keep_lvls)
sd_by_sp_reps <- t(vapply(
  seq_len(nrow(yrep)),
  function(d) tapply(yrep[d, idx_keep], f_keep[idx_keep], sd),
  FUN.VALUE = numeric(G)
))
bayesplot::ppc_dens_overlay(sd_by_sp_obs, sd_by_sp_reps[1:50, ])

message(
  "Random effects kept (after pruning): ",
  paste(design_B$Z_cols, collapse = ", ")
)
print(sessionInfo())
