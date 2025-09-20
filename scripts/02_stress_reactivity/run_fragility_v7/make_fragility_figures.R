# ============================================================
# Figures for psychological fragility (baseline / pre / post)
# Using raw data and model-based posterior predictions
# ============================================================
# This script expects the following objects to already be present in memory:
#   - df_mod : data.frame with one row per observation, including columns
#              y (standardized fragility score), subj (1..S), female,
#              is_pre, is_post, is_pre_bin, is_post_bin, and optional EMA columns.
#   - data_B : list passed to Stan for model B, including elements:
#              N, S, K, X (N x K), Z (N x R), subj (length-N), rows_kept (indices),
#              and any other pieces used by the fitted model.
#   - fit_B  : cmdstanr fit object for the final model
#
# If these objects do not exist, you can load them from RDS files or
# source your modeling script prior to running this file.
#
# Output: figures saved under ./figs/
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(bayesplot)
  library(posterior)
  library(conflicted)
  library(purrr)
})

conflicts_prefer(base::match)

# ------------------------------------------------------------
# Safety checks
# ------------------------------------------------------------
stopifnot(exists("df_mod"), exists("data_B"), exists("fit_B"))
stopifnot(is.list(data_B), is.data.frame(df_mod))

# Make output dir
fig_dir <- file.path(
  here::here("scripts", "02_stress_reactivity", "run_fragility_v7", "figs")
)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ------------------------------------------------------------
# Helper: reconstruct random-effect contributions for each draw
# ------------------------------------------------------------
# b_mat:      nd x (S*R) draws of b_subj[s,r] stacked with r major dimension
# Z:          N x R  random-effects design matrix
# subj:       length-N subject indices in 1..S
# S, R:       subject count and number of random effects
# b_colnames: optional; reorder columns to [r-major, s-minor] if needed
re_contrib_all <- function(b_mat, Z, subj, S, R, b_colnames = NULL) {
  stopifnot(is.matrix(b_mat), is.matrix(Z), length(subj) == nrow(Z))
  if (!is.null(b_colnames)) {
    rx <- regexec("^b_subj\\[(\\d+),(\\d+)\\]$", b_colnames)
    mt <- regmatches(b_colnames, rx)
    ok <- vapply(mt, function(x) length(x) == 3, TRUE)
    if (all(ok)) {
      s_idx <- as.integer(vapply(mt, `[[`, "", 2))
      r_idx <- as.integer(vapply(mt, `[[`, "", 3))
      ord <- order(r_idx, s_idx) # order by r, then s
      b_mat <- b_mat[, ord, drop = FALSE]
    }
  }
  nd <- nrow(b_mat)
  N <- nrow(Z)
  out <- matrix(NA_real_, nd, N)
  for (d in seq_len(nd)) {
    # Fill by columns: first all subjects for r=1, then r=2, ...
    Bsr <- matrix(
      as.numeric(b_mat[d, ]),
      nrow = S,
      ncol = ncol(Z),
      byrow = FALSE
    )
    B_obs <- Bsr[subj, , drop = FALSE] # N x R
    out[d, ] <- rowSums(Z * B_obs) # length N
  }
  out
}

# ------------------------------------------------------------
# Extract a manageable subset of posterior draws
# ------------------------------------------------------------
X <- data_B$X
Z <- data_B$Z
subj <- data_B$subj
N <- length(data_B$y)
S <- max(subj)
R <- ncol(Z)

# Select parameter draws
param_names <- c("beta", "b_subj", "sigma")
if ("sd_delta" %in% dimnames(fit_B$draws())[[3]]) {
  param_names <- c(param_names, "sd_delta")
}
if (any(grepl("^delta_raw\\[", dimnames(fit_B$draws())[[3]]))) {
  param_names <- c(param_names, "delta_raw")
}

dr <- posterior::as_draws_matrix(fit_B$draws(param_names))
set.seed(1)
nd <- min(400L, nrow(dr)) # # of draws to use for figures
i_samp <- sample(seq_len(nrow(dr)), nd)

K <- ncol(X)
beta_draws <- dr[i_samp, paste0("beta[", 1:K, "]"), drop = FALSE]
sigma_draws <- as.numeric(dr[i_samp, "sigma"])

# Subject-specific heteroscedasticity (optional, if present in model)
has_sd_delta <- "sd_delta" %in% colnames(dr)
has_delta_raw <- any(grepl("^delta_raw\\[", colnames(dr)))
if (has_sd_delta && has_delta_raw) {
  sdD_draws <- as.numeric(dr[i_samp, "sd_delta"])
  delta_cols <- grep("^delta_raw\\[", colnames(dr), value = TRUE)
  delta_mat <- dr[i_samp, delta_cols, drop = FALSE] # nd x S
} else {
  sdD_draws <- NULL
  delta_mat <- NULL
}

# Random effects
b_cols <- grep("^b_subj\\[", colnames(dr), value = TRUE)
b_mat <- dr[i_samp, b_cols, drop = FALSE] # nd x (S*R)

# ------------------------------------------------------------
# Build posterior expected means mu for each observation and draw
# ------------------------------------------------------------
mu_fix_all <- beta_draws %*% t(X) # nd x N
re_all <- re_contrib_all(b_mat, Z, subj, S, R, b_colnames = b_cols)
mu_all <- mu_fix_all + re_all # nd x N

# Observation-specific sigma (if heteroscedastic by subject), else constant
if (!is.null(sdD_draws) && !is.null(delta_mat)) {
  sigma_obs <- matrix(NA_real_, nd, N)
  for (d in seq_len(nd)) {
    sig_s <- sigma_draws[d] * exp(sdD_draws[d] * as.numeric(delta_mat[d, 1:S]))
    sigma_obs[d, ] <- sig_s[subj]
  }
} else {
  sigma_obs <- matrix(rep(sigma_draws, each = N), nrow = nd, ncol = N)
}

# Optionally: posterior predictive draws yrep (useful for raw-vs-ppc overlays)
yrep <- matrix(NA_real_, nd, N)
for (d in seq_len(nd)) yrep[d, ] <- rnorm(N, mu_all[d, ], sigma_obs[d, ])

# ------------------------------------------------------------
# Figure 1 — Fragility across periods (raw + model-based)
# ------------------------------------------------------------
# PERIOD vector aligned to data_B rows
per_vec <- with(df_mod, is_pre_bin + 2L * is_post_bin)
per_vec <- per_vec[data_B$rows_kept]
per_lab <- factor(
  per_vec,
  levels = c(0, 1, 2),
  labels = c("Baseline", "Pre-exam", "Post-exam")
)
female_vec <- df_mod$female[data_B$rows_kept]

df_raw <- tibble::tibble(
  y = as.numeric(data_B$y),
  period = per_lab,
  female = factor(female_vec, levels = c(0, 1), labels = c("Male", "Female"))
)

# Model-based marginal means by period: for each draw, average mu_all by period
period_idx <- split(seq_len(N), per_lab)
per_means <- lapply(period_idx, function(idx) {
  rowMeans(mu_all[, idx, drop = FALSE])
})
per_df <- tibble::tibble(
  period = rep(names(period_idx), each = nd),
  draw = rep(seq_len(nd), times = length(period_idx)),
  mean = c(per_means[[1]], per_means[[2]], per_means[[3]])
) %>%
  group_by(period) %>%
  summarise(
    est = median(mean),
    l95 = quantile(mean, 0.025),
    u95 = quantile(mean, 0.975),
    .groups = "drop"
  )

p1 <- ggplot(df_raw, aes(period, y)) +
  geom_violin(aes(fill = period), width = 0.9, alpha = 0.22, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.45) +
  geom_point(
    data = df_raw %>% group_by(period) %>% summarise(y = mean(y)),
    aes(period, y),
    size = 2.2
  ) +
  geom_pointrange(
    data = per_df,
    aes(period, est, ymin = l95, ymax = u95),
    size = 0.7,
    color = "black"
  ) +
  labs(
    title = "Psychological fragility across periods",
    subtitle = "Raw distributions with model-based posterior means (points) and 95% CrI (ranges)",
    x = NULL,
    y = "Fragility (z)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
ggsave(
  file.path(fig_dir, "fig1_fragility_periods.png"),
  p1,
  width = 7,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "fig1_fragility_periods.pdf"),
  p1,
  width = 7,
  height = 4.5
)

# ------------------------------------------------------------
# Figure 2 — Baseline PID-5 Negative Affect vs pre→post change
# ------------------------------------------------------------
# Compute subject-level raw change: mean(Post) - mean(Pre)
subj_ids <- data_B$subj
has_pre <- per_vec == 1L
has_post <- per_vec == 2L

mean_by_sp <- function(yv, s, p) tapply(yv[p], s[p], mean)
y_pre_by_s <- mean_by_sp(data_B$y, subj_ids, has_pre)
y_post_by_s <- mean_by_sp(data_B$y, subj_ids, has_post)

# Keep subjects with both periods observed
common_subj <- intersect(names(y_pre_by_s), names(y_post_by_s))
y_delta_raw <- y_post_by_s[common_subj] - y_pre_by_s[common_subj]
subj_common <- as.integer(common_subj)

# Baseline Negative Affect (z) per subject
# We recover it from df_mod: take first appearance per subject
naff_by_s <- df_mod %>%
  select(subj, z_naff_b) %>%
  distinct(subj, .keep_all = TRUE) %>%
  arrange(subj)

x_naff <- naff_by_s$z_naff_b[subj_common]

df_change_raw <- tibble::tibble(
  subj = subj_common,
  delta_raw = as.numeric(y_delta_raw),
  z_naff_b = x_naff
)

# Model-based subject-level change: for each draw, average mu by subject & period
mu_pre_by_s <- vapply(
  seq_len(nd),
  function(d) tapply(mu_all[d, has_pre], subj_ids[has_pre], mean),
  FUN.VALUE = numeric(S)
)
mu_post_by_s <- vapply(
  seq_len(nd),
  function(d) tapply(mu_all[d, has_post], subj_ids[has_post], mean),
  FUN.VALUE = numeric(S)
)
# Align to common subjects
mu_pre_c <- mu_pre_by_s[, subj_common, drop = FALSE]
mu_post_c <- mu_post_by_s[, subj_common, drop = FALSE]
mu_delta <- mu_post_c - mu_pre_c # nd x |common|

df_change_mod <- tibble::tibble(
  subj = subj_common,
  z_naff_b = x_naff,
  delta_med = apply(mu_delta, 2, median),
  delta_l95 = apply(mu_delta, 2, quantile, 0.025),
  delta_u95 = apply(mu_delta, 2, quantile, 0.975)
)

# Plot: raw scatter + model-based intervals
p2 <- ggplot() +
  geom_point(
    data = df_change_raw,
    aes(z_naff_b, delta_raw),
    alpha = 0.4,
    size = 1.7
  ) +
  geom_errorbar(
    data = df_change_mod,
    aes(z_naff_b, ymin = delta_l95, ymax = delta_u95),
    width = 0
  ) +
  geom_point(
    data = df_change_mod,
    aes(z_naff_b, delta_med),
    size = 1.8,
    shape = 21,
    fill = "white"
  ) +
  geom_smooth(
    data = df_change_raw,
    aes(z_naff_b, delta_raw),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.2
  ) +
  labs(
    title = "Baseline PID-5 Negative Affect predicts pre→post change",
    subtitle = "Points = raw subject-level change;\nhollow points & bars = model-based median and 95% CrI",
    x = "Baseline Negative Affect (z)",
    y = "Δ Fragility (Post − Pre)"
  ) +
  theme_minimal(base_size = 12)
ggsave(
  file.path(fig_dir, "fig2_naff_vs_change.png"),
  p2,
  width = 6.5,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "fig2_naff_vs_change.pdf"),
  p2,
  width = 6.5,
  height = 4.5
)

# ------------------------------------------------------------
# Figure 3 — EMA PID-5 index vs pre→post change (raw & model-based)
# ------------------------------------------------------------
# Build a subject-level EMA index during Pre-exam.
# We compute a *model-weighted* index: X_ema %*% E[beta_ema], averaged per subject.
# If EMA columns are absent, this section will skip gracefully.

# Identify EMA columns by name from df_mod (they were appended after preprocessing)
ema_cols <- setdiff(
  colnames(df_mod),
  c(
    "y",
    "subj",
    "z_naff_b",
    "z_det_b",
    "z_ant_b",
    "z_dis_b",
    "z_psy_b",
    "female",
    "is_pre",
    "is_post",
    "is_pre_bin",
    "is_post_bin"
  )
)
ema_cols <- ema_cols[grepl("^w_", ema_cols)]

fig3_skipped <- FALSE
if (length(ema_cols) == 0) {
  message("No EMA columns found; skipping Figure 3.")
  fig3_skipped <- TRUE
} else {
  # Posterior mean of beta coefficients for EMA columns
  # Note: data_B$X has no column names; we compute betas' means using df_mod alignment.
  beta_means <- colMeans(beta_draws) # K-length
  # Map EMA names to their position in the X design: use x_cols_B order known from pipeline
  # Here we rebuild the X column order as in the modeling script
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
  x_cols_B <- c(x_cols_A, ema_cols)
  # Positions of EMA columns in X:
  ema_pos <- match(ema_cols, x_cols_B)
  stopifnot(all(is.finite(ema_pos)))

  # Rows belonging to data_B after filtering
  rows_kept <- data_B$rows_kept
  pre_rows <- rows_kept[per_vec == 1L]
  subj_pre <- subj[per_vec == 1L]

  # Build a per-row EMA linear predictor using posterior mean betas
  # Take the corresponding EMA columns from df_mod (already standardized there)
  X_ema_pre <- as.matrix(df_mod[pre_rows, ema_cols, drop = TRUE])
  eta_ema_pre <- as.numeric(X_ema_pre %*% beta_means[ema_pos])

  # Aggregate to subject-level mean EMA index during Pre
  ema_index_pre_by_s <- tapply(eta_ema_pre, subj_pre, mean)

  # Align to subjects with both PRE and POST data
  ema_idx_common <- ema_index_pre_by_s[as.character(subj_common)]

  df_change_ema <- tibble::tibble(
    subj = subj_common,
    ema_pre_index = as.numeric(ema_idx_common),
    delta_raw = df_change_raw$delta_raw[match(subj_common, df_change_raw$subj)],
    delta_med = df_change_mod$delta_med[match(subj_common, df_change_mod$subj)],
    delta_l95 = df_change_mod$delta_l95[match(subj_common, df_change_mod$subj)],
    delta_u95 = df_change_mod$delta_u95[match(subj_common, df_change_mod$subj)]
  )

  p3 <- ggplot() +
    geom_point(
      data = df_change_ema,
      aes(ema_pre_index, delta_raw),
      alpha = 0.4,
      size = 1.7
    ) +
    geom_errorbar(
      data = df_change_ema,
      aes(ema_pre_index, ymin = delta_l95, ymax = delta_u95),
      width = 0
    ) +
    geom_point(
      data = df_change_ema,
      aes(ema_pre_index, delta_med),
      size = 1.8,
      shape = 21,
      fill = "white"
    ) +
    geom_smooth(
      data = df_change_ema,
      aes(ema_pre_index, delta_raw),
      method = "lm",
      se = TRUE,
      linewidth = 0.8,
      alpha = 0.2
    ) +
    labs(
      title = "EMA PID-5 index (Pre) vs pre→post change",
      subtitle = "Model-weighted EMA index during Pre;\nraw change (points) and model-based Δ with 95% CrI (bars)",
      x = "EMA risk index (Pre-exam; model-weighted)",
      y = "Δ Fragility (Post − Pre)"
    ) +
    theme_minimal(base_size = 12)
  ggsave(
    file.path(fig_dir, "fig3_ema_vs_change.png"),
    p3,
    width = 6.5,
    height = 4.5,
    dpi = 300
  )
  ggsave(
    file.path(fig_dir, "fig3_ema_vs_change.pdf"),
    p3,
    width = 6.5,
    height = 4.5
  )
}

# ------------------------------------------------------------
# Figure 4 — Gender effects across periods (model-based + raw)
# ------------------------------------------------------------

# 1) Periodo allineato alle righe usate dal modello
per <- df_mod$is_pre_bin + 2L * df_mod$is_post_bin # 0=Baseline,1=Pre,2=Post
per <- per[data_B$rows_kept]
stopifnot(length(per) == data_B$N)
per_lab <- factor(
  per,
  levels = c(0, 1, 2),
  labels = c("Baseline", "Pre-exam", "Post-exam")
)

# 2) Genere BINARIO per OSSERVAZIONE (robusto)
#    Preferisci 'female' (0/1 per soggetto) se esiste e NON è standardizzato.
get_female_obs <- function() {
  # a) usa 'female' 0/1 per soggetto se disponibile
  if (
    exists("female", inherits = TRUE) &&
      is.numeric(female) &&
      length(female) >= max(data_B$subj, na.rm = TRUE) &&
      all(na.omit(unique(female)) %in% c(0, 1))
  ) {
    return(as.integer(female[data_B$subj]))
  }
  # b) prova con 'female_vec' ma forzando a binario nel caso sia stato z-scalato
  if (
    exists("female_vec", inherits = TRUE) &&
      is.numeric(female_vec) &&
      length(female_vec) >= max(data_B$subj, na.rm = TRUE)
  ) {
    uv <- na.omit(sort(unique(female_vec)))
    if (!all(uv %in% c(0, 1))) {
      # sembra standardizzato: riconduci a 0/1 prendendo la mediana come soglia
      thr <- median(female_vec, na.rm = TRUE)
      fv_bin <- as.integer(female_vec > thr)
      return(as.integer(fv_bin[data_B$subj]))
    } else {
      return(as.integer(female_vec[data_B$subj]))
    }
  }
  # c) fallback dai dati grezzi: ricostruisci da 'd' (richiede 'sex_col' e 'keep')
  if (exists("d", inherits = TRUE) && exists("sex_col", inherits = TRUE)) {
    sx <- d[[sex_col]]
    female_by_row <- if (is.numeric(sx)) as.integer(sx != 0) else {
      lv <- tolower(as.character(sx))
      as.integer(
        lv %in% c("f", "female", "donna", "femmina", "woman", "1", "true")
      )
    }
    female_by_row <- female_by_row[keep] # passa a d2
    female_obs <- female_by_row[data_B$rows_kept] # riallinea a data_B
    return(as.integer(female_obs))
  }
  stop("Cannot reconstruct binary female indicator.")
}
female_obs <- get_female_obs()
stopifnot(
  length(female_obs) == data_B$N,
  all(na.omit(unique(female_obs)) %in% c(0, 1))
)
female_fac_obs <- factor(
  female_obs,
  levels = c(0, 1),
  labels = c("Male", "Female")
)

# 3) Indici delle osservazioni per ogni combinazione Gender × Period
gfac <- interaction(female_fac_obs, per_lab, drop = TRUE) # N-level factor
grp_idx <- split(seq_len(N), gfac)

# 4) Riassunto MODEL-BASED: media di mu_all per ciascun gruppo e draw
grp_df <- purrr::map_dfr(names(grp_idx), function(g) {
  idx <- grp_idx[[g]]
  tibble::tibble(
    grp = g,
    draw = seq_len(nd),
    mean = rowMeans(mu_all[, idx, drop = FALSE])
  )
})

# 5) Posterior summaries per gruppo
grp_est <- grp_df %>%
  dplyr::group_by(grp) %>%
  dplyr::summarise(
    est = median(mean),
    l95 = quantile(mean, 0.025),
    u95 = quantile(mean, 0.975),
    .groups = "drop"
  ) %>%
  tidyr::separate(grp, into = c("gender", "period"), sep = "\\.")

# 6) Riassunto RAW (overlay)
raw_df <- tibble::tibble(y = y, grp = gfac) %>%
  dplyr::group_by(grp) %>%
  dplyr::summarise(
    raw_mean = mean(y),
    raw_se = sd(y) / sqrt(dplyr::n()),
    .groups = "drop"
  ) %>%
  tidyr::separate(grp, into = c("gender", "period"), sep = "\\.")

# 7) Plot
pd <- ggplot2::position_dodge(width = 0.35)
p_gender_period <- ggplot2::ggplot(
  grp_est,
  ggplot2::aes(x = period, y = est, color = gender, group = gender)
) +
  ggplot2::geom_point(position = pd, size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = l95, ymax = u95),
    width = 0.12,
    position = pd
  ) +
  ggplot2::geom_point(
    data = raw_df,
    ggplot2::aes(x = period, y = raw_mean, color = gender),
    shape = 1,
    size = 2,
    position = pd,
    inherit.aes = FALSE
  ) +
  ggplot2::labs(
    title = "Gender effects across periods",
    x = NULL,
    y = "Fragility (posterior mean)",
    color = "Gender"
  ) +
  ggplot2::theme_minimal()

print(p_gender_period)

ggsave(
  file.path(fig_dir, "fig_gender_period.png"),
  p_gender_period,
  width = 6.5,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "fig_gender_period.pdf"),
  p_gender_period,
  width = 6.5,
  height = 4.5
)

# ------------------------------------------------------------
# Console notes
# ------------------------------------------------------------
msg <- c(
  sprintf("Saved: %s", file.path(fig_dir, "fig1_fragility_periods.{png,pdf}")),
  sprintf("Saved: %s", file.path(fig_dir, "fig2_naff_vs_change.{png,pdf}")),
  if (fig3_skipped) "Figure 3 skipped (no EMA columns found)." else
    sprintf("Saved: %s", file.path(fig_dir, "fig3_ema_vs_change.{png,pdf}")),
  sprintf("Saved: %s", file.path(fig_dir, "fig4_gender_periods.{png,pdf}"))
)
writeLines(msg)


# ============================================================
# 16) Added predictive value of EMA (Model B vs Model A)
#     — ΔELPD per subject + partial R2 + model weights
# ============================================================

# ------------------------------
# 16.1 ΔELPD per subject (PSIS-LOO)
# ------------------------------
# We use the *subject-level* log-lik matrices and the LOO objects you computed.
stopifnot(!is.null(loo_A_s$pointwise), !is.null(loo_B_s$pointwise))

elpd_A_i <- as.numeric(loo_A_s$pointwise[, "elpd_loo"])
elpd_B_i <- as.numeric(loo_B_s$pointwise[, "elpd_loo"])
stopifnot(length(elpd_A_i) == length(elpd_B_i))
S <- length(elpd_A_i)

delta_elpd_i <- elpd_B_i - elpd_A_i # per-subject improvement
delta_elpd_tot <- sum(delta_elpd_i) # total improvement
delta_elpd_mean <- mean(delta_elpd_i) # average per subject
# Convert to an interpretable multiplicative improvement in average predictive density:
# exp(ΔELPD/S) - 1 gives the *average percentage gain per subject*.
mult_gain <- exp(delta_elpd_mean) - 1

cat(sprintf(
  "\n[ΔELPD] Total ΔELPD = %.1f; mean ΔELPD/subject = %.3f; multiplicative density gain per subject = +%.1f%%\n",
  delta_elpd_tot,
  delta_elpd_mean,
  100 * mult_gain
))

# Plot distribution of ΔELPD per subject (violin + box + zero line)
df_delpd <- tibble::tibble(
  subject = seq_len(S),
  delta_elpd = delta_elpd_i
)

p_delpd <- ggplot(df_delpd, aes(x = "", y = delta_elpd)) +
  geom_violin(fill = "grey85", color = "grey30", width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  annotate(
    "text",
    x = 1.2,
    y = quantile(delta_elpd_i, 0.5),
    label = sprintf("median = %.3f", median(delta_elpd_i)),
    hjust = 0,
    size = 3.5
  ) +
  labs(
    title = "Subject-level out-of-sample improvement: ΔELPD (Model B − Model A)",
    y = "ΔELPD per subject",
    x = NULL,
    caption = sprintf(
      "Mean ΔELPD/subject = %.3f → +%.1f%% predictive density gain on avg.",
      delta_elpd_mean,
      100 * mult_gain
    )
  ) +
  theme_minimal()

ggsave(
  file.path(fig_dir, "fig_delta_elpd_per_subject.png"),
  p_delpd,
  width = 6.5,
  height = 4.5,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "fig_delta_elpd_per_subject.pdf"),
  p_delpd,
  width = 6.5,
  height = 4.5
)


# Optional: a caterpillar to see which subjects benefit most
df_delpd2 <- df_delpd %>% arrange(delta_elpd) %>% mutate(rank = row_number())
p_delpd_cat <- ggplot(df_delpd2, aes(rank, delta_elpd)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  labs(
    title = "ΔELPD by subject (ordered)",
    x = "Subject (ranked)",
    y = "ΔELPD"
  ) +
  theme_minimal()

ggsave(
  file.path(fig_dir, "fig_delta_elpd_caterpillar.png"),
  p_delpd_cat,
  width = 7,
  height = 3.6,
  dpi = 300
)

# ------------------------------
# 16.2 Partial Bayesian R^2 for the EMA block
# ------------------------------
# We compute R2 per draw for each model:
#   R2 = Var(mu) / (Var(mu) + mean(sigma_obs^2))
# where mu = X beta + Z b_subj, and sigma_obs can vary by subject.
# Then we summarize ΔR2 = R2_B - R2_A with median and credible interval.

# Helper: reconstruct mu_all (nd x N) and sigma_obs (nd x N) from a fitted model and its data list
# re_contrib_all <- function(b_mat, Z, subj, S, R, b_colnames = NULL) {
#   stopifnot(is.matrix(b_mat), is.matrix(Z), length(subj) == nrow(Z))
#   if (!is.null(b_colnames)) {
#     rx <- regexec("^b_subj\\[(\\d+),(\\d+)\\]$", b_colnames)
#     mt <- regmatches(b_colnames, rx)
#     ok <- vapply(mt, function(x) length(x) == 3, TRUE)
#     if (all(ok)) {
#       s_idx <- as.integer(vapply(mt, `[[`, "", 2))
#       r_idx <- as.integer(vapply(mt, `[[`, "", 3))
#       ord <- order(r_idx, s_idx) # order by r then s (matches matrix fill below)
#       b_mat <- b_mat[, ord, drop = FALSE]
#     }
#   }
#   nd <- nrow(b_mat); N <- nrow(Z)
#   out <- matrix(NA_real_, nd, N)
#   for (d in seq_len(nd)) {
#     Bsr <- matrix(as.numeric(b_mat[d, ]), nrow = S, ncol = R, byrow = FALSE)
#     B_obs <- Bsr[subj, , drop = FALSE]              # N x R
#     out[d, ] <- rowSums(Z * B_obs)                  # fast rowwise dot
#   }
#   out
# }
#
# get_mu_sigma_draws <- function(fit, data_list, nd = 300) {
#   dr <- posterior::as_draws_matrix(
#     fit$draws(c("beta","b_subj","sigma","sd_delta","delta_raw"))
#   )
#   nd <- min(nd, nrow(dr))
#   set.seed(1)
#   i_samp <- sample(seq_len(nrow(dr)), nd)
#
#   X <- data_list$X; Z <- data_list$Z
#   subj <- data_list$subj
#   N <- nrow(X); K <- ncol(X)
#   R <- ncol(Z); S <- max(subj)
#
#   beta_draws <- dr[i_samp, paste0("beta[", 1:K, "]"), drop = FALSE]
#   sigma_draws <- as.numeric(dr[i_samp, "sigma"])
#   has_sdD <- "sd_delta" %in% colnames(dr)
#   sdD_draws <- if (has_sdD) as.numeric(dr[i_samp, "sd_delta"]) else rep(0, nd)
#   delta_mat <- if (any(grepl("^delta_raw\\[", colnames(dr))))
#     dr[i_samp, grep("^delta_raw\\[", colnames(dr)), drop = FALSE] else
#       matrix(0, nd, S)
#
#   b_cols <- grep("^b_subj\\[", colnames(dr), value = TRUE)
#   b_mat <- dr[i_samp, b_cols, drop = FALSE]
#
#   mu_fix_all <- beta_draws %*% t(X)                                  # nd x N
#   re_all <- re_contrib_all(b_mat, Z, subj, S, R, b_colnames = b_cols) # nd x N
#   mu_all <- mu_fix_all + re_all
#
#   sigma_obs <- matrix(NA_real_, nd, N)
#   for (d in 1:nd) {
#     sigma_s <- sigma_draws[d] * exp(sdD_draws[d] * as.numeric(delta_mat[d, 1:S]))
#     sigma_obs[d, ] <- sigma_s[subj]
#   }
#   list(mu_all = mu_all, sigma_obs = sigma_obs)
# }
#
# # Build draws for A and B and compute R2 per draw
# A_ms <- get_mu_sigma_draws(fit_A, data_A, nd = 400)
# B_ms <- get_mu_sigma_draws(fit_B, data_B, nd = 400)
#
# r2_from_mu_sigma <- function(mu_all, sigma_obs) {
#   apply(mu_all, 1, function(mu_d) {
#     var_mu <- stats::var(mu_d)
#     var_err <- mean(sigma_obs[1, ]^2) # placeholder, replaced below
#   })
# }
#
# r2_of <- function(ms) {
#   nd <- nrow(ms$mu_all)
#   vapply(seq_len(nd), function(d) {
#     var_mu  <- stats::var(ms$mu_all[d, ])
#     var_err <- mean(ms$sigma_obs[d, ]^2)
#     var_mu / (var_mu + var_err)
#   }, numeric(1))
# }
#
# r2_A <- r2_of(A_ms)
# r2_B <- r2_of(B_ms)
# delta_r2 <- r2_B - r2_A
#
# cat(sprintf(
#   "[ΔR2] median ΔR2 = %.003f (95%% CrI [%.003f, %.003f]); Pr(ΔR2>0) = %.3f\n",
#   median(delta_r2), quantile(delta_r2, c(.025,.975))[1], quantile(delta_r2, c(.025,.975))[2],
#   mean(delta_r2 > 0)
# ))
#
# # Plot ΔR2 distribution
# df_r2 <- tibble::tibble(delta_r2 = delta_r2)
# p_dr2 <- ggplot(df_r2, aes(x = delta_r2)) +
#   geom_vline(xintercept = 0, linetype = 2, color = "red") +
#   geom_density(fill = "grey85", color = "grey30") +
#   geom_segment(aes(x = quantile(delta_r2, .025),
#                    xend = quantile(delta_r2, .975), y = 0, yend = 0),
#                linewidth = 1.2) +
#   annotate("text", x = median(delta_r2), y = 0, vjust = -1,
#            label = sprintf("median ΔR2 = %.3f", median(delta_r2))) +
#   labs(title = "Partial Bayesian R² for EMA block (Model B − Model A)",
#        x = expression(Delta*R^2), y = "Density") +
#   theme_minimal()
# ggsave("figs/fig_delta_R2.png", p_dr2, width = 6.6, height = 4, dpi = 300)
# ggsave("figs/fig_delta_R2.pdf",  p_dr2, width = 6.6, height = 4)

# ------------------------------
# 16.3 Model weights (stacking & pseudo-BMA+)
# ------------------------------
# ------------------------------
# 16.3 Model weights (stacking & pseudo-BMA+)
# ------------------------------
# Use the exported helper in 'loo'
stopifnot(all(c("loo_model_weights") %in% getNamespaceExports("loo")))

mw_stack <- loo::loo_model_weights(
  list(A = loo_A_s, B = loo_B_s),
  method = "stacking"
)
mw_pbmaBB <- loo::loo_model_weights(
  list(A = loo_A_s, B = loo_B_s),
  method = "pseudobma",
  BB = TRUE
)

print(mw_stack)
print(mw_pbmaBB)

df_w <- tibble::tibble(
  model = c("A: baseline PID-5", "B: baseline + EMA PID-5"),
  stacking = as.numeric(mw_stack),
  pseudoBMA_BB = as.numeric(mw_pbmaBB)
) |>
  tidyr::pivot_longer(-model, names_to = "method", values_to = "weight")

p_w <- ggplot(df_w, aes(x = method, y = weight, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Model weights (out-of-sample)",
    x = NULL,
    y = "Weight",
    fill = "Model"
  ) +
  theme_minimal()

ggsave(
  file.path(fig_dir, "fig_model_weights.png"),
  p_w,
  width = 6.6,
  height = 4,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "fig_model_weights.pdf"),
  p_w,
  width = 6.6,
  height = 4
)
# Interpretation of model weights: with stacking, if Model B’s weight ≈ 1 and
# Model A’s ≈0, the out-of-sample evidence overwhelmingly favors including EMA.

# ------------------------------
# 16.4 Subject-level ΔELPD vs sample size
# ------------------------------

# 1) Per-subject pointwise ELPD (each column in loo_*_s is a subject)
stopifnot(!is.null(loo_A_s$pointwise), !is.null(loo_B_s$pointwise))
elpd_A_subj <- as.numeric(loo_A_s$pointwise[, "elpd_loo"])
elpd_B_subj <- as.numeric(loo_B_s$pointwise[, "elpd_loo"])
stopifnot(length(elpd_A_subj) == length(elpd_B_subj))

# ΔELPD per subject
dELPD_subj <- elpd_B_subj - elpd_A_subj
S_loo <- length(dELPD_subj)

# 2) Sample size per subject (using the data used in the fit)
N_per_subj <- tabulate(data_B$subj, nbins = max(data_B$subj))
stopifnot(length(N_per_subj) == S_loo)

# 3) Period coverage (how many distinct periods each subject has among kept rows)
#    Rebuild 'per' aligned to data_B if it's not already around:
if (!exists("per") || length(per) != length(data_B$subj)) {
  per_tmp <- df_mod$is_pre_bin + 2L * df_mod$is_post_bin # 0=base,1=pre,2=post
  per <- per_tmp[data_B$rows_kept] # align to data_B
}
stopifnot(length(per) == length(data_B$subj))

coverage_counts <- vapply(
  X = seq_len(S_loo),
  FUN = function(s) dplyr::n_distinct(per[data_B$subj == s]),
  FUN.VALUE = integer(1)
)
coverage_fac <- factor(
  coverage_counts,
  levels = c(1, 2, 3),
  labels = c("1 period", "2 periods", "3 periods")
)

# 4) Assemble dataframe
df_delta <- tibble::tibble(
  subject = seq_len(S_loo),
  dELPD = dELPD_subj,
  N = N_per_subj,
  coverage = coverage_fac
)

# 5) Quick summary to print in console
cat(sprintf(
  "ΔELPD (B−A): median=%.1f (IQR %.1f to %.1f); Spearman ρ(N,ΔELPD)=%.2f\n",
  median(df_delta$dELPD),
  quantile(df_delta$dELPD, .25),
  quantile(df_delta$dELPD, .75),
  suppressWarnings(cor(df_delta$N, df_delta$dELPD, method = "spearman"))
))

# 6) Plot: ΔELPD vs N per subject, colored by period coverage
p_delta <- ggplot(df_delta, aes(x = N, y = dELPD, color = coverage)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  labs(
    title = expression(Delta * "ELPD per subject vs. number of observations"),
    subtitle = paste0(
      "Spearman ρ = ",
      sprintf(
        "%.2f",
        suppressWarnings(cor(df_delta$N, df_delta$dELPD, method = "spearman"))
      ),
      " (higher is better for model B)"
    ),
    x = "Observations per subject",
    y = expression(Delta * "ELPD (" * B - A * ")"),
    color = "Period coverage"
  ) +
  theme_minimal()

# Save the panel

ggsave(
  file.path(fig_dir, "fig_deltaELPD_vs_N.png"),
  p_delta,
  width = 6.6,
  height = 4.4,
  dpi = 300
)
ggsave(
  file.path(fig_dir, "fig_deltaELPD_vs_N.pdf"),
  p_delta,
  width = 6.6,
  height = 4.4
)


# ==== Forest of fixed effects ====

#' 1) Forest plot of standardized fixed effects (Model B)
#' Why: Shows which baseline PID-5 domains and which EMA dimensions carry
#' predictive signal, with uncertainty.
#' What to look for: 95% CrI away from zero → practically relevant effects;
#' add a ROPE band (e.g., ±0.05 SD of y).

# Rebuild the column names in the same order you built X_B
# (base set first, then any EMA-within columns in Wdf, if present)

base_order <- c(
  "z_naff_b",
  "z_det_b",
  "z_ant_b",
  "z_dis_b",
  "z_psy_b",
  "female",
  "is_pre",
  "is_post"
)

ema_order <- character(0)
if (exists("Wdf") && is.data.frame(Wdf) && ncol(Wdf) > 0) {
  ema_order <- names(Wdf) # these were appended after the base columns
} else {
  # fallback: detect EMA columns that start with "w_" in df_mod
  ema_order <- grep("^w_", names(df_mod), value = TRUE)
}

fe_names <- c(base_order, ema_order)
stopifnot(length(fe_names) == ncol(data_B$X))

# (Optional) store names for later use and reattach to X for convenience
data_B$fe_names <- fe_names
colnames(data_B$X) <- fe_names

stopifnot(!is.null(fe_names), length(fe_names) == K)

beta_draws <- posterior::as_draws_matrix(fit_B$draws(paste0("beta[", 1:K, "]")))

post_fe <- tibble::tibble(
  name = fe_names,
  est = apply(beta_draws, 2, median),
  l95 = apply(beta_draws, 2, quantile, .025),
  u95 = apply(beta_draws, 2, quantile, .975)
) %>%
  arrange(desc(abs(est))) %>%
  dplyr::mutate(name = factor(.data$name, levels = rev(.data$name)))

# optional: order by |est| so the largest effects are on top
post_fe <- post_fe %>%
  arrange(desc(abs(est))) %>%
  mutate(term = factor(name, levels = rev(name)))

ggplot(post_fe, aes(x = est, y = term)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_rect(
    aes(xmin = -0.05, xmax = 0.05, ymin = -Inf, ymax = Inf),
    fill = "grey95",
    inherit.aes = FALSE
  ) +
  geom_point() +
  geom_errorbarh(aes(xmin = l95, xmax = u95), height = 0) +
  labs(
    x = "Posterior median (95% CrI), standardized scale",
    y = NULL,
    title = "Model B: Fixed effects"
  ) +
  theme_minimal()

#' 2) Caterpillar plot of subject random effects (intercept & period slopes)
#' Why: Shows heterogeneity across individuals and shrinkage.
#' What to look for: Spread and ordering; whether random period slopes cluster
#' above/below zero (systematic pre/post shifts).

# ==== Caterpillar for random effects ====
# b_subj has S x R random effects (intercept + period slopes + optional EMA slopes)
b_cols <- grep(
  "^b_subj\\[",
  colnames(posterior::as_draws_matrix(fit_B$draws())),
  value = TRUE
)
B <- posterior::as_draws_matrix(fit_B$draws(b_cols))
S <- max(data_B$subj)
R <- ncol(data_B$Z)

# Reconstitute per-RE parameter blocks (one tibble per random effect column r)
make_caterpillar <- function(r) {
  # columns for RE r are every R-th starting at r (Stan uses [s,r] in column names)
  idx <- seq(r, length(b_cols), by = R)
  mat <- B[, idx, drop = FALSE] # draws x S (subjects in columns)
  tibble::tibble(
    subject = 1:S,
    est = apply(mat, 2, median),
    l95 = apply(mat, 2, quantile, .025),
    u95 = apply(mat, 2, quantile, .975),
    re = colnames(data_B$Z)[r]
  )
}
cat_df <- dplyr::bind_rows(lapply(seq_len(R), make_caterpillar))

ggplot(cat_df, aes(y = reorder(subject, est), x = est)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  geom_errorbarh(aes(xmin = l95, xmax = u95), height = 0) +
  facet_wrap(~re, scales = "free_x") +
  labs(
    x = "Random effect (median & 95% CrI)",
    y = "Subject",
    title = "Subject-level random effects (Model B)"
  ) +
  theme_minimal()

#' 3) Period means: raw vs model-based (with uncertainty)
#' Why: Answers “what happens at baseline, pre, post?” using both raw summaries
#' and model-implied means.
#' What to look for: Consistency between raw and model; narrower uncertainty for
#' model-based summaries is common.

# ==== Period means ====
per <- df_mod$is_pre_bin + 2L * df_mod$is_post_bin
per <- per[data_B$rows_kept]
per_lab <- factor(
  per,
  levels = c(0, 1, 2),
  labels = c("Baseline", "Pre-exam", "Post-exam")
)

# Raw
raw_df <- tibble::tibble(y = data_B$y, period = per_lab) %>%
  group_by(period) %>%
  summarise(
    mean = mean(y),
    l95 = mean(y) - 1.96 * sd(y) / sqrt(n()),
    u95 = mean(y) + 1.96 * sd(y) / sqrt(n())
  )

# Model-based (posterior for period means via mu_all already computed in your PPC section)
# If you don't have mu_all: reconstruct for nd draws as in your PPC block.
# Then:
mb_df <- tibble::tibble(
  period = rep(levels(per_lab), each = nrow(mu_all)),
  draw = rep(1:nrow(mu_all), times = length(levels(per_lab))),
  mean = c(
    rowMeans(mu_all[, per == 0, drop = FALSE]),
    rowMeans(mu_all[, per == 1, drop = FALSE]),
    rowMeans(mu_all[, per == 2, drop = FALSE])
  )
) %>%
  group_by(period) %>%
  summarise(
    est = median(mean),
    l95 = quantile(mean, .025),
    u95 = quantile(mean, .975)
  )

ggplot() +
  geom_pointrange(
    data = raw_df,
    aes(x = period, y = mean, ymin = l95, ymax = u95),
    position = position_nudge(x = -0.12),
    shape = 21,
    fill = "white"
  ) +
  geom_pointrange(
    data = mb_df,
    aes(x = period, y = est, ymin = l95, ymax = u95),
    position = position_nudge(x = +0.12),
    shape = 19
  ) +
  labs(
    x = NULL,
    y = "Negative affect (z)",
    title = "Period means: raw vs model-based (Model B)"
  ) +
  theme_minimal()

#' 5) ΔELPD vs baseline severity (or composite PID-5 at baseline)
#' Why: Shows for whom EMA brings the largest gains.
#' What to look for: Positive slope → EMA helps more for high-risk baseline
#' profiles.

# ==== ΔELPD vs baseline composite ====
elpd_A_subj <- as.numeric(loo_A_s$pointwise[, "elpd_loo"])
elpd_B_subj <- as.numeric(loo_B_s$pointwise[, "elpd_loo"])
dELPD <- elpd_B_subj - elpd_A_subj
S <- length(dELPD)

# Subject-level baseline composite (mean of z-scored baseline domains)
Xb_subj <- df_mod %>%
  select(subj, z_naff_b, z_det_b, z_ant_b, z_dis_b, z_psy_b) %>%
  distinct() %>%
  group_by(subj) %>%
  summarise(baseline_comp = mean(c_across(starts_with("z_")), na.rm = TRUE)) %>%
  arrange(subj)

stopifnot(nrow(Xb_subj) >= S)
bl <- Xb_subj$baseline_comp[1:S]

df <- tibble::tibble(subj = 1:S, dELPD = dELPD, baseline = bl)
ggplot(df, aes(baseline, dELPD)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = .7) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  labs(
    x = "Baseline PID-5 composite (z)",
    y = expression(Delta * ELPD ~ "(B − A)"),
    title = "EMA value vs baseline severity"
  ) +
  theme_minimal()

#' 6) Predictive calibration: observed vs posterior predictive mean
#' Why: Overall calibration check; slope ~ 1, intercept ~ 0 is ideal.
#' What to look for: Deviations at tails → discuss discretization artifacts,
#' residual heterogeneity, etc.

# ==== Calibration: observed vs posterior predictive mean ====
# Compute posterior predictive mean \hat{y}_i from mu_all (average over draws)
yhat <- colMeans(mu_all) # N-length
cal_df <- tibble::tibble(y = data_B$y, yhat = yhat)

ggplot(cal_df, aes(yhat, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    x = "Posterior predictive mean",
    y = "Observed",
    title = "Calibration: observed vs posterior predictive mean (Model B)"
  ) +
  theme_minimal()

fit_cal <- lm(y ~ yhat, data = cal_df)
summary(fit_cal)$coef

#' 7) Posterior distribution of within-subject SD (heteroscedasticity)
#' Why: Communicates the extent of subject-level dispersion and its variation
#' (your sigma_s structure).
#' What to look for: Spread of sigma_s shows heterogeneity; optionally stratify
#' by period coverage.

# ==== Distribution of subject-specific sigma_s ====
# sigma_s = sigma * exp(sd_delta * delta_raw[s])
dr <- posterior::as_draws_matrix(fit_B$draws(c(
  "sigma",
  "sd_delta",
  "delta_raw[1]",
  paste0("delta_raw[", 2:max(data_B$subj), "]")
)))
sigma <- dr[, "sigma"]
sd_delta <- dr[, "sd_delta"]
S <- max(data_B$subj)
delta_mat <- dr[, paste0("delta_raw[", 1:S, "]"), drop = FALSE]

# Posterior for sigma_s (take medians per subject to plot)
sigma_s_med <- apply(exp(log(sigma)) * exp(sd_delta)^0 * 0, 1, function(x) NA) # placeholder to avoid confusion
sigma_s_med <- apply(exp(log(sigma) + 0 * 0), 1, function(x) NA) # ignore; use loop below

sigma_s_draws <- array(NA_real_, dim = c(nrow(dr), S))
for (s in 1:S) {
  sigma_s_draws[, s] <- sigma * exp(sd_delta * delta_mat[, s])
}
sigmaS_df <- tibble::tibble(
  subject = rep(1:S, each = nrow(dr)),
  sigma_s = as.numeric(sigma_s_draws)
) %>%
  group_by(subject) %>%
  summarise(
    med = median(sigma_s),
    l95 = quantile(sigma_s, .025),
    u95 = quantile(sigma_s, .975)
  )

ggplot(sigmaS_df, aes(x = med)) +
  geom_histogram(bins = 30) +
  labs(
    x = expression(paste("Subject-specific ", sigma[s])),
    y = "Count",
    title = "Heterogeneity of within-subject SD (Model B)"
  ) +
  theme_minimal()


#' 8) Posterior predictive contrasts: EMA vs baseline contribution
#' Why: Converts coefficients into predictive effect sizes: expected change in
#' NA for a +1 SD increase in a predictor (holding others at typical values).
#' What to look for: Compare distributions of contrasts for baseline vs EMA
#' blocks.

# ==== Predictive contrasts for +1 SD shifts ====
# Choose a reference covariate vector at the median of the data
x_ref <- apply(data_B$X, 2, median)

# Identify columns: baseline block vs EMA block (adjust names as needed)
idx_base <- grep(
  "^z_(naff|det|ant|dis|psy)_b$|^female$|^is_pre$|^is_post$",
  colnames(data_B$X)
)
idx_ema <- setdiff(seq_len(ncol(data_B$X)), idx_base)

# draws for beta
beta_draws <- posterior::as_draws_matrix(fit_B$draws(paste0(
  "beta[",
  1:ncol(data_B$X),
  "]"
)))

# Function: contrast for a +1 SD increase in feature j
contrast <- function(j) {
  # fixed part change: beta_j * 1 (since X are standardized)
  beta_draws[, j]
}

# Summarise blocks
summ_block <- function(idxs, label) {
  # aggregate via average absolute predictive shift
  eff <- rowMeans(abs(beta_draws[, idxs, drop = FALSE]))
  tibble::tibble(
    block = label,
    est = median(eff),
    l95 = quantile(eff, .025),
    u95 = quantile(eff, .975)
  )
}
rbind(
  summ_block(idx_base, "Baseline block"),
  summ_block(idx_ema, "EMA block")
) %>%
  ggplot(aes(x = block, y = est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = l95, ymax = u95), width = .1) +
  labs(
    y = "Average |predictive contrast| for +1 SD",
    x = NULL,
    title = "Predictive contribution: Baseline vs EMA blocks (Model B)"
  ) +
  theme_minimal()
