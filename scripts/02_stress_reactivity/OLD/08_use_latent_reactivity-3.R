# 08_use_latent_reactivity.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(brms)
  library(cmdstanr)
  library(loo)
  library(bayesplot)
  library(conflicted)
})

options(mc.cores = max(1, parallel::detectCores() - 1))

conflict_prefer("filter", "dplyr") # Always use dplyr::filter
conflict_prefer("select", "dplyr") # Always use dplyr::select
conflicts_prefer(stats::sd)

# ---------------------------------------------------------------
# 0) Load files
# ---------------------------------------------------------------
df <- readRDS(here("data", "processed", "df_exam_tagged.rds"))
df <- df %>% filter(exam_period %in% c("baseline", "pre_exam", "post_exam"))

react <- readr::read_csv(
  here("data", "derived", "latent_reactivity_indices.csv"),
  show_col_types = FALSE
)
meas <- readr::read_csv(
  here("data", "derived", "latent_reactivity_meas.csv"),
  show_col_types = FALSE
)
pop <- readr::read_csv(
  here("data", "derived", "latent_reactivity_pop.csv"),
  show_col_types = FALSE
)

# Optional: read fit to extract B (predictor effects on eta)
fit_rds_path <- here("data", "derived", "latent_reactivity_fit.rds")
have_fit <- file.exists(fit_rds_path)
if (have_fit) {
  fit <- readRDS(fit_rds_path)
}

# ---------------------------------------------------------------
# 1) Quick measurement sanity checks
# ---------------------------------------------------------------
message("\n=== Measurement layer (alpha/lambda/sigma) ===")
print(meas)

if (any(meas$param == "lambda_cs_pos")) {
  lam_cs <- meas %>% filter(param == "lambda_cs_pos") %>% pull(mean)
  if (length(lam_cs) == 1 && lam_cs < 0) {
    message(
      "Note: lambda_cs_pos < 0 indicates cs_pos decreases when latent reactivity increases (expected)."
    )
  } else if (length(lam_cs) == 1 && lam_cs > 0) {
    message(
      "Note: lambda_cs_pos > 0 implies cs_pos increases with latent reactivity. ",
      "If your theoretical expectation is the opposite, consider allowing negative loadings ",
      "in Stan by removing '<lower=0>' on lambda_free and using a Normal(0,0.5) prior."
    )
  }
}

message("\n=== Population reactivity (mu0, tau, rho) ===")
print(pop)


# ---------------------------------------------------------------
# 2) Subject-level summaries for external criteria
#    We'll build DASS pre/post means and deltas as external outcomes
#    and also compute EMA affect (neg_aff_ema) as a check.
# ---------------------------------------------------------------
# If neg_aff_ema is not present, try a simple fallback from available items.
if (!"neg_aff_ema" %in% names(df)) {
  if (all(c("sad", "angry") %in% names(df))) {
    message("neg_aff_ema not found; creating fallback as rowMeans(sad, angry).")
    df <- df %>%
      mutate(neg_aff_ema = rowMeans(across(c("sad", "angry")), na.rm = TRUE))
  } else {
    message(
      "neg_aff_ema not found and fallback items not available; skipping neg_aff_ema summaries."
    )
  }
}

summ_cols <- c(
  "dass_sum",
  "dass_stress",
  "dass_depression",
  "dass_anxiety",
  "neg_aff_ema",
  "ucs_neg",
  "cs_pos"
)
avail_cols <- intersect(summ_cols, names(df))

subj_summ <- df %>%
  select(user_id, exam_period, any_of(avail_cols)) %>%
  group_by(user_id, exam_period) %>%
  summarise(
    across(all_of(avail_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = exam_period,
    values_from = all_of(avail_cols),
    names_glue = "{.value}_{exam_period}"
  )

# Compute deltas only if the required columns exist
if (all(c("dass_sum_pre_exam", "dass_sum_baseline") %in% names(subj_summ))) {
  subj_summ <- subj_summ %>%
    mutate(dass_sum_delta = dass_sum_pre_exam - dass_sum_baseline)
}
if (all(c("dass_sum_post_exam", "dass_sum_baseline") %in% names(subj_summ))) {
  subj_summ <- subj_summ %>%
    mutate(dass_sum_postch = dass_sum_post_exam - dass_sum_baseline)
}
if (
  all(c("neg_aff_ema_pre_exam", "neg_aff_ema_baseline") %in% names(subj_summ))
) {
  subj_summ <- subj_summ %>%
    mutate(neg_aff_delta = neg_aff_ema_pre_exam - neg_aff_ema_baseline)
}
if (
  all(c("neg_aff_ema_post_exam", "neg_aff_ema_baseline") %in% names(subj_summ))
) {
  subj_summ <- subj_summ %>%
    mutate(neg_aff_postch = neg_aff_ema_post_exam - neg_aff_ema_baseline)
}

# 2b) Bring static PID-5 (one per subject)
static_cols <- c(
  "domain_negative_affect",
  "domain_detachment",
  "domain_antagonism",
  "domain_disinhibition",
  "domain_psychoticism"
)

first_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) return(NA_real_)
  y[1]
}

static_tbl <- df %>%
  dplyr::select(user_id, dplyr::all_of(static_cols)) %>%
  dplyr::group_by(user_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(static_cols), first_non_na),
    .groups = "drop"
  )

# 2c) Merge latent indices + external summaries + static traits
subj <- react %>%
  dplyr::left_join(subj_summ, by = "user_id") %>%
  dplyr::left_join(static_tbl, by = "user_id") %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(static_cols),
      ~ as.numeric(scale(.x)),
      .names = "z_{.col}"
    ),
    z_magnitude = as.numeric(scale(magnitude)),
    z_asymmetry = as.numeric(scale(asymmetry))
  )
# ---------------------------------------------------------------
# 3) Visualize eta_pre vs eta_post and clusters
# ---------------------------------------------------------------
message("\nSaving figures to data/derived/figs ...")
figdir <- here("data", "derived", "figs")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

p_scatter <- ggplot(subj, aes(x = eta_pre_mean, y = eta_post_mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(alpha = 0.6) +
  coord_equal() +
  labs(
    x = expression(eta[pre] ~ "(mean)"),
    y = expression(eta[post] ~ "(mean)"),
    title = "Latent reactivity: subject means"
  ) +
  theme_minimal(base_size = 12)
ggsave(
  filename = file.path(figdir, "eta_pre_vs_post.png"),
  plot = p_scatter,
  width = 6,
  height = 5,
  dpi = 150
)

# Simple 4-quadrant typology
subj <- subj %>%
  mutate(
    profile = case_when(
      eta_pre_mean <= 0 & eta_post_mean <= 0 ~ "Low-reactive",
      eta_pre_mean > 0 & eta_post_mean <= 0 ~ "Reactive-then-recovers",
      eta_pre_mean > 0 & eta_post_mean > 0 ~ "Prolonged-reactive",
      eta_pre_mean <= 0 & eta_post_mean > 0 ~ "Delayed-reactive"
    )
  )

p_typ <- ggplot(subj, aes(eta_pre_mean, eta_post_mean, color = profile)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(alpha = 0.7) +
  coord_equal() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Idiographic reactivity profiles (quadrants)",
    x = expression(eta[pre]),
    y = expression(eta[post])
  )
ggsave(
  filename = file.path(figdir, "eta_quadrants.png"),
  plot = p_typ,
  width = 6,
  height = 5,
  dpi = 150
)

# ---------------------------------------------------------------
# 4) Incremental validity: do magnitude/asymmetry add beyond static PID-5?
#    Outcome: change in DASS_SUM pre-baseline (dass_sum_delta)
# ---------------------------------------------------------------
pri <- c(
  prior(student_t(3, 0, 1), class = "b"),
  prior(student_t(3, 0, 1), class = "Intercept"),
  prior(exponential(1), class = "sigma")
)

# z-scale predictors for comparability
subj2 <- subj %>%
  mutate(across(
    all_of(static_cols),
    ~ as.numeric(scale(.x)),
    .names = "z_{.col}"
  ))

# same N across models: complete cases on all variables used at least once
pred_static <- paste0("z_", static_cols)
req <- unique(c("dass_sum_delta", pred_static, c("z_magnitude", "z_asymmetry")))
dmod <- subj2 %>% select(user_id, any_of(req)) %>% drop_na()

f_static <- reformulate(pred_static, response = "dass_sum_delta")
f_react <- reformulate(
  c("z_magnitude", "z_asymmetry"),
  response = "dass_sum_delta"
)
f_combined <- reformulate(
  c(pred_static, "z_magnitude", "z_asymmetry"),
  response = "dass_sum_delta"
)

fit_static <- brm(
  f_static,
  data = dmod,
  family = student(),
  prior = pri,
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  seed = 2025,
  refresh = 0
)
fit_react <- brm(
  f_react,
  data = dmod,
  family = student(),
  prior = pri,
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  seed = 2025,
  refresh = 0
)
fit_comb <- brm(
  f_combined,
  data = dmod,
  family = student(),
  prior = pri,
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  seed = 2025,
  refresh = 0
)

loo_s <- loo(fit_static)
loo_r <- loo(fit_react)
loo_c <- loo(fit_comb)
cmp <- loo_compare(loo_s, loo_r, loo_c)
print(cmp)

# Bayes R2
br2_s <- bayes_R2(fit_static)
br2_r <- bayes_R2(fit_react)
br2_c <- bayes_R2(fit_comb)
print(br2_s)
print(br2_r)
print(br2_c)

# Save a concise summary
sink(here("data", "derived", "incremental_validity_dass.txt"))
cat("=== LOO comparison (elpd_loo higher is better) ===\n")
print(cmp)
cat("\n=== Bayes R2 ===\n")
print(br2_s)
print(br2_r)
print(br2_c)
cat("\n=== Coefficients (combined model) ===\n")
print(summary(fit_comb))
sink()

# ---------------------------------------------------------------
# 5) Which predictors in Stan drove eta_pre vs eta_post? (optional)
#    Extract B and save mean/intervals if fit is available.
# ---------------------------------------------------------------
if (have_fit) {
  draws_df <- fit$draws(format = "draws_df")
  # columns like B[1,1], B[2,1], ..., map them back to variable names used in W (see 07_* script)
  # We'll try to reconstruct from 'reactivity_latent_reactivity.R' logic: Z_ names from pred_tbl_z
  # but if not available, we can at least export an index.
  B_cols_pre <- grep("^B\\[[0-9]+,1\\]$", names(draws_df), value = TRUE)
  B_cols_post <- grep("^B\\[[0-9]+,2\\]$", names(draws_df), value = TRUE)
  if (length(B_cols_pre) == length(B_cols_post) && length(B_cols_pre) > 0) {
    # summarise
    summ_fun <- function(x) {
      c(
        mean = mean(x),
        sd = sd(x),
        ll = quantile(x, 0.055),
        ul = quantile(x, 0.945)
      )
    }
    S_pre <- sapply(draws_df[, B_cols_pre, drop = FALSE], summ_fun)
    S_post <- sapply(draws_df[, B_cols_post, drop = FALSE], summ_fun)
    B_idx <- seq_len(ncol(S_pre))
    out_B <- tibble(
      idx = B_idx,
      B_pre_mean = S_pre["mean", ],
      B_pre_ll89 = S_pre["ll", ],
      B_pre_ul89 = S_pre["ul", ],
      B_post_mean = S_post["mean", ],
      B_post_ll89 = S_post["ll", ],
      B_post_ul89 = S_post["ul", ]
    )
    write_csv(out_B, here("data", "derived", "latent_B_coeffs.csv"))
  }
}

message(
  "\nDone. See data/derived/ for outputs and data/derived/figs for plots."
)
