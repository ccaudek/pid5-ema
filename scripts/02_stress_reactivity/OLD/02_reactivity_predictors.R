# 05_reactivity_predictors.R
# From per-subject reactivity indices -> predictors analysis.
# Compares STATIC PID-5 (full) vs DYNAMIC PID-5 EMA in predicting exam reactivity.
#
# Inputs:
#   data/derived/reactivity_indices_ucs_neg.csv
#   data/derived/reactivity_indices_cs_pos.csv
#   data/processed/df_exam_tagged.rds
#
# Outputs:
#   data/derived/reactivity_subject_level.csv  (merged indices + features)
#   data/derived/model_summary_reactivity.txt  (ELPD, stacking weights, Bayes R2)
#   figures/ (optional quick plots)

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(conflicted)
})

options(mc.cores = parallel::detectCores())

conflict_prefer("filter", "dplyr") # Always use dplyr::filter
conflict_prefer("select", "dplyr") # Always use dplyr::select
conflicts_prefer(stats::sd)

# ----------------------- Load reactivity indices -----------------------
read_idx <- function(name) {
  path <- here("data", "derived", paste0("reactivity_indices_", name, ".csv"))
  if (!file.exists(path)) stop("Missing file: ", path)
  readr::read_csv(path, show_col_types = FALSE) %>%
    rename_with(~ paste0(., "_", name), -user_id)
}

idx_ucs <- read_idx("ucs_neg")
idx_cs <- read_idx("cs_pos")

# Merge by user_id
idx <- idx_ucs %>% full_join(idx_cs, by = "user_id")

# ----------------------- Composite indices -----------------------------
# Align direction: "worse reactivity" should be positive.
# For ucs_neg, an increase pre vs baseline is "worse": use +react_pre_ucs_neg
# For cs_pos, a decrease pre vs baseline is "worse": use -react_pre_cs_pos
# We standardize components before averaging.

z <- function(x) as.numeric(scale(x))

idx <- idx %>%
  mutate(
    comp_pre = rowMeans(
      cbind(z(react_pre_ucs_neg), z(-react_pre_cs_pos)),
      na.rm = TRUE
    ),
    comp_post = rowMeans(
      cbind(z(react_post_ucs_neg), z(-react_post_cs_pos)),
      na.rm = TRUE
    ),
    comp_mag = rowMeans(
      cbind(z(magnitude_ucs_neg), z(magnitude_cs_pos)),
      na.rm = TRUE
    ),
    comp_asym = rowMeans(
      cbind(z(asymmetry_ucs_neg), z(asymmetry_cs_pos)),
      na.rm = TRUE
    )
  )

# ----------------------- Visual quick looks ----------------------------
figdir <- here("figures")
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

p1 <- idx %>%
  ggplot(aes(x = comp_pre)) +
  geom_histogram(bins = 40) +
  labs(
    title = "Composite reactivity (pre vs baseline)",
    x = "comp_pre (z-avg)",
    y = "count"
  )
ggsave(
  filename = file.path(figdir, "hist_comp_pre.png"),
  plot = p1,
  width = 6,
  height = 4,
  dpi = 120
)

p2 <- idx %>%
  ggplot(aes(x = comp_post)) +
  geom_histogram(bins = 40) +
  labs(
    title = "Composite reactivity (post vs baseline)",
    x = "comp_post (z-avg)",
    y = "count"
  )
ggsave(
  filename = file.path(figdir, "hist_comp_post.png"),
  plot = p2,
  width = 6,
  height = 4,
  dpi = 120
)

# ----------------------- Aggregate predictors --------------------------
df <- readRDS(here("data", "processed", "df_exam_tagged.rds"))

# STATIC PID-5 (invarianti per riga/soggetto): columns domain_*
static_cols <- c(
  "domain_negative_affect",
  "domain_detachment",
  "domain_antagonism",
  "domain_disinhibition",
  "domain_psychoticism"
)
stopifnot(all(static_cols %in% names(df)))

# DYNAMIC PID-5 EMA domains (3-item sums per ping): pid5_negative_affectivity, ..., pid5_psychoticism
dyn_cols <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)
have_dyn <- dyn_cols %in% names(df)
if (!all(have_dyn)) {
  warning(
    "Missing dynamic PID-5 EMA columns: ",
    paste(dyn_cols[!have_dyn], collapse = ", "),
    "\nWill drop them from dynamic predictors."
  )
  dyn_cols <- dyn_cols[have_dyn]
}

# function to summarize over a period
summ_period <- function(dsub, cols) {
  dsub %>%
    summarise(across(
      all_of(cols),
      list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ))
}

# Build per-subject dynamic summaries in baseline and pre_exam
dyn_by_user <- df %>%
  filter(exam_period %in% c("baseline", "pre_exam")) %>%
  group_by(user_id, exam_period) %>%
  group_modify(~ summ_period(.x, dyn_cols)) %>%
  ungroup() %>%
  pivot_longer(
    -c(user_id, exam_period),
    names_to = "var",
    values_to = "value"
  ) %>%
  separate(var, into = c("basevar", "stat"), sep = "_(?=[^_]+$)") %>%
  unite("var_ep", basevar, exam_period, stat, sep = "_") %>%
  pivot_wider(names_from = var_ep, values_from = value)

# Compute deltas (pre - baseline) for means
for (v in dyn_cols) {
  mn_pre <- paste0(v, "_pre_exam_mean")
  mn_base <- paste0(v, "_baseline_mean")
  dname <- paste0(v, "_delta_mean")
  if (mn_pre %in% names(dyn_by_user) && mn_base %in% names(dyn_by_user)) {
    dyn_by_user[[dname]] <- dyn_by_user[[mn_pre]] - dyn_by_user[[mn_base]]
  }
}

# STATIC features per subject (take first non-missing row)
static_by_user <- df %>%
  select(user_id, all_of(static_cols)) %>%
  group_by(user_id) %>%
  summarise(across(everything(), ~ dplyr::first(na.omit(.x))), .groups = "drop")

# Merge all subject-level features
subj_tbl <- idx %>%
  left_join(static_by_user, by = "user_id") %>%
  left_join(dyn_by_user, by = "user_id")

# Save subject-level table
out_tbl <- here("data", "derived", "reactivity_subject_level.csv")
readr::write_csv(subj_tbl, out_tbl)

# ----------------------- Modeling: Static vs Dynamic -------------------
# Outcome = comp_pre (z-avg; higher = worse). We z-scale predictors for comparability.
# Build three models:
#  M_static:   comp_pre ~ z(static PID-5 domains)
#  M_dynamic:  comp_pre ~ z(dyn baseline mean) + z(dyn delta mean) + z(dyn baseline sd)
#  M_combined: comp_pre ~ both sets

zcols <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  out <- df
  for (c in cols) out[[paste0("z_", c)]] <- as.numeric(scale(out[[c]]))
  out
}

# identify dynamic summary columns
dyn_mean_base <- paste0(dyn_cols, "_baseline_mean")
dyn_mean_base <- dyn_mean_base[dyn_mean_base %in% names(subj_tbl)]
dyn_mean_delta <- paste0(dyn_cols, "_delta_mean")
dyn_mean_delta <- dyn_mean_delta[dyn_mean_delta %in% names(subj_tbl)]
dyn_sd_base <- paste0(dyn_cols, "_baseline_sd")
dyn_sd_base <- dyn_sd_base[dyn_sd_base %in% names(subj_tbl)]

# z-scale all candidate predictors
subj_z <- subj_tbl %>%
  zcols(static_cols) %>%
  zcols(dyn_mean_base) %>%
  zcols(dyn_mean_delta) %>%
  zcols(dyn_sd_base)


# --- ensure identical N across models (same rows) ---

# 1) Define predictor name sets (already z-scored)
pred_static <- paste0(
  "z_",
  c(
    "domain_negative_affect",
    "domain_detachment",
    "domain_antagonism",
    "domain_disinhibition",
    "domain_psychoticism"
  )
)

pred_dyn <- c(
  paste0("z_", paste0(dyn_cols, "_baseline_mean")),
  paste0("z_", paste0(dyn_cols, "_delta_mean")),
  paste0("z_", paste0(dyn_cols, "_baseline_sd"))
)
pred_dyn <- pred_dyn[pred_dyn %in% names(subj_z)] # keep only existing columns

# 2) Build complete-case data on the UNION of all predictors + outcome
req_all <- unique(c("user_id", "comp_pre", pred_static, pred_dyn))
dmod_all <- subj_z %>%
  dplyr::select(any_of(req_all)) %>%
  tidyr::drop_na()

# 3) Build RHS variable vectors actually present in dmod_all
X_static <- intersect(pred_static, names(dmod_all))
X_dynamic <- intersect(pred_dyn, names(dmod_all))
X_combined <- union(X_static, X_dynamic)

# 4) Build formulas programmatically (base R)
f_static <- reformulate(termlabels = X_static, response = "comp_pre") # comp_pre ~ 1 + ...
f_dynamic <- reformulate(termlabels = X_dynamic, response = "comp_pre")
f_combined <- reformulate(termlabels = X_combined, response = "comp_pre")

# (Optional) sanity check: same N for all fits
stopifnot(nrow(dmod_all) > 0)

# 5) Fit the three models on the SAME rows
fit_static <- brm(
  formula = f_static,
  data = dmod_all,
  family = student(),
  prior = pri,
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  seed = 123,
  refresh = 0,
  control = list(adapt_delta = 0.95)
)
pp_check(fit_static)

fit_dynamic <- brm(
  formula = f_dynamic,
  data = dmod_all,
  family = student(),
  prior = pri,
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  seed = 123,
  refresh = 0,
  control = list(adapt_delta = 0.95)
)
pp_check(fit_dynamic)


fit_combined <- brm(
  formula = f_combined,
  data = dmod_all,
  family = student(),
  prior = pri,
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  seed = 123,
  refresh = 0,
  control = list(adapt_delta = 0.95)
)
pp_check(fit_combined)


# 6) Now LOO is comparable (same N)
loo_s <- loo(fit_static)
loo_d <- loo(fit_dynamic)
loo_c <- loo(fit_combined)

cmp <- loo::loo_compare(loo_s, loo_d, loo_c)
print(cmp)

# Stacking weights (same model order as above)
sw <- loo::loo_model_weights(list(loo_s, loo_d, loo_c), method = "stacking")
names(sw) <- c("static", "dynamic", "combined")
print(sw)

# Bayes R2 on the same rows
br2_s <- bayes_R2(fit_static)
br2_d <- bayes_R2(fit_dynamic)
br2_c <- bayes_R2(fit_combined)
print(br2_s)
print(br2_d)
print(br2_c)


# Save summaries
sum_path <- here("data", "derived", "model_summary_reactivity.txt")
sink(sum_path)
cat("=== Model comparison (higher elpd_loo is better) ===\n")
print(cmp)
cat("\n=== Stacking weights ===\n")
print(sw)
cat("\n=== Bayes R2 ===\n")
print(br2_s)
print(br2_d)
print(br2_c)
cat("\n=== Coefficients (posterior means and 89% CI) ===\n")
print(summary(fit_static))
print(summary(fit_dynamic))
print(summary(fit_combined))
sink()

message("Wrote: ", out_tbl, " and ", sum_path)
