# 04_reactivity_sampling_epred.R
# Full NUTS sampling and draws-based extraction of per-subject exam reactivity
# using posterior_linpred() with random effects included.
#
# Model:
#   outcome ~ 0 + exam_period + (0 + exam_period | user_id)
#
# Outputs (per outcome):
#   data/derived/fit_reactivity_<outcome>.rds
#   data/derived/reactivity_indices_<outcome>.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(conflicted)
})

conflict_prefer("filter", "dplyr") # Always use dplyr::filter
conflict_prefer("select", "dplyr") # Always use dplyr::select


# ---------------- Configuration ----------------
min_responses_per_user <- 10
outcomes <- c("ucs_neg", "cs_pos", "neg_aff_ema") # edit as needed
iter <- 4000
chains <- 4
cores <- 4
seed <- 123
adapt_delta <- 0.95
ci_lower <- 0.055 # ~89% CI
ci_upper <- 0.945

# ---------------- Helpers ----------------
coverage_check <- function(df, outcome) {
  stopifnot(outcome %in% names(df))
  x <- df[[outcome]]
  if (is.matrix(x)) x <- as.vector(x)
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  if (length(x) != nrow(df)) stop("Length mismatch for outcome ", outcome)
  tab <- table(df$exam_period, is.na(x))
  message("\nCoverage for ", outcome, " (FALSE = observed, TRUE = NA):\n")
  print(tab)
  invisible(tab)
}

prepare_df <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- readRDS(path)
  stopifnot(all(c("user_id", "exam_period") %in% names(df)))
  if (!is.factor(df$exam_period)) df$exam_period <- factor(df$exam_period)
  df$exam_period <- forcats::fct_relevel(
    df$exam_period,
    "baseline",
    "pre_exam",
    "post_exam"
  )
  valid_ids <- df %>%
    count(user_id, name = "n") %>%
    dplyr::filter(n >= min_responses_per_user) %>%
    pull(user_id)
  df %>% dplyr::filter(user_id %in% valid_ids)
}

fit_reactivity_sampling <- function(dat, outcome) {
  message("Fitting (sampling): ", outcome)
  fml <- as.formula(paste0(
    outcome,
    " ~ 0 + exam_period + (0 + exam_period | user_id)"
  ))
  pri <- c(
    prior(student_t(3, 0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  )
  dat <- dat %>%
    mutate(
      exam_period = forcats::fct_relevel(
        exam_period,
        "baseline",
        "pre_exam",
        "post_exam"
      )
    ) %>%
    dplyr::select(user_id, exam_period, all_of(outcome)) %>%
    tidyr::drop_na(exam_period, all_of(outcome))

  brm(
    formula = fml,
    data = dat,
    family = student(),
    prior = pri,
    chains = chains,
    cores = cores,
    iter = iter,
    seed = seed,
    backend = "cmdstanr",
    algorithm = "sampling",
    control = list(adapt_delta = adapt_delta),
    refresh = 0
  )
}

# Extract subject-specific deltas from posterior linpred including random effects
extract_indices_epred <- function(fit, ci_lo = ci_lower, ci_hi = ci_upper) {
  # Use the user_id levels exactly as in the fit data
  ulevels <- levels(fit$data$user_id)
  if (is.null(ulevels)) ulevels <- levels(factor(fit$data$user_id))

  ep_levels <- c("baseline", "pre_exam", "post_exam")
  newdata <- expand.grid(
    user_id = ulevels,
    exam_period = factor(ep_levels, levels = ep_levels, ordered = FALSE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Draws of linear predictor including random effects
  # re_formula = NULL includes group-level effects
  eta <- posterior_linpred(
    fit,
    newdata = newdata,
    re_formula = NULL,
    transform = FALSE,
    allow_new_levels = FALSE
  )
  # eta: draws x (length(ulevels)*3)

  # Helper to get columns for a given (user, level)
  col_index <- function(user_idx, lvl_idx) {
    # columns are ordered by expand.grid: user major, then level varying slow/fast?
    # expand.grid in R varies the first factor slowest, last fastest
    # So order is: user1-baseline, user1-pre, user1-post, user2-baseline, ...
    (user_idx - 1) * length(ep_levels) + lvl_idx
  }

  res_list <- vector("list", length(ulevels))
  q <- function(x, p) unname(stats::quantile(x, probs = p, na.rm = TRUE))

  for (i in seq_along(ulevels)) {
    col_base <- col_index(i, 1)
    col_pre <- col_index(i, 2)
    col_post <- col_index(i, 3)

    d_pre <- eta[, col_pre] - eta[, col_base]
    d_post <- eta[, col_post] - eta[, col_base]

    magnitude <- sqrt(d_pre^2 + d_post^2)
    asym <- d_post - d_pre

    res_list[[i]] <- tibble::tibble(
      user_id = ulevels[i],
      react_pre = median(d_pre, na.rm = TRUE),
      react_pre_ll = q(d_pre, ci_lo),
      react_pre_ul = q(d_pre, ci_hi),
      react_pre_prob_gt0 = mean(d_pre > 0, na.rm = TRUE),

      react_post = median(d_post, na.rm = TRUE),
      react_post_ll = q(d_post, ci_lo),
      react_post_ul = q(d_post, ci_hi),
      react_post_prob_gt0 = mean(d_post > 0, na.rm = TRUE),

      magnitude = median(magnitude, na.rm = TRUE),
      magnitude_ll = q(magnitude, ci_lo),
      magnitude_ul = q(magnitude, ci_hi),

      asymmetry = median(asym, na.rm = TRUE),
      asymmetry_ll = q(asym, ci_lo),
      asymmetry_ul = q(asym, ci_hi)
    )
  }

  out <- dplyr::bind_rows(res_list)
  out[order(out$user_id), , drop = FALSE]
}

# ---------------- Main ----------------
main <- function() {
  in_path <- here("data", "processed", "df_exam_tagged.rds")
  df <- prepare_df(in_path)

  outdir <- here("data", "derived")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  for (y in outcomes) {
    if (!y %in% names(df)) {
      message("Skipping outcome '", y, "' (not found).")
      next
    }
    coverage_check(df, y)
    dat_y <- df %>%
      dplyr::select(user_id, exam_period, all_of(y)) %>%
      tidyr::drop_na()
    if (nrow(dat_y) == 0) {
      message("No rows for '", y, "'. Skipping.")
      next
    }

    fit <- fit_reactivity_sampling(dat_y, y)
    saveRDS(fit, file = file.path(outdir, paste0("fit_reactivity_", y, ".rds")))

    idx <- extract_indices_epred(fit, ci_lower, ci_upper)
    readr::write_csv(
      idx,
      file.path(outdir, paste0("reactivity_indices_", y, ".csv"))
    )
    message(
      "Saved: ",
      file.path(outdir, paste0("reactivity_indices_", y, ".csv"))
    )
  }
  message("All done.")
}

# Run if sourced directly
if (sys.nframe() == 0) main()
