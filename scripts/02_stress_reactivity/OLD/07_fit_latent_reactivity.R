
# 07_fit_latent_reactivity.R (updated)
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(cmdstanr)
  library(posterior)
  library(conflicted)
  library(tidyr)
  library(readr)
})

conflict_prefer("filter", "dplyr") # Always use dplyr::filter
conflict_prefer("select", "dplyr") # Always use dplyr::select
conflicts_prefer(stats::sd)

# -------------------- Load data --------------------
df <- readRDS(here("data", "processed", "df_exam_tagged.rds"))
df <- df %>% filter(exam_period %in% c("baseline","pre_exam","post_exam"))
df$exam_period <- factor(df$exam_period, levels = c("baseline","pre_exam","post_exam"))

# Outcomes
y_vars <- c("ucs_neg", "cs_pos")
stopifnot(all(y_vars %in% names(df)))

# Subject-level predictors
static_cols <- c("domain_negative_affect","domain_detachment","domain_antagonism",
                 "domain_disinhibition","domain_psychoticism")
dyn_cols <- c("pid5_negative_affectivity","pid5_detachment","pid5_antagonism",
              "pid5_disinhibition","pid5_psychoticism")

# Summaries
summ_period <- function(dsub, cols) {
  dsub %>% summarise(across(all_of(cols),
                            list(mean = ~mean(.x, na.rm = TRUE),
                                 sd   = ~sd(.x,   na.rm = TRUE)),
                            .names = "{.col}_{.fn}"))
}

dyn_by_user <- df %>%
  filter(exam_period %in% c("baseline","pre_exam")) %>%
  group_by(user_id, exam_period) %>%
  group_modify(~summ_period(.x, dyn_cols)) %>%
  ungroup() %>%
  pivot_longer(-c(user_id, exam_period), names_to = "var", values_to = "value") %>%
  separate(var, into = c("basevar","stat"), sep = "_(?=[^_]+$)") %>%
  unite("var_ep", basevar, exam_period, stat, sep = "_") %>%
  pivot_wider(names_from = var_ep, values_from = value)

for (v in dyn_cols) {
  mn_pre  <- paste0(v, "_pre_exam_mean")
  mn_base <- paste0(v, "_baseline_mean")
  dname   <- paste0(v, "_delta_mean")
  if (mn_pre %in% names(dyn_by_user) && mn_base %in% names(dyn_by_user)) {
    dyn_by_user[[dname]] <- dyn_by_user[[mn_pre]] - dyn_by_user[[mn_base]]
  }
}

static_by_user <- df %>%
  select(user_id, all_of(static_cols)) %>%
  group_by(user_id) %>%
  summarise(across(everything(), ~dplyr::first(na.omit(.x))), .groups = "drop")

pred_tbl <- static_by_user %>% left_join(dyn_by_user, by = "user_id")

sel_cols <- c(static_cols,
              paste0(dyn_cols, "_baseline_mean"),
              paste0(dyn_cols, "_baseline_sd"),
              paste0(dyn_cols, "_delta_mean"))
sel_cols <- sel_cols[sel_cols %in% names(pred_tbl)]

pred_tbl_z <- pred_tbl %>%
  mutate(across(all_of(sel_cols), ~as.numeric(scale(.x)), .names = "z_{.col}"))
Z_names <- paste0("z_", sel_cols)

# -------------------- Subject index --------------------
subj_ids <- sort(intersect(unique(df$user_id), pred_tbl_z$user_id))
I <- length(subj_ids)
subj_index <- tibble(user_id = subj_ids, subj = seq_along(subj_ids))

period_map <- c("baseline"=1L, "pre_exam"=2L, "post_exam"=3L)

prep_outcome <- function(var) {
  dd <- df %>% filter(!is.na(.data[[var]])) %>%
    inner_join(subj_index, by = "user_id") %>%
    transmute(y = .data[[var]], subj = subj, period = period_map[as.character(exam_period)])
  list(y = dd$y, subj = dd$subj, period = dd$period, N = nrow(dd))
}

o1 <- prep_outcome(y_vars[1])
o2 <- prep_outcome(y_vars[2])

# Build W (I x Q). If Q==0, create a single zero column and set Q=1.
Z_names <- Z_names[Z_names %in% names(pred_tbl_z)]
W <- pred_tbl_z %>%
  right_join(subj_index, by = "user_id") %>%
  arrange(subj) %>%
  mutate(across(all_of(Z_names), ~replace_na(.x, 0))) %>%
  select(all_of(Z_names))

Q <- ncol(W)
if (Q == 0) {
  W <- matrix(0, nrow = I, ncol = 1)
  Q <- 1
} else {
  W <- as.matrix(W)
}

# -------------------- Compile & sample --------------------
# Adjust the path to where you store the .stan file:
stan_path <- 
  here("script", "02_reactivity_index", "reactivity_latent_eta.stan")  # or here("script","02_reactivity_index","reactivity_latent_eta.stan")
mod <- cmdstan_model(stan_path)

dat <- list(
  I = I, P = 3L, K = 2L,
  N = c(o1$N, o2$N),
  y1 = as.vector(o1$y),
  subj1 = as.integer(o1$subj),
  period1 = as.integer(o1$period),
  y2 = as.vector(o2$y),
  subj2 = as.integer(o2$subj),
  period2 = as.integer(o2$period),
  Q = Q,
  W = W
)

fit <- mod$sample(
  data = dat, seed = 123,
  chains = 4, parallel_chains = 4, iter_warmup = 2000, iter_sampling = 2000,
  refresh = 200
)

saveRDS(fit, here("data","derived","latent_reactivity_fit.rds"))

# -------------------- Posterior summaries --------------------
draws <- as_draws_df(fit$draws())

draw_names <- names(fit$draws(format = "draws_df"))

# helper: summarise a set of variables by name with default q5/q95;
# if those quantiles aren't present, compute 89% (q5.5/q94.5) manually.
sdf <- function(fit, vars) {
  dm <- fit$draws(variables = vars, format = "draws_matrix")
  sum <- posterior::summarise_draws(dm)  # default includes q5 and q95
  if (all(c("q5","q95") %in% names(sum))) {
    sum <- sum %>% transmute(variable, mean, sd, ll = q5, ul = q95)
  } else {
    # manual 89% interval
    ql <- apply(dm, 2, stats::quantile, probs = 0.055, na.rm = TRUE)
    qu <- apply(dm, 2, stats::quantile, probs = 0.945, na.rm = TRUE)
    sum$ll <- as.numeric(ql); sum$ul <- as.numeric(qu)
    sum <- sum %>% transmute(variable, mean, sd, ll, ul)
  }
  sum
}

# --- measurement ---
alpha   <- sdf(fit, c("alpha[1]","alpha[2]"))
sigma   <- sdf(fit, c("sigma[1]","sigma[2]"))
lambda  <- sdf(fit, c("lambda[1]","lambda[2]"))

meas <- tibble::tibble(
  param = c("alpha_ucs_neg","alpha_cs_pos",
            "lambda_ucs_neg","lambda_cs_pos",
            "sigma_ucs_neg","sigma_cs_pos")
) |>
  mutate(
    mean = c(alpha$mean[alpha$variable=="alpha[1]"],
             alpha$mean[alpha$variable=="alpha[2]"],
             lambda$mean[lambda$variable=="lambda[1]"],
             lambda$mean[lambda$variable=="lambda[2]"],
             sigma$mean[sigma$variable=="sigma[1]"],
             sigma$mean[sigma$variable=="sigma[2]"]),
    sd   = c(alpha$sd[alpha$variable=="alpha[1]"],
             alpha$sd[alpha$variable=="alpha[2]"],
             lambda$sd[lambda$variable=="lambda[1]"],
             lambda$sd[lambda$variable=="lambda[2]"],
             sigma$sd[sigma$variable=="sigma[1]"],
             sigma$sd[sigma$variable=="sigma[2]"]),
    ll   = c(alpha$ll[alpha$variable=="alpha[1]"],
             alpha$ll[alpha$variable=="alpha[2]"],
             lambda$ll[lambda$variable=="lambda[1]"],
             lambda$ll[lambda$variable=="lambda[2]"],
             sigma$ll[sigma$variable=="sigma[1]"],
             sigma$ll[sigma$variable=="sigma[2]"]),
    ul   = c(alpha$ul[alpha$variable=="alpha[1]"],
             alpha$ul[alpha$variable=="alpha[2]"],
             lambda$ul[lambda$variable=="lambda[1]"],
             lambda$ul[lambda$variable=="lambda[2]"],
             sigma$ul[sigma$variable=="sigma[1]"],
             sigma$ul[sigma$variable=="sigma[2]"])
  )

readr::write_csv(meas, here("data","derived","latent_reactivity_meas.csv"))

# --- population ---
mu0 <- sdf(fit, c("mu0[1]","mu0[2]"))
tau <- sdf(fit, c("tau_eta[1]","tau_eta[2]"))

if ("Omega[1,2]" %in% draw_names) {
  rho <- sdf(fit, "Omega[1,2]")
  rho_mean <- rho$mean; rho_ll <- rho$ll; rho_ul <- rho$ul
} else {
  LO <- sdf(fit, c("L_Omega[1,1]","L_Omega[2,1]","L_Omega[2,2]"))
  LOm <- setNames(LO$mean, LO$variable)
  rho_mean <- LOm["L_Omega[2,1]"] / sqrt(LOm["L_Omega[2,1]"]^2 + LOm["L_Omega[2,2]"]^2)
  rho_ll <- rho_ul <- NA_real_
}

pop <- tibble::tibble(
  param = c("mu0_pre","mu0_post","tau_pre","tau_post","rho_pre_post"),
  mean  = c(mu0$mean[mu0$variable=="mu0[1]"],
            mu0$mean[mu0$variable=="mu0[2]"],
            tau$mean[tau$variable=="tau_eta[1]"],
            tau$mean[tau$variable=="tau_eta[2]"],
            rho_mean),
  sd    = c(mu0$sd[mu0$variable=="mu0[1]"],
            mu0$sd[mu0$variable=="mu0[2]"],
            tau$sd[tau$variable=="tau_eta[1]"],
            tau$sd[tau$variable=="tau_eta[2]"],
            NA_real_),
  ll    = c(mu0$ll[mu0$variable=="mu0[1]"],
            mu0$ll[mu0$variable=="mu0[2]"],
            tau$ll[tau$variable=="tau_eta[1]"],
            tau$ll[tau$variable=="tau_eta[2]"],
            rho_ll),
  ul    = c(mu0$ul[mu0$variable=="mu0[1]"],
            mu0$ul[mu0$variable=="mu0[2]"],
            tau$ul[tau$variable=="tau_eta[1]"],
            tau$ul[tau$variable=="tau_eta[2]"],
            rho_ul)
)
readr::write_csv(pop, here("data","derived","latent_reactivity_pop.csv"))

# --- per-subject mu_eta (transformed parameters) ---
pre_names  <- grep("^mu_eta\\[[0-9]+,1\\]$", draw_names, value = TRUE)
post_names <- sub(",1\\]", ",2]", pre_names)

mu_pre  <- sdf(fit, pre_names)
mu_post <- sdf(fit, post_names)

# IDs from the earlier 'subj_ids' object if present; else numeric indices
ids <- if (exists("subj_ids")) {
  subj_ids
} else {
  paste0("S", seq_along(pre_names))
}

react <- tibble::tibble(
  user_id       = ids,
  eta_pre_mean  = mu_pre$mean,
  eta_pre_ll89  = mu_pre$ll,
  eta_pre_ul89  = mu_pre$ul,
  eta_post_mean = mu_post$mean,
  eta_post_ll89 = mu_post$ll,
  eta_post_ul89 = mu_post$ul
) %>%
  mutate(
    magnitude = sqrt(eta_pre_mean^2 + eta_post_mean^2),
    asymmetry = eta_post_mean - eta_pre_mean
  )

readr::write_csv(react, here("data","derived","latent_reactivity_indices.csv"))
message("Saved meas/pop/indices CSVs to data/derived/")
'''

patched = head + start_marker + patched_tail
Path("/mnt/data/07_fit_latent_reactivity.R").write_text(patched, encoding="utf-8")
print("Patched 07_fit_latent_reactivity.R")
