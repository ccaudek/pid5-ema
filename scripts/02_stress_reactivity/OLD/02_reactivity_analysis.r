# ============================================================
# Joint Measurement-Structural Model for EMA Negative Affect
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(conflicted)
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("sd", "stats")
conflict_prefer("var", "stats")
conflict_prefer("mad", "stats")


set.seed(20250912)

# ------------------------------------------------------------
# Data preparation for joint model
# ------------------------------------------------------------

# helper per z-score robusto con NA
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m) # fallback: solo centering
  (x - m) / s
}

# ------------------------------------------------------------
# Data preparation for joint model (fix obs_id after filtering)
# ------------------------------------------------------------

# helper per z-score robusto con NA
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m) # fallback: solo centering
  (x - m) / s
}

prepare_joint_ema_data <- function(d_clean) {
  d_clean <- d_clean %>% dplyr::mutate(day = as.Date(day))

  # 1) Item-level -----------------------------------------------------------
  item_data <- d_clean %>%
    dplyr::filter(!is.na(exam_period)) %>%
    dplyr::select(
      user_id,
      day,
      hour,
      exam_period,
      happy,
      sad,
      satisfied,
      angry,
      dass_stress,
      pid5_negative_affect_baseline,
      dass_stress_baseline,
      sex
    ) %>%
    dplyr::mutate(
      # 0–100 -> 1–7
      happy_ord = pmax(1L, pmin(7L, as.integer(round(happy / 100 * 6) + 1))),
      sad_ord = pmax(1L, pmin(7L, as.integer(round(sad / 100 * 6) + 1))),
      satisfied_ord = pmax(
        1L,
        pmin(7L, as.integer(round(satisfied / 100 * 6) + 1))
      ),
      angry_ord = pmax(1L, pmin(7L, as.integer(round(angry / 100 * 6) + 1))),
      period_numeric = dplyr::case_when(
        exam_period == "baseline" ~ 1L,
        exam_period == "pre_exam" ~ 2L,
        exam_period == "post_exam" ~ 3L,
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::filter(!is.na(period_numeric)) %>%
    dplyr::arrange(user_id, day, hour) %>%
    dplyr::mutate(time_numeric = as.integer(day - min(day, na.rm = TRUE)))

  # 2) Mapping soggetti -----------------------------------------------------
  subject_map <- item_data %>%
    dplyr::distinct(user_id) %>%
    dplyr::arrange(user_id) %>%
    dplyr::mutate(subject_numeric = dplyr::row_number())

  # 3) Dati osservazione (FILTRO PRIMA, poi assegna obs_id!) ---------------
  obs_data <- item_data %>%
    dplyr::left_join(subject_map, by = "user_id") %>%
    dplyr::filter(dplyr::if_any(
      c(happy_ord, sad_ord, satisfied_ord, angry_ord),
      ~ !is.na(.)
    )) %>%
    dplyr::arrange(subject_numeric, day, hour) %>% # stabilizza l'ordine
    dplyr::mutate(obs_id = dplyr::row_number()) # ora è 1..N_obs senza buchi

  # 4) Long per item --------------------------------------------------------
  items_long <- obs_data %>%
    dplyr::select(obs_id, happy_ord, sad_ord, satisfied_ord, angry_ord) %>%
    tidyr::pivot_longer(
      cols = tidyselect::ends_with("_ord"),
      names_to = "item_name",
      values_to = "response"
    ) %>%
    dplyr::filter(!is.na(response)) %>%
    dplyr::mutate(
      item_id = dplyr::case_when(
        item_name == "happy_ord" ~ 1L,
        item_name == "sad_ord" ~ 2L,
        item_name == "satisfied_ord" ~ 3L,
        item_name == "angry_ord" ~ 4L
      ),
      response = as.integer(response)
    ) %>%
    dplyr::arrange(obs_id, item_id)

  # 5) Predittori between con imputazione + dummies ------------------------
  between_predictors <- obs_data %>%
    dplyr::group_by(subject_numeric, user_id) %>%
    dplyr::summarise(
      baseline_pid5 = dplyr::first(pid5_negative_affect_baseline),
      baseline_dass = dplyr::first(dass_stress_baseline),
      sex_female = as.numeric(dplyr::first(sex) == "Femmina"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sex_female = dplyr::coalesce(sex_female, 0),
      z_baseline_pid5 = z_(baseline_pid5),
      z_baseline_dass = z_(baseline_dass),
      miss_pid5 = as.integer(is.na(z_baseline_pid5)),
      miss_dass = as.integer(is.na(z_baseline_dass)),
      z_baseline_pid5 = dplyr::coalesce(z_baseline_pid5, 0),
      z_baseline_dass = dplyr::coalesce(z_baseline_dass, 0)
    ) %>%
    dplyr::arrange(subject_numeric)

  # 6) DASS: indice ottimizzato dass_obs -----------------------------------
  dass_data <- obs_data %>%
    dplyr::filter(!is.na(dass_stress)) %>%
    dplyr::select(subject_numeric, period_numeric, dass_stress, obs_id) %>%
    dplyr::arrange(subject_numeric, period_numeric)

  dass_obs <- dass_data %>%
    dplyr::mutate(
      obs_id_match = purrr::pmap_int(
        list(subject_numeric, period_numeric),
        ~ obs_data %>%
          dplyr::filter(subject_numeric == ..1, period_numeric == ..2) %>%
          dplyr::summarise(
            first_obs = dplyr::first(obs_id),
            .groups = "drop"
          ) %>%
          dplyr::pull(first_obs) %>%
          {
            if (length(.) == 0 || is.na(.)) NA_integer_ else as.integer(.)
          }
      )
    ) %>%
    dplyr::pull(obs_id_match)

  if (any(is.na(dass_obs))) {
    keep <- !is.na(dass_obs)
    dass_data <- dass_data[keep, , drop = FALSE]
    dass_obs <- dass_obs[keep]
  }

  # 7) Lookup veloce (i,p) -> first obs ------------------------------------
  I <- max(obs_data$subject_numeric)
  P <- 3L
  first_obs_by_ip <- obs_data %>%
    dplyr::group_by(subject_numeric, period_numeric) %>%
    dplyr::summarise(first_obs = dplyr::first(obs_id), .groups = "drop") %>%
    tidyr::complete(
      subject_numeric = seq_len(I),
      period_numeric = 1:P,
      fill = list(first_obs = NA_integer_)
    ) %>%
    dplyr::arrange(subject_numeric, period_numeric) %>%
    dplyr::pull(first_obs) %>%
    {
      matrix(as.integer(.), nrow = I, ncol = P, byrow = TRUE)
    }
  first_obs_by_ip[is.na(first_obs_by_ip)] <- 0L

  # 8) Matrice W ------------------------------------------------------------
  W_matrix <- as.matrix(between_predictors[, c(
    "z_baseline_pid5",
    "z_baseline_dass",
    "sex_female",
    "miss_pid5",
    "miss_dass"
  )])
  stopifnot(!any(is.na(W_matrix)))

  # 9) Sanity checks sugli indici ------------------------------------------
  stopifnot(max(obs_data$obs_id) == nrow(obs_data))
  stopifnot(all(items_long$obs_id >= 1L & items_long$obs_id <= nrow(obs_data)))
  if (nrow(dass_data) > 0) {
    stopifnot(length(dass_obs) == nrow(dass_data))
    stopifnot(all(dass_obs >= 1L & dass_obs <= nrow(obs_data)))
  }

  # 10) Lista per Stan ------------------------------------------------------
  stan_data <- list(
    I = nrow(between_predictors),
    N_obs = nrow(obs_data),
    K = 4L,
    P = 3L,
    Q = ncol(W_matrix),

    N_items = nrow(items_long),
    y_item = as.integer(items_long$response),
    item_id = as.integer(items_long$item_id),
    obs_id = as.integer(items_long$obs_id),

    subject = as.integer(obs_data$subject_numeric),
    period = as.integer(obs_data$period_numeric),
    time_days = as.numeric(obs_data$time_numeric),

    W = W_matrix,

    N_dass = nrow(dass_data),
    dass_stress = as.numeric(dass_data$dass_stress),
    dass_obs = as.integer(dass_obs),
    first_obs_by_ip = first_obs_by_ip
  )

  list(
    stan_data = stan_data,
    obs_data = obs_data,
    subject_map = subject_map,
    between_predictors = between_predictors,
    items_long = items_long
  )
}


# Prepare data
cat("=== Preparing data for joint model ===\n")
# Carichiamo il dataset completo
d <- rio::import(
  here::here(
    "data",
    "processed",
    "ema_plus_baseline_exam_tagged.csv"
  )
)
joint_data <- prepare_joint_ema_data(d)

cat("Subjects:", joint_data$stan_data$I, "\n")
cat("EMA observations:", joint_data$stan_data$N_obs, "\n")
cat("Item responses:", joint_data$stan_data$N_items, "\n")
cat("DASS observations:", joint_data$stan_data$N_dass, "\n")

# (opzionale) controlla bounds prima di chiamare Stan
with(joint_data$stan_data, {
  stopifnot(all(obs_id >= 1 & obs_id <= N_obs))
  if (N_dass > 0) stopifnot(all(dass_obs >= 1 & dass_obs <= N_obs))
})

# ------------------------------------------------------------
# Fit joint model
# ------------------------------------------------------------

# Compile model
stan_file <- here::here(
  "scripts",
  "02_stress_reactivity",
  "joint_ema_latent_reactivity.stan"
)
# write_stan_file() can be used if you want to create the file directly in R
model_joint <- cmdstan_model(stan_file)

# cat("\n=== Fitting joint model with MCMC ===\n")
# fit_joint <- model_joint$sample(
#   data = joint_data$stan_data,
#   seed = 20250912,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 1000,
#   adapt_delta = 0.95,
#   max_treedepth = 12,
#   refresh = 200
# )
#
# # Check convergence
# fit_joint$cmdstan_diagnose()
# print(fit_joint$diagnostic_summary())

cat("\n=== Fitting joint model with ADVI: fullrank ===\n")
# ADVI: fullrank
fit_joint_vb_mf <- model_joint$variational(
  data = joint_data$stan_data,
  seed = 20250912,
  algorithm = "fullrank", # oppure "fullrank" sotto
  grad_samples = 1, # stima del gradiente per iter
  elbo_samples = 100, # MC per ELBO
  eta = 1, # stepsize (1 è spesso ok)
  adapt_engaged = TRUE, # preconditioning
  tol_rel_obj = 0.001, # criterio arresto (0.1%),
  eval_elbo = 100, # stampa ELBO ogni N iter
  output_samples = 2000, # draws dalla approx variational
  refresh = 200
)

# ------------------------------------------------------------
# Extract and analyze results (ADVI-friendly)
# ------------------------------------------------------------

# Usa l'oggetto fit disponibile: preferisci VB meanfield, poi fullrank, altrimenti HMC
fit <- if (exists("fit_joint_vb_mf")) {
  fit_joint_vb_mf
} else if (exists("fit_joint_vb_fr")) {
  fit_joint_vb_fr
} else if (exists("fit_joint")) {
  fit_joint
} else {
  stop(
    "Nessun oggetto 'fit' trovato (fit_joint_vb_mf / fit_joint_vb_fr / fit_joint)."
  )
}

# ------------------------------------------------------------
# Key parameters summary
# ------------------------------------------------------------
key_params <- c(
  "mu_0",
  "sigma_mu",
  "phi",
  "sigma_eta",
  "gamma",
  "beta",
  "lambda",
  "alpha_dass",
  "lambda_dass",
  "R2_dass"
)
print(fit$summary(variables = key_params), n = Inf)

# ------------------------------------------------------------
# Extract latent reactivity measures (reactivity_pre_post)
# ------------------------------------------------------------
reactivity_draws <- fit$draws("reactivity_pre_post", format = "df") %>%
  # Rimuovi eventuali colonne meta se presenti (in VB potrebbero non esserci)
  select(-any_of(c(".chain", ".iteration", ".draw"))) %>%
  # Riassumi per soggetto (colonne = reactivity_pre_post[1], ..., [I])
  summarise(across(
    everything(),
    list(
      mean = ~ mean(.x),
      sd = ~ sd(.x),
      q025 = ~ quantile(.x, 0.025),
      q975 = ~ quantile(.x, 0.975)
    )
  )) %>%
  pivot_longer(everything()) %>%
  # name è tipo "reactivity_pre_post[12]_mean": separa statistiche
  separate(
    name,
    into = c("variable", "stat"),
    sep = "_(?=[^_]+$)",
    remove = TRUE
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    subject_numeric = as.integer(str_extract(variable, "\\d+"))
  ) %>%
  # Togli soggetti senza stime (p.es. -Inf nei gq)
  filter(is.finite(mean)) %>%
  left_join(joint_data$between_predictors, by = "subject_numeric")

# ------------------------------------------------------------
# Extract eta values by period (eta_by_period)
# ------------------------------------------------------------
eta_period_draws <- fit$draws("eta_by_period", format = "df") %>%
  select(-any_of(c(".chain", ".iteration", ".draw"))) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(everything(), names_to = "name", values_to = "eta_mean") %>%
  mutate(
    # name è tipo "eta_by_period[ i , p ]"
    idx = str_match(name, "eta_by_period\\[(\\d+),(\\d+)\\]"),
    subject_numeric = as.integer(idx[, 2]),
    period = as.integer(idx[, 3])
  ) %>%
  select(-name, -idx) %>%
  # riga per soggetto × periodo
  pivot_wider(
    names_from = period,
    values_from = eta_mean,
    names_prefix = "eta_period_"
  ) %>%
  left_join(joint_data$between_predictors, by = "subject_numeric")

# ------------------------------------------------------------
# (Opzionale) estrai anche i contrasti di periodo
# ------------------------------------------------------------
deltas <- fit$draws(
  c("delta_pre_baseline", "delta_post_baseline", "delta_post_pre"),
  format = "df"
) %>%
  select(-any_of(c(".chain", ".iteration", ".draw"))) %>%
  summarise(across(
    everything(),
    list(
      mean = ~ mean(.x),
      sd = ~ sd(.x),
      q025 = ~ quantile(.x, 0.025),
      q975 = ~ quantile(.x, 0.975)
    )
  )) %>%
  pivot_longer(everything()) %>%
  separate(
    name,
    into = c("variable", "stat"),
    sep = "_(?=[^_]+$)",
    remove = TRUE
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    subject_numeric = as.integer(str_extract(variable, "\\d+")),
    contrast = str_remove(variable, "\\[\\d+\\]")
  ) %>%
  select(subject_numeric, contrast, mean, sd, q025, q975)

# ------------------------------------------------------------
# Analysis: Compare baseline vs EMA predictive power
# ------------------------------------------------------------

# k-fold CV per R^2 predittivo
kfold_predictive_r2 <- function(
  outcome_data,
  predictor_sets,
  k_folds = 10,
  seed = 123
) {
  set.seed(seed)
  n <- nrow(outcome_data)
  folds <- sample(rep(1:k_folds, length.out = n))

  map_dfr(
    predictor_sets,
    function(predictors) {
      r2_folds <- map_dbl(1:k_folds, function(k) {
        train_idx <- which(folds != k)
        test_idx <- which(folds == k)

        train_data <- outcome_data[train_idx, , drop = FALSE]
        test_data <- outcome_data[test_idx, , drop = FALSE]

        if (length(predictors) == 0) {
          y_pred <- rep(mean(train_data$outcome), length(test_idx))
        } else {
          fit_lm <- lm(
            outcome ~ .,
            data = cbind(
              outcome = train_data$outcome,
              train_data[, predictors, drop = FALSE]
            )
          )
          y_pred <- predict(
            fit_lm,
            newdata = test_data[, predictors, drop = FALSE]
          )
        }

        ss_res <- sum((test_data$outcome - y_pred)^2)
        ss_tot <- sum((test_data$outcome - mean(test_data$outcome))^2)
        1 - ss_res / ss_tot
      })

      tibble(
        model = if (length(predictors) == 0) "(Intercept)" else
          paste(predictors, collapse = " + "),
        r2_mean = mean(r2_folds),
        r2_se = sd(r2_folds) / sqrt(k_folds)
      )
    },
    .id = "model_set"
  )
}

# Dati per confronto (usa la mean di reactivity_pre_post per soggetto)
comparison_data <- reactivity_draws %>%
  filter(!is.na(mean), is.finite(mean)) %>%
  select(
    subject_numeric,
    outcome = mean,
    z_baseline_pid5,
    z_baseline_dass,
    sex_female
  )

predictor_sets <- list(
  "Intercept only" = character(0),
  "Baseline PID-5" = "z_baseline_pid5",
  "Baseline DASS" = "z_baseline_dass",
  "Sex" = "sex_female",
  "All baseline" = c("z_baseline_pid5", "z_baseline_dass", "sex_female")
)

cat("\n=== Cross-validation comparison ===\n")
cv_results <- kfold_predictive_r2(
  comparison_data,
  predictor_sets,
  k_folds = 10,
  seed = 2025
)
print(cv_results, n = Inf)

# ------------------------------------------------------------
# Visualizations
# ------------------------------------------------------------

# 1) Reactivity vs baseline PID-5
p1 <- reactivity_draws %>%
  filter(!is.na(mean), is.finite(mean)) %>%
  ggplot(aes(x = z_baseline_pid5, y = mean, color = factor(sex_female))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Latent Reactivity vs Baseline PID-5 Negative Affect",
    x = "Baseline PID-5 (z-scored)",
    y = "Latent Reactivity (max(pre,post) - baseline)",
    color = "Female"
  ) +
  theme_minimal()
print(p1)

# 2) Latent Negative Affect by Period (eta_by_period)
eta_period_summary <- eta_period_draws %>%
  select(subject_numeric, starts_with("eta_period_")) %>%
  pivot_longer(-subject_numeric, names_to = "period", values_to = "eta") %>%
  mutate(period = str_replace(period, "eta_period_", "Period ")) %>%
  filter(is.finite(eta))

p2 <- eta_period_summary %>%
  ggplot(aes(x = period, y = eta)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width = 0.2) +
  labs(
    title = "Latent Negative Affect by Exam Period",
    x = "Exam Period",
    y = "Latent Negative Affect (eta)"
  ) +
  theme_minimal()
print(p2)

# 3) Model validation: DASS R^2
dass_r2 <- fit$draws("R2_dass", format = "df") %>%
  select(-any_of(c(".chain", ".iteration", ".draw"))) %>%
  as_tibble() %>%
  pull(1) %>%
  mean()
cat("\n=== Model validation ===\n")
cat("DASS-Stress R² (mean):", round(dass_r2, 3), "\n")

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------
results_list <- list(
  fit = fit,
  reactivity_estimates = reactivity_draws,
  eta_by_period = eta_period_draws,
  deltas = deltas,
  cv_results = cv_results,
  dass_r2 = dass_r2
)
saveRDS(results_list, "joint_ema_model_results.rds")
cat("\nResults saved to: joint_ema_model_results.rds\n")

# ------------------------------------------------------------
# Answer research question (testo stampato)
# ------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("RESEARCH QUESTION ANALYSIS\n")
cat(strrep("=", 60), "\n", sep = "")

cat(
  "\nQuestion: Are baseline trait measures (PID-5) or repeated EMA measures\n"
)
cat("better at capturing psychological fragility (exam reactivity)?\n\n")

baseline_r2 <- cv_results %>%
  filter(stringr::str_detect(model, "z_baseline_pid5")) %>%
  slice(1) %>%
  pull(r2_mean)

cat("Preliminary Results:\n")
cat("- Baseline PID-5 predictive power (R²):", round(baseline_r2, 3), "\n")
cat(
  "- Joint model successfully estimated latent negative affect (ordinal items + method effects)\n"
)
cat("- DASS validation R²:", round(dass_r2, 3), "\n")

cat("\nNext steps:\n")
cat("1. Compare this joint model to a composite-score approach explicitly.\n")
cat("2. Quantify any increment in predictive accuracy.\n")
cat("3. Assess model stability and measurement invariance across periods.\n")
cat(strrep("=", 60), "\n", sep = "")
