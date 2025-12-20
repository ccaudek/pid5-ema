# ==============================================================================
# Bayesian Moderation Analysis: Composite PID-5 Score × Stress Reactivity
# Testing context-dependent expression of personality pathology
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayestestR)
  library(bayesplot)
  library(tidybayes)
  library(patchwork)
  library(ggdist)
  library(marginaleffects)
  library(purrr)
})

options(brms.backend = "cmdstanr")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS: PID-5 × Stress Reactivity\n")
cat("All 5 PID-5 Domains\n")
cat(rep("=", 70), "\n", sep = "")

# ==============================================================================
# 1. PREPARE DATA FOR MODERATION MODELS
# ==============================================================================

# Check PID-5 distributions at baseline
cat("\n=== PID-5 Baseline Trait Distributions (All 5 Domains) ===\n")
df_analysis %>%
  dplyr::filter(timepoint == "baseline") %>%
  dplyr::select(ends_with("_bl_c")) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    M = mean(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Min = min(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# Create contrast-coded timepoint variable for interaction interpretation
# C1: PRE vs BASELINE (stress reactivity)
# C2: POST vs PRE (stress recovery)

df_analysis <- df_analysis %>%
  mutate(
    # Contrast 1: PRE vs BASELINE (test stress reactivity)
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre" ~ 0.5,
      timepoint == "post" ~ 0
    ),
    # Contrast 2: POST vs PRE (test recovery)
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre" ~ -0.5,
      timepoint == "post" ~ 0.5
    )
  )


df_analysis <- df_analysis %>%
  mutate(
    pid5_mean_z = rowMeans(
      scale(select(
        .,
        pid5_negative_affectivity_bl,
        pid5_detachment_bl,
        pid5_antagonism_bl,
        pid5_disinhibition_bl,
        pid5_psychoticism_bl
      )),
      na.rm = TRUE
    ),
    pid5_mean_z = as.numeric(scale(pid5_mean_z)) # standardizza il composito
  )

df_analysis$f2_std_u <- ifelse(
  df_analysis$f2_std_u == 0,
  12.86300,
  df_analysis$f2_std_u
)

# # PCA sulle 5 dimensioni PID-5 baseline
# pca_pid5 <- df_analysis %>%
#   dplyr::filter(timepoint == "baseline") %>%
#   dplyr::select(ID, ends_with("_bl_c")) %>%
#   distinct()
#
# pca_result <- prcomp(pca_pid5 %>% dplyr::select(-ID), scale. = TRUE)
#
# # Estrai PC1 (general psychopathology)
# pca_pid5$p_factor <- pca_result$x[, 1]
#
# ee <- eigen(cor(pca_pid5[, 2:6]))$values
# ee[1]/sum(ee)
#
# pca_pid5$p_factor <- rowMeans(scale(pca_pid5[, 2:6]))
# # cor(mean_score, pca_pid5$p_factor)
# # [1] -0.9989956
# # la correlazione è molto alta; le due versioni sono praticamente equivalenti.
# # Uso la media perché preserva il segno.
#
# # Aggiungi al dataframe
# df_analysis <- df_analysis %>%
#   left_join(pca_pid5 %>% dplyr::select(ID, p_factor), by = "ID") %>%
#   mutate(p_factor_c = scale(p_factor, center = TRUE, scale = TRUE)[, 1])

# Modello di mediazione con il punteggio composito PID-5
formula_pfactor <- bf(
  f2_std_u ~
    c1_stress *
      pid5_mean_z +
      c2_recovery * pid5_mean_z +
      (1 + c1_stress + c2_recovery | ID)
)

m_pfactor <- brm(
  formula_pfactor,
  data = df_analysis,
  family = lognormal(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4
)
pp_check(m_pfactor)
summary(m_pfactor)
bayes_R2(m_pfactor)

conditional_effects(m_pfactor, "p_factor_c")


# --- 0) impostazioni utili (modifica se vuoi)
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
seed <- 12345
control_list <- list(adapt_delta = 0.95, max_treedepth = 15)

# --- 1) outcome + family (stringhe) in un tibble
outcomes_tbl <- tribble(
  ~var,
  ~family_name,
  "f0_mean_a",
  "gaussian",
  "f0_std_a",
  "lognormal",
  "jitter_a",
  "lognormal",
  "nne_a",
  "gaussian",
  "f2_mean_a",
  "student",
  "f2_std_a",
  "lognormal",
  "f0_mean_i",
  "gaussian",
  "f0_std_i",
  "lognormal",
  "jitter_i",
  "lognormal",
  "nne_i",
  "gaussian",
  "f2_mean_i",
  "student",
  "f2_std_i",
  "lognormal",
  "f0_mean_u",
  "gaussian",
  "f0_std_u",
  "lognormal",
  "jitter_u",
  "lognormal",
  "nne_u",
  "gaussian",
  "f2_mean_u",
  "student",
  "f2_std_u",
  "lognormal"
)

# --- 2) helper: costruisce la bf dinamica per un outcome
make_bf <- function(y) {
  # la formula che volevi
  bf(as.formula(paste0(
    y,
    " ~ c1_stress * pid5_mean_z + c2_recovery * pid5_mean_z + (1 + c1_stress + c2_recovery | ID)"
  )))
}

# --- 3) helper: ritorna l'oggetto family() dato il nome (string)
family_from_name <- function(fname) {
  # mappa name -> funzione family; match.fun funziona se esiste la funzione
  # (es. "gaussian", "student", "lognormal", "skew_normal", ...)
  # qui usiamo match.fun ma proteggiamo con tryCatch
  tryCatch(
    match.fun(fname)(),
    error = function(e)
      stop("Family non trovata: ", fname, " (", e$message, ")")
  )
}

# --- 4) controlla priors_moderation esistente, altrimenti fallback
if (!exists("priors_moderation")) {
  priors_moderation <- c(
    prior(normal(0, 1), class = "b"),
    prior(student_t(3, 0, 2.5), class = "Intercept")
  )
}

# --- 5) funzione che esegue il fit in modo sicuro e ritorna brmsfit o NULL
fit_model_safe <- function(var, family_name, data) {
  fam_obj <- NULL
  try(
    {
      fam_obj <- family_from_name(family_name)
    },
    silent = TRUE
  )
  if (is.null(fam_obj)) {
    message("Skipping ", var, ": family '", family_name, "' non valida.")
    return(NULL)
  }

  f <- make_bf(var)
  message("Fitting ", var, " (family=", family_name, ") ...")
  fit <- tryCatch(
    brm(
      formula = f,
      data = data,
      family = fam_obj,
      prior = priors_moderation,
      iter = iter,
      warmup = warmup,
      chains = chains,
      cores = cores,
      control = control_list,
      seed = seed,
      refresh = 0
    ),
    error = function(e) {
      message("Errore per ", var, ": ", e$message)
      NULL
    },
    warning = function(w) {
      # cattura warning ma non interrompe
      message("Warning (", var, "): ", w$message)
      invokeRestart("muffleWarning")
    }
  )
  fit
}

# --- 6) ciclo di fitting, salva risultati in lista nominata
models <- vector("list", nrow(outcomes_tbl))
names(models) <- outcomes_tbl$var

for (i in seq_len(nrow(outcomes_tbl))) {
  var_i <- outcomes_tbl$var[i]
  fam_i <- outcomes_tbl$family_name[i]
  models[[var_i]] <- fit_model_safe(var_i, fam_i, df_analysis)
}

# --- 7) estrazione parametri d'interesse per ogni modello
# parametri che vuoi riportare
pars_of_interest <- c(
  "c1_stress",
  "pid5_mean_z",
  "c2_recovery",
  "c1_stress:pid5_mean_z",
  "pid5_mean_z:c2_recovery"
)

extract_fixed_summary <- function(fit, var_name, probs = c(0.05, 0.95)) {
  if (is.null(fit))
    return(tibble(
      outcome = var_name,
      parameter = NA,
      Estimate = NA_real_,
      Est.Error = NA_real_,
      Lower = NA_real_,
      Upper = NA_real_
    ))
  fe <- tryCatch(
    fixef(fit, probs = probs) %>% as.data.frame(),
    error = function(e) return(NULL)
  )
  if (is.null(fe))
    return(tibble(
      outcome = var_name,
      parameter = NA,
      Estimate = NA_real_,
      Est.Error = NA_real_,
      Lower = NA_real_,
      Upper = NA_real_
    ))
  # normalizza i nomi delle colonne per Lower/Upper
  colnames(fe) <- colnames(fe) %>% str_replace_all("%", "") %>% make.names()
  # possibili nomi per le quantili: "Q5","Q95" o "Q2.5","Q97.5" o "X2.5." etc.
  # Cerchiamo la colonna di lower/upper automaticamente
  lower_col <- intersect(
    colnames(fe),
    c("Q5", "Q2.5", "X2.5.", "X5.", "Q0.05", "lower", "Lower")
  )[1]
  upper_col <- intersect(
    colnames(fe),
    c("Q95", "Q97.5", "X97.5.", "X95.", "Q0.95", "upper", "Upper")
  )[1]
  if (is.null(lower_col) || is.na(lower_col))
    lower_col <- colnames(fe)[ncol(fe) - 1]
  if (is.null(upper_col) || is.na(upper_col))
    upper_col <- colnames(fe)[ncol(fe)]
  fe_df <- fe %>%
    rownames_to_column("parameter") %>%
    as_tibble() %>%
    filter(parameter %in% pars_of_interest) %>%
    transmute(
      outcome = var_name,
      parameter,
      Estimate = as.numeric(.data$Estimate),
      Est.Error = as.numeric(.data$Est.Error),
      Lower = as.numeric(.data[[lower_col]]),
      Upper = as.numeric(.data[[upper_col]])
    )
  fe_df
}

# costruisco results_all unendo i modelli
results_all <- map2_dfr(models, names(models), extract_fixed_summary)

# --- 8) Aggiungo Bayes R2 per ogni modello (se disponibile)
get_bayesR2 <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  tryCatch(bayes_R2(fit)$Estimate, error = function(e) NA_real_)
}

r2_tbl <- tibble(
  outcome = names(models),
  bayes_R2 = map_dbl(models, get_bayesR2)
)

# unisco R2 a results_all
results_all <- results_all %>%
  left_join(r2_tbl, by = "outcome")

# salva tabella compatta
write_csv(results_all, "brms_fixed_effects_summary.csv")

# --- 9) grafico riepilogativo
# metto valori NA in fondo per ordine coerente
results_all_plot <- results_all %>% dplyr::filter(!is.na(parameter))

ggplot(
  results_all_plot,
  aes(x = Estimate, y = fct_rev(factor(outcome)), xmin = Lower, xmax = Upper)
) +
  geom_point() +
  geom_errorbarh(height = 0.15) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = "Estimate (posterior mean)", y = "Outcome (variabile dipendente)") +
  theme_minimal(base_size = 12)


##################
library(tidyLPA)

# Calcola change scores per ciascun individuo
change_scores <- df_analysis %>%
  select(ID, timepoint, f0_mean_i, f0_std_i, jitter_i, nne_i) %>%
  pivot_wider(names_from = timepoint, values_from = c(f0_mean_i:nne_i)) %>%
  mutate(
    # Reactivity = PRE - BASELINE
    react_f0mean = f0_mean_i_pre - f0_mean_i_baseline,
    react_f0std = f0_std_i_pre - f0_std_i_baseline,
    react_jitter = jitter_i_pre - jitter_i_baseline,
    react_nne = nne_i_pre - nne_i_baseline,
    # Recovery = POST - PRE
    recov_f0mean = f0_mean_i_post - f0_mean_i_pre,
    recov_f0std = f0_std_i_post - f0_std_i_pre,
    recov_jitter = jitter_i_post - jitter_i_pre,
    recov_nne = nne_i_post - nne_i_pre
  )

# LPA sui pattern di reattività
lpa_input <- change_scores %>%
  select(starts_with("react_")) %>%
  scale()

# Testa 2-4 profili
lpa_results <- lpa_input %>%
  estimate_profiles(2:4, variances = "equal", covariances = "zero")

# Scegli modello ottimale
get_fit(lpa_results)

# Estrai cluster membership
best_model <- get_data(lpa_results) # scegli n_profiles ottimale
change_scores$profile <- best_model$Class

# Valida con PID-5
change_scores %>%
  left_join(
    df_analysis %>%
      filter(timepoint == "baseline") %>%
      select(ID, starts_with("pid5_"), -ends_with("_c")),
    by = "ID"
  ) %>%
  group_by(profile) %>%
  summarise(across(starts_with("pid5_"), mean, na.rm = TRUE))

# Plot profili
change_scores %>%
  select(profile, starts_with("react_")) %>%
  pivot_longer(-profile) %>%
  ggplot(aes(x = name, y = value, color = factor(profile), group = profile)) +
  stat_summary(fun = mean, geom = "line", size = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  labs(
    title = "Latent Profiles of Stress Reactivity",
    x = "Acoustic Feature",
    y = "Reactivity (Change Score)",
    color = "Profile"
  ) +
  theme_minimal()
