suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(lubridate)
  library(rio)
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(Matrix)
  library(bayesplot)
  library(conflicted)
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("sd", "stats")
conflict_prefer("var", "stats")
conflict_prefer("mad", "stats")
conflict_prefer("match", "base")
conflict_prefer("lag", "dplyr")

# Funzione z-score sicura
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m)
  (x - m) / s
}


# Carica le funzioni di preparazione dati (assumendo che siano già in memoria)
# source("path_to_data_preparation_script.R")  # se in file separato

set.seed(20250917)

# ---------- Analisi Completa ----------

run_complete_enhanced_analysis <- function(data_path, stan_file_path) {
  cat("=== ANALISI FRAGILITÀ PSICOLOGICA ENHANCED ===\n\n")

  # 1) Carica e prepara i dati
  cat("1. Caricamento e preparazione dati...\n")

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

  frag_data_enh <- prepare_fragility_data_enhanced(d, min_obs_per_period = 1)

  # 2) Compila il modello Stan
  cat("\n2. Compilazione modello Stan...\n")
  model_frag <- cmdstan_model(stan_file_path)

  # 3) Fit Modello A: solo baseline
  cat("\n3. Fit Modello A (baseline-only)...\n")
  data_A <- frag_data_enh$stan_data
  data_A$use_ema <- 0.0

  fit_A <- model_frag$variational(
    data = data_A,
    seed = 20250917,
    algorithm = "meanfield",
    elbo_samples = 100,
    adapt_engaged = TRUE,
    tol_rel_obj = 0.001,
    eval_elbo = 100,
    output_samples = 2000,
    refresh = 0 # silenzioso
  )

  # 4) Fit Modello B: baseline + EMA
  cat("4. Fit Modello B (baseline + EMA)...\n")
  data_B <- frag_data_enh$stan_data
  data_B$use_ema <- 1.0

  fit_B <- model_frag$variational(
    data = data_B,
    seed = 20250917,
    algorithm = "meanfield",
    elbo_samples = 100,
    adapt_engaged = TRUE,
    tol_rel_obj = 0.001,
    eval_elbo = 100,
    output_samples = 2000,
    refresh = 0 # silenzioso
  )

  # 5) Confronto R² Bayesiano
  cat("\n5. Confronto R² Bayesiano...\n")

  R2_A <- as_draws_df(fit_A$draws("R2_frag"))$R2_frag
  R2_B <- as_draws_df(fit_B$draws("R2_frag"))$R2_frag

  cat(sprintf("R² baseline-only:     %.3f (SD = %.3f)\n", mean(R2_A), sd(R2_A)))
  cat(sprintf("R² baseline + EMA:    %.3f (SD = %.3f)\n", mean(R2_B), sd(R2_B)))

  delta_R2 <- R2_B - R2_A
  cat(sprintf(
    "ΔR² (B - A):          %.3f [95%% CrI: %.3f, %.3f]\n",
    mean(delta_R2),
    quantile(delta_R2, 0.025),
    quantile(delta_R2, 0.975)
  ))

  prob_positive <- mean(delta_R2 > 0)
  cat(sprintf("P(ΔR² > 0):           %.3f\n", prob_positive))

  # 6) Analisi coefficienti chiave (solo per modello B)
  cat("\n6. Analisi coefficienti chiave (Modello B)...\n")

  # Estrai coefficienti di interesse
  b_female <- as_draws_df(fit_B$draws("b_female"))$b_female
  c_ema <- as_draws_df(fit_B$draws("c_ema"))
  b_var_ema <- as_draws_df(fit_B$draws("b_var_ema"))$b_var_ema
  b_var_female <- as_draws_df(fit_B$draws("b_var_female"))$b_var_female

  # Effetto genere principale
  cat(sprintf(
    "Effetto genere:       %.3f [%.3f, %.3f] | P(>0) = %.3f\n",
    mean(b_female),
    quantile(b_female, 0.025),
    quantile(b_female, 0.975),
    mean(b_female > 0)
  ))

  # Coefficienti EMA (5 dimensioni)
  cat("Coefficienti EMA principali:\n")
  dim_names <- c(
    "NA_latent",
    "Detachment",
    "Antagonism",
    "Disinhibition",
    "Psychoticism"
  )
  for (i in 1:5) {
    coef_i <- c_ema[[paste0("c_ema[", i, "]")]]
    cat(sprintf(
      "  %s: %.3f [%.3f, %.3f] | P(>0) = %.3f\n",
      dim_names[i],
      mean(coef_i),
      quantile(coef_i, 0.025),
      quantile(coef_i, 0.975),
      mean(coef_i > 0)
    ))
  }

  # Effetto variabilità emotiva
  cat(sprintf(
    "Effetto variabilità EMA: %.3f [%.3f, %.3f] | P(>0) = %.3f\n",
    mean(b_var_ema),
    quantile(b_var_ema, 0.025),
    quantile(b_var_ema, 0.975),
    mean(b_var_ema > 0)
  ))

  cat(sprintf(
    "Moderazione genere × variabilità: %.3f [%.3f, %.3f] | P(>0) = %.3f\n",
    mean(b_var_female),
    quantile(b_var_female, 0.025),
    quantile(b_var_female, 0.975),
    mean(b_var_female > 0)
  ))

  # 7) LOO Comparison
  cat("\n7. Leave-One-Out Cross-Validation...\n")

  ll_A <- fit_A$draws("log_lik_obs", format = "matrix")
  ll_B <- fit_B$draws("log_lik_obs", format = "matrix")

  # LOO a livello osservazione
  loo_A <- loo::loo(ll_A, moment_match = TRUE)
  loo_B <- loo::loo(ll_B, moment_match = TRUE)

  cmp_obs <- loo::loo_compare(list(B = loo_B, A = loo_A))

  cat("LOO Comparison (osservazione-level):\n")
  print(cmp_obs)

  delta_elpd_obs <- as.numeric(cmp_obs["A", "elpd_diff"]) * -1
  se_diff_obs <- as.numeric(cmp_obs["A", "se_diff"])

  cat(sprintf("\nΔELPD (B-A): %.1f ± %.1f\n", delta_elpd_obs, se_diff_obs))
  cat(sprintf(
    "Decisivo se |ΔE| > 2×SE: %s (|%.1f| %s %.1f)\n",
    ifelse(abs(delta_elpd_obs) > 2 * se_diff_obs, "SÌ", "NO"),
    delta_elpd_obs,
    ifelse(abs(delta_elpd_obs) > 2 * se_diff_obs, ">", "≤"),
    2 * se_diff_obs
  ))

  # 8) LOO a livello soggetto (aggregato)
  cat("\n8. LOO Comparison (subject-level)...\n")

  subj_idx <- frag_data_enh$stan_data$subject
  I <- max(subj_idx)
  N_obs <- length(subj_idx)

  # Matrice di raggruppamento sparse
  G <- sparseMatrix(i = subj_idx, j = seq_len(N_obs), x = 1, dims = c(I, N_obs))

  ll_A_subj <- ll_A %*% t(as.matrix(G))
  ll_B_subj <- ll_B %*% t(as.matrix(G))

  loo_A_subj <- loo::loo(ll_A_subj, moment_match = TRUE)
  loo_B_subj <- loo::loo(ll_B_subj, moment_match = TRUE)

  cmp_subj <- loo::loo_compare(list(B = loo_B_subj, A = loo_A_subj))
  print(cmp_subj)

  delta_elpd_subj <- as.numeric(cmp_subj["A", "elpd_diff"]) * -1
  se_diff_subj <- as.numeric(cmp_subj["A", "se_diff"])

  cat(sprintf(
    "\nΔELPD subject-level (B-A): %.1f ± %.1f\n",
    delta_elpd_subj,
    se_diff_subj
  ))
  cat(sprintf(
    "Decisivo: %s (|%.1f| %s %.1f)\n",
    ifelse(abs(delta_elpd_subj) > 2 * se_diff_subj, "SÌ", "NO"),
    delta_elpd_subj,
    ifelse(abs(delta_elpd_subj) > 2 * se_diff_subj, ">", "≤"),
    2 * se_diff_subj
  ))

  # 9) Verdetto finale
  cat("\n=== VERDETTO FINALE ===\n")

  winner_obs <- rownames(cmp_obs)[1]
  winner_subj <- rownames(cmp_subj)[1]

  cat(sprintf(
    "R²: Modello %s ha R² maggiore (Δ = %.3f, P(Δ>0) = %.2f)\n",
    ifelse(mean(delta_R2) > 0, "B", "A"),
    abs(mean(delta_R2)),
    ifelse(mean(delta_R2) > 0, prob_positive, 1 - prob_positive)
  ))

  cat(sprintf(
    "LOO beep-level: Modello %s preferito (ΔELPD = %.1f)\n",
    winner_obs,
    abs(delta_elpd_obs)
  ))

  cat(sprintf(
    "LOO subject-level: Modello %s preferito (ΔELPD = %.1f)\n",
    winner_subj,
    abs(delta_elpd_subj)
  ))

  # Evidenza decisiva?
  decisive_obs <- abs(delta_elpd_obs) > 2 * se_diff_obs
  decisive_subj <- abs(delta_elpd_subj) > 2 * se_diff_subj
  r2_meaningful <- abs(mean(delta_R2)) > 0.02 # soglia arbitraria per rilevanza pratica

  if (decisive_obs || decisive_subj) {
    if (winner_obs == "B" || winner_subj == "B") {
      cat("\n🎯 EVIDENZA A FAVORE DELLE MISURE EMA:\n")
      cat("   Le misure EMA del PID-5 migliorano significativamente\n")
      cat("   la predizione della fragilità psicologica rispetto\n")
      cat("   alle sole misure baseline.\n")
    } else {
      cat("\n⚠️  EVIDENZA CONTRO LE MISURE EMA:\n")
      cat("   Le misure baseline da sole sono sufficienti;\n")
      cat("   le EMA non aggiungono valore predittivo.\n")
    }
  } else {
    cat("\n🤷 EVIDENZA INCONCLUSIVA:\n")
    cat("   La differenza tra i modelli non è statisticamente\n")
    cat("   decisiva. Servono più dati o migliore modellazione.\n")
  }

  # 10) Salva risultati
  results <- list(
    data_summary = list(
      n_subjects = I,
      n_observations = N_obs,
      n_items = frag_data_enh$stan_data$N_items,
      n_female = sum(frag_data_enh$stan_data$female),
      n_male = sum(1 - frag_data_enh$stan_data$female)
    ),

    bayesian_r2 = list(
      baseline_only = list(mean = mean(R2_A), sd = sd(R2_A)),
      baseline_ema = list(mean = mean(R2_B), sd = sd(R2_B)),
      delta = list(
        mean = mean(delta_R2),
        ci_lower = quantile(delta_R2, 0.025),
        ci_upper = quantile(delta_R2, 0.975),
        prob_positive = prob_positive
      )
    ),

    loo_comparison = list(
      obs_level = list(
        delta_elpd = delta_elpd_obs,
        se_diff = se_diff_obs,
        winner = winner_obs,
        decisive = decisive_obs
      ),
      subj_level = list(
        delta_elpd = delta_elpd_subj,
        se_diff = se_diff_subj,
        winner = winner_subj,
        decisive = decisive_subj
      )
    ),

    key_coefficients = list(
      gender_effect = list(
        mean = mean(b_female),
        ci_lower = quantile(b_female, 0.025),
        ci_upper = quantile(b_female, 0.975),
        prob_positive = mean(b_female > 0)
      ),
      variability_effect = list(
        mean = mean(b_var_ema),
        ci_lower = quantile(b_var_ema, 0.025),
        ci_upper = quantile(b_var_ema, 0.975),
        prob_positive = mean(b_var_ema > 0)
      ),
      gender_variability_interaction = list(
        mean = mean(b_var_female),
        ci_lower = quantile(b_var_female, 0.025),
        ci_upper = quantile(b_var_female, 0.975),
        prob_positive = mean(b_var_female > 0)
      )
    ),

    verdict = list(
      ema_beneficial = winner_obs == "B" || winner_subj == "B",
      evidence_strength = ifelse(
        decisive_obs || decisive_subj,
        "strong",
        "weak"
      ),
      recommendation = ifelse(
        (winner_obs == "B" || winner_subj == "B") &&
          (decisive_obs || decisive_subj),
        "Use EMA measures - they improve prediction",
        ifelse(
          winner_obs == "A" &&
            winner_subj == "A" &&
            (decisive_obs || decisive_subj),
          "Baseline measures sufficient - EMA adds no value",
          "Evidence inconclusive - need more data or better model"
        )
      )
    )
  )

  return(list(
    results = results,
    fits = list(baseline_only = fit_A, baseline_ema = fit_B),
    data = frag_data_enh
  ))
}

# Esempio di uso:
analysis_results <- run_complete_enhanced_analysis(
  data_path = here::here("data", "processed", "ema_plus_scales_cleaned.csv"),
  stan_file_path = here::here(
    "scripts",
    "02_stress_reactivity",
    "fragility.stan"
  )
)


sd <- analysis_results$data$stan_data
y_obs <- sd$y_item # vettore lungo N_items
yrep <- posterior::as_draws_matrix(
  analysis_results$fits$baseline_ema$draws("y_item_rep")
) # (ndraws x N_items)

# Confronto delle distribuzioni di categoria (istogrammi)
ppc_bars(y_obs, yrep[1:50, ])

# Oppure rootogram per discreti
ppc_rootogram(y_obs, yrep[1:50, ])


# PPC sugli item ordinali

# osservati
y_item_obs <- analysis_results$data$stan_data$y_item

# repliche (matrice: draws x N_items)
y_item_rep_mat <- posterior::as_draws_matrix(
  analysis_results$fits$baseline_ema$draws("y_item_rep")
)

# Istogrammi per categorie
ppc_bars(y_item_obs, y_item_rep_mat[1:200, ])

# (opzionale) rootogram per discreti
ppc_rootogram(y_item_obs, y_item_rep_mat[1:200, ])


# PPC sulla fragilità pre–post per soggetto
# --- INPUT ---
sd <- analysis_results$data$stan_data

# Draws degli item replicati: (ndraws x N_items)
Yrep <- posterior::as_draws_matrix(
  analysis_results$fits$baseline_ema$draws("y_item_rep")
)

# Allinea gli item verso NA: inverti happy (1) e satisfied (3) -> 8 - y
invert <- sd$item_id %in% c(1L, 3L)
Yrep_na <- Yrep
Yrep_na[, invert] <- 8 - Yrep_na[, invert]

# 1) Media per osservazione (beep): ogni osservazione ha 4 item
#    rowsum somma per gruppi di righe (qui: gli item, raggruppati per obs_id)
#    t(Yrep_na) è (N_items x ndraws); il risultato è (N_obs x ndraws)
NA_obs_t <- rowsum(t(Yrep_na), group = sd$obs_id, reorder = FALSE) / 4

# 2) Medie per soggetto in PRE (period == 2) e POST (period == 3)
I <- sd$I
is_pre <- sd$period == 2L
is_post <- sd$period == 3L
ndraws <- ncol(NA_obs_t)

avg_by_subj <- function(rows_flag) {
  out <- matrix(NA_real_, nrow = I, ncol = ndraws)
  idx <- sd$subject[rows_flag] # soggetti corrispondenti alle righe selezionate
  sums <- rowsum(
    NA_obs_t[rows_flag, , drop = FALSE],
    group = idx,
    reorder = FALSE
  )
  cnt <- tabulate(idx, nbins = I)
  present <- which(cnt > 0)
  out[present, ] <- sums[present, , drop = FALSE] / cnt[present]
  out # (I x ndraws)
}

PRE_t <- avg_by_subj(is_pre) # (I x ndraws)
POST_t <- avg_by_subj(is_post) # (I x ndraws)

# 3) Differenza pre - post per soggetto e draw
DIFF_t <- PRE_t - POST_t # (I x ndraws)
diff_rep_mat <- t(DIFF_t) # (ndraws x I) -> pronto per bayesplot

# 4) Costruisci la proxy osservata y_obs (se non l'hai già)
library(dplyr)
library(tidyr)

df <- tibble(
  n = sd$obs_id,
  subj = sd$subject[sd$obs_id],
  per = sd$period[sd$obs_id],
  item = sd$item_id,
  y = sd$y_item
) |>
  mutate(y_na = if_else(item %in% c(1, 3), 8L - y, y))

m_obs <- df |>
  group_by(n, subj, per) |>
  summarise(na_mean = mean(y_na), .groups = "drop")

y_obs <- m_obs |>
  group_by(subj) |>
  summarise(
    pre = if (any(per == 2)) mean(na_mean[per == 2]) else NA_real_,
    post = if (any(per == 3)) mean(na_mean[per == 3]) else NA_real_
  ) |>
  mutate(diff = pre - post) |>
  arrange(subj) |>
  pull(diff)

keep <- is.finite(y_obs)
y_obs <- y_obs[keep]
diff_rep_mat <- diff_rep_mat[, keep, drop = FALSE]

# 5) PPC: scatter e intervalli predittivi
bayesplot::ppc_scatter_avg(y = y_obs, yrep = diff_rep_mat)
bayesplot::ppc_intervals(
  y = y_obs,
  yrep = diff_rep_mat,
  prob = 0.5,
  prob_outer = 0.9
)
