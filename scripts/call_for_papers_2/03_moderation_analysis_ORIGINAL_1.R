# ==============================================================================
# Bayesian Moderation Analysis: PID-5 Traits × Stress Reactivity
# Testing context-dependent expression of personality pathology
# ==============================================================================
# VERSIONE CORRETTA: Ripristina l'analisi originale (senza log-transform)
# Modifiche rispetto alla versione modificata:
#   1. RIMOSSA trasformazione logaritmica su f0_std, f2_std, jitter
#   2. CAMBIATO f2_mean: asym_laplace() -> gaussian()
#   3. AGGIORNATE le prior per usare valori sulla scala originale (non log)
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
  library(tidyr)
  library(stringr)
})

options(brms.backend = "cmdstanr")
options(max.print = 5000)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS: PID-5 × Stress Reactivity\n")
cat("All 5 PID-5 Domains\n")
cat("VERSIONE ORIGINALE (senza log-transform)\n")
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

# Gestione valori zero per f2_std_i (se necessario)
df_analysis$f2_std_i <- ifelse(
  df_analysis$f2_std_i == 0,
  min(df_analysis$f2_std_i[df_analysis$f2_std_i > 0]),
  df_analysis$f2_std_i
)

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

# ==============================================================================
# NOTA: RIMOSSA la trasformazione logaritmica che era presente nella versione
# modificata. I dati rimangono sulla scala originale (Hz, %, etc.)
# ==============================================================================
# Le seguenti righe sono state RIMOSSE:
# df_analysis <- df_analysis %>%
#   mutate(across(
#     matches("^(f0_std|f2_std|jitter)_"),
#     ~ log(.x)
#   ))
# ==============================================================================

# --- Gestione valori zero/negativi per outcome lognormal ----------------------
# Lognormal richiede valori > 0. Sostituiamo zeri con il minimo valore positivo.
lognormal_vars <- c("f0_std", "f2_std", "jitter")
vowels_all <- c("a", "i", "u")

for (var in lognormal_vars) {
  for (v in vowels_all) {
    col_name <- paste0(var, "_", v)
    if (col_name %in% names(df_analysis)) {
      # Trova valori <= 0
      bad_idx <- which(df_analysis[[col_name]] <= 0)
      if (length(bad_idx) > 0) {
        # Sostituisci con il minimo valore positivo
        min_pos <- min(
          df_analysis[[col_name]][df_analysis[[col_name]] > 0],
          na.rm = TRUE
        )
        df_analysis[[col_name]][bad_idx] <- min_pos
        cat(
          "Nota:",
          length(bad_idx),
          "valori <= 0 in",
          col_name,
          "sostituiti con",
          round(min_pos, 4),
          "\n"
        )
      }
    }
  }
}

# --- Parametri empirici calcolati dinamicamente dai dati ---------------------
compute_vowel_means <- function(data) {
  outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
  vowels <- c("a", "i", "u")

  # Outcome che usano lognormal (calcoliamo media geometrica, devono essere > 0)
  lognormal_outcomes <- c("f0_std", "f2_std", "jitter")

  vowel_means <- list()

  for (v in vowels) {
    vowel_means[[v]] <- list()

    for (out in outcomes) {
      col_name <- paste0(out, "_", v)

      if (col_name %in% names(data)) {
        values <- data[[col_name]]
        values <- values[!is.na(values)] # rimuovi solo NA

        if (out %in% lognormal_outcomes) {
          # Media geometrica per outcome lognormal (richiede valori > 0)
          values <- values[values > 0]
          vowel_means[[v]][[out]] <- exp(mean(log(values)))
        } else {
          # Media aritmetica per outcome gaussian/student (può includere negativi)
          vowel_means[[v]][[out]] <- mean(values)
        }
      }
    }
  }

  return(vowel_means)
}

# Calcola i parametri dai dati
vowel_means <- compute_vowel_means(df_analysis)

# Verifica
cat("\n=== Computed Vowel Means for Priors ===\n")
for (v in names(vowel_means)) {
  cat("\nVowel /", v, "/:\n", sep = "")
  for (out in names(vowel_means[[v]])) {
    cat("  ", out, ": ", round(vowel_means[[v]][[out]], 3), "\n", sep = "")
  }
}

make_priors_vowel <- function(outcome, vowel) {
  vm <- vowel_means[[vowel]]

  fmt <- function(x, digits = 6) {
    formatC(x, format = "f", digits = digits)
  }

  # -------------------------------
  # Intercept priors
  # For gaussian: on original scale
  # For lognormal: on LOG scale (because brms models log(y))
  # -------------------------------
  intercept_prior <- switch(
    outcome,

    # GAUSSIAN MODELS - intercept on original scale
    "f0_mean" = prior_string(
      paste0("student_t(3, ", fmt(vm$f0_mean), ", 30)"),
      class = "Intercept"
    ),

    "nne" = prior_string(
      paste0("student_t(3, ", fmt(vm$nne), ", 5)"),
      class = "Intercept"
    ),

    "f2_mean" = prior_string(
      paste0("student_t(3, ", fmt(vm$f2_mean), ", 150)"),
      class = "Intercept"
    ),

    # LOGNORMAL MODELS - intercept on LOG scale
    # Widened SD from 0.5 to 1 for better convergence
    "f0_std" = prior_string(
      paste0("student_t(3, ", fmt(log(vm$f0_std)), ", 1)"),
      class = "Intercept"
    ),

    "jitter" = prior_string(
      paste0("student_t(3, ", fmt(log(vm$jitter)), ", 1)"),
      class = "Intercept"
    ),

    "f2_std" = prior_string(
      paste0("student_t(3, ", fmt(log(vm$f2_std)), ", 1)"),
      class = "Intercept"
    ),

    stop("Outcome non riconosciuto")
  )

  # -----------------------------------
  # Fixed effects priors
  # Note: brms doesn't support prefix matching with coef=
  # So we use a single regularizing prior for all fixed effects
  # Adjust scale based on whether outcome uses lognormal
  # -----------------------------------
  is_lognormal <- outcome %in% c("f0_std", "f2_std", "jitter")

  if (is_lognormal) {
    # For lognormal: coefficients on log scale
    # Widened from 0.5 to 1 for better convergence
    priors_fixed <- prior_string("normal(0, 1)", class = "b")
  } else {
    # For gaussian/student: coefficients on original scale
    priors_fixed <- prior_string("normal(0, 10)", class = "b")
  }

  # -------------------------
  # Random effects + residual
  # Widened for better convergence
  # -------------------------
  priors_random <- c(
    prior_string("exponential(0.3)", class = "sd"),
    prior_string("exponential(0.3)", class = "sigma")
  )

  # -------------------------
  # Final prior set
  # -------------------------
  c(
    intercept_prior,
    priors_fixed,
    priors_random
  )
}


# --- Mappa outcome -> family -------------------------------------------------
# Family selection based on variable characteristics:
#   - gaussian(): f0_mean, nne (symmetric, no outliers)
#   - student(): f2_mean (symmetric but with outliers - robust)
#   - lognormal(): f0_std, f2_std, jitter (positive, right-skewed)
outcome_families <- list(
  f0_mean = gaussian(),
  f0_std = lognormal(),
  jitter = lognormal(),
  nne = gaussian(),
  f2_mean = student(),
  f2_std = lognormal()
)

# --- Lista degli outcome di interesse e vocali ------------------------------
outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
vowels <- c("a", "i", "u")

# --- Opzioni di sampling (usa i tuoi valori) -------------------------------
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
seed <- 123
# Aumentato adapt_delta per ridurre transizioni divergenti
control <- list(adapt_delta = 0.995, max_treedepth = 18)

# --- Cartella per salvare modelli (opzionale) ------------------------------
dir.create("models_original", showWarnings = FALSE)

# --- Flag: se TRUE lancia i modelli; se FALSE prepara tutto e stampa il piano -
run_models <- TRUE

# --- Loop per costruire formule e (opzionalmente) eseguire brm ---------------
fitted_models <- list()
for (v in vowels) {
  message("\n--- Vocale: /", v, "/ -------------------------------")
  for (out in outcomes) {
    # nome colonna effettivo nel df (es. f0_mean_a)
    colname <- paste0(out, "_", v)
    fam <- outcome_families[[out]]
    # safety: controlla che la col esista
    if (!colname %in% names(df_analysis)) {
      message("  Skipping ", colname, " (colonna non trovata nel dataset).")
      next
    }

    # formula: come nel tuo script, interaction con tutti i 5 domini e random slopes
    fmla <- bf(
      as.formula(
        paste0(
          colname,
          " ~ c1_stress * (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c + pid5_antagonism_bl_c + pid5_disinhibition_bl_c + pid5_psychoticism_bl_c) + ",
          "c2_recovery * (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c + pid5_antagonism_bl_c + pid5_disinhibition_bl_c + pid5_psychoticism_bl_c) + ",
          "(1 + c1_stress + c2_recovery | ID)"
        )
      )
    )

    priors_here <- make_priors_vowel(out, v)

    model_name <- paste0("m_", out, "_", v)

    # estrai nome della family in modo sicuro
    fam_name <- NULL
    if (!is.null(fam$family)) {
      fam_name <- fam$family
    } else if (!is.null(attr(fam, "family"))) {
      fam_name <- attr(fam, "family")
    } else {
      fam_name <- paste0(class(fam)[1])
    }

    message("  Model: ", model_name, "  (family: ", fam_name, ")")

    if (run_models) {
      # fitta il modello e salva l'oggetto
      fit <- brm(
        formula = fmla,
        data = df_analysis,
        family = fam,
        prior = priors_here,
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores,
        seed = seed,
        control = control,
        file = file.path("models_original", model_name) # salva su disco per riuso
      )
      fitted_models[[model_name]] <- fit

      # rapido check diagnostico (opzionale)
      print(summary(fit, waic = FALSE, odds = FALSE))
      pp_check(fit) # puoi commentare se vuoi meno output
    } else {
      # se non si esegue, mostriamo i priors che useremmo
      message("    (run_models = FALSE) - priors che verrebbero usati:")
      print(priors_here)
    }
  } # fine loop outcomes
} # fine loop vowels


# Ri-fitta m_jitter_u con più iterazioni
cat("Refitting m_jitter_u with more iterations...\n")
fitted_models[["m_jitter_u"]] <- brm(
  formula = bf(
    jitter_u ~
      c1_stress *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
        c2_recovery *
          (pid5_negative_affectivity_bl_c +
            pid5_detachment_bl_c +
            pid5_antagonism_bl_c +
            pid5_disinhibition_bl_c +
            pid5_psychoticism_bl_c) +
        (1 + c1_stress + c2_recovery | ID)
  ),
  data = df_analysis,
  family = lognormal(),
  prior = c(
    prior(student_t(3, -0.3, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(0.3), class = "sd"),
    prior(exponential(0.3), class = "sigma")
  ),
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.995, max_treedepth = 18),
  file = "models_original/m_jitter_u_refit"
)

# Ri-fitta m_nne_u con più iterazioni
cat("Refitting m_nne_u with more iterations...\n")
fitted_models[["m_nne_u"]] <- brm(
  formula = bf(
    nne_u ~
      c1_stress *
        (pid5_negative_affectivity_bl_c +
          pid5_detachment_bl_c +
          pid5_antagonism_bl_c +
          pid5_disinhibition_bl_c +
          pid5_psychoticism_bl_c) +
        c2_recovery *
          (pid5_negative_affectivity_bl_c +
            pid5_detachment_bl_c +
            pid5_antagonism_bl_c +
            pid5_disinhibition_bl_c +
            pid5_psychoticism_bl_c) +
        (1 + c1_stress + c2_recovery | ID)
  ),
  data = df_analysis,
  family = gaussian(),
  prior = c(
    prior(student_t(3, -28, 5), class = "Intercept"),
    prior(normal(0, 10), class = "b"),
    prior(exponential(0.3), class = "sd"),
    prior(exponential(0.3), class = "sigma")
  ),
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.995, max_treedepth = 18),
  file = "models_original/m_nne_u_refit"
)

# Verifica convergenza
cat("\nChecking convergence after refit...\n")
for (name in c("m_jitter_u", "m_nne_u")) {
  r <- brms::rhat(fitted_models[[name]])
  if (any(r > 1.05, na.rm = TRUE)) {
    cat("⚠️", name, "still has issues\n")
  } else {
    cat(
      "✓",
      name,
      "converged (max Rhat:",
      round(max(r, na.rm = TRUE), 3),
      ")\n"
    )
  }
}


# --- Uso degli oggetti fitted_models: se run_models = TRUE -> fitted_models popolato
# --- Fine

# --- 1. Estrarre summary da tutti i modelli -------------------------
results_table <- lapply(names(fitted_models), function(mn) {
  fit <- fitted_models[[mn]]
  fe <- fixef(fit, summary = TRUE) # posterior mean, sd, 95% CI

  # determina outcome e vocale dal nome del modello
  # es. m_f0_mean_a -> outcome = f0_mean, vowel = a
  parts <- str_split(mn, "_", simplify = TRUE)
  outcome <- paste(parts[2], parts[3], sep = "_")
  vowel <- parts[4]

  # identifica se il termine è main o interazione
  terms <- rownames(fe)
  type <- ifelse(str_detect(terms, ":"), "Interaction", "Main")

  tibble(
    outcome = outcome,
    vowel = vowel,
    parameter = terms,
    type = type,
    estimate = fe[, "Estimate"],
    se = fe[, "Est.Error"],
    ci_lower = fe[, "Q2.5"],
    ci_upper = fe[, "Q97.5"],
    significant = ifelse(fe[, "Q2.5"] > 0 | fe[, "Q97.5"] < 0, TRUE, FALSE)
  )
}) %>%
  bind_rows() %>%
  arrange(outcome, vowel, type, parameter)

# Visualizza tabella completa
cat("\n=== RESULTS TABLE ===\n")
print(data.frame(results_table), row.names = FALSE)

# --- 2. Estrai solo le interazioni significative -------------------------
sig_interactions <- results_table %>%
  filter(significant == TRUE, type == "Interaction")

cat("\n=== SIGNIFICANT INTERACTIONS (95% CI excludes zero) ===\n")
if (nrow(sig_interactions) > 0) {
  print(data.frame(sig_interactions), row.names = FALSE)
} else {
  cat("Nessuna interazione significativa trovata.\n")
}

# --- 3. Forest plot bayesiano -------------------------
p_forest <- ggplot(
  results_table %>% filter(type == "Interaction"),
  aes(
    x = estimate,
    y = parameter,
    xmin = ci_lower,
    xmax = ci_upper,
    color = outcome,
    shape = significant
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(outcome ~ vowel, scales = "free") +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
  theme_bw(base_size = 11) +
  labs(
    x = "Posterior Estimate",
    y = "Parametro",
    color = "Outcome",
    shape = "Significativo (CI esclude 0)",
    title = "Bayesian Forest Plot: Interazioni PID-5 × Stress",
    subtitle = "gaussian: scala originale | lognormal: scala log"
  ) +
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom"
  )

ggsave(
  "figures/forest_plot_interactions_ORIGINAL.png",
  p_forest,
  width = 14,
  height = 16,
  dpi = 300
)

cat(
  "\n✓ Forest plot salvato in figures/forest_plot_interactions_ORIGINAL.png\n"
)

# ==============================================================================
# 4. HELPER FUNCTIONS FOR LOGNORMAL COEFFICIENT INTERPRETATION
# ==============================================================================

#' Transform lognormal moderation coefficients to interpretable metrics
#'
#' For lognormal models, coefficients are on log scale. This function provides:
#' - Multiplicative effect (ratio): exp(beta)
#' - Percent change: (exp(beta) - 1) * 100
#' - Absolute change at geometric mean: geom_mean * (exp(beta) - 1)

interpret_lognormal_moderation <- function(model, outcome_var, data) {
  # Get geometric mean of outcome
  geom_mean <- exp(mean(log(data[[outcome_var]]), na.rm = TRUE))

  # Extract fixed effects
  fe <- fixef(model, summary = TRUE)

  # Create interpretation table
  tibble(
    parameter = rownames(fe),
    beta_log = fe[, "Estimate"],
    se_log = fe[, "Est.Error"],
    ci_lower_log = fe[, "Q2.5"],
    ci_upper_log = fe[, "Q97.5"],
    # Multiplicative effect (ratio)
    ratio = exp(beta_log),
    ratio_ci_lower = exp(ci_lower_log),
    ratio_ci_upper = exp(ci_upper_log),
    # Percent change
    pct_change = (exp(beta_log) - 1) * 100,
    pct_ci_lower = (exp(ci_lower_log) - 1) * 100,
    pct_ci_upper = (exp(ci_upper_log) - 1) * 100,
    # Absolute change at geometric mean
    abs_change_at_geomean = geom_mean * (exp(beta_log) - 1),
    abs_ci_lower = geom_mean * (exp(ci_lower_log) - 1),
    abs_ci_upper = geom_mean * (exp(ci_upper_log) - 1),
    # Reference
    geometric_mean = geom_mean,
    # Significance
    significant = (ci_lower_log > 0) | (ci_upper_log < 0)
  )
}

# --- 5. Create interpretable results for lognormal models -------------------------
cat("\n=== LOGNORMAL MODEL INTERPRETATIONS ===\n")
cat("Converting log-scale coefficients to original scale for:\n")
cat("  - f0_std (all vowels)\n")
cat("  - f2_std (all vowels)\n")
cat("  - jitter (all vowels)\n\n")

lognormal_outcomes <- c("f0_std", "f2_std", "jitter")
lognormal_interpretations <- list()

for (v in vowels) {
  for (out in lognormal_outcomes) {
    model_name <- paste0("m_", out, "_", v)
    col_name <- paste0(out, "_", v)

    if (model_name %in% names(fitted_models)) {
      cat("Processing:", model_name, "\n")
      interp <- interpret_lognormal_moderation(
        fitted_models[[model_name]],
        col_name,
        df_analysis
      )
      lognormal_interpretations[[model_name]] <- interp

      # Show significant interactions
      sig_interp <- interp %>%
        dplyr::filter(significant, grepl(":", parameter))

      if (nrow(sig_interp) > 0) {
        cat("  Significant interactions:\n")
        print(
          sig_interp %>%
            dplyr::select(parameter, pct_change, pct_ci_lower, pct_ci_upper) %>%
            mutate(across(where(is.numeric), ~ round(., 2)))
        )
      }
    }
  }
}

# --- 6. Salva risultati -------------------------
saveRDS(fitted_models, "results/fitted_models_ORIGINAL.rds")
saveRDS(results_table, "results/results_table_ORIGINAL.rds")
saveRDS(
  lognormal_interpretations,
  "results/lognormal_moderation_interpretations.rds"
)
save.image("results/moderation_analysis_ORIGINAL.RData")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("ANALISI COMPLETATA\n")
cat(rep("=", 70), "\n", sep = "")

# ==============================================================================
# NOTA INTERPRETATIVA
# ==============================================================================
#
# MODELLI GAUSSIAN (f0_mean, f2_mean, nne):
# - Coefficienti direttamente interpretabili nella scala originale
# - Es: β = 5.03 per f0_mean significa +5.03 Hz per +1 SD nel predittore
#
# MODELLI LOGNORMAL (f0_std, f2_std, jitter):
# - Coefficienti su scala LOG
# - Per interpretare:
#   * exp(β) = effetto moltiplicativo (ratio)
#   * (exp(β) - 1) * 100 = variazione percentuale
#   * geom_mean * (exp(β) - 1) = variazione assoluta
#
# Esempio: se β = 0.15 per un'interazione su jitter:
#   * exp(0.15) = 1.162 → 16.2% di aumento
#   * Se geom_mean(jitter) = 0.9% → aumento di ~0.15 punti percentuali
#
# Per le interazioni PID-5:
# - I coefficienti rappresentano come il tratto di personalità MODERA
#   l'effetto dello stress (c1) o del recovery (c2)
# - Un β positivo per c1_stress:pid5_X significa che livelli più alti
#   del tratto X amplificano la risposta allo stress
# ==============================================================================
