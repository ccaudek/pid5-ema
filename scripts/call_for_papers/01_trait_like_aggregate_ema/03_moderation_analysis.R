# ==============================================================================
# Bayesian Moderation Analysis: PID-5 Traits × Stress Reactivity
# Robust/stable version with special handling for jitter models
# ==============================================================================
# PRINCIPALI MODIFICHE vs versione precedente:
# 1) Random effects: usa "||" (no correlazioni) per stabilità generale
# 2) JITTER: priors più regolarizzanti + random effects più parsimoniosi
# 3) Opzionale fallback per JITTER: Gamma(link="log") se lognormal continua a non convergere
# 4) CORRETTO: rimosso prior LKJ quando si usa || (no correlazioni)
# 5) CORRETTO: gestione appropriata di sigma/shape per diverse famiglie
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayesplot)
  library(tidybayes)
  library(ggdist)
  library(tidyr)
  library(stringr)
})

options(brms.backend = "cmdstanr")
options(max.print = 5000)

dir.create("models", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("MODERATION ANALYSIS: PID-5 × Stress Reactivity (ROBUST)\n")
cat("Special handling for JITTER models\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

# ==============================================================================
# 0. LOAD DATA
# ==============================================================================
if (!exists("df_analysis")) {
  if (file.exists("results/df_analysis.rds")) {
    df_analysis <- readRDS("results/df_analysis.rds")
    message("Loaded df_analysis from results/df_analysis.rds")
  } else {
    stop("df_analysis non trovato e results/df_analysis.rds non esiste.")
  }
}

df_analysis <- df_analysis %>%
  mutate(timepoint = factor(timepoint, levels = c("baseline", "pre", "post")))

# ==============================================================================
# 1. CONTRAST CODING + LOGNORMAL CLEANING
# ==============================================================================
df_analysis <- df_analysis %>%
  mutate(
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre" ~ 0.5,
      timepoint == "post" ~ 0
    ),
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre" ~ -0.5,
      timepoint == "post" ~ 0.5
    )
  )

# Ensure strictly positive for lognormal/gamma
logpositive_vars <- c("f0_std", "f2_std", "jitter")
vowels_all <- c("a", "i", "u")

cat(
  "\n=== Checking for <=0 values in positive outcomes (lognormal/gamma) ===\n"
)
for (var in logpositive_vars) {
  for (v in vowels_all) {
    col <- paste0(var, "_", v)
    if (!col %in% names(df_analysis)) next

    bad <- which(df_analysis[[col]] <= 0)
    if (length(bad) > 0) {
      min_pos <- min(df_analysis[[col]][df_analysis[[col]] > 0], na.rm = TRUE)
      df_analysis[[col]][bad] <- min_pos
      cat(
        "Replaced",
        length(bad),
        "values <=0 in",
        col,
        "with",
        signif(min_pos, 6),
        "\n"
      )
    }
  }
}

# ==============================================================================
# 2. PRIORS (DATA-DRIVEN) + FAMILIES
# ==============================================================================
compute_vowel_means <- function(data) {
  outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
  vowels <- c("a", "i", "u")
  pos_outcomes <- c("f0_std", "f2_std", "jitter")

  out <- list()
  for (v in vowels) {
    out[[v]] <- list()
    for (o in outcomes) {
      nm <- paste0(o, "_", v)
      if (!nm %in% names(data)) next
      x <- data[[nm]]
      x <- x[!is.na(x)]
      if (o %in% pos_outcomes) {
        x <- x[x > 0]
        out[[v]][[o]] <- exp(mean(log(x)))
      } else {
        out[[v]][[o]] <- mean(x)
      }
    }
  }
  out
}

vowel_means <- compute_vowel_means(df_analysis)

fmt <- function(x, digits = 6) formatC(x, format = "f", digits = digits)

# ------------------------------------------------------------------------------
# Base priors builder (gaussian/lognormal/student/gamma)
# CORRETTO: niente prior LKJ quando use_correlations = FALSE
# CORRETTO: gestione appropriata di sigma/shape per diverse famiglie
# ------------------------------------------------------------------------------
make_priors_vowel <- function(
  outcome,
  vowel,
  family_name = c("gaussian", "lognormal", "student", "gamma"),
  use_correlations = FALSE
) {
  family_name <- match.arg(family_name)
  vm <- vowel_means[[vowel]]

  # --- Intercept prior ---
  intercept_prior <- switch(
    outcome,
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

  # --- Fixed effects prior ---
  is_pos <- outcome %in% c("f0_std", "f2_std", "jitter")
  is_log_scale <- family_name %in% c("lognormal", "gamma") && is_pos

  b_prior <- if (outcome == "jitter") {
    # jitter: più shrinkage per stabilità
    prior_string("normal(0, 0.5)", class = "b")
  } else if (is_log_scale) {
    prior_string("normal(0, 1)", class = "b")
  } else {
    prior_string("normal(0, 10)", class = "b")
  }

  # --- Random effects prior ---
  # Solo sd, NO cor quando si usa || (use_correlations = FALSE)
  sd_rate <- if (outcome == "jitter") 1 else 0.3
  re_priors <- prior_string(paste0("exponential(", sd_rate, ")"), class = "sd")

  # Aggiungi prior LKJ solo se il modello usa correlazioni (|)
  if (use_correlations) {
    re_priors <- c(re_priors, prior_string("lkj(2)", class = "cor"))
  }

  # --- Residual/dispersion prior ---
  # - gaussian e student: hanno sigma

  # - lognormal: ha sigma (sd della distribuzione normale sottostante)
  # - gamma: ha shape (non sigma)
  resid_prior <- switch(
    family_name,
    "gaussian" = prior_string("exponential(0.3)", class = "sigma"),
    "student" = prior_string("exponential(0.3)", class = "sigma"),
    "lognormal" = prior_string("exponential(1)", class = "sigma"),
    "gamma" = prior_string("gamma(0.01, 0.01)", class = "shape"),
    NULL
  )

  # --- Student nu prior ---
  nu_prior <- if (family_name == "student") {
    prior_string("gamma(2, 0.1)", class = "nu")
  } else {
    NULL
  }

  # Combina tutti i priors (rimuove NULL automaticamente)
  c(intercept_prior, b_prior, re_priors, resid_prior, nu_prior)
}

# Families (coerenti con main script; jitter può essere forzato a Gamma come fallback)
family_map <- list(
  f0_mean = gaussian(),
  f0_std = lognormal(),
  jitter = lognormal(), # <-- default; puoi cambiare a Gamma(link="log") sotto
  nne = gaussian(),
  f2_mean = student(),
  f2_std = lognormal()
)

# ==============================================================================
# 3. MODEL SPECIFICATIONS (con gestione speciale per jitter)
# ==============================================================================
outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
vowels <- c("a", "i", "u")

# Sampling options
iter <- 5000
warmup <- 2500
chains <- 4
cores <- 4
seed <- 123
control <- list(adapt_delta = 0.995, max_treedepth = 18)

# Se TRUE: per jitter usa Gamma(log) invece di lognormal
use_gamma_for_jitter <- FALSE

# Random effects strategy:
# Usa || (no correlazioni) per maggiore stabilità
# Se ancora non converge per jitter, prova strutture più parsimoniose
jitter_re_structure <- "(1 + c1_stress + c2_recovery || ID)"
other_re_structure <- "(1 + c1_stress + c2_recovery || ID)"

# Build formula strings
traits <- "pid5_negative_affectivity_bl_c + pid5_detachment_bl_c + pid5_antagonism_bl_c + pid5_disinhibition_bl_c + pid5_psychoticism_bl_c"

build_formula <- function(colname, re_struct) {
  as.formula(
    paste0(
      colname,
      " ~ c1_stress * (",
      traits,
      ") + ",
      "c2_recovery * (",
      traits,
      ") + ",
      re_struct
    )
  )
}

# ==============================================================================
# 4. FIT MODELS
# ==============================================================================
fitted_models <- list()

for (v in vowels) {
  message("\n--- Vocale: /", v, "/ -------------------------------")
  for (out in outcomes) {
    colname <- paste0(out, "_", v)
    if (!colname %in% names(df_analysis)) {
      message("  Skipping ", colname, " (colonna non trovata).")
      next
    }

    # Select family
    fam <- family_map[[out]]
    fam_name <- tryCatch(fam$family, error = function(e) class(fam)[1])

    if (out == "jitter" && use_gamma_for_jitter) {
      fam <- Gamma(link = "log")
      fam_name <- "Gamma(log)"
    }

    # Random effects structure
    re_struct <- if (out == "jitter") jitter_re_structure else
      other_re_structure

    # Formula
    fmla <- bf(build_formula(colname, re_struct))

    # Priors (family-aware)
    # Determina la chiave per la famiglia
    fam_key <- if (out == "f2_mean") {
      "student"
    } else if (out %in% c("f0_std", "f2_std", "jitter")) {
      if (out == "jitter" && use_gamma_for_jitter) "gamma" else "lognormal"
    } else {
      "gaussian"
    }

    # use_correlations = FALSE perché usiamo ||
    priors_here <- make_priors_vowel(
      out,
      v,
      family_name = fam_key,
      use_correlations = FALSE
    )

    model_name <- paste0("m_", out, "_", v)
    message(
      "  Model: ",
      model_name,
      "  (family: ",
      fam_name,
      ", RE: ",
      re_struct,
      ")"
    )

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
      file = file.path("models", model_name)
    )

    fitted_models[[model_name]] <- fit
    print(summary(fit, waic = FALSE, loo = FALSE))
  }
}

# ==============================================================================
# 5. CONVERGENCE CHECK
# ==============================================================================
cat("\n=== Convergence check (max Rhat per modello) ===\n")
for (mn in names(fitted_models)) {
  r <- brms::rhat(fitted_models[[mn]])
  cat(mn, " max Rhat: ", round(max(r, na.rm = TRUE), 3), "\n", sep = "")
}

# ==============================================================================
# 6. RESULTS TABLE (fixed effects)
# ==============================================================================
.parse_model_name <- function(mn) {
  stopifnot(startsWith(mn, "m_"))
  rest <- sub("^m_", "", mn)
  parts <- str_split(rest, "_", simplify = TRUE)
  vowel <- parts[ncol(parts)]
  outcome <- paste(parts[1:(ncol(parts) - 1)], collapse = "_")
  list(outcome = outcome, vowel = vowel)
}

results_table <- lapply(names(fitted_models), function(mn) {
  fit <- fitted_models[[mn]]
  fe <- fixef(fit, summary = TRUE)

  info <- .parse_model_name(mn)
  outcome <- info$outcome
  vowel <- info$vowel

  terms <- rownames(fe)
  type <- ifelse(str_detect(terms, ":"), "Interaction", "Main")

  tibble(
    model = mn,
    outcome = outcome,
    vowel = vowel,
    parameter = terms,
    type = type,
    estimate = fe[, "Estimate"],
    se = fe[, "Est.Error"],
    ci_lower = fe[, "Q2.5"],
    ci_upper = fe[, "Q97.5"],
    significant = (fe[, "Q2.5"] > 0) | (fe[, "Q97.5"] < 0)
  )
}) %>%
  bind_rows() %>%
  arrange(outcome, vowel, type, parameter)

cat("\n=== RESULTS TABLE (fixed effects) ===\n")
print(results_table, n = min(nrow(results_table), 200))

sig_interactions <- results_table %>%
  filter(significant, type == "Interaction")

cat("\n=== SIGNIFICANT INTERACTIONS (95% CI excludes zero) ===\n")
if (nrow(sig_interactions) > 0) {
  print(sig_interactions, n = min(nrow(sig_interactions), 200))
} else {
  cat("Nessuna interazione significativa trovata.\n")
}

# ==============================================================================
# 7. FOREST PLOT (Interactions)
# ==============================================================================
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
  geom_point(size = 2.7) +
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
    subtitle = "gaussian/student: scala originale | lognormal/gamma: coefficienti su scala log"
  ) +
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom"
  )

ggsave(
  "figures/forest_plot_interactions.png",
  p_forest,
  width = 14,
  height = 16,
  dpi = 300
)

cat("\n✓ Forest plot salvato in figures/forest_plot_interactions.png\n")

# ==============================================================================
# 8. LOGNORMAL INTERPRETATION (jitter/f0_std/f2_std)
# ==============================================================================
interpret_lognormal_moderation <- function(model, outcome_var, data) {
  geom_mean <- exp(mean(log(data[[outcome_var]]), na.rm = TRUE))
  fe <- fixef(model, summary = TRUE)

  tibble(
    parameter = rownames(fe),
    beta_log = fe[, "Estimate"],
    se_log = fe[, "Est.Error"],
    ci_lower_log = fe[, "Q2.5"],
    ci_upper_log = fe[, "Q97.5"],
    ratio = exp(beta_log),
    ratio_ci_lower = exp(ci_lower_log),
    ratio_ci_upper = exp(ci_upper_log),
    pct_change = (exp(beta_log) - 1) * 100,
    pct_ci_lower = (exp(ci_lower_log) - 1) * 100,
    pct_ci_upper = (exp(ci_upper_log) - 1) * 100,
    geometric_mean = geom_mean,
    significant = (ci_lower_log > 0) | (ci_upper_log < 0)
  )
}

cat("\n=== LOGNORMAL/GAMMA COEFFICIENT INTERPRETATIONS (percent change) ===\n")
pos_outcomes <- c("f0_std", "f2_std", "jitter")
lognormal_interpretations <- list()

for (v in vowels) {
  for (out in pos_outcomes) {
    mn <- paste0("m_", out, "_", v)
    col <- paste0(out, "_", v)

    if (mn %in% names(fitted_models)) {
      message("Processing: ", mn)
      lognormal_interpretations[[mn]] <- interpret_lognormal_moderation(
        fitted_models[[mn]],
        col,
        df_analysis
      )
    }
  }
}

# ==============================================================================
# 9. SAVE
# ==============================================================================
saveRDS(fitted_models, "results/fitted_models.rds")
saveRDS(results_table, "results/results_table.rds")
saveRDS(
  lognormal_interpretations,
  "results/lognormal_moderation_interpretations.rds"
)
save.image("results/moderation_analysis.RData")

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("ANALISI COMPLETATA (ROBUST)\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
