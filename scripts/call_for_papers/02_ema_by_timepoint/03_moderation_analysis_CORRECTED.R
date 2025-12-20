# ==============================================================================
# Bayesian Moderation Analysis: PID-5 Traits × Stress Reactivity
# CORRECTED VERSION - Coerente con abstract
# ==============================================================================
# CORREZIONI:
#   1. Usa i predittori PID-5 preparati dallo script 02_CORRECTED
#   2. Nomi variabili coerenti: pid5_*_c (centrati, between-person)
#   3. Contrast coding per stress (c1) e recovery (c2)
#   4. Random effects stabili (||, no correlazioni)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayesplot)
  library(tidybayes)
  library(ggdist)
  library(stringr)
})

options(brms.backend = "cmdstanr")
options(max.print = 5000)

dir.create("models", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("MODERATION ANALYSIS: PID-5 × Stress Reactivity (CORRECTED)\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

# ==============================================================================
# 0. LOAD DATA
# ==============================================================================

if (!exists("df_analysis")) {
  if (file.exists("results/df_analysis.rds")) {
    df_analysis <- readRDS("results/df_analysis.rds")
    message("Loaded df_analysis from results/df_analysis.rds")
  } else {
    stop("Esegui prima 02_voice_personality_analysis_CORRECTED.R")
  }
}

# Verifica che le variabili PID-5 centrate esistano
pid5_vars_c <- c(
  "pid5_negative_affectivity_c",
  "pid5_detachment_c",
  "pid5_antagonism_c",
  "pid5_disinhibition_c",
  "pid5_psychoticism_c"
)

if (!all(pid5_vars_c %in% names(df_analysis))) {
  stop("Variabili PID-5 centrate non trovate. Esegui 02_CORRECTED.R prima.")
}

df_analysis <- df_analysis %>%
  mutate(timepoint = factor(timepoint, levels = c("baseline", "pre", "post")))

# ==============================================================================
# 1. CONTRAST CODING
# ==============================================================================

# c1_stress: PRE vs BASELINE (effetto dello stress)
# c2_recovery: POST vs PRE (effetto del recupero)
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

# Ensure strictly positive for lognormal
logpositive_vars <- c("f0_std", "f2_std", "jitter")
vowels_all <- c("a", "i", "u")

for (var in logpositive_vars) {
  for (v in vowels_all) {
    col <- paste0(var, "_", v)
    if (!col %in% names(df_analysis)) next
    bad <- which(df_analysis[[col]] <= 0)
    if (length(bad) > 0) {
      min_pos <- min(df_analysis[[col]][df_analysis[[col]] > 0], na.rm = TRUE)
      df_analysis[[col]][bad] <- min_pos
    }
  }
}

# ==============================================================================
# 2. PRIORS E FAMILIES
# ==============================================================================

# Geometric means per prior informativi
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

# Prior builder
make_priors_vowel <- function(outcome, vowel, family_name = c("gaussian", "lognormal", "student")) {
  family_name <- match.arg(family_name)
  vm <- vowel_means[[vowel]]
  
  # Intercept prior
  intercept_prior <- switch(outcome,
    "f0_mean" = prior_string(paste0("student_t(3, ", fmt(vm$f0_mean), ", 30)"), class = "Intercept"),
    "nne" = prior_string(paste0("student_t(3, ", fmt(vm$nne), ", 5)"), class = "Intercept"),
    "f2_mean" = prior_string(paste0("student_t(3, ", fmt(vm$f2_mean), ", 150)"), class = "Intercept"),
    "f0_std" = prior_string(paste0("student_t(3, ", fmt(log(vm$f0_std)), ", 1)"), class = "Intercept"),
    "jitter" = prior_string(paste0("student_t(3, ", fmt(log(vm$jitter)), ", 1)"), class = "Intercept"),
    "f2_std" = prior_string(paste0("student_t(3, ", fmt(log(vm$f2_std)), ", 1)"), class = "Intercept")
  )
  
  # Fixed effects prior
  is_log_scale <- family_name %in% c("lognormal") && outcome %in% c("f0_std", "f2_std", "jitter")
  b_prior <- if (is_log_scale) {
    prior_string("normal(0, 0.5)", class = "b")  # Più shrinkage per log-scale
  } else {
    prior_string("normal(0, 10)", class = "b")
  }
  
  # Random effects
  re_prior <- prior_string("exponential(0.5)", class = "sd")
  
  # Residual
  resid_prior <- switch(family_name,
    "gaussian" = prior_string("exponential(0.3)", class = "sigma"),
    "student" = prior_string("exponential(0.3)", class = "sigma"),
    "lognormal" = prior_string("exponential(1)", class = "sigma")
  )
  
  # Nu for student
  nu_prior <- if (family_name == "student") {
    prior_string("gamma(2, 0.1)", class = "nu")
  } else NULL
  
  c(intercept_prior, b_prior, re_prior, resid_prior, nu_prior)
}

# Families
family_map <- list(
  f0_mean = gaussian(),
  f0_std = lognormal(),
  jitter = lognormal(),
  nne = gaussian(),
  f2_mean = student(),
  f2_std = lognormal()
)

# ==============================================================================
# 3. MODEL SPECIFICATIONS
# ==============================================================================

outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
vowels <- c("a", "i", "u")

# Sampling settings
iter <- 5000
warmup <- 2500
chains <- 4
cores <- 4
seed <- 123
control <- list(adapt_delta = 0.995, max_treedepth = 18)

# Formula con interazioni PID-5 × contrasti
# Nota: usiamo le variabili _c (centrate)
traits <- "pid5_negative_affectivity_c + pid5_detachment_c + pid5_antagonism_c + pid5_disinhibition_c + pid5_psychoticism_c"

build_formula <- function(colname) {
  as.formula(
    paste0(
      colname, " ~ c1_stress * (", traits, ") + c2_recovery * (", traits, ") + (1 + c1_stress + c2_recovery || ID)"
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
    
    fam <- family_map[[out]]
    fam_key <- if (out == "f2_mean") "student" else if (out %in% c("f0_std", "f2_std", "jitter")) "lognormal" else "gaussian"
    
    fmla <- bf(build_formula(colname))
    priors_here <- make_priors_vowel(out, v, family_name = fam_key)
    
    model_name <- paste0("m_", out, "_", v)
    message("  Fitting: ", model_name)
    
    fit <- brm(
      formula = fmla,
      data = df_analysis,
      family = fam,
      prior = priors_here,
      iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
      control = control,
      file = file.path("models", model_name)
    )
    
    fitted_models[[model_name]] <- fit
  }
}

# ==============================================================================
# 5. CONVERGENCE CHECK
# ==============================================================================

cat("\n=== Convergence check ===\n")
for (mn in names(fitted_models)) {
  r <- brms::rhat(fitted_models[[mn]])
  max_r <- max(r, na.rm = TRUE)
  status <- if (max_r < 1.01) "✓" else "⚠"
  cat(status, " ", mn, " max Rhat: ", round(max_r, 3), "\n", sep = "")
}

# ==============================================================================
# 6. RESULTS TABLE
# ==============================================================================

.parse_model_name <- function(mn) {
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
  
  tibble(
    model = mn,
    outcome = info$outcome,
    vowel = info$vowel,
    parameter = rownames(fe),
    type = ifelse(str_detect(rownames(fe), ":"), "Interaction", "Main"),
    estimate = fe[, "Estimate"],
    se = fe[, "Est.Error"],
    ci_lower = fe[, "Q2.5"],
    ci_upper = fe[, "Q97.5"],
    significant = (fe[, "Q2.5"] > 0) | (fe[, "Q97.5"] < 0)
  )
}) %>%
  bind_rows() %>%
  arrange(outcome, vowel, type, parameter)

# ==============================================================================
# 7. KEY RESULTS - Matching Abstract
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("KEY RESULTS - Verifica coerenza con abstract\n")
cat(rep("=", 70), "\n", sep = "")

cat("\n--- Abstract claims to verify ---\n")
cat("1. Negative Affectivity amplifies stress-induced pitch increases (β=5.03 Hz)\n")
cat("2. Detachment impairs post-stress recovery (β=-4.37)\n")
cat("3. Antagonism facilitates recovery (β=3.44-4.17 across vowels)\n")
cat("4. Psychoticism shows pervasive pitch elevation and voice quality degradation\n")
cat("5. Disinhibition shows F2 alterations (β=63 Hz)\n\n")

# 1. NA × stress on F0
na_stress_f0 <- results_table %>%
  filter(str_detect(parameter, "c1_stress:pid5_negative_affectivity"),
         outcome == "f0_mean")
cat("\n[1] Negative Affectivity × Stress → F0 Mean:\n")
print(na_stress_f0 %>% select(vowel, estimate, ci_lower, ci_upper, significant))

# 2. Detachment × recovery
det_recovery_f0 <- results_table %>%
  filter(str_detect(parameter, "c2_recovery:pid5_detachment"),
         outcome == "f0_mean")
cat("\n[2] Detachment × Recovery → F0 Mean:\n")
print(det_recovery_f0 %>% select(vowel, estimate, ci_lower, ci_upper, significant))

# 3. Antagonism × recovery
ant_recovery <- results_table %>%
  filter(str_detect(parameter, "c2_recovery:pid5_antagonism"))
cat("\n[3] Antagonism × Recovery (all outcomes):\n")
print(ant_recovery %>% filter(significant) %>% select(outcome, vowel, estimate, ci_lower, ci_upper))

# 4. Psychoticism main effects
psy_effects <- results_table %>%
  filter(str_detect(parameter, "pid5_psychoticism"),
         type == "Main" | str_detect(parameter, "c1_stress"))
cat("\n[4] Psychoticism effects:\n")
print(psy_effects %>% filter(significant) %>% select(outcome, vowel, parameter, estimate, ci_lower, ci_upper))

# 5. Disinhibition × F2
dis_f2 <- results_table %>%
  filter(str_detect(parameter, "pid5_disinhibition"),
         outcome == "f2_mean")
cat("\n[5] Disinhibition × F2 Mean:\n")
print(dis_f2 %>% select(vowel, parameter, estimate, ci_lower, ci_upper, significant))

# ==============================================================================
# 8. SIGNIFICANT INTERACTIONS SUMMARY
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("ALL SIGNIFICANT INTERACTIONS (95% CI excludes zero)\n")
cat(rep("=", 70), "\n", sep = "")

sig_interactions <- results_table %>%
  filter(significant, type == "Interaction") %>%
  arrange(outcome, vowel, parameter)

if (nrow(sig_interactions) > 0) {
  print(sig_interactions %>% select(outcome, vowel, parameter, estimate, ci_lower, ci_upper), n = 100)
} else {
  cat("Nessuna interazione con 95% CI che esclude zero.\n")
  cat("\nInterazioni con 90% CI che esclude zero:\n")
  
  # Ricalcola con 90% CI
  sig_90 <- lapply(names(fitted_models), function(mn) {
    fit <- fitted_models[[mn]]
    post <- as_draws_df(fit)
    
    info <- .parse_model_name(mn)
    
    # Estrai solo interazioni
    int_vars <- grep("^b_c[12].*:pid5", names(post), value = TRUE)
    
    lapply(int_vars, function(v) {
      samples <- post[[v]]
      tibble(
        model = mn,
        outcome = info$outcome,
        vowel = info$vowel,
        parameter = sub("^b_", "", v),
        estimate = mean(samples),
        ci90_lower = quantile(samples, 0.05),
        ci90_upper = quantile(samples, 0.95),
        pd = mean(samples > 0)
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  sig_90_filtered <- sig_90 %>%
    mutate(significant_90 = (ci90_lower > 0) | (ci90_upper < 0)) %>%
    filter(significant_90)
  
  if (nrow(sig_90_filtered) > 0) {
    print(sig_90_filtered, n = 50)
  }
}

# ==============================================================================
# 9. FOREST PLOT
# ==============================================================================

p_forest <- ggplot(
  results_table %>% filter(type == "Interaction"),
  aes(x = estimate, y = parameter,
      xmin = ci_lower, xmax = ci_upper,
      color = outcome, shape = significant)
) +
  geom_point(size = 2.5) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(outcome ~ vowel, scales = "free") +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
  theme_bw(base_size = 10) +
  labs(
    x = "Posterior Estimate",
    y = "Parameter",
    title = "PID-5 × Stress Interactions",
    subtitle = "95% CI | Filled = significant"
  ) +
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom"
  )

ggsave("figures/forest_plot_interactions.png", p_forest, width = 14, height = 16, dpi = 300)

# ==============================================================================
# 10. SAVE
# ==============================================================================

saveRDS(fitted_models, "results/fitted_models.rds")
saveRDS(results_table, "results/results_table.rds")
save.image("results/moderation_analysis.RData")

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")
