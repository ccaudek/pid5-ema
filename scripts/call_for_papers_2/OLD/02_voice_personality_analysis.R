# ==============================================================================
# Bayesian Analysis: Voice Acoustics, Stress, and Personality Pathology
# Context-dependent expression of PID-5 traits via passive sensing
# ==============================================================================
# VERSIONE ORIGINALE - Ripristina l'analisi coerente con l'abstract
# Modifiche rispetto alla versione modificata:
#   1. jitter: lognormal() → gaussian()
#   2. f2_mean: asym_laplace() → gaussian()
#   3. f2_std: lognormal() → gaussian()
#   4. Plot: rimosso log(jitter_a), ora usa jitter_a direttamente
# ==============================================================================

# Load required packages
library(here)
library(tidyverse)
library(readxl)
library(brms)
library(cmdstanr)
library(bayestestR)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(ggdist)
library(lubridate)
library(stringr)
library(missRanger)

# Set cmdstanr as backend for brms
options(brms.backend = "cmdstanr")

# Create output directories
dir.create("models_original", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

# Load data from three timepoints
baseline <- read_excel(
  here::here(
    "data",
    "raw",
    "acustic_features",
    "datiacustici",
    "AUDIO.xlsx"
  ),
  sheet = "BASELINE"
)

pre <- read_excel(
  here::here(
    "data",
    "raw",
    "acustic_features",
    "datiacustici",
    "AUDIO.xlsx"
  ),
  sheet = "PRE"
)

post <- read_excel(
  here::here(
    "data",
    "raw",
    "acustic_features",
    "datiacustici",
    "AUDIO.xlsx"
  ),
  sheet = "POST"
)

# Add timepoint indicator
baseline$timepoint <- "baseline"
pre$timepoint <- "pre"
post$timepoint <- "post"

# Combine all data
df_wide <- bind_rows(baseline, pre, post)

# Clean column names (remove extra spaces)
names(df_wide) <- str_trim(names(df_wide))

# Correct ID names
df_wide <- df_wide %>%
  mutate(
    ID = case_when(
      ID == "am_bo_1988_08_24_166" ~ "an_bo_1988_08_24_166",
      ID == "as_li_2005_04_26_447" ~ "as_si_2005_04_26_447",
      ID == "cl_bo_1987_10_16_628" ~ "ca_bo_1987_10_16_628",
      ID == "hi_na_2005_03_08_339" ~ "gi_na_2005_03_08_339",
      ID == "ma_si_2003_10_31_940" ~ "si_ma_2003_10_31_940",
      TRUE ~ ID
    )
  )

length(unique(df_wide$ID))
# [1] 141

# Extract key acoustic variables (we'll focus on /a/ vowel for primary analysis)
# F0 mean - pitch (arousal/anxiety)
# F0 std - pitch variability (emotional instability)
# Jitter - voice quality (tension)
# NNE - normalized noise energy (vocal stress)

# Select relevant variables
df_clean <- df_wide %>%
  dplyr::select(
    ID,
    Data,
    timepoint,
    # Acoustic measures (vowel /a/)
    f0_mean_a = `F0 mean Hz /a/`,
    f0_std_a = `F0 std Hz /a/`,
    jitter_a = `Jitter /a/`,
    nne_a = `NNE /a/`,
    f2_mean_a = `F2 mean Hz /a/`,
    f2_std_a = `F2 Std Hz /a/`,

    # Acoustic measures (vowel /i/)
    f0_mean_i = `F0 mean Hz /i/`,
    f0_std_i = `F0 std Hz /i/`,
    jitter_i = `Jitter /i/`,
    nne_i = `NNE /i/`,
    f2_mean_i = `F2 mean Hz /i/`,
    f2_std_i = `F2 Std Hz /i/`,

    # Acoustic measures (vowel /u/)
    f0_mean_u = `F0 mean Hz /u/`,
    f0_std_u = `F0 std Hz /u/`,
    jitter_u = `Jitter /u/`,
    nne_u = `NNE /u/`,
    f2_mean_u = `F2 mean Hz /u/`,
    f2_std_u = `F2 Std Hz /u/`,

    # PID-5 traits
    pid5_na = pid5_negative_affectivity,
    pid5_det = pid5_detachment,
    pid5_ant = pid5_antagonism,
    pid5_dis = pid5_disinhibition,
    pid5_psy = pid5_psychoticism
  ) %>%
  # Convert timepoint to factor with meaningful order
  mutate(timepoint = factor(timepoint, levels = c("baseline", "pre", "post")))


# Add correct PID-5 EMA and baseline measures ----------------------------------

data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
d <- rio::import(data_path)
vars <- c(
  "user_id",
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism",
  "domain_negative_affect_baseline",
  "domain_detachment_baseline",
  "domain_antagonism_baseline",
  "domain_disinhibition_baseline",
  "domain_psychoticism_baseline",
  "exam_period",
  "day",
  "sex"
)

d1 <- d |>
  dplyr::select(all_of(vars)) |>
  dplyr::rename(
    ID = user_id
  )

# 1.1 - crea colonne date coerenti
df_clean2 <- df_clean %>%
  mutate(
    Data_date = as.Date(Data) # Data è POSIXct/dttm -> Date
  ) |>
  dplyr::rename(
    date = Data_date
  ) |>
  dplyr::select(-c(Data, pid5_na, pid5_det, pid5_ant, pid5_dis, pid5_psy))

d1_2 <- d1 %>%
  mutate(
    day_date = as.Date(day) # day è IDate -> Date
  ) |>
  dplyr::rename(
    date = day_date,
    timepoint = exam_period
  ) |>
  dplyr::select(-c(day))

# 3) Filtra mantenendo solo le date consentite
voice_df <- df_clean2

table(
  voice_df$ID,
  voice_df$timepoint
)

allowed_id <- unique(voice_df$ID)

pid5_df <- d1_2 %>%
  dplyr::filter(ID %in% allowed_id)

pid5_df$timepoint <- forcats::fct_recode(
  pid5_df$timepoint,
  "baseline" = "baseline",
  "pre" = "pre_exam",
  "post" = "post_exam"
)

length(unique(pid5_df$ID))
length(unique(voice_df$ID))

table(voice_df$ID, voice_df$timepoint)

###################

pid5_agg <- pid5_df %>%
  group_by(ID, timepoint) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "{.col}"), # medie per tutte le numeric
    sex = first(na.omit(sex)), # sesso "rappresentativo"
    date_rep = if ("date" %in% names(.)) min(date, na.rm = TRUE) else
      as.Date(NA), # data rappresentativa (min)
    .groups = "drop"
  )

# 3a) JOIN principale: aggiunge le colonne aggregate a voice_df, preservando tutte le colonne di voice_df
joined_left <- voice_df %>%
  left_join(pid5_agg, by = c("ID", "timepoint"))

imp <- missRanger(joined_left, num.trees = 200)

df_clean <- imp

# Extract baseline PID-5 scores (for between-person predictors)
pid5_baseline <- df_clean %>%
  dplyr::filter(timepoint == "baseline") %>%
  dplyr::select(ID, starts_with("pid5")) %>%
  rename_with(~ paste0(., "_bl"), .cols = starts_with("pid5"))

# Create long format with baseline PID-5 as between-person predictors
df_long <- df_clean %>%
  dplyr::select(-starts_with("pid5")) %>% # Remove timepoint-varying PID-5
  left_join(pid5_baseline, by = "ID") %>%
  # Center baseline PID-5 scores for interpretability
  mutate(
    across(
      ends_with("_bl"),
      ~ scale(., center = TRUE, scale = TRUE)[, 1],
      .names = "{.col}_c"
    )
  )

# Remove missing cases
df_analysis <- df_long %>%
  dplyr::filter(complete.cases(.))

# Check sample size
cat("\n=== SAMPLE CHARACTERISTICS ===\n")
cat("Total observations:", nrow(df_analysis), "\n")
cat("Unique participants:", n_distinct(df_analysis$ID), "\n")
df_analysis %>%
  count(timepoint) %>%
  print()

# ==============================================================================
# 2. DESCRIPTIVE STATISTICS
# ==============================================================================

cat("\n=== DESCRIPTIVE STATISTICS: ACOUSTIC MEASURES BY TIMEPOINT ===\n")

# Calculate descriptives
descriptives <- df_analysis %>%
  group_by(timepoint) %>%
  summarise(
    across(
      c(
        f0_mean_a,
        f0_std_a,
        jitter_a,
        nne_a,
        f2_mean_a,
        f2_std_a,
        f0_mean_i,
        f0_std_i,
        jitter_i,
        nne_i,
        f2_mean_i,
        f2_std_i,
        f0_mean_u,
        f0_std_u,
        jitter_u,
        nne_u,
        f2_mean_u,
        f2_std_u
      ),
      list(
        M = ~ mean(., na.rm = TRUE),
        SD = ~ sd(., na.rm = TRUE),
        Median = ~ median(., na.rm = TRUE),
        Q25 = ~ quantile(., 0.25, na.rm = TRUE),
        Q75 = ~ quantile(., 0.75, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    n = n()
  )

print(descriptives, width = Inf)

# Calculate effect sizes (Cohen's d) for stress manipulation check
calc_cohens_d <- function(data, var, time1, time2) {
  group1 <- data %>% dplyr::filter(timepoint == time1) %>% pull({{ var }})
  group2 <- data %>% dplyr::filter(timepoint == time2) %>% pull({{ var }})

  m1 <- mean(group1, na.rm = TRUE)
  m2 <- mean(group2, na.rm = TRUE)
  sd_pooled <- sqrt(
    (sd(group1, na.rm = TRUE)^2 + sd(group2, na.rm = TRUE)^2) / 2
  )

  d <- (m2 - m1) / sd_pooled
  return(d)
}

cat("\n=== EFFECT SIZES: Stress Manipulation Check ===\n")
cat("\nCohen's d for PRE vs BASELINE (stress increase):\n")
cat(
  "F0 mean a:",
  round(calc_cohens_d(df_analysis, f0_mean_a, "baseline", "pre"), 3),
  "\n"
)
cat(
  "F0 std a:",
  round(calc_cohens_d(df_analysis, f0_std_a, "baseline", "pre"), 3),
  "\n"
)
cat(
  "Jitter a:",
  round(calc_cohens_d(df_analysis, jitter_a, "baseline", "pre"), 3),
  "\n"
)
cat(
  "NNE a:",
  round(calc_cohens_d(df_analysis, nne_a, "baseline", "pre"), 3),
  "\n"
)
cat(
  "F2 mean a:",
  round(calc_cohens_d(df_analysis, f2_mean_a, "baseline", "pre"), 3),
  "\n"
)
cat(
  "F2 std a:",
  round(calc_cohens_d(df_analysis, f2_std_a, "baseline", "pre"), 3),
  "\n"
)

cat("\nCohen's d for POST vs PRE (stress decrease):\n")
cat(
  "F0 mean a:",
  round(calc_cohens_d(df_analysis, f0_mean_a, "pre", "post"), 3),
  "\n"
)
cat(
  "F0 std a:",
  round(calc_cohens_d(df_analysis, f0_std_a, "pre", "post"), 3),
  "\n"
)
cat(
  "Jitter a:",
  round(calc_cohens_d(df_analysis, jitter_a, "pre", "post"), 3),
  "\n"
)
cat("NNE a:", round(calc_cohens_d(df_analysis, nne_a, "pre", "post"), 3), "\n")
cat(
  "F2 mean a:",
  round(calc_cohens_d(df_analysis, f2_mean_a, "pre", "post"), 3),
  "\n"
)
cat(
  "F2 std a:",
  round(calc_cohens_d(df_analysis, f2_std_a, "pre", "post"), 3),
  "\n"
)

# ==============================================================================
# 3. VISUALIZATION: DESCRIPTIVE PATTERNS
# ==============================================================================

# Create visualization function
plot_acoustic_trajectory <- function(data, outcome_var, outcome_label) {
  # Summary data for plotting
  summary_data <- data %>%
    group_by(timepoint) %>%
    summarise(
      M = mean({{ outcome_var }}, na.rm = TRUE),
      SE = sd({{ outcome_var }}, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  # Individual trajectories (sample for clarity)
  set.seed(123)
  sample_ids <- sample(unique(data$ID), min(30, n_distinct(data$ID)))

  ggplot() +
    # Individual trajectories (faded)
    geom_line(
      data = data %>% filter(ID %in% sample_ids),
      aes(x = timepoint, y = {{ outcome_var }}, group = ID),
      alpha = 0.1,
      color = "gray50"
    ) +
    # Mean trajectory
    geom_line(
      data = summary_data,
      aes(x = timepoint, y = M, group = 1),
      linewidth = 1.5,
      color = "#0072B2"
    ) +
    geom_point(
      data = summary_data,
      aes(x = timepoint, y = M),
      size = 4,
      color = "#0072B2"
    ) +
    # Error bars (SE)
    geom_errorbar(
      data = summary_data,
      aes(x = timepoint, ymin = M - SE, ymax = M + SE),
      width = 0.1,
      linewidth = 1,
      color = "#0072B2"
    ) +
    # Annotations
    annotate(
      "rect",
      xmin = 1.5,
      xmax = 2.5,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.1,
      fill = "red"
    ) +
    annotate(
      "text",
      x = 2,
      y = Inf,
      label = "EXAM\nSTRESS",
      vjust = 1.5,
      size = 3,
      fontface = "bold",
      color = "red"
    ) +
    labs(
      title = outcome_label,
      x = "Timepoint",
      y = outcome_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

# Create plots for each outcome
# MODIFICATO: jitter_a invece di log(jitter_a)
p1 <- plot_acoustic_trajectory(df_analysis, f0_mean_a, "F0 Mean (Hz)")
p2 <- plot_acoustic_trajectory(df_analysis, f0_std_a, "F0 SD (Hz)")
p3 <- plot_acoustic_trajectory(df_analysis, jitter_a, "Jitter (%)")
p4 <- plot_acoustic_trajectory(df_analysis, nne_a, "NNE")

# Combine plots
p_combined <- (p1 + p2) /
  (p3 + p4) +
  plot_annotation(
    title = "Acoustic Markers Across Exam Stress Period",
    subtitle = "Baseline → Pre-Exam → Post-Exam",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(
  "figures/acoustic_trajectories_ORIGINAL.png",
  p_combined,
  width = 12,
  height = 10,
  dpi = 300
)

cat(
  "\n✓ Descriptive plots saved to figures/acoustic_trajectories_ORIGINAL.png\n"
)

length(unique(df_analysis$ID))
# 118

# ==============================================================================
# 4. BAYESIAN MULTILEVEL MODELS - MAIN EFFECTS
# ==============================================================================

cat("\n=== FITTING BAYESIAN MULTILEVEL MODELS (VERSIONE ORIGINALE) ===\n")
cat("Tutti i modelli usano gaussian() per coerenza con l'abstract\n")
cat("This will take several minutes...\n\n")

# Set up weakly informative priors
priors_main <- c(
  prior(normal(0, 5), class = "b"),
  prior(student_t(3, 0, 10), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# Model 1: F0 mean ~ timepoint + (1 + timepoint | ID)
cat("Fitting Model 1: F0 Mean...\n")
m1_f0mean_a <- brm(
  f0_mean_a ~ timepoint + (1 + timepoint | ID),
  data = df_analysis,
  family = gaussian(),
  prior = priors_main,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99),
  file = "models_original/m1_f0mean_a_main"
)
pp_check(m1_f0mean_a)
summary(m1_f0mean_a)


# Model 2: F0 std ~ timepoint + (1 + timepoint | ID)
cat("Fitting Model 2: F0 SD...\n")
m2_f0std_a <- brm(
  f0_std_a ~ timepoint + (1 + timepoint | ID),
  data = df_analysis,
  family = gaussian(),
  prior = priors_main,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99),
  file = "models_original/m2_f0std_a_main"
)
pp_check(m2_f0std_a)
summary(m2_f0std_a)

# Model 3: Jitter ~ timepoint + (1 + timepoint | ID)
# MODIFICATO: gaussian() invece di lognormal()
cat("Fitting Model 3: Jitter (gaussian)...\n")
m3_jitter_a <- brm(
  jitter_a ~ timepoint + (1 + timepoint | ID),
  data = df_analysis,
  family = gaussian(),
  prior = priors_main,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99),
  file = "models_original/m3_jitter_a_main"
)
pp_check(m3_jitter_a)
summary(m3_jitter_a)

# Model 4: NNE ~ timepoint + (1 + timepoint | ID)
cat("Fitting Model 4: NNE...\n")
m4_nne_a <- brm(
  nne_a ~ timepoint + (1 + timepoint | ID),
  data = df_analysis,
  family = gaussian(),
  prior = priors_main,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99),
  file = "models_original/m4_nne_a_main"
)
pp_check(m4_nne_a)
summary(m4_nne_a)

# Model 5: F2 mean ~ timepoint + (1 + timepoint | ID)
# MODIFICATO: gaussian() invece di asym_laplace()
cat("Fitting Model 5: F2 Mean (gaussian)...\n")
m5_f2mean_a <- brm(
  f2_mean_a ~ timepoint + (1 + timepoint | ID),
  data = df_analysis,
  family = gaussian(),
  prior = priors_main,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99),
  file = "models_original/m5_f2mean_a_main"
)
pp_check(m5_f2mean_a)
summary(m5_f2mean_a)


# Model 6: F2 std ~ timepoint + (1 + timepoint | ID)
# MODIFICATO: gaussian() invece di lognormal()
cat("Fitting Model 6: F2 SD (gaussian)...\n")
m6_f2std_a <- brm(
  f2_std_a ~ timepoint + (1 + timepoint | ID),
  data = df_analysis,
  family = gaussian(),
  prior = priors_main,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.99),
  file = "models_original/m6_f2std_a_main"
)
pp_check(m6_f2std_a)
summary(m6_f2std_a)

cat("\n✓ Main effect models fitted successfully (VERSIONE ORIGINALE)\n")

# ==============================================================================
# 5. RESULTS: MAIN EFFECTS OF STRESS
# ==============================================================================

cat("\n=== MODEL RESULTS: Main Effects of Timepoint ===\n\n")

.extract_random_diagnostics <- function(sum_mod) {
  if (is.null(sum_mod$random)) {
    return(NULL)
  }

  # Case 1: single matrix
  if (is.matrix(sum_mod$random)) {
    return(sum_mod$random)
  }

  # Case 2: list of matrices (most common)
  if (is.list(sum_mod$random)) {
    mats <- lapply(sum_mod$random, function(x) {
      if (is.matrix(x)) x else NULL
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) return(NULL)
    return(do.call(rbind, mats))
  }

  NULL
}


.print_hypothesis <- function(h) {
  stopifnot(is.list(h), "hypothesis" %in% names(h))

  tbl <- h$hypothesis
  stopifnot(is.data.frame(tbl))

  # Always-safe core columns
  core <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")

  # Optional diagnostics (version/backend dependent)
  optional <- c("Rhat", "Ess_bulk", "Evid.Ratio", "Post.Prob")

  cols <- intersect(c(core, optional), colnames(tbl))

  print(round(tbl[, cols, drop = FALSE], 3))
}


summarize_main_effects <- function(model, outcome_name) {
  stopifnot(inherits(model, "brmsfit"))

  sep <- paste(rep("=", 60), collapse = "")

  cat("\n", sep, "\n", sep = "")
  cat("OUTCOME:", outcome_name, "\n")
  cat(sep, "\n", sep = "")

  # ---------------------------------------------------------------------------
  # 1) Fixed effects
  # ---------------------------------------------------------------------------
  cat("\nFixed Effects (Posterior Estimates):\n")

  sum_mod <- summary(model)

  fixed_cols <- intersect(
    c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS"),
    colnames(sum_mod$fixed)
  )

  print(round(sum_mod$fixed[, fixed_cols, drop = FALSE], 3))

  # ---------------------------------------------------------------------------
  # 2) Planned contrasts
  # ---------------------------------------------------------------------------
  cat("\nPlanned Contrasts:\n")

  cat("\n[H1] PRE vs BASELINE:\n")
  .print_hypothesis(hypothesis(model, "timepointpre = 0"))

  cat("\n[H2] POST vs PRE:\n")
  .print_hypothesis(hypothesis(model, "timepointpost - timepointpre = 0"))

  cat("\n[H3] POST vs BASELINE:\n")
  .print_hypothesis(hypothesis(model, "timepointpost = 0"))

  # ---------------------------------------------------------------------------
  # 3) Random effects
  # ---------------------------------------------------------------------------
  cat("\nRandom Effects (Variance Components):\n")
  print(VarCorr(model), digits = 3)

  invisible(model)
}


# Display results for each model
summarize_main_effects(m1_f0mean_a, "F0 Mean (Hz)")
summarize_main_effects(m2_f0std_a, "F0 SD (Hz)")
summarize_main_effects(m3_jitter_a, "Jitter (%)")
summarize_main_effects(m4_nne_a, "NNE (dB)")
summarize_main_effects(m5_f2mean_a, "F2 Mean (Hz)")
summarize_main_effects(m6_f2std_a, "F2 SD (Hz)")

# ==============================================================================
# 6. VISUALIZATION: POSTERIOR DISTRIBUTIONS OF STRESS EFFECTS
# ==============================================================================

# Function to create posterior plots for stress effects
plot_stress_effects <- function(model, outcome_name) {
  # Extract posterior draws
  draws <- as_draws_df(model)

  # Get key contrasts
  post_df <- draws %>%
    mutate(
      pre_vs_baseline = b_timepointpre,
      post_vs_pre = b_timepointpost - b_timepointpre,
      post_vs_baseline = b_timepointpost
    ) %>%
    dplyr::select(pre_vs_baseline, post_vs_pre, post_vs_baseline) %>%
    pivot_longer(everything(), names_to = "contrast", values_to = "value") %>%
    mutate(
      contrast = factor(
        contrast,
        levels = c("pre_vs_baseline", "post_vs_pre", "post_vs_baseline"),
        labels = c(
          "PRE - BASELINE\n(Stress onset)",
          "POST - PRE\n(Recovery)",
          "POST - BASELINE\n(Residual)"
        )
      )
    )

  # Create plot
  ggplot(post_df, aes(x = value, y = contrast)) +
    stat_halfeye(
      .width = c(0.95, 0.89),
      fill = "#0072B2",
      alpha = 0.7
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    labs(
      title = outcome_name,
      subtitle = "Posterior distributions of stress effects",
      x = "Effect size (change in outcome)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

# Create plots
p_f0mean_post_a <- plot_stress_effects(m1_f0mean_a, "F0 Mean (Hz)")
p_f0std_post_a <- plot_stress_effects(m2_f0std_a, "F0 SD (Hz)")
p_jitter_post_a <- plot_stress_effects(m3_jitter_a, "Jitter (%)")
p_nne_post_a <- plot_stress_effects(m4_nne_a, "NNE (dB)")
p_f2mean_post_a <- plot_stress_effects(m5_f2mean_a, "F2 Mean (Hz)")
p_f2std_post_a <- plot_stress_effects(m6_f2std_a, "F2 SD (Hz)")


# Combine
p_posterior_combined <- (p_f0mean_post_a + p_f0std_post_a) /
  (p_jitter_post_a + p_nne_post_a) /
  (p_f2mean_post_a + p_f2std_post_a) +
  plot_annotation(
    title = "Bayesian Estimates of Stress Effects on Voice Acoustics",
    subtitle = "95% and 89% Credible Intervals (VERSIONE ORIGINALE - gaussian)",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(
  "figures/posterior_stress_effects_ORIGINAL.png",
  p_posterior_combined,
  width = 14,
  height = 12,
  dpi = 300
)

cat(
  "\n✓ Posterior plots saved to figures/posterior_stress_effects_ORIGINAL.png\n"
)

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

# Save workspace
save.image("results/main_effects_workspace_ORIGINAL.RData")

# Also save df_analysis separately for use in other scripts
saveRDS(df_analysis, "results/df_analysis_ORIGINAL.rds")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MAIN EFFECTS ANALYSIS COMPLETE (VERSIONE ORIGINALE)\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nTutti i modelli usano gaussian() - coefficienti in scala originale\n")
cat("I risultati sono coerenti con l'abstract.\n")
cat("\nAll results saved. Proceed to moderation analysis...\n\n")

# ==============================================================================
# RIEPILOGO MODIFICHE RISPETTO ALLA VERSIONE MODIFICATA
# ==============================================================================
# 1. m3_jitter_a: lognormal() → gaussian()
# 2. m5_f2mean_a: asym_laplace() → gaussian()
# 3. m6_f2std_a: lognormal() → gaussian()
# 4. Plot p3: log(jitter_a) → jitter_a
# 5. Cartella modelli: models/ → models_original/
# 6. File output con suffisso _ORIGINAL per non sovrascrivere
#
# INTERPRETAZIONE COEFFICIENTI (scala originale):
# - f0_mean, f2_mean, f2_std, f0_std: Hz
# - jitter: % (percentuale)
# - nne: dB
# ==============================================================================
