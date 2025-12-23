# ==============================================================================
# 05_calculate_within_subject_variability.R
# Calculate within-subject variability for PID-5 from ALL EMA data
# ==============================================================================
# PURPOSE:
#   Compute robust measures of emotional lability/instability using all
#   available EMA assessments (not just the 3 timepoints with voice)
#
# OUTPUT:
#   Dataset with person-level variability measures to use as moderators
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(rio)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("WITHIN-SUBJECT VARIABILITY: PID-5 EMA\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD ALL EMA DATA
# ==============================================================================

cat("Loading complete EMA data...\n")

data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
if (!file.exists(data_path)) {
  stop("Run 01_clean_after_merge_FINAL.R first")
}

d_ema <- rio::import(data_path)

# Rename for consistency
d_ema <- d_ema %>%
  rename(ID = any_of(c("user_id", "UserID", "id_anon")))

# PID-5 variables
pid5_vars <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

cat("✓ EMA data loaded: N =", nrow(d_ema), "observations\n")
cat("  N subjects:", n_distinct(d_ema$ID), "\n\n")

# ==============================================================================
# 2. COMPUTE WITHIN-SUBJECT VARIABILITY MEASURES
# ==============================================================================

cat("Computing within-subject variability measures...\n\n")

# For each subject, calculate:
# 1. SD - Standard deviation (overall variability)
# 2. IQR - Interquartile range (robust variability)
# 3. RMSSD - Root mean square of successive differences (instability/lability)
# 4. n_obs - Number of observations (for reliability check)

variability_measures <- d_ema %>%
  # Keep only valid PID-5 observations
  filter(!is.na(pid5_negative_affectivity)) %>%
  # Sort by subject and date
  arrange(ID, day) %>%
  # Group by subject
  group_by(ID) %>%
  # Calculate variability for each PID-5 domain
  summarise(
    # Sample size
    n_ema_total = n(),

    # === NEGATIVE AFFECTIVITY ===
    pid5_negative_affectivity_mean = mean(
      pid5_negative_affectivity,
      na.rm = TRUE
    ),
    pid5_negative_affectivity_sd = sd(pid5_negative_affectivity, na.rm = TRUE),
    pid5_negative_affectivity_iqr = IQR(
      pid5_negative_affectivity,
      na.rm = TRUE
    ),
    pid5_negative_affectivity_rmssd = sqrt(mean(
      diff(pid5_negative_affectivity)^2,
      na.rm = TRUE
    )),

    # === DETACHMENT ===
    pid5_detachment_mean = mean(pid5_detachment, na.rm = TRUE),
    pid5_detachment_sd = sd(pid5_detachment, na.rm = TRUE),
    pid5_detachment_iqr = IQR(pid5_detachment, na.rm = TRUE),
    pid5_detachment_rmssd = sqrt(mean(
      diff(pid5_detachment)^2,
      na.rm = TRUE
    )),

    # === ANTAGONISM ===
    pid5_antagonism_mean = mean(pid5_antagonism, na.rm = TRUE),
    pid5_antagonism_sd = sd(pid5_antagonism, na.rm = TRUE),
    pid5_antagonism_iqr = IQR(pid5_antagonism, na.rm = TRUE),
    pid5_antagonism_rmssd = sqrt(mean(
      diff(pid5_antagonism)^2,
      na.rm = TRUE
    )),

    # === DISINHIBITION ===
    pid5_disinhibition_mean = mean(pid5_disinhibition, na.rm = TRUE),
    pid5_disinhibition_sd = sd(pid5_disinhibition, na.rm = TRUE),
    pid5_disinhibition_iqr = IQR(pid5_disinhibition, na.rm = TRUE),
    pid5_disinhibition_rmssd = sqrt(mean(
      diff(pid5_disinhibition)^2,
      na.rm = TRUE
    )),

    # === PSYCHOTICISM ===
    pid5_psychoticism_mean = mean(pid5_psychoticism, na.rm = TRUE),
    pid5_psychoticism_sd = sd(pid5_psychoticism, na.rm = TRUE),
    pid5_psychoticism_iqr = IQR(pid5_psychoticism, na.rm = TRUE),
    pid5_psychoticism_rmssd = sqrt(mean(
      diff(pid5_psychoticism)^2,
      na.rm = TRUE
    )),

    .groups = "drop"
  )

cat("✓ Variability computed for", nrow(variability_measures), "subjects\n\n")

# ==============================================================================
# 3. DATA QUALITY CHECKS
# ==============================================================================

cat(rep("=", 70), "\n")
cat("DATA QUALITY CHECKS\n")
cat(rep("=", 70), "\n\n")

# Check minimum observations needed for reliable variability
cat("EMA observations per subject:\n")
cat("  Min:", min(variability_measures$n_ema_total), "\n")
cat("  Median:", median(variability_measures$n_ema_total), "\n")
cat("  Mean:", round(mean(variability_measures$n_ema_total), 1), "\n")
cat("  Max:", max(variability_measures$n_ema_total), "\n\n")

# Flag subjects with too few observations
min_obs_threshold <- 10
n_low_reliability <- sum(variability_measures$n_ema_total < min_obs_threshold)

if (n_low_reliability > 0) {
  cat(
    "⚠ WARNING:",
    n_low_reliability,
    "subjects have <",
    min_obs_threshold,
    "EMA observations\n"
  )
  cat("  Variability estimates may be unreliable for these subjects\n")
  cat("  Consider excluding them or flagging in analysis\n\n")

  variability_measures <- variability_measures %>%
    mutate(reliable_variability = n_ema_total >= min_obs_threshold)
} else {
  cat("✓ All subjects have ≥", min_obs_threshold, "observations\n\n")
  variability_measures$reliable_variability <- TRUE
}

# Summary statistics
cat("Variability measures summary:\n\n")

for (domain in c(
  "negative_affectivity",
  "detachment",
  "antagonism",
  "disinhibition",
  "psychoticism"
)) {
  var_sd <- paste0("pid5_", domain, "_sd")
  var_rmssd <- paste0("pid5_", domain, "_rmssd")

  cat(str_to_title(str_replace_all(domain, "_", " ")), ":\n")
  cat(sprintf(
    "  SD:    M = %.3f, SD = %.3f, Range = [%.3f, %.3f]\n",
    mean(variability_measures[[var_sd]], na.rm = TRUE),
    sd(variability_measures[[var_sd]], na.rm = TRUE),
    min(variability_measures[[var_sd]], na.rm = TRUE),
    max(variability_measures[[var_sd]], na.rm = TRUE)
  ))
  cat(sprintf(
    "  RMSSD: M = %.3f, SD = %.3f, Range = [%.3f, %.3f]\n\n",
    mean(variability_measures[[var_rmssd]], na.rm = TRUE),
    sd(variability_measures[[var_rmssd]], na.rm = TRUE),
    min(variability_measures[[var_rmssd]], na.rm = TRUE),
    max(variability_measures[[var_rmssd]], na.rm = TRUE)
  ))
}

# ==============================================================================
# 4. CORRELATIONS: MEAN VS VARIABILITY
# ==============================================================================

cat(rep("=", 70), "\n")
cat("CORRELATIONS: Mean vs Variability\n")
cat(rep("=", 70), "\n\n")

cat("Are mean levels correlated with variability?\n")
cat("(High correlations suggest multicollinearity in models)\n\n")

for (domain in c(
  "negative_affectivity",
  "detachment",
  "antagonism",
  "disinhibition",
  "psychoticism"
)) {
  var_mean <- paste0("pid5_", domain, "_mean")
  var_sd <- paste0("pid5_", domain, "_sd")

  cor_val <- cor(
    variability_measures[[var_mean]],
    variability_measures[[var_sd]],
    use = "complete.obs"
  )

  domain_name <- str_to_title(str_replace_all(domain, "_", " "))
  cat(sprintf("  %s: r = %.3f", domain_name, cor_val))

  if (abs(cor_val) > 0.70) {
    cat(" ← HIGH (multicollinearity concern)")
  } else if (abs(cor_val) > 0.50) {
    cat(" ← MODERATE")
  }
  cat("\n")
}

cat("\n")

# ==============================================================================
# 5. MERGE WITH VOICE DATA
# ==============================================================================

cat(rep("=", 70), "\n")
cat("MERGING WITH VOICE DATA\n")
cat(rep("=", 70), "\n\n")

# Load voice data
if (!file.exists("results/df_analysis.rds")) {
  stop("Run 02_voice_personality_analysis_FINAL.R first")
}

df_voice <- readRDS("results/df_analysis.rds")

cat("✓ Voice data loaded: N =", nrow(df_voice), "observations\n")

# Merge - variability is constant per subject
df_variability <- df_voice %>%
  left_join(variability_measures, by = "ID")

# Check merge success
n_missing <- sum(is.na(df_variability$pid5_negative_affectivity_mean))
cat("✓ Merge complete\n")
cat(
  "  Observations with variability data:",
  nrow(df_variability) - n_missing,
  "\n"
)

if (n_missing > 0) {
  cat("  ⚠ Missing variability data:", n_missing, "observations\n")
  cat("    (Subjects without sufficient EMA data)\n")
}

cat("\n")

# ==============================================================================
# 6. CENTER PREDICTORS
# ==============================================================================

cat("Centering variability predictors...\n")

# Center all mean and variability measures
for (domain in c(
  "negative_affectivity",
  "detachment",
  "antagonism",
  "disinhibition",
  "psychoticism"
)) {
  # Mean (already have this from trait analysis, but recompute for consistency)
  var_mean <- paste0("pid5_", domain, "_mean")
  var_mean_c <- paste0("pid5_", domain, "_mean_c")
  df_variability[[var_mean_c]] <- scale(
    df_variability[[var_mean]],
    center = TRUE,
    scale = FALSE
  )[, 1]

  # SD variability
  var_sd <- paste0("pid5_", domain, "_sd")
  var_sd_c <- paste0("pid5_", domain, "_sd_c")
  df_variability[[var_sd_c]] <- scale(
    df_variability[[var_sd]],
    center = TRUE,
    scale = FALSE
  )[, 1]

  # RMSSD variability
  var_rmssd <- paste0("pid5_", domain, "_rmssd")
  var_rmssd_c <- paste0("pid5_", domain, "_rmssd_c")
  df_variability[[var_rmssd_c]] <- scale(
    df_variability[[var_rmssd]],
    center = TRUE,
    scale = FALSE
  )[, 1]
}

cat("✓ Predictors centered\n\n")

# ==============================================================================
# 7. SAVE OUTPUT
# ==============================================================================

cat("Saving output...\n")

dir.create("results/variability", showWarnings = FALSE, recursive = TRUE)

# Save full dataset
saveRDS(df_variability, "results/variability/df_with_variability.rds")
rio::export(df_variability, "results/variability/df_with_variability.csv")

# Save variability measures only (for reference)
saveRDS(variability_measures, "results/variability/variability_measures.rds")
rio::export(
  variability_measures,
  "results/variability/variability_measures.csv"
)

cat("✓ Saved: df_with_variability.rds\n")
cat("✓ Saved: variability_measures.csv\n\n")

cat(rep("=", 70), "\n")
cat("VARIABILITY COMPUTATION COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Next steps:\n")
cat(
  "1. Run 06_moderation_with_variability.R to test variability as moderator\n"
)
cat("2. Compare models: Mean-only vs Mean+SD vs SD-only\n")
cat("3. Visualize how lability/instability relates to stress reactivity\n\n")

cat("Key interpretations:\n")
cat("- SD = Overall variability (are they variable across all timepoints?)\n")
cat("- RMSSD = Instability (how much do they fluctuate day-to-day?)\n")
cat("- High SD/RMSSD = Emotional lability, dysregulation\n\n")
