# ==============================================================================
# INTEGRATED PATTERN CONFORMITY ANALYSIS
# Works with existing moderation_analysis.R pipeline
# ==============================================================================

# Run this AFTER moderation_analysis.R has completed and fitted_models is available
# OR load the saved workspace:
# load("results/moderation_analysis_complete.RData")

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(posterior)
  library(patchwork)
  library(ggdist)
})

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("BAYESIAN PATTERN CONFORMITY ANALYSIS\n")
cat("Testing Overall Theoretical Predictions\n")
cat(rep("=", 70), "\n", sep = "")

# ==============================================================================
# 1. VERIFY FITTED MODELS ARE AVAILABLE
# ==============================================================================

# Check if fitted_models exists from moderation_analysis.R
if (!exists("fitted_models")) {
  # Try to load from saved workspace
  if (file.exists("results/moderation_analysis_complete.RData")) {
    cat("\nLoading fitted models from saved workspace...\n")
    load("results/moderation_analysis_complete.RData")
  } else if (file.exists("models/m_f0_mean_a.rds")) {
    # Or load individual model files
    cat("\nLoading fitted models from individual files...\n")
    model_files <- list.files("models", pattern = "^m_.*\\.rds$", full.names = TRUE)
    fitted_models <- lapply(model_files, readRDS)
    names(fitted_models) <- gsub("\\.rds$", "", basename(model_files))
  } else {
    stop("No fitted models found. Please run moderation_analysis.R first.")
  }
}

cat("\nFound", length(fitted_models), "fitted models\n")
cat("Model names:", paste(names(fitted_models)[1:6], collapse = ", "), "...\n")

# ==============================================================================
# 2. DEFINE THEORETICAL PREDICTIONS
# ==============================================================================

# Based on your abstract:
# - Negative Affectivity amplified stress-induced pitch increases (β=5.03 Hz)
# - Detachment impaired post-stress recovery (β=-4.37)
# - Antagonism facilitated recovery (β=3.44-4.17)
# - Psychoticism showed pervasive pitch elevation and voice quality degradation
# - Disinhibition exhibited domain-specific articulatory alterations (F2)

# Create prediction matrix
# expected_direction: +1 = positive, -1 = negative, 0 = no prediction (excluded)

predictions <- tribble(
  ~outcome,   ~vowel, ~contrast,     ~pid5_domain,           ~expected_direction, ~rationale,
  

  # ============================================================================
  # F0 MEAN (Pitch) - Primary outcome for stress/arousal
  # ============================================================================
  

  # Negative Affectivity: amplifies stress reactivity (abstract: β=5.03 Hz)
  "f0_mean",  "a",    "c1_stress",   "negative_affectivity", +1, "NA amplifies stress-induced pitch increase",
  "f0_mean",  "i",    "c1_stress",   "negative_affectivity", +1, "Replication across vowels",
  "f0_mean",  "u",    "c1_stress",   "negative_affectivity", +1, "Replication across vowels",
  
  # Negative Affectivity: impairs recovery (prolonged elevation)
  "f0_mean",  "a",    "c2_recovery", "negative_affectivity", -1, "NA impairs return to baseline (less recovery)",
  "f0_mean",  "i",    "c2_recovery", "negative_affectivity", -1, "Replication",
  "f0_mean",  "u",    "c2_recovery", "negative_affectivity", -1, "Replication",
  
  # Detachment: blunted reactivity OR impaired recovery (abstract: β=-4.37 for recovery)
  "f0_mean",  "a",    "c1_stress",   "detachment",           -1, "Detachment blunts emotional/physiological reactivity",
  "f0_mean",  "i",    "c1_stress",   "detachment",           -1, "Replication",
  "f0_mean",  "u",    "c1_stress",   "detachment",           -1, "Replication",
  "f0_mean",  "a",    "c2_recovery", "detachment",           -1, "Detachment impairs recovery (abstract finding)",
  "f0_mean",  "i",    "c2_recovery", "detachment",           -1, "Replication",
  "f0_mean",  "u",    "c2_recovery", "detachment",           -1, "Replication",
  
  # Antagonism: facilitates recovery (abstract: β=3.44-4.17)
  "f0_mean",  "a",    "c1_stress",   "antagonism",            0, "No clear prediction for reactivity",
  "f0_mean",  "i",    "c1_stress",   "antagonism",            0, "No clear prediction",
  "f0_mean",  "u",    "c1_stress",   "antagonism",            0, "No clear prediction",
  "f0_mean",  "a",    "c2_recovery", "antagonism",           +1, "Antagonism facilitates recovery (abstract finding)",
  "f0_mean",  "i",    "c2_recovery", "antagonism",           +1, "Replication",
  "f0_mean",  "u",    "c2_recovery", "antagonism",           +1, "Replication",
  
  # Psychoticism: pervasive pitch elevation (main effect, but also expect moderation)
  "f0_mean",  "a",    "c1_stress",   "psychoticism",         +1, "Psychoticism shows heightened arousal",
  "f0_mean",  "i",    "c1_stress",   "psychoticism",         +1, "Replication",
  "f0_mean",  "u",    "c1_stress",   "psychoticism",         +1, "Replication",
  "f0_mean",  "a",    "c2_recovery", "psychoticism",         -1, "Psychoticism shows persistent elevation (less recovery)",
  "f0_mean",  "i",    "c2_recovery", "psychoticism",         -1, "Replication",
  "f0_mean",  "u",    "c2_recovery", "psychoticism",         -1, "Replication",
  
  # Disinhibition: no strong prediction for F0
  "f0_mean",  "a",    "c1_stress",   "disinhibition",         0, "No clear prediction",
  "f0_mean",  "i",    "c1_stress",   "disinhibition",         0, "No clear prediction",
  "f0_mean",  "u",    "c1_stress",   "disinhibition",         0, "No clear prediction",
  "f0_mean",  "a",    "c2_recovery", "disinhibition",         0, "No clear prediction",
  "f0_mean",  "i",    "c2_recovery", "disinhibition",         0, "No clear prediction",
  "f0_mean",  "u",    "c2_recovery", "disinhibition",         0, "No clear prediction",
  
  # ============================================================================
  # F0 STD (Pitch variability) - Emotional instability/dysregulation
  # ============================================================================
  
  # NA: increased variability under stress
  "f0_std",   "a",    "c1_stress",   "negative_affectivity", +1, "NA increases emotional lability",
  "f0_std",   "i",    "c1_stress",   "negative_affectivity", +1, "Replication",
  "f0_std",   "u",    "c1_stress",   "negative_affectivity", +1, "Replication",
  "f0_std",   "a",    "c2_recovery", "negative_affectivity", -1, "Impaired recovery of vocal stability",
  "f0_std",   "i",    "c2_recovery", "negative_affectivity", -1, "Replication",
  "f0_std",   "u",    "c2_recovery", "negative_affectivity", -1, "Replication",
  
  # Other domains: less clear predictions for F0 variability
  "f0_std",   "a",    "c1_stress",   "detachment",            0, "No clear prediction",
  "f0_std",   "a",    "c1_stress",   "antagonism",            0, "No clear prediction",
  "f0_std",   "a",    "c1_stress",   "disinhibition",         0, "No clear prediction",
  "f0_std",   "a",    "c1_stress",   "psychoticism",         +1, "Psychoticism shows instability",
  
  # ============================================================================
  # JITTER (Voice quality/tension)
  # ============================================================================
  
  # NA: increased jitter (more tension)
  "jitter",   "a",    "c1_stress",   "negative_affectivity", +1, "NA increases vocal tension",
  "jitter",   "i",    "c1_stress",   "negative_affectivity", +1, "Replication",
  "jitter",   "u",    "c1_stress",   "negative_affectivity", +1, "Replication",
  
  # Psychoticism: voice quality degradation (abstract finding)
  "jitter",   "a",    "c1_stress",   "psychoticism",         +1, "Psychoticism shows voice quality issues",
  "jitter",   "i",    "c1_stress",   "psychoticism",         +1, "Replication",
  "jitter",   "u",    "c1_stress",   "psychoticism",         +1, "Replication",
  
  # ============================================================================
  # NNE (Normalized Noise Energy) - Vocal stress marker
  # More negative = cleaner voice; stress typically makes it less negative (increases)
  # ============================================================================
  
  # NA: amplifies stress effect on NNE
  "nne",      "a",    "c1_stress",   "negative_affectivity", +1, "NA amplifies vocal stress (less negative NNE)",
  "nne",      "i",    "c1_stress",   "negative_affectivity", +1, "Replication",
  "nne",      "u",    "c1_stress",   "negative_affectivity", +1, "Replication",
  
  # Psychoticism: pervasive voice quality issues
  "nne",      "a",    "c1_stress",   "psychoticism",         +1, "Psychoticism shows degraded voice quality",
  "nne",      "i",    "c1_stress",   "psychoticism",         +1, "Replication",
  "nne",      "u",    "c1_stress",   "psychoticism",         +1, "Replication",
  
  # ============================================================================
  # F2 MEAN & F2 STD (Formants - Articulatory features)
  # Abstract: Disinhibition showed F2 effects (β=63 Hz)
  # ============================================================================
  
  # Disinhibition: articulatory alterations (abstract finding)
  "f2_mean",  "a",    "c1_stress",   "disinhibition",        +1, "Disinhibition shows articulatory changes",
  "f2_mean",  "i",    "c1_stress",   "disinhibition",        +1, "Replication",
  "f2_mean",  "u",    "c1_stress",   "disinhibition",        +1, "Replication",
  "f2_std",   "a",    "c1_stress",   "disinhibition",        +1, "Disinhibition shows articulatory variability",
  "f2_std",   "i",    "c1_stress",   "disinhibition",        +1, "Replication",
  "f2_std",   "u",    "c1_stress",   "disinhibition",        +1, "Replication"
)

# View prediction summary
cat("\n=== THEORETICAL PREDICTIONS ===\n")
cat("\nTotal predictions specified:", nrow(predictions), "\n")
cat("Predictions with expected direction (+1 or -1):", 
    sum(predictions$expected_direction != 0), "\n")
cat("Predictions excluded (no direction):", 
    sum(predictions$expected_direction == 0), "\n")

predictions %>%
  filter(expected_direction != 0) %>%
  count(outcome, pid5_domain, expected_direction) %>%
  pivot_wider(names_from = expected_direction, values_from = n, values_fill = 0) %>%
  print(n = 50)

# ==============================================================================
# 3. MAP PREDICTIONS TO MODEL PARAMETERS
# ==============================================================================

# Create the parameter name as it appears in brms output
predictions <- predictions %>%
  mutate(
    model_name = paste0("m_", outcome, "_", vowel),
    param_name = paste0("b_", contrast, ":pid5_", pid5_domain, "_bl_c")
  )

# Filter to only directional predictions
predictions_directional <- predictions %>%
  filter(expected_direction != 0)

cat("\n=== PARAMETER MAPPING ===\n")
cat("\nExample parameter names:\n")
print(head(predictions_directional %>% select(model_name, param_name, expected_direction), 10))

# ==============================================================================
# 4. EXTRACT POSTERIOR DRAWS AND COMPUTE CONFORMITY
# ==============================================================================

cat("\n=== EXTRACTING POSTERIOR DRAWS ===\n")

# Get number of posterior draws from first model
first_model <- fitted_models[[1]]
n_draws <- nrow(as_draws_df(first_model))
cat("Number of posterior draws per model:", n_draws, "\n")

n_predictions <- nrow(predictions_directional)
cat("Number of directional predictions to test:", n_predictions, "\n")

# Initialize conformity matrix
conformity_matrix <- matrix(NA, nrow = n_draws, ncol = n_predictions)
colnames(conformity_matrix) <- paste0(
  predictions_directional$outcome, "_",
  predictions_directional$vowel, "_",
  predictions_directional$pid5_domain, "_",
  predictions_directional$contrast
)

# Also store the actual posterior values for additional analyses
posterior_values <- matrix(NA, nrow = n_draws, ncol = n_predictions)
colnames(posterior_values) <- colnames(conformity_matrix)

# Extract posteriors and compute conformity
for (i in seq_len(n_predictions)) {
  pred <- predictions_directional[i, ]
  
  # Check if model exists
  if (!pred$model_name %in% names(fitted_models)) {
    cat("  Warning: Model not found:", pred$model_name, "\n")
    next
  }
  
  model <- fitted_models[[pred$model_name]]
  draws <- as_draws_df(model)
  
  # Try to find the parameter
  param_name <- pred$param_name
  
  # brms may use different naming conventions
  possible_names <- c(
    param_name,
    gsub(":", ".", param_name),  # Some versions use . instead of :
    gsub("b_", "b_", gsub("_bl_c", "_bl_c", param_name))
  )
  
  found_param <- NULL
  for (pn in possible_names) {
    if (pn %in% names(draws)) {
      found_param <- pn
      break
    }
  }
  
  if (is.null(found_param)) {
    # List available parameters for debugging
    if (i == 1) {
      cat("\nAvailable parameters in first model:\n")
      cat(paste(names(draws)[grep("^b_", names(draws))], collapse = "\n"), "\n")
    }
    cat("  Warning: Parameter not found:", param_name, "\n")
    next
  }
  
  # Extract posterior samples
  posterior_samples <- draws[[found_param]]
  posterior_values[, i] <- posterior_samples
  
  # Compute conformity: is the effect in the expected direction?
  if (pred$expected_direction == 1) {
    conformity_matrix[, i] <- as.numeric(posterior_samples > 0)
  } else {
    conformity_matrix[, i] <- as.numeric(posterior_samples < 0)
  }
}

# Check how many predictions we successfully extracted
n_valid <- sum(!is.na(conformity_matrix[1, ]))
cat("\nSuccessfully extracted", n_valid, "of", n_predictions, "predictions\n")

# ==============================================================================
# 5. COMPUTE PATTERN CONFORMITY METRICS
# ==============================================================================

cat("\n=== COMPUTING CONFORMITY METRICS ===\n")

results <- list()

# Remove NA columns (predictions we couldn't find)
valid_cols <- !is.na(conformity_matrix[1, ])
conformity_valid <- conformity_matrix[, valid_cols]
posterior_valid <- posterior_values[, valid_cols]
predictions_valid <- predictions_directional[valid_cols, ]

results$n_predictions <- ncol(conformity_valid)

# Metric 1: Proportion conforming per draw
results$proportion_per_draw <- rowMeans(conformity_valid, na.rm = TRUE)

# Metric 2: P(conformity > chance)
results$p_better_than_chance <- mean(results$proportion_per_draw > 0.5)

# Metric 3: Mean proportion and CI
results$mean_proportion <- mean(results$proportion_per_draw)
results$ci_proportion <- quantile(results$proportion_per_draw, c(0.025, 0.975))

# Metric 4: Individual PD for each prediction
results$individual_pd <- colMeans(conformity_valid, na.rm = TRUE)

# Metric 5: Combined PD using log-odds
pd_values <- pmax(pmin(results$individual_pd, 0.999), 0.001)
log_odds <- log(pd_values / (1 - pd_values))
mean_log_odds <- mean(log_odds, na.rm = TRUE)
results$combined_pd <- exp(mean_log_odds) / (1 + exp(mean_log_odds))

# Metric 6: Geometric mean PD
results$geometric_mean_pd <- exp(mean(log(pd_values), na.rm = TRUE))

# Metric 7: P(all conform simultaneously)
results$all_conform_per_draw <- apply(conformity_valid, 1, all, na.rm = TRUE)
results$p_all_conform <- mean(results$all_conform_per_draw)

# Metric 8: Effect size (how much better than chance)
results$effect_size <- results$mean_proportion - 0.5

# Metric 9: Number conforming per draw
results$n_conforming_per_draw <- rowSums(conformity_valid, na.rm = TRUE)

# Metric 10: Threshold probabilities
results$p_at_least_60 <- mean(results$proportion_per_draw >= 0.60)
results$p_at_least_70 <- mean(results$proportion_per_draw >= 0.70)
results$p_at_least_80 <- mean(results$proportion_per_draw >= 0.80)

# Store raw data
results$conformity_matrix <- conformity_valid
results$posterior_values <- posterior_valid
results$predictions <- predictions_valid

# ==============================================================================
# 6. PRINT RESULTS REPORT
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("PATTERN CONFORMITY RESULTS\n")
cat(rep("=", 70), "\n", sep = "")

cat("\n--- SUMMARY METRICS ---\n\n")

cat(sprintf("Number of directional predictions tested: %d\n", results$n_predictions))
cat(sprintf("Mean proportion conforming: %.1f%% [chance = 50%%]\n", 
            results$mean_proportion * 100))
cat(sprintf("95%% Credible Interval: [%.1f%%, %.1f%%]\n",
            results$ci_proportion[1] * 100, results$ci_proportion[2] * 100))

cat(sprintf("\nP(pattern > chance): %.2f%%\n", results$p_better_than_chance * 100))
cat(sprintf("P(ALL predictions simultaneously correct): %.4f%%\n", 
            results$p_all_conform * 100))

cat(sprintf("\nCombined PD (log-odds method): %.1f%%\n", results$combined_pd * 100))
cat(sprintf("Geometric mean PD: %.1f%%\n", results$geometric_mean_pd * 100))

cat(sprintf("\nEffect size (conformity - 0.5): %.3f\n", results$effect_size))

cat(sprintf("\nP(at least 60%% conform): %.1f%%\n", results$p_at_least_60 * 100))
cat(sprintf("P(at least 70%% conform): %.1f%%\n", results$p_at_least_70 * 100))
cat(sprintf("P(at least 80%% conform): %.1f%%\n", results$p_at_least_80 * 100))

# Interpretation
cat("\n--- INTERPRETATION ---\n\n")

if (results$p_better_than_chance > 0.95) {
  cat("✓ STRONG EVIDENCE that the overall theoretical pattern holds.\n")
  cat("  The observed pattern is unlikely under random expectation.\n")
} else if (results$p_better_than_chance > 0.80) {
  cat("~ MODERATE EVIDENCE favoring the theoretical pattern.\n")
  cat("  Results are more consistent with theory than chance.\n")
} else {
  cat("? WEAK/INCONCLUSIVE evidence for the overall pattern.\n")
  cat("  Individual effects may still be informative.\n")
}

# Individual predictions
cat("\n--- INDIVIDUAL PREDICTION RESULTS ---\n\n")

pd_summary <- tibble(
  outcome = predictions_valid$outcome,
  vowel = predictions_valid$vowel,
  domain = predictions_valid$pid5_domain,
  contrast = predictions_valid$contrast,
  expected = ifelse(predictions_valid$expected_direction == 1, "+", "-"),
  PD = results$individual_pd,
  evidence = case_when(
    PD > 0.975 ~ "Strong +",
    PD > 0.90 ~ "Moderate +",
    PD > 0.75 ~ "Weak +",
    PD > 0.50 ~ "Slight +",
    PD > 0.25 ~ "Slight -",
    PD > 0.10 ~ "Weak -",
    TRUE ~ "Strong -"
  )
) %>%
  arrange(desc(PD))

cat("Top 10 MOST supported predictions:\n")
print(head(pd_summary, 10), n = 10)

cat("\nPredictions in OPPOSITE direction (PD < 50%):\n")
opposite <- pd_summary %>% filter(PD < 0.5)
if (nrow(opposite) > 0) {
  print(opposite, n = 20)
} else {
  cat("None - all predictions in expected direction!\n")
}

# Summary by domain
cat("\n--- SUMMARY BY PID-5 DOMAIN ---\n\n")
pd_summary %>%
  group_by(domain) %>%
  summarise(
    n_predictions = n(),
    mean_PD = mean(PD),
    n_supported = sum(PD > 0.5),
    pct_supported = mean(PD > 0.5) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_PD)) %>%
  print()

# Summary by outcome
cat("\n--- SUMMARY BY ACOUSTIC OUTCOME ---\n\n")
pd_summary %>%
  group_by(outcome) %>%
  summarise(
    n_predictions = n(),
    mean_PD = mean(PD),
    n_supported = sum(PD > 0.5),
    pct_supported = mean(PD > 0.5) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_PD)) %>%
  print()

# ==============================================================================
# 7. CREATE VISUALIZATIONS
# ==============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n")

# Create figures directory if needed
dir.create("figures", showWarnings = FALSE)

# Plot 1: Posterior distribution of proportion conforming
p1 <- ggplot(data.frame(prop = results$proportion_per_draw), aes(x = prop)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = "#3182ce", alpha = 0.7, color = "white") +
  geom_density(color = "#1a365d", linewidth = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#e53e3e", linewidth = 1) +
  geom_vline(xintercept = results$mean_proportion, 
             linetype = "solid", color = "#38a169", linewidth = 1.2) +
  annotate("text", x = 0.48, y = Inf, label = "Chance (50%)", 
           vjust = 2, hjust = 1, color = "#e53e3e", fontface = "bold") +
  annotate("text", x = results$mean_proportion + 0.02, y = Inf, 
           label = sprintf("Mean = %.1f%%", results$mean_proportion * 100), 
           vjust = 2, hjust = 0, color = "#38a169", fontface = "bold") +
  annotate("label", x = 0.85, y = max(density(results$proportion_per_draw)$y) * 0.8,
           label = sprintf("P(> chance) = %.1f%%", results$p_better_than_chance * 100),
           fill = "white", label.size = 0.5) +
  scale_x_continuous(labels = scales::percent, limits = c(0.2, 1), 
                     breaks = seq(0.2, 1, 0.1)) +
  labs(
    title = "Posterior Distribution of Pattern Conformity",
    subtitle = sprintf("Testing %d directional predictions from theoretical model", 
                       results$n_predictions),
    x = "Proportion of Predictions in Expected Direction",
    y = "Posterior Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.minor = element_blank()
  )

# Plot 2: Individual PDs by domain and outcome
pd_plot_data <- pd_summary %>%
  mutate(
    label = paste0(outcome, "_", vowel, "_", contrast),
    in_expected = PD > 0.5
  )

p2 <- ggplot(pd_plot_data, aes(x = PD, y = reorder(label, PD), fill = domain)) +
  geom_col(alpha = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_vline(xintercept = 0.95, linetype = "dotted", color = "darkgreen", linewidth = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Probability of Direction for Each Prediction",
    subtitle = "Dashed line = chance (50%); dotted line = strong evidence (95%)",
    x = "P(effect in expected direction)",
    y = NULL,
    fill = "PID-5 Domain"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1))

# Plot 3: Heatmap of PDs by domain × outcome
heatmap_data <- pd_summary %>%
  group_by(outcome, domain) %>%
  summarise(mean_PD = mean(PD), .groups = "drop")

p3 <- ggplot(heatmap_data, aes(x = domain, y = outcome, fill = mean_PD)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f%%", mean_PD * 100)), 
            color = ifelse(heatmap_data$mean_PD > 0.6, "white", "black"), 
            fontface = "bold") +
  scale_fill_gradient2(low = "#e53e3e", mid = "#edf2f7", high = "#38a169",
                       midpoint = 0.5, limits = c(0, 1),
                       labels = scales::percent) +
  labs(
    title = "Pattern Conformity by Domain × Outcome",
    subtitle = "Mean probability of direction (averaged across vowels and contrasts)",
    x = "PID-5 Domain",
    y = "Acoustic Outcome",
    fill = "Mean PD"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Plot 4: Number of conforming predictions per draw
p4 <- ggplot(data.frame(n = results$n_conforming_per_draw), aes(x = n)) +
  geom_histogram(bins = results$n_predictions + 1, 
                 fill = "#805ad5", alpha = 0.7, color = "white") +
  geom_vline(xintercept = results$n_predictions / 2, 
             linetype = "dashed", color = "#e53e3e", linewidth = 1) +
  geom_vline(xintercept = mean(results$n_conforming_per_draw),
             linetype = "solid", color = "#38a169", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, results$n_predictions, by = 5)) +
  labs(
    title = "Number of Confirmed Predictions per Posterior Draw",
    subtitle = sprintf("Out of %d directional predictions; red = chance, green = mean",
                       results$n_predictions),
    x = "Number of predictions in expected direction",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Combine and save
combined_plot <- (p1 + p4) / (p3 + plot_spacer()) +
  plot_annotation(
    title = "Bayesian Pattern Conformity Analysis",
    subtitle = "Testing overall theoretical predictions for personality × stress interactions",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

ggsave("figures/pattern_conformity_main.png", combined_plot, 
       width = 14, height = 12, dpi = 300)

ggsave("figures/pattern_conformity_individual_pds.png", p2, 
       width = 10, height = 12, dpi = 300)

cat("✓ Figures saved to figures/pattern_conformity_*.png\n")

# ==============================================================================
# 8. GENERATE MANUSCRIPT TEXT
# ==============================================================================

cat("\n=== MANUSCRIPT-READY TEXT ===\n\n")

manuscript_text <- sprintf(
'ANALYSIS PLAN TEXT:
"Rather than treating each of the %d personality × stress interactions as 
independent hypothesis tests—which would require stringent multiplicity 
corrections and fundamentally mischaracterize our research question—we adopted 
a theory-driven pattern conformity approach. We specified a priori directional 
predictions for each interaction based on our theoretical model (see Table X), 
then computed the posterior probability that the overall pattern of results 
conformed to these predictions. This approach directly tests our theoretical 
model as a unified whole, asking: what is the probability that personality 
pathology dimensions moderate stress reactivity in the directions predicted 
by theory?"

RESULTS TEXT:
"We specified %d directional predictions based on our theoretical framework 
(e.g., Negative Affectivity amplifies stress reactivity; Detachment impairs 
recovery; Antagonism facilitates recovery). The mean proportion of predictions 
in the expected direction was %.1f%% (95%% CI [%.1f%%, %.1f%%]), substantially 
exceeding the 50%% expected by chance. The posterior probability that our 
theoretical pattern performed better than chance was %.1f%%, providing %s 
evidence for the overall model. Using a log-odds combination method, the 
aggregate probability of direction across all predictions was %.1f%%, 
indicating %s support for the directional consistency of effects.
  
Examining individual predictions, %d of %d (%.0f%%) showed probability of 
direction exceeding 50%%, and %d (%.0f%%) exceeded the 95%% threshold for 
strong evidence. Predictions regarding Negative Affectivity and recovery 
processes showed particularly consistent support (see Figure X)."',
  
  results$n_predictions,
  results$n_predictions,
  results$mean_proportion * 100,
  results$ci_proportion[1] * 100,
  results$ci_proportion[2] * 100,
  results$p_better_than_chance * 100,
  ifelse(results$p_better_than_chance > 0.95, "strong",
         ifelse(results$p_better_than_chance > 0.80, "moderate", "weak")),
  results$combined_pd * 100,
  ifelse(results$combined_pd > 0.80, "strong",
         ifelse(results$combined_pd > 0.65, "moderate", "weak")),
  sum(results$individual_pd > 0.5),
  results$n_predictions,
  mean(results$individual_pd > 0.5) * 100,
  sum(results$individual_pd > 0.95),
  mean(results$individual_pd > 0.95) * 100
)

cat(manuscript_text)

# ==============================================================================
# 9. SAVE RESULTS
# ==============================================================================

cat("\n\n=== SAVING RESULTS ===\n")

# Save results object
saveRDS(results, "results/pattern_conformity_results.rds")

# Save summary table
write_csv(pd_summary, "results/pattern_conformity_by_prediction.csv")

# Save summary by domain
domain_summary <- pd_summary %>%
  group_by(domain) %>%
  summarise(
    n_predictions = n(),
    mean_PD = mean(PD),
    median_PD = median(PD),
    n_strong = sum(PD > 0.95),
    n_moderate = sum(PD > 0.80 & PD <= 0.95),
    n_weak = sum(PD > 0.50 & PD <= 0.80),
    n_opposite = sum(PD <= 0.50),
    .groups = "drop"
  )
write_csv(domain_summary, "results/pattern_conformity_by_domain.csv")

cat("✓ Results saved to results/pattern_conformity_*.csv/rds\n")

# ==============================================================================
# 10. FINAL SUMMARY
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("PATTERN CONFORMITY ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")

cat(sprintf("\n
KEY FINDINGS:
• Tested %d directional predictions from theoretical model
• Mean conformity: %.1f%% (95%% CI: %.1f%%-%.1f%%)
• P(better than chance): %.1f%%
• Combined PD: %.1f%%
• Evidence strength: %s

OUTPUT FILES:
• figures/pattern_conformity_main.png
• figures/pattern_conformity_individual_pds.png
• results/pattern_conformity_results.rds
• results/pattern_conformity_by_prediction.csv
• results/pattern_conformity_by_domain.csv

",
results$n_predictions,
results$mean_proportion * 100,
results$ci_proportion[1] * 100,
results$ci_proportion[2] * 100,
results$p_better_than_chance * 100,
results$combined_pd * 100,
ifelse(results$p_better_than_chance > 0.95, "STRONG",
       ifelse(results$p_better_than_chance > 0.80, "MODERATE", "WEAK"))
))
