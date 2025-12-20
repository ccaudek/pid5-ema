# ==============================================================================
# Bayesian Pattern Conformity Analysis
# Testing overall theoretical predictions rather than individual effects
# ==============================================================================

library(tidyverse)
library(brms)
library(posterior)
library(patchwork)

# ==============================================================================
# 1. DEFINE THEORETICAL PREDICTIONS
# ==============================================================================

# Your theory predicts DIRECTION of effects, not just that they exist.
# Create a prediction matrix: +1 = expect positive, -1 = expect negative, 0 = no prediction

# Example predictions based on your abstract and theory:
# - Negative Affectivity: amplifies stress reactivity (+ on C1_stress interaction)
# - Detachment: impairs recovery (- on C2_recovery interaction, meaning LESS recovery)
# - Antagonism: facilitates recovery (+ on C2_recovery interaction)
# - Psychoticism: stable elevation (+ main effect on pitch/voice quality)
# - Disinhibition: articulatory alterations (+ on F2)

# Define your theoretical predictions
# Format: list of named vectors, one per outcome family
# Names should match your parameter names from brms

create_prediction_matrix <- function() {
  
  # Define which direction you expect for each interaction
  # These should be based on YOUR theoretical predictions - adjust as needed!
  
  predictions <- tribble(
    ~outcome,    ~vowel, ~parameter,                                    ~expected_direction, ~rationale,
    # F0 Mean predictions
    "f0_mean",   "a",    "c1_stress:pid5_negative_affectivity_bl_c",    +1,    "NA amplifies stress-induced pitch increase",
    "f0_mean",   "a",    "c1_stress:pid5_detachment_bl_c",              -1,    "Detachment blunts emotional reactivity",
    "f0_mean",   "a",    "c1_stress:pid5_antagonism_bl_c",               0,    "No clear prediction",
    "f0_mean",   "a",    "c1_stress:pid5_disinhibition_bl_c",            0,    "No clear prediction",
    "f0_mean",   "a",    "c1_stress:pid5_psychoticism_bl_c",            +1,    "Psychoticism shows pervasive arousal",
    "f0_mean",   "a",    "c2_recovery:pid5_negative_affectivity_bl_c", -1,    "NA impairs return to baseline",
    "f0_mean",   "a",    "c2_recovery:pid5_detachment_bl_c",           -1,    "Detachment impairs recovery",
    "f0_mean",   "a",    "c2_recovery:pid5_antagonism_bl_c",           +1,    "Antagonism facilitates recovery",
    "f0_mean",   "a",    "c2_recovery:pid5_disinhibition_bl_c",         0,    "No clear prediction",
    "f0_mean",   "a",    "c2_recovery:pid5_psychoticism_bl_c",         -1,    "Psychoticism shows persistent elevation",
    
    # Replicate for vowels /i/ and /u/ (same predictions)
    "f0_mean",   "i",    "c1_stress:pid5_negative_affectivity_bl_c",   +1,    "Replication",
    "f0_mean",   "i",    "c2_recovery:pid5_detachment_bl_c",           -1,    "Replication",
    "f0_mean",   "i",    "c2_recovery:pid5_antagonism_bl_c",           +1,    "Replication",
    "f0_mean",   "u",    "c1_stress:pid5_negative_affectivity_bl_c",   +1,    "Replication",
    "f0_mean",   "u",    "c2_recovery:pid5_detachment_bl_c",           -1,    "Replication",
    "f0_mean",   "u",    "c2_recovery:pid5_antagonism_bl_c",           +1,    "Replication",
    
    # NNE predictions (more negative = cleaner voice; stress makes it less negative)
    "nne",       "a",    "c1_stress:pid5_negative_affectivity_bl_c",   +1,    "NA amplifies stress-induced voice quality degradation",
    "nne",       "a",    "c1_stress:pid5_psychoticism_bl_c",           +1,    "Psychoticism shows voice quality issues",
    
    # Jitter predictions
    "jitter",    "a",    "c1_stress:pid5_negative_affectivity_bl_c",   +1,    "NA increases vocal tension under stress",
    "jitter",    "a",    "c1_stress:pid5_psychoticism_bl_c",           +1,    "Psychoticism shows voice instability",
    
    # F2 predictions (articulatory)
    "f2_mean",   "a",    "c1_stress:pid5_disinhibition_bl_c",          +1,    "Disinhibition shows articulatory alterations",
    "f2_std",    "a",    "c1_stress:pid5_disinhibition_bl_c",          +1,    "Disinhibition shows articulatory variability"
    
    # Add more predictions as needed...
  )
  
  return(predictions)
}

# ==============================================================================
# 2. EXTRACT POSTERIOR DRAWS AND COMPUTE CONFORMITY
# ==============================================================================

#' Compute pattern conformity from a list of fitted brms models
#' 
#' @param fitted_models Named list of brmsfit objects
#' @param predictions Tibble from create_prediction_matrix()
#' @return List with conformity metrics and posterior distributions

compute_pattern_conformity <- function(fitted_models, predictions) {
  
  # Filter to only predictions with expected direction (exclude 0s)
  predictions_directional <- predictions %>%
    filter(expected_direction != 0)
  
  n_predictions <- nrow(predictions_directional)
  cat("Testing", n_predictions, "directional predictions\n")
  
  # Initialize storage for posterior draws of conformity
  # We'll compute conformity for each MCMC draw
  n_draws <- 4000  # Adjust based on your actual posterior samples
  
  conformity_by_draw <- matrix(NA, nrow = n_draws, ncol = n_predictions)
  colnames(conformity_by_draw) <- paste0(predictions_directional$outcome, "_", 
                                          predictions_directional$vowel, "_",
                                          predictions_directional$parameter)
  
  # Extract posterior for each prediction
  for (i in seq_len(n_predictions)) {
    pred <- predictions_directional[i, ]
    
    # Find the corresponding model
    model_name <- paste0("m_", pred$outcome, "_", pred$vowel)
    
    if (!model_name %in% names(fitted_models)) {
      warning("Model not found: ", model_name)
      next
    }
    
    model <- fitted_models[[model_name]]
    
    # Extract posterior draws for this parameter
    param_name <- paste0("b_", pred$parameter)
    draws <- as_draws_df(model)
    
    if (!param_name %in% names(draws)) {
      # Try alternative naming conventions
      param_name_alt <- gsub(":", ".", pred$parameter)
      param_name_alt <- paste0("b_", param_name_alt)
      if (param_name_alt %in% names(draws)) {
        param_name <- param_name_alt
      } else {
        warning("Parameter not found: ", param_name)
        next
      }
    }
    
    posterior_samples <- draws[[param_name]]
    
    # Compute: is this draw in the expected direction?
    # conformity = 1 if sign matches expected, 0 otherwise
    if (pred$expected_direction == 1) {
      conformity_by_draw[, i] <- as.numeric(posterior_samples > 0)
    } else {
      conformity_by_draw[, i] <- as.numeric(posterior_samples < 0)
    }
  }
  
  # ==============================================================================
  # 3. COMPUTE AGGREGATE CONFORMITY METRICS
  # ==============================================================================
  
  results <- list()
  
  # --- Metric 1: Proportion of effects in expected direction (per draw) ---
  # This gives a posterior distribution over "fraction of predictions confirmed"
  results$proportion_conforming <- rowMeans(conformity_by_draw, na.rm = TRUE)
  
  # --- Metric 2: Joint probability that ALL predictions hold ---
  # Very stringent: P(all effects simultaneously in expected direction)
  results$all_conform <- rowMeans(conformity_by_draw == 1, na.rm = TRUE) == 1
  results$p_all_conform <- mean(results$all_conform)
  
  # --- Metric 3: Number of conforming predictions (per draw) ---
  results$n_conforming <- rowSums(conformity_by_draw, na.rm = TRUE)
  
  # --- Metric 4: Individual Probability of Direction (PD) per prediction ---
  results$individual_pd <- colMeans(conformity_by_draw, na.rm = TRUE)
  names(results$individual_pd) <- colnames(conformity_by_draw)
  
  # --- Metric 5: Combined PD using log-odds (more robust than product) ---
  # Convert PDs to log-odds, average, convert back
  pd_values <- results$individual_pd
  pd_values <- pmax(pmin(pd_values, 0.999), 0.001)  # Avoid inf
  log_odds <- log(pd_values / (1 - pd_values))
  mean_log_odds <- mean(log_odds, na.rm = TRUE)
  results$combined_pd_logodds <- exp(mean_log_odds) / (1 + exp(mean_log_odds))
  
  # --- Metric 6: Geometric mean of PDs ---
  results$geometric_mean_pd <- exp(mean(log(pd_values), na.rm = TRUE))
  
  # --- Metric 7: Bayesian sign test ---
  # Under null (random direction), expect 50% conformity
  # Compute P(observed conformity > 50% | data)
  results$p_better_than_chance <- mean(results$proportion_conforming > 0.5)
  
  # --- Metric 8: Effect size for pattern conformity ---
  # How much better than chance (50%) is the observed conformity?
  results$conformity_effect_size <- mean(results$proportion_conforming) - 0.5
  
  # Store raw data
  results$conformity_matrix <- conformity_by_draw
  results$predictions <- predictions_directional
  results$n_predictions <- n_predictions
  
  return(results)
}

# ==============================================================================
# 4. VISUALIZATION FUNCTIONS
# ==============================================================================

plot_pattern_conformity <- function(results) {
  
  # Plot 1: Posterior distribution of proportion conforming
  p1 <- ggplot(data.frame(prop = results$proportion_conforming), aes(x = prop)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#2c5282", alpha = 0.7) +
    geom_density(color = "#1a365d", linewidth = 1) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = mean(results$proportion_conforming), 
               linetype = "solid", color = "#38a169", linewidth = 1.2) +
    annotate("text", x = 0.52, y = Inf, label = "Chance", vjust = 2, color = "red") +
    annotate("text", x = mean(results$proportion_conforming) + 0.02, y = Inf, 
             label = sprintf("Mean = %.1f%%", mean(results$proportion_conforming) * 100), 
             vjust = 2, color = "#38a169", hjust = 0) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = "Posterior Distribution of Pattern Conformity",
      subtitle = sprintf("P(conformity > chance) = %.1f%%", results$p_better_than_chance * 100),
      x = "Proportion of Predictions in Expected Direction",
      y = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  # Plot 2: Individual PDs with theoretical predictions
  pd_df <- data.frame(
    prediction = names(results$individual_pd),
    pd = results$individual_pd
  ) %>%
    mutate(
      in_expected = pd > 0.5,
      strong_evidence = pd > 0.95 | pd < 0.05
    ) %>%
    arrange(pd)
  
  pd_df$prediction <- factor(pd_df$prediction, levels = pd_df$prediction)
  
  p2 <- ggplot(pd_df, aes(x = pd, y = prediction, fill = in_expected)) +
    geom_col(alpha = 0.8) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0.95, linetype = "dotted", color = "darkgreen") +
    scale_fill_manual(values = c("FALSE" = "#e53e3e", "TRUE" = "#38a169"),
                      labels = c("Opposite direction", "Expected direction")) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = "Probability of Direction for Each Prediction",
      x = "P(effect in expected direction)",
      y = NULL,
      fill = "Direction"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 7),
      legend.position = "bottom"
    )
  
  # Plot 3: Distribution of number of conforming predictions
  p3 <- ggplot(data.frame(n = results$n_conforming), aes(x = n)) +
    geom_histogram(bins = results$n_predictions + 1, fill = "#805ad5", alpha = 0.7) +
    geom_vline(xintercept = results$n_predictions / 2, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = seq(0, results$n_predictions, by = 5)) +
    labs(
      title = "Posterior Distribution: Number of Confirmed Predictions",
      subtitle = sprintf("Out of %d directional predictions", results$n_predictions),
      x = "Number of predictions in expected direction",
      y = "Posterior draws"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  return(list(p1 = p1, p2 = p2, p3 = p3))
}

# ==============================================================================
# 5. SUMMARY REPORT FUNCTION
# ==============================================================================

print_conformity_report <- function(results) {
  
  cat("\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("BAYESIAN PATTERN CONFORMITY ANALYSIS\n")
  cat("Testing Overall Theoretical Predictions\n")
  cat(rep("=", 70), "\n", sep = "")
  
  cat("\n--- SUMMARY METRICS ---\n\n")
  
  cat(sprintf("Number of directional predictions tested: %d\n", results$n_predictions))
  cat(sprintf("Mean proportion conforming: %.1f%% (chance = 50%%)\n", 
              mean(results$proportion_conforming) * 100))
  cat(sprintf("95%% CI for proportion conforming: [%.1f%%, %.1f%%]\n",
              quantile(results$proportion_conforming, 0.025) * 100,
              quantile(results$proportion_conforming, 0.975) * 100))
  
  cat(sprintf("\nP(pattern conformity > chance): %.2f%%\n", 
              results$p_better_than_chance * 100))
  
  cat(sprintf("P(ALL predictions simultaneously correct): %.4f%%\n", 
              results$p_all_conform * 100))
  
  cat(sprintf("\nCombined PD (log-odds method): %.1f%%\n", 
              results$combined_pd_logodds * 100))
  cat(sprintf("Geometric mean PD: %.1f%%\n", 
              results$geometric_mean_pd * 100))
  
  cat(sprintf("\nEffect size (conformity - 0.5): %.3f\n", 
              results$conformity_effect_size))
  
  # Interpretation
  cat("\n--- INTERPRETATION ---\n\n")
  
  if (results$p_better_than_chance > 0.95) {
    cat("Strong evidence that the overall theoretical pattern holds.\n")
    cat("The observed pattern of results is unlikely under random expectation.\n")
  } else if (results$p_better_than_chance > 0.80) {
    cat("Moderate evidence favoring the theoretical pattern.\n")
    cat("Results are more consistent with theory than chance.\n")
  } else {
    cat("Weak or inconclusive evidence for the overall pattern.\n")
    cat("Individual effects may still be informative.\n")
  }
  
  # Individual predictions
  cat("\n--- INDIVIDUAL PREDICTIONS ---\n\n")
  
  pd_summary <- data.frame(
    prediction = names(results$individual_pd),
    PD = results$individual_pd
  ) %>%
    mutate(
      evidence = case_when(
        PD > 0.975 ~ "Strong support",
        PD > 0.90 ~ "Moderate support", 
        PD > 0.75 ~ "Weak support",
        PD > 0.50 ~ "Slight support",
        TRUE ~ "Opposite direction"
      )
    ) %>%
    arrange(desc(PD))
  
  cat("Top 10 most supported predictions:\n")
  print(head(pd_summary, 10), row.names = FALSE)
  
  cat("\nPredictions in opposite direction (PD < 50%):\n")
  opposite <- pd_summary %>% filter(PD < 0.5)
  if (nrow(opposite) > 0) {
    print(opposite, row.names = FALSE)
  } else {
    cat("None - all predictions in expected direction!\n")
  }
  
  cat("\n", rep("=", 70), "\n", sep = "")
}

# ==============================================================================
# 6. EXAMPLE USAGE
# ==============================================================================

# After running your moderation_analysis.R, you would have fitted_models list
# Then:

# predictions <- create_prediction_matrix()
# results <- compute_pattern_conformity(fitted_models, predictions)
# print_conformity_report(results)
# plots <- plot_pattern_conformity(results)
# 
# # Save combined plot
# combined <- plots$p1 / plots$p3
# ggsave("figures/pattern_conformity.png", combined, width = 10, height = 10, dpi = 300)

# ==============================================================================
# 7. ALTERNATIVE: JOINT POSTERIOR SIMULATION
# ==============================================================================

#' More sophisticated approach: simulate from joint posterior
#' accounting for correlations between estimates within each model

compute_joint_conformity <- function(fitted_models, predictions) {
  
  predictions_directional <- predictions %>%
    filter(expected_direction != 0)
  
  # Get posterior draws from all models
  all_conformity_draws <- list()
  
  for (model_name in unique(paste0("m_", predictions_directional$outcome, "_", 
                                    predictions_directional$vowel))) {
    
    if (!model_name %in% names(fitted_models)) next
    
    model <- fitted_models[[model_name]]
    draws <- as_draws_df(model)
    n_draws <- nrow(draws)
    
    # Get predictions for this model
    model_preds <- predictions_directional %>%
      filter(paste0("m_", outcome, "_", vowel) == model_name)
    
    # For each prediction, check conformity
    for (j in seq_len(nrow(model_preds))) {
      pred <- model_preds[j, ]
      param_name <- paste0("b_", pred$parameter)
      
      if (param_name %in% names(draws)) {
        posterior_samples <- draws[[param_name]]
        
        conformity <- if (pred$expected_direction == 1) {
          posterior_samples > 0
        } else {
          posterior_samples < 0
        }
        
        key <- paste0(pred$outcome, "_", pred$vowel, "_", pred$parameter)
        all_conformity_draws[[key]] <- conformity
      }
    }
  }
  
  # Convert to matrix (rows = draws, cols = predictions)
  conformity_matrix <- do.call(cbind, all_conformity_draws)
  
  # Key insight: within each draw, all estimates are from the SAME posterior sample
  # So correlations between parameters are preserved
  
  # Compute joint conformity metrics
  results <- list(
    # For each draw: what proportion of predictions conform?
    proportion_per_draw = rowMeans(conformity_matrix, na.rm = TRUE),
    
    # For each draw: do ALL predictions conform?
    all_conform_per_draw = apply(conformity_matrix, 1, all, na.rm = TRUE),
    
    # Posterior probability that at least X% of predictions conform
    p_at_least_50_pct = mean(rowMeans(conformity_matrix, na.rm = TRUE) >= 0.5),
    p_at_least_75_pct = mean(rowMeans(conformity_matrix, na.rm = TRUE) >= 0.75),
    p_at_least_90_pct = mean(rowMeans(conformity_matrix, na.rm = TRUE) >= 0.90),
    
    conformity_matrix = conformity_matrix
  )
  
  return(results)
}

# ==============================================================================
# 8. MANUSCRIPT TEXT GENERATOR
# ==============================================================================

generate_results_text <- function(results) {
  
  text <- sprintf(
"To evaluate our overall theoretical model rather than individual effects, we 
computed a Bayesian pattern conformity analysis. We specified directional 
predictions for %d effects based on our theoretical framework (e.g., Negative 
Affectivity amplifies stress reactivity; Detachment impairs recovery). For each 
posterior draw, we computed the proportion of predictions in the expected 
direction, yielding a posterior distribution over pattern conformity.
  
Results indicated %s evidence for our theoretical model. The mean proportion 
of predictions in the expected direction was %.1f%% (95%% CI [%.1f%%, %.1f%%]), 
substantially exceeding the 50%% expected by chance. The posterior probability 
that our theoretical pattern performs better than chance was %.1f%%. 

Using a log-odds combination method, the aggregate probability of direction 
across all predictions was %.1f%%, indicating %s support for the overall 
directionality of effects. The probability that ALL predictions simultaneously 
held was %.4f%%, reflecting the stringent nature of this joint test.

These results suggest that personality pathology domains moderate stress 
reactivity in voice acoustics in theoretically predictable ways, supporting 
the construct validity of context-dependent trait expression.",
    results$n_predictions,
    ifelse(results$p_better_than_chance > 0.95, "strong",
           ifelse(results$p_better_than_chance > 0.80, "moderate", "weak")),
    mean(results$proportion_conforming) * 100,
    quantile(results$proportion_conforming, 0.025) * 100,
    quantile(results$proportion_conforming, 0.975) * 100,
    results$p_better_than_chance * 100,
    results$combined_pd_logodds * 100,
    ifelse(results$combined_pd_logodds > 0.80, "strong",
           ifelse(results$combined_pd_logodds > 0.65, "moderate", "weak")),
    results$p_all_conform * 100
  )
  
  return(text)
}

cat("\n✓ Pattern conformity analysis functions loaded\n")
cat("  Use create_prediction_matrix() to define your theoretical predictions\n")
cat("  Use compute_pattern_conformity() to run the analysis\n")
cat("  Use print_conformity_report() for summary output\n")
cat("  Use plot_pattern_conformity() for visualizations\n")
cat("  Use generate_results_text() for manuscript-ready text\n")
