# ==============================================================================
# 14_visualize_between_within.R
# Visualize between-within decomposition results
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(ggdist)
  library(patchwork)
  library(ggtext)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("BETWEEN-WITHIN VISUALIZATION\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND RESULTS
# ==============================================================================

dir.create("figures/between_within", showWarnings = FALSE, recursive = TRUE)

df <- readRDS("results/between_within/df_between_within.rds")
comparison <- readRDS("results/between_within/f0_mean_a_comparison.rds")
model_comp <- readRDS("results/between_within/model_comparison_summary.rds")
fit_bw <- readRDS("results/between_within/fit_f0_mean_a_bw.rds")

cat("Data and results loaded\n\n")

# ==============================================================================
# 2. VARIANCE DECOMPOSITION PLOT
# ==============================================================================

cat("Creating variance decomposition plot...\n")

# Calculate ICC for each domain
pid5_domains <- c(
  "negative_affectivity",
  "detachment",
  "antagonism",
  "disinhibition",
  "psychoticism"
)

icc_data <- map_df(pid5_domains, function(domain) {
  var_between <- paste0("pid5_", domain, "_between")
  var_within <- paste0("pid5_", domain, "_within")

  v_between <- var(df[[var_between]], na.rm = TRUE)
  v_within <- var(df[[var_within]], na.rm = TRUE)

  icc <- v_between / (v_between + v_within)

  tibble(
    domain = str_to_title(str_replace_all(domain, "_", " ")),
    variance_between = v_between,
    variance_within = v_within,
    variance_total = v_between + v_within,
    icc = icc,
    pct_between = 100 * icc,
    pct_within = 100 * (1 - icc)
  )
})

# Reshape for plotting
icc_long <- icc_data %>%
  select(domain, pct_between, pct_within) %>%
  pivot_longer(
    cols = c(pct_between, pct_within),
    names_to = "component",
    values_to = "percentage"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = c("pct_between", "pct_within"),
      labels = c("Between-person (Trait)", "Within-person (State)")
    )
  )

p_variance <- ggplot(
  icc_long,
  aes(x = domain, y = percentage, fill = component)
) +
  geom_col(position = "stack", alpha = 0.8) +
  geom_text(
    aes(label = sprintf("%.0f%%", percentage)),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Between-person (Trait)" = "#2E86AB",
      "Within-person (State)" = "#A23B72"
    ),
    name = "Variance Component"
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Variance Decomposition: Trait vs State",
    subtitle = "Proportion of total PID-5 variance at each level",
    x = "PID-5 Domain",
    y = "Percentage of Total Variance (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/between_within/variance_decomposition.png",
  p_variance,
  width = 10,
  height = 6,
  dpi = 300
)

cat("✓ Saved: variance_decomposition.png\n")

# ==============================================================================
# 3. FOREST PLOT: BETWEEN VS WITHIN EFFECTS
# ==============================================================================

cat("Creating forest plot comparing between vs within effects...\n")

# Prepare data for forest plot
forest_data <- comparison %>%
  select(
    domain,
    starts_with("between_main"),
    starts_with("within_main"),
    starts_with("between_stress"),
    starts_with("within_stress"),
    starts_with("between_recovery"),
    starts_with("within_recovery")
  ) %>%
  pivot_longer(
    cols = -domain,
    names_to = c("level", "effect", "stat"),
    names_pattern = "(between|within)_(.+)_(mean|ci_lower|ci_upper|credible)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    domain = str_to_title(str_replace_all(domain, "_", " ")),
    level = str_to_title(level),
    effect = str_to_title(str_replace_all(effect, "_", " ")),
    effect = factor(effect, levels = c("Main", "Stress", "Recovery"))
  )

# Main effects only for clarity
forest_main <- forest_data %>%
  filter(effect == "Main", !is.na(mean)) %>%
  mutate(
    level = factor(level, levels = c("Within", "Between")),
    credible = as.logical(credible)
  )

p_forest <- ggplot(
  forest_main,
  aes(x = mean, y = domain, color = level, shape = credible)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_linerange(
    aes(xmin = ci_lower, xmax = ci_upper),
    position = position_dodge(width = 0.5),
    linewidth = 1
  ) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  scale_color_manual(
    values = c("Between" = "#2E86AB", "Within" = "#A23B72"),
    name = "Level"
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    name = "Credible",
    labels = c(
      "FALSE" = "No (95% CI includes 0)",
      "TRUE" = "Yes (95% CI excludes 0)"
    )
  ) +
  labs(
    title = "Between-Person vs Within-Person Main Effects on F0",
    subtitle = "Effect sizes with 95% credible intervals",
    x = "Effect Size (β)",
    y = "PID-5 Domain"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  "figures/between_within/forest_plot_main_effects.png",
  p_forest,
  width = 10,
  height = 7,
  dpi = 300
)

cat("✓ Saved: forest_plot_main_effects.png\n")

# ==============================================================================
# 4. INTERACTION PLOTS FOR CREDIBLE EFFECTS
# ==============================================================================

cat("Creating interaction plots for credible effects...\n")

# Initialize list to track created plots
interaction_plots <- list()

# Identify credible within-person interactions
credible_within <- comparison %>%
  filter(within_stress_credible == TRUE | within_recovery_credible == TRUE) %>%
  select(domain, within_stress_credible, within_recovery_credible)

if (nrow(credible_within) > 0) {
  # Function to create interaction plot
  plot_interaction <- function(domain_name, effect_type) {
    # Get variable names
    domain_clean <- tolower(str_replace_all(domain_name, " ", "_"))
    var_within <- paste0("pid5_", domain_clean, "_within_c")

    # Create prediction grid - FIXED: use correct timepoint levels
    pred_grid <- expand.grid(
      timepoint = factor(
        c("baseline", "pre_exam", "post_exam"),
        levels = c("baseline", "pre_exam", "post_exam")
      ),
      state_level = c(-1, 0, 1) # Low, Mean, High state
    )

    pred_grid[[var_within]] <- pred_grid$state_level

    # Add other predictors at mean (0 for centered)
    other_domains <- setdiff(pid5_domains, domain_clean)
    for (d in other_domains) {
      pred_grid[[paste0("pid5_", d, "_between_c")]] <- 0
      pred_grid[[paste0("pid5_", d, "_within_c")]] <- 0
    }

    # This domain's between component at mean
    pred_grid[[paste0("pid5_", domain_clean, "_between_c")]] <- 0

    # Add contrast codes
    pred_grid <- pred_grid %>%
      mutate(
        c1_stress = case_when(
          timepoint == "baseline" ~ -0.5,
          timepoint == "pre_exam" ~ 0.5,
          timepoint == "post_exam" ~ 0
        ),
        c2_recovery = case_when(
          timepoint == "baseline" ~ 0,
          timepoint == "pre_exam" ~ -0.5,
          timepoint == "post_exam" ~ 0.5
        ),
        state_label = factor(
          state_level,
          levels = c(-1, 0, 1),
          labels = c("Low (-1 SD)", "Mean", "High (+1 SD)")
        ),
        timepoint_label = factor(
          timepoint,
          levels = c("baseline", "pre_exam", "post_exam"),
          labels = c("Baseline", "Pre-Exam", "Post-Exam")
        )
      )

    # Get predictions
    preds <- predict(fit_bw, newdata = pred_grid, summary = TRUE) %>%
      as_tibble() %>%
      bind_cols(pred_grid)

    # Plot
    p <- ggplot(
      preds,
      aes(
        x = timepoint_label,
        y = Estimate,
        color = state_label,
        group = state_label
      )
    ) +
      geom_ribbon(
        aes(ymin = Q2.5, ymax = Q97.5, fill = state_label),
        alpha = 0.2,
        color = NA
      ) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_color_viridis_d(option = "plasma", end = 0.9) +
      scale_fill_viridis_d(option = "plasma", end = 0.9) +
      labs(
        title = paste("State-dependent trajectory:", domain_name),
        subtitle = "How within-person fluctuations predict F0 across stress",
        x = "Timepoint",
        y = "Predicted F0 (Hz)",
        color = "State Level",
        fill = "State Level"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )

    return(p)
  }

  # Create plots for each credible effect
  for (i in 1:nrow(credible_within)) {
    domain_name <- str_to_title(str_replace_all(
      credible_within$domain[i],
      "_",
      " "
    ))

    if (
      credible_within$within_stress_credible[i] |
        credible_within$within_recovery_credible[i]
    ) {
      tryCatch(
        {
          p <- plot_interaction(domain_name, "combined")
          interaction_plots[[domain_name]] <- p

          filename <- paste0(
            "figures/between_within/interaction_",
            tolower(str_replace_all(domain_name, " ", "_")),
            ".png"
          )
          ggsave(filename, p, width = 10, height = 6, dpi = 300)
          cat("✓ Saved:", basename(filename), "\n")
        },
        error = function(e) {
          cat("✗ Error creating plot for", domain_name, ":", e$message, "\n")
        }
      )
    }
  }
} else {
  cat("No credible within-person interactions to plot\n")
}

# ==============================================================================
# 5. MODEL COMPARISON R² PLOT
# ==============================================================================

cat("Creating model comparison summary figure...\n")

# Extract key metrics
r2_data <- model_comp %>%
  filter(metric %in% c("R2_trait", "R2_trait_state")) %>%
  mutate(
    model = c("Trait-only", "Trait + State"),
    model = factor(model, levels = c("Trait-only", "Trait + State"))
  )

p_r2 <- ggplot(r2_data, aes(x = model, y = value)) +
  geom_col(aes(fill = model), alpha = 0.7, width = 0.6) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2,
    linewidth = 1
  ) +
  geom_text(
    aes(label = sprintf("R² = %.3f", value)),
    vjust = -2,
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(
    values = c("Trait-only" = "#2E86AB", "Trait + State" = "#A23B72")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Model Comparison: Variance Explained",
    subtitle = "Bayesian R² with 95% credible intervals",
    x = NULL,
    y = "Variance Explained (R²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(
  "figures/between_within/model_comparison_r2.png",
  p_r2,
  width = 8,
  height = 6,
  dpi = 300
)

cat("✓ Saved: model_comparison_r2.png\n")

# ==============================================================================
# 6. COMBINED SUMMARY FIGURE
# ==============================================================================

cat("Creating combined summary figure...\n")

# Combine variance decomposition + forest plot + R² comparison
p_combined <- (p_variance / p_forest / p_r2) +
  plot_layout(heights = c(1, 1.2, 1)) +
  plot_annotation(
    title = "Between-Person vs Within-Person Analysis: Complete Summary",
    subtitle = "Decomposition of trait and state contributions to voice-personality associations",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave(
  "figures/between_within/combined_summary.png",
  p_combined,
  width = 14,
  height = 16,
  dpi = 300
)

cat("✓ Saved: combined_summary.png\n")

# ==============================================================================
# 7. PUBLICATION TABLE
# ==============================================================================

cat("Creating summary table...\n")

# Clean comparison table for publication
pub_table <- comparison %>%
  mutate(
    domain = str_to_title(str_replace_all(domain, "_", " ")),

    # Between main
    between_main = sprintf(
      "%.2f [%.2f, %.2f]%s",
      between_main_mean,
      between_main_ci_lower,
      between_main_ci_upper,
      ifelse(between_main_credible, "*", "")
    ),

    # Within main
    within_main = sprintf(
      "%.2f [%.2f, %.2f]%s",
      within_main_mean,
      within_main_ci_lower,
      within_main_ci_upper,
      ifelse(within_main_credible, "*", "")
    ),

    # Between stress
    between_stress = sprintf(
      "%.2f [%.2f, %.2f]%s",
      between_stress_mean,
      between_stress_ci_lower,
      between_stress_ci_upper,
      ifelse(between_stress_credible, "*", "")
    ),

    # Within stress
    within_stress = sprintf(
      "%.2f [%.2f, %.2f]%s",
      within_stress_mean,
      within_stress_ci_lower,
      within_stress_ci_upper,
      ifelse(within_stress_credible, "*", "")
    )
  ) %>%
  select(
    Domain = domain,
    `Between Main` = between_main,
    `Within Main` = within_main,
    `Between × Stress` = between_stress,
    `Within × Stress` = within_stress
  )

rio::export(pub_table, "results/between_within/publication_table.csv")
cat("✓ Saved: publication_table.csv\n\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat(rep("=", 70), "\n")
cat("VISUALIZATION COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Figures created:\n")
cat("  1. variance_decomposition.png - ICC for each domain\n")
cat("  2. forest_plot_main_effects.png - Between vs Within comparison\n")
cat("  3. model_comparison_r2.png - Model fit comparison\n")
cat("  4. combined_summary.png - All panels together\n")

if (length(interaction_plots) > 0) {
  cat(sprintf(
    "  5. interaction_*.png - %d state × timepoint plots\n",
    length(interaction_plots)
  ))
}

cat("\nTables created:\n")
cat("  1. publication_table.csv - Formatted parameter estimates\n\n")

cat("All outputs in: figures/between_within/ and results/between_within/\n\n")
