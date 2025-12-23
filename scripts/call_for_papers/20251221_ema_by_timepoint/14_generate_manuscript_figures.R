# ==============================================================================
# 14_generate_manuscript_figures.R
# Generate publication-ready figures
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggdist)
  library(ggtext)
  library(scales)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("GENERATING MANUSCRIPT FIGURES\n")
cat(rep("=", 70), "\n\n")

dir.create("figures", showWarnings = FALSE)

# Load data
stress_effects <- read_csv(
  "results/stress_main_effects.csv",
  show_col_types = FALSE
)
df <- readRDS("results/df_analysis.rds")

# ==============================================================================
# FIGURE 1: Study Design
# ==============================================================================

cat("Creating Figure 1: Study Design...\n")

# This is a conceptual diagram - create manually or with ggplot annotations
# For now, I'll create a timeline visualization

timeline_data <- tibble(
  timepoint = factor(
    c("Baseline", "Pre-Exam", "Post-Exam"),
    levels = c("Baseline", "Pre-Exam", "Post-Exam")
  ),
  day = c(-60, 0, 1),
  label = c(
    "Baseline\n(~2 months before exam)",
    "Pre-Exam\n(day before exam)",
    "Post-Exam\n(day after exam)"
  ),
  ema_period = c("EMA Period", "EMA Period", "End")
)

p_timeline <- ggplot(timeline_data, aes(x = day, y = 1)) +
  # Timeline arrow
  geom_segment(
    aes(x = -70, xend = 10, y = 1, yend = 1),
    arrow = arrow(length = unit(0.3, "cm")),
    linewidth = 1,
    color = "gray40"
  ) +

  # Timepoints
  geom_point(size = 8, color = "#2E86AB") +
  geom_text(aes(label = timepoint), vjust = -2, size = 4, fontface = "bold") +
  geom_text(aes(label = label), vjust = 3, size = 3.5, lineheight = 0.9) +

  # EMA period
  annotate(
    "rect",
    xmin = -65,
    xmax = 5,
    ymin = 0.7,
    ymax = 0.8,
    fill = "#A23B72",
    alpha = 0.3
  ) +
  annotate(
    "text",
    x = -30,
    y = 0.6,
    label = "EMA Period (~2.5 months, ~30 assessments)",
    size = 3.5,
    fontface = "italic"
  ) +

  # Voice recordings
  geom_point(
    aes(y = 1.2),
    size = 6,
    shape = 21,
    fill = "#D32F2F",
    color = "white"
  ) +
  annotate(
    "text",
    x = -30,
    y = 1.35,
    label = "Voice Recordings (3 timepoints)",
    size = 3.5,
    fontface = "italic"
  ) +

  scale_x_continuous(expand = c(0.05, 0.05)) +
  coord_cartesian(ylim = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(
  "figures/figure1_study_design.png",
  p_timeline,
  width = 12,
  height = 4,
  dpi = 300
)

cat("✓ Figure 1 saved\n")

# ==============================================================================
# FIGURE 2: Stress Main Effects
# ==============================================================================

cat("Creating Figure 2: Stress Main Effects...\n")

# Panel A: Forest plot
plot_data_stress <- stress_effects %>%
  filter(!is.na(stress_beta)) %>%
  mutate(
    outcome_label = str_replace(outcome, " /", "\n/"),
    credible = stress_credible
  )

p2a <- ggplot(
  plot_data_stress,
  aes(x = stress_d, y = fct_reorder(outcome_label, stress_d))
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray70") +
  geom_linerange(
    aes(
      xmin = stress_ci_lower / sd(df$f0_mean_a, na.rm = TRUE),
      xmax = stress_ci_upper / sd(df$f0_mean_a, na.rm = TRUE),
      color = credible
    ),
    linewidth = 1.2
  ) +
  geom_point(aes(color = credible, shape = credible), size = 4) +
  scale_color_manual(
    values = c("TRUE" = "#D32F2F", "FALSE" = "gray60"),
    name = "Credible",
    labels = c("No", "Yes")
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    name = "Credible",
    labels = c("No", "Yes")
  ) +
  labs(
    title = "A. Stress Effects (Baseline → Pre-Exam)",
    x = "Cohen's d",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10)
  )

# Panel B: Effect size magnitudes
plot_magnitude <- stress_effects %>%
  filter(!is.na(stress_d)) %>%
  mutate(
    effect_size = abs(stress_d),
    magnitude = case_when(
      effect_size < 0.2 ~ "Negligible",
      effect_size < 0.5 ~ "Small",
      effect_size < 0.8 ~ "Medium",
      TRUE ~ "Large"
    ),
    magnitude = factor(
      magnitude,
      levels = c("Negligible", "Small", "Medium", "Large")
    )
  )

p2b <- ggplot(plot_magnitude, aes(x = magnitude, fill = magnitude)) +
  geom_bar(alpha = 0.8, color = "black", linewidth = 0.5) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  labs(
    title = "B. Distribution of Effect Sizes",
    x = "Effect Size Category",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Combine panels
p_figure2 <- p2a +
  p2b +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = "Figure 2. Main Effects of Academic Stress on Acoustic Features",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave(
  "figures/figure2_stress_main_effects.png",
  p_figure2,
  width = 12,
  height = 6,
  dpi = 300
)

cat("✓ Figure 2 saved\n")

# ==============================================================================
# FIGURE 3: Personality Moderation (Key Interactions)
# ==============================================================================

cat("Creating Figure 3: Personality Moderation...\n")

# Load models for predictions
m_f0_i <- readRDS("models/m_f0_mean_i.rds")

# Create prediction grid for Negative Affectivity × Stress
pred_grid <- expand.grid(
  pid5_negative_affectivity_c = seq(-2, 2, 0.5),
  c1_stress = c(-0.5, 0.5),
  c2_recovery = 0,
  # Set other predictors to 0
  pid5_detachment_c = 0,
  pid5_antagonism_c = 0,
  pid5_disinhibition_c = 0,
  pid5_psychoticism_c = 0
) %>%
  mutate(
    timepoint = factor(
      ifelse(c1_stress == -0.5, "Baseline", "Pre-Exam"),
      levels = c("Baseline", "Pre-Exam")
    ),
    trait_level = case_when(
      pid5_negative_affectivity_c < -1 ~ "Low (-1 SD)",
      abs(pid5_negative_affectivity_c) <= 1 ~ "Average",
      pid5_negative_affectivity_c > 1 ~ "High (+1 SD)"
    ),
    trait_level = factor(
      trait_level,
      levels = c("Low (-1 SD)", "Average", "High (+1 SD)")
    )
  )

# Get predictions
preds <- predict(
  m_f0_i,
  newdata = pred_grid,
  summary = TRUE,
  re_formula = NA
) %>%
  as_tibble() %>%
  bind_cols(pred_grid)

# Plot
p_figure3 <- ggplot(
  preds %>% filter(trait_level != "Average"),
  aes(
    x = timepoint,
    y = Estimate,
    group = interaction(pid5_negative_affectivity_c, trait_level),
    color = trait_level
  )
) +
  geom_ribbon(
    aes(ymin = Q2.5, ymax = Q97.5, fill = trait_level),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("Low (-1 SD)" = "#2E86AB", "High (+1 SD)" = "#D32F2F"),
    name = "Negative Affectivity"
  ) +
  scale_fill_manual(
    values = c("Low (-1 SD)" = "#2E86AB", "High (+1 SD)" = "#D32F2F"),
    name = "Negative Affectivity"
  ) +
  labs(
    title = "Figure 3. Negative Affectivity Moderates Stress-Related Pitch Increase",
    subtitle = "F0 Mean /i/ across baseline and pre-exam timepoints",
    x = "Timepoint",
    y = "Predicted F0 (Hz)",
    caption = "Shaded areas represent 95% credible intervals. Other PID-5 dimensions held at mean."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/figure3_personality_moderation.png",
  p_figure3,
  width = 10,
  height = 6,
  dpi = 300
)

cat("✓ Figure 3 saved\n")

# ==============================================================================
# FIGURE 4: EMA vs Baseline Comparison
# ==============================================================================

cat("Creating Figure 4: EMA vs Baseline Comparison...\n")

comparison <- readRDS("results/ema_vs_baseline_comparison.rds")

plot_data_comp <- comparison$r2_comparison %>%
  pivot_longer(
    cols = c(R2_EMA, R2_Baseline),
    names_to = "Model",
    values_to = "R2"
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = c("R2_EMA", "R2_Baseline"),
      labels = c(
        "EMA\n(15 items, ~30 assessments)",
        "Baseline\n(220 items, single)"
      )
    )
  )

p_figure4 <- ggplot(plot_data_comp, aes(x = Model, y = R2, fill = Model)) +
  geom_col(alpha = 0.8, width = 0.6, color = "black", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("R² = %.3f", R2)),
    vjust = -0.5,
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "EMA\n(15 items, ~30 assessments)" = "#2E86AB",
      "Baseline\n(220 items, single)" = "#A23B72"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = percent_format(accuracy = 1)
  ) +
  facet_wrap(~Outcome, nrow = 1) +
  labs(
    title = "Figure 4. Predictive Accuracy: EMA vs Baseline PID-5 Assessment",
    subtitle = "Bayesian R² comparison shows equivalent performance",
    x = NULL,
    y = "Variance Explained (R²)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "figures/figure4_ema_vs_baseline.png",
  p_figure4,
  width = 10,
  height = 6,
  dpi = 300
)

cat("✓ Figure 4 saved\n")

# ==============================================================================
# FIGURE 5: Integrated Summary (The "Everything" Figure)
# ==============================================================================

cat("Creating Figure 5: Integrated Summary...\n")

# Panel A: Stress main effects (simplified)
p5a <- p2a +
  labs(title = "A. Stress Main Effects") +
  theme(legend.position = "none")

# Panel B: Number of credible interactions by domain
if (file.exists("results/moderation_results_summary.csv")) {
  mod_results <- read_csv(
    "results/moderation_results_summary.csv",
    show_col_types = FALSE
  )

  interaction_counts <- mod_results %>%
    filter(grepl(":", parameter), significant == TRUE) %>%
    mutate(
      domain = str_extract(
        parameter,
        "negative_affectivity|detachment|antagonism|disinhibition|psychoticism"
      ) %>%
        str_replace("negative_affectivity", "Neg. Affect.") %>%
        str_replace("detachment", "Detach.") %>%
        str_replace("antagonism", "Antag.") %>%
        str_replace("disinhibition", "Disinhib.") %>%
        str_replace("psychoticism", "Psychot."),
      type = ifelse(grepl("c1_stress", parameter), "Stress", "Recovery")
    ) %>%
    count(domain, type)

  p5b <- ggplot(interaction_counts, aes(x = domain, y = n, fill = type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(
      values = c("Stress" = "#D32F2F", "Recovery" = "#1976D2"),
      name = "Interaction Type"
    ) +
    labs(
      title = "B. Personality Moderation",
      subtitle = "Number of credible interactions",
      x = "PID-5 Domain",
      y = "Count"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
} else {
  p5b <- ggplot() + theme_void()
}

# Panel C: EMA vs Baseline (simplified)
p5c <- plot_data_comp %>%
  ggplot(aes(x = Outcome, y = R2, fill = Model)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(
    values = c(
      "EMA\n(15 items, ~30 assessments)" = "#2E86AB",
      "Baseline\n(220 items, single)" = "#A23B72"
    ),
    name = NULL
  ) +
  labs(
    title = "C. Assessment Validity",
    subtitle = "EMA ≈ Baseline",
    x = NULL,
    y = "R²"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8)
  )

# Panel D: Variability contribution (if available)
if (file.exists("results/variability/version_a_results.csv")) {
  var_results <- read_csv(
    "results/variability/version_a_results.csv",
    show_col_types = FALSE
  )

  var_counts <- tibble(
    Type = c("Mean Effects", "Variability Effects"),
    Count = c(
      sum(
        var_results$mean_stress_credible | var_results$mean_recovery_credible,
        na.rm = TRUE
      ),
      sum(
        var_results$sd_stress_credible | var_results$sd_recovery_credible,
        na.rm = TRUE
      )
    )
  )

  p5d <- ggplot(var_counts, aes(x = Type, y = Count, fill = Type)) +
    geom_col(alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("#2E86AB", "#A23B72")) +
    geom_text(aes(label = Count), vjust = -0.5, fontface = "bold") +
    labs(
      title = "D. Variability Contribution",
      x = NULL,
      y = "Credible Effects"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold")
    )
} else {
  p5d <- ggplot() + theme_void()
}

# Combine all panels
p_figure5 <- (p5a | (p5b / p5c / p5d)) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = "Figure 5. Integrated Summary: Stress, Personality, and Voice",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(
  "figures/figure5_integrated_summary.png",
  p_figure5,
  width = 14,
  height = 10,
  dpi = 300
)

cat("✓ Figure 5 saved\n")

cat("\n", rep("=", 70), "\n")
cat("FIGURE GENERATION COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Figures saved in figures/ directory:\n")
cat("  - figure1_study_design.png\n")
cat("  - figure2_stress_main_effects.png\n")
cat("  - figure3_personality_moderation.png\n")
cat("  - figure4_ema_vs_baseline.png\n")
cat("  - figure5_integrated_summary.png\n\n")
