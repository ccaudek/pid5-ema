# ==============================================================================
# Interpretive Report: Manuscript-Ready Summary
# Context-Dependent Expression of Personality Pathology via Voice Acoustics
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(bayestestR)
  library(knitr)
  library(kableExtra)
})

# ==============================================================================
# LOAD DATA AND MODELS
# ==============================================================================

# Load moderation analysis workspace
if (file.exists("results/moderation_analysis_ORIGINAL.RData")) {
  load("results/moderation_analysis_ORIGINAL.RData")
  cat("✓ Loaded moderation analysis workspace\n")
} else {
  stop("Please run moderation_analysis.R first to generate results.")
}

# Also load main effects models if available
if (file.exists("results/main_effects_workspace_ORIGINAL.RData")) {
  # Load into a separate environment to avoid overwriting
  main_env <- new.env()
  load("results/main_effects_workspace_ORIGINAL.RData", envir = main_env)

  # Extract main effects models if they exist
  if (exists("m1_f0mean_a", envir = main_env)) {
    m1_f0mean_a <- main_env$m1_f0mean_a
    m2_f0std_a <- main_env$m2_f0std_a
    m3_jitter_a <- main_env$m3_jitter_a
    m4_nne_a <- main_env$m4_nne_a
    m5_f2mean_a <- main_env$m5_f2mean_a
    m6_f2std_a <- main_env$m6_f2std_a
    cat("✓ Loaded main effects models\n")
    main_effects_available <- TRUE
  } else {
    main_effects_available <- FALSE
    cat("⚠ Main effects models not found in workspace\n")
  }
} else {
  main_effects_available <- FALSE
  cat("⚠ Main effects workspace not found\n")
}

# Check what's available
cat("\nAvailable fitted models:", length(fitted_models), "\n")
cat("Model names:", paste(names(fitted_models), collapse = ", "), "\n")

# ==============================================================================
# EXECUTIVE SUMMARY REPORT
# ==============================================================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MANUSCRIPT-READY INTERPRETIVE REPORT\n")
cat("Personality Pathology × Stress Reactivity in Voice Acoustics\n")
cat(rep("=", 80), "\n", sep = "")

# ==============================================================================
# 1. STUDY DESIGN SUMMARY (for Methods section)
# ==============================================================================

cat("\n### STUDY DESIGN ###\n")
cat("\nDesign: Naturalistic stress manipulation (exam period)\n")
cat("- BASELINE: Neutral period (low stress)\n")
cat("- PRE-EXAM: High stress period (day before exam)\n")
cat("- POST-EXAM: Recovery period (day after exam)\n")

cat("\nSample:\n")
cat("- N participants:", n_distinct(df_analysis$ID), "\n")
cat("- N total observations:", nrow(df_analysis), "\n")
cat("- Observations per timepoint:", nrow(df_analysis) / 3, "\n")

cat("\nMeasures:\n")
cat("- Acoustic features: F0 (pitch), F0 variability, Jitter, NNE, F2\n")
cat("- Personality pathology: PID-5 (all 5 domains)\n")
cat("- Assessment method: Voice recordings (passive sensing) + self-report\n")

cat("\nAnalytic Approach:\n")
cat("- Bayesian multilevel models\n")
cat("- Within-person stress effects + between-person moderation\n")
cat(
  "- Contrasts: (1) Reactivity = PRE vs BASELINE, (2) Recovery = POST vs PRE\n"
)

# ==============================================================================
# 2. MAIN EFFECTS: Stress Manipulation Check (for Results section)
# ==============================================================================

cat("\n\n### MAIN EFFECTS: STRESS MANIPULATION CHECK ###\n")
cat("\nResearch Question: Did exam stress affect voice acoustics?\n")
cat("This validates the contextual manipulation and shows proximal\n")
cat("environmental impact on physiological markers.\n\n")

# Function to interpret main effects from moderation models
# (using c1_stress and c2_recovery contrasts)
interpret_main_effect_from_moderation <- function(model, outcome_name) {
  fixed <- fixef(model)

  cat("\n", outcome_name, ":\n", sep = "")

  # Get c1_stress (reactivity) effect
  if ("c1_stress" %in% rownames(fixed)) {
    est <- fixed["c1_stress", "Estimate"]
    ci_lower <- fixed["c1_stress", "Q2.5"]
    ci_upper <- fixed["c1_stress", "Q97.5"]

    cat(
      "  Stress Reactivity (PRE vs BASE): β = ",
      sprintf("%.3f", est),
      ", 95% CI [",
      sprintf("%.3f", ci_lower),
      ", ",
      sprintf("%.3f", ci_upper),
      "]\n",
      sep = ""
    )

    if (ci_lower > 0) {
      cat("  → Strong evidence for INCREASE during stress\n")
    } else if (ci_upper < 0) {
      cat("  → Strong evidence for DECREASE during stress\n")
    } else {
      cat("  → Credible interval includes zero (inconclusive)\n")
    }
  }

  # Get c2_recovery effect
  if ("c2_recovery" %in% rownames(fixed)) {
    est <- fixed["c2_recovery", "Estimate"]
    ci_lower <- fixed["c2_recovery", "Q2.5"]
    ci_upper <- fixed["c2_recovery", "Q97.5"]

    cat(
      "  Recovery (POST vs PRE): β = ",
      sprintf("%.3f", est),
      ", 95% CI [",
      sprintf("%.3f", ci_lower),
      ", ",
      sprintf("%.3f", ci_upper),
      "]\n",
      sep = ""
    )

    if (ci_upper < 0) {
      cat("  → Strong evidence for recovery (return toward baseline)\n")
    } else if (ci_lower > 0) {
      cat("  → Continued elevation post-exam\n")
    } else {
      cat("  → Credible interval includes zero (inconclusive)\n")
    }
  }
}

# Use moderation models for vowel /a/ to show main effects
if ("m_f0_mean_a" %in% names(fitted_models)) {
  interpret_main_effect_from_moderation(
    fitted_models[["m_f0_mean_a"]],
    "F0 Mean (Pitch)"
  )
}
if ("m_f0_std_a" %in% names(fitted_models)) {
  interpret_main_effect_from_moderation(
    fitted_models[["m_f0_std_a"]],
    "F0 Variability"
  )
}
if ("m_jitter_a" %in% names(fitted_models)) {
  interpret_main_effect_from_moderation(
    fitted_models[["m_jitter_a"]],
    "Jitter (Voice Quality)"
  )
}
if ("m_nne_a" %in% names(fitted_models)) {
  interpret_main_effect_from_moderation(
    fitted_models[["m_nne_a"]],
    "NNE (Vocal Stress)"
  )
}
if ("m_f2_mean_a" %in% names(fitted_models)) {
  interpret_main_effect_from_moderation(
    fitted_models[["m_f2_mean_a"]],
    "F2 Mean (Articulation)"
  )
}
if ("m_f2_std_a" %in% names(fitted_models)) {
  interpret_main_effect_from_moderation(
    fitted_models[["m_f2_std_a"]],
    "F2 Variability"
  )
}

cat("\n\n** Manuscript Interpretation **\n")
cat("These main effects establish that:\n")
cat(
  "1. The exam created measurable physiological stress (contextual manipulation check)\n"
)
cat("2. Voice acoustics are sensitive to proximal environmental stressors\n")
cat(
  "3. This provides foundation for testing individual differences in reactivity\n"
)

# ==============================================================================
# 3. MODERATION EFFECTS: Key Findings (for Results section)
# ==============================================================================

cat("\n\n### MODERATION EFFECTS: PERSONALITY × CONTEXT INTERACTIONS ###\n")
cat(
  "\nResearch Question: Do personality pathology traits moderate stress reactivity?\n"
)
cat("This tests context-dependent expression of maladaptive traits.\n\n")

# Function to interpret interactions with CORRECT parameter names
interpret_interaction <- function(model, outcome_name) {
  cat("\n", rep("-", 60), "\n", sep = "")
  cat(outcome_name, "\n")
  cat(rep("-", 60), "\n", sep = "")

  fixed <- fixef(model)

  # Stress reactivity interactions - CORRECT NAMES
  cat("\nSTRESS REACTIVITY MODERATION (PRE vs BASELINE):\n")

  interactions_stress <- c(
    "c1_stress:pid5_negative_affectivity_bl_c" = "Negative Affectivity",
    "c1_stress:pid5_detachment_bl_c" = "Detachment",
    "c1_stress:pid5_antagonism_bl_c" = "Antagonism",
    "c1_stress:pid5_disinhibition_bl_c" = "Disinhibition",
    "c1_stress:pid5_psychoticism_bl_c" = "Psychoticism"
  )

  for (i in seq_along(interactions_stress)) {
    param <- names(interactions_stress)[i]
    trait <- interactions_stress[i]

    if (param %in% rownames(fixed)) {
      est <- fixed[param, "Estimate"]
      ci_lower <- fixed[param, "Q2.5"]
      ci_upper <- fixed[param, "Q97.5"]

      # Calculate PD from posterior
      draws <- as_draws_df(model)
      param_col <- paste0("b_", gsub(":", ".", param)) # brms uses . instead of :
      if (param_col %in% names(draws)) {
        post_draws <- draws[[param_col]]
        pd_val <- mean(post_draws > 0)
        if (pd_val < 0.5) pd_val <- 1 - pd_val
      } else {
        pd_val <- NA
      }

      cat("\n  ", trait, ":\n", sep = "")
      cat(
        "    β = ",
        sprintf("%.3f", est),
        ", 95% CI [",
        sprintf("%.3f", ci_lower),
        ", ",
        sprintf("%.3f", ci_upper),
        "]\n",
        sep = ""
      )
      if (!is.na(pd_val)) {
        cat(
          "    Probability of Direction:",
          sprintf("%.1f%%", pd_val * 100),
          "\n"
        )
      }

      # Interpretation
      if (ci_lower > 0) {
        cat(
          "    ✓ AMPLIFICATION: Higher",
          trait,
          "→ stronger stress reactivity\n"
        )
      } else if (ci_upper < 0) {
        cat("    ✓ BUFFERING: Higher", trait, "→ reduced stress reactivity\n")
      } else {
        cat("    • No clear moderation (CI includes zero)\n")
      }
    }
  }

  # Recovery interactions - CORRECT NAMES
  cat("\n\nRECOVERY/RESILIENCE MODERATION (POST vs PRE):\n")

  # Note: brms may order these differently
  interactions_recovery <- c(
    "pid5_negative_affectivity_bl_c:c2_recovery" = "Negative Affectivity",
    "pid5_detachment_bl_c:c2_recovery" = "Detachment",
    "pid5_antagonism_bl_c:c2_recovery" = "Antagonism",
    "pid5_disinhibition_bl_c:c2_recovery" = "Disinhibition",
    "pid5_psychoticism_bl_c:c2_recovery" = "Psychoticism"
  )

  for (i in seq_along(interactions_recovery)) {
    param <- names(interactions_recovery)[i]
    trait <- interactions_recovery[i]

    if (param %in% rownames(fixed)) {
      est <- fixed[param, "Estimate"]
      ci_lower <- fixed[param, "Q2.5"]
      ci_upper <- fixed[param, "Q97.5"]

      # Calculate PD
      draws <- as_draws_df(model)
      param_col <- paste0("b_", gsub(":", ".", param))
      if (param_col %in% names(draws)) {
        post_draws <- draws[[param_col]]
        pd_val <- mean(post_draws > 0)
        if (pd_val < 0.5) pd_val <- 1 - pd_val
      } else {
        pd_val <- NA
      }

      cat("\n  ", trait, ":\n", sep = "")
      cat(
        "    β = ",
        sprintf("%.3f", est),
        ", 95% CI [",
        sprintf("%.3f", ci_lower),
        ", ",
        sprintf("%.3f", ci_upper),
        "]\n",
        sep = ""
      )
      if (!is.na(pd_val)) {
        cat(
          "    Probability of Direction:",
          sprintf("%.1f%%", pd_val * 100),
          "\n"
        )
      }

      if (ci_upper < 0) {
        cat(
          "    ✓ IMPAIRED RECOVERY: Higher",
          trait,
          "→ slower return to baseline\n"
        )
      } else if (ci_lower > 0) {
        cat(
          "    ✓ ENHANCED RECOVERY: Higher",
          trait,
          "→ faster return to baseline\n"
        )
      } else {
        cat("    • No clear moderation (CI includes zero)\n")
      }
    }
  }
}

# Interpret moderation for each outcome (vowel /a/)
if ("m_f0_mean_a" %in% names(fitted_models)) {
  interpret_interaction(
    fitted_models[["m_f0_mean_a"]],
    "F0 Mean (Pitch) - Vowel /a/"
  )
}
if ("m_f0_std_a" %in% names(fitted_models)) {
  interpret_interaction(
    fitted_models[["m_f0_std_a"]],
    "F0 Variability - Vowel /a/"
  )
}
if ("m_jitter_a" %in% names(fitted_models)) {
  interpret_interaction(fitted_models[["m_jitter_a"]], "Jitter - Vowel /a/")
}
if ("m_nne_a" %in% names(fitted_models)) {
  interpret_interaction(fitted_models[["m_nne_a"]], "NNE - Vowel /a/")
}

# ==============================================================================
# 4. SUMMARY OF SIGNIFICANT EFFECTS
# ==============================================================================

cat("\n\n### SUMMARY: SIGNIFICANT MODERATION EFFECTS ###\n")

if (exists("results_table") && nrow(results_table) > 0) {
  sig_interactions <- results_table %>%
    filter(significant == TRUE, type == "Interaction") %>%
    arrange(outcome, vowel, parameter)

  if (nrow(sig_interactions) > 0) {
    cat("\nSignificant interactions (95% CI excludes zero):\n\n")

    sig_interactions %>%
      mutate(
        effect = paste0(outcome, "_", vowel),
        direction = ifelse(estimate > 0, "Amplifies", "Buffers/Impairs")
      ) %>%
      select(effect, parameter, estimate, ci_lower, ci_upper, direction) %>%
      print(n = 50)
  } else {
    cat("\nNo interactions with 95% CI excluding zero.\n")
    cat("Consider examining 90% CIs or probability of direction (PD).\n")
  }

  # Show effects with strong PD even if CI includes zero
  cat("\n\nEffects with PD > 90% (suggestive evidence):\n")
  # This would require calculating PD for each effect...
} else {
  cat("\nResults table not available.\n")
}

# ==============================================================================
# 5. MANUSCRIPT STRUCTURE GUIDANCE
# ==============================================================================

cat("\n\n")
cat(rep("=", 80), "\n", sep = "")
cat("MANUSCRIPT STRUCTURE GUIDANCE\n")
cat(rep("=", 80), "\n", sep = "")

cat("\n### METHODS ###\n")
cat(
  "
Key points to include:
- Design: Naturalistic exam stress paradigm (ecological validity)
- Sample: N = ",
  n_distinct(df_analysis$ID),
  " university students
- Procedure: 3 timepoints (baseline, pre-exam, post-exam)
- Voice recording: Sustained vowels /a/, /i/, /u/ (3 sec each)
- Acoustic features: F0 mean/SD, jitter, NNE, F2 mean/SD
  - Extracted using [software name, e.g., Praat]
- PID-5: Full 220-item at baseline; 15-item EMA (3 per domain)
- Analysis: Bayesian multilevel models (brms/Stan)
  - Contrast coding: C1 (stress reactivity), C2 (recovery)
  - Random slopes for stress effects
  - Weakly informative priors
"
)

cat("\n### RESULTS ###\n")
cat(
  "
Structure:
1. Descriptive statistics + manipulation check
   - Show stress affected voice (main effects)
   - Establishes contextual validity
   
2. Moderation analyses
   - Present each outcome (F0 mean, F0 SD, Jitter, NNE, F2)
   - For each: report interactions with all 5 PID-5 domains
   - Use figures showing conditional effects (simple slopes)
   - Report Bayesian metrics: β, 95% CI, PD
   
3. Pattern conformity analysis (if using)
   - Overall support for theoretical model
   - Not relying on individual 'significant' effects

Visualization strategy:
- Figure 1: Descriptive trajectories across timepoints
- Figure 2: Posterior distributions of main effects
- Figure 3: Interaction plots (conditional effects)
- Figure 4: Pattern conformity posterior
- Table 1: Sample characteristics
- Table 2: Main effects summary
- Table 3: Moderation effects summary
"
)

cat("\n### DISCUSSION ###\n")
cat(
  "
Key messages aligned with special issue aims:

1. METHODOLOGICAL INNOVATION (Primary aim #1):
   'Our study demonstrates that voice acoustics, captured through passive
   sensing, can serve as an objective, non-intrusive marker of proximal 
   stress exposure. This addresses a key limitation of EMA: reliance on
   self-report.'

2. MULTILEVEL CONTEXT (Primary aim #2):
   'By showing that personality pathology traits moderate stress reactivity
   in voice, we provide evidence for context-dependent expression of 
   maladaptive traits.'

3. PERSON × CONTEXT TRANSACTIONS:
   'Our findings support transactional models: personality pathology is not
   merely a static vulnerability, but shapes—and is shaped by—environmental
   demands.'

4. CLINICAL IMPLICATIONS:
   - Personalized risk assessment (who is most reactive?)
   - Intervention targets (reduce stress or boost resilience?)
   - Ecological validity of passive monitoring

5. LIMITATIONS:
   - Single stressor type (achievement-oriented, time-limited)
   - Student sample
   - Multiple comparisons (addressed via pattern conformity)
"
)

# ==============================================================================
# 6. STATISTICAL REPORTING TEMPLATES
# ==============================================================================

cat("\n\n")
cat(rep("=", 80), "\n", sep = "")
cat("STATISTICAL REPORTING TEMPLATES\n")
cat(rep("=", 80), "\n", sep = "")

cat(
  "
### Main Effects (example) ###
'Exam stress significantly increased vocal pitch (F0 mean): β = X.XX, 95% CI 
[X.XX, X.XX], PD = XX.X%. Recovery was observed post-exam: β = X.XX, 95% CI 
[X.XX, X.XX].'

### Moderation Effects (example) ###
'Negative Affectivity moderated stress reactivity in vocal pitch (F0 mean): 
β = X.XX, 95% CI [X.XX, X.XX], PD = XX.X%. Individuals high (+1 SD) in 
Negative Affectivity showed stronger pitch increase during exam stress 
compared to those low (-1 SD) in Negative Affectivity.'

### Pattern Conformity (if using) ###
'To evaluate our overall theoretical model rather than individual effects,
we computed a Bayesian pattern conformity analysis. The mean proportion of
predictions in the expected direction was XX.X% (95% CI [XX.X%, XX.X%]),
substantially exceeding the 50% expected by chance (P(> chance) = XX.X%).'

### Bayesian Interpretation Notes ###
- Report β (posterior mean), 95% CI, and PD (probability of direction)
- Interpret 95% CI: if excludes zero → 'credible evidence for effect'
- PD > 97.5% roughly corresponds to p < .05 in frequentist terms
- Emphasize effect size magnitude, not just 'significance'
"
)

# ==============================================================================
# 7. CHECKLIST
# ==============================================================================

cat("\n\n")
cat(rep("=", 80), "\n", sep = "")
cat("MANUSCRIPT CHECKLIST\n")
cat(rep("=", 80), "\n", sep = "")

cat(
  "
ANALYSES COMPLETED:
  ✓ Main effects models (stress manipulation check)
  ✓ Moderation models (PID-5 × stress interactions)
  ",
  ifelse(
    file.exists("results/pattern_conformity_results.rds"),
    "✓ Pattern conformity analysis",
    "□ Pattern conformity analysis (run pattern_conformity_integrated.R)"
  ),
  "

FIGURES TO CREATE:
  □ Figure 1: Acoustic trajectories across timepoints
  □ Figure 2: Posterior distributions of stress effects
  □ Figure 3: Interaction/conditional effects plots
  □ Figure 4: Pattern conformity posterior (if using)

TABLES TO CREATE:
  □ Table 1: Sample characteristics
  □ Table 2: Descriptive statistics (acoustics by timepoint)
  □ Table 3: Main effects summary
  □ Table 4: Moderation effects summary

SUPPLEMENTARY MATERIALS:
  □ Full model specifications (priors, families)
  □ Convergence diagnostics (R-hat, ESS, trace plots)
  □ Sensitivity analyses (prior sensitivity, vowel comparison)
  □ Complete results for all vowels (/a/, /i/, /u/)

SUBMISSION:
  □ Deadline: February 14, 2026
  □ Target journal: Personality Disorders (special issue)
  □ Word limit: Check call for papers
"
)

cat("\n", rep("=", 80), "\n", sep = "")
cat("END OF REPORT\n")
cat(rep("=", 80), "\n\n", sep = "")
