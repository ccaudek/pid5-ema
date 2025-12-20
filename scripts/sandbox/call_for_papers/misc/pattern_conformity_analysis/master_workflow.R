# ==============================================================================
# MASTER WORKFLOW: Voice Acoustics × Personality Pathology Analysis
# ==============================================================================
#
# This script shows the complete analysis pipeline and execution order.
# Run each section in order, or source individual scripts.
#
# ==============================================================================

# ==============================================================================
# STEP 0: SETUP
# ==============================================================================

# Required packages
required_packages <- c(
  "tidyverse",
  "brms",
  "cmdstanr",
  "bayestestR",
  "bayesplot",
  "tidybayes",
  "patchwork",
  "ggdist",
  "marginaleffects",
  "posterior",
  "here",
  "readxl",
  "missRanger",
  "knitr",
  "kableExtra"
)

# Install missing packages
missing <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]
if (length(missing) > 0) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set options
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)

cat("Setup complete.\n")

# ==============================================================================
# STEP 1: DATA PREPARATION & MAIN EFFECTS
# ==============================================================================
#
# Run: voice_personality_analysis.R
#
# This script:
# - Loads and cleans acoustic data from Excel files
# - Merges with PID-5 baseline and EMA data
# - Handles missing data with missRanger
# - Fits main effects models (timepoint effect on acoustics)
# - Creates descriptive visualizations
#
# Output:
# - df_analysis: cleaned dataset for modeling
# - Main effects models (m1_f0mean_a, m2_f0std_a, etc.)
# - results/main_effects_workspace.RData
# - figures/acoustic_trajectories.png
# - figures/posterior_stress_effects.png

cat("\n=== STEP 1: Data Preparation & Main Effects ===\n")
source(
  here::here("scripts", "call_for_papers", "voice_personality_analysis.R")
)

# ==============================================================================
# STEP 2: MODERATION ANALYSIS (All PID-5 Domains × All Outcomes × All Vowels)
# ==============================================================================
#
# Run: moderation_analysis.R
#
# This script:
# - Creates contrast codes (C1: stress reactivity, C2: recovery)
# - Sets up weakly informative priors calibrated to acoustic scales
# - Fits 18 moderation models (6 outcomes × 3 vowels)
# - Each model includes all 5 PID-5 domains as moderators
# - Extracts and summarizes all fixed effects
#
# Output:
# - fitted_models: list of 18 brmsfit objects
# - results_table: summary of all effects
# - models/m_[outcome]_[vowel].rds: individual model files
# - results/moderation_analysis_complete.RData

cat("\n=== STEP 2: Moderation Analysis ===\n")
source(here::here("scripts", "call_for_papers", "moderation_analysis_fixed.R"))

# ==============================================================================
# STEP 3: PATTERN CONFORMITY ANALYSIS (NEW)
# ==============================================================================
#
# Run: pattern_conformity_integrated.R
#
# This script:
# - Defines theoretical predictions (expected direction for each effect)
# - Extracts posterior draws from all fitted models
# - Computes pattern conformity metrics:
#   * Proportion of predictions in expected direction (per draw)
#   * P(pattern > chance)
#   * Combined probability of direction
#   * Individual PDs for each prediction
# - Creates visualizations of overall pattern support
# - Generates manuscript-ready text
#
# Output:
# - results/pattern_conformity_results.rds
# - results/pattern_conformity_by_prediction.csv
# - results/pattern_conformity_by_domain.csv
# - figures/pattern_conformity_main.png
# - figures/pattern_conformity_individual_pds.png

cat("\n=== STEP 3: Pattern Conformity Analysis ===\n")
source(here::here(
  "scripts",
  "call_for_papers",
  "pattern_conformity_analysis",
  "pattern_conformity_integrated.R"
))

# ==============================================================================
# STEP 4: SUPPLEMENTARY ANALYSES
# ==============================================================================
#
# Run: utilities.R
#
# This script:
# - Creates publication-ready tables (sample characteristics, effects)
# - Generates diagnostic plots (trace, R-hat, ESS, pp_check)
# - Conducts sensitivity analyses:
#   * Prior sensitivity (skeptical vs. default priors)
#   * Outlier influence
# - Calculates standardized effect sizes
# - Creates supplementary figures
#
# Output:
# - results/table1_sample_characteristics.csv
# - results/table2_main_effects.csv
# - results/table3_moderation_effects.csv
# - results/sensitivity_*.csv
# - results/effect_sizes.csv
# - figures/diagnostics_*.png
# - figures/supplementary_*.png

cat("\n=== STEP 4: Supplementary Analyses ===\n")
# source(here::here("scripts","call_for_papers","utilities.R"))

# ==============================================================================
# STEP 5: COMPOSITE PID-5 ANALYSIS (Optional/Supplementary)
# ==============================================================================
#
# Run: moderation_composite_pid5_score.R
#
# This script:
# - Creates composite PID-5 score (mean of standardized domains)
# - Fits models with composite as single moderator
# - Useful for showing general personality pathology effect
# - Optional: Latent Profile Analysis of stress reactivity patterns
#
# Output:
# - models with composite PID-5 moderator
# - brms_fixed_effects_summary.csv

cat("\n=== STEP 5: Composite PID-5 Analysis (Optional) ===\n")
source(here::here(
  "scripts",
  "call_for_papers",
  "moderation_composite_pid5_score.R"
))

# ==============================================================================
# STEP 6: GENERATE MANUSCRIPT REPORT
# ==============================================================================
#
# Run: manuscript_report.R
#
# This script:
# - Loads all results
# - Generates interpretive summaries for each section
# - Provides templates for Results and Discussion
# - Creates checklist of figures and tables
#
# Output:
# - Console output with manuscript guidance

cat("\n=== STEP 6: Manuscript Report ===\n")
source(here::here("scripts", "call_for_papers", "manuscript_report.R"))

# ==============================================================================
# QUICK EXECUTION (after models are fitted)
# ==============================================================================

# If models are already fitted and saved, you can run just the new analysis:

run_pattern_conformity_only <- function() {
  cat("\n=== Running Pattern Conformity Analysis ===\n")

  # Load fitted models
  if (file.exists("results/moderation_analysis_complete.RData")) {
    load("results/moderation_analysis_complete.RData")
  } else {
    stop("Please run moderation_analysis.R first to fit the models.")
  }

  # Run pattern conformity
  source("pattern_conformity_integrated.R")

  cat("\n✓ Pattern conformity analysis complete.\n")
  cat("Check results/ and figures/ for outputs.\n")
}

# Uncomment to run:
# run_pattern_conformity_only()

# ==============================================================================
# FILE STRUCTURE
# ==============================================================================
#
# Your project should have this structure:
#
# project/
# ├── data/
# │   ├── raw/
# │   │   └── acustic_features/datiacustici/AUDIO.xlsx
# │   └── processed/
# │       └── ema_plus_scales_cleaned.csv
# ├── models/
# │   └── m_[outcome]_[vowel].rds (saved brms models)
# ├── results/
# │   ├── moderation_analysis_complete.RData
# │   ├── pattern_conformity_results.rds
# │   ├── pattern_conformity_by_prediction.csv
# │   ├── pattern_conformity_by_domain.csv
# │   └── [other tables]
# ├── figures/
# │   ├── pattern_conformity_main.png
# │   ├── pattern_conformity_individual_pds.png
# │   └── [other figures]
# └── scripts/
#     ├── voice_personality_analysis.R
#     ├── moderation_analysis.R
#     ├── pattern_conformity_integrated.R   # NEW
#     ├── utilities.R
#     ├── moderation_composite_pid5_score.R
#     ├── manuscript_report.R
#     └── master_workflow.R                 # THIS FILE
#
# ==============================================================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("WORKFLOW GUIDE LOADED\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nTo run the complete pipeline:\n")
cat("1. source('voice_personality_analysis.R')\n")
cat("2. source('moderation_analysis.R')\n")
cat("3. source('pattern_conformity_integrated.R')  # NEW\n")
cat("4. source('utilities.R')\n")
cat("5. source('manuscript_report.R')\n")
cat("\nOr, if models already fitted:\n")
cat("run_pattern_conformity_only()\n\n")
