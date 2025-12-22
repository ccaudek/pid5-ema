# VOICE ACOUSTICS & PERSONALITY PATHOLOGY - ANALYSIS PIPELINE
# Final Version - Female Sample Only (N = 141)
# ==============================================================================

## OVERVIEW

This pipeline analyzes the relationship between personality pathology (PID-5)
and voice acoustic changes during exam stress using Bayesian multilevel models.

**Sample**: 141 female university students
**Design**: Naturalistic stress manipulation (exam period)
**Timepoints**: Baseline → Pre-exam → Post-exam
**Measures**: 
  - Voice: F0, formants, jitter, NNE (vowels /a/, /i/, /u/)
  - Personality: PID-5 full (220 items) + EMA (15 items over 2.5 months)

## EXECUTION ORDER

### CORE ANALYSIS (Required)

**1. Data Cleaning**
```r
source("01_clean_after_merge_FINAL.R")
```
- Filters: careless responding, n_ema range, quality metrics
- **CRITICAL**: Filters to FEMALE ONLY (sex == "Femmina")
- Output: `data/processed/ema_plus_scales_cleaned.csv`

**2. Main Effects & EMA Validation**
```r
source("02_voice_personality_analysis_FINAL.R")
```
- Aggregates EMA PID-5 (mean across 2.5 months)
- **NEW**: Validates EMA vs PID-5 full (convergent correlations)
- Fits main effects models (stress manipulation check)
- Output: `results/df_analysis.rds`, `results/ema_validation_correlations.rds`

**3. Moderation Analysis**
```r
source("03_moderation_analysis_FINAL.R")
```
- Tests PID-5 × Stress interactions (primary hypotheses)
- Contrast coding: c1_stress (reactivity), c2_recovery
- Fits models for all outcomes × vowels
- Output: `results/fitted_models.rds`, `results/results_table.rds`

### SUPPLEMENTARY ANALYSES (Strongly Recommended)

**4. Temporal Covariation**
```r
source("08_temporal_covariation_analysis.R")
```
- **NEW**: Tests if EMA trait fluctuations co-vary with acoustic changes
- Within-person vs between-person decomposition
- Supports abstract claim: "co-varied temporally with acoustic changes"
- Output: `results/temporal_covariation_results.rds`

**5. Model Comparison: EMA vs Baseline**
```r
source("09_model_comparison_EMA_vs_Baseline.R")
```
- **NEW**: Compares EMA (15 items) vs PID-5 full (220 items)
- Uses R² and LOO-IC for model comparison
- Supports abstract claim: "EMA added predictive value beyond baseline"
- Output: `results/ema_vs_baseline_comparison.rds`

**6. Alternative Models for Disinhibition**
```r
source("07_disinhibition_alternative_models.R")
```
- Explores non-moderation patterns for Disinhibition
- Tests: variability, non-linearity, composite outcomes
- **KEY FINDING**: Disinhibition → higher baseline F0 (main effect)
- Output: `results/disinhibition_alternative_models.rds`

## KEY RESULTS TO VERIFY

After running the pipeline, check these results against the abstract:

### Primary Findings (Abstract Claims)

1. **Negative Affectivity × Stress → F0**
   - Expected: β ≈ 5 Hz [1.26, 8.76]
   - Check: `results/results_table.rds` → filter for NA × c1_stress

2. **Detachment × Recovery → F0**
   - Expected: β ≈ -4.37 [-8.30, -0.46]
   - Check: `results/results_table.rds` → filter for Det × c2_recovery

3. **Antagonism × Recovery**
   - Expected: β ≈ 3.44-4.17 across vowels
   - Check: `results/results_table.rds` → filter for Ant × c2_recovery

4. **Psychoticism**: Pervasive effects
   - Check: Multiple significant main effects + interactions

5. **Disinhibition × F2**
   - Expected: β ≈ 63 Hz [26, 99]
   - Check: `results/results_table.rds` → filter for Dis × F2

### New Findings (From Alternative Models)

6. **Disinhibition Main Effect on F0**
   - Finding: β ≈ 3.87 [0.53, 7.27]
   - Location: `results/disinhibition_alternative_models.rds`
   - Interpretation: Higher baseline pitch, not stress moderation

7. **Psychoticism → Reduced Composite Reactivity**
   - Finding: β ≈ -0.228 [-0.34, -0.12]
   - Location: `results/disinhibition_alternative_models.rds`
   - Interpretation: Blunted emotional reactivity

### Validation Analyses

8. **EMA Convergent Validity**
   - Check: `results/ema_validation_correlations.rds`
   - Expected: r > 0.50 for all domains (ideally > 0.70)

9. **Temporal Covariation**
   - Check: `results/temporal_covariation_results.rds`
   - Look for: Significant within-person effects

10. **EMA vs Baseline Comparison**
    - Check: `results/ema_vs_baseline_comparison.rds`
    - Look for: Positive R² difference or negative LOO difference

## MODEL DIAGNOSTICS

All models should show:
- Rhat < 1.01 (convergence)
- ESS > 400 (effective sample size)
- No divergent transitions

Check diagnostics in console output during fitting.

## EXPECTED RUNTIME

- Script 01: ~30 seconds
- Script 02: ~30-45 minutes (6 main effects models)
- Script 03: ~3-5 hours (18 moderation models: 6 outcomes × 3 vowels)
- Script 08: ~2-3 hours (temporal covariation models)
- Script 09: ~1-2 hours (comparison models)
- Script 07: ~2-3 hours (alternative models)

**Total**: ~8-12 hours of computation

## OUTPUT FILES

### Data
- `data/processed/ema_plus_scales_cleaned.csv` - Clean EMA + scales
- `results/df_analysis.rds` - Analysis-ready dataset

### Models
- `models/m1_f0mean_a.rds` - Main effect: F0 mean
- `models/m_f0_mean_a.rds` - Moderation: F0 mean (vowel /a/)
- `models/m_temporal_*.rds` - Temporal covariation models
- `models/m_baseline_*.rds` - Baseline PID-5 comparison models
- ... (many more)

### Results
- `results/ema_validation_correlations.rds` - EMA validity
- `results/results_table.rds` - All moderation effects
- `results/temporal_covariation_results.rds` - Within-person effects
- `results/ema_vs_baseline_comparison.rds` - Model comparison
- `results/disinhibition_alternative_models.rds` - Alternative patterns

### Figures
- `figures/ema_validation.png` - Convergent validity plot
- `figures/temporal_covariation_f0.png` - Within-person effects
- `figures/forest_plot_interactions.png` - All moderation effects

## NOTES FOR MANUSCRIPT

### Sample Description
- **CRITICAL**: Sample is FEMALE ONLY (N = 141)
- Rationale: Control for sex-related variation in vocal pitch
- Limitation: Findings may not generalize to males

### Abstract Updates Required
1. Change "111 university students" → "141 female university students"
2. Verify all effect sizes replicate with N = 141
3. Ensure claims about "temporal covariation" and "predictive value" 
   are supported by new analyses (Scripts 08-09)

### Methods Section
Must include:
- Sample composition (female only) + justification
- EMA validation (convergent correlations with full PID-5)
- Temporal resolution (twice-weekly over 2.5 months)
- Model comparison strategy (EMA vs baseline)

### Results Section
1. Main Effects (manipulation check)
2. Moderation Effects (primary hypotheses)
3. Supplementary Analyses:
   - EMA validation
   - Temporal covariation
   - Alternative patterns (Disinhibition, Psychoticism)
   - Model comparison

## TROUBLESHOOTING

**Problem**: Models don't converge (Rhat > 1.01)
- **Solution**: Increase adapt_delta to 0.999, increase iterations

**Problem**: "Variables not found" errors
- **Solution**: Ensure scripts run in order (01 → 02 → 03)

**Problem**: Missing baseline PID-5 variables
- **Solution**: Check that `ema_plus_scales_merged.RDS` includes 
  `domain_*_baseline` columns

**Problem**: Very long runtime
- **Solution**: Use `cores = parallel::detectCores() - 1` for faster fitting

## CONTACT

For questions about the analysis pipeline, contact Corrado Caudek.

## VERSION HISTORY

- **v1.0** (Dec 2024): Initial pipeline with N = 111
- **v2.0** (Dec 2024): Updated with N = 141, female only, added validation analyses

---
Last updated: December 2024
