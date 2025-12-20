# Pattern Conformity Analysis - Quick Start Guide

## Overview

This analysis tests whether your theoretical model of personality × stress interactions is supported by the overall pattern of results, rather than testing each of 180 effects individually.

## Files

1. **`pattern_conformity_integrated.R`** - Main analysis script (run after moderation_analysis.R)
2. **`master_workflow.R`** - Shows how all scripts fit together

## How to Run

```r
# After running your existing moderation_analysis.R:

# Option 1: If workspace is saved
load("results/moderation_analysis_complete.RData")
source("pattern_conformity_integrated.R")

# Option 2: If fitted_models list exists in environment
source("pattern_conformity_integrated.R")
```

## Customizing Your Predictions (IMPORTANT)

The key step is defining your theoretical predictions in the `predictions` tibble (around line 40 of `pattern_conformity_integrated.R`).

### Structure of Predictions Table

| Column | Description |
|--------|-------------|
| `outcome` | Acoustic variable: f0_mean, f0_std, jitter, nne, f2_mean, f2_std |
| `vowel` | Which vowel: a, i, u |
| `contrast` | Which stress contrast: c1_stress (reactivity) or c2_recovery |
| `pid5_domain` | PID-5 domain: negative_affectivity, detachment, antagonism, disinhibition, psychoticism |
| `expected_direction` | +1 (positive), -1 (negative), or 0 (no prediction/exclude) |
| `rationale` | Brief explanation (for documentation) |

### Example Predictions

```r
predictions <- tribble(
  ~outcome,   ~vowel, ~contrast,     ~pid5_domain,           ~expected_direction, ~rationale,
  
  # Negative Affectivity amplifies stress reactivity
  "f0_mean",  "a",    "c1_stress",   "negative_affectivity", +1, "NA increases pitch under stress",
  
  # Detachment impairs recovery
  "f0_mean",  "a",    "c2_recovery", "detachment",           -1, "Detachment slows return to baseline",
  
  # No prediction for this combination (will be excluded)
  "jitter",   "a",    "c1_stress",   "antagonism",            0, "No theoretical prediction"
)
```

### Prediction Guidelines

1. **Only include predictions you have theoretical justification for**
   - Effects with `expected_direction = 0` are excluded from conformity calculation
   - This is a feature, not a bug—it keeps you honest about what you're testing

2. **Direction coding:**
   - `+1` = You expect the interaction to be POSITIVE
   - `-1` = You expect the interaction to be NEGATIVE
   - `0` = No prediction (excluded from analysis)

3. **For recovery contrast (c2_recovery):**
   - Positive coefficient = faster/better recovery
   - Negative coefficient = impaired/slower recovery
   - Be careful about the direction interpretation!

4. **Replication across vowels:**
   - I've included predictions for all three vowels where you expect the effect to replicate
   - You could also restrict to /a/ only for a more focused test

### Current Predictions (Based on Your Abstract)

The script currently includes predictions based on your abstract findings:

| Domain | × Stress Reactivity | × Recovery | Rationale |
|--------|---------------------|------------|-----------|
| Negative Affectivity | + (amplifies) | − (impairs) | Heightened emotional reactivity |
| Detachment | − (blunts) | − (impairs) | Emotional withdrawal |
| Antagonism | 0 (no prediction) | + (facilitates) | Callousness aids rapid recovery |
| Disinhibition | + on F2 only | 0 | Articulatory/motor effects |
| Psychoticism | + | − (persistent) | Stable arousal/disorganization |

**Please review these and adjust based on your actual theoretical predictions!**

## Output

The script produces:

### Console Output
- Summary metrics (mean conformity, P(>chance), combined PD)
- Individual prediction results
- Summary by domain and outcome
- Manuscript-ready text

### Files
- `results/pattern_conformity_results.rds` - Full results object
- `results/pattern_conformity_by_prediction.csv` - Per-prediction PDs
- `results/pattern_conformity_by_domain.csv` - Summary by PID-5 domain
- `figures/pattern_conformity_main.png` - Main figure (posterior of conformity)
- `figures/pattern_conformity_individual_pds.png` - Individual PDs bar chart

## Interpreting Results

### Key Metrics

| Metric | Interpretation |
|--------|----------------|
| **Mean proportion conforming** | What % of predictions are in the expected direction? (chance = 50%) |
| **P(>chance)** | Posterior probability that conformity exceeds 50%. Values >95% = strong evidence |
| **Combined PD** | Aggregate probability of direction. >80% = strong directional consistency |
| **P(all conform)** | Very stringent: probability ALL predictions are simultaneously correct |

### Evidence Thresholds

- **P(>chance) > 95%**: Strong evidence the pattern holds
- **P(>chance) 80-95%**: Moderate evidence
- **P(>chance) < 80%**: Weak/inconclusive

## Troubleshooting

### "Parameter not found" warnings

The script tries to match parameter names from your brms models. If you see warnings:

1. Check the parameter names printed in the console output
2. The expected format is: `b_c1_stress:pid5_negative_affectivity_bl_c`
3. Adjust the `param_name` construction in the script if needed

### "Model not found" warnings

Make sure your model names follow the pattern: `m_[outcome]_[vowel]`
(e.g., `m_f0_mean_a`, `m_jitter_i`, `m_nne_u`)

## Questions?

The key things to customize:
1. Your theoretical predictions in the `predictions` tribble
2. Whether to include all vowels or focus on /a/
3. Which outcomes have strong enough theoretical predictions to include

Everything else should work automatically with your existing fitted models.
