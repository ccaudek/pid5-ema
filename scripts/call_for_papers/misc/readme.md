# Bayesian Analysis: Voice Acoustics, Stress, and Personality Pathology

## Overview

This analysis package implements a Bayesian multilevel modeling approach to examine how personality pathology traits (PID-5) moderate voice acoustic responses to exam stress. The study addresses the special issue call on "Context in Personality Pathology Assessment" by demonstrating:

1. **Innovative passive sensing**: Voice acoustics as objective stress markers
2. **Context-dependent expression**: Personality × situation interactions
3. **Multilevel modeling**: Within-person stress effects moderated by between-person traits

---

## File Structure

```
project/
├── AUDIO.xlsx                           # Input data (3 sheets: BASELINE, PRE, POST)
├── voice_personality_analysis.R         # STEP 1: Main effects & descriptives
├── moderation_analysis.R                # STEP 2: Moderation by PID-5 traits
├── manuscript_report.R                  # STEP 3: Interpretive report
├── utilities.R                          # STEP 4: Tables, diagnostics, sensitivity
├── master_script.R                      # RUN ALL: Complete pipeline
├── README.md                            # This file
│
├── models/                              # Fitted brms models (auto-saved)
│   ├── m1_f0mean_main.rds
│   ├── m1_f0mean_moderation.rds
│   └── ...
│
├── figures/                             # Generated plots
│   ├── acoustic_trajectories.png
│   ├── posterior_stress_effects.png
│   ├── moderation_interactions.png
│   ├── interaction_posteriors.png
│   ├── diagnostics_*.png
│   └── supplementary_*.png
│
└── results/                             # Output tables & workspaces
    ├── main_effects_workspace.RData
    ├── moderation_analysis_complete.RData
    ├── table1_sample_characteristics.csv
    ├── table2_main_effects.csv
    ├── table3_moderation_effects.csv
    ├── diagnostic_summary.csv
    └── effect_sizes.csv
```

---

## Quick Start

### Prerequisites

```r
# Required packages
install.packages(c(
  "tidyverse", "readxl", "brms", "cmdstanr",
  "bayestestR", "bayesplot", "tidybayes",
  "patchwork", "ggdist", "knitr", "kableExtra",
  "marginaleffects"
))

# Install cmdstan (if not already installed)
cmdstanr::install_cmdstan()
```

### Running the Full Pipeline

**Option 1: Run everything at once** (recommended for first run)
```r
source("master_script.R")
```
⏱️ Expected runtime: 2-3 hours (parallel chains)

**Option 2: Run step-by-step** (recommended for exploration)
```r
source("voice_personality_analysis.R")      # ~30 min
source("moderation_analysis.R")             # ~90 min
source("manuscript_report.R")               # ~5 min
source("utilities.R")                       # ~20 min
```

---

## Analysis Pipeline Details

### STEP 1: Main Effects & Descriptives (`voice_personality_analysis.R`)

**Purpose**: Establish that exam stress affects voice acoustics (manipulation check)

**What it does**:
- Loads and cleans data from three timepoints
- Computes descriptive statistics and effect sizes
- Fits 4 Bayesian models testing stress effects:
  - F0 mean ~ timepoint + (1 + timepoint | ID)
  - F0 SD ~ timepoint + (1 + timepoint | ID)
  - Jitter ~ timepoint + (1 + timepoint | ID)
  - NNE ~ timepoint + (1 + timepoint | ID)
- Creates visualization of stress manipulation

**Key outputs**:
- `figures/acoustic_trajectories.png` - Shows mean changes across timepoints
- `figures/posterior_stress_effects.png` - Bayesian estimates with credible intervals
- `results/main_effects_workspace.RData` - Workspace for next steps

**Interpretation focus**:
- Do voice features change from BASELINE → PRE (stress onset)?
- Do they recover from PRE → POST (stress offset)?
- Effect sizes: Cohen's d for within-person changes

---

### STEP 2: Moderation Analysis (`moderation_analysis.R`)

**Purpose**: Test if PID-5 traits moderate stress reactivity (key research question)

**What it does**:
- Creates contrast codes:
  - C1: PRE vs BASELINE (reactivity to stress)
  - C2: POST vs PRE (recovery from stress)
- Fits 4 moderation models with cross-level interactions:
  - Outcome ~ C1 × PID5_NA + C1 × PID5_DIS + C1 × PID5_PSY +
              C2 × PID5_NA + C2 × PID5_DIS + C2 × PID5_PSY +
              (1 + C1 + C2 | ID)
- Tests 6 interactions per outcome (3 traits × 2 contrasts)
- Compares moderation vs main-effects-only models (LOO-CV)

**Key outputs**:
- `figures/moderation_interactions.png` - Conditional effects plots (simple slopes)
- `figures/interaction_posteriors.png` - Posterior distributions of interactions
- `results/moderation_summary_table.csv` - Publication-ready coefficients
- `results/moderation_analysis_complete.RData` - Full workspace

**Interpretation focus**:
- Which traits amplify (or buffer) stress reactivity?
- Which traits predict impaired recovery (or resilience)?
- Example: "High NA individuals show β = +X.XX stronger pitch increase"

**Statistical metrics reported**:
- β: Posterior mean of interaction
- 95% CI: Credible interval
- PD: Probability of direction (similar to p-value)
- % in ROPE: Practical equivalence region (±0.1)

---

### STEP 3: Interpretive Report (`manuscript_report.R`)

**Purpose**: Generate manuscript-ready summaries and interpretations

**What it does**:
- Produces executive summary of all findings
- Interprets each interaction in plain language
- Provides manuscript structure recommendations
- Suggests key literatures to cite
- Creates statistical reporting templates

**Key outputs**:
- Console output with comprehensive interpretations
- Guidance for Introduction, Methods, Results, Discussion
- Templates for APA-style statistical reporting

**Use this for**:
- Understanding what the results mean
- Writing the manuscript efficiently
- Ensuring alignment with special issue aims

---

### STEP 4: Utilities (`utilities.R`)

**Purpose**: Create publication materials and conduct robustness checks

**What it does**:
- Generates formatted tables (Table 1, 2, 3)
- Creates diagnostic plots (trace, R-hat, ESS, pp-checks)
- Runs sensitivity analyses:
  - Alternative priors (skeptical vs original)
  - Outlier exclusion
- Calculates standardized effect sizes
- Creates supplementary figures

**Key outputs**:
- `results/table1_sample_characteristics.csv` - Demographics & PID-5
- `results/table2_main_effects.csv` - Stress manipulation results
- `results/table3_moderation_effects.csv` - Interaction effects
- `figures/diagnostics_*.png` - Model convergence checks
- `results/diagnostic_summary.csv` - All models R-hat/ESS
- `results/sensitivity_*.csv` - Robustness checks
- `results/effect_sizes.csv` - Cohen's d for all contrasts

**Use this for**:
- Creating manuscript tables
- Verifying model quality
- Responding to reviewer concerns about robustness

---

## Understanding the Models

### Model Structure

**Level 1 (Within-Person)**:
```
Voice_it = β0i + β1i*(Stress_Contrast) + β2i*(Recovery_Contrast) + eit
```
Each person has their own intercept and slopes for stress/recovery.

**Level 2 (Between-Person)**:
```
β0i = γ00 + γ01*(PID5_NA) + γ02*(PID5_DIS) + γ03*(PID5_PSY) + u0i
β1i = γ10 + γ11*(PID5_NA) + γ12*(PID5_DIS) + γ13*(PID5_PSY) + u1i
β2i = γ20 + γ21*(PID5_NA) + γ22*(PID5_DIS) + γ23*(PID5_PSY) + u2i
```
Individual differences in reactivity/recovery predicted by traits.

**Key Parameters**:
- γ11, γ12, γ13: How NA, DIS, PSY moderate stress reactivity
- γ21, γ22, γ23: How NA, DIS, PSY moderate recovery

---

## Interpreting Results

### Bayesian Metrics

1. **β (Posterior Mean)**
   - Point estimate of effect size
   - Interpreted like a regression coefficient
   - In standardized units (PID-5 centered, so β = effect per 1 SD)

2. **95% Credible Interval**
   - Range containing true value with 95% probability
   - If excludes 0 → "credible evidence" for effect
   - More intuitive than frequentist CI

3. **Probability of Direction (PD)**
   - % of posterior on same side of zero
   - PD > 97.5% ≈ p < .05 (two-tailed)
   - PD > 95% = "moderate evidence"
   - PD > 99.9% = "very strong evidence"

4. **Region of Practical Equivalence (ROPE)**
   - % of posterior in [-0.1, +0.1] range
   - If < 5% in ROPE → effect is "practically significant"
   - Distinguishes statistical from practical significance

### Effect Size Interpretation

For **interaction coefficients**:
- Small: |β| ≈ 0.1 (10% SD change per 1 SD of moderator)
- Medium: |β| ≈ 0.3
- Large: |β| ≈ 0.5

**Example interpretation**:
> "Negative Affectivity moderated stress reactivity in F0 mean, β = 0.25, 
> 95% CI [0.10, 0.40], PD = 99.8%. Individuals 1 SD above the mean in NA 
> showed a 0.25 SD larger increase in vocal pitch during exam stress, 
> representing a meaningful amplification of physiological stress response."

---

## Figures for Manuscript

### Main Figures

1. **Figure 1**: `acoustic_trajectories.png`
   - Shows raw trajectories (spaghetti + means)
   - Caption: "Voice acoustic features across exam stress period"

2. **Figure 2**: `moderation_interactions.png`
   - Conditional effects (simple slopes at ±1 SD)
   - Caption: "Personality pathology moderates stress reactivity"

3. **Figure 3**: `interaction_posteriors.png`
   - Forest plot of all interactions
   - Caption: "Posterior distributions of moderation effects"

### Supplementary Figures

4. **SF1**: `posterior_stress_effects.png` - Main effects
5. **SF2**: `diagnostics_*.png` - Model convergence
6. **SF3**: `supplementary_trajectories_*.png` - Individual differences

---

## Troubleshooting

### Issue: Models won't converge (R-hat > 1.01)

**Solutions**:
1. Increase `adapt_delta`: Change from 0.95 to 0.99
   ```r
   control = list(adapt_delta = 0.99)
   ```

2. Increase `max_treedepth`: Change from 10 to 15
   ```r
   control = list(max_treedepth = 15)
   ```

3. Run longer chains:
   ```r
   iter = 6000, warmup = 3000
   ```

4. Reparameterize random effects:
   ```r
   (1 + timepoint || ID)  # Independent random effects
   ```

### Issue: Out of memory

**Solutions**:
- Reduce parallel chains: `chains = 2, cores = 2`
- Fit models sequentially rather than all at once
- Reduce `iter` (but keep warmup at 50%)

### Issue: Taking too long

**Expected runtimes** (4 cores, default settings):
- Main effects models: ~5 min each
- Moderation models: ~15 min each
- Full pipeline: ~2-3 hours

**Speed up**:
- Use only 2 chains (less reliable but faster)
- Reduce `iter` to 2000 (1000 warmup)
- Fit only key models (e.g., F0 mean + F0 SD)

---

## Extending the Analysis

### Adding Negative Affect EMA as Mediator

If you have momentary negative affect:

```r
# Test if voice changes mediate stress → NA
m_mediation <- brm(
  na_ema ~ timepoint + f0_mean + (1 | ID),
  data = df_analysis,
  ...
)

# Use marginaleffects for indirect effects
library(marginaleffects)
med_effects <- avg_slopes(m_mediation, 
                         variables = "timepoint",
                         by = "f0_mean")
```

### Testing Other PID-5 Domains

Currently focuses on NA, Disinhibition, Psychoticism. To add Detachment/Antagonism:

```r
# In moderation models, add:
+ c1_stress * pid5_det_bl_c + 
+ c1_stress * pid5_ant_bl_c + 
+ c2_recovery * pid5_det_bl_c + 
+ c2_recovery * pid5_ant_bl_c
```

### Analyzing Other Vowels

Data includes /a/, /i/, /u/. Currently uses /a/. To compare:

```r
# Reshape to long format with vowel as factor
df_long_vowel <- df_wide %>%
  pivot_longer(
    cols = starts_with("F0 mean Hz /"),
    names_to = "vowel",
    values_to = "f0_mean"
  ) %>%
  mutate(vowel = str_extract(vowel, "/[aiu]/"))

# Add vowel as predictor
m_vowel <- brm(
  f0_mean ~ timepoint * vowel * pid5_na_bl_c + 
            (1 + timepoint | ID),
  ...
)
```

---

## Manuscript Checklist

### Before Submission

- [ ] All models converged (R-hat < 1.01)
- [ ] All ESS ratios > 0.1
- [ ] Posterior predictive checks look reasonable
- [ ] Tables 1-3 formatted
- [ ] All figures have clear captions
- [ ] Effect sizes calculated and reported
- [ ] Sensitivity analyses show robustness
- [ ] Supplementary materials prepared:
  - [ ] Model specifications (priors, formulas)
  - [ ] Full posterior summaries
  - [ ] Diagnostic plots
  - [ ] Additional vowel analyses (if applicable)

### Key Messages for Special Issue

**Methodological Innovation** (Primary Contribution):
- Voice acoustics as passive, objective stress marker
- Goes beyond traditional self-report EMA
- Captures implicit physiological processes
- Scalable and non-intrusive

**Theoretical Advancement**:
- Demonstrates context-dependent expression of personality pathology
- Validates transactional models (person × situation)
- Shows individual differences in stress reactivity
- Identifies moderators of adaptive vs. maladaptive response

**Clinical Implications**:
- Personalizes risk assessment (who is most vulnerable?)
- Informs intervention targets (reduce stress vs. build resilience)
- Enables real-world monitoring in natural contexts

---

## Support & Contact

Questions about:
- **Statistical approach**: Check `manuscript_report.R` interpretations
- **Code issues**: Review this README's troubleshooting section
- **Conceptual interpretation**: See discussion templates in report script

Good luck with the submission! 🎯

---

## Citation

If you use or adapt this code:

```
[Your citation here once published]

Analysis code: Bayesian multilevel modeling of voice acoustics and 
personality pathology. R with brms (Bürkner, 2017) and cmdstanr 
(Stan Development Team, 2024).
```

---

## Version History

- v1.0 (2025-XX-XX): Initial analysis pipeline
  - Main effects models
  - Moderation by NA, Disinhibition, Psychoticism
  - Visualization suite
  - Diagnostic & sensitivity analyses
  