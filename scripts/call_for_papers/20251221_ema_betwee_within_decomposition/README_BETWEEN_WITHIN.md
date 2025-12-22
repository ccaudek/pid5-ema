# BETWEEN-WITHIN DECOMPOSITION ANALYSIS
## Voice-Personality Study: Trait vs State Components

---

## OVERVIEW

This analysis suite tests whether **within-person state fluctuations** in personality pathology predict voice characteristics **beyond stable between-person traits**.

**Research Question:**
> Do momentary personality states (measured via EMA at voice recording timepoints) contribute to vocal stress responses beyond what stable traits can explain?

**Model:**
```
Voice[t] = Trait_between[constant] + State_within[varying] + Timepoint + Interactions
```

---

## CONCEPTUAL FRAMEWORK

### Current Analysis (Trait-Only)
- **Predictors:** EMA PID-5 aggregated across ~27 assessments (stable trait)
- **Outcomes:** Voice features at 3 timepoints (baseline, pre-exam, post-exam)
- **Question:** Do stable traits moderate stress reactivity?

### New Analysis (Trait + State)  
- **Between-person (Trait):** Person's average across all EMA assessments
- **Within-person (State):** Deviation from person mean at each timepoint
- **Question:** Do momentary states add predictive value beyond traits?

### Why This Matters
If within-person effects are credible:
1. **Theoretical:** Demonstrates state-trait dissociation in personality-physiology coupling
2. **Practical:** Suggests real-time assessment could improve prediction
3. **Clinical:** Momentary interventions could target state fluctuations

---

## FILE STRUCTURE

```
between_within_analysis/
├── 11_prepare_between_within_data.R      # Data preparation & decomposition
├── 12_fit_between_within_models.R        # Fit Bayesian models (F0 /a/ first)
├── 13_compare_between_within.R           # Model comparison (trait vs trait+state)
├── 14_visualize_between_within.R         # Create publication figures
├── 15_extend_to_all_outcomes.R           # Extend to all 18 outcomes (if promising)
└── README.md                             # This file
```

**Outputs:**
```
results/between_within/
├── df_between_within.rds                 # Prepared data with decomposition
├── f0_mean_a_comparison.rds              # Between vs within parameter estimates
├── model_comparison_summary.rds          # R², LOO-IC, decision
├── all_outcomes_comparison.rds           # Results from all 18 models
├── within_effects_summary.csv            # Summary by outcome
└── publication_table.csv                 # Formatted table for manuscript

figures/between_within/
├── variance_decomposition.png            # ICC for each domain
├── forest_plot_main_effects.png          # Between vs within comparison
├── model_comparison_r2.png               # R² comparison
├── combined_summary.png                  # All panels together
└── interaction_*.png                     # Credible state × time effects

models/between_within/
├── f0_mean_a_bw.rds                      # Main model
└── all_outcomes/                         # All 18 models (if extended)
```

---

## PIPELINE

### STEP 1: Data Preparation (Script 11)

**What it does:**
1. Loads EMA data with PID-5 scores at all timepoints
2. Matches EMA assessments to voice recording dates
3. Decomposes PID-5 into between-person and within-person components:
   - **Between:** Mean across all ~27 EMA assessments (stable trait)
   - **Within:** Deviation from person mean at each timepoint (state)
4. Verifies decomposition (within-person deviations should average to 0)
5. Centers predictors and adds contrast codes

**Output:** `df_between_within.rds` with variables:
- `pid5_*_between_c`: Grand-mean centered traits (constant per person)
- `pid5_*_within_c`: Person-mean centered states (vary across timepoints)

**Run:**
```r
source("11_prepare_between_within_data.R")
```

**Expected time:** 2-3 minutes

**Key checks:**
- N observations = 423 (141 subjects × 3 timepoints)
- ICC values: High ICC = mostly trait, Low ICC = mostly state
- Decomposition verification: Max deviation should be < 0.001

---

### STEP 2: Fit Models (Script 12)

**What it does:**
1. Fits Bayesian multilevel model for F0 mean /a/ with:
   - Between-person main effects + interactions
   - Within-person main effects + interactions
2. Extracts posterior distributions for all effects
3. Compares between vs within effects for each domain
4. Identifies credible within-person effects (95% CI excludes 0)
5. Provides decision criteria for manuscript inclusion

**Model specification:**
```r
f0_mean_a ~ 
  # Main effects
  pid5_*_between_c + pid5_*_within_c + c1_stress + c2_recovery +
  
  # Between × Timepoint interactions
  pid5_*_between_c:c1_stress + pid5_*_between_c:c2_recovery +
  
  # Within × Timepoint interactions (KEY!)
  pid5_*_within_c:c1_stress + pid5_*_within_c:c2_recovery +
  
  # Random effects
  (1 + c1_stress + c2_recovery || ID)
```

**Output:** 
- Model: `models/between_within/f0_mean_a_bw.rds`
- Comparison: `results/between_within/f0_mean_a_comparison.csv`

**Run:**
```r
source("12_fit_between_within_models.R")
```

**Expected time:** 10-15 minutes

**Decision criteria:**
- **≥2 credible within-person effects:** Include in main manuscript
- **1 credible effect:** Include in supplementary materials
- **0 credible effects:** Supplementary only, trait model sufficient

---

### STEP 3: Model Comparison (Script 13)

**What it does:**
1. Compares trait-only model (from Script 03_FINAL) vs trait+state model
2. Metrics:
   - **Bayesian R²:** Variance explained
   - **LOO-IC:** Out-of-sample predictive accuracy
   - **WAIC:** Alternative information criterion
   - **Parameter efficiency:** Effective number of parameters
3. Provides formal recommendation based on:
   - Statistical evidence (R², LOO)
   - Number of credible within-person effects
   - Parsimony considerations

**Output:** `model_comparison_summary.rds` with decision

**Run:**
```r
source("13_compare_between_within.R")
```

**Expected time:** 3-5 minutes

**Interpretation:**
- **ΔR² > 0.01 AND LOO favors trait+state:** Strong evidence for state
- **ΔR² < 0.005 OR LOO equivalent:** Trait-only sufficient
- **Intermediate:** Consider theoretical value

---

### STEP 4: Visualization (Script 14)

**What it does:**
Creates publication-quality figures:
1. **Variance decomposition:** ICC for each domain (trait vs state %)
2. **Forest plot:** Between vs within effects with CIs
3. **Interaction plots:** Credible state × timepoint effects
4. **Model comparison:** R² comparison bar chart
5. **Combined summary:** All panels in one figure

**Output:** Multiple PNG files in `figures/between_within/`

**Run:**
```r
source("14_visualize_between_within.R")
```

**Expected time:** 2-3 minutes

**Figures suitable for:**
- Main manuscript (if results are strong)
- Supplementary materials (if results are moderate)
- Presentations/talks (always useful!)

---

### STEP 5: Extension to All Outcomes (Script 15) - OPTIONAL

**Only run if:**
- F0 results showed ≥1 credible within-person effect, OR
- Model comparison favored trait+state model

**What it does:**
1. Fits trait+state models for all 18 outcomes in parallel:
   - 6 features: f0_mean, f0_std, f2_mean, f2_std, jitter, nne
   - 3 vowels: /a/, /i/, /u/
2. Extracts within-person effects from all models
3. Identifies outcomes with most robust state effects
4. Provides final recommendation for manuscript

**Output:** 
- 18 models in `models/between_within/all_outcomes/`
- Summary: `all_outcomes_comparison.csv`
- Master table: `all_credible_within_effects.csv`

**Run:**
```r
source("15_extend_to_all_outcomes.R")
```

**Expected time:** 40-60 minutes (parallel on 4 cores)

**Decision:**
- **≥10 total within-person effects:** Strong case for main manuscript
- **5-9 effects:** Moderate case, main or extensive supplement
- **<5 effects:** Supplementary only

---

## INTERPRETATION GUIDE

### Between-Person Effects (Trait)
- **What it means:** Individuals who are generally high in this trait (averaged across all assessments) show this vocal pattern
- **Example:** "People high in Negative Affectivity have higher baseline F0"
- **Level:** Individual differences (stable)

### Within-Person Effects (State)
- **What it means:** When an individual's personality state is elevated (relative to their own average), their voice changes
- **Example:** "On days when someone feels more detached than usual, their F0 recovery is impaired"
- **Level:** Intra-individual fluctuation (dynamic)

### Why Both Matter
**Between-person (trait):**
- Identifies vulnerability factors (who is at risk)
- Stable individual differences
- Already well-studied

**Within-person (state):**
- Identifies dynamic processes (when/how risk manifests)
- Malleable targets for intervention
- Often overlooked in voice research

---

## DECISION TREE FOR MANUSCRIPT INCLUSION

```
                    START
                      |
                      v
          Fit F0 model (Script 12)
                      |
                      v
        ┌─────────────┴─────────────┐
        |                           |
    ≥2 credible              0-1 credible
    within effects           within effects
        |                           |
        v                           v
   Run Script 13                Run Script 13
   (model comp)                 (model comp)
        |                           |
        v                           v
  R² improvement            No R² improvement
    ΔR² > 0.01                  ΔR² < 0.005
        |                           |
        v                           v
   Run Script 15              STOP: Trait-only
   (extend all)               sufficient
        |                           |
        v                           v
  ≥10 total effects          Supplementary
        |                    materials only
        v
   MAIN MANUSCRIPT
   - Methods: Full decomposition
   - Results: Trait+State section
   - Discussion: State dynamics
        |
        v
   Create figures (Script 14)
```

---

## EXPECTED RESULTS SCENARIOS

### Scenario A: Strong State Effects (Ideal)
- **F0 /a/:** 3-4 credible within-person effects
- **Model comparison:** ΔR² = 0.03, LOO favors trait+state
- **Extension:** 12+ total effects across outcomes
- **Action:** Full integration in main manuscript

**Implications:**
- State fluctuations are meaningful beyond traits
- Momentary assessment adds substantial value
- Theoretical contribution: Dynamic personality processes

### Scenario B: Moderate State Effects
- **F0 /a/:** 1-2 credible within-person effects  
- **Model comparison:** ΔR² = 0.01, LOO equivalent
- **Extension:** 5-8 total effects
- **Action:** Supplementary materials with main text mention

**Implications:**
- Some state variance exists but limited
- Trait model tells main story
- State analysis is interesting but not essential

### Scenario C: Null State Effects
- **F0 /a/:** 0 credible within-person effects
- **Model comparison:** ΔR² < 0.005, LOO favors trait-only
- **Extension:** <3 total effects
- **Action:** Brief supplementary note only

**Implications:**
- Personality-voice coupling is primarily trait-based
- EMA states do not add predictive value
- Trait aggregation across timepoints is appropriate

---

## TECHNICAL NOTES

### Data Requirements
1. **EMA PID-5 at voice timepoints:** Must have personality assessments on or near (±3 days) voice recording dates
2. **Minimum observations:** At least 2 timepoints per person (prefer 3)
3. **Between-person variance:** Need sufficient variability in person means

### Potential Issues

**Issue 1: Low within-person variance**
- **Symptom:** ICC > 0.95 for all domains
- **Cause:** Personality very stable during exam period
- **Solution:** This is a substantive finding! Report high stability, don't force within-person analysis

**Issue 2: Date matching failures**
- **Symptom:** Many missing EMA-voice matches
- **Cause:** EMA assessments not collected near voice recordings
- **Solution:** Script 11 tries both date-matching and exam_period approach

**Issue 3: Convergence problems**
- **Symptom:** Rhat > 1.01, low ESS
- **Cause:** Complex model with limited within-person data
- **Solution:** Increase adapt_delta to 0.995, iterations to 6000

**Issue 4: Overfitting**
- **Symptom:** p_loo >> number of parameters
- **Cause:** Model too complex for data
- **Solution:** Use trait-only model, report in supplement

### Assumptions
1. **Linearity:** Effects of personality on voice are linear
2. **Independence:** Within-person deviations independent across timepoints
3. **No carry-over:** State at T1 doesn't affect state at T2
4. **Measurement equivalence:** PID-5 items measure same construct across time

---

## MANUSCRIPT INTEGRATION

### If Results Support Inclusion

**Methods additions:**
```markdown
## Between-Person and Within-Person Decomposition

To test whether momentary personality states contribute beyond stable 
traits, we decomposed EMA PID-5 scores into between-person (person's 
average across all assessments) and within-person (deviation from 
person mean at each timepoint) components. This allows simultaneous 
estimation of:

1. Between-person effects: Whether individuals high in a trait show 
   different vocal patterns
2. Within-person effects: Whether fluctuations around one's typical 
   level predict vocal changes

The model was:
Voice[t] = β_between × Trait + β_within × State + ...
```

**Results additions:**
```markdown
## State-Trait Decomposition

Beyond stable trait effects, within-person state fluctuations showed 
[N] credible associations with vocal features. Specifically, when 
individuals experienced elevated [Domain] relative to their own 
average, [vocal effect]. This suggests that personality-voice coupling 
reflects both stable individual differences and dynamic state processes.
```

**Discussion additions:**
```markdown
The significant within-person effects demonstrate that acoustic 
biomarkers capture not only trait-level vulnerabilities but also 
real-time psychological states. This has implications for:

1. Ambulatory monitoring: Real-time personality assessment feasible
2. Intervention timing: Target moments of elevated vulnerability
3. Theory: Personality pathology involves dynamic state processes
```

### If Results Don't Support Inclusion

**Brief supplementary note:**
```markdown
## Supplementary Note: State-Trait Decomposition

We tested whether within-person state fluctuations in PID-5 domains 
(deviations from person mean) contributed beyond between-person trait 
differences. Model comparison (ΔR² = [X], LOO-IC difference = [Y]) 
showed [no substantial improvement/equivalent fit]. Trait aggregation 
across EMA assessments appears sufficient for capturing personality-
voice relationships in this sample. Full results available upon request.
```

---

## CONTACT & SUPPORT

**Questions about:**
- **Conceptual framework:** Review Curran & Bauer (2011) on within-between models
- **Bayesian implementation:** See brms documentation on multilevel models
- **Interpretation:** Email study PI or consult with methodologist

**Common problems:**
- Check TROUBLESHOOTING.txt (if created)
- Review convergence diagnostics
- Consult with Bayesian statistics expert if needed

---

## REFERENCES

**Methodological:**
- Curran, P. J., & Bauer, D. J. (2011). The disaggregation of within-person and between-person effects in longitudinal models of change. *Annual Review of Psychology, 62*, 583-619.
- Hoffman, L. (2015). *Longitudinal analysis: Modeling within-person fluctuation and change*. Routledge.
- McNeish, D., & Hamaker, E. L. (2020). A primer on two-level dynamic structural equation models for intensive longitudinal data in Mplus. *Psychological Methods, 25*(5), 610-635.

**Software:**
- Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software, 80*(1), 1-28.
- Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing, 27*(5), 1413-1432.

---

## VERSION HISTORY

- **v1.0** (2025-12-21): Initial release
  - Scripts 11-15 created
  - F0 mean /a/ as test case
  - Extension to all outcomes

---

**GOOD LUCK WITH YOUR ANALYSIS!**

Remember: 
- Null results are still results
- Trait-only model is perfectly valid
- State analysis is exploratory - let data guide decisions
