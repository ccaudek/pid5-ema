# Predictive value of EMA-based PID-5 components for psychological fragility

*Manuscript-ready methods & results (concise), followed by a detailed Supplement.*

---

## 1) Manuscript-ready summary (concise)

### Research question
Do Ecological Momentary Assessment (EMA) measures aligned with the PID‑5 domains **improve the prediction of psychological fragility** during the exam period **beyond** a single baseline administration of the full PID‑5?

### Outcome (“psychological fragility”)
At each EMA prompt participants reported four affective items on 1–7 scales: *happy, sad, satisfied, angry*. We constructed a momentary fragility index as
\[
\text{Fragility}_{t} = (7-\texttt{happy}_{t}) + (7-\texttt{satisfied}_{t}) + \texttt{sad}_{t} + \texttt{angry}_{t},
\]
then standardized it (z-score) across all observations. Larger values denote higher negative affect / fragility.

### Predictors
- **Baseline (between-person)**: PID‑5 *Negative Affect, Detachment, Antagonism, Disinhibition, Psychoticism* measured once at study start (first non-missing per subject); each domain was standardized across subjects. Sex (female vs other) was included as a binary covariate. Exam **period** was coded with two dummies (*pre* and *post*, baseline as reference).
- **EMA (within-person)**: Short EMA measures aligned to PID‑5 components (Detachment, Antagonism, Disinhibition, Psychoticism) collected at each prompt and standardized at the observation level. These were included **only** in Model **B** (see below).

### Models
We estimated hierarchical Gaussian mixed-effects models with subject-level random effects and subject-specific residual scale:
\[
y_n \sim \mathcal{N}\big(x_n^\top \beta + z_n^\top b_{s(n)}, \ \sigma_{s(n)}^2\big),\quad
b_s \sim \mathcal{N}\!\left(0,\ \mathrm{diag}(\tau_1^2,\dots,\tau_R^2)\right),\quad
\log \sigma_s = \log \sigma + \delta_s,\ \ \delta_s \sim \mathcal{N}(0,\sigma_\delta^2).
\]
Here \(x_n\) are fixed effects and \(z_n\) are random-effects design vectors (subject-specific intercept and random slopes for the exam-period dummies, plus EMA slopes in Model B). Predictors were centered and scaled. Weakly informative priors were used (Normal(0,1) for standardized betas; half-Normals for scale parameters).

Two nested models were compared:
- **Model A (baseline-only)**: baseline PID‑5 domains + sex + exam-period dummies.
- **Model B (baseline + EMA)**: Model A + EMA PID‑5 components (within-person). Random slopes were specified for the period indicators and EMA terms.

### Estimation and out-of-sample comparison
We used **CmdStanR** with 4 chains (1000 warmup + 1000 sampling per chain). The primary metric was **subject-wise PSIS-LOO** (leave-one-subject-out) computed from the **marginal subject likelihood** (random effects integrated in closed form) with **moment matching** correction. Reliability was checked via Pareto-\(k\) diagnostics; we also computed \(r_\mathrm{eff}\) from chain IDs when appropriate.

### Results (model comparison)
The analysis included **N = 6,226** observations from **S = 229** participants. Model **B** showed a **large and precise** improvement in subject-wise predictive accuracy over Model **A**:
\[
\Delta \text{ELPD}_{\text{subject}} = \text{ELPD}_B - \text{ELPD}_A = \mathbf{134.1} \pm \mathbf{9.8}\ \text{SE} \quad (\text{PSIS-LOO with moment matching}),
\]
with **all Pareto-\(k\) < 0.7** for both models. This corresponds to a decrease in LOOIC of ≈ **268** points in favor of Model B.

Posterior predictive checks (PPC) indicated good calibration (ECDF and QQ overlays) and, critically, the model with subject-specific residual scales reproduced the **distribution of within-subject standard deviations**. A visible discrepancy in the global density overlay reflects the discretized 1–7 construction of the outcome rather than poor fit.

### Conclusion
EMA-based PID‑5 components provide **substantial out-of-sample predictive value** for momentary psychological fragility **over and above** a one‑off baseline PID‑5 assessment. Given the magnitude and robustness of the LOO improvement, EMA information meaningfully enhances prediction at the **subject level**, which is aligned with the study’s primary inferential target.

---

## 2) Supplementary Methods (detailed)

### 2.1 Data processing
- **Scales.** EMA items (*happy, satisfied, sad, angry*) were coerced to 1–7 integers where needed (strings parsed to numbers; 0–100 scales linearly mapped to 1–7 when present).
- **Fragility index.** For each observation \(t\): \((7-\text{happy}) + (7-\text{satisfied}) + \text{sad} + \text{angry}\). The composite was **z-scored** over all observations.
- **Baseline PID‑5 domains.** For each subject, we took the **first non-missing** value per domain at study start and then standardized across subjects.
- **EMA PID‑5 components.** Detachment, Antagonism, Disinhibition, Psychoticism short EMA measures were standardized **within the full sample** at the observation level.
- **Period coding.** Three phases (baseline/pre/post) were coded as two indicator variables (*pre*, *post*; baseline is reference). These dummies entered both as fixed and random slopes.
- **Sex.** A binary female indicator was created from available information.
- **Alignment.** All predictors were constructed in a single long-format table (one row per observation), then split into fixed-effects \(X\) and random-effects \(Z\) matrices to ensure **perfect alignment** across models A and B.

### 2.2 Statistical model
For observation \(n=1,\dots,N\) with subject \(s=s(n)\):
\[
\begin{aligned}
y_n \mid b_{s(n)}, \sigma_{s(n)} &\sim \mathcal{N}\!\left(\mu_n,\ \sigma_{s(n)}^2\right),\quad
\mu_n = x_n^\top \beta + z_n^\top b_{s(n)}, \\
b_s &\sim \mathcal{N}\!\left(0,\ \mathrm{diag}(\tau_1^2,\dots,\tau_R^2)\right), \\
\log \sigma_s &= \log \sigma + \delta_s,\quad \delta_s \sim \mathcal{N}(0,\sigma_\delta^2).
\end{aligned}
\]
- \(x_n\) (fixed effects) included baseline PID‑5 domains, sex, and period dummies; in Model B, EMA components were added.
- \(z_n\) (random effects) included an intercept and random slopes for period; in Model B we additionally allowed random slopes for the EMA regressors (to capture between-subject heterogeneity in EMA–fragility coupling).
- **Priors.** For standardized predictors, \(\beta_j \sim \mathcal{N}(0,1)\). Random-effect SDs \(\tau_r \sim \mathrm{half\text{-}Normal}(0,1)\). Overall residual SD \(\sigma \sim \mathrm{half\text{-}Normal}(0,1)\). Heteroscedasticity parameter \(\sigma_\delta \sim \mathrm{half\text{-}Normal}(0,1)\). These weakly-informative priors regularize the model while keeping estimates primarily likelihood-driven.

> **Note.** The outcome combines ordinal inputs but is treated as approximately continuous. This choice is justified by the composite’s near-continuous support and the focus on **predictive** performance at the subject level; nevertheless, we examined posterior predictive discrepancies (Section 2.6).

### 2.3 Subject-wise marginal likelihood
Let \(y_s\) be the vector of outcomes for subject \(s\) and \(X_s, Z_s\) be the design blocks. Conditional on \(\sigma_s\), integrating the Gaussian random effects yields the exact marginal:
\[
y_s \mid \beta,\sigma_s,\tau \ \sim\ \mathcal{N}\!\left(X_s\beta,\ \Sigma_s\right), \quad
\Sigma_s \equiv Z_s \,\mathrm{diag}(\tau_1^2,\dots,\tau_R^2)\, Z_s^\top \;+\; \sigma_s^2 I_{n_s}.
\]
Hence the **subject-level log-likelihood** is
\[
\log p(y_s\mid\theta)= -\tfrac12\!\left[n_s\log(2\pi)+\log|\Sigma_s|+(y_s-X_s\beta)^\top \Sigma_s^{-1}(y_s-X_s\beta)\right],
\]
which we compute in Stan for all \(s=1,\dots,S\). This **marginal** form is used both in the target density (through the conditional representation) and to report the **per-subject log-likelihoods** used by PSIS-LOO.

### 2.4 Estimation
- **Engine.** CmdStanR, 4 chains, 1000 warmup + 1000 sampling iterations, `adapt_delta = 0.95`.
- **Data scaling.** All continuous predictors were centered and scaled. The two period indicators were also used unscaled in the heteroscedastic submodel (binary 0/1) to define subject-specific \(\sigma_s\) where needed.
- **Random effects.** We used independent (diagonal) subject-level random-effect SDs \(\tau_r\); empirical checks suggested this captured the dominant heterogeneity without the additional cost/instability of a full correlation matrix at this sample size.

### 2.5 Out-of-sample comparison
We targeted the **subject** as the unit of generalization. We computed **PSIS-LOO** on the \(S\)-vector of per-subject log-likelihoods \(\{\log p(y_s\mid\theta^{(m)})\}_{m=1}^M\), using **moment matching** to stabilize any influential terms. We passed chain-based `r_eff` to `loo::loo()` to reflect MCMC autocorrelation. We report the difference in expected log predictive density (ELPD) and its standard error:
\[
\Delta \mathrm{ELPD} = \mathrm{ELPD}_B - \mathrm{ELPD}_A.
\]
Positive values favor Model B. We also confirm diagnostic reliability by inspecting the distribution of **Pareto-\(k\)** values.

### 2.6 Posterior predictive checks
We conducted PPCs focusing on distributional calibration and within-subject dynamics:
1. **Global density overlay** of \(y\) vs replicated \(y^\mathrm{rep}\). Small discrepancies persisted—as expected—because \(y\) is a discretized composite of ordinal inputs.
2. **ECDF overlay** and **QQ overlays** showed good overall calibration.
3. **Distribution of within-subject SDs** \(\{\mathrm{sd}(y_{s\cdot})\}\) was well reproduced **only** when allowing subject-specific residual scales \(\sigma_s\) (heteroscedastic model), supporting this component of the specification.
4. **Subject × period SDs** were also examined by computing SDs within each subject-period cell (retaining cells with ≥2 observations); replicated distributions were close to observed ones.

### 2.7 Results (detailed)
- Data: **N = 6,226** observations; **S = 229** subjects.
- **Model comparison.** PSIS-LOO at the subject level strongly favored **Model B**:
  - \(\Delta \mathrm{ELPD} = \mathbf{134.1} \pm \mathbf{9.8}\) SE (Model B – Model A).
  - All **Pareto-\(k < 0.7\)** for both models (with moment matching), indicating reliable approximations.
- **Diagnostics.** Chains mixed well; effective sample sizes were adequate for the key parameters (details available from `posterior::summarise_draws`). PPCs supported distributional adequacy for the main targets (Section 2.6).
- **Interpretation.** Incorporating **EMA-based PID‑5 components** meaningfully improves **subject-level** out-of-sample prediction of momentary fragility beyond baseline-only information, consistent with the hypothesis that **dynamic, within-person personality-affect signals** carry incremental predictive value during the exam period.

### 2.8 Sensitivity & robustness
We explored several alternatives during model development:
- **Residual distribution.** A Student‑t model improved tail fit but left the main conclusion unchanged (Model B > Model A).
- **Inference scheme.** Variational inference (meanfield) and short MCMC both yielded **consistent** model rankings, though final inferences rely on MCMC.
- **Random-effects structure.** Adding random slopes (period and EMA) and subject-specific \(\sigma_s\) improved PPCs (notably within-subject SDs) and PSIS diagnostics.
Across these variations, **Model B consistently dominated** in subject-wise predictive accuracy.

### 2.9 Limitations
- The outcome aggregates ordinal items and is treated as continuous; a full ordinal-IRT measurement model could be pursued in future work.
- We used independent random-effect SDs for computational stability; a full covariance may capture additional structure but at higher cost.
- Heteroscedasticity was modeled at the subject (not time-varying) level to target the LOO-subject estimand; extensions to time-varying scale are possible.

### 2.10 Reproducibility
- Data file: `data/processed/ema_plus_scales_cleaned.csv`.
- Script (R): `fragility_subject_loo_rslope_sigmaS.R` (provided).
- Model: `subject_loo_normal_rslope_sigmaS.stan`.
- Seed: `20250919`. CmdStanR 4 chains × (1000 warmup + 1000 sampling), `adapt_delta = 0.95`.
- All preprocessing ensures that **Model A** and **Model B** are built from the **same aligned rows** and differ only by the inclusion of EMA terms.

---

*End of document.*