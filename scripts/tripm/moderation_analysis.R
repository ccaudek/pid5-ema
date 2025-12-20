# ==============================================================================
# Bayesian Moderation Analysis: PID-5 Traits × Stress Reactivity
# Testing context-dependent expression of personality pathology
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayestestR)
  library(bayesplot)
  library(tidybayes)
  library(patchwork)
  library(ggdist)
  library(marginaleffects)
  library(tidyr)
  library(stringr)
})

options(brms.backend = "cmdstanr")
options(max.print = 5000)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODERATION ANALYSIS: PID-5 × Stress Reactivity\n")
cat("All 5 PID-5 Domains\n")
cat(rep("=", 70), "\n", sep = "")

# ==============================================================================
# 1. PREPARE DATA FOR MODERATION MODELS
# ==============================================================================

# Check PID-5 distributions at baseline
cat("\n=== PID-5 Baseline Trait Distributions (All 5 Domains) ===\n")
df_analysis %>%
  dplyr::filter(timepoint == "baseline") %>%
  dplyr::select(ends_with("_bl_c")) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    M = mean(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Min = min(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

df_analysis$f2_std_i <- ifelse(
  df_analysis$f2_std_i == 0,
  min(df_analysis$f2_std_i[df_analysis$f2_std_i > 0]),
  df_analysis$f2_std_i
)

# Create contrast-coded timepoint variable for interaction interpretation
# C1: PRE vs BASELINE (stress reactivity)
# C2: POST vs PRE (stress recovery)

df_analysis <- df_analysis %>%
  mutate(
    # Contrast 1: PRE vs BASELINE (test stress reactivity)
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre" ~ 0.5,
      timepoint == "post" ~ 0
    ),
    # Contrast 2: POST vs PRE (test recovery)
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre" ~ -0.5,
      timepoint == "post" ~ 0.5
    )
  )


# --- Parametri empirici ricavati dal summary
vowel_means <- list(
  a = list(
    f0_mean = 189.9,
    f0_std = 6.230,
    jitter = 0.8366,
    nne = -24.42,
    f2_mean = 1207.4,
    f2_std = 59.37
  ),
  i = list(
    f0_mean = 194.4,
    f0_std = 8.888,
    jitter = 1.058,
    nne = -27.89,
    f2_mean = 2189.0,
    f2_std = 141.04
  ),
  u = list(
    f0_mean = 195.1,
    f0_std = 9.724,
    jitter = 1.286,
    nne = -28.429,
    f2_mean = 1123.2,
    f2_std = 161.57
  )
)

# ---------------------------------------------------------------------
# make_priors_vowel: versione sicura che usa prior_string() con numeri
# ---------------------------------------------------------------------
make_priors_vowel <- function(outcome, vowel) {
  vm <- vowel_means[[vowel]]
  # helper per formattare numeri (evita notazione scientifica spiacevole)
  fmt <- function(x, digits = 6) formatC(x, format = "f", digits = digits)

  switch(
    outcome,

    "f0_mean" = c(
      prior_string(
        paste0("student_t(3, ", fmt(vm$f0_mean), ", 30)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 10)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "f0_std" = c(
      # intercept on log scale -> log(mean)
      prior_string(
        paste0("student_t(3, ", fmt(log(vm$f0_std)), ", 0.5)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 0.3)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "jitter" = c(
      prior_string(
        paste0("student_t(3, ", fmt(log(vm$jitter)), ", 0.5)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 0.3)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "nne" = c(
      prior_string(
        paste0("student_t(3, ", fmt(vm$nne), ", 5)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 2)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "f2_mean" = c(
      prior_string(
        paste0("student_t(3, ", fmt(vm$f2_mean), ", 150)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 50)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    "f2_std" = c(
      prior_string(
        paste0("student_t(3, ", fmt(log(vm$f2_std)), ", 0.6)"),
        class = "Intercept"
      ),
      prior_string("normal(0, 0.4)", class = "b"),
      prior_string("exponential(1)", class = "sd"),
      prior_string("exponential(1)", class = "sigma")
    ),

    stop("Outcome non riconosciuto")
  )
}

# --- Mappa outcome -> family (senza suffisso vocale) -----------------------
outcome_families <- list(
  f0_mean = gaussian(),
  f0_std = lognormal(),
  jitter = lognormal(),
  nne = gaussian(),
  f2_mean = asym_laplace(),
  f2_std = lognormal()
)

# --- Lista degli outcome di interesse e vocali ------------------------------
outcomes <- c("f0_mean", "f0_std", "jitter", "nne", "f2_mean", "f2_std")
vowels <- c("a", "i", "u")

# --- Opzioni di sampling (usa i tuoi valori) -------------------------------
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
seed <- 123
control <- list(adapt_delta = 0.95, max_treedepth = 12)

# --- Cartella per salvare modelli (opzionale) ------------------------------
dir.create("models", showWarnings = FALSE)

# --- Flag: se TRUE lancia i modelli; se FALSE prepara tutto e stampa il piano -
run_models <- TRUE

# --- Loop per costruire formule e (opzionalmente) eseguire brm ---------------
# --- Loop corretto per costruire formule ed eseguire brm ---------------
fitted_models <- list()
for (v in vowels) {
  message("\n--- Vocale: /", v, "/ -------------------------------")
  for (out in outcomes) {
    # nome colonna effettivo nel df (es. f0_mean_a)
    colname <- paste0(out, "_", v)
    fam <- outcome_families[[out]]
    # safety: controlla che la col esista
    if (!colname %in% names(df_analysis)) {
      message("  Skipping ", colname, " (colonna non trovata nel dataset).")
      next
    }

    # formula: come nel tuo script, interaction con tutti i 5 domini e random slopes
    fmla <- bf(
      as.formula(
        paste0(
          colname,
          " ~ c1_stress * (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c + pid5_antagonism_bl_c + pid5_disinhibition_bl_c + pid5_psychoticism_bl_c) + ",
          "c2_recovery * (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c + pid5_antagonism_bl_c + pid5_disinhibition_bl_c + pid5_psychoticism_bl_c) + ",
          "(1 + c1_stress + c2_recovery | ID)"
        )
      )
    )

    priors_here <- make_priors_vowel(out, v)

    model_name <- paste0("m_", out, "_", v)

    # estrai nome della family in modo sicuro
    fam_name <- NULL
    if (!is.null(fam$family)) {
      fam_name <- fam$family
    } else if (!is.null(attr(fam, "family"))) {
      fam_name <- attr(fam, "family")
    } else {
      fam_name <- paste0(class(fam)[1])
    }

    message("  Model: ", model_name, "  (family: ", fam_name, ")")

    if (run_models) {
      # fitta il modello e salva l'oggetto
      fit <- brm(
        formula = fmla,
        data = df_analysis,
        family = fam,
        prior = priors_here,
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores,
        seed = seed,
        control = control,
        file = file.path("models", model_name) # salva su disco per riuso
      )
      fitted_models[[model_name]] <- fit

      # rapido check diagnostico (opzionale)
      print(summary(fit, waic = FALSE, odds = FALSE))
      pp_check(fit) # puoi commentare se vuoi meno output
    } else {
      # se non si esegue, mostriamo i priors che useremmo
      message("    (run_models = FALSE) - priors che verrebbero usati:")
      print(priors_here)
    }
  } # fine loop outcomes
} # fine loop vowels


# --- Uso degli oggetti fitted_models: se run_models = TRUE -> fitted_models popolato
# --- Fine

# --- 1. Estrarre summary da tutti i modelli -------------------------
results_table <- lapply(names(fitted_models), function(mn) {
  fit <- fitted_models[[mn]]
  fe <- fixef(fit, summary = TRUE) # posterior mean, sd, 95% CI

  # determina outcome e vocale dal nome del modello
  # es. m_f0_mean_a -> outcome = f0_mean, vowel = a
  parts <- str_split(mn, "_", simplify = TRUE)
  outcome <- parts[2]
  vowel <- parts[3]

  # identifica se il termine è main o interazione
  terms <- rownames(fe)
  type <- ifelse(str_detect(terms, ":"), "Interaction", "Main")

  tibble(
    Outcome = outcome,
    Vowel = vowel,
    Term = terms,
    Type = type,
    Estimate = fe[, "Estimate"],
    SE = fe[, "Est.Error"],
    CI_low = fe[, "Q2.5"],
    CI_high = fe[, "Q97.5"],
    Credible = ifelse(fe[, "Q2.5"] > 0 | fe[, "Q97.5"] < 0, TRUE, FALSE)
  )
}) %>%
  bind_rows() %>%
  arrange(Outcome, Vowel, Type, Term)

# Visualizza tabella completa
data.frame(results_table)


# --- 2. Forest plot bayesiano -------------------------
ggplot(
  results_table,
  aes(
    x = Estimate,
    y = Term,
    xmin = CI_low,
    xmax = CI_high,
    color = Type,
    shape = Credible
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  facet_grid(Outcome ~ Vowel, scales = "free_y") +
  scale_color_manual(values = c("Main" = "blue", "Interaction" = "red")) +
  scale_shape_manual(values = c(FALSE = 1, TRUE = 16)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Posterior Estimate",
    y = "Term",
    color = "Type",
    shape = "Credible (CI not including 0)",
    title = "Bayesian Forest Plot: Effects and Interactions"
  ) +
  theme(strip.text.y = element_text(angle = 0))


data.frame(results_table)
# Outcome Vowel                                       Term        Type      Estimate          SE        CI_low       CI_high Credible
# 1        f0  mean             c1_stress:pid5_antagonism_bl_c Interaction  1.294138e+00  1.46430754 -1.549582e+00  4.151366e+00    FALSE
# 2        f0  mean             c1_stress:pid5_antagonism_bl_c Interaction -4.659519e-01  1.75855417 -3.966726e+00  2.997478e+00    FALSE
# 3        f0  mean             c1_stress:pid5_antagonism_bl_c Interaction  1.248467e+00  1.69135151 -2.052353e+00  4.644937e+00    FALSE
# 4        f0  mean             c1_stress:pid5_detachment_bl_c Interaction -1.006429e+00  1.63481214 -4.211004e+00  2.210809e+00    FALSE
# 5        f0  mean             c1_stress:pid5_detachment_bl_c Interaction  2.297715e-01  1.99687954 -3.714312e+00  4.152842e+00    FALSE
# 6        f0  mean             c1_stress:pid5_detachment_bl_c Interaction -1.199222e+00  1.84843540 -4.864950e+00  2.433974e+00    FALSE
# 7        f0  mean          c1_stress:pid5_disinhibition_bl_c Interaction -3.512273e-01  1.79282704 -3.913807e+00  3.189504e+00    FALSE
# 8        f0  mean          c1_stress:pid5_disinhibition_bl_c Interaction -1.239840e+00  2.11716036 -5.469902e+00  2.864737e+00    FALSE
# 9        f0  mean          c1_stress:pid5_disinhibition_bl_c Interaction -9.791346e-02  2.04556007 -4.167068e+00  3.844431e+00    FALSE
# 10       f0  mean   c1_stress:pid5_negative_affectivity_bl_c Interaction  1.126961e+00  1.60488143 -2.098278e+00  4.302168e+00    FALSE
# 11       f0  mean   c1_stress:pid5_negative_affectivity_bl_c Interaction  5.032815e+00  1.92072258  1.259305e+00  8.762684e+00     TRUE
# 12       f0  mean   c1_stress:pid5_negative_affectivity_bl_c Interaction  2.000881e+00  1.83207814 -1.617558e+00  5.524989e+00    FALSE
# 13       f0  mean           c1_stress:pid5_psychoticism_bl_c Interaction  4.270257e-01  1.71764962 -2.984029e+00  3.727359e+00    FALSE
# 14       f0  mean           c1_stress:pid5_psychoticism_bl_c Interaction -1.877573e+00  2.00487274 -5.713690e+00  1.995476e+00    FALSE
# 15       f0  mean           c1_stress:pid5_psychoticism_bl_c Interaction -1.254012e+00  1.92786118 -5.025094e+00  2.505021e+00    FALSE
# 16       f0  mean           pid5_antagonism_bl_c:c2_recovery Interaction  3.435216e+00  1.50409233  4.883394e-01  6.422963e+00     TRUE
# 17       f0  mean           pid5_antagonism_bl_c:c2_recovery Interaction  3.305100e+00  1.77765752 -1.986298e-01  6.738806e+00    FALSE
# 18       f0  mean           pid5_antagonism_bl_c:c2_recovery Interaction  4.172447e+00  1.67809935  8.419006e-01  7.441151e+00     TRUE
# 19       f0  mean           pid5_detachment_bl_c:c2_recovery Interaction -8.379729e-01  1.65364964 -4.068156e+00  2.441557e+00    FALSE
# 20       f0  mean           pid5_detachment_bl_c:c2_recovery Interaction -4.365976e+00  2.00236764 -8.297237e+00 -4.556315e-01     TRUE
# 21       f0  mean           pid5_detachment_bl_c:c2_recovery Interaction -1.472385e+00  1.81448302 -5.015600e+00  2.072842e+00    FALSE
# 22       f0  mean        pid5_disinhibition_bl_c:c2_recovery Interaction -9.173861e-01  1.81180481 -4.511608e+00  2.676038e+00    FALSE
# 23       f0  mean        pid5_disinhibition_bl_c:c2_recovery Interaction -7.631532e-01  2.13440849 -4.999948e+00  3.466006e+00    FALSE
# 24       f0  mean        pid5_disinhibition_bl_c:c2_recovery Interaction -6.606141e-01  2.01752211 -4.596291e+00  3.365854e+00    FALSE
# 25       f0  mean pid5_negative_affectivity_bl_c:c2_recovery Interaction -6.359743e-01  1.62048481 -3.794251e+00  2.552436e+00    FALSE
# 26       f0  mean pid5_negative_affectivity_bl_c:c2_recovery Interaction  1.145929e+00  1.90924623 -2.645072e+00  4.917208e+00    FALSE
# 27       f0  mean pid5_negative_affectivity_bl_c:c2_recovery Interaction -3.868473e-01  1.82107352 -4.033826e+00  3.137309e+00    FALSE
# 28       f0  mean         pid5_psychoticism_bl_c:c2_recovery Interaction -1.427037e+00  1.69824236 -4.761669e+00  1.873480e+00    FALSE
# 29       f0  mean         pid5_psychoticism_bl_c:c2_recovery Interaction -9.557453e-01  2.04211114 -4.996326e+00  2.972449e+00    FALSE
# 30       f0  mean         pid5_psychoticism_bl_c:c2_recovery Interaction -2.736859e+00  1.89525427 -6.473222e+00  9.657059e-01    FALSE
# 31       f0  mean                                  Intercept        Main  1.898886e+02  1.81286342  1.863152e+02  1.934339e+02     TRUE
# 32       f0  mean                                  Intercept        Main  1.942461e+02  1.78745057  1.907134e+02  1.977226e+02     TRUE
# 33       f0  mean                                  Intercept        Main  1.949964e+02  1.83614846  1.913622e+02  1.986032e+02     TRUE
# 34       f0  mean                                  c1_stress        Main  3.211738e+00  1.40987030  3.935737e-01  5.954814e+00     TRUE
# 35       f0  mean                                  c1_stress        Main  2.619957e+00  1.67879445 -6.615051e-01  5.939920e+00    FALSE
# 36       f0  mean                                  c1_stress        Main  3.684628e+00  1.61999140  5.032481e-01  6.841354e+00     TRUE
# 37       f0  mean                                c2_recovery        Main  2.617845e-02  1.42466665 -2.776343e+00  2.824029e+00    FALSE
# 38       f0  mean                                c2_recovery        Main -3.130228e+00  1.68006843 -6.385817e+00  1.974881e-01    FALSE
# 39       f0  mean                                c2_recovery        Main  1.537488e-01  1.59783064 -2.951853e+00  3.251020e+00    FALSE
# 40       f0  mean                       pid5_antagonism_bl_c        Main  8.860396e-01  1.92861780 -2.904145e+00  4.571349e+00    FALSE
# 41       f0  mean                       pid5_antagonism_bl_c        Main  1.642594e+00  1.83991566 -1.939310e+00  5.256168e+00    FALSE
# 42       f0  mean                       pid5_antagonism_bl_c        Main  1.215620e+00  1.89635555 -2.530956e+00  4.879516e+00    FALSE
# 43       f0  mean                       pid5_detachment_bl_c        Main  1.333539e-01  2.06801377 -3.934403e+00  4.152322e+00    FALSE
# 44       f0  mean                       pid5_detachment_bl_c        Main  8.446231e-01  1.98552307 -3.007055e+00  4.783629e+00    FALSE
# 45       f0  mean                       pid5_detachment_bl_c        Main -1.281993e-01  2.04828094 -4.152257e+00  3.903979e+00    FALSE
# 46       f0  mean                    pid5_disinhibition_bl_c        Main -1.271109e+00  2.33016283 -5.930822e+00  3.247451e+00    FALSE
# 47       f0  mean                    pid5_disinhibition_bl_c        Main -2.575605e+00  2.21167563 -6.930158e+00  1.872581e+00    FALSE
# 48       f0  mean                    pid5_disinhibition_bl_c        Main -2.039401e+00  2.18350860 -6.397024e+00  2.180870e+00    FALSE
# 49       f0  mean             pid5_negative_affectivity_bl_c        Main -1.985536e+00  2.08641574 -6.037475e+00  2.212636e+00    FALSE
# 50       f0  mean             pid5_negative_affectivity_bl_c        Main -1.555106e+00  2.03054052 -5.566415e+00  2.438296e+00    FALSE
# 51       f0  mean             pid5_negative_affectivity_bl_c        Main -1.778000e+00  1.97553354 -5.605899e+00  2.101692e+00    FALSE
# 52       f0  mean                     pid5_psychoticism_bl_c        Main  6.712931e+00  2.18658591  2.457660e+00  1.106226e+01     TRUE
# 53       f0  mean                     pid5_psychoticism_bl_c        Main  6.873989e+00  2.07624905  2.784480e+00  1.080840e+01     TRUE
# 54       f0  mean                     pid5_psychoticism_bl_c        Main  7.213986e+00  2.09150216  3.156904e+00  1.132339e+01     TRUE
# 55       f0   std             c1_stress:pid5_antagonism_bl_c Interaction  8.076617e-02  0.10549380 -1.263542e-01  2.857114e-01    FALSE
# 56       f0   std             c1_stress:pid5_antagonism_bl_c Interaction  3.149906e-02  0.11390867 -1.903878e-01  2.538599e-01    FALSE
# 57       f0   std             c1_stress:pid5_antagonism_bl_c Interaction  7.496488e-02  0.12154467 -1.626890e-01  3.142556e-01    FALSE
# 58       f0   std             c1_stress:pid5_detachment_bl_c Interaction -4.217637e-02  0.11630115 -2.707145e-01  1.857270e-01    FALSE
# 59       f0   std             c1_stress:pid5_detachment_bl_c Interaction -1.487454e-02  0.12204146 -2.592041e-01  2.210976e-01    FALSE
# 60       f0   std             c1_stress:pid5_detachment_bl_c Interaction -3.943773e-02  0.12710534 -2.847985e-01  2.156320e-01    FALSE
# 61       f0   std          c1_stress:pid5_disinhibition_bl_c Interaction -1.300121e-01  0.12375993 -3.783728e-01  1.164964e-01    FALSE
# 62       f0   std          c1_stress:pid5_disinhibition_bl_c Interaction -1.325320e-01  0.13005869 -3.925580e-01  1.266108e-01    FALSE
# 63       f0   std          c1_stress:pid5_disinhibition_bl_c Interaction  8.953569e-02  0.13537066 -1.751533e-01  3.593409e-01    FALSE
# 64       f0   std   c1_stress:pid5_negative_affectivity_bl_c Interaction  7.357600e-02  0.11419506 -1.514760e-01  2.969743e-01    FALSE
# 65       f0   std   c1_stress:pid5_negative_affectivity_bl_c Interaction  1.585219e-01  0.12173748 -8.243155e-02  3.991897e-01    FALSE
# 66       f0   std   c1_stress:pid5_negative_affectivity_bl_c Interaction  6.026472e-02  0.12510675 -1.874151e-01  3.081773e-01    FALSE
# 67       f0   std           c1_stress:pid5_psychoticism_bl_c Interaction  6.767537e-02  0.11815422 -1.585328e-01  2.977500e-01    FALSE
# 68       f0   std           c1_stress:pid5_psychoticism_bl_c Interaction -8.979592e-02  0.12738006 -3.358076e-01  1.573834e-01    FALSE
# 69       f0   std           c1_stress:pid5_psychoticism_bl_c Interaction -6.823232e-02  0.13327802 -3.303423e-01  1.873280e-01    FALSE
# 70       f0   std           pid5_antagonism_bl_c:c2_recovery Interaction -7.588927e-02  0.10739869 -2.865811e-01  1.292507e-01    FALSE
# 71       f0   std           pid5_antagonism_bl_c:c2_recovery Interaction -4.047882e-02  0.11529177 -2.679538e-01  1.859204e-01    FALSE
# 72       f0   std           pid5_antagonism_bl_c:c2_recovery Interaction  5.307562e-02  0.10542900 -1.537817e-01  2.624158e-01    FALSE
# 73       f0   std           pid5_detachment_bl_c:c2_recovery Interaction -5.974303e-02  0.11626194 -2.869957e-01  1.692939e-01    FALSE
# 74       f0   std           pid5_detachment_bl_c:c2_recovery Interaction -3.761835e-02  0.12294158 -2.750121e-01  2.079222e-01    FALSE
# 75       f0   std           pid5_detachment_bl_c:c2_recovery Interaction -5.024425e-02  0.11128667 -2.671021e-01  1.738401e-01    FALSE
# 76       f0   std        pid5_disinhibition_bl_c:c2_recovery Interaction -2.142601e-02  0.12391395 -2.644730e-01  2.212849e-01    FALSE
# 77       f0   std        pid5_disinhibition_bl_c:c2_recovery Interaction  6.739135e-02  0.13264621 -1.921027e-01  3.294883e-01    FALSE
# 78       f0   std        pid5_disinhibition_bl_c:c2_recovery Interaction -1.099372e-01  0.11996604 -3.452808e-01  1.253133e-01    FALSE
# 79       f0   std pid5_negative_affectivity_bl_c:c2_recovery Interaction -8.847876e-02  0.11537581 -3.161672e-01  1.376970e-01    FALSE
# 80       f0   std pid5_negative_affectivity_bl_c:c2_recovery Interaction  5.301172e-03  0.12080689 -2.297206e-01  2.425691e-01    FALSE
# 81       f0   std pid5_negative_affectivity_bl_c:c2_recovery Interaction  1.385179e-01  0.10981335 -7.909924e-02  3.595915e-01    FALSE
# 82       f0   std         pid5_psychoticism_bl_c:c2_recovery Interaction  1.336824e-01  0.12236348 -1.083452e-01  3.724430e-01    FALSE
# 83       f0   std         pid5_psychoticism_bl_c:c2_recovery Interaction  1.902167e-03  0.12804412 -2.442761e-01  2.541125e-01    FALSE
# 84       f0   std         pid5_psychoticism_bl_c:c2_recovery Interaction  1.182866e-01  0.11561689 -1.093742e-01  3.485819e-01    FALSE
# 85       f0   std                                  Intercept        Main  1.512476e+00  0.04639521  1.421699e+00  1.603883e+00     TRUE
# 86       f0   std                                  Intercept        Main  1.762635e+00  0.05820224  1.648342e+00  1.877625e+00     TRUE
# 87       f0   std                                  Intercept        Main  1.822518e+00  0.06710575  1.693626e+00  1.955016e+00     TRUE
# 88       f0   std                                  c1_stress        Main -1.516473e-01  0.10195368 -3.529591e-01  4.873063e-02    FALSE
# 89       f0   std                                  c1_stress        Main -1.303922e-01  0.11177541 -3.477341e-01  9.372580e-02    FALSE
# 90       f0   std                                  c1_stress        Main -9.596299e-02  0.11853264 -3.328188e-01  1.312097e-01    FALSE
# 91       f0   std                                c2_recovery        Main -1.342866e-01  0.10317605 -3.361481e-01  6.591235e-02    FALSE
# 92       f0   std                                c2_recovery        Main -6.940950e-02  0.11114167 -2.910503e-01  1.479897e-01    FALSE
# 93       f0   std                                c2_recovery        Main -1.465477e-01  0.10259253 -3.478508e-01  5.566283e-02    FALSE
# 94       f0   std                       pid5_antagonism_bl_c        Main  9.525815e-04  0.04787711 -9.330596e-02  9.417538e-02    FALSE
# 95       f0   std                       pid5_antagonism_bl_c        Main -7.060336e-02  0.05978601 -1.900459e-01  4.276532e-02    FALSE
# 96       f0   std                       pid5_antagonism_bl_c        Main -4.946473e-02  0.06952927 -1.878512e-01  8.586350e-02    FALSE
# 97       f0   std                       pid5_detachment_bl_c        Main  6.160272e-02  0.05119323 -3.843009e-02  1.613605e-01    FALSE
# 98       f0   std                       pid5_detachment_bl_c        Main  9.764140e-02  0.06626297 -3.118926e-02  2.307738e-01    FALSE
# 99       f0   std                       pid5_detachment_bl_c        Main  2.630672e-02  0.07369356 -1.189395e-01  1.729608e-01    FALSE
# 100      f0   std                    pid5_disinhibition_bl_c        Main  3.227912e-02  0.05686042 -7.766733e-02  1.450418e-01    FALSE
# 101      f0   std                    pid5_disinhibition_bl_c        Main  7.722565e-03  0.07339299 -1.377658e-01  1.515870e-01    FALSE
# 102      f0   std                    pid5_disinhibition_bl_c        Main  6.222286e-02  0.08637608 -1.034897e-01  2.316705e-01    FALSE
# 103      f0   std             pid5_negative_affectivity_bl_c        Main -5.601089e-02  0.05125838 -1.589504e-01  4.456822e-02    FALSE
# 104      f0   std             pid5_negative_affectivity_bl_c        Main  1.074047e-04  0.06615604 -1.309111e-01  1.299057e-01    FALSE
# 105      f0   std             pid5_negative_affectivity_bl_c        Main -4.731237e-02  0.07346526 -1.939295e-01  9.061103e-02    FALSE
# 106      f0   std                     pid5_psychoticism_bl_c        Main -5.021465e-02  0.05480650 -1.563244e-01  5.793145e-02    FALSE
# 107      f0   std                     pid5_psychoticism_bl_c        Main -8.317651e-03  0.06873705 -1.421752e-01  1.265255e-01    FALSE
# 108      f0   std                     pid5_psychoticism_bl_c        Main -9.733773e-02  0.07743395 -2.504956e-01  5.337643e-02    FALSE
# 109      f2  mean             c1_stress:pid5_antagonism_bl_c Interaction -6.368006e-01 13.41043970 -2.598744e+01  2.653609e+01    FALSE
# 110      f2  mean             c1_stress:pid5_antagonism_bl_c Interaction -4.137619e+01 35.78203125 -1.086165e+02  3.147617e+01    FALSE
# 111      f2  mean             c1_stress:pid5_antagonism_bl_c Interaction  1.133067e+01 11.39136060 -1.021018e+01  3.486975e+01    FALSE
# 112      f2  mean             c1_stress:pid5_detachment_bl_c Interaction -1.645255e+01 13.95103014 -4.351030e+01  1.118241e+01    FALSE
# 113      f2  mean             c1_stress:pid5_detachment_bl_c Interaction -7.245299e+01 41.96997937 -1.545831e+02  1.101561e+01    FALSE
# 114      f2  mean             c1_stress:pid5_detachment_bl_c Interaction -1.478618e+00 13.07329326 -2.743508e+01  2.385009e+01    FALSE
# 115      f2  mean          c1_stress:pid5_disinhibition_bl_c Interaction  7.442064e+00 16.08475592 -2.348242e+01  3.857215e+01    FALSE
# 116      f2  mean          c1_stress:pid5_disinhibition_bl_c Interaction -6.223337e+01 38.14279249 -1.352062e+02  1.278004e+01    FALSE
# 117      f2  mean          c1_stress:pid5_disinhibition_bl_c Interaction -7.066993e+00 13.63628222 -3.319947e+01  1.982108e+01    FALSE
# 118      f2  mean   c1_stress:pid5_negative_affectivity_bl_c Interaction -1.589780e+01 13.32521083 -4.250967e+01  9.855632e+00    FALSE
# 119      f2  mean   c1_stress:pid5_negative_affectivity_bl_c Interaction  1.657200e+01 41.35156509 -6.519042e+01  9.609439e+01    FALSE
# 120      f2  mean   c1_stress:pid5_negative_affectivity_bl_c Interaction  8.417163e+00 11.56496697 -1.402924e+01  3.086333e+01    FALSE
# 121      f2  mean           c1_stress:pid5_psychoticism_bl_c Interaction  2.259753e+01 15.53689859 -7.360707e+00  5.340924e+01    FALSE
# 122      f2  mean           c1_stress:pid5_psychoticism_bl_c Interaction -8.339904e+00 31.85768753 -7.033735e+01  5.435498e+01    FALSE
# 123      f2  mean           c1_stress:pid5_psychoticism_bl_c Interaction -1.053275e+01 14.36599785 -4.003921e+01  1.645328e+01    FALSE
# 124      f2  mean           pid5_antagonism_bl_c:c2_recovery Interaction -2.635467e+00 16.63323714 -3.409736e+01  2.966325e+01    FALSE
# 125      f2  mean           pid5_antagonism_bl_c:c2_recovery Interaction  6.318461e+01 29.71817925  2.453617e+00  1.198689e+02     TRUE
# 126      f2  mean           pid5_antagonism_bl_c:c2_recovery Interaction  6.268118e+00 11.56782106 -1.539416e+01  3.005137e+01    FALSE
# 127      f2  mean           pid5_detachment_bl_c:c2_recovery Interaction -1.546860e+01 14.08663181 -4.328827e+01  1.270248e+01    FALSE
# 128      f2  mean           pid5_detachment_bl_c:c2_recovery Interaction  2.936920e+01 31.01500003 -3.214802e+01  8.964844e+01    FALSE
# 129      f2  mean           pid5_detachment_bl_c:c2_recovery Interaction -2.046368e+01 13.85814070 -4.663958e+01  7.172836e+00    FALSE
# 130      f2  mean        pid5_disinhibition_bl_c:c2_recovery Interaction  1.984811e+01 13.90968622 -7.536217e+00  4.697815e+01    FALSE
# 131      f2  mean        pid5_disinhibition_bl_c:c2_recovery Interaction -1.589317e+01 36.74148700 -8.580746e+01  5.768073e+01    FALSE
# 132      f2  mean        pid5_disinhibition_bl_c:c2_recovery Interaction  1.126745e+01 14.52571942 -1.524944e+01  4.104526e+01    FALSE
# 133      f2  mean pid5_negative_affectivity_bl_c:c2_recovery Interaction  1.150396e-01 12.32877774 -2.327611e+01  2.422225e+01    FALSE
# 134      f2  mean pid5_negative_affectivity_bl_c:c2_recovery Interaction  4.702960e+01 28.61787114 -8.256720e+00  1.033488e+02    FALSE
# 135      f2  mean pid5_negative_affectivity_bl_c:c2_recovery Interaction  1.540669e+01 10.57856517 -6.258826e+00  3.514136e+01    FALSE
# 136      f2  mean         pid5_psychoticism_bl_c:c2_recovery Interaction -7.723505e+00 14.66821579 -3.746821e+01  2.083261e+01    FALSE
# 137      f2  mean         pid5_psychoticism_bl_c:c2_recovery Interaction -1.173792e+02 29.09910488 -1.711910e+02 -5.882386e+01     TRUE
# 138      f2  mean         pid5_psychoticism_bl_c:c2_recovery Interaction -1.050294e+01 14.97974757 -4.183978e+01  1.597560e+01    FALSE
# 139      f2  mean                                  Intercept        Main  1.053858e+03 13.06711670  1.024897e+03  1.077242e+03     TRUE
# 140      f2  mean                                  Intercept        Main  1.408541e+03 11.90186313  1.383902e+03  1.430136e+03     TRUE
# 141      f2  mean                                  Intercept        Main  9.712852e+02  8.57587360  9.556392e+02  9.884473e+02     TRUE
# 142      f2  mean                                  c1_stress        Main -1.694343e+01 11.95738250 -3.931634e+01  7.174596e+00    FALSE
# 143      f2  mean                                  c1_stress        Main -3.992303e+01 31.04796674 -9.968799e+01  2.176479e+01    FALSE
# 144      f2  mean                                  c1_stress        Main -1.166943e+01 10.27120085 -3.071244e+01  8.843516e+00    FALSE
# 145      f2  mean                                c2_recovery        Main  9.189637e-01 11.48860228 -2.249955e+01  2.310300e+01    FALSE
# 146      f2  mean                                c2_recovery        Main -1.301835e+02 24.40120025 -1.765759e+02 -8.168476e+01     TRUE
# 147      f2  mean                                c2_recovery        Main -1.231276e+01 11.89112427 -3.543106e+01  1.148436e+01    FALSE
# 148      f2  mean                       pid5_antagonism_bl_c        Main  8.519361e+00  5.51991949 -1.982733e+00  1.953250e+01    FALSE
# 149      f2  mean                       pid5_antagonism_bl_c        Main -3.821342e+00 14.80898668 -3.415880e+01  2.377932e+01    FALSE
# 150      f2  mean                       pid5_antagonism_bl_c        Main -4.539485e-01  4.60996213 -9.781298e+00  8.369596e+00    FALSE
# 151      f2  mean                       pid5_detachment_bl_c        Main -1.795514e+00  6.00425752 -1.381914e+01  9.542964e+00    FALSE
# 152      f2  mean                       pid5_detachment_bl_c        Main -7.135139e+00 16.10092444 -3.913511e+01  2.313434e+01    FALSE
# 153      f2  mean                       pid5_detachment_bl_c        Main  2.426548e+00  5.75019815 -8.431167e+00  1.395338e+01    FALSE
# 154      f2  mean                    pid5_disinhibition_bl_c        Main  2.293661e+00  5.81616613 -9.331241e+00  1.375442e+01    FALSE
# 155      f2  mean                    pid5_disinhibition_bl_c        Main  6.299064e+01 18.53537025  2.611378e+01  9.885974e+01     TRUE
# 156      f2  mean                    pid5_disinhibition_bl_c        Main  5.718244e+00  5.23339106 -4.505542e+00  1.619462e+01    FALSE
# 157      f2  mean             pid5_negative_affectivity_bl_c        Main -6.779744e+00  4.55702629 -1.573455e+01  2.047044e+00    FALSE
# 158      f2  mean             pid5_negative_affectivity_bl_c        Main -1.945714e+01 16.33900438 -5.135107e+01  1.216891e+01    FALSE
# 159      f2  mean             pid5_negative_affectivity_bl_c        Main -9.193221e+00  3.82240272 -1.692842e+01 -1.644015e+00     TRUE
# 160      f2  mean                     pid5_psychoticism_bl_c        Main -4.721375e+00  6.15248726 -1.694772e+01  6.839710e+00    FALSE
# 161      f2  mean                     pid5_psychoticism_bl_c        Main -7.056496e+01 14.37467715 -9.916890e+01 -4.287246e+01     TRUE
# 162      f2  mean                     pid5_psychoticism_bl_c        Main -1.932322e+00  5.22374290 -1.241144e+01  8.153576e+00    FALSE
# 163      f2   std             c1_stress:pid5_antagonism_bl_c Interaction  1.650895e-03  0.05194209 -1.013145e-01  1.050374e-01    FALSE
# 164      f2   std             c1_stress:pid5_antagonism_bl_c Interaction  1.161080e-01  0.10915730 -1.019964e-01  3.302171e-01    FALSE
# 165      f2   std             c1_stress:pid5_antagonism_bl_c Interaction  4.870020e-02  0.11739203 -1.820542e-01  2.834508e-01    FALSE
# 166      f2   std             c1_stress:pid5_detachment_bl_c Interaction  2.230647e-03  0.05751813 -1.107600e-01  1.149524e-01    FALSE
# 167      f2   std             c1_stress:pid5_detachment_bl_c Interaction -1.479948e-01  0.11878367 -3.777223e-01  8.798004e-02    FALSE
# 168      f2   std             c1_stress:pid5_detachment_bl_c Interaction -7.594856e-02  0.12900102 -3.259567e-01  1.791405e-01    FALSE
# 169      f2   std          c1_stress:pid5_disinhibition_bl_c Interaction -4.064155e-02  0.06350725 -1.642782e-01  8.264055e-02    FALSE
# 170      f2   std          c1_stress:pid5_disinhibition_bl_c Interaction  1.494912e-02  0.12632892 -2.307384e-01  2.625443e-01    FALSE
# 171      f2   std          c1_stress:pid5_disinhibition_bl_c Interaction  1.021137e-01  0.13929907 -1.742720e-01  3.753856e-01    FALSE
# 172      f2   std   c1_stress:pid5_negative_affectivity_bl_c Interaction -2.096494e-01  0.05731336 -3.217830e-01 -9.650887e-02     TRUE
# 173      f2   std   c1_stress:pid5_negative_affectivity_bl_c Interaction  2.843419e-01  0.11782636  5.443838e-02  5.116346e-01     TRUE
# 174      f2   std   c1_stress:pid5_negative_affectivity_bl_c Interaction  2.674435e-02  0.12639941 -2.210490e-01  2.764040e-01    FALSE
# 175      f2   std           c1_stress:pid5_psychoticism_bl_c Interaction  1.029040e-01  0.05914100 -1.206224e-02  2.185343e-01    FALSE
# 176      f2   std           c1_stress:pid5_psychoticism_bl_c Interaction -2.701583e-01  0.12607542 -5.250820e-01 -2.617595e-02     TRUE
# 177      f2   std           c1_stress:pid5_psychoticism_bl_c Interaction -2.479053e-01  0.13357014 -5.112263e-01  1.383848e-02    FALSE
# 178      f2   std           pid5_antagonism_bl_c:c2_recovery Interaction -5.956532e-03  0.05851446 -1.211947e-01  1.081112e-01    FALSE
# 179      f2   std           pid5_antagonism_bl_c:c2_recovery Interaction  8.055534e-02  0.10950760 -1.343421e-01  2.992582e-01    FALSE
# 180      f2   std           pid5_antagonism_bl_c:c2_recovery Interaction -7.699939e-02  0.13591834 -3.433275e-01  1.952523e-01    FALSE
# 181      f2   std           pid5_detachment_bl_c:c2_recovery Interaction -9.721318e-02  0.06410169 -2.233048e-01  2.631737e-02    FALSE
# 182      f2   std           pid5_detachment_bl_c:c2_recovery Interaction -2.122229e-01  0.12117292 -4.463754e-01  2.655067e-02    FALSE
# 183      f2   std           pid5_detachment_bl_c:c2_recovery Interaction  1.415746e-01  0.14612152 -1.506290e-01  4.263749e-01    FALSE
# 184      f2   std        pid5_disinhibition_bl_c:c2_recovery Interaction  1.044750e-01  0.07063223 -3.464280e-02  2.428235e-01    FALSE
# 185      f2   std        pid5_disinhibition_bl_c:c2_recovery Interaction  9.167585e-02  0.12856882 -1.588466e-01  3.483753e-01    FALSE
# 186      f2   std        pid5_disinhibition_bl_c:c2_recovery Interaction -1.040198e-02  0.15563642 -3.208059e-01  2.958683e-01    FALSE
# 187      f2   std pid5_negative_affectivity_bl_c:c2_recovery Interaction -1.278797e-01  0.06419114 -2.555408e-01 -3.998769e-03     TRUE
# 188      f2   std pid5_negative_affectivity_bl_c:c2_recovery Interaction  1.772370e-01  0.11671359 -5.383442e-02  4.052798e-01    FALSE
# 189      f2   std pid5_negative_affectivity_bl_c:c2_recovery Interaction  6.436767e-02  0.14379368 -2.208560e-01  3.479308e-01    FALSE
# 190      f2   std         pid5_psychoticism_bl_c:c2_recovery Interaction  2.890662e-02  0.06843144 -1.052898e-01  1.638020e-01    FALSE
# 191      f2   std         pid5_psychoticism_bl_c:c2_recovery Interaction -1.038312e-01  0.12600511 -3.487518e-01  1.441295e-01    FALSE
# 192      f2   std         pid5_psychoticism_bl_c:c2_recovery Interaction -1.747192e-01  0.15354834 -4.817651e-01  1.245746e-01    FALSE
# 193      f2   std                                  Intercept        Main  3.999850e+00  0.02772694  3.946248e+00  4.053862e+00     TRUE
# 194      f2   std                                  Intercept        Main  4.706050e+00  0.04694456  4.613095e+00  4.799830e+00     TRUE
# 195      f2   std                                  Intercept        Main  4.251209e+00  0.05086596  4.151896e+00  4.353687e+00     TRUE
# 196      f2   std                                  c1_stress        Main -8.145949e-02  0.04863092 -1.761955e-01  1.541260e-02    FALSE
# 197      f2   std                                  c1_stress        Main  1.436008e-01  0.10557661 -6.565959e-02  3.482302e-01    FALSE
# 198      f2   std                                  c1_stress        Main  8.698071e-02  0.11453048 -1.334099e-01  3.114380e-01    FALSE
# 199      f2   std                                c2_recovery        Main -7.601671e-02  0.05495625 -1.849480e-01  3.136028e-02    FALSE
# 200      f2   std                                c2_recovery        Main  9.303037e-02  0.10368433 -1.115410e-01  2.966201e-01    FALSE
# 201      f2   std                                c2_recovery        Main  7.871486e-02  0.13201171 -1.845654e-01  3.413856e-01    FALSE
# 202      f2   std                       pid5_antagonism_bl_c        Main  2.753450e-03  0.02959971 -5.534267e-02  6.194452e-02    FALSE
# 203      f2   std                       pid5_antagonism_bl_c        Main  2.102833e-02  0.04959476 -7.635486e-02  1.181312e-01    FALSE
# 204      f2   std                       pid5_antagonism_bl_c        Main -5.150948e-02  0.05188662 -1.526771e-01  5.199411e-02    FALSE
# 205      f2   std                       pid5_detachment_bl_c        Main  5.563161e-02  0.03240634 -9.473002e-03  1.179686e-01    FALSE
# 206      f2   std                       pid5_detachment_bl_c        Main -1.949615e-02  0.05379538 -1.245596e-01  8.550827e-02    FALSE
# 207      f2   std                       pid5_detachment_bl_c        Main  9.172191e-02  0.05592246 -1.805555e-02  2.030442e-01    FALSE
# 208      f2   std                    pid5_disinhibition_bl_c        Main -2.196333e-02  0.03536245 -9.199075e-02  4.495031e-02    FALSE
# 209      f2   std                    pid5_disinhibition_bl_c        Main -1.254620e-02  0.05880967 -1.282695e-01  1.039333e-01    FALSE
# 210      f2   std                    pid5_disinhibition_bl_c        Main -5.628636e-02  0.06245418 -1.778354e-01  6.498220e-02    FALSE
# 211      f2   std             pid5_negative_affectivity_bl_c        Main  1.797662e-02  0.03207978 -4.620802e-02  7.992583e-02    FALSE
# 212      f2   std             pid5_negative_affectivity_bl_c        Main -2.690954e-02  0.05305925 -1.320197e-01  7.687892e-02    FALSE
# 213      f2   std             pid5_negative_affectivity_bl_c        Main -1.020362e-01  0.05535436 -2.107998e-01  7.013042e-03    FALSE
# 214      f2   std                     pid5_psychoticism_bl_c        Main -7.961758e-03  0.03420969 -7.369189e-02  5.857468e-02    FALSE
# 215      f2   std                     pid5_psychoticism_bl_c        Main  1.141929e-01  0.05616726  5.018629e-03  2.263354e-01     TRUE
# 216      f2   std                     pid5_psychoticism_bl_c        Main  4.121774e-02  0.05951186 -7.455901e-02  1.599363e-01    FALSE
# 217  jitter     a             c1_stress:pid5_antagonism_bl_c Interaction  7.108331e-02  0.08038414 -8.907745e-02  2.252683e-01    FALSE
# 218  jitter     a             c1_stress:pid5_detachment_bl_c Interaction  4.827653e-02  0.08795486 -1.242687e-01  2.212288e-01    FALSE
# 219  jitter     a          c1_stress:pid5_disinhibition_bl_c Interaction -1.273839e-01  0.09608120 -3.135344e-01  6.127762e-02    FALSE
# 220  jitter     a   c1_stress:pid5_negative_affectivity_bl_c Interaction  5.894936e-02  0.08825726 -1.111721e-01  2.309242e-01    FALSE
# 221  jitter     a           c1_stress:pid5_psychoticism_bl_c Interaction  2.620051e-02  0.09101798 -1.541773e-01  2.038844e-01    FALSE
# 222  jitter     a           pid5_antagonism_bl_c:c2_recovery Interaction -6.424005e-02  0.08138759 -2.221849e-01  9.796710e-02    FALSE
# 223  jitter     a           pid5_detachment_bl_c:c2_recovery Interaction -3.084417e-02  0.08828822 -2.057247e-01  1.393446e-01    FALSE
# 224  jitter     a        pid5_disinhibition_bl_c:c2_recovery Interaction  2.021429e-02  0.09620543 -1.713577e-01  2.096604e-01    FALSE
# 225  jitter     a pid5_negative_affectivity_bl_c:c2_recovery Interaction -1.057480e-01  0.08716523 -2.720065e-01  6.657547e-02    FALSE
# 226  jitter     a         pid5_psychoticism_bl_c:c2_recovery Interaction  1.016141e-01  0.09121669 -7.818922e-02  2.825419e-01    FALSE
# 227  jitter     a                                  Intercept        Main -4.455586e-01  0.04426741 -5.317236e-01 -3.593284e-01     TRUE
# 228  jitter     a                                  c1_stress        Main -1.724069e-01  0.07853224 -3.218645e-01 -1.629371e-02     TRUE
# 229  jitter     a                                c2_recovery        Main -1.005783e-01  0.07785258 -2.519338e-01  4.904537e-02    FALSE
# 230  jitter     a                       pid5_antagonism_bl_c        Main  7.265148e-03  0.04704202 -8.694956e-02  9.889274e-02    FALSE
# 231  jitter     a                       pid5_detachment_bl_c        Main  6.571251e-02  0.05102318 -3.360957e-02  1.675951e-01    FALSE
# 232  jitter     a                    pid5_disinhibition_bl_c        Main  1.030740e-03  0.05565132 -1.069423e-01  1.090629e-01    FALSE
# 233  jitter     a             pid5_negative_affectivity_bl_c        Main -3.747156e-02  0.04929024 -1.336702e-01  6.236363e-02    FALSE
# 234  jitter     a                     pid5_psychoticism_bl_c        Main -8.251522e-02  0.05410268 -1.886206e-01  2.473022e-02    FALSE
# 235  jitter     i             c1_stress:pid5_antagonism_bl_c Interaction  3.793794e-02  0.10639760 -1.734801e-01  2.494901e-01    FALSE
# 236  jitter     i             c1_stress:pid5_detachment_bl_c Interaction  1.060888e-01  0.11387238 -1.170516e-01  3.299477e-01    FALSE
# 237  jitter     i          c1_stress:pid5_disinhibition_bl_c Interaction -1.187150e-01  0.12349982 -3.626516e-01  1.240125e-01    FALSE
# 238  jitter     i   c1_stress:pid5_negative_affectivity_bl_c Interaction  5.828255e-02  0.11359561 -1.648302e-01  2.769237e-01    FALSE
# 239  jitter     i           c1_stress:pid5_psychoticism_bl_c Interaction -1.562110e-01  0.12001560 -3.874392e-01  7.932256e-02    FALSE
# 240  jitter     i           pid5_antagonism_bl_c:c2_recovery Interaction  2.614399e-02  0.09723831 -1.660608e-01  2.168701e-01    FALSE
# 241  jitter     i           pid5_detachment_bl_c:c2_recovery Interaction  5.872888e-02  0.10296457 -1.457123e-01  2.578286e-01    FALSE
# 242  jitter     i        pid5_disinhibition_bl_c:c2_recovery Interaction  9.785055e-03  0.11094214 -2.037670e-01  2.268496e-01    FALSE
# 243  jitter     i pid5_negative_affectivity_bl_c:c2_recovery Interaction  3.416306e-02  0.10205809 -1.708981e-01  2.319369e-01    FALSE
# 244  jitter     i         pid5_psychoticism_bl_c:c2_recovery Interaction -5.757892e-02  0.10759807 -2.653640e-01  1.575481e-01    FALSE
# 245  jitter     i                                  Intercept        Main -3.475328e-01  0.05377292 -4.554308e-01 -2.408323e-01     TRUE
# 246  jitter     i                                  c1_stress        Main -1.420508e-01  0.10497614 -3.436256e-01  7.118672e-02    FALSE
# 247  jitter     i                                c2_recovery        Main -7.344440e-02  0.09468919 -2.588136e-01  1.118276e-01    FALSE
# 248  jitter     i                       pid5_antagonism_bl_c        Main -5.816019e-02  0.05493120 -1.684563e-01  4.734910e-02    FALSE
# 249  jitter     i                       pid5_detachment_bl_c        Main  1.101245e-01  0.06051880 -1.019753e-02  2.282159e-01    FALSE
# 250  jitter     i                    pid5_disinhibition_bl_c        Main  2.785394e-02  0.06677582 -1.038540e-01  1.643875e-01    FALSE
# 251  jitter     i             pid5_negative_affectivity_bl_c        Main -3.130140e-02  0.05976760 -1.503725e-01  8.213375e-02    FALSE
# 252  jitter     i                     pid5_psychoticism_bl_c        Main -3.560203e-02  0.06404164 -1.652429e-01  8.579853e-02    FALSE
# 253  jitter     u             c1_stress:pid5_antagonism_bl_c Interaction  1.119663e-01  0.11732839 -1.182148e-01  3.446699e-01    FALSE
# 254  jitter     u             c1_stress:pid5_detachment_bl_c Interaction -2.215273e-02  0.12582865 -2.724181e-01  2.240818e-01    FALSE
# 255  jitter     u          c1_stress:pid5_disinhibition_bl_c Interaction  4.732995e-02  0.13494387 -2.183455e-01  3.065850e-01    FALSE
# 256  jitter     u   c1_stress:pid5_negative_affectivity_bl_c Interaction  1.343165e-01  0.12282499 -1.091155e-01  3.748340e-01    FALSE
# 257  jitter     u           c1_stress:pid5_psychoticism_bl_c Interaction -1.752622e-01  0.12874848 -4.254449e-01  7.964234e-02    FALSE
# 258  jitter     u           pid5_antagonism_bl_c:c2_recovery Interaction  5.250207e-02  0.08934300 -1.258983e-01  2.307380e-01    FALSE
# 259  jitter     u           pid5_detachment_bl_c:c2_recovery Interaction -4.707672e-02  0.09573845 -2.359154e-01  1.401717e-01    FALSE
# 260  jitter     u        pid5_disinhibition_bl_c:c2_recovery Interaction -1.126410e-01  0.10492962 -3.148242e-01  9.320347e-02    FALSE
# 261  jitter     u pid5_negative_affectivity_bl_c:c2_recovery Interaction  1.520489e-01  0.09602983 -3.701214e-02  3.410650e-01    FALSE
# 262  jitter     u         pid5_psychoticism_bl_c:c2_recovery Interaction  4.694004e-02  0.10129335 -1.516889e-01  2.469365e-01    FALSE
# 263  jitter     u                                  Intercept        Main -2.857360e-01  0.06399234 -4.084288e-01 -1.587654e-01     TRUE
# 264  jitter     u                                  c1_stress        Main -1.227199e-01  0.11391685 -3.479232e-01  9.897548e-02    FALSE
# 265  jitter     u                                c2_recovery        Main -1.234193e-01  0.08640645 -2.932768e-01  4.484114e-02    FALSE
# 266  jitter     u                       pid5_antagonism_bl_c        Main -3.197404e-02  0.06455193 -1.584286e-01  9.276564e-02    FALSE
# 267  jitter     u                       pid5_detachment_bl_c        Main  6.209557e-02  0.07237866 -7.893968e-02  2.045238e-01    FALSE
# 268  jitter     u                    pid5_disinhibition_bl_c        Main  4.323133e-02  0.07590482 -1.032994e-01  1.973288e-01    FALSE
# 269  jitter     u             pid5_negative_affectivity_bl_c        Main -1.162723e-02  0.06941949 -1.481872e-01  1.226861e-01    FALSE
# 270  jitter     u                     pid5_psychoticism_bl_c        Main -1.315692e-01  0.07633800 -2.848609e-01  1.465873e-02    FALSE
# 271     nne     a             c1_stress:pid5_antagonism_bl_c Interaction -4.068014e-01  0.48012813 -1.359189e+00  5.341058e-01    FALSE
# 272     nne     a             c1_stress:pid5_detachment_bl_c Interaction  2.905842e-01  0.51589285 -7.219582e-01  1.313250e+00    FALSE
# 273     nne     a          c1_stress:pid5_disinhibition_bl_c Interaction  9.129720e-02  0.56792391 -1.016519e+00  1.197250e+00    FALSE
# 274     nne     a   c1_stress:pid5_negative_affectivity_bl_c Interaction -4.197508e-02  0.51370052 -1.049974e+00  9.439225e-01    FALSE
# 275     nne     a           c1_stress:pid5_psychoticism_bl_c Interaction  3.403435e-01  0.53250429 -7.037047e-01  1.388744e+00    FALSE
# 276     nne     a           pid5_antagonism_bl_c:c2_recovery Interaction -5.848014e-01  0.46870888 -1.511187e+00  3.433352e-01    FALSE
# 277     nne     a           pid5_detachment_bl_c:c2_recovery Interaction  1.704537e-01  0.49976631 -8.258911e-01  1.157091e+00    FALSE
# 278     nne     a        pid5_disinhibition_bl_c:c2_recovery Interaction  6.415194e-02  0.55082067 -1.002104e+00  1.116285e+00    FALSE
# 279     nne     a pid5_negative_affectivity_bl_c:c2_recovery Interaction -2.472533e-01  0.50220963 -1.242850e+00  7.291988e-01    FALSE
# 280     nne     a         pid5_psychoticism_bl_c:c2_recovery Interaction  5.132109e-01  0.52763053 -5.144537e-01  1.539359e+00    FALSE
# 281     nne     a                                  Intercept        Main -2.436671e+01  0.29883814 -2.495840e+01 -2.379161e+01     TRUE
# 282     nne     a                                  c1_stress        Main -1.263863e+00  0.44231964 -2.148526e+00 -3.816885e-01     TRUE
# 283     nne     a                                c2_recovery        Main -5.555404e-01  0.43856644 -1.407629e+00  3.273418e-01    FALSE
# 284     nne     a                       pid5_antagonism_bl_c        Main -1.522165e-01  0.30899590 -7.636353e-01  4.521150e-01    FALSE
# 285     nne     a                       pid5_detachment_bl_c        Main -1.661763e-01  0.34321243 -8.399683e-01  5.068904e-01    FALSE
# 286     nne     a                    pid5_disinhibition_bl_c        Main  3.996403e-01  0.37568948 -3.414221e-01  1.134557e+00    FALSE
# 287     nne     a             pid5_negative_affectivity_bl_c        Main  3.538848e-01  0.34059723 -3.247348e-01  1.012177e+00    FALSE
# 288     nne     a                     pid5_psychoticism_bl_c        Main -7.934598e-01  0.35455452 -1.497900e+00 -9.377364e-02     TRUE
# 289     nne     i             c1_stress:pid5_antagonism_bl_c Interaction  4.844081e-02  0.40885550 -7.433355e-01  8.431950e-01    FALSE
# 290     nne     i             c1_stress:pid5_detachment_bl_c Interaction  2.370561e-01  0.44052451 -6.238408e-01  1.092450e+00    FALSE
# 291     nne     i          c1_stress:pid5_disinhibition_bl_c Interaction -1.488630e-01  0.49653607 -1.118962e+00  8.099853e-01    FALSE
# 292     nne     i   c1_stress:pid5_negative_affectivity_bl_c Interaction  1.140988e-02  0.44395820 -8.758261e-01  8.705917e-01    FALSE
# 293     nne     i           c1_stress:pid5_psychoticism_bl_c Interaction -7.824342e-02  0.46725091 -9.960636e-01  8.270054e-01    FALSE
# 294     nne     i           pid5_antagonism_bl_c:c2_recovery Interaction -1.573365e-01  0.43381213 -9.841818e-01  7.051827e-01    FALSE
# 295     nne     i           pid5_detachment_bl_c:c2_recovery Interaction -2.478633e-02  0.46858824 -9.344615e-01  8.898950e-01    FALSE
# 296     nne     i        pid5_disinhibition_bl_c:c2_recovery Interaction -5.260341e-02  0.52177834 -1.091849e+00  9.714626e-01    FALSE
# 297     nne     i pid5_negative_affectivity_bl_c:c2_recovery Interaction -2.721817e-02  0.47539057 -9.600925e-01  8.870038e-01    FALSE
# 298     nne     i         pid5_psychoticism_bl_c:c2_recovery Interaction  3.909753e-01  0.49040722 -5.903001e-01  1.362933e+00    FALSE
# 299     nne     i                                  Intercept        Main -2.787844e+01  0.25346301 -2.838254e+01 -2.738132e+01     TRUE
# 300     nne     i                                  c1_stress        Main -3.300022e-01  0.38063885 -1.075745e+00  4.284524e-01    FALSE
# 301     nne     i                                c2_recovery        Main  3.823585e-01  0.39631188 -4.024278e-01  1.148870e+00    FALSE
# 302     nne     i                       pid5_antagonism_bl_c        Main  1.073587e-02  0.27493256 -5.444321e-01  5.363737e-01    FALSE
# 303     nne     i                       pid5_detachment_bl_c        Main  2.230644e-01  0.30268438 -3.547518e-01  8.154063e-01    FALSE
# 304     nne     i                    pid5_disinhibition_bl_c        Main -3.989422e-02  0.31668239 -6.555632e-01  5.948967e-01    FALSE
# 305     nne     i             pid5_negative_affectivity_bl_c        Main  2.537861e-01  0.29022345 -3.247584e-01  8.176235e-01    FALSE
# 306     nne     i                     pid5_psychoticism_bl_c        Main -6.206394e-01  0.30804408 -1.219488e+00 -1.388641e-02     TRUE
# 307     nne     u             c1_stress:pid5_antagonism_bl_c Interaction  1.382761e-01  0.42190037 -6.834820e-01  9.818131e-01    FALSE
# 308     nne     u             c1_stress:pid5_detachment_bl_c Interaction  1.924514e-01  0.46123487 -7.298015e-01  1.092210e+00    FALSE
# 309     nne     u          c1_stress:pid5_disinhibition_bl_c Interaction  2.235015e-01  0.50112135 -7.633798e-01  1.218567e+00    FALSE
# 310     nne     u   c1_stress:pid5_negative_affectivity_bl_c Interaction -3.851728e-01  0.45491437 -1.265414e+00  5.308303e-01    FALSE
# 311     nne     u           c1_stress:pid5_psychoticism_bl_c Interaction -3.609780e-01  0.47920084 -1.300999e+00  5.817405e-01    FALSE
# 312     nne     u           pid5_antagonism_bl_c:c2_recovery Interaction -2.705601e-01  0.39402647 -1.025887e+00  4.984782e-01    FALSE
# 313     nne     u           pid5_detachment_bl_c:c2_recovery Interaction  3.343792e-01  0.42720251 -4.858700e-01  1.166485e+00    FALSE
# 314     nne     u        pid5_disinhibition_bl_c:c2_recovery Interaction -1.094524e-01  0.46314313 -1.005351e+00  7.866441e-01    FALSE
# 315     nne     u pid5_negative_affectivity_bl_c:c2_recovery Interaction -2.578975e-01  0.42376841 -1.079011e+00  5.698824e-01    FALSE
# 316     nne     u         pid5_psychoticism_bl_c:c2_recovery Interaction  3.458444e-01  0.44991307 -5.292957e-01  1.219891e+00    FALSE
# 317     nne     u                                  Intercept        Main -2.841481e+01  0.25260464 -2.890133e+01 -2.791552e+01     TRUE
# 318     nne     u                                  c1_stress        Main -9.361528e-01  0.40720492 -1.733858e+00 -1.198827e-01     TRUE
# 319     nne     u                                c2_recovery        Main -2.513920e-01  0.37278605 -9.717288e-01  4.848716e-01    FALSE
# 320     nne     u                       pid5_antagonism_bl_c        Main -1.179981e-01  0.26431061 -6.360230e-01  4.134542e-01    FALSE
# 321     nne     u                       pid5_detachment_bl_c        Main  3.316347e-01  0.29009503 -2.360970e-01  9.096477e-01    FALSE
# 322     nne     u                    pid5_disinhibition_bl_c        Main  2.506848e-01  0.31921882 -3.538397e-01  8.828269e-01    FALSE
# 323     nne     u             pid5_negative_affectivity_bl_c        Main  1.984905e-02  0.28097959 -5.270628e-01  5.793302e-01    FALSE
# 324     nne     u                     pid5_psychoticism_bl_c        Main -6.071220e-01  0.30033709 -1.195607e+00 -1.699462e-02     TRUE
