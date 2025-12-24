# ==============================================================================
# MODELLI ALTERNATIVI: Come può Disinhibition influenzare la voce?
# ==============================================================================
# Esploriamo approcci diversi dalla moderazione lineare standard
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(bayestestR)
  library(tidybayes)
})

options(brms.backend = "cmdstanr")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("MODELLI ALTERNATIVI PER DISINHIBITION\n")
cat(rep("=", 70), "\n\n", sep = "")

# ==============================================================================
# 1. CARICA DATI
# ==============================================================================

df <- readRDS("results/df_analysis.rds")

# Assicurati che i contrasti siano presenti
df <- df %>%
  mutate(
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre" ~ 0.5,
      timepoint == "post" ~ 0
    ),
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre" ~ -0.5,
      timepoint == "post" ~ 0.5
    )
  )

cat("N osservazioni:", nrow(df), "\n")
cat("N soggetti:", n_distinct(df$ID), "\n\n")

# ==============================================================================
# APPROCCIO 1: VARIABILITÀ INTRA-INDIVIDUALE
# ==============================================================================
# Teoria: Disinhibition = instabilità, disregolazione
# Ipotesi: Chi ha alta Disinhibition mostra MAGGIORE VARIABILITÀ nella risposta
#          vocale attraverso i timepoint (non un livello medio diverso)

cat(rep("=", 70), "\n", sep = "")
cat("APPROCCIO 1: Disinhibition predice VARIABILITÀ intra-individuale?\n")
cat(rep("=", 70), "\n\n", sep = "")

# Calcola variabilità intra-soggetto per ogni outcome
variability_df <- df %>%
  group_by(ID) %>%
  summarise(
    # SD intra-individuale per ogni outcome
    sd_f0_mean_a = sd(f0_mean_a, na.rm = TRUE),
    sd_f0_std_a = sd(f0_std_a, na.rm = TRUE),
    sd_jitter_a = sd(jitter_a, na.rm = TRUE),
    sd_nne_a = sd(nne_a, na.rm = TRUE),
    sd_f2_mean_a = sd(f2_mean_a, na.rm = TRUE),
    sd_f2_std_a = sd(f2_std_a, na.rm = TRUE),

    # Range (max - min) come misura alternativa
    range_f0_mean_a = diff(range(f0_mean_a, na.rm = TRUE)),
    range_f2_mean_a = diff(range(f2_mean_a, na.rm = TRUE)),

    # Coefficiente di variazione (CV = SD/Mean)
    cv_f0_mean_a = sd(f0_mean_a, na.rm = TRUE) / mean(f0_mean_a, na.rm = TRUE),
    cv_f2_mean_a = sd(f2_mean_a, na.rm = TRUE) / mean(f2_mean_a, na.rm = TRUE),

    # PID-5 (prendi il primo valore - sono costanti per soggetto)
    pid5_dis = first(pid5_disinhibition),
    pid5_na = first(pid5_negative_affectivity),
    pid5_det = first(pid5_detachment),
    pid5_ant = first(pid5_antagonism),
    pid5_psy = first(pid5_psychoticism),

    .groups = "drop"
  ) %>%
  mutate(
    pid5_dis_c = scale(pid5_dis)[, 1],
    pid5_na_c = scale(pid5_na)[, 1],
    pid5_det_c = scale(pid5_det)[, 1],
    pid5_ant_c = scale(pid5_ant)[, 1],
    pid5_psy_c = scale(pid5_psy)[, 1]
  )

cat("Correlazioni Disinhibition con variabilità intra-individuale:\n\n")

var_outcomes <- c(
  "sd_f0_mean_a",
  "sd_f0_std_a",
  "sd_jitter_a",
  "sd_nne_a",
  "sd_f2_mean_a",
  "sd_f2_std_a",
  "range_f0_mean_a",
  "range_f2_mean_a",
  "cv_f0_mean_a",
  "cv_f2_mean_a"
)

cor_results <- sapply(var_outcomes, function(v) {
  cor(variability_df$pid5_dis, variability_df[[v]], use = "complete.obs")
})

cor_df <- tibble(
  outcome = var_outcomes,
  r_disinhibition = round(cor_results, 3)
) %>%
  arrange(desc(abs(r_disinhibition)))

print(cor_df)

# Test formale con modello bayesiano
cat("\n--- Modello: SD(F0 mean) ~ Disinhibition + altri tratti ---\n")

# Rimuovi NA e valori infiniti
variability_clean <- variability_df %>%
  filter(is.finite(sd_f0_mean_a), sd_f0_mean_a > 0)

m_var_f0 <- brm(
  sd_f0_mean_a ~ pid5_dis_c + pid5_na_c + pid5_det_c + pid5_ant_c + pid5_psy_c,
  data = variability_clean,
  family = lognormal(), # SD è sempre positivo
  prior = c(
    prior(normal(0, 1), class = b),
    prior(student_t(3, 2, 1), class = Intercept)
  ),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  file = "models/m_variability_f0"
)

cat("\nRisultati:\n")
print(summary(m_var_f0)$fixed)

cat("\n--- Modello: SD(F2 mean) ~ Disinhibition + altri tratti ---\n")

variability_clean2 <- variability_df %>%
  filter(is.finite(sd_f2_mean_a), sd_f2_mean_a > 0)

m_var_f2 <- brm(
  sd_f2_mean_a ~ pid5_dis_c + pid5_na_c + pid5_det_c + pid5_ant_c + pid5_psy_c,
  data = variability_clean2,
  family = lognormal(),
  prior = c(
    prior(normal(0, 1), class = b),
    prior(student_t(3, 4, 1), class = Intercept)
  ),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  file = "models/m_variability_f2"
)

cat("\nRisultati:\n")
print(summary(m_var_f2)$fixed)

# ==============================================================================
# APPROCCIO 2: EFFETTI NON-LINEARI (SOGLIA / QUADRATICI)
# ==============================================================================
# Teoria: Forse solo chi ha Disinhibition MOLTO ALTA mostra effetti
# Ipotesi: Effetto quadratico o effetto soglia

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("APPROCCIO 2: Effetti NON-LINEARI di Disinhibition?\n")
cat(rep("=", 70), "\n\n", sep = "")

# Crea termine quadratico e gruppi
df <- df %>%
  mutate(
    pid5_dis_c = scale(pid5_disinhibition)[, 1],
    pid5_dis_c2 = pid5_dis_c^2, # termine quadratico
    pid5_dis_high = ifelse(pid5_dis_c > 1, 1, 0), # sopra +1 SD
    pid5_dis_group = case_when(
      pid5_dis_c < -0.5 ~ "low",
      pid5_dis_c > 0.5 ~ "high",
      TRUE ~ "medium"
    )
  )

cat("Distribuzione gruppi Disinhibition:\n")
print(table(df$pid5_dis_group) / 3) # /3 perché ogni soggetto ha 3 obs

# Modello con termine quadratico
cat("\n--- Modello F0 mean con effetto QUADRATICO di Disinhibition ---\n")

m_quad <- brm(
  f0_mean_a ~
    c1_stress *
      (pid5_dis_c + pid5_dis_c2) +
      c2_recovery * (pid5_dis_c + pid5_dis_c2) +
      (1 + c1_stress + c2_recovery || ID),
  data = df,
  family = gaussian(),
  prior = c(
    prior(student_t(3, 220, 30), class = Intercept),
    prior(normal(0, 10), class = b),
    prior(exponential(0.3), class = sigma),
    prior(exponential(0.3), class = sd)
  ),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95),
  file = "models/m_f0_disinhibition_quadratic"
)

cat("\nEffetti fissi (focus su Disinhibition):\n")
fe_quad <- fixef(m_quad)
dis_rows <- grep("dis", rownames(fe_quad))
print(round(fe_quad[dis_rows, ], 3))

# Modello con F2 (l'outcome dell'abstract)
cat("\n--- Modello F2 mean con effetto QUADRATICO di Disinhibition ---\n")

m_quad_f2 <- brm(
  f2_mean_a ~
    c1_stress *
      (pid5_dis_c + pid5_dis_c2) +
      c2_recovery * (pid5_dis_c + pid5_dis_c2) +
      (1 + c1_stress + c2_recovery || ID),
  data = df,
  family = student(),
  prior = c(
    prior(student_t(3, 1200, 150), class = Intercept),
    prior(normal(0, 50), class = b),
    prior(exponential(0.01), class = sigma),
    prior(exponential(0.01), class = sd),
    prior(gamma(2, 0.1), class = nu)
  ),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95),
  file = "models/m_f2_disinhibition_quadratic"
)

cat("\nEffetti fissi (focus su Disinhibition):\n")
fe_quad_f2 <- fixef(m_quad_f2)
dis_rows_f2 <- grep("dis", rownames(fe_quad_f2))
print(round(fe_quad_f2[dis_rows_f2, ], 3))

# ==============================================================================
# APPROCCIO 3: INTERAZIONI TRATTO × TRATTO
# ==============================================================================
# Teoria: Disinhibition potrebbe MODERARE l'effetto di altri tratti
# Es: NA causa stress vocale, ma solo se anche Disinhibition è alta

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("APPROCCIO 3: Disinhibition MODERA l'effetto di altri tratti?\n")
cat(rep("=", 70), "\n\n", sep = "")

# Prepara variabili
df <- df %>%
  mutate(
    pid5_na_c = scale(pid5_negative_affectivity)[, 1],
    pid5_psy_c = scale(pid5_psychoticism)[, 1]
  )

# NA × Disinhibition interaction
cat("--- Modello: F0 mean ~ Stress × NA × Disinhibition ---\n")

m_3way_na <- brm(
  f0_mean_a ~
    c1_stress *
      pid5_na_c *
      pid5_dis_c +
      c2_recovery * pid5_na_c * pid5_dis_c +
      (1 + c1_stress + c2_recovery || ID),
  data = df,
  family = gaussian(),
  prior = c(
    prior(student_t(3, 220, 30), class = Intercept),
    prior(normal(0, 10), class = b),
    prior(exponential(0.3), class = sigma),
    prior(exponential(0.3), class = sd)
  ),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95),
  file = "models/m_f0_NA_x_Dis_interaction"
)

cat("\nEffetti fissi (focus su interazioni con Disinhibition):\n")
fe_3way <- fixef(m_3way_na)
dis_rows_3way <- grep("dis", rownames(fe_3way), ignore.case = TRUE)
print(round(fe_3way[dis_rows_3way, ], 3))

# Test specifico per l'interazione a 3 vie
cat("\n--- Hypothesis test: Stress × NA × Disinhibition ---\n")
h_3way <- hypothesis(m_3way_na, "c1_stress:pid5_na_c:pid5_dis_c = 0")
print(h_3way)

# ==============================================================================
# APPROCCIO 4: COMPOSITE OUTCOME - Reattività vocale globale
# ==============================================================================
# Idea: Invece di singoli outcome, creare un indice composito

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("APPROCCIO 4: Indice COMPOSITO di reattività vocale\n")
cat(rep("=", 70), "\n\n", sep = "")

# Calcola il cambiamento pre-baseline per ogni outcome (stress reactivity)
reactivity_df <- df %>%
  select(
    ID,
    timepoint,
    f0_mean_a,
    f0_std_a,
    jitter_a,
    nne_a,
    f2_mean_a,
    f2_std_a,
    pid5_disinhibition,
    pid5_negative_affectivity,
    pid5_detachment,
    pid5_antagonism,
    pid5_psychoticism
  ) %>%
  pivot_wider(
    id_cols = c(ID, starts_with("pid5")),
    names_from = timepoint,
    values_from = c(f0_mean_a, f0_std_a, jitter_a, nne_a, f2_mean_a, f2_std_a)
  ) %>%
  mutate(
    # Stress reactivity (pre - baseline), standardizzato
    react_f0_mean = scale(f0_mean_a_pre - f0_mean_a_baseline)[, 1],
    react_f0_std = scale(f0_std_a_pre - f0_std_a_baseline)[, 1],
    react_jitter = scale(jitter_a_pre - jitter_a_baseline)[, 1],
    react_nne = scale(nne_a_pre - nne_a_baseline)[, 1],
    react_f2_mean = scale(f2_mean_a_pre - f2_mean_a_baseline)[, 1],
    react_f2_std = scale(f2_std_a_pre - f2_std_a_baseline)[, 1],

    # Composite: media dei cambiamenti standardizzati (in valore assoluto)
    # Cattura "quanto cambia la voce" indipendentemente dalla direzione
    composite_reactivity_abs = (abs(react_f0_mean) +
      abs(react_f0_std) +
      abs(react_jitter) +
      abs(react_nne) +
      abs(react_f2_mean) +
      abs(react_f2_std)) /
      6,

    # Composite direzionale (aumenti = positivo)
    composite_reactivity_dir = (react_f0_mean +
      react_f0_std +
      react_jitter +
      react_nne) /
      4,

    # Standardizza predittori
    pid5_dis_c = scale(pid5_disinhibition)[, 1],
    pid5_na_c = scale(pid5_negative_affectivity)[, 1],
    pid5_det_c = scale(pid5_detachment)[, 1],
    pid5_ant_c = scale(pid5_antagonism)[, 1],
    pid5_psy_c = scale(pid5_psychoticism)[, 1]
  )

cat("Correlazioni con composite reactivity (absolute):\n")
cor_composite <- c(
  Disinhibition = cor(
    reactivity_df$pid5_dis_c,
    reactivity_df$composite_reactivity_abs,
    use = "complete.obs"
  ),
  `NA` = cor(
    # Modifica qui: NA tra backticks
    reactivity_df$pid5_na_c,
    reactivity_df$composite_reactivity_abs,
    use = "complete.obs"
  ),
  Detachment = cor(
    reactivity_df$pid5_det_c,
    reactivity_df$composite_reactivity_abs,
    use = "complete.obs"
  ),
  Antagonism = cor(
    reactivity_df$pid5_ant_c,
    reactivity_df$composite_reactivity_abs,
    use = "complete.obs"
  ),
  Psychoticism = cor(
    reactivity_df$pid5_psy_c,
    reactivity_df$composite_reactivity_abs,
    use = "complete.obs"
  )
)
print(round(sort(cor_composite, decreasing = TRUE), 3))

# Modello formale
cat("\n--- Modello: Composite Reactivity ~ tutti i tratti PID-5 ---\n")

reactivity_clean <- reactivity_df %>%
  filter(is.finite(composite_reactivity_abs))

m_composite <- brm(
  composite_reactivity_abs ~
    pid5_dis_c + pid5_na_c + pid5_det_c + pid5_ant_c + pid5_psy_c,
  data = reactivity_clean,
  family = lognormal(), # sempre positivo
  prior = c(
    prior(normal(0, 1), class = b),
    prior(student_t(3, 0, 1), class = Intercept)
  ),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  file = "models/m_composite_reactivity"
)

cat("\nRisultati:\n")
print(summary(m_composite)$fixed)

# ==============================================================================
# APPROCCIO 5: CONFRONTO HIGH vs LOW Disinhibition
# ==============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("APPROCCIO 5: Confronto HIGH vs LOW Disinhibition\n")
cat(rep("=", 70), "\n\n", sep = "")

# Dividi in terzili
df <- df %>%
  mutate(
    dis_tertile = ntile(pid5_disinhibition, 3),
    dis_group = factor(dis_tertile, labels = c("Low", "Medium", "High"))
  )

# Statistiche descrittive per gruppo
cat("Reattività F0 per gruppo Disinhibition:\n\n")

f0_by_group <- df %>%
  group_by(dis_group, timepoint) %>%
  summarise(
    M_f0 = mean(f0_mean_a, na.rm = TRUE),
    SD_f0 = sd(f0_mean_a, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = timepoint, values_from = c(M_f0, SD_f0))

print(f0_by_group)

# Calcola stress reactivity per gruppo
cat("\n\nDelta F0 (pre - baseline) per gruppo:\n")

delta_by_group <- df %>%
  select(ID, dis_group, timepoint, f0_mean_a, f2_mean_a) %>%
  pivot_wider(names_from = timepoint, values_from = c(f0_mean_a, f2_mean_a)) %>%
  mutate(
    delta_f0 = f0_mean_a_pre - f0_mean_a_baseline,
    delta_f2 = f2_mean_a_pre - f2_mean_a_baseline
  ) %>%
  group_by(dis_group) %>%
  summarise(
    M_delta_f0 = mean(delta_f0, na.rm = TRUE),
    SD_delta_f0 = sd(delta_f0, na.rm = TRUE),
    M_delta_f2 = mean(delta_f2, na.rm = TRUE),
    SD_delta_f2 = sd(delta_f2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(delta_by_group)

# ==============================================================================
# RIEPILOGO FINALE
# ==============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("RIEPILOGO: Quali approcci mostrano effetti di Disinhibition?\n")
cat(rep("=", 70), "\n", sep = "")

cat(
  "
APPROCCIO 1 - Variabilità intra-individuale:
→ Disinhibition predice maggiore INSTABILITÀ nella voce attraverso i timepoint?

APPROCCIO 2 - Effetti non-lineari:
→ Solo chi ha Disinhibition ALTA mostra effetti (soglia/quadratico)?

APPROCCIO 3 - Interazione tratto × tratto:
→ Disinhibition MODERA l'effetto di NA o Psychoticism?

APPROCCIO 4 - Composite outcome:
→ Disinhibition predice reattività vocale GLOBALE?

APPROCCIO 5 - Confronto gruppi:
→ Pattern diversi tra LOW/MEDIUM/HIGH Disinhibition?

Vedi i risultati sopra per determinare quali approcci sono promettenti.

Se NESSUN approccio funziona:
- Disinhibition potrebbe semplicemente non influenzare la voce in contesti di stress da esame
- Potrebbe richiedere contesti interpersonali/impulsivi per emergere
- L'effetto originale potrebbe essere stato un falso positivo
"
)

# Salva risultati
results_alternative <- list(
  variability_correlations = cor_df,
  composite_correlations = cor_composite,
  delta_by_group = delta_by_group
)
saveRDS(results_alternative, "results/disinhibition_alternative_models.rds")

cat("\n✓ Analisi alternative completate.\n\n")
