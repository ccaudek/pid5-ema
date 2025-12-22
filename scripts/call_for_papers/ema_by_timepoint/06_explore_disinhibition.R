# ==============================================================================
# ESPLORAZIONE COMPLETA: Dove emerge Disinhibition?
# ==============================================================================
# L'abstract riporta "Disinhibition F2: β=63 Hz [26,99]"
# Ma non specifica: mean vs std? quale vocale? stress vs recovery?
# Questo script esplora sistematicamente tutte le possibilità
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(cmdstanr)
  library(rio)
  library(here)
})

options(brms.backend = "cmdstanr")

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("ESPLORAZIONE: Effetti Disinhibition su TUTTI gli outcome\n")
cat(rep("=", 70), "\n\n", sep = "")

# ==============================================================================
# 1. CARICA MODELLI GIÀ FITTATI
# ==============================================================================

if (file.exists("results/fitted_models.rds")) {
  fitted_models <- readRDS("results/fitted_models.rds")
  cat("✓ Caricati", length(fitted_models), "modelli\n")
} else {
  stop("Esegui prima 03_moderation_analysis_CORRECTED.R")
}

# ==============================================================================
# 2. ESTRAI TUTTI GLI EFFETTI DISINHIBITION
# ==============================================================================

cat("\n=== Estrazione effetti Disinhibition ===\n")

extract_disinhibition_effects <- function(model, model_name) {
  fe <- fixef(model, summary = TRUE)

  # Trova tutte le righe con "disinhibition"
  dis_rows <- grep("disinhibition", rownames(fe), ignore.case = TRUE)

  if (length(dis_rows) == 0) return(NULL)

  dis_fe <- fe[dis_rows, , drop = FALSE]

  # Parse model name per outcome e vocale
  parts <- strsplit(gsub("^m_", "", model_name), "_")[[1]]
  vowel <- parts[length(parts)]
  outcome <- paste(parts[1:(length(parts) - 1)], collapse = "_")

  tibble(
    model = model_name,
    outcome = outcome,
    vowel = vowel,
    parameter = rownames(dis_fe),
    estimate = dis_fe[, "Estimate"],
    se = dis_fe[, "Est.Error"],
    ci_lower = dis_fe[, "Q2.5"],
    ci_upper = dis_fe[, "Q97.5"],
    significant = (ci_lower > 0) | (ci_upper < 0),
    # Calcola probability of direction
    effect_type = case_when(
      str_detect(parameter, "c1_stress") ~ "stress_moderation",
      str_detect(parameter, "c2_recovery") ~ "recovery_moderation",
      TRUE ~ "main_effect"
    )
  )
}

dis_effects <- map_dfr(names(fitted_models), function(mn) {
  extract_disinhibition_effects(fitted_models[[mn]], mn)
})

# ==============================================================================
# 3. RIEPILOGO EFFETTI DISINHIBITION
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("TUTTI GLI EFFETTI DISINHIBITION\n")
cat(rep("=", 70), "\n\n", sep = "")

# Ordina per dimensione dell'effetto (valore assoluto)
dis_effects <- dis_effects %>%
  mutate(abs_estimate = abs(estimate)) %>%
  arrange(desc(abs_estimate))

cat("--- Ordinati per dimensione effetto (|β|) ---\n\n")
print(
  dis_effects %>%
    select(
      outcome,
      vowel,
      effect_type,
      estimate,
      ci_lower,
      ci_upper,
      significant
    ) %>%
    head(20),
  n = 20
)

# ==============================================================================
# 4. EFFETTI SIGNIFICATIVI O QUASI-SIGNIFICATIVI
# ==============================================================================

cat("\n\n--- Effetti SIGNIFICATIVI (95% CI esclude zero) ---\n")
sig_dis <- dis_effects %>% filter(significant)
if (nrow(sig_dis) > 0) {
  print(
    sig_dis %>%
      select(outcome, vowel, effect_type, estimate, ci_lower, ci_upper)
  )
} else {
  cat("Nessun effetto significativo.\n")
}

cat("\n--- Effetti con 90% CI che esclude zero ---\n")
# Ricalcola con 90% CI
dis_90 <- map_dfr(names(fitted_models), function(mn) {
  model <- fitted_models[[mn]]
  post <- as_draws_df(model)

  # Trova colonne Disinhibition
  dis_cols <- grep("disinhibition", names(post), value = TRUE)
  dis_cols <- dis_cols[startsWith(dis_cols, "b_")]

  if (length(dis_cols) == 0) return(NULL)

  parts <- strsplit(gsub("^m_", "", mn), "_")[[1]]
  vowel <- parts[length(parts)]
  outcome <- paste(parts[1:(length(parts) - 1)], collapse = "_")

  map_dfr(dis_cols, function(col) {
    samples <- post[[col]]
    tibble(
      model = mn,
      outcome = outcome,
      vowel = vowel,
      parameter = gsub("^b_", "", col),
      estimate = mean(samples),
      ci90_lower = quantile(samples, 0.05),
      ci90_upper = quantile(samples, 0.95),
      pd = max(mean(samples > 0), mean(samples < 0)),
      sig_90 = (ci90_lower > 0) | (ci90_upper < 0)
    )
  })
})

sig_90 <- dis_90 %>% filter(sig_90)
if (nrow(sig_90) > 0) {
  print(
    sig_90 %>%
      select(outcome, vowel, parameter, estimate, ci90_lower, ci90_upper, pd)
  )
} else {
  cat("Nessun effetto con 90% CI che esclude zero.\n")
}

cat("\n--- Effetti con PD > 90% (direzione consistente) ---\n")
high_pd <- dis_90 %>% filter(pd > 0.90) %>% arrange(desc(pd))
if (nrow(high_pd) > 0) {
  print(high_pd %>% select(outcome, vowel, parameter, estimate, pd), n = 20)
} else {
  cat("Nessun effetto con PD > 90%.\n")
}

# ==============================================================================
# 5. FOCUS SU F2 (mean E std)
# ==============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("FOCUS: Disinhibition × F2 (mean E std)\n")
cat(rep("=", 70), "\n\n", sep = "")

f2_effects <- dis_effects %>%
  filter(str_detect(outcome, "f2"))

cat("--- F2 mean ---\n")
print(
  f2_effects %>%
    filter(outcome == "f2_mean") %>%
    select(vowel, effect_type, estimate, ci_lower, ci_upper, significant)
)

cat("\n--- F2 std ---\n")
print(
  f2_effects %>%
    filter(outcome == "f2_std") %>%
    select(vowel, effect_type, estimate, ci_lower, ci_upper, significant)
)

# ==============================================================================
# 6. CONFRONTO CON ALTRI TRATTI
# ==============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("CONFRONTO: Disinhibition vs altri tratti PID-5\n")
cat(rep("=", 70), "\n\n", sep = "")

# Estrai tutti gli effetti di moderazione per confronto
extract_all_moderation <- function(model, model_name) {
  fe <- fixef(model, summary = TRUE)

  # Trova interazioni (contengono ":")
  int_rows <- grep(":", rownames(fe))
  if (length(int_rows) == 0) return(NULL)

  int_fe <- fe[int_rows, , drop = FALSE]

  parts <- strsplit(gsub("^m_", "", model_name), "_")[[1]]
  vowel <- parts[length(parts)]
  outcome <- paste(parts[1:(length(parts) - 1)], collapse = "_")

  tibble(
    model = model_name,
    outcome = outcome,
    vowel = vowel,
    parameter = rownames(int_fe),
    estimate = int_fe[, "Estimate"],
    ci_lower = int_fe[, "Q2.5"],
    ci_upper = int_fe[, "Q97.5"],
    significant = (ci_lower > 0) | (ci_upper < 0)
  ) %>%
    mutate(
      trait = case_when(
        str_detect(parameter, "negative_aff") ~ "NA",
        str_detect(parameter, "detachment") ~ "Det",
        str_detect(parameter, "antagonism") ~ "Ant",
        str_detect(parameter, "disinhibition") ~ "Dis",
        str_detect(parameter, "psychoticism") ~ "Psy"
      ),
      contrast = case_when(
        str_detect(parameter, "c1_stress") ~ "stress",
        str_detect(parameter, "c2_recovery") ~ "recovery"
      )
    )
}

all_moderation <- map_dfr(names(fitted_models), function(mn) {
  extract_all_moderation(fitted_models[[mn]], mn)
})

# Conta effetti significativi per tratto
cat("Numero di effetti significativi per tratto PID-5:\n\n")
sig_by_trait <- all_moderation %>%
  filter(significant) %>%
  count(trait, contrast) %>%
  pivot_wider(names_from = contrast, values_from = n, values_fill = 0) %>%
  mutate(total = stress + recovery) %>%
  arrange(desc(total))

print(sig_by_trait)

cat("\n\nEffetti significativi per Disinhibition:\n")
dis_sig_details <- all_moderation %>%
  filter(significant, trait == "Dis")
if (nrow(dis_sig_details) > 0) {
  print(
    dis_sig_details %>%
      select(outcome, vowel, contrast, estimate, ci_lower, ci_upper)
  )
} else {
  cat("Nessuno.\n")
}

# ==============================================================================
# 7. PATTERN ANALYSIS: Disinhibition ha effetti consistenti?
# ==============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("PATTERN ANALYSIS: Consistenza effetti Disinhibition\n")
cat(rep("=", 70), "\n\n", sep = "")

# Per ogni outcome, guarda se la direzione è consistente across vowels
pattern_analysis <- dis_effects %>%
  filter(effect_type != "main_effect") %>%
  group_by(outcome, effect_type) %>%
  summarise(
    n_vowels = n(),
    n_positive = sum(estimate > 0),
    n_negative = sum(estimate < 0),
    mean_estimate = mean(estimate),
    min_estimate = min(estimate),
    max_estimate = max(estimate),
    consistent_direction = (n_positive == n_vowels) | (n_negative == n_vowels),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(mean_estimate)))

cat("Pattern di direzione attraverso le vocali:\n\n")
print(pattern_analysis)

cat(
  "\n\nOutcome con direzione CONSISTENTE (tutte vocali nella stessa direzione):\n"
)
consistent <- pattern_analysis %>% filter(consistent_direction)
if (nrow(consistent) > 0) {
  print(consistent)

  cat("\n\nDettagli per outcome con pattern consistente:\n")
  for (i in 1:nrow(consistent)) {
    out <- consistent$outcome[i]
    eff <- consistent$effect_type[i]

    cat("\n", out, " - ", eff, ":\n", sep = "")
    details <- dis_effects %>%
      filter(outcome == out, effect_type == eff) %>%
      select(vowel, estimate, ci_lower, ci_upper)
    print(details)
  }
} else {
  cat("Nessun pattern consistente.\n")
}

# ==============================================================================
# 8. CONCLUSIONI
# ==============================================================================

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("CONCLUSIONI\n")
cat(rep("=", 70), "\n", sep = "")

cat(
  "
L'abstract riporta: 'Disinhibition exhibited domain-specific articulatory 
alterations (F2: β=63 Hz [26,99])'

Possibili spiegazioni per la mancata replica:

1. EFFETTO SPECIFICO PER SOTTOCAMPIONE
   - L'abstract aveva N=111, ora hai più soggetti
   - L'effetto potrebbe essere stato un falso positivo
   - Oppure specifico per un sottogruppo (es. solo un corso)

2. DIFFERENZE METODOLOGICHE
   - Diversa aggregazione PID-5 (EMA vs Full)
   - Diversa gestione outlier
   - Diversi prior/modelli

3. L'EFFETTO ERA SU UNA VARIABILE DIVERSA
   - 'F2' potrebbe riferirsi a un composito
   - Oppure a una specifica vocale nella lettura della frase
   - O a F1 invece di F2 (typo nell'abstract?)

4. DISINHIBITION HA EFFETTI PIÙ DEBOLI DEGLI ALTRI TRATTI
   - Come mostrato dal confronto, NA/Det/Ant/Psy hanno più effetti significativi
   - Disinhibition potrebbe avere effetti più piccoli o più variabili

RACCOMANDAZIONE PER IL MANOSCRITTO:
- Riporta onestamente che l'effetto F2 non si replica
- Enfatizza gli effetti che SI replicano (NA, Det, Ant, Psy)
- Discuti le possibili ragioni nella sezione limitazioni
- Considera di rimuovere la claim su Disinhibition dall'abstract rivisto
"
)

# Salva risultati
saveRDS(dis_effects, "results/disinhibition_exploration.rds")
saveRDS(all_moderation, "results/all_moderation_effects.rds")

cat("\n✓ Esplorazione completata.\n\n")
