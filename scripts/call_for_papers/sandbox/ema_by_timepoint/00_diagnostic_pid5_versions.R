# ==============================================================================
# DIAGNOSTIC SCRIPT: Confronto versioni PID-5 per replicare abstract
# ==============================================================================
# Questo script aiuta a determinare quale versione dei predittori PID-5
# produce risultati più vicini all'abstract originale
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
cat(rep("=", 80), "\n", sep = "")
cat("DIAGNOSTIC: Confronto versioni PID-5\n")
cat(rep("=", 80), "\n\n", sep = "")

# ==============================================================================
# 1. Carica e ispeziona dati EMA
# ==============================================================================

data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
d <- rio::import(data_path)

cat("=== Struttura dati EMA ===\n")
cat("N righe totali:", nrow(d), "\n")
cat("N soggetti:", n_distinct(d$user_id), "\n\n")

cat("Distribuzione exam_period:\n")
print(table(d$exam_period, useNA = "ifany"))

cat("\n\nVariabili PID-5 disponibili:\n")
pid5_cols <- grep("pid5|domain.*baseline", names(d), value = TRUE)
print(pid5_cols)

# ==============================================================================
# 2. Confronta aggregazioni
# ==============================================================================

cat("\n\n=== Confronto metodi di aggregazione PID-5 ===\n\n")

# Versione A: Media di TUTTI i valori EMA
ema_all <- d %>%
  group_by(user_id) %>%
  summarise(
    na_ema_all = mean(pid5_negative_affectivity, na.rm = TRUE),
    det_ema_all = mean(pid5_detachment, na.rm = TRUE),
    ant_ema_all = mean(pid5_antagonism, na.rm = TRUE),
    dis_ema_all = mean(pid5_disinhibition, na.rm = TRUE),
    psy_ema_all = mean(pid5_psychoticism, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

# Versione B: Solo EMA del periodo baseline
ema_baseline <- d %>%
  filter(exam_period == "baseline") %>%
  group_by(user_id) %>%
  summarise(
    na_ema_bl = mean(pid5_negative_affectivity, na.rm = TRUE),
    det_ema_bl = mean(pid5_detachment, na.rm = TRUE),
    ant_ema_bl = mean(pid5_antagonism, na.rm = TRUE),
    dis_ema_bl = mean(pid5_disinhibition, na.rm = TRUE),
    psy_ema_bl = mean(pid5_psychoticism, na.rm = TRUE),
    n_obs_bl = n(),
    .groups = "drop"
  )

# Versione C: PID-5 completo baseline (domain_*_baseline)
pid5_full <- d %>%
  group_by(user_id) %>%
  summarise(
    na_full = first(na.omit(domain_negative_affect_baseline)),
    det_full = first(na.omit(domain_detachment_baseline)),
    ant_full = first(na.omit(domain_antagonism_baseline)),
    dis_full = first(na.omit(domain_disinhibition_baseline)),
    psy_full = first(na.omit(domain_psychoticism_baseline)),
    .groups = "drop"
  )

# Merge per confronto
comparison <- ema_all %>%
  left_join(ema_baseline, by = "user_id") %>%
  left_join(pid5_full, by = "user_id")

cat("Statistiche descrittive per versione:\n\n")

cat("--- Negative Affectivity ---\n")
cat(
  "EMA All:     M =",
  round(mean(comparison$na_ema_all, na.rm = TRUE), 2),
  ", SD =",
  round(sd(comparison$na_ema_all, na.rm = TRUE), 2),
  "\n"
)
cat(
  "EMA Baseline: M =",
  round(mean(comparison$na_ema_bl, na.rm = TRUE), 2),
  ", SD =",
  round(sd(comparison$na_ema_bl, na.rm = TRUE), 2),
  "\n"
)
cat(
  "PID-5 Full:  M =",
  round(mean(comparison$na_full, na.rm = TRUE), 2),
  ", SD =",
  round(sd(comparison$na_full, na.rm = TRUE), 2),
  "\n"
)

cat("\n--- Correlazioni tra versioni (NA) ---\n")
cor_na <- cor(
  comparison[, c("na_ema_all", "na_ema_bl", "na_full")],
  use = "complete.obs"
)
print(round(cor_na, 3))

cat("\n--- Detachment ---\n")
cat(
  "EMA All:     M =",
  round(mean(comparison$det_ema_all, na.rm = TRUE), 2),
  ", SD =",
  round(sd(comparison$det_ema_all, na.rm = TRUE), 2),
  "\n"
)
cat(
  "EMA Baseline: M =",
  round(mean(comparison$det_ema_bl, na.rm = TRUE), 2),
  ", SD =",
  round(sd(comparison$det_ema_bl, na.rm = TRUE), 2),
  "\n"
)
cat(
  "PID-5 Full:  M =",
  round(mean(comparison$det_full, na.rm = TRUE), 2),
  ", SD =",
  round(sd(comparison$det_full, na.rm = TRUE), 2),
  "\n"
)

# ==============================================================================
# 3. Verifica osservazioni per periodo
# ==============================================================================

cat("\n\n=== Osservazioni EMA per periodo ===\n")

obs_by_period <- d %>%
  group_by(user_id, exam_period) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = exam_period, values_from = n, values_fill = 0)

cat("\nDistribuzione N osservazioni per periodo:\n")
print(summary(obs_by_period[, -1]))

cat("\nSoggetti con almeno 1 osservazione per periodo:\n")
obs_by_period %>%
  mutate(
    has_baseline = baseline > 0,
    has_pre = pre_exam > 0,
    has_post = post_exam > 0,
    has_all = has_baseline & has_pre & has_post
  ) %>%
  summarise(
    n_baseline = sum(has_baseline),
    n_pre = sum(has_pre),
    n_post = sum(has_post),
    n_all_three = sum(has_all)
  ) %>%
  print()

# ==============================================================================
# 4. Raccomandazioni
# ==============================================================================

cat("\n\n")
cat(rep("=", 80), "\n", sep = "")
cat("RACCOMANDAZIONI\n")
cat(rep("=", 80), "\n", sep = "")

cat(
  "
L'abstract dice:
- 'PID-5 at baseline, then conducted intensive EMA over 2.5 months'
- 'each EMA PID-5 domain moderated stress reactivity'

Questo suggerisce che i PID-5 EMA dovrebbero essere usati come predittori
trait-like (between-person), aggregando su TUTTO il periodo EMA.

VERSIONE RACCOMANDATA: 'ema_all'
- Aggrega tutti i valori EMA per soggetto (non solo baseline)
- Crea misure trait-like stabili basate su 2.5 mesi di dati
- Coerente con la logica dell'abstract

VERSIONE PROBLEMATICA: 'ema_baseline' (attuale negli script)
- Usa solo i valori EMA del periodo 'baseline' (prima dell'esame)
- Ignora i dati EMA raccolti pre e post esame
- NON coerente con l'abstract

VERSIONE ALTERNATIVA: 'pid5_full'
- Usa i punteggi PID-5 completo (220 item) del baseline
- L'abstract menziona analisi supplementari che mostrano 'EMA measures 
  added predictive value beyond baseline assessment'
- Quindi questa versione dovrebbe essere usata per confronto, non come primaria

VERSIONE TIME-VARYING: 'ema_timevarying'
- Valori EMA diversi per ogni timepoint vocale
- Potrebbe aggiungere rumore perché ci sono solo 3 timepoint vocali
- Da testare come analisi di sensibilità

"
)

cat(
  "\n\n--- Per procedere, modifica PID5_VERSION in 02_..._CORRECTED.R ---\n\n"
)
