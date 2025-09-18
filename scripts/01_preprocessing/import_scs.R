# --- Import & selezione
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(rio)
})

raw <- rio::import(here::here(
  "data",
  "raw",
  "quest",
  "2025_EMA_quest2.xlsx"
)) %>%
  dplyr::rename(
    user_id = `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)`
  )

# Colonne SCS: 47..73 esclusa 66 (careless)
scs_raw_cols <- setdiff(47:73, 66)

scs_df <- raw %>%
  dplyr::select(user_id, dplyr::all_of(scs_raw_cols))

# Rinomina nell’ordine in cui compaiono → scs_1:scs_26
names(scs_df)[2:27] <- paste0("scs_", 1:26)

# (Facoltativo ma consigliato) Coercizione numerica sicura e controllo range
scs_num <- scs_df %>%
  mutate(across(starts_with("scs_"), ~ suppressWarnings(as.numeric(.x))))

# Check veloce: valori fuori 1..5?
rng_chk <- scs_num %>%
  pivot_longer(starts_with("scs_")) %>%
  summarise(
    min_val = min(value, na.rm = TRUE),
    max_val = max(value, na.rm = TRUE),
    n_out = sum(!is.na(value) & (value < 1 | value > 5))
  )
print(rng_chk) # dovrebbe dare min>=1, max<=5, n_out=0

# --- Scoring (stessa logica tua)
scoring_scs <- function(d) {
  stopifnot("user_id" %in% names(d))

  # Positivi
  d$self_kindness <- d$scs_5 + d$scs_12 + d$scs_19 + d$scs_23 + d$scs_26
  d$common_humanity <- d$scs_3 + d$scs_7 + d$scs_10 + d$scs_15
  d$mindfulness <- d$scs_9 + d$scs_14 + d$scs_17 + d$scs_22

  # Negativi (reverse: 6 - x ≡ abs(x-6))
  d$self_judgment <- abs(d$scs_1 - 6) +
    abs(d$scs_8 - 6) +
    abs(d$scs_11 - 6) +
    abs(d$scs_16 - 6) +
    abs(d$scs_21 - 6)
  d$isolation <- abs(d$scs_4 - 6) +
    abs(d$scs_13 - 6) +
    abs(d$scs_18 - 6) +
    abs(d$scs_25 - 6)
  d$over_identification <- abs(d$scs_2 - 6) +
    abs(d$scs_6 - 6) +
    abs(d$scs_20 - 6) +
    abs(d$scs_24 - 6)

  d$neg_self_compassion <- d$self_judgment + d$isolation + d$over_identification
  d$pos_self_compassion <- d$self_kindness + d$common_humanity + d$mindfulness
  d$scs_ts <- d$neg_self_compassion + d$pos_self_compassion

  # Subscale NON-reversed per output (come richiesto nei tuoi commenti)
  d$sk <- d$self_kindness
  d$ch <- d$common_humanity
  d$mi <- d$mindfulness
  d$sj <- d$scs_1 + d$scs_8 + d$scs_11 + d$scs_16 + d$scs_21 # not reversed
  d$is <- d$scs_4 + d$scs_13 + d$scs_18 + d$scs_25 # not reversed
  d$oi <- d$scs_2 + d$scs_6 + d$scs_20 + d$scs_24 # not reversed

  tibble(
    user_id = d$user_id,
    self_kindness_baseline = d$sk,
    common_humanity_baseline = d$ch,
    mindfulness_baseline = d$mi,
    self_judgment_baseline = d$sj,
    isolation_baseline = d$is,
    over_identification_baseline = d$oi,
    scs_total_score_baseline = d$scs_ts
  )
}

scs_scores <- scoring_scs(scs_num)

# # (Opzionale) Soglia di validità sul totale (es. ≥ 20 item risposti)
# valid_counts <- scs_num %>%
#   transmute(user_id, n_valid = rowSums(!is.na(across(starts_with("scs_")))))
#
# scs_scores <- scs_scores %>%
#   left_join(valid_counts, by = "user_id") %>%
#   mutate(scs_total_score = ifelse(n_valid >= 20, scs_total_score, NA_real_)) %>%
#   select(-n_valid)

# Export
rio::export(scs_scores, here::here("data", "processed", "scs_scores.csv"))
