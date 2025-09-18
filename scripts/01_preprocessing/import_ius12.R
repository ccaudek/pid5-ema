# Script: import_ius12.R
# Scopo: Importare e calcolare i punteggi dello IUS-12

library(tidyverse)
library(here)
library(rio)
library(stringr)

# Importa dati grezzi
raw_data <- rio::import(
  here::here(
    "data",
    "raw",
    "quest",
    "2025_EMA_quest2.xlsx"
  )
)

# Rinomina colonna user_id e seleziona items IUS-12 (colonne 25-36)
ius_df <- raw_data %>%
  rename(
    "user_id" = `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)`
  ) %>%
  dplyr::select(user_id, 25:36)

# Rinomina le colonne degli items
names(ius_df)[2:13] <- paste0("ius_", 1:12)

# Pulisci user_id (coerenza con altri script)
ius_df <- ius_df %>%
  mutate(
    user_id = user_id %>%
      str_to_lower() %>%
      str_replace_all("__+", "_") %>%
      str_replace("_$", "") %>%
      str_replace("^fla_me_", "fl_me_")
  )

# Funzione di scoring IUS-12
score_ius12 <- function(d) {
  # Prospective IU items (1-7)
  prospective_iu <- d %>%
    dplyr::select(ius_1, ius_2, ius_3, ius_4, ius_5, ius_6, ius_7) %>%
    mutate(across(everything(), as.numeric))

  d$prospective_iu <- rowSums(prospective_iu, na.rm = TRUE)

  # Inhibitory IU items (8-12)
  inhibitory_iu <- d %>%
    select(ius_8, ius_9, ius_10, ius_11, ius_12) %>%
    mutate(across(everything(), as.numeric))

  d$inhibitory_iu <- rowSums(inhibitory_iu, na.rm = TRUE)

  d$ius_tot <- d$prospective_iu + d$inhibitory_iu

  # Seleziona solo user_id e punteggi
  ius_scale <- d %>%
    select(user_id, inhibitory_iu, prospective_iu, ius_tot)

  return(ius_scale)
}

# Applica lo scoring
ius_scores <- score_ius12(ius_df)

ius_scores <- ius_scores |>
  dplyr::rename(
    inhibitory_ius_baseline = inhibitory_iu,
    prospective_ius_baseline = prospective_iu,
    ius_tot_baseline = ius_tot
  )

# Esporta
rio::export(
  ius_scores,
  here::here("data", "processed", "ius12_scores.csv")
)

# eof ---
