# Script: import_dass21.R
# Scopo: Importare e score il questionario DASS-21 dai dati grezzi

library(tidyverse)
library(here)
library(rio)
library(stringr)

# Importa il file grezzo
raw_data <- rio::import(
  here::here(
    "data",
    "raw",
    "quest",
    "2025_EMA_quest2.xlsx"
  )
)

# Rinomina la colonna dell'user_id
dass_df <- raw_data %>%
  rename(
    "user_id" = `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)`
  ) %>%
  dplyr::select(user_id, everything())

# Seleziona le colonne del DASS-21 (da 3 a 20 e 22 a 24)
dass_items <- dass_df[, c(3:20, 22:24)]

# Funzione per mappare le risposte testuali in numeri
map_dass21_response <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x0 <- tolower(str_squish(as.character(x)))
  dict <- c(
    "non mi è mai accaduto" = 0,
    "mi è capitato qualche volta" = 1,
    "mi è capitato con una certa frequenza" = 2,
    "mi è capitato quasi sempre" = 3
  )
  out <- unname(dict[x0])
  bad <- is.na(out) & !is.na(x0)
  if (any(bad)) out[bad] <- suppressWarnings(as.numeric(x[bad]))
  as.numeric(out)
}

# Applica la funzione a tutte le colonne degli item
dass_items_numeric <- dass_items %>%
  mutate(across(everything(), map_dass21_response))

# Assegna nomi alle colonne per gli item
names(dass_items_numeric) <- paste0("dass21_", 1:21)

# Definisci gli indici delle sottoscale
stress_idx <- c(1, 6, 8, 11, 12, 14, 18)
anxiety_idx <- c(2, 4, 7, 9, 15, 19, 20)
depression_idx <- c(3, 5, 10, 13, 16, 17, 21)

# Funzione per calcolare il punteggio di una sottoscale
calculate_scale_score <- function(data, indices) {
  rowSums(data[, indices, drop = FALSE], na.rm = FALSE)
}

# Calcola i punteggi
dass_scores <- tibble(
  user_id = dass_df$user_id,
  dass_stress_baseline = calculate_scale_score(dass_items_numeric, stress_idx),
  dass_anxiety_baseline = calculate_scale_score(
    dass_items_numeric,
    anxiety_idx
  ),
  dass_depression_baseline = calculate_scale_score(
    dass_items_numeric,
    depression_idx
  )
)

# Pulisci l'user_id (per coerenza con gli altri script)
dass_scores <- dass_scores %>%
  mutate(
    user_id = user_id %>%
      str_to_lower() %>%
      str_replace_all("__+", "_") %>%
      str_replace("_$", "") %>%
      str_replace("^fla_me_", "fl_me_")
  )

# Salva il file
rio::export(
  dass_scores,
  here::here("data", "processed", "dass21_scores.csv")
)

# eof ---
