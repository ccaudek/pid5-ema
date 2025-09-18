# Script aggiornato: import_rosenberg.R (versione robusta)
library(tidyverse)
library(here)
library(rio)
library(stringr)

# -------------------------------
# Importa dati grezzi
# -------------------------------
raw_data <- rio::import(
  here::here("data", "raw", "quest", "2025_EMA_quest2.xlsx")
)

# Selezione colonna user_id (adatta la stringa se diversa)
raw_data <- raw_data %>%
  rename(
    user_id = `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)`
  )

# -------------------------------
# Individua colonne Rosenberg
# Due opzioni: per posizione (37:46) oppure cercando i nomi che iniziano con "1)" .. "10)".
# Controlla quale funziona meglio nel tuo file.
# -------------------------------

# Opzione A: per posizione (se sei sicuro che siano 37:46)
ros_cols_pos <- 37:46

# Opzione B: individuazione via regex sui nomi (se le colonne hanno "1) " all'inizio)
ros_idx_by_name <- grep("^\\s*\\d+\\)\\s*", names(raw_data))

# Scegli quale usare: preferisci la ricerca per nome se ros_idx_by_name non è vuoto
if (length(ros_idx_by_name) == 10) {
  ros_cols <- ros_idx_by_name
} else {
  ros_cols <- ros_cols_pos
}

# Estrai il dataframe con user_id e items Rosenberg
rosenberg_df <- raw_data %>%
  dplyr::select(user_id, all_of(ros_cols))

# Rinomina le colonne items in ros_1 ... ros_10
names(rosenberg_df)[2:11] <- paste0("ros_", 1:10)

# -------------------------------
# Normalizzazione e mappatura delle etichette testuali -> numerico 1..4
# Possibili etichette incontrate: "Fortemente d'accordo", "D'accordo",
# "In disaccordo", "Fortemente in disaccordo"
# Ma possono esserci variazioni di apostrofo, spazi, maiuscole.
# -------------------------------

# --- Funzione vettoriale: etichette -> 1..4
label_to_num <- function(x) {
  s <- as.character(x)
  # normalizza: trim, apostrofi tipografici, minuscole, spazi multipli, punto finale
  s <- s %>%
    stringr::str_squish() %>%
    stringr::str_replace_all("[\u2019\u2018]", "'") %>% # apostrofi tipografici -> '
    stringr::str_to_lower(locale = "it") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_remove("\\.$") # rimuovi eventuale punto finale

  dplyr::case_when(
    s == "fortemente d'accordo" ~ 4,
    s == "d'accordo" ~ 3,
    s == "in disaccordo" ~ 2,
    s == "fortemente in disaccordo" ~ 1,
    # eventuali abbreviazioni (se le usi davvero)
    s %in% c("fd", "forte d'accordo") ~ 4,
    s %in% c("d") ~ 3,
    s %in% c("id") ~ 2,
    s %in% c("fid") ~ 1,
    TRUE ~ NA_real_
  )
}

# Applica la funzione a tutte le colonne ros_
rosenberg_df <- rosenberg_df %>%
  mutate(across(starts_with("ros_"), ~ label_to_num(.x)))

# -------------------------------
# Funzione di scoring Rosenberg
# Reverse items: 3,5,8,9,10 (come nel tuo script originale)
# Nota: con scala 1..4 il reverse è 5 - x
# Impongo una soglia minima di item validi (es. almeno 8/10) per restituire punteggio,
# altrimenti NA; modifica 'min_items' se preferisci.
# -------------------------------
score_rosenberg <- function(d, min_items = 8L) {
  d2 <- d %>%
    # reverse coding più compatto
    mutate(across(
      c(ros_3, ros_5, ros_8, ros_9, ros_10),
      ~ if_else(is.na(.x), NA_real_, 5 - .x)
    )) %>%
    rowwise() %>%
    mutate(
      n_valid = sum(!is.na(c_across(ros_1:ros_10))),
      rosenberg_score_raw = if (n_valid >= min_items)
        sum(c_across(ros_1:ros_10), na.rm = TRUE) else NA_real_
    ) %>%
    ungroup() %>%
    dplyr::select(-n_valid)

  # Rinomina colonna risultato
  d2 <- d2 %>% rename(rosenberg_score = rosenberg_score_raw)
  # Se vuoi avere anche media per item: rosenberg_mean = rosenberg_score / number_valid_items (ma qui richiede contare di nuovo)
  return(d2 %>% dplyr::select(user_id, rosenberg_score, everything()))
}

# Applica lo scoring
rosenberg_scored <- score_rosenberg(rosenberg_df, min_items = 8)

# Se preferisci salvare solo user_id e rosenberg_score:
rosenberg_scores <- rosenberg_scored %>% dplyr::select(user_id, rosenberg_score)

rosenberg_scores <- rosenberg_scores |>
  rename(
    rosenberg_score_baseline = rosenberg_score
  )

# Esporta
rio::export(
  rosenberg_scores,
  here::here("data", "processed", "rosenberg_scores.csv")
)

# eof ---
