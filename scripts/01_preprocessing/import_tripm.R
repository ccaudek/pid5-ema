#' Dai nomi in italiano, l’ordine delle colonne 74:96 e 98:132 corrisponde
#' all’ordine standard degli item 1..58 del TriPM (esempi di match: 1 “Il più
#' delle volte sono ottimista” → I’m optimistic more often than not; 2 “Per me
#' è importante come si sentono gli altri” → How other people feel is important
#' to me; 3 “Spesso agisco in base a bisogni immediati” → I often act on
#' immediate needs). Le chiavi di sottoscala e degli item “reverse” sono quelle
#' canoniche nei manuali/scoring key ufficiali.

#' Lo script import_tripm.R:
#' - seleziona le colonne 74:96, 98:132 (escludendo 97),
#' - le rinomina tripm_1:tripm_58 nell’ordine originale,
#' - normalizza la codifica alle 4 categorie 0..3 (se il file usa 1..4 sottrae 1;
#' se è già 0..3 lascia così),
#' - applica la reversione sugli item [F] standard,
#' - calcola Boldness (19 item) e Meanness (19 item),
#'- esporta in data/processed/tripm_scores.csv.
#'
#'Reverse-coded [F]: 2, 4, 10, 11, 16, 21, 25, 30, 33, 35, 39, 41, 44, 47, 50, 52, 57
#'Boldness: 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 32, 35, 38, 41, 44, 47, 50, 54, 57
#'Meanness: 2, 6, 8, 11, 14, 17, 20, 23, 26, 29, 33, 36, 39, 40, 42, 45, 48, 52, 55.

# import_tripm.R — TriPM (Boldness & Meanness), escluso careless
# - Seleziona colonne 74:96 e 98:132 (97 = careless, esclusa)
# - Rinomina in tripm_1:tripm_58
# - Ricodifica IT -> 0..3 ("vera"=3, "abbastanza vera"=2, "abbastanza falsa"=1, "falsa"=0)
# - Reverse sugli item [F]
# - Calcola Boldness e Meanness + n_valid_* nella stessa pipeline
# - (Opzionale) tripm_total = boldness + meanness
# - Esporta CSV

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(rio)
  library(here)
})

# 0) Import e user_id
raw <- rio::import(here::here(
  "data",
  "raw",
  "quest",
  "2025_EMA_quest2.xlsx"
)) %>%
  dplyr::rename(
    user_id = `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)`
  )

# 1) Selezione item TriPM: 74..96 e 98..132 (97 = careless escluso)
tripm_cols <- c(74:96, 98:132)
tripm <- raw %>% dplyr::select(user_id, dplyr::all_of(tripm_cols))

# 2) Rinomina in tripm_1 .. tripm_58 nell'ordine del file
names(tripm)[2:59] <- paste0("tripm_", 1:58)

# 3) Ricodifica IT -> 0..3 (vettoriale, senza warning)
#    Mappa: "vera"=3; "abbastanza vera"=2; "abbastanza falsa"=1; "falsa"=0
#    Riconosce prefisso "questa affermazione è ..." e numeri come stringhe ("0..3" o "1..4")
to_03_it <- function(x) {
  s <- as.character(x) |>
    stringr::str_squish() |>
    stringr::str_replace_all("[\u2019\u2018]", "'") |>
    stringr::str_to_lower(locale = "it") |>
    stringr::str_remove("^questa affermazione è\\s+")
  out <- rep(NA_real_, length(s))
  out[s == "vera"] <- 3
  out[s == "abbastanza vera"] <- 2
  out[s == "abbastanza falsa"] <- 1
  out[s == "falsa"] <- 0
  idx03 <- s %in% c("0", "1", "2", "3")
  out[idx03] <- as.numeric(s[idx03])
  idx14 <- s %in% c("1", "2", "3", "4")
  out[idx14] <- as.numeric(s[idx14]) - 1
  out
}

# 4) Applica ricodifica a tutte le colonne tripm_*
tripm_num <- tripm %>%
  dplyr::mutate(across(starts_with("tripm_"), to_03_it))

# (Facoltativo) Controllo di qualità: valori non mappati (lista deve risultare vuota)
# tripm %>%
#   tidyr::pivot_longer(starts_with("tripm_")) %>%
#   dplyr::mutate(mapped = to_03_it(value)) %>%
#   dplyr::filter(is.na(mapped) & !is.na(value)) %>%
#   dplyr::count(value, sort = TRUE)

# 5) Reverse sugli item [F] (codifica 0..3 -> 3..0)
rev_items <- c(2, 4, 10, 11, 16, 21, 25, 30, 33, 35, 39, 41, 44, 47, 50, 52, 57)
tripm_rev <- tripm_num %>%
  dplyr::mutate(across(
    paste0("tripm_", rev_items),
    ~ ifelse(is.na(.x), NA_real_, 3 - .x)
  ))

# 6) Indici sottoscale (solo Boldness e Meanness)
bold_items <- c(
  1,
  4,
  7,
  10,
  13,
  16,
  19,
  22,
  25,
  28,
  32,
  35,
  38,
  41,
  44,
  47,
  50,
  54,
  57
) # 19 item
mean_items <- c(
  2,
  6,
  8,
  11,
  14,
  17,
  20,
  23,
  26,
  29,
  33,
  36,
  39,
  40,
  42,
  45,
  48,
  52,
  55
) # 19 item

row_sum_idx <- function(df, idx)
  rowSums(df[paste0("tripm_", idx)], na.rm = TRUE)

# 7) Calcolo score + conteggio item validi (tutto in un'unica transmute: no join, no duplicati)
tripm_scores <- tripm_rev %>%
  dplyr::transmute(
    user_id,
    boldness_baseline = row_sum_idx(cur_data(), bold_items),
    meanness_baseline = row_sum_idx(cur_data(), mean_items)
    # n_valid_bold = rowSums(!is.na(across(paste0("tripm_", bold_items)))),
    # n_valid_mean = rowSums(!is.na(across(paste0("tripm_", mean_items))))
  )
# %>%
# (Opzionale) totale riassuntivo senza Disinhibition
# dplyr::mutate(tripm_total = boldness + meanness)

# (Opzionale) Soglie minime per validità (esempio: ≥ 15/19 item)
# tripm_scores <- tripm_scores %>%
#   dplyr::mutate(
#     boldness = ifelse(n_valid_bold >= 15, boldness, NA_real_),
#     meanness = ifelse(n_valid_mean >= 15, meanness, NA_real_)
#   )

# (Opzionale) Se desideri forzare una riga per user_id (in caso di duplicati nel raw)
# tripm_scores <- tripm_scores %>% dplyr::distinct(user_id, .keep_all = TRUE)

# 8) Export
rio::export(tripm_scores, here::here("data", "processed", "tripm_scores.csv"))

# Fine script
