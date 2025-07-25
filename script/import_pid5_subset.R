# Overview ----------------------------------------------------------------
# Associated project:
# Script purpose: Import the data and code the 5 domains of the PID-5, after
#   removing the 3 items for each domain that are used in the EMA procedure.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-07-25
# Last update:
# Status: In progress
# Notes:
# Remove these items: "pid5_13"                   "pid5_15"
# "pid5_11"                   "pid5_3"                    "pid5_2"
# "pid5_7"                    "pid5_14"                   "pid5_6"
# "pid5_4"                    "pid5_12"                   "pid5_1"
# "pid5_9"                    "pid5_5"                    "pid5_8"
# "pid5_10"

# Load necessary libraries ------------------------------------------------

library(tidyverse)
library(here)
library(rio)
library(stringr)


temp1 <- rio::import(
  here::here(
    "data",
    "raw",
    "quest",
    "PID-5_2025_EMA.xlsx"
  )
)

temp <- temp1 %>%
  mutate(
    `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)` = `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)` %>%
      str_to_lower() %>% # tutto in minuscolo
      str_replace_all("__+", "_") %>% # "__" → "_"
      str_replace("_$", "") %>% # rimuove "_" finale
      str_replace("^fla_me_", "fl_me_") # fla_me_ → fl_me_
  )

# indici delle colonne item (escludo le due domande di attenzione 68, 161)
item_idx <- setdiff(12:233, c(68, 161))

df_pid5 <- temp %>%
  dplyr::select(
    date = 1,
    user_id = 3,
    age = 4,
    sex = 5,
    all_of(item_idx)
  ) %>%
  rename_with(~ paste0("i_", seq_along(.)), -(1:4))


# 1. Define items to remove
to_be_removed <- c(
  "i_13",
  "i_15",
  "i_11",
  "i_3",
  "i_2",
  "i_7",
  "i_14",
  "i_6",
  "i_4",
  "i_12",
  "i_1",
  "i_9",
  "i_5",
  "i_8",
  "i_10"
)

# 2. Reverse-score specified items
rev_items <- c(
  7,
  30,
  35,
  58,
  87,
  90,
  96,
  97,
  98,
  131,
  142,
  155,
  164,
  177,
  210,
  215
)
df_pid5 <- df_pid5 %>%
  mutate(across(
    all_of(paste0("i_", rev_items)),
    ~ ifelse(is.na(.x), NA_real_, 3 - as.numeric(.x))
  ))

# 3. Define original facet and domain mappings (as before)
facet_items <- list(
  anhedonia = c(1, 23, 26, 30, 124, 155, 157, 189),
  anxiousness = c(79, 93, 95, 96, 109, 110, 130, 141, 174),
  attention_seeking = c(14, 43, 74, 111, 113, 173, 191, 211),
  callousness = c(
    11,
    13,
    19,
    54,
    72,
    73,
    90,
    153,
    166,
    183,
    198,
    200,
    207,
    208
  ),
  deceitfulness = c(41, 53, 56, 76, 126, 134, 142, 206, 214, 218),
  depressivity = c(
    27,
    61,
    66,
    81,
    86,
    104,
    119,
    148,
    151,
    163,
    168,
    169,
    178,
    212
  ),
  distractibility = c(6, 29, 47, 68, 88, 118, 132, 144, 199),
  eccentricity = c(5, 21, 24, 25, 33, 52, 55, 70, 71, 152, 172, 185, 205),
  emotional_lability = c(18, 62, 102, 122, 138, 165, 181),
  grandiosity = c(40, 65, 114, 179, 187, 197),
  hostility = c(28, 32, 38, 85, 92, 116, 158, 170, 188, 216),
  impulsivity = c(4, 16, 17, 22, 58, 204),
  intimacy_avoidance = c(89, 97, 108, 120, 145, 203),
  irresponsibility = c(31, 129, 156, 160, 171, 201, 210),
  manipulativeness = c(107, 125, 162, 180, 219),
  perceptual_dysregulation = c(
    36,
    37,
    42,
    44,
    59,
    77,
    83,
    154,
    192,
    193,
    213,
    217
  ),
  perseveration = c(46, 51, 60, 78, 80, 100, 121, 128, 137),
  restricted_affectivity = c(8, 45, 84, 91, 101, 167, 184),
  rigid_perfectionism = c(34, 49, 105, 115, 123, 135, 140, 176, 196, 220),
  risk_taking = c(3, 7, 35, 39, 48, 67, 69, 87, 98, 112, 159, 164, 195, 215),
  separation_insecurity = c(12, 50, 57, 64, 127, 149, 175),
  submissiveness = c(9, 15, 63, 202),
  suspiciousness = c(2, 103, 117, 131, 133, 177, 190),
  unusual_beliefs_exp = c(94, 99, 106, 139, 143, 150, 194, 209),
  withdrawal = c(10, 20, 75, 82, 136, 146, 147, 161, 182, 186)
)

domain_facets <- list(
  negative_affect = c(
    "emotional_lability",
    "anxiousness",
    "separation_insecurity"
  ),
  detachment = c("withdrawal", "anhedonia", "intimacy_avoidance"),
  antagonism = c("manipulativeness", "deceitfulness", "grandiosity"),
  disinhibition = c("irresponsibility", "impulsivity", "distractibility"),
  psychoticism = c(
    "unusual_beliefs_exp",
    "eccentricity",
    "perceptual_dysregulation"
  )
)

# 4. Remove the EMA items from all facets
to_be_removed_nums <- as.integer(gsub("i_", "", to_be_removed))
facet_items_filtered <- lapply(facet_items, function(items) {
  setdiff(items, to_be_removed_nums)
})

# 5. Recompute the scores excluding the EMA items
df_scores <- df_pid5

# Facet scores
for (facet in names(facet_items_filtered)) {
  item_names <- paste0("i_", facet_items_filtered[[facet]])
  df_scores[[paste0("facet_", facet)]] <- rowSums(
    df_scores[, item_names],
    na.rm = TRUE
  )
}

# Domain scores
for (dom in names(domain_facets)) {
  facet_names <- paste0("facet_", domain_facets[[dom]])
  df_scores[[paste0("domain_", dom)]] <- rowSums(
    df_scores[, facet_names],
    na.rm = TRUE
  )
}

# $ domain_negative_affect         <dbl> 24, 17, 56, 44, 38, 10, 34, 19, 50, 15, 30, 41, 24,…
# $ domain_detachment              <dbl> 13, 12, 11, 12, 12, 1, 3, 28, 23, 0, 33, 13, 12, 17…
# $ domain_antagonism              <dbl> 8, 7, 11, 6, 22, 4, 24, 5, 1, 3, 7, 9, 7, 5, 26, 8,…
# $ domain_disinhibition           <dbl> 8, 5, 12, 7, 14, 10, 40, 15, 6, 6, 33, 17, 14, 14, …
# $ domain_psychoticism            <dbl> 14, 6, 57, 24, 40, 9, 34, 26, 30, 9, 22, 23, 27, 5,…

#--- 3. (facoltativo) Seleziona solo id + punteggi ----------------------------
df_scores <- df_scores %>%
  dplyr::select(date, user_id, starts_with("facet_"), starts_with("domain_"))


# Careless responding ---------------------------------------------------------

catch_idx <- c(68, 161) # posizioni dei due catch‑item
catch_key <- c(1, 2) # risposte corrette

att1 <- as.numeric(temp[[catch_idx[1]]])
att2 <- as.numeric(temp[[catch_idx[2]]])

wrong_rows <- which(
  (is.na(att1) | att1 != catch_key[1]) &
    (is.na(att2) | att2 != catch_key[2])
)

wrong_ids <- df_pid5$user_id[wrong_rows]

user_id_with_careless_responding <- c(
  "ma_se_2005_11_14_490",
  "reve20041021036",
  "di_ma_2005_10_20_756",
  "pa_sc_2005_09_10_468",
  "il_re_2006_01_18_645",
  "so_ma_2003_10_13_804",
  "lo_ca_2005_05_07_05_437",
  "va_ma_2005_05_31_567",
  "no_un_2005_06_29_880",
  "an_bo_1988_08_24_166",
  "st_ma_2004_04_21_426",
  "an_st_2005_10_16_052",
  "vi_de_2002_12_30_067",
  "gi_ru_2005_03_08_033",
  "al_mi_2005_03_05_844",
  "la_ma_2006_01_31_787",
  "gi_lo_2004_06_27_237",
  "ch_bi_2001_01_28_407",
  "al_pe_2001_04_20_079",
  "le_de_2003_09_05_067",
  "fe_gr_2002_02_19_434",
  "ma_ba_2002_09_09_052",
  "ca_gi_2003_09_16_737",
  "an_to_2003_08_06_114",
  "al_se_2003_07_28_277",
  "ja_tr_2002_10_06_487",
  "el_ci_2002_02_15_057",
  "se_ti_2000_03_04_975",
  "co_ga_2003_10_29_614",
  "al_ba_2003_18_07_905",
  "bi_ro_2003_09_07_934",
  "an_va_2004_04_08_527",
  "ev_cr_2003_01_27_573"
)


rio::export(
  df_scores,
  here::here(
    "data",
    "processed",
    "pid5.csv"
  )
)
