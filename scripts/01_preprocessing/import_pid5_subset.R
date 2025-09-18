# PID-5 domains excluding EMA items (and excluding catch items from scoring)
# - questionnaire columns 12..233 -> i_1..i_222
# - do NOT drop subjects
# - exclude from scoring: catch at 68 and 161; EMA set: i_1..i_15

library(tidyverse)
library(here)
library(rio)
library(stringr)

# ---------------------------
# 0) Import & user_id cleaning
# ---------------------------
temp1 <- rio::import(here::here("data", "raw", "quest", "PID-5_2025_EMA.xlsx"))

temp <- temp1 %>%
  mutate(
    `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)` = `Inserisci il tuo codice anonimo (esempio: ma_ro_1997_05_04_174)` %>%
      str_to_lower() %>%
      str_replace_all("__+", "_") %>%
      str_replace("_$", "") %>%
      str_replace("^fla_me_", "fl_me_")
  )

# ---------------------------
# 1) Selezione colonne e rinomina
#    - Item: 12..233 (inclusi catch 68 e 161)
#    - Rinominati come i_1..i_222 (ordine del questionario)
# ---------------------------
item_idx <- 12:233
non_item_pos <- c(date = 1, user_id = 3, age = 4, sex = 5)

df_pid5 <- temp %>%
  dplyr::select(
    date = all_of(names(temp)[non_item_pos["date"]]),
    user_id = all_of(names(temp)[non_item_pos["user_id"]]),
    age = all_of(names(temp)[non_item_pos["age"]]),
    sex = all_of(names(temp)[non_item_pos["sex"]]),
    all_of(names(temp)[item_idx])
  )

# i_1..i_222 (i_1 = col 12, ..., i_222 = col 233)
new_item_names <- paste0("i_", seq_along(item_idx))
names(df_pid5)[
  (ncol(df_pid5) - length(new_item_names) + 1):ncol(df_pid5)
] <- new_item_names

# coerzione numerica sicura
df_pid5 <- df_pid5 %>%
  mutate(across(all_of(new_item_names), ~ suppressWarnings(as.numeric(.x))))

# ---------------------------
# 2) Reverse coding 0..3 -> 3..0 (posizioni del questionario)
# ---------------------------
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
  mutate(across(paste0("i_", rev_items), ~ ifelse(is.na(.x), NA_real_, 3 - .x)))

# ---------------------------
# 3) Chiavi facet (posizioni del questionario i_1..i_222)
# ---------------------------
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
  distractibility = c(6, 29, 47, 68, 88, 118, 132, 144, 199), # 68 = catch
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
  unusual_beliefs_exp = c(94, 99, 106, 139, 143, 150, 194, 209), # (nota: qui 150 NON è catch; i catch sono 68,161)
  withdrawal = c(10, 20, 75, 82, 136, 146, 147, 161, 182, 186) # 161 = catch
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

# ---------------------------
# 4) Esclusioni dal calcolo
#    - careless nelle posizioni 68 e 161
#    - EMA items da rimuovere (pid5_1..pid5_15 -> i_1..i_15)
# ---------------------------
careless_qpos <- c(68, 161)
ema_exclude_qpos <- 1:15 # dal tuo elenco: i_1..i_15

exclude_qpos <- sort(unique(c(careless_qpos, ema_exclude_qpos)))

# ---------------------------
# 5) Calcolo facet e domini (somma), con n_valid per QA
# ---------------------------
df_scores <- df_pid5

for (facet in names(facet_items)) {
  idx_q <- setdiff(facet_items[[facet]], exclude_qpos)
  cols <- paste0("i_", idx_q)
  nm <- paste0("facet_", facet)

  if (length(cols) == 0) {
    df_scores[[nm]] <- NA_real_
    df_scores[[paste0(nm, "_nvalid")]] <- 0L
  } else {
    df_scores[[nm]] <- rowSums(df_scores[, cols, drop = FALSE], na.rm = TRUE)
    df_scores[[paste0(nm, "_nvalid")]] <- rowSums(
      !is.na(df_scores[, cols, drop = FALSE])
    )
  }
}

for (dom in names(domain_facets)) {
  fcols <- paste0("facet_", domain_facets[[dom]])
  df_scores[[paste0("domain_", dom)]] <-
    rowSums(df_scores[, fcols, drop = FALSE], na.rm = TRUE)
}

# ---------------------------
# 6) Output "snello": id + scores
# ---------------------------
df_out <- df_scores %>%
  dplyr::select(user_id, starts_with("domain_")) |>
  dplyr::rename(
    domain_negative_affect_baseline = domain_negative_affect,
    domain_detachment_baseline = domain_detachment,
    domain_antagonism_baseline = domain_antagonism,
    domain_disinhibition_baseline = domain_disinhibition,
    domain_psychoticism_baseline = domain_psychoticism
  )

# ---------------------------
# 7) (Opzionale) Flag careless (report, NON filtro)
# ---------------------------
att1 <- suppressWarnings(as.numeric(df_pid5[["i_68"]]))
att2 <- suppressWarnings(as.numeric(df_pid5[["i_161"]]))
catch_key <- c(1, 2)
careless_flag <- (is.na(att1) | att1 != catch_key[1]) &
  (is.na(att2) | att2 != catch_key[2])
careless_ids <- df_pid5$user_id[which(careless_flag)]
# non usati per escludere soggetti

# ---------------------------
# 8) Export
# ---------------------------
rio::export(
  df_out,
  here::here("data", "processed", "pid5_scores_noEMA.csv")
)
