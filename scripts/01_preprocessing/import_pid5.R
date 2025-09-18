# PID-5 scoring (questionario 12..233, 2 catch dentro ma esclusi dai punteggi)

library(tidyverse)
library(here)
library(rio)
library(stringr)

# ---------------------------
# 0) Import & cleaning user_id
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
# 1) Selezione colonne
#    - Item: 12..233 (inclusi i due catch alle posizioni 68 e 161)
#    - Rinomina item come i_1..i_222 (ordine del questionario)
# ---------------------------
item_idx <- 12:233
non_item_pos <- c(date = 1, user_id = 3, age = 4, sex = 5)

date_col <- names(temp)[non_item_pos["date"]]
user_col <- names(temp)[non_item_pos["user_id"]]
age_col <- names(temp)[non_item_pos["age"]]
sex_col <- names(temp)[non_item_pos["sex"]]

df_pid5 <- temp %>%
  dplyr::select(
    date = all_of(date_col),
    user_id = all_of(user_col),
    age = all_of(age_col),
    sex = all_of(sex_col),
    all_of(names(temp)[item_idx])
  )

# Rinomina sequenziale: i_1 .. i_222 (i_1 = col12, i_222 = col233)
new_item_names <- paste0("i_", seq_along(item_idx))
names(df_pid5)[
  (ncol(df_pid5) - length(new_item_names) + 1):ncol(df_pid5)
] <- new_item_names

# Coerzione numerica (atteso 0..3); non errore se ci sono NA
df_pid5 <- df_pid5 %>%
  mutate(across(all_of(new_item_names), ~ suppressWarnings(as.numeric(.x))))

# ---------------------------
# 2) Reverse coding 0..3 -> 3..0
#    (numerazione nell'ordine del questionario: i_1..i_222)
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
rev_cols <- paste0("i_", rev_items)

df_pid5 <- df_pid5 %>%
  mutate(across(any_of(rev_cols), ~ ifelse(is.na(.x), NA_real_, 3 - .x)))

# ---------------------------
# 3) Chiavi facet (ordine del questionario i_1..i_222)
#    NB: includono le posizioni 68 e 161 (careless) perché sono nel file,
#        ma verranno escluse dal calcolo
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
  unusual_beliefs_exp = c(94, 99, 106, 139, 143, 150, 194, 209), # 150 = catch
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
# 4) Calcolo facet escludendo i careless (posizioni nel questionario)
#    careless alle posizioni 68 e 161 => corrispondono a i_57 e i_150
#    MA qui lavoriamo in numerazione del questionario, quindi:
#    -> escludiamo semplicemente 68 e 161 dai vettori della facet al volo
# ---------------------------
careless_qpos <- c(68, 161)

df_scores <- df_pid5

# (opzionale) salvare numero di item validi per facet
for (facet in names(facet_items)) {
  idx <- setdiff(facet_items[[facet]], careless_qpos) # rimuovi catch
  cols <- paste0("i_", idx)
  facet_name <- paste0("facet_", facet)

  df_scores[[facet_name]] <- rowSums(
    df_scores[, cols, drop = FALSE],
    na.rm = TRUE
  )

  # conteggio item validi (QA opzionale)
  df_scores[[paste0(facet_name, "_nvalid")]] <-
    rowSums(!is.na(df_scores[, cols, drop = FALSE]))
}

# ---------------------------
# 5) Calcolo domini (somma delle facet)
# ---------------------------
for (dom in names(domain_facets)) {
  facet_cols <- paste0("facet_", domain_facets[[dom]])
  df_scores[[paste0("domain_", dom)]] <-
    rowSums(df_scores[, facet_cols, drop = FALSE], na.rm = TRUE)
}

# ---------------------------
# 6) (Opzionale) keep solo id + punteggi
# ---------------------------
df_out <- df_scores %>%
  select(date, user_id, starts_with("facet_"), starts_with("domain_"))

# ---------------------------
# 7) (Opzionale) flag careless (report, NON filtro)
# ---------------------------
# In numerazione del questionario: i_68 e i_161 sono i catch
att1 <- suppressWarnings(as.numeric(df_pid5[["i_68"]]))
att2 <- suppressWarnings(as.numeric(df_pid5[["i_161"]]))

catch_key <- c(1, 2) # risposte attese
careless_flag <- (is.na(att1) | att1 != catch_key[1]) &
  (is.na(att2) | att2 != catch_key[2])

careless_ids <- df_pid5$user_id[which(careless_flag)]
# (solo per QA: non usato per escludere soggetti)

# Modifica i nomi per coerenza con gli altri questionari
df1_out <- df_out |>
  dplyr::select(user_id, all_of(contains("domain_"))) |>
  dplyr::rename(
    domain_negative_affect_baseline = domain_negative_affect,
    domain_detachment_baseline = domain_detachment,
    domain_antagonism_baseline = domain_antagonism,
    domain_disinhibition_baseline = domain_disinhibition,
    domain_psychoticism_baseline = domain_psychoticism
  )

# ---------------------------
# 8) Export
# ---------------------------
rio::export(
  df_out,
  here::here("data", "processed", "pid5_scores.csv")
)
