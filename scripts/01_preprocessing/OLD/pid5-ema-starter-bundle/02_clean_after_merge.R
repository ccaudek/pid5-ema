# 02_clean_after_merge.R
# Data cleaning post-merge EMA + baseline
# - Esclusione careless responding (lista a priori)
# - Filtro numerosità EMA: 5 <= n_ema <= 40
# - Import metadati (genere) + tag esami (exam_period)
# - QA essenziale + export dataset finale e log esclusioni

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(lubridate)
  library(readr)
  library(rio)
})

# ============================
# 0) PERCORSI / INPUT
# ============================

# File "raw totale" già salvato dallo script precedente (RDS/CSV)
# -> Se esistono entrambi, preferisco RDS
INPUT_RDS <- here::here("data", "processed", "ema_plus_scales_merged.RDS")
INPUT_CSV <- here::here("data", "processed", "ema_plus_scales_merged.csv")

# Metadati (almeno: user_id, genere) — adatta il percorso/nome file
META_PATH <- here::here("data", "raw", "meta", "all_combined_sex_NEW_1.xlsx")
# deve contenere colonne: user_id, genere

# Tag esami (finestre) — adatta il percorso/nome file
# formato atteso: start_date, end_date, exam_period (YYYY-MM-DD)
EXAM_TAGS_PATH <- here::here("data", "raw", "meta", "exam_periods.csv")

# Output
OUT_CLEAN_RDS <- here::here("data", "processed", "ema_plus_scales_cleaned.rds")
OUT_CLEAN_CSV <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
OUT_LOG_EXCL <- here::here("data", "processed", "ema_exclusion_log.csv")
OUT_QA_TXT <- here::here("data", "processed", "ema_cleaning_QA.txt")

# ============================
# 1) LETTURA RAW TOTALE
# ============================

dat_raw <- if (file.exists(INPUT_RDS)) {
  readRDS(INPUT_RDS)
} else if (file.exists(INPUT_CSV)) {
  suppressMessages(readr::read_csv(INPUT_CSV, show_col_types = FALSE))
} else {
  stop(
    "Non trovo il file raw totale (né RDS né CSV). Aggiorna INPUT_RDS / INPUT_CSV."
  )
}

# Normalizzo alcuni nomi tipici (se presenti)
dat_raw <- dat_raw %>%
  rename(
    user_id = any_of(c("user_id", "UserID", "id_anon")),
    date = any_of(c("date", "Date", "baseline_date")),
    ema_time = any_of(c("ema_time", "timestamp", "time", "datetime"))
  )

# ============================
# 2) LISTA CARELESS A PRIORI
#    (presa dallo script allegato)
# ============================

careless_ids <- c(
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

# ============================
# 3) CALCOLO n_ema (robusto)
# ============================

# ============================
# 3) CALCOLO n_ema (robusto)
# ============================

if ("n_ema" %in% names(dat_raw)) {
  n_ema_tbl <- dat_raw %>% dplyr::distinct(user_id, n_ema)
} else {
  # colonne candidate a indicare una riga EMA
  ind_candidates <- c("ema_time", "day", "beep", "ema_id", "ema_wave")
  present_inds <- intersect(ind_candidates, names(dat_raw))

  # matrice logica: per ogni colonna presente TRUE se non NA
  flag_mat <- NULL
  if (length(present_inds) > 0) {
    flag_mat <- sapply(present_inds, function(nm) !is.na(dat_raw[[nm]]))
    # sapply con 1 sola colonna restituisce un vettore: forziamo matrice
    if (!is.matrix(flag_mat)) flag_mat <- matrix(flag_mat, ncol = 1)
  }

  # flag da "source" se esiste ed è testuale con "ema" dentro
  source_flag <- if ("source" %in% names(dat_raw)) {
    grepl("ema", dat_raw[["source"]], ignore.case = TRUE)
  } else {
    rep(FALSE, nrow(dat_raw))
  }

  # combina gli indicatori senza mai referenziare colonne mancanti
  if (is.null(flag_mat)) {
    is_ema <- source_flag
  } else {
    is_ema <- (rowSums(flag_mat) > 0) | source_flag
  }

  dat_raw$.is_ema <- is_ema

  n_ema_tbl <- dat_raw %>%
    dplyr::group_by(user_id) %>%
    dplyr::summarise(n_ema = sum(.is_ema, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(n_ema = as.integer(n_ema))
}

# Attacca n_ema al dataset principale (se non c'è già)
if (!"n_ema" %in% names(dat_raw)) {
  dat_raw <- dat_raw %>% dplyr::left_join(n_ema_tbl, by = "user_id")
}

# ============================
# 4) METADATI (genere) + TAG ESAMI (exam_period)
# ============================

# ============================
# 4) METADATI (genere, corso) + TAG ESAMI (exam_period)
# ============================

# -- 4a) Metadati: user_id, course, sex (data in meta è "dd_mm_yyyy", qui NON usata per il tag)
meta_df <- suppressMessages(rio::import(META_PATH)) %>%
  dplyr::rename(
    user_id = any_of(c("subj_code", "user_id", "id_anon")),
    date_meta = any_of(c("date", "Date", "data"))
  ) %>%
  dplyr::mutate(
    # parsiamo la data se presente (non è necessaria per il tag esami, ma utile per QA)
    date_meta = suppressWarnings(as.Date(date_meta, format = "%d_%m_%Y")),
    # normalizziamo il nome corso
    course = dplyr::case_when(
      stringr::str_detect(course, regex("^psico", ignore_case = TRUE)) ~
        "Psicometria",
      stringr::str_detect(course, regex("^test", ignore_case = TRUE)) ~
        "Testing",
      stringr::str_detect(course, regex("^inter", ignore_case = TRUE)) ~
        "Interventi",
      stringr::str_detect(course, regex("^clin", ignore_case = TRUE)) ~
        "Clinica",
      TRUE ~ course
    ),
    sex = dplyr::case_when(
      stringr::str_detect(sex, regex("^f", ignore_case = TRUE)) ~ "Femmina",
      stringr::str_detect(sex, regex("^m", ignore_case = TRUE)) ~ "Maschio",
      TRUE ~ sex
    )
  ) %>%
  dplyr::select(user_id, course, sex, date_meta) %>%
  dplyr::distinct()

# -- 4b) Join metadati
# NB: dat_final deve esistere già dal passo precedente (dopo esclusioni careless + n_ema)
dat_final <- dat_raw %>%
  dplyr::left_join(meta_df, by = "user_id")

# -- 4c) Definizione date d’esame per corso (baseline altrimenti)
psico_pre <- as.Date(c("2025-04-14", "2025-05-21"))
psico_post <- as.Date(c("2025-04-15", "2025-05-22"))

test_pre <- as.Date(c("2025-04-14", "2025-05-25"))
test_post <- as.Date(c("2025-04-15", "2025-05-26"))

interv_pre <- as.Date("2025-05-12")
interv_post <- as.Date("2025-05-13")

# -- 4d) Tag exam_period in base al corso e alla data EMA (colonna 'day')
dat_final <- dat_final %>%
  dplyr::mutate(
    day = as.Date(day),
    exam_period = dplyr::case_when(
      course == "Clinica" ~ "baseline",

      course == "Psicometria" & day %in% psico_pre ~ "pre_exam",
      course == "Psicometria" & day %in% psico_post ~ "post_exam",

      course == "Testing" & day %in% test_pre ~ "pre_exam",
      course == "Testing" & day %in% test_post ~ "post_exam",

      course == "Interventi" & day %in% interv_pre ~ "pre_exam",
      course == "Interventi" & day %in% interv_post ~ "post_exam",

      TRUE ~ "baseline"
    ),
    exam_period = factor(
      exam_period,
      levels = c("baseline", "pre_exam", "post_exam")
    )
  )

# -- 4e) QA rapido (opzionale)
dplyr::count(dat_final, course, exam_period)
table(dat_final$exam_period)

# ============================
# 5) ESCLUSIONI
#   A) careless responding (lista a priori)
#   B) numerosità EMA: 5 <= n_ema <= 40
# ============================

# (A) careless
excl_careless <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(user_id %in% careless_ids) %>%
  mutate(reason = "careless_responding (a priori)")

# (B) n_ema fuori range
excl_nema <- n_ema_tbl %>%
  dplyr::filter(is.na(n_ema) | n_ema < 5 | n_ema > 40) %>%
  mutate(reason = "n_ema_out_of_range (must be 5..40)")

# Log esclusioni unificato
excl_log <- bind_rows(excl_careless, excl_nema) %>%
  dplyr::distinct(user_id, reason)

# Soggetti da mantenere
keep_ids <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(!(user_id %in% excl_log$user_id)) %>%
  pull(user_id)

dat_clean <- dat_final %>%
  dplyr::filter(user_id %in% keep_ids)

# ============================
# 6) QA ESSENZIALE (stampata e salvata su file)
# ============================

qa_lines <- c()

n_subj_raw <- n_distinct(dat_raw$user_id)
n_subj_clean <- n_distinct(dat_clean$user_id)
n_excl <- n_distinct(excl_log$user_id)

qa_lines <- c(
  qa_lines,
  sprintf("Totale soggetti (raw): %d", n_subj_raw),
  sprintf("Esclusi (tot): %d", n_excl),
  sprintf("   - Careless: %d", n_distinct(excl_careless$user_id)),
  sprintf("   - n_ema fuori [5,40]: %d", n_distinct(excl_nema$user_id)),
  sprintf("Soggetti inclusi (finale): %d", n_subj_clean),
  ""
)

# Distribuzione n_ema
nema_keep <- n_ema_tbl %>% filter(user_id %in% keep_ids)
qa_lines <- c(
  qa_lines,
  sprintf(
    "n_ema (included) — min: %d, q1: %d, median: %d, mean: %.2f, q3: %d, max: %d",
    min(nema_keep$n_ema, na.rm = TRUE),
    quantile(nema_keep$n_ema, 0.25, na.rm = TRUE) %>% as.integer(),
    median(nema_keep$n_ema, na.rm = TRUE) %>% as.integer(),
    mean(nema_keep$n_ema, na.rm = TRUE),
    quantile(nema_keep$n_ema, 0.75, na.rm = TRUE) %>% as.integer(),
    max(nema_keep$n_ema, na.rm = TRUE)
  )
)

# Genere (se presente)
if ("sex" %in% names(dat_clean)) {
  gen_tab <- dat_clean %>%
    distinct(user_id, sex) %>%
    count(sex, sort = TRUE)
  qa_lines <- c(qa_lines, "", "Distribuzione genere (soggetti unici):")
  qa_lines <- c(qa_lines, paste0(" - ", gen_tab$sex, ": ", gen_tab$n))
}

# Exam period (se disponibile)
if ("exam_period" %in% names(dat_clean)) {
  ex_tab <- dat_clean %>%
    distinct(user_id, exam_period) %>%
    count(exam_period, sort = TRUE)
  qa_lines <- c(qa_lines, "", "Distribuzione exam_period (soggetti unici):")
  qa_lines <- c(
    qa_lines,
    paste0(
      " - ",
      ifelse(is.na(ex_tab$exam_period), "NA", ex_tab$exam_period),
      ": ",
      ex_tab$n
    )
  )
}

# (A) careless
excl_careless <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(user_id %in% careless_ids) %>%
  dplyr::left_join(n_ema_tbl, by = "user_id") %>%
  dplyr::mutate(reason = "careless_responding (a priori)")

# (B) n_ema fuori range
excl_nema <- n_ema_tbl %>%
  dplyr::filter(is.na(n_ema) | n_ema < 5 | n_ema > 40) %>%
  dplyr::mutate(reason = "n_ema_out_of_range (must be 5..40)")

# Log
excl_log <- dplyr::bind_rows(excl_careless, excl_nema) %>%
  dplyr::distinct(user_id, n_ema, reason)

# Stampa a console
cat(paste(qa_lines, collapse = "\n"), "\n")

# Salva QA su file
writeLines(qa_lines, con = OUT_QA_TXT)

# ============================
# 7) EXPORT: dataset finale + log esclusioni
# ============================

# Dataset finale
rio::export(dat_clean, OUT_CLEAN_CSV)
saveRDS(dat_clean, OUT_CLEAN_RDS)

# Log esclusioni (uno per riga con motivazione)
excl_log %>%
  arrange(reason, user_id) %>%
  rio::export(OUT_LOG_EXCL)

message(
  "Pulizia completata.\n- Dataset: ",
  OUT_CLEAN_CSV,
  " | ",
  OUT_CLEAN_RDS,
  "\n- Log esclusioni: ",
  OUT_LOG_EXCL,
  "\n- QA: ",
  OUT_QA_TXT
)

# eof ---
