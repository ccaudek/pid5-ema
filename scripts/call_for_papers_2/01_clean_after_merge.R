# 02_clean_after_merge_UPDATED.R
# Data cleaning post-merge EMA + baseline (versione con criteri qualità soggetto)
# - Esclusione careless responding (lista a priori)
# - Filtro numerosità EMA: 5 <= n_ema <= 40
# - Import metadati (genere/corso) + tag esami (exam_period)
# - Nuovi criteri soggetto (pre-specificati, NON basati su outcome):
#     * copertura minima per periodo (>= min_pre, >= min_post)
#     * variabilità intra-soggetto minima su NA (scala 1..7)
#     * preferenza eccessiva per cifre tonde (multipli di 10 su 0..100)
# - QA esteso + export dataset finale e log esclusioni

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(lubridate)
  library(readr)
  library(rio)
})

# ============================
# 0) PARAMETRI / PERCORSI
# ============================

# ---- Parametri soglie (modifica liberamente) ----
MIN_PRE <- 1L # beep minimi in PRE per soggetto (consigliato 2 o 3)
MIN_POST <- 1L # beep minimi in POST per soggetto (consigliato 2 o 3)
MIN_SD_NA <- 0.30 # variabilità minima su NA (scala 1..7)
MAX_MULT10 <- 0.80 # soglia massima frazione cifre tonde su 0..100

# File "raw totale" già salvato dallo script precedente (RDS/CSV)
INPUT_RDS <- here::here("data", "processed", "ema_plus_scales_merged.RDS")
INPUT_CSV <- here::here("data", "processed", "ema_plus_scales_merged.csv")

# Metadati (almeno: user_id, course, sex) — adatta il percorso/nome file
META_PATH <- here::here("data", "raw", "meta", "all_combined_sex_NEW_1.xlsx")

# Tag esami per corsi (qui implementati come finestre pre/post fisse per corso)
# Se hai un file esterno, puoi ignorare queste costanti e leggere EXAM_TAGS_PATH
EXAM_TAGS_PATH <- here::here("data", "raw", "meta", "exam_periods.csv")

# Output
OUT_CLEAN_RDS <- here::here(
  "data",
  "processed",
  "ema_plus_scales_cleaned_NO_PREPOST_FILTER.rds"
)

OUT_CLEAN_CSV <- here::here(
  "data",
  "processed",
  "ema_plus_scales_cleaned_NO_PREPOST_FILTER.csv"
)

OUT_LOG_EXCL <- here::here(
  "data",
  "processed",
  "ema_exclusion_log_NO_PREPOST_FILTER.csv"
)

OUT_QA_TXT <- here::here(
  "data",
  "processed",
  "ema_cleaning_QA_NO_PREPOST_FILTER.txt"
)

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
    day = any_of(c("day", "date", "Date", "baseline_date")),
    ema_time = any_of(c("ema_time", "timestamp", "time", "datetime"))
  )

# ============================
# 2) LISTA CARELESS A PRIORI
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

if ("n_ema" %in% names(dat_raw)) {
  n_ema_tbl <- dat_raw %>% dplyr::distinct(user_id, n_ema)
} else {
  ind_candidates <- c("ema_time", "day", "beep", "ema_id", "ema_wave")
  present_inds <- intersect(ind_candidates, names(dat_raw))
  flag_mat <- NULL
  if (length(present_inds) > 0) {
    flag_mat <- sapply(present_inds, function(nm) !is.na(dat_raw[[nm]]))
    if (!is.matrix(flag_mat)) flag_mat <- matrix(flag_mat, ncol = 1)
  }
  source_flag <- if ("source" %in% names(dat_raw))
    grepl("ema", dat_raw[["source"]], ignore.case = TRUE) else
    rep(FALSE, nrow(dat_raw))
  is_ema <- if (is.null(flag_mat)) source_flag else
    (rowSums(flag_mat) > 0) | source_flag
  dat_raw$.is_ema <- is_ema

  n_ema_tbl <- dat_raw %>%
    dplyr::group_by(user_id) %>%
    dplyr::summarise(
      n_ema = as.integer(sum(.is_ema, na.rm = TRUE)),
      .groups = "drop"
    )
}
if (!"n_ema" %in% names(dat_raw)) {
  dat_raw <- dat_raw %>% dplyr::left_join(n_ema_tbl, by = "user_id")
}

# ============================
# 4) METADATI (course/sex) + TAG ESAMI (exam_period)
# ============================

meta_df <- suppressMessages(rio::import(META_PATH)) %>%
  dplyr::rename(
    user_id = any_of(c("subj_code", "user_id", "id_anon")),
    date_meta = any_of(c("date", "Date", "data")),
    course = any_of(c("course", "Corso", "insegnamento")),
    sex = any_of(c("sex", "Sex", "genere", "gender"))
  ) %>%
  dplyr::mutate(
    date_meta = suppressWarnings(as.Date(date_meta, format = "%d_%m_%Y")),
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

dat_final <- dat_raw %>%
  dplyr::left_join(meta_df, by = "user_id") %>%
  dplyr::mutate(day = as.Date(day))

# Finestre PRE/POST per corso (adatta se necessario)
psico_pre <- as.Date(c("2025-04-14", "2025-05-21"))
psico_post <- as.Date(c("2025-04-15", "2025-05-22"))
test_pre <- as.Date(c("2025-04-14", "2025-05-25"))
test_post <- as.Date(c("2025-04-15", "2025-05-26"))
interv_pre <- as.Date("2025-05-12")
interv_post <- as.Date("2025-05-13")

dat_final <- dat_final %>%
  dplyr::mutate(
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

# ============================
# 5) ESCLUSIONI — step 1 (careless + n_ema)
# ============================

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

excl_log_step1 <- dplyr::bind_rows(excl_careless, excl_nema) %>%
  dplyr::distinct(user_id, n_ema, reason)

keep_step1 <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(!(user_id %in% excl_log_step1$user_id)) %>%
  dplyr::pull(user_id)

dat_step1 <- dat_final %>% dplyr::filter(user_id %in% keep_step1)

# ============================
# 6) COSTRUZIONE METRICHE QUALITÀ per soggetto
# ============================

# Helper per mappare 0–100 -> 1..7 (per NA medio su scala item)
as_item_to_1_7 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(as.integer(x))
  if (all(x[is.finite(x)] %in% 1:7)) return(as.integer(x))
  xmin <- suppressWarnings(min(x, na.rm = TRUE))
  xmax <- suppressWarnings(max(x, na.rm = TRUE))
  if (is.finite(xmin) && is.finite(xmax) && xmin >= 0 && xmax <= 100) {
    brk <- seq(0, 100, length.out = 8) # 7 classi
    return(as.integer(findInterval(
      x,
      brk,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )))
  }
  y <- 1 + 6 * (x - xmin) / (xmax - xmin) # fallback lineare
  as.integer(round(pmax(1, pmin(7, y))))
}

# Long degli item EMA con:
# - per: 1=baseline, 2=pre, 3=post
# - item in 0–100 e in 1..7
# - NA medio orientato a "negative affect" (happy/satisfied invertiti)
items_long <- dat_step1 %>%
  transmute(
    user_id,
    per = dplyr::case_when(
      exam_period == "baseline" ~ 1L,
      exam_period == "pre_exam" ~ 2L,
      exam_period == "post_exam" ~ 3L,
      TRUE ~ 1L
    ),
    # 0–100
    happy_100 = happy,
    sad_100 = sad,
    satis_100 = satisfied,
    angry_100 = angry,
    # 1..7 (per PPC/qualità)
    happy = as_item_to_1_7(happy),
    sad = as_item_to_1_7(sad),
    satis = as_item_to_1_7(satisfied),
    angry = as_item_to_1_7(angry)
  ) %>%
  tidyr::drop_na(happy, sad, satis, angry) %>%
  mutate(
    # Allinea a NA: inverti happy e satisfied
    na_item1_7 = (8L - happy + sad + (8L - satis) + angry) / 4
  )

# Metriche per soggetto (nessun join: i conteggi sono calcolati qui)
quality_by_user <- items_long %>%
  dplyr::group_by(user_id) %>%
  dplyr::summarise(
    # quantità di dati
    n_total = dplyr::n(),
    n_per_1 = sum(per == 1L, na.rm = TRUE),
    n_per_2 = sum(per == 2L, na.rm = TRUE),
    n_per_3 = sum(per == 3L, na.rm = TRUE),

    # variabilità/qualità su scala 1..7
    sd_within_1_7 = sd(na_item1_7, na.rm = TRUE),

    pct_same_consec_1_7 = {
      v <- na_item1_7
      if (length(v) < 2) NA_real_ else mean(diff(v) == 0, na.rm = TRUE)
    },

    max_run_len_1_7 = {
      v <- na_item1_7
      if (length(v) == 0) NA_integer_ else {
        r <- rle(v)
        max(r$lengths, na.rm = TRUE)
      }
    },

    pct_extremes_1_7 = mean(na_item1_7 %in% c(1, 7), na.rm = TRUE),

    entropy_1_7 = {
      v <- na_item1_7
      p <- prop.table(table(v))
      -sum(p * log(p + 1e-12))
    },

    # varianza intra-soggetto sulle scale 0–100 (per ciascun item)
    sd_happy_100 = sd(happy_100, na.rm = TRUE),
    sd_sad_100 = sd(sad_100, na.rm = TRUE),
    sd_satis_100 = sd(satis_100, na.rm = TRUE),
    sd_angry_100 = sd(angry_100, na.rm = TRUE),

    # preferenza per cifre “tonde” (0–100)
    pct_mult10 = {
      v <- c(happy_100, sad_100, satis_100, angry_100)
      mean(v %% 10 == 0, na.rm = TRUE)
    },

    .groups = "drop"
  )

# ============================
# 7) ESCLUSIONI — step 2 (criteri qualità pre-specificati)
# ============================

excl_quality <- quality_by_user %>%
  transmute(
    user_id,
    reason_low_sd = ifelse(
      sd_within_1_7 < MIN_SD_NA,
      sprintf("low_within_sd(<%.2f)", MIN_SD_NA),
      NA_character_
    ),
    reason_mult10 = ifelse(
      pct_mult10 > MAX_MULT10,
      sprintf("too_many_round_digits(>%.2f)", MAX_MULT10),
      NA_character_
    )
  ) %>%
  tidyr::pivot_longer(
    -user_id,
    names_to = "rname",
    values_to = "reason"
  ) %>%
  dplyr::filter(!is.na(reason)) %>%
  dplyr::select(user_id, reason)

# Unifica log step1 + step2
excl_log <- bind_rows(excl_log_step1, excl_quality) %>%
  distinct(user_id, .keep_all = TRUE)

keep_ids <- dat_step1 %>%
  distinct(user_id) %>%
  filter(!(user_id %in% excl_log$user_id)) %>%
  pull(user_id)

dat_clean <- dat_step1 %>%
  filter(user_id %in% keep_ids)

# ============================
# 8) QA ESSENZIALE
# ============================

qa_lines <- c()

n_subj_raw <- n_distinct(dat_raw$user_id)
n_subj_step1 <- n_distinct(dat_step1$user_id)
n_subj_final <- n_distinct(dat_clean$user_id)

qa_lines <- c(
  qa_lines,
  sprintf("Totale soggetti (raw): %d", n_subj_raw),
  sprintf(
    "Esclusi step1 (careless + n_ema out of [5,40]): %d",
    n_distinct(excl_log_step1$user_id)
  ),
  sprintf("  - Careless: %d", n_distinct(excl_careless$user_id)),
  sprintf("  - n_ema fuori [5,40]: %d", n_distinct(excl_nema$user_id)),
  sprintf("Soggetti dopo step1: %d", n_subj_step1),
  sprintf(
    "Esclusi step2 (quality rules): %d",
    n_distinct(setdiff(excl_log$user_id, excl_log_step1$user_id))
  ),
  sprintf(
    "  - Copertura bassa (PRE<%d o POST<%d): %d",
    MIN_PRE,
    MIN_POST,
    sum(grepl("^low_coverage", excl_quality$reason))
  ),
  sprintf(
    "  - SD(NA)<%.2f: %d",
    MIN_SD_NA,
    sum(grepl("^low_within_sd", excl_quality$reason))
  ),
  sprintf(
    "  - Cifre tonde >%.2f: %d",
    MAX_MULT10,
    sum(grepl("^too_many_round_digits", excl_quality$reason))
  ),
  sprintf("Soggetti inclusi (finale): %d", n_subj_final),
  ""
)

# Distribuzione n_ema (inclusi finali)
nema_keep <- n_ema_tbl %>% filter(user_id %in% dat_clean$user_id)
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

# Genere
if ("sex" %in% names(dat_clean)) {
  gen_tab <- dat_clean %>% distinct(user_id, sex) %>% count(sex, sort = TRUE)
  qa_lines <- c(qa_lines, "", "Distribuzione genere (soggetti unici):")
  qa_lines <- c(qa_lines, paste0(" - ", gen_tab$sex, ": ", gen_tab$n))
}

# Exam period (conteggio beep per periodo sui soggetti inclusi)
if ("exam_period" %in% names(dat_clean)) {
  ex_tab <- dat_clean %>% count(exam_period)
  qa_lines <- c(qa_lines, "", "Distribuzione exam_period (righe/beep):")
  qa_lines <- c(qa_lines, paste0(" - ", ex_tab$exam_period, ": ", ex_tab$n))
}

cat(paste(qa_lines, collapse = "\n"), "\n")
writeLines(qa_lines, con = OUT_QA_TXT)

# Remove baseline-only participants
dat2_clean <- dat_clean |>
  dplyr::filter(
    course != "Clinica" &
      sex == "Femmina"
  )

length(unique(dat2_clean$user_id))

# ============================
# 9) EXPORT
# ============================

rio::export(dat2_clean, OUT_CLEAN_CSV)
saveRDS(dat_clean, OUT_CLEAN_RDS)

excl_log %>%
  arrange(user_id) %>%
  rio::export(OUT_LOG_EXCL)

message(
  "Pulizia completata.\n- Dataset: ",
  OUT_CLEAN_CSV,
  " | ",
  OUT_CLEAN_RDS,
  "\n- Log esclusioni: ",
  OUT_LOG_EXCL,
  "\n- QA: ",
  OUT_QA_TXT,
  "\nParametri: MIN_PRE=",
  MIN_PRE,
  " MIN_POST=",
  MIN_POST,
  " MIN_SD_NA=",
  MIN_SD_NA,
  " MAX_MULT10=",
  MAX_MULT10
)

# eof ---
