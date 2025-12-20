# ==============================================================================
# 01_clean_after_merge_IMPROVED.R
# Data cleaning post-merge EMA + baseline
#
# CHANGES FROM ORIGINAL:
# 1. Added acoustic outlier detection (physiologically implausible values)
# 2. Increased MIN_SD_NA from 0.30 to 0.50 (stricter response quality)
# 3. Sample is all-female (documented, no sex covariate needed)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(lubridate)
  library(readr)
  library(rio)
})

# ==============================================================================
# 0) PARAMETERS
# ==============================================================================

# ---- EMA Quality Thresholds ----
MIN_PRE <- 1L # Minimum beeps in PRE period
MIN_POST <- 1L # Minimum beeps in POST period
MIN_SD_NA <- 0.50 # CHANGED: was 0.30, now stricter (scale 1..7)
MAX_MULT10 <- 0.80 # Maximum fraction of round numbers (0..100 scale)

# ---- Acoustic Outlier Thresholds (for female speakers) ----
# These are standard thresholds from voice research literature

# F0 (fundamental frequency) - female typical range 150-300 Hz
# Using wider bounds to be conservative
F0_MIN <- 100 # Hz - below this is physiologically implausible for females
F0_MAX <- 400 # Hz - above this is physiologically implausible for females

# F0 SD (pitch variability) - typically 10-50 Hz for sustained vowels
F0_SD_MIN <- 1 # Hz - essentially monotone, likely measurement error
F0_SD_MAX <- 100 # Hz - extremely high variability, likely artifact

# Jitter (cycle-to-cycle pitch perturbation) - healthy voice < 1%, acceptable < 2%
JITTER_MAX <- 5 # % - above this suggests recording artifact or pathology

# NNE (Normalized Noise Energy) - typical range -35 to -15 dB
NNE_MIN <- -50 # dB - extremely low, possibly measurement issue
NNE_MAX <- 0 # dB - very noisy signal, unreliable

# F2 (second formant) - varies by vowel
# /a/: ~1000-1400 Hz, /i/: ~2000-2800 Hz, /u/: ~700-1200 Hz for females
F2_MIN <- 500 # Hz - below this is implausible
F2_MAX <- 3500 # Hz - above this is implausible

# F2 SD (formant variability)
F2_SD_MIN <- 5 # Hz - essentially no variation
F2_SD_MAX <- 500 # Hz - extreme variation, likely artifact

# ---- File Paths ----
INPUT_RDS <- here::here("data", "processed", "ema_plus_scales_merged.RDS")
INPUT_CSV <- here::here("data", "processed", "ema_plus_scales_merged.csv")
META_PATH <- here::here("data", "raw", "meta", "all_combined_sex_NEW_1.xlsx")
EXAM_TAGS_PATH <- here::here("data", "raw", "meta", "exam_periods.csv")

# Output
OUT_CLEAN_RDS <- here::here("data", "processed", "ema_plus_scales_cleaned.rds")
OUT_CLEAN_CSV <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
OUT_LOG_EXCL <- here::here("data", "processed", "ema_exclusion_log.csv")
OUT_QA_TXT <- here::here("data", "processed", "ema_cleaning_QA.txt")
OUT_ACOUSTIC_LOG <- here::here("data", "processed", "acoustic_outliers_log.csv")

# ==============================================================================
# 1) LOAD RAW DATA
# ==============================================================================

dat_raw <- if (file.exists(INPUT_RDS)) {
  readRDS(INPUT_RDS)
} else if (file.exists(INPUT_CSV)) {
  suppressMessages(readr::read_csv(INPUT_CSV, show_col_types = FALSE))
} else {
  stop(
    "Cannot find raw data file (neither RDS nor CSV). Update INPUT_RDS/INPUT_CSV."
  )
}

# Normalize column names
dat_raw <- dat_raw %>%
  rename(
    user_id = any_of(c("user_id", "UserID", "id_anon")),
    day = any_of(c("day", "date", "Date", "baseline_date")),
    ema_time = any_of(c("ema_time", "timestamp", "time", "datetime"))
  )

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("DATA CLEANING WITH ACOUSTIC OUTLIER DETECTION\n")
cat(rep("=", 70), "\n", sep = "")
cat(
  "\nRaw data loaded:",
  nrow(dat_raw),
  "rows,",
  n_distinct(dat_raw$user_id),
  "participants\n"
)

# ==============================================================================
# 2) CARELESS RESPONDING LIST (A PRIORI)
# ==============================================================================

careless_ids <- c(
  # "ma_se_2005_11_14_490",
  # "reve20041021036",
  # "di_ma_2005_10_20_756",
  # "pa_sc_2005_09_10_468",
  # "il_re_2006_01_18_645",
  # "so_ma_2003_10_13_804",
  # "lo_ca_2005_05_07_05_437",
  # "va_ma_2005_05_31_567",
  # "no_un_2005_06_29_880",
  # "an_bo_1988_08_24_166",
  # "st_ma_2004_04_21_426",
  # "an_st_2005_10_16_052",
  # "vi_de_2002_12_30_067",
  # "gi_ru_2005_03_08_033",
  # "al_mi_2005_03_05_844",
  # "la_ma_2006_01_31_787",
  # "gi_lo_2004_06_27_237",
  # "ch_bi_2001_01_28_407",
  # "al_pe_2001_04_20_079",
  # "le_de_2003_09_05_067",
  # "fe_gr_2002_02_19_434",
  # "ma_ba_2002_09_09_052",
  # "ca_gi_2003_09_16_737",
  # "an_to_2003_08_06_114",
  # "al_se_2003_07_28_277",
  # "ja_tr_2002_10_06_487",
  # "el_ci_2002_02_15_057",
  # "se_ti_2000_03_04_975",
  # "co_ga_2003_10_29_614",
  # "al_ba_2003_18_07_905",
  # "bi_ro_2003_09_07_934",
  # "an_va_2004_04_08_527",
  # "ev_cr_2003_01_27_573"
)

# ==============================================================================
# 3) COMPUTE n_ema
# ==============================================================================

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

# ==============================================================================
# 4) METADATA (course/sex) + EXAM PERIOD TAGS
# ==============================================================================

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

# Exam period windows by course
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

# ==============================================================================
# 5) EXCLUSIONS - STEP 1: Careless + n_ema
# ==============================================================================

cat("\n--- Step 1: Careless responding + n_ema range ---\n")

# (A) Careless
excl_careless <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(user_id %in% careless_ids) %>%
  dplyr::left_join(n_ema_tbl, by = "user_id") %>%
  dplyr::mutate(reason = "careless_responding (a priori)")

# (B) n_ema out of range
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

cat("  Excluded (careless):", n_distinct(excl_careless$user_id), "\n")
cat("  Excluded (n_ema):", n_distinct(excl_nema$user_id), "\n")
cat("  Remaining:", n_distinct(dat_step1$user_id), "participants\n")

# ==============================================================================
# 6) EXCLUSIONS - STEP 2: ACOUSTIC OUTLIER DETECTION (NEW)
# ==============================================================================

cat("\n--- Step 2: Acoustic outlier detection ---\n")

# Initialize acoustic outlier counters (for QA safety)
n_acoustic_outliers <- 0L
n_participants_with_outliers <- 0L

# Identify acoustic variable columns (for all vowels)
# Expected naming: f0_mean_a, f0_std_a, jitter_a, nne_a, f2_mean_a, f2_std_a, etc.

acoustic_vars <- list(
  f0_mean = c("f0_mean_a", "f0_mean_i", "f0_mean_u"),
  f0_std = c("f0_std_a", "f0_std_i", "f0_std_u"),
  jitter = c("jitter_a", "jitter_i", "jitter_u"),
  nne = c("nne_a", "nne_i", "nne_u"),
  f2_mean = c("f2_mean_a", "f2_mean_i", "f2_mean_u"),
  f2_std = c("f2_std_a", "f2_std_i", "f2_std_u")
)

# Check which acoustic variables exist in the data
available_acoustic <- unlist(acoustic_vars)
available_acoustic <- available_acoustic[
  available_acoustic %in% names(dat_step1)
]

if (length(available_acoustic) > 0) {
  cat(
    "  Found acoustic variables:",
    paste(available_acoustic, collapse = ", "),
    "\n"
  )

  # Create a function to flag outliers for each observation
  flag_acoustic_outliers <- function(df) {
    df <- df %>%
      mutate(
        acoustic_outlier = FALSE,
        acoustic_outlier_reason = NA_character_
      )

    # F0 mean outliers
    for (v in intersect(acoustic_vars$f0_mean, names(df))) {
      outlier_flag <- !is.na(df[[v]]) & (df[[v]] < F0_MIN | df[[v]] > F0_MAX)
      df$acoustic_outlier <- df$acoustic_outlier | outlier_flag
      df$acoustic_outlier_reason <- ifelse(
        outlier_flag & is.na(df$acoustic_outlier_reason),
        paste0(v, "_out_of_range"),
        df$acoustic_outlier_reason
      )
    }

    # F0 SD outliers
    for (v in intersect(acoustic_vars$f0_std, names(df))) {
      outlier_flag <- !is.na(df[[v]]) &
        (df[[v]] < F0_SD_MIN | df[[v]] > F0_SD_MAX)
      df$acoustic_outlier <- df$acoustic_outlier | outlier_flag
      df$acoustic_outlier_reason <- ifelse(
        outlier_flag & is.na(df$acoustic_outlier_reason),
        paste0(v, "_out_of_range"),
        df$acoustic_outlier_reason
      )
    }

    # Jitter outliers (assuming jitter is in %)
    for (v in intersect(acoustic_vars$jitter, names(df))) {
      outlier_flag <- !is.na(df[[v]]) & (df[[v]] > JITTER_MAX)
      df$acoustic_outlier <- df$acoustic_outlier | outlier_flag
      df$acoustic_outlier_reason <- ifelse(
        outlier_flag & is.na(df$acoustic_outlier_reason),
        paste0(v, "_too_high"),
        df$acoustic_outlier_reason
      )
    }

    # NNE outliers
    for (v in intersect(acoustic_vars$nne, names(df))) {
      outlier_flag <- !is.na(df[[v]]) & (df[[v]] < NNE_MIN | df[[v]] > NNE_MAX)
      df$acoustic_outlier <- df$acoustic_outlier | outlier_flag
      df$acoustic_outlier_reason <- ifelse(
        outlier_flag & is.na(df$acoustic_outlier_reason),
        paste0(v, "_out_of_range"),
        df$acoustic_outlier_reason
      )
    }

    # F2 mean outliers
    for (v in intersect(acoustic_vars$f2_mean, names(df))) {
      outlier_flag <- !is.na(df[[v]]) & (df[[v]] < F2_MIN | df[[v]] > F2_MAX)
      df$acoustic_outlier <- df$acoustic_outlier | outlier_flag
      df$acoustic_outlier_reason <- ifelse(
        outlier_flag & is.na(df$acoustic_outlier_reason),
        paste0(v, "_out_of_range"),
        df$acoustic_outlier_reason
      )
    }

    # F2 SD outliers
    for (v in intersect(acoustic_vars$f2_std, names(df))) {
      outlier_flag <- !is.na(df[[v]]) &
        (df[[v]] < F2_SD_MIN | df[[v]] > F2_SD_MAX)
      df$acoustic_outlier <- df$acoustic_outlier | outlier_flag
      df$acoustic_outlier_reason <- ifelse(
        outlier_flag & is.na(df$acoustic_outlier_reason),
        paste0(v, "_out_of_range"),
        df$acoustic_outlier_reason
      )
    }

    return(df)
  }

  # Apply outlier detection
  dat_step1 <- flag_acoustic_outliers(dat_step1)

  # Summary of acoustic outliers
  n_acoustic_outliers <- sum(dat_step1$acoustic_outlier, na.rm = TRUE)
  n_participants_with_outliers <- dat_step1 %>%
    filter(acoustic_outlier) %>%
    distinct(user_id) %>%
    nrow()

  cat("  Observations flagged as acoustic outliers:", n_acoustic_outliers, "\n")
  cat(
    "  Participants with at least one outlier:",
    n_participants_with_outliers,
    "\n"
  )

  # Log acoustic outliers (for inspection)
  acoustic_outlier_log <- dat_step1 %>%
    filter(acoustic_outlier) %>%
    select(
      user_id,
      day,
      exam_period,
      acoustic_outlier_reason,
      any_of(available_acoustic)
    )

  if (nrow(acoustic_outlier_log) > 0) {
    rio::export(acoustic_outlier_log, OUT_ACOUSTIC_LOG)
    cat("  Acoustic outlier log saved to:", OUT_ACOUSTIC_LOG, "\n")
  }

  # DECISION: Set acoustic outlier values to NA (rather than excluding entire participant)
  # This preserves the participant but removes unreliable measurements
  cat(
    "\n  Strategy: Setting outlier acoustic values to NA (preserving participants)\n"
  )

  for (v in available_acoustic) {
    if (v %in% names(dat_step1)) {
      # Determine which threshold to apply
      if (grepl("f0_mean", v)) {
        dat_step1[[v]] <- ifelse(
          dat_step1[[v]] < F0_MIN | dat_step1[[v]] > F0_MAX,
          NA_real_,
          dat_step1[[v]]
        )
      } else if (grepl("f0_std", v)) {
        dat_step1[[v]] <- ifelse(
          dat_step1[[v]] < F0_SD_MIN | dat_step1[[v]] > F0_SD_MAX,
          NA_real_,
          dat_step1[[v]]
        )
      } else if (grepl("jitter", v)) {
        dat_step1[[v]] <- ifelse(
          dat_step1[[v]] > JITTER_MAX,
          NA_real_,
          dat_step1[[v]]
        )
      } else if (grepl("nne", v)) {
        dat_step1[[v]] <- ifelse(
          dat_step1[[v]] < NNE_MIN | dat_step1[[v]] > NNE_MAX,
          NA_real_,
          dat_step1[[v]]
        )
      } else if (grepl("f2_mean", v)) {
        dat_step1[[v]] <- ifelse(
          dat_step1[[v]] < F2_MIN | dat_step1[[v]] > F2_MAX,
          NA_real_,
          dat_step1[[v]]
        )
      } else if (grepl("f2_std", v)) {
        dat_step1[[v]] <- ifelse(
          dat_step1[[v]] < F2_SD_MIN | dat_step1[[v]] > F2_SD_MAX,
          NA_real_,
          dat_step1[[v]]
        )
      }
    }
  }

  # Report how many values were set to NA
  cat("  Acoustic values set to NA due to outliers:\n")
  for (v in available_acoustic) {
    if (v %in% names(dat_step1)) {
      n_na <- sum(is.na(dat_step1[[v]]))
      cat("    ", v, ":", n_na, "NA values\n")
    }
  }
} else {
  cat(
    "  No acoustic variables found in dataset. Skipping acoustic outlier detection.\n"
  )
  dat_step1$acoustic_outlier <- FALSE
  dat_step1$acoustic_outlier_reason <- NA_character_
}

# ==============================================================================
# 7) COMPUTE QUALITY METRICS BY PARTICIPANT
# ==============================================================================

cat("\n--- Step 3: EMA response quality metrics ---\n")

# Helper to map 0-100 -> 1..7
as_item_to_1_7 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(as.integer(x))
  if (all(x[is.finite(x)] %in% 1:7)) return(as.integer(x))
  xmin <- suppressWarnings(min(x, na.rm = TRUE))
  xmax <- suppressWarnings(max(x, na.rm = TRUE))
  if (is.finite(xmin) && is.finite(xmax) && xmin >= 0 && xmax <= 100) {
    brk <- seq(0, 100, length.out = 8)
    return(as.integer(findInterval(
      x,
      brk,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )))
  }
  y <- 1 + 6 * (x - xmin) / (xmax - xmin)
  as.integer(round(pmax(1, pmin(7, y))))
}

# Create long format for quality metrics
items_long <- dat_step1 %>%
  transmute(
    user_id,
    per = dplyr::case_when(
      exam_period == "baseline" ~ 1L,
      exam_period == "pre_exam" ~ 2L,
      exam_period == "post_exam" ~ 3L,
      TRUE ~ 1L
    ),
    happy_100 = happy,
    sad_100 = sad,
    satis_100 = satisfied,
    angry_100 = angry,
    happy = as_item_to_1_7(happy),
    sad = as_item_to_1_7(sad),
    satis = as_item_to_1_7(satisfied),
    angry = as_item_to_1_7(angry)
  ) %>%
  tidyr::drop_na(happy, sad, satis, angry) %>%
  mutate(
    na_item1_7 = (8L - happy + sad + (8L - satis) + angry) / 4
  )

# Quality metrics by participant
quality_by_user <- items_long %>%
  dplyr::group_by(user_id) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_per_1 = sum(per == 1L, na.rm = TRUE),
    n_per_2 = sum(per == 2L, na.rm = TRUE),
    n_per_3 = sum(per == 3L, na.rm = TRUE),
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
    sd_happy_100 = sd(happy_100, na.rm = TRUE),
    sd_sad_100 = sd(sad_100, na.rm = TRUE),
    sd_satis_100 = sd(satis_100, na.rm = TRUE),
    sd_angry_100 = sd(angry_100, na.rm = TRUE),
    pct_mult10 = {
      v <- c(happy_100, sad_100, satis_100, angry_100)
      mean(v %% 10 == 0, na.rm = TRUE)
    },
    .groups = "drop"
  )

# ==============================================================================
# 8) EXCLUSIONS - STEP 3: Quality-based exclusions
# ==============================================================================

excl_quality <- quality_by_user %>%
  transmute(
    user_id,
    reason_coverage = ifelse(
      n_per_2 < MIN_PRE | n_per_3 < MIN_POST,
      sprintf("low_coverage(pre<%d or post<%d)", MIN_PRE, MIN_POST),
      NA_character_
    ),
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
  pivot_longer(-user_id, names_to = "rname", values_to = "reason") %>%
  filter(!is.na(reason)) %>%
  select(user_id, reason)

cat(
  "  Excluded (low coverage):",
  sum(grepl("^low_coverage", excl_quality$reason)),
  "\n"
)
cat(
  "  Excluded (low SD, <",
  MIN_SD_NA,
  "):",
  sum(grepl("^low_within_sd", excl_quality$reason)),
  "\n"
)
cat(
  "  Excluded (round digits):",
  sum(grepl("^too_many_round_digits", excl_quality$reason)),
  "\n"
)

# Combine all exclusion logs
excl_log <- bind_rows(excl_log_step1, excl_quality) %>%
  distinct(user_id, .keep_all = TRUE)

keep_ids <- dat_step1 %>%
  distinct(user_id) %>%
  filter(!(user_id %in% excl_log$user_id)) %>%
  pull(user_id)

dat_clean <- dat_step1 %>%
  filter(user_id %in% keep_ids)

# ==============================================================================
# 9) QA SUMMARY
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("QUALITY ASSURANCE SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")

qa_lines <- c()

n_subj_raw <- n_distinct(dat_raw$user_id)
n_subj_step1 <- n_distinct(dat_step1$user_id)
n_subj_final <- n_distinct(dat_clean$user_id)

qa_lines <- c(
  qa_lines,
  sprintf("Total participants (raw): %d", n_subj_raw),
  sprintf(
    "Excluded step 1 (careless + n_ema out of [5,40]): %d",
    n_distinct(excl_log_step1$user_id)
  ),
  sprintf("  - Careless: %d", n_distinct(excl_careless$user_id)),
  sprintf("  - n_ema out of [5,40]: %d", n_distinct(excl_nema$user_id)),
  sprintf("Participants after step 1: %d", n_subj_step1),
  "",
  sprintf(
    "Acoustic outliers: %d observations flagged (values set to NA)",
    n_acoustic_outliers
  ),
  "",
  sprintf(
    "Excluded step 2 (quality rules): %d",
    n_distinct(setdiff(excl_log$user_id, excl_log_step1$user_id))
  ),
  sprintf(
    "  - Low coverage (PRE<%d or POST<%d): %d",
    MIN_PRE,
    MIN_POST,
    sum(grepl("^low_coverage", excl_quality$reason))
  ),
  sprintf(
    "  - SD(NA) < %.2f: %d",
    MIN_SD_NA,
    sum(grepl("^low_within_sd", excl_quality$reason))
  ),
  sprintf(
    "  - Round digits > %.2f: %d",
    MAX_MULT10,
    sum(grepl("^too_many_round_digits", excl_quality$reason))
  ),
  "",
  sprintf("FINAL SAMPLE: %d participants", n_subj_final),
  ""
)

# n_ema distribution
nema_keep <- n_ema_tbl %>% filter(user_id %in% dat_clean$user_id)
qa_lines <- c(
  qa_lines,
  sprintf(
    "n_ema (included) — min: %d, median: %d, mean: %.2f, max: %d",
    min(nema_keep$n_ema, na.rm = TRUE),
    median(nema_keep$n_ema, na.rm = TRUE) %>% as.integer(),
    mean(nema_keep$n_ema, na.rm = TRUE),
    max(nema_keep$n_ema, na.rm = TRUE)
  )
)

# Sex distribution (should be all female)
if ("sex" %in% names(dat_clean)) {
  gen_tab <- dat_clean %>% distinct(user_id, sex) %>% count(sex, sort = TRUE)
  qa_lines <- c(qa_lines, "", "Sex distribution (should be all female):")
  qa_lines <- c(qa_lines, paste0("  - ", gen_tab$sex, ": ", gen_tab$n))
}

# Exam period distribution
if ("exam_period" %in% names(dat_clean)) {
  ex_tab <- dat_clean %>% count(exam_period)
  qa_lines <- c(qa_lines, "", "Exam period distribution (observations):")
  qa_lines <- c(qa_lines, paste0("  - ", ex_tab$exam_period, ": ", ex_tab$n))
}

# Print and save QA
cat(paste(qa_lines, collapse = "\n"), "\n")
writeLines(qa_lines, con = OUT_QA_TXT)

# ==============================================================================
# 10) EXPORT
# ==============================================================================

cat("\n--- Saving outputs ---\n")

# Remove temporary columns before export
dat_clean <- dat_clean %>%
  select(-any_of(c("acoustic_outlier", "acoustic_outlier_reason", ".is_ema")))

rio::export(dat_clean, OUT_CLEAN_CSV)
saveRDS(dat_clean, OUT_CLEAN_RDS)

excl_log %>%
  arrange(user_id) %>%
  rio::export(OUT_LOG_EXCL)

cat("  Dataset:", OUT_CLEAN_CSV, "\n")
cat("  Dataset:", OUT_CLEAN_RDS, "\n")
cat("  Exclusion log:", OUT_LOG_EXCL, "\n")
cat("  QA report:", OUT_QA_TXT, "\n")
if (exists("acoustic_outlier_log") && nrow(acoustic_outlier_log) > 0) {
  cat("  Acoustic outlier log:", OUT_ACOUSTIC_LOG, "\n")
}

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("CLEANING COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nParameters used:\n")
cat("  MIN_PRE =", MIN_PRE, "\n")
cat("  MIN_POST =", MIN_POST, "\n")
cat("  MIN_SD_NA =", MIN_SD_NA, "(INCREASED from 0.30)\n")
cat("  MAX_MULT10 =", MAX_MULT10, "\n")
cat("\nAcoustic thresholds (for female speakers):\n")
cat("  F0: ", F0_MIN, "-", F0_MAX, " Hz\n")
cat("  F0 SD: ", F0_SD_MIN, "-", F0_SD_MAX, " Hz\n")
cat("  Jitter: <", JITTER_MAX, "%\n")
cat("  NNE: ", NNE_MIN, "-", NNE_MAX, " dB\n")
cat("  F2: ", F2_MIN, "-", F2_MAX, " Hz\n")
cat("  F2 SD: ", F2_SD_MIN, "-", F2_SD_MAX, " Hz\n")

# eof ---
