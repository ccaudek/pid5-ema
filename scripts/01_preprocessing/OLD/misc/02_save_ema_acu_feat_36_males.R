# Overview ----------------------------------------------------------------
# Associated project: PID-5, EMA, and acustic features
# Script purpose: Save 2 csv files with complete EMA variables and acustic
#   features, by considering only the 36 male participants.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-09-14
# Last update:
# Status: In progress
# Notes:

# Load necessary libraries ------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(rio)
  library(lubridate)
})

d <- rio::import(
  here::here(
    "data",
    "processed",
    "ema_plus_baseline_exam_tagged.csv"
  )
)

# 1) Tieni solo pre/post e ordina temporalmente
d_pp <- d %>%
  dplyr::filter(
    exam_period %in% c("pre_exam", "post_exam"),
    !is.na(neg_affect_ema)
  ) %>%
  mutate(
    # costruisco un timestamp; se hour mancante, metto mezzogiorno
    dt = as.POSIXct(day) + hours(ifelse(is.na(hour), 12L, hour))
  ) %>%
  arrange(user_id, dt)

# 2) Crea un indice d'esame per soggetto: 1° esame, 2° esame, ...
#    L'idea: per le righe "pre" l'indice è (#post visti finora) + 1,
#    per le righe "post" l'indice è (#post fino a quella riga).
d_pp <- d_pp %>%
  group_by(user_id) %>%
  arrange(dt, .by_group = TRUE) %>%
  mutate(
    cum_post = cumsum(exam_period == "post_exam"),
    exam_idx = if_else(exam_period == "pre_exam", cum_post + 1L, cum_post)
  ) %>%
  ungroup()

# 3) Riassumi per soggetto × esame: media di neg_affect_ema in pre e in post
prepost_diff <- d_pp %>%
  group_by(user_id, exam_idx) %>%
  summarise(
    pre_mean = mean(neg_affect_ema[exam_period == "pre_exam"], na.rm = TRUE),
    post_mean = mean(neg_affect_ema[exam_period == "post_exam"], na.rm = TRUE),
    n_pre = sum(exam_period == "pre_exam"),
    n_post = sum(exam_period == "post_exam"),
    pre_first = suppressWarnings(min(
      dt[exam_period == "pre_exam"],
      na.rm = TRUE
    )),
    post_last = suppressWarnings(max(
      dt[exam_period == "post_exam"],
      na.rm = TRUE
    )),
    .groups = "drop"
  ) %>%
  # Tieni solo le coppie con entrambi i periodi presenti
  dplyr::filter(is.finite(pre_mean), is.finite(post_mean)) %>%
  mutate(
    diff_post_minus_pre = post_mean - pre_mean
  ) %>%
  arrange(user_id, exam_idx)

prepost_diff$diff_post_minus_pre |> hist()

length(unique(prepost_diff$user_id))

bysubj_delta_df <- prepost_diff |>
  group_by(user_id) |>
  summarize(
    delta_neg_aff = mean(diff_post_minus_pre)
  )
# valori negativi grandi: grande differenza tra pre e post, dove pre è più alto di post

# Media DASS per soggetto considerando SOLO il pre_exam
dass_pre_means <- d %>%
  dplyr::filter(exam_period == "pre_exam") %>%
  group_by(user_id) %>%
  summarise(
    dass_stress_mean = mean(dass_stress, na.rm = TRUE),
    dass_depression_mean = mean(dass_depression, na.rm = TRUE),
    dass_anxiety_mean = mean(dass_anxiety, na.rm = TRUE),
    # utili per QC:
    n_pre_rows = n(),
    n_stress_nonmiss = sum(!is.na(dass_stress)),
    n_depr_nonmiss = sum(!is.na(dass_depression)),
    n_anx_nonmiss = sum(!is.na(dass_anxiety)),
    .groups = "drop"
  )

dass_pre_means

delta_df <- left_join(bysubj_delta_df, dass_pre_means, by = "user_id")

temp <- delta_df[, 2:5]
round(cor(temp), 2)

plot(temp$delta_neg_aff, temp$dass_stress_mean)

pid5_mean_df <- d |>
  group_by(user_id) |>
  summarize(
    pid5_negative_affect_baseline_m = mean(
      pid5_negative_affect_baseline,
      na.rm = T
    ),
    pid5_detachment_baseline_m = mean(pid5_detachment_baseline, na.rm = T),
    pid5_antagonism_baseline_m = mean(pid5_antagonism_baseline, na.rm = T),
    pid5_disinhibition_baseline_m = mean(
      pid5_disinhibition_baseline,
      na.rm = T
    ),
    pid5_psychoticism_baseline_m = mean(pid5_psychoticism_baseline, na.rm = T)
  )

delta2_df <- left_join(bysubj_delta_df, pid5_mean_df, by = "user_id")

round(cor(delta2_df[, 2:6], use = "complete.obs"), 2)

# Import Excel sheet pre
acustic_pre <- read_excel(
  here("data", "raw", "acustic_features", "database_acfeat.xlsx"),
  sheet = "Pre"
)

acustic_pre <- acustic_pre %>%
  tidyr::extract(
    Case,
    into = c("day_raw", "user_id"),
    regex = "^\\s*(\\d{2}_\\d{2}_\\d{4})\\s*-\\s*(.+)\\s*$", # "gg_mm_aaaa - resto"
    remove = FALSE
  ) %>%
  mutate(
    day = as.Date(day_raw, format = "%d_%m_%Y")
  ) %>%
  dplyr::select(-day_raw) %>%
  relocate(day, user_id, .after = Case)

# Import Excel sheet post
acustic_post <- read_excel(
  here("data", "raw", "acustic_features", "database_acfeat.xlsx"),
  sheet = "Post"
)

acustic_post <- acustic_post %>%
  tidyr::extract(
    Case,
    into = c("day_raw", "user_id"),
    regex = "^\\s*(\\d{2}_\\d{2}_\\d{4})\\s*-\\s*(.+)\\s*$", # "gg_mm_aaaa - resto"
    remove = FALSE
  ) %>%
  mutate(
    day = as.Date(day_raw, format = "%d_%m_%Y")
  ) %>%
  dplyr::select(-day_raw) %>%
  relocate(day, user_id, .after = Case)

# Selct the ids of the 36 subjects that have the acustic features
selected_ids <- acustic_pre$user_id

# From the complete psychological dataset, select only the 36 subject that
# can be matched with the acustic features data
d1 <- d[d$user_id %in% selected_ids, ]

d2 <- d1 |>
  dplyr::select(
    user_id,
    day,
    pid5_negative_affect_baseline,
    pid5_detachment_baseline,
    pid5_antagonism_baseline,
    pid5_disinhibition_baseline,
    pid5_psychoticism_baseline,
    dass_stress_baseline,
    dass_anxiety_baseline,
    dass_depression_baseline,
    neg_affect_ema,
    dass_stress,
    dass_depression,
    dass_anxiety
  )

d2 <- d2 %>% mutate(day = as.Date(day))

pre_tot_df <- left_join(acustic_pre, d2, by = c("user_id", "day"))
post_tot_df <- left_join(acustic_post, d2, by = c("user_id", "day"))


hz2mel <- function(f) 1127 * log(1 + f / 700)
tri_area <- function(x, y)
  0.5 * abs(x[1] * (y[2] - y[3]) + x[2] * (y[3] - y[1]) + x[3] * (y[1] - y[2]))

pre_tot_df <- pre_tot_df %>%
  # rinomina comodo (adatta ai tuoi nomi esatti: occhio agli spazi!)
  mutate(
    F1a = `F1 mean Hz /a/`,
    F2a = `F2 mean Hz /a/`,
    F1i = `F1 mean Hz  /i/`,
    F2i = `F2 mean Hz  /i/`,
    F1u = `F1 mean Hz   /u/`,
    F2u = `F2 mean Hz   /u/`
  ) %>%
  mutate(
    # Vowel Space Area su scala Mel
    VSA_mel = tri_area(hz2mel(c(F1a, F1i, F1u)), hz2mel(c(F2a, F2i, F2u))),
    # Tilt spettrale (proxy)
    tilt_mfcc1 = `MFCC1 mean`,
    # Dispersione cepstrale (minore = inviluppo più "compresso")
    cep_var = rowSums(
      dplyr::select(cur_data_all(), matches("^MFCC\\d+ std$"))^2,
      na.rm = TRUE
    ),
    # Rumore medio sulle tre vocali (controllo)
    nne_mean = rowMeans(
      dplyr::select(cur_data_all(), matches("^NNE")),
      na.rm = TRUE
    ),
    jitter_mean = rowMeans(
      dplyr::select(cur_data_all(), matches("^Jitter")),
      na.rm = TRUE
    ),
    F0_mean_a = `F0 mean Hz /a/`
  )

post_tot_df <- post_tot_df %>%
  # rinomina comodo (adatta ai tuoi nomi esatti: occhio agli spazi!)
  mutate(
    F1a = `F1 mean Hz /a/`,
    F2a = `F2 mean Hz /a/`,
    F1i = `F1 mean Hz  /i/`,
    F2i = `F2 mean Hz  /i/`,
    F1u = `F1 mean Hz   /u/`,
    F2u = `F2 mean Hz   /u/`
  ) %>%
  mutate(
    # Vowel Space Area su scala Mel
    VSA_mel = tri_area(hz2mel(c(F1a, F1i, F1u)), hz2mel(c(F2a, F2i, F2u))),
    # Tilt spettrale (proxy)
    tilt_mfcc1 = `MFCC1 mean`,
    # Dispersione cepstrale (minore = inviluppo più "compresso")
    cep_var = rowSums(
      dplyr::select(cur_data_all(), matches("^MFCC\\d+ std$"))^2,
      na.rm = TRUE
    ),
    # Rumore medio sulle tre vocali (controllo)
    nne_mean = rowMeans(
      dplyr::select(cur_data_all(), matches("^NNE")),
      na.rm = TRUE
    ),
    jitter_mean = rowMeans(
      dplyr::select(cur_data_all(), matches("^Jitter")),
      na.rm = TRUE
    ),
    F0_mean_a = `F0 mean Hz /a/`
  )

delta <- post_tot_df$cep_var - pre_tot_df$cep_var
hist(delta)

delta_na_ema <- post_tot_df$neg_affect_ema - pre_tot_df$neg_affect_ema

plot(delta_nne, delta_na_ema)

rio::export(
  pre_tot_df,
  "pre_exam.csv"
)
rio::export(
  post_tot_df,
  "post_exam.csv"
)


# eof ---
