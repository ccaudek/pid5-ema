# ==============================================================================
# Bayesian Analysis: Voice Acoustics, Stress, and Personality Pathology
# Context-dependent expression of PID-5 traits via passive sensing
# ==============================================================================

# Load required packages
library(here)
library(tidyverse)
library(readxl)
library(brms)
library(cmdstanr)
library(bayestestR)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(ggdist)
library(lubridate)
library(stringr)
library(missRanger)

# Set cmdstanr as backend for brms
options(brms.backend = "cmdstanr")

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

# Load data from three timepoints
baseline <- read_excel(
  here::here(
    "data",
    "raw",
    "acustic_features",
    "datiacustici",
    "AUDIO.xlsx"
  ),
  sheet = "BASELINE"
)

pre <- read_excel(
  here::here(
    "data",
    "raw",
    "acustic_features",
    "datiacustici",
    "AUDIO.xlsx"
  ),
  sheet = "PRE"
)

post <- read_excel(
  here::here(
    "data",
    "raw",
    "acustic_features",
    "datiacustici",
    "AUDIO.xlsx"
  ),
  sheet = "POST"
)

# Add timepoint indicator
baseline$timepoint <- "baseline"
pre$timepoint <- "pre"
post$timepoint <- "post"

# Combine all data
df_wide <- bind_rows(baseline, pre, post)

# Clean column names (remove extra spaces)
names(df_wide) <- str_trim(names(df_wide))

# Extract key acoustic variables (we'll focus on /a/ vowel for primary analysis)
# F0 mean - pitch (arousal/anxiety)
# F0 std - pitch variability (emotional instability)
# Jitter - voice quality (tension)
# NNE - normalized noise energy (vocal stress)

# Select relevant variables
df_clean <- df_wide %>%
  dplyr::select(
    ID,
    Data,
    timepoint,
    # Acoustic measures (vowel /a/)
    f0_mean_a = `F0 mean Hz /a/`,
    f0_std_a = `F0 std Hz /a/`,
    jitter_a = `Jitter /a/`,
    nne_a = `NNE /a/`,
    f2_mean_a = `F2 mean Hz /a/`,
    f2_std_a = `F2 Std Hz /a/`,

    # Acoustic measures (vowel /i/)
    f0_mean_i = `F0 mean Hz /i/`,
    f0_std_i = `F0 std Hz /i/`,
    jitter_i = `Jitter /i/`,
    nne_i = `NNE /i/`,
    f2_mean_i = `F2 mean Hz /i/`,
    f2_std_i = `F2 Std Hz /i/`,

    # Acoustic measures (vowel /u/)
    f0_mean_u = `F0 mean Hz /u/`,
    f0_std_u = `F0 std Hz /u/`,
    jitter_u = `Jitter /u/`,
    nne_u = `NNE /u/`,
    f2_mean_u = `F2 mean Hz /u/`,
    f2_std_u = `F2 Std Hz /u/`,

    # PID-5 traits
    pid5_na = pid5_negative_affectivity,
    pid5_det = pid5_detachment,
    pid5_ant = pid5_antagonism,
    pid5_dis = pid5_disinhibition,
    pid5_psy = pid5_psychoticism
  ) %>%
  # Convert timepoint to factor with meaningful order
  mutate(timepoint = factor(timepoint, levels = c("baseline", "pre", "post")))


# Add correct PID-5 EMA and baseline measures ----------------------------------

data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
d <- rio::import(data_path)


vars <- c(
  "user_id",
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism",
  "domain_negative_affect_baseline",
  "domain_detachment_baseline",
  "domain_antagonism_baseline",
  "domain_disinhibition_baseline",
  "domain_psychoticism_baseline",
  "tripm_1",
  "tripm_2",
  "tripm_3",
  "tripm_4",
  "tripm_4_rev",
  "tripm_sum",
  "tripm_boldness",
  "tripm_meanness",
  "boldness_baseline",
  "meanness_baseline",
  "domain_disinhibition_baseline",
  "esi_bf_baseline",
  "exam_period",
  "day",
  "sex"
)

d1 <- d |>
  dplyr::select(all_of(vars)) |>
  dplyr::rename(
    ID = user_id
  )

# 1.1 - crea colonne date coerenti
df_clean2 <- df_clean %>%
  mutate(
    Data_date = as.Date(Data) # Data è POSIXct/dttm -> Date
  ) |>
  dplyr::rename(
    date = Data_date
  ) |>
  dplyr::select(-c(Data, pid5_na, pid5_det, pid5_ant, pid5_dis, pid5_psy))

d1_2 <- d1 %>%
  mutate(
    day_date = as.Date(day) # day è IDate -> Date
  ) |>
  dplyr::rename(
    date = day_date,
    timepoint = exam_period
  ) |>
  dplyr::select(-c(day))

# 1) vettore di date consentite (uniche)
# allowed_dates <- c(
#   "2005-04-06", "2025-03-01", "2025-03-04", "2025-03-29", "2025-03-30", "2025-03-31", "2025-04-01",
#   "2025-04-02", "2025-04-03", "2025-04-04", "2025-04-05", "2025-04-06", "2025-04-07", "2025-04-09",
#   "2025-04-11", "2025-04-14", "2025-04-15", "2025-05-16", "2025-05-23"
# )

# 3) Filtra mantenendo solo le date consentite
voice_df <- df_clean2

table(
  voice_df$ID,
  voice_df$timepoint
)

allowed_id <- unique(voice_df$ID)

pid5_df <- d1_2 %>%
  dplyr::filter(ID %in% allowed_id)

pid5_df$timepoint <- forcats::fct_recode(
  pid5_df$timepoint,
  "baseline" = "baseline",
  "pre" = "pre_exam",
  "post" = "post_exam"
)

length(unique(pid5_df$ID))
length(unique(voice_df$ID))

table(voice_df$ID, voice_df$timepoint)

###################

pid5_agg <- pid5_df %>%
  group_by(ID, timepoint) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "{.col}"), # medie per tutte le numeric
    sex = first(na.omit(sex)), # sesso "rappresentativo"
    date_rep = if ("date" %in% names(.)) min(date, na.rm = TRUE) else
      as.Date(NA), # data rappresentativa (min)
    .groups = "drop"
  )

# 3a) JOIN principale: aggiunge le colonne aggregate a voice_df, preservando tutte le colonne di voice_df
joined_left <- voice_df %>%
  left_join(pid5_agg, by = c("ID", "timepoint"))

imp <- missRanger(joined_left, num.trees = 100)

df_clean <- imp


# Extract baseline PID-5 scores (for between-person predictors)
pid5_baseline <- df_clean %>%
  dplyr::filter(timepoint == "baseline") %>%
  dplyr::select(ID, starts_with("pid5")) %>%
  rename_with(~ paste0(., "_bl"), .cols = starts_with("pid5"))

# Create long format with baseline PID-5 as between-person predictors
df_long <- df_clean %>%
  dplyr::select(-starts_with("pid5")) %>% # Remove timepoint-varying PID-5
  left_join(pid5_baseline, by = "ID") %>%
  # Center baseline PID-5 scores for interpretability
  mutate(
    across(
      ends_with("_bl"),
      ~ scale(., center = TRUE, scale = TRUE)[, 1],
      .names = "{.col}_c"
    )
  )

# Remove missing cases
df_analysis <- df_long %>%
  dplyr::filter(complete.cases(.))
