# 00_data_dictionary.R
# Build a reproducible data dictionary (CSV) for the PID-5 EMA project.
# - Reads the main processed file `ema_plus_baseline_exam_tagged.csv` if present,
#   otherwise constructs a dictionary template based on known variable name patterns.
# - Outputs: data/processed/data_dictionary.csv
#
# Run from project root (so that here::here() resolves correctly).

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(readr)
})

out_path <- here::here("data", "processed", "data_dictionary.csv")

# Helper: classify variables by regex/prefix
classify_var <- function(v) {
  case_when(
    v %in% c("user_id", "course", "sex") ~ "id/meta",
    v %in% c("day", "hour", "bysubj_day", "exam_period") ~ "time/exam",
    str_detect(v, "^pid5_.*_baseline$") ~ "baseline: PID-5 domains",
    v == "esi_bf_baseline" ~ "baseline: ESI-BF",
    str_detect(v, "^dass_(stress|anxiety|depression)_baseline$") ~
      "baseline: DASS-21",
    str_detect(
      v,
      "^pid5_(negative_affectivity|detachment|antagonism|disinhibition|psychoticism)$"
    ) ~
      "EMA: PID-5 (subset)",
    str_detect(v, "^cope_") ~ "EMA: COPE",
    str_detect(v, "^tripm_") ~ "EMA: Tri-PM",
    str_detect(v, "^ipv_") ~ "EMA: IPV",
    v %in% c("cs_pos", "ucs_neg") ~ "EMA: SCS (state composites)",
    v %in% c("happy", "sad", "satisfied", "angry", "neg_affect_ema") ~
      "EMA: mood / negative affect",
    str_detect(v, "^context_(quality|control|support|threat)$") ~
      "EMA: context",
    TRUE ~ "other"
  )
}

# Try to read the processed dataset; otherwise create a skeleton list from known names
dat_path <- here::here("data", "processed", "ema_plus_baseline_exam_tagged.csv")
if (file.exists(dat_path)) {
  d <- suppressMessages(readr::read_csv(
    dat_path,
    show_col_types = FALSE,
    n_max = 100
  ))
  vars <- names(d)
} else {
  vars <- c(
    "user_id",
    "course",
    "sex",
    "day",
    "hour",
    "bysubj_day",
    "exam_period",
    "esi_bf_baseline",
    "pid5_negative_affect_baseline",
    "pid5_detachment_baseline",
    "pid5_antagonism_baseline",
    "pid5_disinhibition_baseline",
    "pid5_psychoticism_baseline",
    "dass_stress_baseline",
    "dass_anxiety_baseline",
    "dass_depression_baseline",
    "pid5_negative_affectivity",
    "pid5_detachment",
    "pid5_antagonism",
    "pid5_disinhibition",
    "pid5_psychoticism",
    "cope_avoid",
    "cope_prob_or",
    "cope_social_support",
    "cope_positive_att",
    "cope_trascendent_or",
    "tripm_boldness",
    "tripm_meanness",
    "ipv_sum",
    "context_quality",
    "context_control",
    "context_support",
    "context_threat",
    "happy",
    "sad",
    "satisfied",
    "angry",
    "cs_pos",
    "ucs_neg",
    "neg_affect_ema"
  )
}

# Build dictionary with short descriptions (editable)
dict <- tibble::tibble(
  variable = vars,
  group = purrr::map_chr(vars, classify_var),
  description = dplyr::case_when(
    variable == "user_id" ~ "Anonymous subject code",
    variable == "course" ~
      "Course the subject belongs to (Psicometria, Testing, Interventi, Clinica)",
    variable == "sex" ~ "Biological sex as reported (F/M/other coding)",
    variable == "day" ~ "Date of EMA response (YYYY-MM-DD)",
    variable == "hour" ~ "Hour of EMA response (0-23)",
    variable == "bysubj_day" ~
      "Sequential EMA index within subject (if available)",
    variable == "exam_period" ~
      "Exam proximity tag: baseline / pre_exam / post_exam",
    variable == "esi_bf_baseline" ~ "ESI-BF total score at baseline",
    variable == "pid5_negative_affect_baseline" ~
      "PID-5 domain Negative Affect (baseline, EMA items removed)",
    variable == "pid5_detachment_baseline" ~
      "PID-5 domain Detachment (baseline, EMA items removed)",
    variable == "pid5_antagonism_baseline" ~
      "PID-5 domain Antagonism (baseline, EMA items removed)",
    variable == "pid5_disinhibition_baseline" ~
      "PID-5 domain Disinhibition (baseline, EMA items removed)",
    variable == "pid5_psychoticism_baseline" ~
      "PID-5 domain Psychoticism (baseline, EMA items removed)",
    variable == "dass_stress_baseline" ~ "DASS-21 Stress (baseline)",
    variable == "dass_anxiety_baseline" ~ "DASS-21 Anxiety (baseline)",
    variable == "dass_depression_baseline" ~ "DASS-21 Depression (baseline)",
    variable == "pid5_negative_affectivity" ~
      "EMA PID-5 Negative Affectivity composite (subset items)",
    variable == "pid5_detachment" ~
      "EMA PID-5 Detachment composite (subset items)",
    variable == "pid5_antagonism" ~
      "EMA PID-5 Antagonism composite (subset items)",
    variable == "pid5_disinhibition" ~
      "EMA PID-5 Disinhibition composite (subset items)",
    variable == "pid5_psychoticism" ~
      "EMA PID-5 Psychoticism composite (subset items)",
    variable == "cope_avoid" ~ "EMA COPE: Avoidance",
    variable == "cope_prob_or" ~ "EMA COPE: Problem-oriented coping",
    variable == "cope_social_support" ~ "EMA COPE: Social support",
    variable == "cope_positive_att" ~ "EMA COPE: Positive attitude",
    variable == "cope_trascendent_or" ~ "EMA COPE: Transcendent orientation",
    variable == "tripm_boldness" ~ "EMA Tri-PM: Boldness",
    variable == "tripm_meanness" ~ "EMA Tri-PM: Meanness",
    variable == "ipv_sum" ~ "EMA Intimate Partner Violence (sum)",
    variable == "context_quality" ~ "EMA context appraisal: quality (1-5)",
    variable == "context_control" ~ "EMA context appraisal: control (1-5)",
    variable == "context_support" ~ "EMA context appraisal: support (1-5)",
    variable == "context_threat" ~ "EMA context appraisal: threat (1-5)",
    variable == "happy" ~ "EMA mood (0-100, higher = more happy)",
    variable == "sad" ~ "EMA mood (0-100, higher = more sad)",
    variable == "satisfied" ~ "EMA mood (0-100, higher = more satisfied)",
    variable == "angry" ~ "EMA mood (0-100, higher = more angry)",
    variable == "cs_pos" ~
      "State Self-Compassion: positive component (sum of SCS*_pos items)",
    variable == "ucs_neg" ~
      "State Self-Compassion: uncompassionate/negative component (sum of SCS*_neg items)",
    variable == "neg_affect_ema" ~
      "Standardized composite negative affect: reverse(happy,satisfied) + sad + angry (z-score)",
    TRUE ~ ""
  )
) %>%
  arrange(group, variable)

readr::write_csv(dict, out_path)
message("Wrote data dictionary to: ", out_path)
print(dict, n = 200)
