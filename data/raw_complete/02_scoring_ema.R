suppressPackageStartupMessages({
  library(psych)
  library(ggplot2)
  library(sjPlot)
  library(dplyr)
  library(rio)
  library(here)
  library(purrr)
  library(lubridate)
  library(readxl)
  library(janitor)
  library(stringr)
})


df <- rio::import(here(
  "data",
  "raw_complete",
  "data_interim",
  "ema_data_clean.RDS"
))


# SCORING PID-5 -------------------------------------------------------------------

# PID5
df$pid5_sum <- rowSums(df[, 5:19])
df$pid5_negative_affectivity <- rowSums(df[, c(9, 15, 8)])
df$pid5_detachment <- rowSums(df[, c(13, 17, 12)])
df$pid5_antagonism <- rowSums(df[, c(10, 16, 18)])
df$pid5_disinhibition <- rowSums(df[, c(19, 7, 14)])
df$pid5_psychoticism <- rowSums(df[, c(5, 6, 11)])

# IPV
df$ipv_sum <- rowSums(df[, 38:41])

# TRIPM
df$tripm_4_rev <- 5 - df$tripm_4

df$tripm_sum <- rowSums(df[, 20:23])
df$tripm_boldness <- rowSums(df[, c(20, 22)])
df$tripm_meanness <- rowSums(df[, c(21, 68)])

# DASS21
df$dass_sum <- rowSums(df[, 24:29])
df$dass_stress <- rowSums(df[, c(29, 24)])
df$dass_depression <- rowSums(df[, c(26, 27)])
df$dass_anxiety <- rowSums(df[, c(25, 28)])

# COPE
df$cope_10_rev <- 5 - df$cope_nvi_10

df$cope_avoid <- rowSums(df[, c(42, 48)])
df$cope_prob_or <- rowSums(df[, c(49, 50)])
df$cope_social_support <- rowSums(df[, c(44, 51)])
df$cope_positive_att <- rowSums(df[, c(43, 45)])
df$cope_trascendent_or <- rowSums(df[, c(46, 76)])

# SCS

df$cs_pos <- rowSums(df[, c(30, 32, 35, 37)])
df$ucs_neg <- rowSums(df[, c(31, 33, 34, 36)])

rio::export(
  df,
  here("data", "raw_complete", "data_interim", "ema_data_scoring.RDS")
)
