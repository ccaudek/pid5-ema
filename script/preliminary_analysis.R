library(tidyverse)
library(here)
library(rio)
library(brms)


esi_bf <- rio::import(
  here::here(
    "data",
    "processed",
    "esi_bf.csv"
  )
) |>
  dplyr::distinct(user_id, .keep_all = TRUE) |>
  dplyr::select(user_id, esi_bf)

pid5 <- rio::import(
  here::here(
    "data",
    "processed",
    "pid5.csv"
  )
) |>
  dplyr::distinct(user_id, .keep_all = TRUE) |>
  dplyr::select(user_id, starts_with("domain_"))


df <- left_join(esi_bf, pid5, by = "user_id")

fit <- brm(
  esi_bf ~
    domain_negative_affect +
      domain_detachment +
      domain_antagonism +
      domain_disinhibition +
      domain_psychoticism,
  family = skew_normal(),
  data = df,
  backend = "cmdstanr"
)
pp_check(fit)
summary(fit)

bayes_R2(fit)


hist(df$esi_bf)

conditional_effects(fit, "domain_negative_affect")
conditional_effects(fit, "domain_disinhibition")
conditional_effects(fit, "domain_psychoticism")
