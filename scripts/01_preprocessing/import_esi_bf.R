library(tidyverse)
library(here)
library(rio)
library(stringr)

# Abbiamo usato il reverse-scoring per gli item i_15, i_17 e i_18.

temp <- rio::import(
  here::here(
    "data",
    "raw",
    "quest",
    "2025_EMA_quest3.xlsx"
  )
)

esi_bdf <- temp[, c(1, 2, 99:116)]


df_recoded <- esi_bdf %>%
  # ---- (1) rinomina colonne --------------------------
  dplyr::rename(
    date = 1,
    user_id = 2
  ) %>%
  rename_with(
    .cols = 3:20,
    .fn = ~ sprintf("i_%02d", seq_along(.))
  ) %>%
  # ---- (2) ricodifica Vero/Falso in 1/0 --------------
  mutate(
    across(
      starts_with("i_"),
      ~ case_when(
        .x == "Vero" ~ 1,
        .x == "Falso" ~ 0,
        TRUE ~ NA_real_
      )
    )
  ) %>%
  # ---- (3) reverse‑scoring specifico -----------------
  mutate(
    across(
      c(i_15, i_17, i_18),
      ~ ifelse(is.na(.x), NA_real_, 1 - .x) # 1→0, 0→1, NA rimane NA
    )
  )


df_scored <- df_recoded %>% # <-- il tuo df già rinominato/ricodificato
  mutate(
    # Somma delle colonne che iniziano con "i_"
    esi_bf = rowSums(
      dplyr::select(., starts_with("i_")), # seleziona i_01 … i_18
      na.rm = TRUE # ignora eventuali NA
    )
  )


# Correct user_id

temp <- rio::import(
  here::here(
    "data",
    "raw",
    "quest",
    "PID-5_2025_EMA.xlsx"
  )
)

temp[, 3]


# "gi_st_07_22_311" # teniamo così
#
# "lo_ca_2005_05_07_05_437" # teniamo così
#
# "ch_ar_00_08_8_388" # ok
#
# "an_bi_2001_02_01" # ok
#
# "gi_mi_2005_12_854" # ok
#
# "be_mo_2005_29_815" # ok
#
# "vi_da_1999_11_499" # ok
#
# "gi_zo_1997_04_358_" # ok
#
# "gi_ri_2003_24_820" # ok
#
# "ma_ta_2005_07_21_866_"
# "el_ci_2002_02_15_057_"
# "ma_pa_2005_06_22_205_"
# "fla_me_2004_02_17_358"
# "el_ci_2002_02_15_057_"
# "an_va_2004_04_08_527_"
# "ca_mo_2004_09_02_700_"

df_scored <- df_scored %>%
  mutate(
    user_id = user_id %>%
      str_to_lower() %>% # tutto in minuscolo
      str_replace_all("__+", "_") %>% # "__" → "_"
      str_replace("_$", "") %>% # rimuove "_" finale
      str_replace("^fla_me_", "fl_me_") # fla_me_ → fl_me_
  )

esi_bf_scores <- df_scored |>
  dplyr::select(user_id, esi_bf) |>
  dplyr::rename(
    esi_bf_baseline = esi_bf
  )

rio::export(
  esi_bf_scores,
  here::here(
    "data",
    "processed",
    "esi_bf_scores.csv"
  )
)
