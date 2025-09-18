# Load necessary libraries
library(tidyverse)
library(here)
library(rio)
library(brms)
library(stringr)
library(purrr)
library(cmdstanr)
options(mc.cores = parallel::detectCores())
library(loo)
library(ppcor)
library(tidyr)
library(broom)
library(tibble)
library(mice)
library(corrplot)


# Read and process 'esi_bf' data
esi_bf <- rio::import(
  here::here(
    "data",
    "processed",
    "esi_bf.csv"
  )
) |>
  dplyr::distinct(user_id, .keep_all = TRUE) |> # Keep only distinct user_id
  dplyr::select(user_id, esi_bf) # Select relevant columns

# Read and process 'pid5' data
pid5 <- rio::import(
  here::here(
    "data",
    "processed",
    "pid5.csv"
  )
) |>
  dplyr::distinct(user_id, .keep_all = TRUE) |> # Keep only distinct user_id
  dplyr::select(user_id, starts_with("domain_")) # Select domain variables

# Merge 'esi_bf' and 'pid5' data by user_id
df <- left_join(esi_bf, pid5, by = "user_id")


# Define list of user IDs with careless responding
user_id_with_careless_responding <- c(
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

# Filter out users with careless responses
df1 <- df[!(df$user_id %in% user_id_with_careless_responding), ]

# Read EMA data and rename 'subj_code' to 'user_id'
ema_raw <- readRDS(
  here::here(
    "data",
    "raw",
    "ema",
    "ema_data_scoring.RDS"
  )
) |>
  dplyr::rename(
    user_id = subj_code
  )

# Merge EMA data with filtered main data
df2 <- left_join(df1, ema_raw, by = "user_id")

# Verify number of unique users
length(unique(df2$user_id))

### Scopo: Analizzare la moderazione esercitata dai tratti PID-5 dinamici (EMA) nella relazione tra affetto negativo momentaneo e UCS (Uncompassionate Self-Responding)

# Calcola la correlazione tra CS e UCS per ogni partecipante
cor_df <- df2 %>%
  group_by(user_id) %>%
  filter(!is.na(cs_pos) & !is.na(ucs_neg)) %>%
  filter(n() >= 2) %>% # Almeno due osservazioni complete
  summarise(corr_cs_ucs = cor(cs_pos, ucs_neg, use = "complete.obs")) %>%
  ungroup()

# Calcola le deviazioni standard per CS e UCS per ogni soggetto
var_df <- df2 %>%
  group_by(user_id) %>%
  summarise(
    sd_cs = sd(cs_pos, na.rm = TRUE),
    sd_ucs = sd(ucs_neg, na.rm = TRUE)
  ) %>%
  ungroup()

# Costruisce una misura media dell'affetto negativo momentaneo
# NOTA: assumiamo che le variabili "sad", "angry" siano affetti negativi,
# "happy" positivo quindi lo escludiamo dal calcolo.
# df2 <- df2 %>%
#   mutate(
#     neg_aff_ema = rowMeans(dplyr::select(., sad, angry), na.rm = TRUE)
#   )

md.pattern(df2[, c("sad", "angry", "happy", "satisfied")], plot = TRUE)

# Seleziona solo le colonne rilevanti (per velocità)
items <- c("sad", "angry", "happy", "satisfied")

# Imputa i missing (1 solo imputazione, dato che i NA sono pochi)
imputed <- mice(df2[, items], m = 1, maxit = 10, seed = 123)

# Estrai il dataset imputato e sostituisci le colonne originali
df2_imputed <- complete(imputed)
df2[, items] <- df2_imputed[, items]

df2 <- df2 %>%
  mutate(
    # Calcola differenze attese (es. sad dovrebbe essere opposto a happy)
    incoherence = abs(sad - (100 - happy)) + abs(angry - (100 - satisfied)),
    # Standardizza l'indice (valori alti = potenziali careless responders)
    incoherence_z = scale(incoherence)
  )

careless_subjects <- df2 %>%
  group_by(user_id) %>%
  summarise(
    n_borderline = sum(incoherence_z > quantile(incoherence_z, 0.95)), # Conteggio osservazioni "sospette"
    mean_incoherence = mean(incoherence_z)
  ) %>%
  dplyr::filter(
    n_borderline >= 3 | mean_incoherence > quantile(mean_incoherence, 0.95)
  ) %>%
  pull(user_id)

# Numero di soggetti da escludere (5% del totale)
n_esclusi <- length(careless_subjects)
n_totali <- n_distinct(df2$user_id)
percentuale_esclusi <- round(n_esclusi / n_totali * 100, 1)
cat(
  "Soggetti esclusi:",
  n_esclusi,
  "/",
  n_totali,
  "(",
  percentuale_esclusi,
  "%)"
)

df2 %>%
  dplyr::filter(user_id %in% careless_subjects) %>%
  summarise(
    mean_sad = mean(sad),
    mean_happy = mean(happy),
    mean_incoherence = mean(incoherence_z),
    prop_extreme = mean(
      sad == 100 | happy == 100 | angry == 100 | satisfied == 100
    )
  )


df2_clean <- df2 %>%
  dplyr::filter(!user_id %in% careless_subjects)

# Verifica le dimensioni finali
length(unique(df2_clean$user_id)) # Dovresti vedere ~407 soggetti rimanenti

# Criteri usati: ≥3 osservazioni con incoherence_z > 95° percentile oppure
# mean_incoherence > 95° percentile tra i soggetti.
# Percentuale esclusa: ~5.4%.

df2 <- df2_clean %>%
  mutate(
    happy_reversed = 100 - happy, # Scala 0-100
    satisfied_reversed = 100 - satisfied,
    neg_aff_ema = rowMeans(
      cbind(sad, angry, happy_reversed, satisfied_reversed),
      na.rm = TRUE
    )
  )


# Stima la reattività affettiva UCS ~ NA per ogni soggetto (slopes idiografici)
slopes_df <- df2 %>%
  group_by(user_id) %>%
  dplyr::filter(n() >= 10) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(ucs_neg ~ neg_aff_ema, data = .x)),
    tidied = map(model, ~ broom::tidy(.x))
  ) %>%
  unnest(tidied) %>%
  dplyr::filter(term == "neg_aff_ema") %>%
  dplyr::select(user_id, slope_ucs_neg_aff = estimate)

# Unisce correlazioni, varianze e slope in un unico data frame
all_df <- reduce(list(cor_df, var_df, slopes_df), full_join, by = "user_id")

# Standardizza affetto negativo all'interno di ciascun soggetto
# per ridurre la collinearità con i tratti
df2 <- df2 %>%
  group_by(user_id) %>%
  mutate(neg_aff_ema_c = scale(neg_aff_ema)) %>%
  ungroup()

# Prepara dataset finale con slope e predittori PID-5 EMA
final_df <- df2 %>%
  dplyr::select(user_id, starts_with("pid5_")) %>%
  dplyr::distinct() %>%
  left_join(all_df, by = "user_id") %>%
  mutate(across(where(is.numeric) & !matches("user_id"), scale))

# Modello lineare: i tratti EMA predicono la reattività affettiva idiografica
model_df <- final_df %>%
  dplyr::filter(
    complete.cases(
      dplyr::across(
        c(slope_ucs_neg_aff, dplyr::starts_with("pid5_"))
      )
    )
  )

model <- lm(
  slope_ucs_neg_aff ~
    pid5_negative_affectivity +
      pid5_detachment +
      pid5_antagonism +
      pid5_disinhibition +
      pid5_psychoticism,
  data = model_df
)
summary(model)

# Modello bayesiano di base: UCS ~ affetto negativo momentaneo
model_base <- brm(
  ucs_neg ~ neg_aff_ema_c + (1 + neg_aff_ema_c | user_id),
  data = df2,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  cores = 4,
  iter = 2000,
  seed = 123,
  backend = "cmdstanr"
)

# Modello alternativo: effetto moderatore dei tratti EMA
model_alt <- brm(
  ucs_neg ~
    neg_aff_ema_c *
      (pid5_negative_affectivity +
        pid5_detachment +
        pid5_antagonism +
        pid5_disinhibition +
        pid5_psychoticism) +
      (1 + neg_aff_ema_c | user_id),
  data = df2,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  cores = 4,
  iter = 2000,
  seed = 123,
  backend = "cmdstanr"
)

pp_check(model_base)
pp_check(model_alt)

# Confronto predittivo tramite ELPD
loo_base <- loo(model_base)
loo_alt <- loo(model_alt)
loo_compare(loo_base, loo_alt)

# R2 bayesiano per valutare la varianza spiegata
bayes_R2(model_base)
bayes_R2(model_alt)

# Esplora interazioni: visualizza effetti condizionali
conditional_effects(model_alt, effects = "neg_aff_ema_c:pid5_detachment")
conditional_effects(model_alt, effects = "neg_aff_ema_c:pid5_disinhibition")

# Analizza slope stimati per soggetto e relazione con i tratti
slope_post <- ranef(model_alt)$user_id[,, "neg_aff_ema_c"]
slope_df <- tibble(
  user_id = rownames(slope_post),
  neg_aff_slope = slope_post[, "Estimate"]
)
subject_traits <- df2 %>%
  dplyr::select(user_id, starts_with("pid5_")) %>%
  dplyr::distinct()
left_join(slope_df, subject_traits, by = "user_id") %>%
  ggplot(aes(x = pid5_negative_affectivity, y = neg_aff_slope)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Does trait NA predict reactivity (UCS ~ neg_aff)?")


print(model_alt)
