# ============================================================
# Preprocessing PID-5 & feature acustiche: media delle imputazioni
# Output:
#   - acoustic_pid5_long_avg.csv (medio per imputazione, formato lungo)
#   - acoustic_pid5_delta.csv    (wide a livello soggetto, con Δ post–pre)
# ============================================================

# ============================================================
# 0) Setup
# ============================================================

suppressPackageStartupMessages({
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
  # library(tidyr)
  library(broom)
  library(tibble)
  library(mice)
  library(lubridate)
  library(readxl)
  library(janitor)
  library(mice)
  library(miceadds)
  library(posterior)
  # library(tidybayes)
  library(lmtest)
  library(sandwich)
  library(metafor)
  library(robustbase)
  library(lme4)
  library(broom.mixed)
  library(conflicted)
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("var", "stats")
conflict_prefer("mad", "stats")
conflict_prefer("ar", "brms")
conflict_prefer("chisq.test", "stats")
conflict_prefer("dstudent_t", "brms")
conflict_prefer("epilepsy", "brms")
conflict_prefer("expand", "Matrix")
conflict_prefer("factorize", "lme4")
conflict_prefer("fisher.test", "stats")
conflict_prefer("lag", "stats")
conflict_prefer("match", "base")
conflict_prefer("milk", "loo")
conflict_prefer("ngrps", "brms")
conflict_prefer("pack", "Matrix")
conflict_prefer("sd", "stats")
conflict_prefer("unpack", "Matrix")

set.seed(1234)

# ------------------------------------------------------------
# 1) Lettura questionari (ESI-BF + PID-5) e pulizia
# ------------------------------------------------------------
esi_bf <- rio::import(here::here("data", "processed", "esi_bf.csv")) |>
  dplyr::distinct(user_id, .keep_all = TRUE) |>
  dplyr::select(user_id, esi_bf)

pid5 <- rio::import(here::here("data", "processed", "pid5.csv")) |>
  dplyr::distinct(user_id, .keep_all = TRUE) |>
  dplyr::select(user_id, starts_with("domain_"))

quest_df <- left_join(esi_bf, pid5, by = "user_id")

# eventuale filtro risposte "careless"
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

quest_df <- quest_df |> filter(!(user_id %in% user_id_with_careless_responding))

# ------------------------------------------------------------
# 2) Lettura Excel con feature acustiche (Baseline/Pre/Post)
#    + normalizzazione nomi + parsing Case -> day + user_id
# ------------------------------------------------------------
file_path <- here::here(
  "data",
  "raw",
  "acustic_features",
  "database_acfeat.xlsx"
)

fix_names <- function(nm) {
  nm %>%
    str_replace_all("[\\s\\u00A0]+", " ") %>%
    str_trim() %>%
    str_replace_all("\\s*/\\s*", "_") %>%
    str_replace_all("/+", "_") %>%
    str_squish()
}

wanted_sheets <- c("Baseline", "Pre", "Post")
sheets <- intersect(readxl::excel_sheets(file_path), wanted_sheets)

df1 <- purrr::map_dfr(
  sheets,
  ~ readxl::read_excel(file_path, sheet = .x) |>
    dplyr::rename_with(fix_names) |> # pulisce gli header (spazi, slash, NBSP)
    janitor::clean_names(case = "snake") |> # snake_case finale
    mutate(sheet = .x)
)

# La colonna 'case' deve esistere: contiene data e user_id separati da " - "
stopifnot("case" %in% names(df1))

df1 <- df1 %>%
  tidyr::separate_wider_regex(
    case,
    patterns = c(
      day = "\\d{2}[-_]\\d{2}[-_]\\d{4}",
      "\\s*-\\s*",
      user_id = ".+"
    )
  ) %>%
  mutate(
    day = lubridate::dmy(stringr::str_replace_all(day, "[-_]", "/")),
    phase = factor(
      sheet,
      levels = c("Baseline", "Pre", "Post"),
      labels = c("neutral", "pre", "post")
    )
  ) %>%
  dplyr::select(-sheet)

# ------------------------------------------------------------
# 3) Merge con questionari
# ------------------------------------------------------------
df2 <- left_join(df1, quest_df, by = "user_id")

# ------------------------------------------------------------
# 4) Definizione colonne di interesse
# ------------------------------------------------------------
pid5_cols <- c(
  "domain_negative_affect",
  "domain_detachment",
  "domain_antagonism",
  "domain_disinhibition",
  "domain_psychoticism"
)

acoustic_feats <- c(
  "f0_min_hz_a",
  "f0_median_hz_i",
  "f0_mean_hz_u",
  "f0_median_hz_u",
  "f1_min_hz_u",
  "mfcc3_std",
  "mfcc9_std",
  "mfcc11_skewness",
  "mfcc11_median",
  "mfcc3_iqr",
  "mfcc9_iqr",
  "mfcc9_p25th"
)

# Check che le colonne esistano
missing_pid <- setdiff(pid5_cols, names(df2))
missing_aco <- setdiff(acoustic_feats, names(df2))
if (length(missing_pid) > 0)
  warning("Mancano PID-5: ", paste(missing_pid, collapse = ", "))
if (length(missing_aco) > 0)
  warning("Mancano feature acustiche: ", paste(missing_aco, collapse = ", "))

# ------------------------------------------------------------
# 5) IMPUTAZIONE SOLO PID-5 a livello soggetto (m semplice)
#    -> media per soggetto (collasso imputazioni)
# ------------------------------------------------------------
pid_subject <- df2 %>%
  dplyr::select(user_id, esi_bf, all_of(pid5_cols)) %>%
  dplyr::distinct(user_id, .keep_all = TRUE)

# Schema metodi: PMM per domini + ESI-BF (se serve)
meth <- make.method(pid_subject)
meth[pid5_cols] <- "pmm"
meth["esi_bf"] <- "pmm" # opzionale

# PredMatrix: non usare user_id come predittore
pred <- make.predictorMatrix(pid_subject)
pred[, "user_id"] <- 0
pred["user_id", ] <- 0

imp_pid <- mice(
  pid_subject,
  m = 20,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  printFlag = TRUE
)

# Lista dei dataset imputati
pid_imputed_list <- lapply(
  1:imp_pid$m,
  function(k) complete(imp_pid, action = k)
)

# Media per soggetto (collasso imputazioni)
pid_avg <- bind_rows(pid_imputed_list, .id = ".imp") |>
  group_by(user_id) |>
  summarise(
    across(c(esi_bf, all_of(pid5_cols)), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Arrotonda i 5 domini se desideri restare su interi
pid_avg <- pid_avg |> mutate(across(all_of(pid5_cols), ~ round(.x)))

# ------------------------------------------------------------
# 6) Rientro nel dataset lungo con le registrazioni vocali
#    sostituendo i PID-5 con la versione media imputata
# ------------------------------------------------------------
df_long <- df2 %>%
  dplyr::select(-all_of(pid5_cols), -esi_bf) %>%
  left_join(pid_avg, by = "user_id")

# ------------------------------------------------------------
# 7) IMPUTAZIONE PRELIMINARE DELLE FEATURE ACUSTICHE
#    (una sola passata, semplice; non imputiamo id/phase/day)
# ------------------------------------------------------------
# Costruiamo un DF minimale per imputare le feature acustiche
impute_df <- df_long %>%
  dplyr::select(
    user_id,
    day,
    phase,
    all_of(pid5_cols),
    esi_bf,
    all_of(acoustic_feats)
  )

# mice richiede fattori come tali
impute_df <- impute_df %>%
  mutate(
    phase = factor(phase, levels = c("neutral", "pre", "post"))
  )

# Metodi: PMM per tutte le feature acustiche (robusto), niente imputazione per id/day/phase
meth2 <- make.method(impute_df)
meth2[c("user_id", "day", "phase")] <- "" # non imputare
meth2[pid5_cols] <- "" # sono già stati imputati/mediati
meth2["esi_bf"] <- "" # idem
meth2[acoustic_feats] <- "pmm" # imputiamo solo le feature vocali

pred2 <- make.predictorMatrix(impute_df)
# non usare user_id come predittore
pred2[, "user_id"] <- 0
pred2["user_id", ] <- 0
# possiamo usare day e phase come predittori
pred2["day", ] <- 0
pred2[, "day"] <- 1
pred2["day", "day"] <- 0 # day solo come predittore
pred2["phase", ] <- 0
pred2[, "phase"] <- 1 # phase solo come predittore

imp_aco <- mice(
  impute_df,
  m = 20,
  maxit = 15,
  method = meth2,
  predictorMatrix = pred2,
  printFlag = TRUE
)

# Collassiamo le imputazioni acustiche con la media per (user_id, phase, day)
aco_avg <- bind_rows(
  lapply(1:imp_aco$m, \(k) complete(imp_aco, k)),
  .id = ".imp"
) |>
  group_by(user_id, phase, day) |>
  summarise(
    across(
      all_of(c(pid5_cols, "esi_bf", acoustic_feats)),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 8) Dataset "lungo" medio sulle imputazioni + z-score PID-5
# ------------------------------------------------------------
acoustic_pid5_long_avg <- aco_avg |>
  mutate(across(
    all_of(pid5_cols),
    \(x) as.numeric(scale(x)),
    .names = "{.col}_z"
  ))

readr::write_csv(acoustic_pid5_long_avg, "acoustic_pid5_long_avg.csv")
message(">> Salvato: acoustic_pid5_long_avg.csv")

# ------------------------------------------------------------
# 9) Dataset "wide" per soggetto con Δ(post - pre) delle feature
# ------------------------------------------------------------
# restringo a una riga per (user_id, phase), mediando su day se ci sono più registrazioni
long_phase <- acoustic_pid5_long_avg |>
  group_by(user_id, phase) |>
  summarise(
    across(
      all_of(c(pid5_cols, paste0(pid5_cols, "_z"), "esi_bf", acoustic_feats)),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

wide_acoustic <- long_phase |>
  dplyr::select(user_id, phase, all_of(acoustic_feats)) |>
  pivot_wider(
    id_cols = user_id,
    names_from = phase,
    values_from = all_of(acoustic_feats),
    names_sep = "_"
  )

# Δ post–pre per ogni feature acustica
for (v in acoustic_feats) {
  post_name <- paste0(v, "_post")
  pre_name <- paste0(v, "_pre")
  delta_name <- paste0(v, "_delta_post_pre")
  if (all(c(post_name, pre_name) %in% names(wide_acoustic))) {
    wide_acoustic[[delta_name]] <- wide_acoustic[[post_name]] -
      wide_acoustic[[pre_name]]
  }
}

# PID-5 medi per soggetto (dalla tabella long_phase)
pid5_by_id <- long_phase |>
  group_by(user_id) |>
  summarise(
    across(all_of(pid5_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
    across(
      paste0(pid5_cols, "_z"),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    esi_bf_mean = mean(esi_bf, na.rm = TRUE),
    .groups = "drop"
  )

acoustic_pid5_delta <- wide_acoustic |>
  left_join(pid5_by_id, by = "user_id")

readr::write_csv(acoustic_pid5_delta, "acoustic_pid5_delta.csv")
message(">> Salvato: acoustic_pid5_delta.csv")

# ------------------------------------------------------------
# 10) Check rapidi
# ------------------------------------------------------------
message(">> Varianze PID-5 (z) sul dataset lungo medio:")
print(
  acoustic_pid5_long_avg |>
    summarise(across(paste0(pid5_cols, "_z"), ~ var(.x, na.rm = TRUE)))
)

delta_cols <- grep("_delta_post_pre$", names(acoustic_pid5_delta), value = TRUE)
message(">> Riepilogo Δ(post-pre) delle feature acustiche:")
print(summary(acoustic_pid5_delta[delta_cols]))

# 1) Completezza per fase (neutral/pre/post) per ogni soggetto
tab_fasi <- acoustic_pid5_long_avg |>
  count(user_id, phase) |>
  tidyr::pivot_wider(names_from = phase, values_from = n, values_fill = 0)

# chi ha sia pre che post (necessario per Δ)
with_prepost <- tab_fasi |> filter(pre > 0 & post > 0)
nrow(with_prepost)
n_distinct(acoustic_pid5_long_avg$user_id)


# 2) Hai calcolato Δ solo quando esistono sia pre che post?

delta_cols <- grep("_delta_post_pre$", names(acoustic_pid5_delta), value = TRUE)
colSums(is.na(acoustic_pid5_delta[delta_cols]))

# Se compaiono NA in massa per una feature, significa che per quella feature il
# pre o il post mancano spesso.

# 3) Duplicati di (user_id, day, phase)

dup <- acoustic_pid5_long_avg |>
  count(user_id, day, phase) |>
  dplyr::filter(n > 1)
dup

# Se ci sono duplicati, valuta se mediare all’interno di (user_id, phase) (lo
# hai già fatto più avanti, bene) o scegliere un criterio (p.es. migliore SNR).

# 4) Outliers nei Δ (soprattutto F0/formanti)

flag_out <- function(x)
  abs(x - median(x, na.rm = TRUE)) > 3 * mad(x, na.rm = TRUE)

out_df <- acoustic_pid5_delta |>
  transmute(
    user_id,
    across(all_of(delta_cols), flag_out, .names = "{.col}_out")
  )
colSums(out_df[-1], na.rm = TRUE)

# Se alcuni Δ hanno molti outlier, conviene fare:
# * un controllo qualità sulle registrazioni corrispondenti;
# * in analisi preliminare, usare robust regression o winsorize (solo per esplorazione).

# 5) Coerenza unità/scale delle MFCC

# Controlla che MFCC siano comparabili (stesso front-end, stessa normalizzazione).
# In caso di scale molto diverse:
# z-score solo per le analisi (non per salvare il "dati grezzi")
delta_z <- acoustic_pid5_delta |>
  mutate(across(
    all_of(delta_cols),
    ~ as.numeric(scale(.x)),
    .names = "{.col}_z"
  ))


# 6) Varianza PID-5 (già ok) + varianza delle feature per fase

# Varianza feature per fase (diagnostica)
acoustic_pid5_long_avg |>
  group_by(phase) |>
  summarise(across(
    c(
      f0_min_hz_a,
      f0_median_hz_i,
      f0_mean_hz_u,
      f0_median_hz_u,
      f1_min_hz_u,
      mfcc3_std,
      mfcc9_std,
      mfcc11_skewness,
      mfcc11_median,
      mfcc3_iqr,
      mfcc9_iqr,
      mfcc9_p25th
    ),
    ~ var(.x, na.rm = TRUE)
  ))

# A) Variazione pre–post (mixed model semplice per ogni feature) ---------------

feat <- acoustic_feats[1] # cambia la feature

df_phase <- acoustic_pid5_long_avg |>
  group_by(user_id, phase) |>
  summarise(val = mean(.data[[feat]], na.rm = TRUE), .groups = "drop") |>
  filter(phase %in% c("pre", "post")) |>
  mutate(phase = relevel(phase, ref = "pre"))

m <- lmer(val ~ phase + (1 | user_id), data = df_phase)
broom.mixed::tidy(m)

# Interpretazione: il coefficiente di phasepost ≈ Δ(post–pre)
# stimato controllando per soggetto.

# B) Δ(post–pre) ~ PID-5 (ipotesi “più vulnerabilità ⇒ Δ maggiore”) ------------
### 6, 7, 10, 11
# queste quattro feature appartengono alla stessa “famiglia” → sono indici di 
# variabilità temporale di MFCC di basso e medio ordine, quindi tutte misurano 
# quanto cambia nel tempo l’inviluppo spettrale della voce, a diversi livelli di 
# dettaglio.


feat <- acoustic_feats[13] # cambia la feature
y <- paste0(feat, "_delta_post_pre")

df_delta <- acoustic_pid5_delta |>
  select(user_id, all_of(y), ends_with("_z_mean")) |> # i PID-5 z medi per soggetto
  rename_with(~ sub("_z_mean$", "", .x), ends_with("_z_mean")) # nomi puliti

# regressione lineare (preliminare)
fm <- reformulate(
  c(
    "domain_negative_affect",
    "domain_detachment",
    "domain_antagonism",
    "domain_disinhibition",
    "domain_psychoticism"
  ),
  response = y
)

fit_lm <- lm(fm, data = df_delta)
summary(fit_lm)


# Bayes con brms + LOO, priors deboli e regolarizzazione
pri <- c(
  prior(normal(0, 0.2), class = "b"),
  prior(student_t(3, 0, 0.5), class = "sigma")
)
fit_brms <- brm(
  fm,
  data = df_delta,
  family = gaussian(),
  prior = pri,
  iter = 4000,
  chains = 4,
  seed = 123,
  backend = "cmdstanr"
)
summary(fit_brms)
pp_check(fit_brms)



################################################################################
