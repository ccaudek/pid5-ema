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


df1 <- df[!(df$user_id %in% user_id_with_careless_responding), ]


fit <- brm(
  esi_bf ~
    domain_negative_affect +
      domain_detachment +
      domain_antagonism +
      domain_disinhibition +
      domain_psychoticism,
  family = asym_laplace(),
  data = df1,
  backend = "cmdstanr"
)
pp_check(fit)
summary(fit)

bayes_R2(fit)


hist(df$esi_bf)

conditional_effects(fit, "domain_negative_affect")
conditional_effects(fit, "domain_disinhibition")
conditional_effects(fit, "domain_psychoticism")


# Import EMA data -------------------------------------------------------------

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

df2 <- left_join(df1, ema_raw, by = "user_id")
length(unique(df2$user_id))


# Scegli le variabili temporali da decomporsi
vars <- c("pid5_negative_affectivity", "context_threat", "happy") # puoi aggiungere

# Centro le variabili su tre livelli
ema_long <- df2 %>%
  group_by(user_id) %>%
  mutate(across(
    all_of(vars),
    .fns = list(
      person = ~ mean(.x, na.rm = TRUE)
    ),
    .names = "{.col}_person"
  )) %>%
  group_by(user_id, day) %>%
  mutate(across(
    all_of(vars),
    .fns = list(
      day = ~ mean(.x, na.rm = TRUE)
    ),
    .names = "{.col}_day"
  )) %>%
  mutate(across(
    all_of(vars),
    .fns = list(
      moment = ~ .x - get(paste0(cur_column(), "_day"))
    ),
    .names = "{.col}_moment"
  )) %>%
  mutate(across(
    paste0(vars, "_day"),
    .fns = list(
      centered = ~ .x - get(sub("_day", "_person", cur_column()))
    ),
    .names = "{.col}_centered"
  )) %>%
  ungroup()


form <- bf(
  I(esi_bf + .01) ~
    pid5_negative_affectivity_moment +
      pid5_negative_affectivity_day_centered +
      pid5_negative_affectivity_person +
      context_threat_moment +
      context_threat_day_centered +
      context_threat_person +
      happy_moment +
      happy_day_centered +
      happy_person +
      (1 +
        pid5_negative_affectivity_moment +
        context_threat_moment ||
        user_id) +
      (1 | user_id:day)
)

fit_ema <- brm(
  formula = form,
  data = ema_long,
  family = skew_normal(),
  backend = "cmdstanr",
  # control = list(adapt_delta = 0.99),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  # seed = 123,
  algorithm = "fullrank"
)

print(fit_ema)


cor_df <- df2 %>%
  group_by(user_id) %>%
  dplyr::filter(!is.na(cs_pos) & !is.na(ucs_neg)) %>%
  dplyr::filter(n() >= 2) %>% # almeno due coppie complete
  summarise(corr_cs_ucs = cor(cs_pos, ucs_neg, use = "complete.obs")) %>%
  ungroup()

var_df <- df2 %>%
  group_by(user_id) %>%
  summarise(
    sd_cs = sd(cs_pos, na.rm = TRUE),
    sd_ucs = sd(ucs_neg, na.rm = TRUE)
  ) %>%
  ungroup()

df2 <- df2 %>%
  mutate(
    neg_aff_ema = rowMeans(
      dplyr::select(., sad, angry, happy = -happy),
      na.rm = TRUE
    )
  )


# Per ogni soggetto, stima UCS ~ neg_aff
slopes_df <- df2 %>%
  group_by(user_id) %>%
  dplyr::filter(n() >= 10) %>% # solo soggetti con abbastanza dati
  nest() %>%
  mutate(
    model = map(data, ~ lm(ucs_neg ~ neg_aff_ema, data = .x)),
    tidied = map(model, ~ tidy(.x))
  ) %>%
  unnest(tidied) %>%
  dplyr::filter(term == "neg_aff_ema") %>%
  dplyr::select(user_id, slope_ucs_neg_aff = estimate)


all_df <- reduce(
  list(cor_df, var_df, slopes_df),
  full_join,
  by = "user_id"
)

df2 <- df2 %>%
  group_by(user_id) %>%
  mutate(neg_aff_ema_c = scale(neg_aff_ema, scale = FALSE)) %>%
  ungroup()

final_df <- df2 %>%
  dplyr::select(user_id, starts_with("pid5_")) %>%
  dplyr::distinct() %>%
  left_join(all_df, by = "user_id")

model_df <- final_df %>%
  dplyr::filter(complete.cases(
    slope_ucs_neg_aff,
    pid5_negative_affectivity,
    pid5_detachment,
    pid5_antagonism
  ))

model <- lm(
  slope_ucs_neg_aff ~
    pid5_negative_affectivity + pid5_detachment + pid5_antagonism,
  data = model_df
)

summary(model)


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

pp_check(model_base)

print(model_base)
bayes_R2(model_base)


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

pp_check(model_alt)
print(model_alt)
bayes_R2(model_alt)


loo0 <- loo(model_base) # baseline model
loo1 <- loo(model_alt) # augmented model
loo_compare(loo0, loo1)


###############################

ema_sd <- df2 %>% # df2 = dati EMA “long”
  group_by(user_id) %>%
  summarise(
    ema_neg_aff_sd = sd(pid5_negative_affectivity, na.rm = TRUE),
    ema_det_sd = sd(pid5_detachment, na.rm = TRUE),
    ema_ant_sd = sd(pid5_antagonism, na.rm = TRUE),
    ema_dis_sd = sd(pid5_disinhibition, na.rm = TRUE),
    ema_psy_sd = sd(pid5_psychoticism, na.rm = TRUE),
    .groups = "drop"
  )

df_ml <- df2 %>% # 1 riga per soggetto
  distinct(
    user_id,
    esi_bf,
    domain_negative_affect,
    domain_detachment,
    domain_antagonism,
    domain_disinhibition,
    domain_psychoticism
  ) %>%
  left_join(ema_sd, by = "user_id")


#––– formula principale ––––––––––––––––––––––––––––––––––––––––––––––––––––––

form <-
  bf(
    esi_bf ~
      domain_negative_affect +
        domain_detachment +
        domain_antagonism +
        domain_disinhibition +
        domain_psychoticism +
        mi(ema_neg_aff_sd) +
        mi(ema_det_sd) +
        mi(ema_ant_sd) +
        mi(ema_dis_sd) +
        mi(ema_psy_sd),
    family = asym_laplace()
  ) +

  bf(ema_neg_aff_sd | mi() ~ 1, family = gaussian()) +
  bf(ema_det_sd | mi() ~ 1, family = gaussian()) +
  bf(ema_ant_sd | mi() ~ 1, family = gaussian()) +
  bf(ema_dis_sd | mi() ~ 1, family = gaussian()) +
  bf(ema_psy_sd | mi() ~ 1, family = gaussian()) +

  set_rescor(FALSE) # niente correlazioni residue fra gli “outcome”


fit_var <- brm(
  formula = form,
  data = df_ml, # una riga per participant
  chains = 4,
  iter = 4000,
  warmup = 1000,
  algorithm = control = list(adapt_delta = 0.95)
)

summary(fit_var)
bayes_R2(fit_var)


df2 <- df2 %>%
  group_by(user_id) %>%
  mutate(across(starts_with("pid5_"), ~ .x - mean(.x, na.rm = TRUE))) %>%
  ungroup()

df2 <- df2 %>%
  group_by(user_id) %>%
  mutate(
    across(
      starts_with("pid5_"),
      ~ {
        mu <- mean(.x, na.rm = TRUE)
        sd_ <- sd(.x, na.rm = TRUE)
        z <- (.x - mu) / sd_
        replace(z, is.na(z) | is.infinite(z), 0) # se sd = 0 o NA → 0
      }
    )
  ) %>%
  ungroup()

form <- bf(
  esi_bf ~
    domain_negative_affect +
      domain_detachment +
      domain_antagonism +
      domain_disinhibition +
      domain_psychoticism + # tratti (livello‑2)

      pid5_negative_affectivity +
      pid5_detachment +
      pid5_antagonism +
      pid5_disinhibition +
      pid5_psychoticism + # fluttuazioni EMA (livello‑1)

      (1 | user_id) # intercetta casuale
)

library(dplyr)

df2 <- df2 %>%
  # standardizzazione WITHIN‑person per le 5 variabili EMA ripetute
  group_by(user_id) %>%
  mutate(
    across(
      c(
        pid5_negative_affectivity,
        pid5_detachment,
        pid5_antagonism,
        pid5_disinhibition,
        pid5_psychoticism
      ),
      ~ {
        mu <- mean(.x, na.rm = TRUE)
        sd_ <- sd(.x, na.rm = TRUE)
        z <- (.x - mu) / sd_
        replace(z, is.na(z) | is.infinite(z), 0)
      }
    )
  ) %>%
  ungroup() %>%

  # standardizzazione GRAND‑mean per outcome + 5 domini trait
  mutate(
    across(
      c(
        esi_bf,
        domain_negative_affect,
        domain_detachment,
        domain_antagonism,
        domain_disinhibition,
        domain_psychoticism
      ),
      ~ as.numeric(scale(.x))
    )
  )


fit_w <- brm(
  formula = form,
  data = df2,
  family = skew_normal(),
  chains = 4,
  iter = 400000,
  warmup = 1000,
  # control = list(adapt_delta = .95),
  backend = "cmdstanr",
  algorithm = "meanfield"
)

summary(fit_w)

bayes_R2(fit_w)


ema_domains <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

# flag exam windows -------------------------------------------
stress_map <- df2 %>% # define once per study
  group_by(user_id) %>%
  mutate(
    stress_t = calendar_day >= 30 & calendar_day <= 37, # exam week
    poststress = calendar_day > 37 & calendar_day <= 44
  ) %>%
  ungroup()


# ────────────────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
# 1.  Helper to compute AR(1) safely  ───────────────────────────
safe_ar1 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  x1 <- head(x, -1)
  x2 <- tail(x, -1)
  if (sd(x1) == 0 || sd(x2) == 0) return(NA_real_)
  cor(x1, x2)
}

# ────────────────────────────────────────────────────────────────
# 2.  Build EMA features  (μ, σ, φ)  ────────────────────────────
ema_domains <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

stress_map <- df2 %>% # example cut‑offs
  group_by(user_id) %>%
  mutate(
    stress_t = calendar_day %in% 30:37,
    poststress = calendar_day %in% 38:44
  ) %>%
  ungroup()

ema_feats <- stress_map %>%
  group_by(user_id) %>%
  summarise(
    across(
      all_of(ema_domains),
      list(
        mu = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        phi = ~ safe_ar1(.x)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# ────────────────────────────────────────────────────────────────
# 3.  Baseline table + merge  ───────────────────────────────────
baseline <- df2 %>%
  select(
    user_id,
    esi_bf,
    domain_negative_affect,
    domain_detachment,
    domain_antagonism,
    domain_disinhibition,
    domain_psychoticism
  ) %>%
  distinct() %>%
  mutate(
    across(starts_with("domain_"), ~ replace_na(.x, mean(.x, na.rm = TRUE))),
    esibf = esi_bf
  ) # drop underscore

dat <- baseline %>%
  select(-esi_bf) %>% # use esibf only
  left_join(ema_feats, by = "user_id")

# ────────────────────────────────────────────────────────────────
# 4.  Drop *_beta,   rename predictors to remove '_'  ───────────
dyn_terms_old <- names(dat) %>% str_subset("_(mu|sd|phi)$")
dyn_terms_new <- str_replace_all(dyn_terms_old, "_", "")
names(dat)[match(dyn_terms_old, names(dat))] <- dyn_terms_new
dyn_terms <- dyn_terms_new

# ────────────────────────────────────────────────────────────────
# 5.  Z‑score all continuous predictors  ────────────────────────
scale_cols <- function(d, cols)
  d %>%
    mutate(across(all_of(cols), ~ as.vector(scale(.x)), .names = "{.col}_z"))
baseline_z <- paste0(
  c(
    "domain_negative_affect",
    "domain_detachment",
    "domain_antagonism",
    "domain_disinhibition",
    "domain_psychoticism"
  ),
  "_z"
)
dyn_terms_z <- paste0(dyn_terms, "_z")

dat_z <- dat %>%
  scale_cols(
    c(baseline_z, dyn_terms_z) %>% str_remove("_z$")
  )

# ────────────────────────────────────────────────────────────────
# 6.  Build brms formula with mi()  ─────────────────────────────
with_mi <- \(v) paste0("mi(", v, ")", collapse = " + ")

main <- bf(
  paste("esibf ~", with_mi(baseline_z), "+", with_mi(dyn_terms_z))
)

imp_forms <- lapply(
  c(baseline_z, dyn_terms_z),
  \(v) bf(paste0(v, " | mi() ~ 1"))
)

full_f <- purrr::reduce(imp_forms, `+`, .init = main) +
  set_rescor(FALSE)

# ────────────────────────────────────────────────────────────────
# 7.  Priors  (only on slopes of esibf)  ────────────────────────
priors <- prior(normal(0, 0.5), class = "b", resp = "esibf")


# ────────────────────────────────────────────────────────────────
# 8.  Short HMC using VB means as inits  ────────────────────────
fit_mcmc <- brm(
  full_f,
  data = dat_z,
  prior = priors,
  backend = "cmdstanr",
  chains = 2,
  cores = 2,
  iter = 2000,
  warmup = 1000,
  adapt_delta = 0.9,
  seed = 456
)

print(fit_mcmc, pars = "^b_esibf") # slopes for esibf only
bayes_R2(fit_mcmc)


# ─── already in memory: stress_map with stress_t / poststress flags ───
ema_domains <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

# ────────────────────────────────────────────────────────────────
# 1.  Build per‑domain μ, SD, φ  (as you did)  -------------------
pool_feats <- stress_map %>%
  group_by(user_id) %>%
  summarise(
    across(
      all_of(ema_domains),
      list(
        mu = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        phi = ~ safe_ar1(.x)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%

  # ─── ADD the pooled indices here  ──────────────────────────────
  mutate(
    mu_bar = rowMeans(select(., ends_with("_mu")), na.rm = TRUE),
    sd_bar = rowMeans(select(., ends_with("_sd")), na.rm = TRUE),
    phi_bar = rowMeans(select(., ends_with("_phi")), na.rm = TRUE)
  )

# 2.  Long table
long_dyn <- pool_feats %>%
  pivot_longer(
    cols = matches("_(mu|sd|phi)$"), # β dropped: all NA
    names_to = c("domain", "metric"),
    names_pattern = "^(.*)_(mu|sd|phi)$",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("mu", "sd", "phi")))

# 3.  Merge baseline + pooled indices (no suffixes)
long_dat <- long_dyn %>%
  left_join(baseline, by = "user_id") %>%
  left_join(
    pool_feats %>% select(user_id, mu_bar, sd_bar, phi_bar),
    by = "user_id",
    keep = FALSE
  )

# 4.  (optional) scale numeric columns
long_dat <- long_dat %>%
  mutate(across(where(is.numeric) & !c(user_id), ~ as.vector(scale(.x))))


# ─── after long_dat is ready  ──────────────────────────────────
with_mi <- \(v) paste0("mi(", v, ")", collapse = " + ")

full_f <- bf(
  paste(
    "esibf ~",
    # static PID‐5 domains (already z‑scaled + imputed)
    "mi(domain_negative_affect) + mi(domain_detachment) +",
    "mi(domain_antagonism) + mi(domain_disinhibition) +",
    "mi(domain_psychoticism) +",
    # pooled EMA indices actually present in long_dat
    "mi(mu_bar) + mi(sd_bar) + mi(phi_bar) +",
    # random deviation for each domain × metric
    "(0 + mi(value) | metric / domain)"
  )
) +
  set_rescor(FALSE)

# (no beta_bar now)

# ─── fit  ──────────────────────────────────────────────────────
fit_dyn <- brm(
  formula = full_f,
  data = long_dat,
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  adapt_delta = 0.95,
  backend = "cmdstanr",
  seed = 987
)


# 5.  Fit
fit_dyn <- brm(
  full_f,
  data = long_dat,
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  adapt_delta = 0.95,
  backend = "cmdstanr",
  seed = 987
)

################################

safe_ar1 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  x1 <- head(x, -1)
  x2 <- tail(x, -1)
  if (sd(x1) == 0 || sd(x2) == 0) return(NA_real_)
  cor(x1, x2)
}

#─────────────────────────────────────────────────────────────────
# 1.  Per‑participant EMA indices  -------------------------------
ema_domains <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

pool_feats <- stress_map %>% # <- you already have this tibble
  group_by(user_id) %>%
  summarise(
    across(
      all_of(ema_domains),
      list(
        mu = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        phi = ~ safe_ar1(.x)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    mu_bar = rowMeans(select(., ends_with("_mu")), na.rm = TRUE),
    sd_bar = rowMeans(select(., ends_with("_sd")), na.rm = TRUE),
    phi_bar = rowMeans(select(., ends_with("_phi")), na.rm = TRUE)
  ) %>%
  dplyr::select(user_id, mu_bar, sd_bar, phi_bar) # keep only what we need

# sanity
pool_feats %>% slice_head(n = 3)

#─────────────────────────────────────────────────────────────────
# 2.  Baseline static PID‑5 domains  -----------------------------
baseline <- df2 %>%
  dplyr::select(
    user_id,
    esi_bf,
    domain_negative_affect,
    domain_detachment,
    domain_antagonism,
    domain_disinhibition,
    domain_psychoticism
  ) %>%
  distinct() %>%
  # mean‑impute 64 missing domain scores
  mutate(across(
    starts_with("domain_"),
    ~ replace_na(.x, mean(.x, na.rm = TRUE))
  ))

#─────────────────────────────────────────────────────────────────
# 3.  Merge & z‑scale predictors  -------------------------------
dat <- baseline %>%
  left_join(pool_feats, by = "user_id") %>%
  # drop participants without EMA indices (rare)
  filter(!is.na(mu_bar)) %>%
  # z‑scale every numeric column *except* user_id
  mutate(across(where(is.numeric) & !c(user_id), ~ as.vector(scale(.x))))

# sanity
names(dat)
summary(dat$mu_bar)

#─────────────────────────────────────────────────────────────────
# 4.  Model formulas  -------------------------------------------
form_baseline <- bf(
  esi_bf ~
    domain_negative_affect +
      domain_detachment +
      domain_antagonism +
      domain_disinhibition +
      domain_psychoticism
)

# form_augmented <- bf(
#   esi_bf ~
#     domain_negative_affect +
#       domain_detachment +
#       domain_antagonism +
#       domain_disinhibition +
#       domain_psychoticism +
#       mu_bar +
#       sd_bar +
#       phi_bar
# )

form_augmented <-
  bf(
    esi_bf ~
      domain_negative_affect +
        domain_detachment +
        domain_antagonism +
        domain_disinhibition +
        domain_psychoticism +
        mi(mu_bar) +
        mi(sd_bar) +
        mi(phi_bar)
  ) +

  bf(mu_bar | mi() ~ 1) +
  bf(sd_bar | mi() ~ 1) +
  bf(phi_bar | mi() ~ 1) +
  set_rescor(FALSE)


#─────────────────────────────────────────────────────────────────
# 5.  Priors (only one line is really needed)  ------------------
priors <- prior(normal(0, 0.4), class = "b") # weak, but not huge

#─────────────────────────────────────────────────────────────────
# 6.  Fit both models  ------------------------------------------
fit0 <- brm(
  form_baseline,
  data = dat,
  prior = set_prior("horseshoe(0.5)", class = "b"),
  chains = 4,
  cores = 4,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.995, max_treedepth = 15), # key change
  iter = 2000,
  warmup = 1000,
  seed = 123
)

fit1 <- brm(
  form_augmented,
  data = dat,
  prior = priors,
  chains = 4,
  cores = 4,
  backend = "cmdstanr",
  iter = 2000,
  warmup = 1000,
  seed = 456
)

#─────────────────────────────────────────────────────────────────
# 7.  Compare predictive performance  ---------------------------

loo0 <- loo(fit0, resp = "esibf") # baseline model
loo1 <- loo(fit1, resp = "esibf") # augmented model
loo_compare(loo0, loo1)

bayes_R2(fit0)
bayes_R2(fit1)

#─────────────────────────────────────────────────────────────────
# 8.  Inspect EMA slopes  ---------------------------------------
print(fit1, pars = c("b_mimu_bar", "b_misd_bar", "b_miphi_bar"))


mood_vars <- c("happy", "sad", "satisfied", "angry")

mood_sd <- stress_map %>%
  group_by(user_id) %>%
  summarise(
    # rowwise = each participant separately
    rowwise() %>%
      mutate(
        mood_sd = mean(
          c_across(all_of(mood_vars)) |>
            sapply(sd, na.rm = TRUE),
          na.rm = TRUE
        )
      ) %>%
      ungroup() %>%
      distinct(user_id, mood_sd)
  )

mood_vars <- c("happy", "sad", "satisfied", "angry")
mood_sd <- stress_map %>%
  group_by(user_id) %>%
  summarise(
    mood_sd = mean(
      across(all_of(mood_vars), ~ sd(.x, na.rm = TRUE)),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

dat2 <- dat %>%
  left_join(mood_sd, by = "user_id") %>%
  mutate(mood_sd = as.vector(scale(mood_sd)))

fit2 <- brm(
  esi_bf ~ domain_ * +mi(mood_sd) + mi(sd_bar) + mi(phi_bar),
  data = dat2,
  prior = priors,
  chains = 4,
  cores = 4,
  backend = "cmdstanr"
)

print(fit2, pars = "b_mimood_sd") # expect a stronger positive slope


# 1. Replace plain SD with RMSSD (volatility of change)
pool_feats <- stress_map %>%
  group_by(user_id) %>%
  summarise(
    across(
      all_of(ema_domains),
      list(
        mu = ~ mean(.x, na.rm = TRUE),
        rmssd = ~ sqrt(mean(diff(.x)^2, na.rm = TRUE)),
        phi = ~ safe_ar1(.x)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    mu_bar = rowMeans(select(., ends_with("_mu")), na.rm = TRUE),
    rmssd_bar = rowMeans(select(., ends_with("_rmssd")), na.rm = TRUE),
    phi_bar = rowMeans(select(., ends_with("_phi")), na.rm = TRUE)
  )

# 2. Merge & z‑scale
dat <- baseline %>%
  left_join(pool_feats, by = "user_id") %>%
  mutate(across(where(is.numeric) & !c(user_id), ~ as.vector(scale(.))))

# 3. Fit augmented model with horseshoe prior
fit_hs <- brm(
  esi_bf ~
    domain_negative_affect +
      domain_detachment +
      domain_antagonism +
      domain_disinhibition +
      domain_psychoticism +
      mu_bar +
      rmssd_bar +
      phi_bar,
  data = dat,
  prior = set_prior("horseshoe(0.5)", class = "b"), # slightly stronger shrinkage
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1500,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.995, max_treedepth = 15), # key change
  seed = 789
)


print(fit_hs) # expect a stronger positive slope

bayes_R2(fit0)
bayes_R2(fit_hs)

dat_380 <- dat # already only the 380 ids with EMA
fit0_380 <- update(fit0, newdata = dat_380, recompile = FALSE)
loo0 <- loo(fit0_380, resp = "esi_bf")
loo1 <- loo(fit_hs, resp = "esi_bf")
loo_compare(loo0, loo1) # now legitimate

mood_vars <- c("happy", "sad", "satisfied", "angry")
mood_rmssd <- stress_map %>%
  group_by(user_id) %>%
  summarise(
    mood_rmssd = sqrt(mean(
      rowMeans(sapply(mood_vars, \(v) diff(.data[[v]])^2), na.rm = TRUE)
    )),
    .groups = "drop"
  )
dat2 <- dat_380 %>%
  left_join(mood_rmssd, by = "user_id") %>%
  mutate(mood_rmssd = scale(mood_rmssd)[, 1])


# °°°°°°°
vol_vars <- names(dyn_screen)[grep("_rmssd$", names(dyn_screen))]

dat_full <- dat %>%
  left_join(select(dyn_screen, user_id, all_of(vol_vars)), by = "user_id") %>%
  mutate(across(all_of(vol_vars), scale))

rhs <- paste(
  c(
    "domain_negative_affect",
    "domain_detachment",
    "domain_antagonism",
    "domain_disinhibition",
    "domain_psychoticism",
    "rmssd_bar",
    "phi_bar",
    vol_vars
  ),
  collapse = " + "
)

fit_full <- brm(
  bf(paste("esi_bf ~", rhs)),
  data = dat_full,
  prior = set_prior("horseshoe(0.5)", class = "b"),
  chains = 4,
  cores = 4,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.995)
)
