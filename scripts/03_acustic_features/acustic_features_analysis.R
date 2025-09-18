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
library(lubridate)
library(readxl)
library(janitor)
library(mice)
library(miceadds)
library(posterior)
library(tidybayes)
library(broom)
library(lmtest)
library(sandwich)
library(metafor)


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
quest_df <- left_join(esi_bf, pid5, by = "user_id")


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
quest_df <- quest_df[
  !(quest_df$user_id %in% user_id_with_careless_responding),
]

file_path <- here::here(
  "data",
  "raw",
  "acustic_features",
  "database_acfeat.xlsx"
)

# 1) Funzione per ripulire spazi (anche NBSP) e slash
fix_names <- function(nm) {
  nm %>%
    str_replace_all("[\\s\\u00A0]+", " ") %>% # comprime spazi e NBSP
    str_trim() %>% # toglie spazi ai bordi
    str_replace_all("\\s*/\\s*", "_") %>% # " / a / " -> "_a_"
    str_replace_all("/+", "_") %>% # eventuali "/" residui -> "_"
    str_squish()
}

wanted_sheets <- c("Baseline", "Pre", "Post")
sheets <- intersect(excel_sheets(file_path), wanted_sheets)

# 2) Import + normalizzazione nomi (fix + snake_case)
df1 <- map_dfr(
  sheets,
  ~ read_excel(file_path, sheet = .x) %>%
    {
      setNames(., fix_names(names(.)))
    } %>%
    clean_names(case = "snake") %>% # snake_case definitivo
    mutate(sheet = .x)
)

# 3) Separa Case in day e user_id (se serve)
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
    sheet = factor(sheet, levels = c("Baseline", "Pre", "Post"))
  )


df2 <- left_join(df1, quest_df, by = "user_id")


acoustic_feats <- c(
  "f0_min_hz_a",
  "f0_median_hz_i",
  "signalduration_s_i",
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

# check rapido: quali mancano?
setdiff(acoustic_feats, names(df2))


# colonne PID-5 (come da tuo glimpse)
pid5_cols <- c(
  "domain_negative_affect",
  "domain_detachment",
  "domain_antagonism",
  "domain_disinhibition",
  "domain_psychoticism"
)

# tabella a livello soggetto: una riga per user_id
pid_subject <- df2 %>%
  dplyr::select(user_id, esi_bf, all_of(pid5_cols)) %>%
  dplyr::distinct(user_id, .keep_all = TRUE)


set.seed(1234)

# costruisci lo schema di metodi: PMM per i 5 domini (sono punteggi interi ma PMM va benissimo)
meth <- make.method(pid_subject)
meth[pid5_cols] <- "pmm" # predittive mean matching
meth["esi_bf"] <- "pmm" # opzionale, se vuoi permettere l'imputazione di esi_bf

# matrice di predittori: non usare user_id come predittore
pred <- make.predictorMatrix(pid_subject)
pred[, "user_id"] <- 0
pred["user_id", ] <- 0

# avvia imputazione (m=20 consigliato)
imp_pid <- mice(
  pid_subject,
  m = 20,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  printFlag = TRUE
)


# ottieni i m data.frame imputati a livello soggetto
pid_imputed_list <- lapply(
  1:imp_pid$m,
  function(k) complete(imp_pid, action = k)
)

# (opzionale) arrotonda ai numeri interi i 5 domini
pid_imputed_list <- lapply(pid_imputed_list, function(d) {
  d %>%
    mutate(across(all_of(pid5_cols), ~ round(.x)))
})

# crea i m dataset analitici completi (join con df2)
analysis_list <- lapply(pid_imputed_list, function(pid_imp) {
  df2 %>%
    dplyr::select(-all_of(pid5_cols)) %>% # togli eventuali versioni con NA
    left_join(pid_imp, by = "user_id") # innesta i domini imputati
})

acoustic_feats <- c(
  "f0_min_hz_a",
  "f0_median_hz_i",
  # "signalduration_s_i",  # esclusa come richiesto
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

build_deltas <- function(dat) {
  ac_long <- dat %>%
    dplyr::select(
      user_id,
      sheet,
      dplyr::all_of(acoustic_feats),
      dplyr::all_of(pid5_cols),
      dplyr::any_of("esi_bf") # <-- stringa, non simbolo
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(acoustic_feats),
      names_to = "feature",
      values_to = "value"
    )

  deltas <- ac_long %>%
    dplyr::arrange(user_id, feature, sheet) %>%
    dplyr::group_by(user_id, feature) %>%
    dplyr::summarize(
      delta_pre_minus_base = if (all(c("Baseline", "Pre") %in% sheet)) {
        value[sheet == "Pre"] - value[sheet == "Baseline"]
      } else NA_real_,
      delta_post_minus_pre = if (all(c("Pre", "Post") %in% sheet)) {
        value[sheet == "Post"] - value[sheet == "Pre"]
      } else NA_real_,
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      c(delta_pre_minus_base, delta_post_minus_pre),
      names_to = "contrast",
      values_to = "delta"
    ) %>%
    dplyr::mutate(
      contrast = factor(
        contrast,
        levels = c("delta_pre_minus_base", "delta_post_minus_pre"),
        labels = c("Pre - Baseline", "Post - Pre")
      )
    ) %>%
    dplyr::left_join(
      dat %>%
        dplyr::select(
          user_id,
          dplyr::all_of(pid5_cols),
          dplyr::any_of("esi_bf")
        ) %>%
        dplyr::distinct(),
      by = "user_id"
    )

  deltas
}

fit_one_imputation <- function(deltas_df) {
  dd <- deltas_df %>%
    dplyr::mutate(
      dplyr::across(all_of(pid5_cols), ~ as.numeric(scale(.x))),
      pid5_composite = rowMeans(dplyr::across(all_of(pid5_cols)), na.rm = TRUE)
    ) %>%
    dplyr::filter(!is.na(delta))

  res <- dd %>%
    dplyr::group_by(feature, contrast) %>%
    dplyr::group_modify(
      ~ {
        fit <- stats::lm(delta ~ pid5_composite, data = .x)
        broom::tidy(fit)
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(term == "pid5_composite") %>%
    dplyr::select(feature, contrast, estimate, std.error)

  res
}

# STEP 1 — costruisci i delta per ciascuna imputazione
deltas_list <- lapply(analysis_list, build_deltas)

# controllo rapido: quante righe per la prima imputazione?
dplyr::glimpse(deltas_list[[1]])

# check: percentuale di NA nei delta (prima imputazione)
deltas_list[[1]] %>%
  dplyr::summarise(prop_na = mean(is.na(delta)))

# STEP 2 — stima su UNA imputazione (ad es. la prima)
ests_1 <- fit_one_imputation(deltas_list[[1]])
dplyr::glimpse(ests_1)

# (facoltativo) p-value "di controllo" per ordinare i risultati
ests_1 %>%
  dplyr::mutate(
    stat = estimate / std.error,
    p = 2 * pt(abs(stat), df = Inf, lower.tail = FALSE)
  ) %>%
  dplyr::arrange(p) %>%
  head(10)


# 3a) stima per ciascuna imputazione
ests_list <- lapply(deltas_list, fit_one_imputation)

# 3b) attacca l'indice di imputazione e unisci
ests_long <- dplyr::bind_rows(
  Map(
    function(df, k) dplyr::mutate(df, .imp = k),
    ests_list,
    seq_along(ests_list)
  )
)

# 3c) pooling Rubin per ogni (feature, contrast)
# helper per pooling Rubin
rubin_pool <- function(est, se) {
  m <- length(est)
  Qb <- mean(est) # Q-bar
  W <- mean(se^2) # within-imputation variance
  B <- stats::var(est) # between-imputation variance
  Tt <- W + (1 + 1 / m) * B # total variance
  se_pooled <- sqrt(Tt)

  # frazione d'incertezza da imputazione
  lambda <- if (Tt > 0) ((1 + 1 / m) * B) / Tt else 0
  df <- if (lambda > 0) (m - 1) / (lambda^2) else Inf

  tval <- Qb / se_pooled
  pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
  ci_low <- Qb - qt(.975, df = df) * se_pooled
  ci_high <- Qb + qt(.975, df = df) * se_pooled

  tibble::tibble(
    effect = Qb,
    se = se_pooled,
    df = df,
    t = tval,
    p = pval,
    ci_low = ci_low,
    ci_high = ci_high
  )
}

# 3c) pooling Rubin per ogni (feature, contrast)
pooled <- ests_long %>%
  dplyr::group_by(feature, contrast) %>%
  dplyr::summarise(
    rubin_pool(estimate, std.error),
    .groups = "drop"
  ) %>%
  dplyr::mutate(p_bh = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(contrast, p)

pooled %>% dplyr::slice_head(n = 12)

pooled %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice_min(order_by = p, n = 5) %>%
  dplyr::ungroup()

ggplot(
  pooled,
  aes(x = effect, y = feature, xmin = ci_low, xmax = ci_high, color = contrast)
) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~contrast, scales = "free_y") +
  labs(
    x = "Effetto del PID-5 composito sul delta acustico",
    y = "Feature acustica"
  ) +
  theme_minimal()


# ricostruisci un dataset long per mfcc9_iqr
example_df <- deltas_list[[1]] %>%
  dplyr::filter(feature == "mfcc9_iqr", contrast == "Pre - Baseline") %>%
  dplyr::mutate(
    pid5_composite = scale(rowMeans(
      dplyr::across(all_of(pid5_cols)),
      na.rm = TRUE
    ))
  )

ggplot(example_df, aes(x = pid5_composite, y = delta)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "PID-5 composito (z-score)", y = "Δ Pre - Baseline (mfcc9_iqr)") +
  theme_minimal()


pid5_cols <- c(
  "domain_negative_affect",
  "domain_detachment",
  "domain_antagonism",
  "domain_disinhibition",
  "domain_psychoticism"
)

prep_for_brms <- function(dd) {
  dd %>%
    # z-score dei domini e composito (media z)
    mutate(
      across(all_of(pid5_cols), ~ as.numeric(scale(.x))),
      pid5_composite = rowMeans(across(all_of(pid5_cols)), na.rm = TRUE)
    ) %>%
    dplyr::filter(!is.na(delta), !is.na(pid5_composite)) %>%
    mutate(
      delta_z = as.numeric(scale(delta)), # target standardizzato (stabilizza)
      feature = factor(feature),
      contrast = factor(contrast, levels = c("Pre - Baseline", "Post - Pre")),
      user_id = factor(user_id)
    )
}


bfml <- bf(
  delta_z ~
    1 +
      pid5_composite * contrast +
      (1 + pid5_composite | feature) +
      (1 | user_id)
)

priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sd"), # varianza effetti casuali
  prior(student_t(3, 0, 1), class = "sigma") # robustezza residui
)


# 1) Pooled dataset (media su imputazioni)
deltas_pooled <- purrr::map_dfr(deltas_list, prep_for_brms, .id = "imp")

deltas_mean <- deltas_pooled %>%
  group_by(user_id, feature, contrast) %>%
  summarise(
    delta_z = mean(delta_z, na.rm = TRUE),
    pid5_composite = mean(pid5_composite, na.rm = TRUE),
    .groups = "drop"
  )


# usa deltas_mean costruito prima
bfml_A <- bf(
  delta_z ~ 1 + pid5_composite * contrast + (1 | feature) + (1 | user_id)
)

priors_A <- c(
  prior(normal(0, 0.5), class = "Intercept"),
  prior(normal(0, 0.3), class = "b"),
  prior(exponential(2), class = "sd"), # RE più regolarizzati
  prior(exponential(2), class = "sigma")
)

fit_A <- brm(
  bfml_A,
  data = deltas_mean,
  family = gaussian(),
  prior = priors_A,
  chains = 4,
  iter = 3000,
  warmup = 1000,
  seed = 2025,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 0
)
summary(fit_A)


#### TODO -----------------

# quante imputazioni usare? (puoi mettere length(deltas_list))
M_to_fit <- length(deltas_list)

fits <- vector("list", M_to_fit)

for (k in seq_len(M_to_fit)) {
  dd_k <- prep_for_brms(deltas_list[[k]])
  fits[[k]] <- brm(
    formula = bfml,
    data = dd_k,
    family = student(),
    prior = priors,
    chains = 4,
    iter = 3000,
    warmup = 1000,
    seed = 2025 + k,
    cores = max(1, parallel::detectCores() - 1),
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    silent = 2,
    refresh = 0
  )
}


# estrai e concatena draws sulle stesse colonne
draws_list <- lapply(fits, as_draws_df)
common <- Reduce(intersect, lapply(draws_list, colnames))
draws <- do.call(rbind, lapply(draws_list, function(x) x[, common]))

# Effetto globale Pre−Baseline
glob_pre <- draws %>%
  transmute(beta_global_pre = b_pid5_composite)

# Effetto globale Post−Pre = main + interaction
glob_post <- draws %>%
  transmute(
    beta_global_post = b_pid5_composite + `b_pid5_composite:contrastPost - Pre`
  )

# Riassunti
summ_pre <- glob_pre %>% median_qi(beta_global_pre, .width = c(.80, .95))
summ_post <- glob_post %>% median_qi(beta_global_post, .width = c(.80, .95))

summ_pre
summ_post


#' Facciamo una PCA sui Δ (separatamente per Pre − Baseline e Post − Pre),
#' usando un dataset pooled (media sulle imputazioni), e poi testiamo la
#' relazione tra PC1 e PID-5 composito (lineare e, opzionale, quadratica).

# deltas_list: lista dei m dataset (già creata nei passaggi precedenti)
# pid5_cols: i 5 domini (già definito)

# unisci tutte le imputazioni, marcando l'indice
pooled0 <- purrr::map_dfr(
  seq_along(deltas_list),
  ~ dplyr::mutate(deltas_list[[.x]], .imp = .x)
)

# media del delta per (user_id, contrast, feature) sulle m imputazioni
deltas_mean <- pooled0 %>%
  group_by(user_id, contrast, feature) %>%
  summarise(delta = mean(delta, na.rm = TRUE), .groups = "drop")

# PID-5 composito pooled per soggetto (media sui m set, poi z-score)
pid5_mean <- pooled0 %>%
  group_by(.imp, user_id) %>%
  summarise(across(all_of(pid5_cols), ~ first(.x)), .groups = "drop") %>% # i domini sono costanti per riga (replicati su feature)
  group_by(user_id) %>%
  summarise(across(all_of(pid5_cols), mean), .groups = "drop") %>%
  mutate(
    pid5_composite = as.numeric(scale(rowMeans(
      across(all_of(pid5_cols)),
      na.rm = TRUE
    )))
  )


run_pca_and_reg <- function(contrast_label) {
  # matrice larga: colonne = feature, righe = soggetti
  wide_dat <- deltas_mean %>%
    dplyr::filter(contrast == contrast_label) %>%
    dplyr::select(user_id, feature, delta) %>%
    tidyr::pivot_wider(names_from = feature, values_from = delta)

  # rimuovi righe con NA
  wide_cc <- wide_dat %>% tidyr::drop_na()

  # separa X (numerica) e id
  X <- wide_cc %>%
    dplyr::select(-user_id) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.matrix()

  # PCA (centra e scala)
  pca_fit <- prcomp(X, center = TRUE, scale. = TRUE)

  # varianza spiegata
  var_prop <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)

  # PC1/PC2 per soggetto
  scores <- tibble::tibble(
    user_id = wide_cc$user_id,
    PC1 = pca_fit$x[, 1],
    PC2 = pca_fit$x[, 2]
  ) %>%
    dplyr::left_join(
      pid5_mean %>% dplyr::select(user_id, pid5_composite),
      by = "user_id"
    )

  # modelli: lineare e quadratico (screening)
  fit_lin <- stats::lm(PC1 ~ pid5_composite, data = scores)
  fit_quad <- stats::lm(
    PC1 ~ poly(pid5_composite, 2, raw = TRUE),
    data = scores
  )

  list(
    contrast = contrast_label,
    pca = pca_fit,
    var_prop = var_prop,
    scores = scores,
    fit_lin = broom::tidy(fit_lin),
    fit_glance = broom::glance(fit_lin),
    fit_quad = broom::tidy(fit_quad),
    quad_glance = broom::glance(fit_quad)
  )
}


res_pre <- run_pca_and_reg("Pre - Baseline")
res_post <- run_pca_and_reg("Post - Pre")

# Varianza spiegata dalle prime componenti
round(res_pre$var_prop[1:3], 3)
round(res_post$var_prop[1:3], 3)

# Effetto lineare PID-5 → PC1 (Pre−Baseline)
res_pre$fit_lin
res_pre$fit_glance %>% dplyr::select(r.squared, adj.r.squared, p.value)

# Effetto lineare PID-5 → PC1 (Post−Pre)
res_post$fit_lin
res_post$fit_glance %>% dplyr::select(r.squared, adj.r.squared, p.value)

# (Opzionale) Verifica curvatura
res_pre$fit_quad
res_pre$quad_glance %>% dplyr::select(r.squared, adj.r.squared, p.value)
res_post$fit_quad
res_post$quad_glance %>% dplyr::select(r.squared, adj.r.squared, p.value)


# helper: stima beta e SE per una singola feature
slope_one_feature <- function(df, robust = TRUE) {
  fit <- lm(delta ~ pid5_composite, data = df)
  if (robust) {
    V <- sandwich::vcovHC(fit, type = "HC3")
    ct <- lmtest::coeftest(fit, vcov = V)
    tibble(
      beta = unname(ct["pid5_composite", "Estimate"]),
      se = unname(ct["pid5_composite", "Std. Error"]),
      n = nobs(fit)
    )
  } else {
    s <- summary(fit)$coef["pid5_composite", ]
    tibble(
      beta = unname(s["Estimate"]),
      se = unname(s["Std. Error"]),
      n = nobs(fit)
    )
  }
}

# pendenze per tutte le feature, dentro un contrasto
slopes_by_contrast <- function(contrast_label, robust = TRUE) {
  deltas_mean %>%
    dplyr::filter(contrast == contrast_label) %>%
    left_join(
      pid5_mean %>% dplyr::select(user_id, pid5_composite),
      by = "user_id"
    ) %>%
    group_by(feature) %>%
    group_modify(~ slope_one_feature(.x, robust = robust)) %>%
    ungroup() %>%
    mutate(contrast = contrast_label)
}

slopes_pre <- slopes_by_contrast("Pre - Baseline", robust = TRUE)
slopes_post <- slopes_by_contrast("Post - Pre", robust = TRUE)

slopes_pre
slopes_post

### Non collasso le 5 dimensioni del PID-5 -------------------------------------

# Stima beta e SE robusto (HC3) per una feature con UN SOLO predittore "pred"
slope_one_feature_pred <- function(df, pred) {
  fml <- stats::as.formula(paste("delta ~", pred))
  fit <- stats::lm(fml, data = df)
  V <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov = V)
  tibble::tibble(
    beta = unname(ct[pred, "Estimate"]),
    se = unname(ct[pred, "Std. Error"]),
    n = stats::nobs(fit)
  )
}

# Pendenze per tutte le feature, per un dato CONTRASTO e un dato DOMINIO (pred)
slopes_by_contrast_and_domain <- function(contrast_label, domain_var) {
  deltas_mean %>%
    dplyr::filter(contrast == contrast_label) %>%
    dplyr::left_join(
      pid5_mean %>% dplyr::select(user_id, dplyr::all_of(domain_var)),
      by = "user_id"
    ) %>%
    dplyr::rename(pred = !!domain_var) %>% # pred sarà il nome usato nella formula
    dplyr::group_by(feature) %>%
    dplyr::group_modify(~ slope_one_feature_pred(.x, "pred")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      contrast = contrast_label,
      domain = domain_var
    )
}

# calcola pendenze per tutti i domini e contrasti
slopes_all <- purrr::map_dfr(
  pid5_cols,
  ~ bind_rows(
    slopes_by_contrast_and_domain("Pre - Baseline", .x),
    slopes_by_contrast_and_domain("Post - Pre", .x)
  )
)

# occhiata sintetica: media (non pesata) e proporzione di segni per dominio/contrasto
slopes_all %>%
  group_by(domain, contrast) %>%
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    prop_pos = mean(beta > 0, na.rm = TRUE),
    k = dplyr::n(),
    .groups = "drop"
  )


meta_one <- function(df_dom_con) {
  # yi = stima; sei = errore standard; RE model (DL come default)
  res <- metafor::rma.uni(
    yi = beta,
    sei = se,
    data = df_dom_con,
    method = "REML"
  )
  tibble::tibble(
    k = res$k,
    mu_hat = as.numeric(res$b), # effetto medio pooled
    se_mu = res$se,
    ci_lb = res$ci.lb,
    ci_ub = res$ci.ub,
    tau2 = res$tau2, # eterogeneità
    Q = res$QE, # test eterogeneità
    Q_p = res$QEp,
    p_mu = res$pval # p-value su mu_hat
  )
}

meta_table <- slopes_all %>%
  group_by(domain, contrast) %>%
  group_modify(~ meta_one(.x)) %>%
  ungroup() %>%
  arrange(contrast, domain)

meta_table


forest_one <- function(domain_var, contrast_label) {
  df <- slopes_all %>%
    dplyr::filter(domain == domain_var, contrast == contrast_label) %>%
    mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se)

  meta_row <- meta_table %>%
    filter(domain == domain_var, contrast == contrast_label) %>%
    mutate(feature = "Pooled (REML)", beta = mu_hat, lo = ci_lb, hi = ci_ub)

  bind_rows(
    df %>% mutate(feature = as.character(feature)),
    meta_row %>% dplyr::select(feature, beta, lo, hi)
  ) %>%
    mutate(feature = factor(feature, levels = rev(unique(feature)))) %>%
    ggplot(aes(x = beta, y = feature)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste0("Forest: ", domain_var, " — ", contrast_label),
      x = "Pendenza (Δ per 1 unità dominio PID-5)",
      y = NULL
    ) +
    theme_minimal()
}

# esempi:
forest_one("domain_negative_affect", "Pre - Baseline")
forest_one("domain_negative_affect", "Post - Pre")


# pendenze "parziali" da modello multivariato per feature: delta ~ 5 domini
slope_multivar_one_feature <- function(df) {
  fit <- lm(
    delta ~
      domain_negative_affect +
        domain_detachment +
        domain_antagonism +
        domain_disinhibition +
        domain_psychoticism,
    data = df
  )
  V <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov = V)
  tibble::tibble(
    term = rownames(ct),
    beta = ct[, "Estimate"],
    se = ct[, "Std. Error"]
  ) %>%
    dplyr::filter(term != "(Intercept)")
}

# wrapper per contrasto: restituisce beta,se per OGNI dominio e feature
slopes_multivar_by_contrast <- function(contrast_label) {
  deltas_mean %>%
    dplyr::filter(contrast == contrast_label) %>%
    left_join(
      pid5_mean %>% dplyr::select(user_id, all_of(pid5_cols)),
      by = "user_id"
    ) %>%
    group_by(feature) %>%
    group_modify(~ slope_multivar_one_feature(.x)) %>%
    ungroup() %>%
    mutate(
      domain = term,
      contrast = contrast_label
    ) %>%
    dplyr::select(feature, domain, beta, se, contrast)
}

slopes_mv_pre <- slopes_multivar_by_contrast("Pre - Baseline")
slopes_mv_post <- slopes_multivar_by_contrast("Post - Pre")

slopes_mv_pre |> as.data.frame()
slopes_mv_post |> as.data.frame()


# 1) Combina i due contrasti
slopes_mv_both <- bind_rows(slopes_mv_pre, slopes_mv_post) %>%
  mutate(
    # etichetta più compatta per i domini (opzionale)
    domain_short = recode(
      domain,
      "domain_negative_affect" = "NegAff",
      "domain_detachment" = "Detach",
      "domain_antagonism" = "Antag",
      "domain_disinhibition" = "Disinh",
      "domain_psychoticism" = "Psych"
    ),
    # 2) Classifica le feature in famiglie
    family = case_when(
      str_detect(feature, "^f0_") ~ "F0",
      str_detect(feature, "^f1_|^f2_|^f3_") ~ "Formants",
      str_detect(feature, "^mfcc") ~ "MFCC",
      TRUE ~ "Other"
    )
  )


# meta per un singolo data.frame (una cella: es. un dominio in un contrasto)
meta_one <- function(df) {
  # rma.uni calcola anche I2 e test di eterogeneità (QE)
  fit <- rma.uni(yi = beta, sei = se, data = df, method = "REML")
  tibble(
    k = fit$k,
    mu_hat = as.numeric(fit$b),
    se_mu = fit$se,
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    tau2 = fit$tau2,
    I2 = fit$I2,
    Q = fit$QE,
    Q_p = fit$QEp,
    p_mu = fit$pval
  )
}

# 1) Pooled per dominio × contrasto
meta_domain <- slopes_mv_both %>%
  group_by(contrast, domain_short, domain) %>%
  group_modify(~ meta_one(.x)) %>%
  ungroup() %>%
  arrange(contrast, domain_short)

# 2) Pooled per dominio × contrasto × famiglia (F0, Formants, MFCC)
meta_domain_family <- slopes_mv_both %>%
  group_by(contrast, domain_short, domain, family) %>%
  group_modify(~ meta_one(.x)) %>%
  ungroup() %>%
  arrange(contrast, domain_short, family)

# stampa sintetica
meta_domain %>%
  dplyr::select(contrast, domain_short, k, mu_hat, ci_lb, ci_ub, p_mu, tau2, I2)
meta_domain_family %>%
  dplyr::select(
    contrast,
    domain_short,
    family,
    k,
    mu_hat,
    ci_lb,
    ci_ub,
    p_mu,
    tau2,
    I2
  )


forest_one <- function(domain_short_label, contrast_label) {
  df <- slopes_mv_both %>%
    dplyr::filter(
      domain_short == domain_short_label,
      contrast == contrast_label
    ) %>%
    mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se)

  pooled_row <- meta_domain %>%
    dplyr::filter(
      domain_short == domain_short_label,
      contrast == contrast_label
    ) %>%
    transmute(feature = "Pooled (REML)", beta = mu_hat, lo = ci_lb, hi = ci_ub)

  bind_rows(
    df %>% dplyr::select(feature, beta, lo, hi),
    pooled_row
  ) %>%
    mutate(feature = factor(feature, levels = rev(unique(feature)))) %>%
    ggplot(aes(x = beta, y = feature)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste0("Forest — ", domain_short_label, " — ", contrast_label),
      x = "Pendenza (Δ per 1 unità del dominio PID-5)",
      y = NULL
    ) +
    theme_minimal(base_size = 12)
}

# esempi:
forest_one("NegAff", "Pre - Baseline")
forest_one("NegAff", "Post - Pre")


forest_by_family <- function(domain_short_label, contrast_label) {
  df <- slopes_mv_both %>%
    dplyr::filter(
      domain_short == domain_short_label,
      contrast == contrast_label
    ) %>%
    mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se)

  pooled_rows <- meta_domain_family %>%
    dplyr::filter(
      domain_short == domain_short_label,
      contrast == contrast_label
    ) %>%
    transmute(
      feature = paste0("Pooled (", family, ")"),
      beta = mu_hat,
      lo = ci_lb,
      hi = ci_ub
    )

  bind_rows(
    df %>% dplyr::select(feature, beta, lo, hi),
    pooled_rows
  ) %>%
    mutate(feature = factor(feature, levels = rev(unique(feature)))) %>%
    ggplot(aes(x = beta, y = feature)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste0(
        "Forest — ",
        domain_short_label,
        " — ",
        contrast_label,
        " (per famiglia)"
      ),
      x = "Pendenza (Δ per 1 unità del dominio PID-5)",
      y = NULL
    ) +
    theme_minimal(base_size = 12)
}

# esempio:
forest_by_family("NegAff", "Pre - Baseline")
forest_by_family("Detach", "Pre - Baseline")
forest_by_family("Antag", "Pre - Baseline")
forest_by_family("Disinh", "Pre - Baseline")
forest_by_family("Psych", "Pre - Baseline")


nice_tbl <- meta_domain %>%
  mutate(
    estimate_CI = sprintf("%.3f [%.3f, %.3f]", mu_hat, ci_lb, ci_ub),
    I2_fmt = sprintf("%.1f%%", I2),
    p_fmt = format.pval(p_mu, digits = 3, eps = .001)
  ) %>%
  dplyr::select(
    contrast,
    domain = domain_short,
    k,
    estimate_CI,
    p = p_fmt,
    tau2,
    I2 = I2_fmt
  )

nice_tbl


# Filtra solo Post - Pre e i due domini interessanti
df_focus <- slopes_mv_post %>%
  dplyr::filter(domain %in% c("domain_antagonism", "domain_disinhibition"))

# Uniforma i nomi
df_focus <- df_focus %>%
  mutate(
    domain_clean = recode(
      domain,
      "domain_antagonism" = "Antagonism",
      "domain_disinhibition" = "Disinhibition"
    )
  )

# Funzione per meta-analisi per singolo dominio
run_meta <- function(dat, dom) {
  res <- rma(
    yi = beta,
    sei = se,
    data = filter(dat, domain_clean == dom),
    method = "REML"
  )
  list(domain = dom, res = res)
}

meta_antag <- run_meta(df_focus, "Antagonism")
meta_disinh <- run_meta(df_focus, "Disinhibition")

summary(meta_antag$res)
summary(meta_disinh$res)

forest(
  meta_antag$res,
  slab = df_focus$feature[df_focus$domain_clean == "Antagonism"],
  xlab = "Beta",
  main = "Post - Pre: Antagonism"
)

forest(
  meta_disinh$res,
  slab = df_focus$feature[df_focus$domain_clean == "Disinhibition"],
  xlab = "Beta",
  main = "Post - Pre: Disinhibition"
)


df_focus |>
  dplyr::filter(
    contrast == "Post - Pre",
    domain %in% c("Antagonism", "Disinhibition")
  ) |>
  arrange(beta)
