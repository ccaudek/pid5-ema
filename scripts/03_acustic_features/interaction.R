# ============================================================
# Δ stress nel lungo: phase × PID-5 con random (1 + phase | user_id)
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(loo)
  library(tidybayes)
  library(posterior)
  library(stringr)
  library(conflicted)
  options(mc.cores = parallel::detectCores())
})

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("var", "stats")
conflict_prefer("sd", "stats")
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
conflict_prefer("pstudent_t", "brms")
conflict_prefer("qstudent_t", "brms")
conflict_prefer("rstudent_t", "brms")
conflict_prefer("lmer", "lme4")

set.seed(123)

# ------------------------------------------------------------
# 1) Lettura e preparazione dati
# ------------------------------------------------------------
dat_long <- readr::read_csv(
  "acoustic_pid5_long_avg.csv",
  show_col_types = FALSE
)

# teniamo solo pre/post per testare l'effetto dello stress
dat_pp <- dat_long %>%
  filter(phase %in% c("pre", "post")) %>%
  mutate(
    phase = factor(phase, levels = c("pre", "post")) # ref = pre
  )

# soggetti con entrambe le fasi
have_prepost <- dat_pp %>%
  count(user_id, phase) %>%
  tidyr::pivot_wider(names_from = phase, values_from = n, values_fill = 0) %>%
  filter(pre > 0 & post > 0) %>%
  pull(user_id)

dat_pp <- dat_pp %>% filter(user_id %in% have_prepost)
length(unique(dat_pp$user_id))

# PID-5 su scala z disponibili nel file lungo:
pid5_z <- c(
  "domain_negative_affect_z",
  "domain_detachment_z",
  "domain_antagonism_z",
  "domain_disinhibition_z",
  "domain_psychoticism_z"
)
stopifnot(all(pid5_z %in% names(dat_pp)))

# 4 feature MFCC di interesse (livello del segnale, non Δ):
feat_list <- c("mfcc3_std", "mfcc3_iqr", "mfcc9_std", "mfcc9_iqr")
missing_feat <- setdiff(feat_list, names(dat_pp))
if (length(missing_feat) > 0) {
  stop(
    "Mancano queste feature in acoustic_pid5_long_avg.csv: ",
    paste(missing_feat, collapse = ", ")
  )
}


# ================== MODELLO SINGOLO ==================

# Scegli la feature e il dominio PID-5
# yvar   <- "mfcc3_std"                 # <- cambia qui se vuoi mfcc3_iqr / mfcc9_std / mfcc9_iqr
# pidvar <- "domain_disinhibition_z"    # <- cambia qui il dominio PID-5

yvar <- "mfcc3_iqr"
pidvar <- "domain_antagonism_z"

# 1) Costruisci il sotto-dataset per questo modello
d <- dat_pp %>%
  transmute(
    user_id,
    phase = factor(phase, levels = c("pre", "post")), # ref = pre
    y = as.numeric(scale(.data[[yvar]])), # z-score della feature
    pid_z = .data[[pidvar]]
  ) %>%
  drop_na(user_id, phase, y, pid_z)

# Check minimi
stopifnot(length(unique(d$phase)) == 2)
if (sd(d$pid_z) == 0 || is.na(sd(d$pid_z)))
  stop("Varianza nulla in pid_z: impossibile stimare interazione.")
# almeno 2 righe per id (pre+post)
if ((d %>% count(user_id) %>% summarise(min(n)) %>% pull()) < 2) {
  stop("Alcuni id hanno una sola fase osservata in questa feature.")
}

# 2) Specifica formule
f_main <- bf(y ~ phase + pid_z + (1 + phase | user_id))
f_int <- bf(y ~ phase * pid_z + (1 + phase | user_id))

# 3) Priors "sicuri" via get_prior (evita duplicati)
make_priors <- function(formula, data, family = student()) {
  pri_tbl <- get_prior(formula = formula, data = data, family = family)
  pri_tbl <- pri_tbl %>%
    mutate(
      prior = case_when(
        class == "b" ~ "normal(0, 0.3)", # slope piccoli, regolarizzanti
        class == "Intercept" ~ "student_t(3, 0, 1)",
        class == "sigma" ~ "student_t(3, 0, 1)",
        class == "nu" ~ "gamma(2, 0.1)", # solo se Student-t
        class == "sd" ~ "exponential(1)", # RE SD
        class == "cor" ~ "lkj(2)", # correlazioni RE
        TRUE ~ prior
      )
    )
  do.call(
    c,
    apply(pri_tbl, 1, function(r) {
      args <- list(prior = r[["prior"]], class = r[["class"]])
      if (nzchar(r[["resp"]])) args$resp <- r[["resp"]]
      if (nzchar(r[["coef"]])) args$coef <- r[["coef"]]
      if (nzchar(r[["group"]])) args$group <- r[["group"]]
      do.call(prior, args)
    })
  )
}

pri_main <- make_priors(f_main, d, family = student())
pri_int <- make_priors(f_int, d, family = student())

# 4) Fit (prima interaction, poi main per confronto LOO)
m_int <- brm(
  formula = f_int,
  data = d,
  family = student(),
  prior = pri_int,
  chains = 4,
  iter = 4000,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  backend = "cmdstanr"
)

m_main <- brm(
  formula = f_main,
  data = d,
  family = student(),
  prior = pri_main,
  chains = 4,
  iter = 4000,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  backend = "cmdstanr"
)

# 5) LOO: interaction vs main-effects
m_int <- add_criterion(m_int, "loo")
m_main <- add_criterion(m_main, "loo")
print(loo::loo_compare(m_main$criteria$loo, m_int$criteria$loo))

# 6) Estrai l'interazione "phasepost:pid_z" in modo sicuro
coef_name <- "b_phasepost:pid_z"
dr <- posterior::as_draws(m_int)
if (!coef_name %in% posterior::variables(dr)) {
  stop(
    "Il coefficiente di interazione non è presente. Controlla livelli di phase e varianza di pid_z."
  )
}
dr_sub <- posterior::subset_draws(dr, variables = coef_name)
est_tbl <- posterior::summarise_draws(dr_sub) %>% as_tibble()
v <- posterior::as_draws_df(dr_sub)[[coef_name]]
pd_val <- max(mean(v > 0), mean(v < 0))

cat("\n=== Interazione phase:PID ===\n")
cat("feature:", yvar, "  pid:", pidvar, "\n")
print(est_tbl %>% dplyr::select(mean, q5, q95, sd))
cat("PD (probabilità direzionale) =", round(pd_val, 3), "\n\n")

# 7) Posterior predictive check e effetti marginali
pp_check(m_int, type = "dens_overlay")

# Effetti condizionali: come cambia y (z) a -2..+2 SD di pid_z per pre/post
print(conditional_effects(m_int, effects = "pid_z:phase"))

# 8) Contrasto post-pre a pid_z basso/alto (es. -1 e +1 SD)
newd <- tidyr::expand_grid(
  phase = factor(c("pre", "post"), levels = c("pre", "post")),
  pid_z = c(-1, 1)
)
pred <- fitted(m_int, newdata = newd, re_formula = NA, summary = TRUE) %>%
  as_tibble() %>%
  bind_cols(newd)

# calcola Δ(post-pre) a pid_z = -1 e +1
delta_tbl <- pred %>%
  select(phase, pid_z, Estimate) %>%
  pivot_wider(names_from = phase, values_from = Estimate) %>%
  mutate(delta_post_pre = post - pre)

cat("\n=== Δ(post-pre) previsto a pid_z = -1 e +1 ===\n")
print(delta_tbl)


##---------------------------------------------

# SCEGLI la coppia
yvar <- "mfcc3_std"
pidvar <- "domain_disinhibition_z"

# Dati per il modello (come prima)
d <- dat_pp %>%
  transmute(
    user_id,
    phase = factor(phase, levels = c("pre", "post")),
    y = as.numeric(scale(.data[[yvar]])),
    pid_z = .data[[pidvar]]
  ) %>%
  drop_na(user_id, phase, y, pid_z)
stopifnot(length(unique(d$phase)) == 2, sd(d$pid_z) > 0)

# Formula distribuzionale:
# - media:   y ~ phase * pid_z + (1 + phase | user_id)
# - scala σ: sigma ~ phase * pid_z     (log-link implicito)
f_dist <- bf(
  y ~ phase * pid_z + (1 + phase | user_id),
  sigma ~ phase * pid_z
)

mkpri_strong <- function(formula, data, family = student()) {
  tab <- get_prior(formula, data = data, family = family) |>
    mutate(
      prior = case_when(
        class == "b" ~ "normal(0, 0.25)", # slope più regolarizzati
        class == "Intercept" ~ "student_t(3, 0, 1)",
        class == "sigma" ~ "student_t(3, 0, 1)",
        class == "nu" ~ "gamma(2, 0.1)",
        class == "sd" & coef == "phasepost" ~ "exponential(2)", # shrink extra sullo slope di phase
        class == "sd" ~ "exponential(1.5)",
        class == "cor" ~ "lkj(2)",
        TRUE ~ prior
      )
    )
  do.call(
    c,
    apply(tab, 1, \(r) {
      args <- list(prior = r[["prior"]], class = r[["class"]])
      if (nzchar(r[["dpar"]])) args$dpar <- r[["dpar"]]
      if (nzchar(r[["coef"]])) args$coef <- r[["coef"]]
      if (nzchar(r[["group"]])) args$group <- r[["group"]]
      do.call(prior, args)
    })
  )
}

f_dist <- bf(y ~ phase * pid_z + (1 + phase | user_id), sigma ~ phase * pid_z)

pri_str <- mkpri_strong(f_dist, d, student())

fit_dist2 <- brm(
  f_dist,
  data = d,
  family = student(),
  prior = pri_str,
  chains = 4,
  iter = 6000,
  warmup = 3000,
  seed = 123,
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  backend = "cmdstanr"
)
summary(fit_dist2)


# Interazione su mu
coef_mu <- "b_phasepost:pid_z"
dr <- posterior::as_draws(fit_dist2)
v_mu <- posterior::as_draws_df(posterior::subset_draws(
  dr,
  variables = coef_mu
))[[coef_mu]]
est_mu <- tibble(
  mean = mean(v_mu),
  sd = sd(v_mu),
  q2.5 = quantile(v_mu, 0.025),
  q97.5 = quantile(v_mu, 0.975),
  pd = max(mean(v_mu > 0), mean(v_mu < 0))
)
est_mu

# Interazione su sigma (ricorda: è su log(sigma))
coef_sigma <- "b_sigma_phasepost:pid_z"
v_si <- posterior::as_draws_df(posterior::subset_draws(
  dr,
  variables = coef_sigma
))[[coef_sigma]]
est_si <- tibble(
  mean = mean(v_si),
  sd = sd(v_si),
  q2.5 = quantile(v_si, 0.025),
  q97.5 = quantile(v_si, 0.975),
  pd = max(mean(v_si > 0), mean(v_si < 0))
)
est_si


# stesse 'd', yvar, pidvar usati per fit_dist2
f_fullRE <- bf(y ~ phase * pid_z + (1 + phase | user_id), sigma ~ phase * pid_z)
f_intRE <- bf(y ~ phase * pid_z + (1 | user_id), sigma ~ phase * pid_z)

pri_full <- mkpri_strong(f_fullRE, d, student()) # la funzione mkpri_strong che hai già definito
pri_int <- mkpri_strong(f_intRE, d, student())

fit_fullRE <- brm(
  f_fullRE,
  data = d,
  family = student(),
  prior = pri_full,
  chains = 4,
  iter = 6000,
  warmup = 3000,
  seed = 123,
  control = list(adapt_delta = .995, max_treedepth = 15),
  backend = "cmdstanr"
) |>
  add_criterion("loo")

fit_intRE <- brm(
  f_intRE,
  data = d,
  family = student(),
  prior = pri_int,
  chains = 4,
  iter = 6000,
  warmup = 3000,
  seed = 123,
  control = list(adapt_delta = .995, max_treedepth = 15),
  backend = "cmdstanr"
) |>
  add_criterion("loo")

loo::loo_compare(fit_intRE$criteria$loo, fit_fullRE$criteria$loo)

###########

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(performance)

# scegli la coppia
yvar <- "mfcc3_std" # oppure mfcc3_iqr, mfcc9_std, mfcc9_iqr
pidvar <- "domain_disinhibition_z" # oppure un altro dominio PID-5

# dati (come già preparato: dat_pp = pre/post, PID z)
d <- dat_pp %>%
  transmute(
    user_id,
    phase = factor(phase, levels = c("pre", "post")),
    y = as.numeric(scale(.data[[yvar]])),
    pid_z = .data[[pidvar]]
  ) %>%
  drop_na()

# modello: interazione sulla media, random intercept per soggetto
m_lmer <- lmer(y ~ phase * pid_z + (1 | user_id), data = d, REML = TRUE)

cat("\n=== Riassunto (lmerTest) ===\n")
print(summary(m_lmer)) # stime + p-val Satterthwaite
cat("\n=== ANOVA Type III ===\n")
print(anova(m_lmer, type = 3)) # F-test per l'interazione

# CI 95% per l'interazione
cat("\nCI 95% (Wald) per phasepost:pid_z\n")
print(confint(m_lmer, parm = "phasepost:pid_z", method = "Wald"))

# R2 di Nakagawa
cat("\nR2 Nakagawa (marginale/condizionato):\n")
print(performance::r2(m_lmer))

# Contrasto Δ(post-pre) a pid_z = -1 e +1 SD
em <- emmeans(m_lmer, ~ phase | pid_z, at = list(pid_z = c(-1, 1)))
cat("\n=== EMM per fase a pid_z = -1 e +1 ===\n")
print(em)
cat("\n=== Δ(post - pre) a pid_z = -1 e +1 ===\n")
print(contrast(em, method = list("post - pre" = c(-1, 1))))

# Tabella compatta della sola interazione
cat("\n=== Interazione (tabella compatta) ===\n")
print(
  broom.mixed::tidy(m_lmer, effects = "fixed") %>%
    dplyr::filter(term == "phasepost:pid_z") %>%
    dplyr::select(term, estimate, std.error, statistic)
)


##########

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(performance)

fit_lmer_one <- function(
  dat_pp,
  yvar,
  pidvar,
  pid_grid = seq(-2, 2, by = 0.25)
) {
  d <- dat_pp %>%
    transmute(
      user_id,
      phase = factor(phase, levels = c("pre", "post")),
      y = as.numeric(scale(.data[[yvar]])),
      pid_z = .data[[pidvar]]
    ) %>%
    drop_na()
  stopifnot(length(unique(d$phase)) == 2, sd(d$pid_z) > 0)

  m <- lmer(y ~ phase * pid_z + (1 | user_id), data = d, REML = TRUE)

  # Coefficienti fissi (con p di lmerTest)
  fixed_tab <- broom.mixed::tidy(m, effects = "fixed")

  # CI 95% (Wald) per l'interazione
  ci_int <- confint(m, parm = "phasepost:pid_z", method = "Wald")

  # ANOVA (Kenward–Roger via lmerTest)
  anova3 <- anova(m)

  # *** SIMPLE SLOPES ***  (fix: serve 'var = "pid_z"')
  slopes <- emtrends(m, specs = ~phase, var = "pid_z")
  slopes_diff <- contrast(slopes, method = "revpairwise") # (slope_post - slope_pre)

  # Δ(post-pre) su una griglia di PID (per visualizzare l'interazione)
  emm_grid <- emmeans(m, ~ phase | pid_z, at = list(pid_z = pid_grid))
  delta_grid <- contrast(emm_grid, method = list("post - pre" = c(-1, 1)))

  # Riassunto compatto dell'interazione
  inter_row <- fixed_tab %>%
    filter(term == "phasepost:pid_z") %>%
    transmute(term, estimate, std.error, statistic)

  list(
    model = m,
    fixed = fixed_tab,
    anova = anova3,
    ci_interaction = ci_int,
    slopes = slopes,
    slopes_diff = slopes_diff,
    delta_grid = delta_grid,
    interaction_row = inter_row
  )
}


out1 <- fit_lmer_one(
  dat_pp,
  yvar = "mfcc3_std",
  pidvar = "domain_disinhibition_z"
)
out1$interaction_row # β_interazione, SE, t, p
out1$slopes # slope di PID in PRE e in POST
out1$slopes_diff # differenza tra slope (≈ interazione)
as_tibble(out1$delta_grid) %>% slice(1:5) # Δ(post-pre) lungo PID
plot(out1$delta_grid) # CI dell’effetto post-pre come funzione di PID_z


out1 <- fit_lmer_one(
  dat_pp,
  yvar = "mfcc3_iqr",
  pidvar = "domain_disinhibition_z"
)
out1$interaction_row # β_interazione, SE, t, p
out1$slopes # slope di PID in PRE e in POST
out1$slopes_diff # differenza tra slope (≈ interazione)
as_tibble(out1$delta_grid) %>% slice(1:5) # Δ(post-pre) lungo PID
plot(out1$delta_grid) # CI dell’effetto post-pre come funzione di PID_z


out1 <- fit_lmer_one(
  dat_pp,
  yvar = "mfcc9_std",
  pidvar = "domain_disinhibition_z"
)
out1$interaction_row # β_interazione, SE, t, p
out1$slopes # slope di PID in PRE e in POST
out1$slopes_diff # differenza tra slope (≈ interazione)
as_tibble(out1$delta_grid) %>% slice(1:5) # Δ(post-pre) lungo PID
plot(out1$delta_grid) # CI dell’effetto post-pre come funzione di PID_z


out1 <- fit_lmer_one(
  dat_pp,
  yvar = "mfcc9_iqr",
  pidvar = "domain_disinhibition_z"
)
out1$interaction_row # β_interazione, SE, t, p
out1$slopes # slope di PID in PRE e in POST
out1$slopes_diff # differenza tra slope (≈ interazione)
as_tibble(out1$delta_grid) %>% slice(1:5) # Δ(post-pre) lungo PID
plot(out1$delta_grid) # CI dell’effetto post-pre come funzione di PID_z


############## composito

# ===== 1) Composite MFCC (PC1) =====
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

dat_long <- readr::read_csv(
  "acoustic_pid5_long_avg.csv",
  show_col_types = FALSE
)

mfcc_feats <- c("mfcc3_std", "mfcc3_iqr", "mfcc9_std", "mfcc9_iqr")
stopifnot(all(mfcc_feats %in% names(dat_long)))

# Z-score delle MFCC (sul lungo) e PCA su complete cases delle 4 MFCC
X <- dat_long %>%
  select(all_of(mfcc_feats)) %>%
  mutate(across(everything(), scale)) %>%
  as_tibble()

cc <- complete.cases(X)
pca_fit <- prcomp(X[cc, ], center = FALSE, scale. = FALSE) # già z-scored

# PC1 scores per le righe complete; metti NA altrimenti
pc1 <- rep(NA_real_, nrow(X))
pc1[cc] <- as.numeric(predict(pca_fit, newdata = X[cc, ])[, "PC1"])

# Orienta il segno di PC1 verso cariche positive su mfcc3_iqr (più “recupero” = più alto)
sgn <- sign(cor(pc1[cc], X$mfcc3_iqr[cc], use = "pairwise"))
pc1 <- sgn * pc1

dat_long <- dat_long %>%
  mutate(mfcc_pc1 = as.numeric(scale(pc1)))

# Salva dataset con PC1
readr::write_csv(dat_long, "acoustic_pid5_long_with_pc1.csv")
message(">> Salvato: acoustic_pid5_long_with_pc1.csv")

# Info PCA rapida (carichi)
loadings <- as_tibble(
  pca_fit$rotation[, 1, drop = FALSE],
  rownames = "feature"
) %>%
  rename(loading_PC1 = PC1) %>%
  mutate(loading_PC1 = sgn * loading_PC1)
print(loadings)


# ===== 2) Tre fasi con C1/C2 e interazione con PID =====
# ===== Passo 2: Tre fasi con C1/C2 e interazione con PID =====
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(broom.mixed)
  library(performance)
})

# Usare Kenward–Roger (se preferisci Satterthwaite, scommenta la riga sotto)
emm_options(lmer.df = "kenward-roger")
# emm_options(lmer.df = "satterthwaite")

# --- 2.1) Lettura dati e definizione dei contrasti ortogonali ---
dat3 <- readr::read_csv(
  "acoustic_pid5_long_with_pc1.csv",
  show_col_types = FALSE
) %>%
  filter(phase %in% c("neutral", "pre", "post")) %>%
  mutate(
    # C1: recovery = post - pre (neutral=0; pre=-.5; post=+.5)
    C1 = case_when(
      phase == "pre" ~ -0.5,
      phase == "post" ~ 0.5,
      TRUE ~ 0
    ),
    # C2: attivazione = exam - neutral (neutral=-1; pre=+.5; post=+.5)
    C2 = case_when(
      phase == "neutral" ~ -1,
      phase %in% c("pre", "post") ~ 0.5,
      TRUE ~ NA_real_
    )
  )

# Outcome composito e dominio PID-5 primario (puoi cambiare qui)
yvar <- "mfcc_pc1"
pidvar <- "domain_disinhibition_z"

stopifnot(all(c(yvar, pidvar, "user_id", "C1", "C2") %in% names(dat3)))

# Dataset per il modello (y standardizzata per interpretazione in SD)
d3 <- dat3 %>%
  transmute(
    user_id,
    C1,
    C2,
    y = as.numeric(scale(.data[[yvar]])),
    pid_z = .data[[pidvar]]
  ) %>%
  drop_na()

# --- 2.2) Modello mixed con intercept casuale per soggetto ---
form3 <- y ~ C1 * pid_z + C2 * pid_z + (1 | user_id)
m3 <- lmer(form3, data = d3, REML = TRUE)

cat("\n=== Summary (lmerTest) ===\n")
print(summary(m3))

cat("\n=== ANOVA (KR/Satt) ===\n")
print(anova(m3)) # df via emm_options(lmer.df=...)

# CI 95% (Wald) per le interazioni principali
coef_names <- names(lme4::fixef(m3))
parm_ok <- intersect(
  c("C1:pid_z", "pid_z:C1", "C2:pid_z", "pid_z:C2"),
  coef_names
)
cat("\n=== CI 95% (Wald) per le interazioni disponibili ===\n")
print(confint(m3, parm = parm_ok, method = "Wald"))

# --- 2.3) Interpretazione con emmeans/emtrends ---

# (A) Δ(post–pre) come funzione di PID: slope su C1
#     (C1 è codificato in modo che lo slope = post - pre)
emm_C1 <- emtrends(
  m3,
  specs = ~ 1 | pid_z,
  var = "C1",
  at = list(pid_z = seq(-2, 2, by = 0.25))
)
cat("\n=== Δ(post–pre) = slope(C1) lungo PID (−2..+2) ===\n")
print(summary(emm_C1))

# Valori a PID = −1 e +1 SD
emm_C1_pts <- emtrends(
  m3,
  specs = ~ 1 | pid_z,
  var = "C1",
  at = list(pid_z = c(-1, 1))
)
cat("\n=== Δ(post–pre) a PID = −1 e +1 SD ===\n")
print(emm_C1_pts)

# Differenza di Δ tra +1 e −1 SD (coincide con l'interazione C1:pid_z)
# NB: serve by = NULL per confrontare tra livelli di 'pid_z'
cat("\n=== Confronto Δ@+1 − Δ@−1 (atteso = coefficiente C1:pid_z) ===\n")
print(contrast(emm_C1_pts, method = list("Δ@+1 − Δ@−1" = c(-1, 1)), by = NULL))

# (B) Attivazione exam–neutral come funzione di PID: slope su C2
emm_C2 <- emtrends(
  m3,
  specs = ~ 1 | pid_z,
  var = "C2",
  at = list(pid_z = seq(-2, 2, by = 0.25))
)
cat("\n=== Exam − Neutral = slope(C2) lungo PID (−2..+2) ===\n")
print(summary(emm_C2))

# --- 2.4) (Opzionale) Ricostruzione PRE/POST e contrasto post − pre ---
# PRE  = C1 = −0.5, C2 = +0.5;  POST = C1 = +0.5, C2 = +0.5
emm_PP <- emmeans(
  m3,
  ~ C1 | pid_z,
  at = list(
    C1 = c(-0.5, 0.5), # pre, post
    C2 = 0.5, # esame (media di pre e post)
    pid_z = c(-1, 1)
  )
)
cat("\n=== PRE vs POST a PID = −1 e +1 (emmeans sui valori di C1/C2) ===\n")
print(emm_PP)

cat("\n=== Contrasto (POST − PRE) a PID = −1 e +1 ===\n")
print(contrast(emm_PP, method = list("post − pre" = c(-1, 1)), by = "pid_z"))

# --- 2.5) Sanity check analitico: Δ(post–pre)(pid) = b_C1 + b_(C1:pid_z)*pid ---
b <- fixef(m3)
delta_at <- function(pid) {
  b["C1"] +
    (if ("C1:pid_z" %in% names(b)) b["C1:pid_z"] else b["pid_z:C1"]) * pid
}
cat("\n=== Δ(post–pre) atteso dal modello (analitico) a PID = −1 e +1 ===\n")
print(
  tibble(pid_z = c(-1, 1), Delta_post_pre = c(delta_at(-1), delta_at(1)))
)

# --- 2.6) (Opzionale) Tabella compatta dei fissi + R2 ---
cat("\n=== Fissi (righe chiave) ===\n")
print(
  broom.mixed::tidy(m3, effects = "fixed") %>%
    filter(
      term %in% c("(Intercept)", "C1", "pid_z", "C2", "C1:pid_z", "pid_z:C2")
    )
)

cat("\n=== R2 Nakagawa (marginale / condizionato) ===\n")
print(performance::r2(m3))

# --- 2.7) (Opzionale) dati per plot in ggplot ---
# Δ(post–pre) vs PID con CI
emmC1_df <- as_tibble(summary(emm_C1)) %>%
  rename(pid_z = pid_z, delta_post_pre = C1.trend, SE = SE) %>%
  mutate(lwr = lower.CL, upr = upper.CL)
# Esempio di plot:
# ggplot(emmC1_df, aes(pid_z, delta_post_pre)) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .2) +
#   geom_line() + geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "PID-5 (z)", y = "Δ (post − pre) su MFCC_PC1 (SD)")

# ===== 3) LOSO sui coefficienti C1:pid_z e C2:pid_z =====
loso_res <- lapply(unique(d3$user_id), function(id_drop) {
  m <- lmer(
    y ~ C1 * pid_z + C2 * pid_z + (1 | user_id),
    data = d3 %>% filter(user_id != id_drop),
    REML = TRUE
  )
  cf <- broom.mixed::tidy(m, effects = "fixed") %>%
    filter(term %in% c("C1:pid_z", "C2:pid_z")) %>%
    select(term, estimate, std.error, statistic) %>%
    mutate(left_out = id_drop)
  cf
}) %>%
  bind_rows()

cat("\n=== LOSO summary ===\n")
print(
  loso_res %>%
    group_by(term) %>%
    summarise(
      n = n(),
      mean = mean(estimate),
      sd = sd(estimate),
      min = min(estimate),
      max = max(estimate),
      sign_consistency = mean(sign(estimate) == sign(mean(estimate))),
      .groups = "drop"
    )
)


# ===== 4) Cluster bootstrap (per soggetto) =====
set.seed(123)
B <- 200 # aumenta se vuoi, ma con cautela per i tempi
ids <- unique(d3$user_id)

boot_fun <- function() {
  samp_ids <- sample(ids, size = length(ids), replace = TRUE)
  d_boot <- map_dfr(
    samp_ids,
    ~ d3 %>% filter(user_id == .x) %>% mutate(user_id = paste0(.x, "_b"))
  )
  m <- lmer(
    y ~ C1 * pid_z + C2 * pid_z + (1 | user_id),
    data = d_boot,
    REML = TRUE
  )
  broom.mixed::tidy(m, effects = "fixed") %>%
    filter(term %in% c("C1:pid_z", "C2:pid_z")) %>%
    select(term, estimate)
}

boot_res <- replicate(B, boot_fun(), simplify = FALSE) %>% bind_rows()
boot_ci <- boot_res %>%
  group_by(term) %>%
  summarise(
    mean = mean(estimate),
    l95 = quantile(estimate, 0.025),
    u95 = quantile(estimate, 0.975),
    .groups = "drop"
  )
cat("\n=== Cluster bootstrap 95% CI ===\n")
print(boot_ci)


# ===== 5a) Covariate (se presenti) =====
covs <- intersect(c("esi_bf", "sex", "age"), names(dat3))
message("Covariate trovate: ", paste(covs, collapse = ", "))

form_cov <- as.formula(
  paste(
    "y ~ C1*pid_z + C2*pid_z",
    if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else "",
    "+ (1 | user_id)"
  )
)

m3_cov <- lmer(
  form_cov,
  data = dat3 %>%
    transmute(
      user_id,
      C1,
      C2,
      y = as.numeric(scale(.data[[yvar]])),
      pid_z = .data[[pidvar]],
      across(all_of(covs), identity)
    ) %>%
    drop_na(),
  REML = TRUE
)
summary(m3_cov)

# ===== 5b) Power (simr) per l'interazione C1:pid_z =====
# install.packages("simr")
library(simr)

conflict_prefer("fixed", "simr")

# ===== Power analysis con simr per C1:pid_z =====

set.seed(123)

# 1) Modello base (già stimato prima come m3, ma lo rifaccio per chiarezza)
m_base <- lmer(
  y ~ C1 * pid_z + C2 * pid_z + (1 | user_id),
  data = d3,
  REML = TRUE
)

# 2) Imposta esplicitamente l'effetto da testare (puoi cambiare 'target' se vuoi pianificare un effetto minimo rilevante)
target <- fixef(m_base)["C1:pid_z"] # usa la stima osservata
fixef(m_base)["C1:pid_z"] <- as.numeric(target)

# 3) Potenza a N corrente (36 soggetti; 108 righe = 3 fasi x 36)
ps <- powerSim(m_base, nsim = 200, test = fixed("C1:pid_z"))
ps

# 4) Power curve aumentando i soggetti (mantiene 3 osservazioni per soggetto)
#    *** QUI STA IL FIX: specifica along = "user_id" ***
breaks_N <- seq(36, 80, by = 8) # prova anche fino a 100 se vuoi
pc <- powerCurve(
  m_base,
  along = "user_id",
  breaks = breaks_N,
  nsim = 100,
  test = fixed("C1:pid_z")
)

pc # tabella con N e potenza stimata
plot(pc) # curva di potenza

# 5) (Opzionale) tabella in tibble per report
pc_tbl <- as.data.frame(pc) %>%
  rename(N_subjects = x, power = y, power_lwr = ymin, power_upr = ymax)
pc_tbl
