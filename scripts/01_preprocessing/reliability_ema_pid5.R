# Carica le librerie necessarie
library(lavaan)
library(psych)
library(MASS)
library(dplyr)
library(lavaan)
library(semTools)

#  R function multilevel_alpha(), a wrapper for computing
#  the reliability indices discussed in
#  Lai, M. H. C. (2021). Composite reliability of multilevel data:
#      It’s about observed scores and construct meanings.
#      Psychological Methods, 26(1), 90–102.
#      https://doi.org/10.1037/met0000287
#  Copyright (C) 2021 Lai, Hok Chio (Mark)
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

multilevel_alpha <- function(
  data,
  id,
  nsim = 5000,
  conf_level = .95,
  se = "robust.huber.white",
  pa_args = list(fa = "pc"),
  ...
) {
  if (!require(lavaan)) stop("The lavaan package needs to be installed.")
  if (!require(psych)) stop("The psych package needs to be installed.")
  nitem <- ncol(data)
  ynames <- paste0("y", seq_len(nitem))
  colnames(data) <- ynames
  data <- cbind(data, id = id)
  tab_id <- table(id)
  hmean_cs <- 1 / mean(1 / tab_id[tab_id > 0])
  # Alpha
  # Generate syntax for saturated model
  sat_syntax <- (function(y) {
    if (length(y) <= 1) {
      return(NULL)
    }
    paste(
      c(paste(y[1], "~~", paste(y[-1], collapse = " + ")), Recall(y[-1])),
      collapse = "\n  "
    )
  })(ynames)
  msat <- paste0("level: 1\n  ", sat_syntax, "\nlevel: 2\n  ", sat_syntax)
  msat_fit <- cfa(
    msat,
    data = data,
    cluster = "id",
    se = se,
    test = "none",
    h1 = FALSE,
    baseline = FALSE,
    ...
  )
  coef_msat <- coef(msat_fit, type = "user")
  vcov_msat <- vcov(msat_fit)
  vw <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs == rhs & level == 1))
  ]
  cvw <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs != rhs & level == 1))
  ]
  vb <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs == rhs & level == 2))
  ]
  cvb <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs != rhs & level == 2))
  ]
  Sw <- sum(coef_msat[vw], 2 * coef_msat[cvw])
  Sb <- sum(coef_msat[vb], 2 * coef_msat[cvb])
  alpha_const <- nitem / (nitem - 1)
  alphaw <- alpha_const * sum(2 * coef_msat[cvw]) / Sw
  alpha2l <- alpha_const * sum(2 * coef_msat[c(cvw, cvb)]) / (Sb + Sw)
  alphab <- alpha_const * sum(2 * coef_msat[cvb]) / (Sb + Sw / hmean_cs)
  sim_coef_msat <- MASS::mvrnorm(
    nsim,
    mu = coef_msat[c(vw, vb, cvw, cvb)],
    Sigma = vcov_msat[c(vw, vb, cvw, cvb), c(vw, vb, cvw, cvb)]
  )
  sim_Sw <- rowSums(cbind(sim_coef_msat[, vw], 2 * sim_coef_msat[, cvw]))
  sim_Sb <- rowSums(cbind(sim_coef_msat[, vb], 2 * sim_coef_msat[, cvb]))
  sim_alphaw <- alpha_const * rowSums(2 * sim_coef_msat[, cvw]) / sim_Sw
  sim_alpha2l <- alpha_const *
    rowSums(2 * sim_coef_msat[, c(cvw, cvb)]) /
    (sim_Sb + sim_Sw)
  sim_alphab <- alpha_const *
    rowSums(2 * sim_coef_msat[, cvb]) /
    (sim_Sb + sim_Sw / hmean_cs)
  sim_alpha_cis <- lapply(
    list(alpha2l = sim_alpha2l, alphab = sim_alphab, alphaw = sim_alphaw),
    quantile,
    probs = .5 + c(-conf_level, conf_level) / 2
  )
  # Omega
  loading_labels <- paste0("l", seq_len(nitem))
  g_syntax <- paste(loading_labels, "*", ynames, collapse = " + ")
  mcfa <- paste0(
    "level: 1\n  fw =~ ",
    g_syntax,
    "\nlevel: 2\n  fb =~ ",
    g_syntax
  )
  mcfa_fit <- cfa(
    mcfa,
    data = data,
    cluster = "id",
    se = se,
    test = "none",
    h1 = TRUE,
    baseline = FALSE,
    ...
  )
  mcfa_pt <- partable(mcfa_fit)
  coef_mcfa <- coef(mcfa_fit)
  vcov_mcfa <- vcov(mcfa_fit)
  coef_mcfa[loading_labels]
  ld <- names(coef_mcfa)[with(
    mcfa_pt,
    free[which(
      op == "=~" &
        level == 1
    )]
  )]
  evw <- names(coef_mcfa)[with(
    mcfa_pt,
    free[which(
      op == "~~" &
        lhs == rhs &
        lhs != "fw" &
        level == 1
    )]
  )]
  fvw <- names(coef_mcfa)[with(
    mcfa_pt,
    free[which(
      op == "~~" &
        lhs == "fw"
    )]
  )]
  evb <- names(coef_mcfa)[with(
    mcfa_pt,
    free[which(
      op == "~~" &
        lhs == rhs &
        lhs != "fb" &
        level == 2
    )]
  )]
  fvb <- names(coef_mcfa)[with(
    mcfa_pt,
    free[which(
      op == "~~" &
        lhs == "fb"
    )]
  )]
  sumldsq <- sum(1, coef_mcfa[ld])^2
  sumevw <- sum(coef_mcfa[evw])
  sumevb <- sum(coef_mcfa[evb])
  omegaw <- sumldsq * coef_mcfa[[fvw]] / (sumldsq * coef_mcfa[[fvw]] + sumevw)
  omega2l <- sum(sumldsq * coef_mcfa[c(fvw, fvb)]) /
    sum(sumldsq * coef_mcfa[c(fvw, fvb)], sumevw, sumevb)
  omegab <- sumldsq *
    coef_mcfa[[fvb]] /
    (sumldsq *
      (coef_mcfa[[fvb]] + coef_mcfa[[fvw]] / hmean_cs) +
      sumevb +
      sumevw / hmean_cs)
  sim_coef_mcfa <- MASS::mvrnorm(
    nsim,
    mu = coef_mcfa[c(ld, fvw, fvb, evw, evb)],
    Sigma = vcov_mcfa[c(ld, fvw, fvb, evw, evb), c(ld, fvw, fvb, evw, evb)]
  )
  sim_sumldsq <- (1 + rowSums(sim_coef_mcfa[, ld]))^2
  sim_sumevw <- rowSums(sim_coef_mcfa[, evw])
  sim_sumevb <- rowSums(sim_coef_mcfa[, evb])
  sim_omegaw <- sim_sumldsq *
    sim_coef_mcfa[, fvw] /
    (sim_sumldsq * sim_coef_mcfa[, fvw] + sim_sumevw)
  sim_omega2l <- rowSums(sim_sumldsq * sim_coef_mcfa[, c(fvw, fvb)]) /
    (rowSums(sim_sumldsq * sim_coef_mcfa[, c(fvw, fvb)]) +
      sim_sumevw +
      sim_sumevb)
  sim_omegab <- sim_sumldsq *
    sim_coef_mcfa[, fvb] /
    (sim_sumldsq *
      (sim_coef_mcfa[, fvb] +
        sim_coef_mcfa[, fvw] / hmean_cs) +
      sim_sumevb +
      sim_sumevw / hmean_cs)
  sim_omega_cis <- lapply(
    list(omega2l = sim_omega2l, omegab = sim_omegab, omegaw = sim_omegaw),
    quantile,
    probs = .5 + c(-conf_level, conf_level) / 2
  )
  # resid_corb <- resid(mcfa_fit, type = "cor")$id$cov
  # diag(resid_corb) <- 1
  # psych::KMO(resid_corb)$MSA
  # psych::fa.parallel(resid_corb, fm = "pa", fa = "fa",
  #                    n.obs = lavTech(msat_fit, "nclusters")[[1]],
  #                    n.iter = 30 * nitem,
  #                    plot = FALSE)$nfact
  # Dimensionality
  corw <- lavTech(msat_fit, what = "cor.ov")$within
  corb <- lavTech(msat_fit, what = "cor.ov")$id
  paw <- do.call(
    fa.parallel,
    args = c(
      list(
        x = corw,
        n.obs = nobs(msat_fit) -
          lavTech(msat_fit, "nclusters")[[1]],
        n.iter = 30 * nitem,
        plot = FALSE
      ),
      pa_args
    )
  )
  pab <- do.call(
    fa.parallel,
    args = c(
      list(
        x = corb,
        n.obs = lavTech(msat_fit, "nclusters")[[1]],
        n.iter = 30 * nitem,
        plot = FALSE
      ),
      pa_args
    )
  )
  if (pa_args$fa == "pc") {
    ncompw <- paw$ncomp
    ncompb <- pab$ncomp
  } else if (pa_args$fa == "fa") {
    ncompw <- paw$nfact
    ncompb <- pab$nfact
  }
  list(
    alpha = c(alpha2l = alpha2l, alphab = alphab, alphaw = alphaw),
    alpha_ci = do.call(rbind, sim_alpha_cis),
    omega = c(omega2l = omega2l, omegab = omegab, omegaw = omegaw),
    omega_ci = do.call(rbind, sim_omega_cis),
    ncomp = c(within = ncompw, between = ncompb)
  )
}

# ==============================================================================
# 1. LOAD EMA DATA
# ==============================================================================

data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
d <- rio::import(data_path)


# PID5 EMA
# df$pid5_sum <- rowSums(df[, 5:19])
# df$pid5_negative_affectivity  # pid5_1, pid5_2, pid5_3
# df$pid5_detachment # pid5_4, pid5_5, pid5_6
# df$pid5_antagonism  # pid5_7, pid5_8, pid5_9
# df$pid5_disinhibition # pid5_10, pid5_11, pid5_12
# df$pid5_psychoticism  # pid5_13, pid5_14, pid5_15
# 15 item, 3 per dominio
ema_pid5_items <- c(
  "pid5_2",
  "pid5_11",
  "pid5_13",
  "pid5_15",
  "pid5_3",
  "pid5_7",
  "pid5_14",
  "pid5_6",
  "pid5_4",
  "pid5_12",
  "pid5_1",
  "pid5_9",
  "pid5_5",
  "pid5_8",
  "pid5_10"
)

d1 <- d %>%
  dplyr::select(
    user_id,
    day,
    all_of(ema_pid5_items),
    # all_of(baseline_pid5_vars)
  ) %>%
  dplyr::rename(ID = user_id)

glimpse(d1)


# Distribuzione osservazioni per soggetto
obs_per_subject <- d1 %>%
  group_by(ID) %>%
  summarise(
    n_obs = n(),
    .groups = "drop"
  )

cat("=== DISTRIBUZIONE OSSERVAZIONI PER SOGGETTO ===\n")
summary(obs_per_subject$n_obs)
cat("\nSoggetti con < 5 osservazioni:", sum(obs_per_subject$n_obs < 5), "\n")
cat("Soggetti con < 10 osservazioni:", sum(obs_per_subject$n_obs < 10), "\n")

# Visualizza distribuzione
hist(
  obs_per_subject$n_obs,
  main = "Distribuzione osservazioni per soggetto",
  xlab = "N osservazioni",
  ylab = "Frequenza",
  breaks = 30
)
abline(v = 5, col = "red", lwd = 2, lty = 2)

# Identifica soggetti con poche osservazioni
low_obs_subjects <- obs_per_subject %>%
  filter(n_obs < 10 | n_obs > 35) %>%
  pull(ID)

if (length(low_obs_subjects) > 0) {
  cat("\nSoggetti da considerare per esclusione:\n")
  print(obs_per_subject %>% filter(ID %in% low_obs_subjects))
}

bad_id <- c(
  "al_be_2002_09_10_525",
  "al_pe_2004_02_28_769",
  "el_to_2005_09_02_232",
  "la_wh_2004_12_13_861",
  "ma_ca_2001_03_13_546"
)
d2 <- d1[!(d1$ID %in% bad_id), ]


# Definisci la struttura degli item per dimensione
dim_structure <- list(
  Dim1 = c("pid5_2", "pid5_11", "pid5_13"),
  Dim2 = c("pid5_15", "pid5_3", "pid5_7"),
  Dim3 = c("pid5_14", "pid5_6", "pid5_4"),
  Dim4 = c("pid5_12", "pid5_1", "pid5_9"),
  Dim5 = c("pid5_5", "pid5_8", "pid5_10")
)

# Prepara i dati
all_items <- unlist(dim_structure)
model_data <- d2[, c("ID", all_items)]
colnames(model_data) <- c("ID", paste0("y", 1:15))
# Crea etichette per item e loading
item_labels <- paste0("y", 1:15)
loading_labels <- paste0("l", 1:15)

model_data2 <- model_data

# Se sono factor/ordered: converti correttamente preservando i valori 0-4
model_data2[items] <- lapply(model_data2[items], function(x) {
  if (is.factor(x) || is.ordered(x)) {
    as.numeric(as.character(x))
  } else {
    as.numeric(x)
  }
})

# Controllo veloce: devono essere numeric e NON ordered/factor
str(model_data2[items])

items <- paste0("y", 1:15)

# conta quanti item hanno varianza within = 0 per ogni soggetto
zero_within_by_id <- aggregate(
  model_data2[items],
  by = list(ID = model_data2$ID),
  FUN = function(x) var(x, na.rm = TRUE)
)

# var==0 -> item costante
zero_within_by_id$zero_items <- rowSums(
  zero_within_by_id[items] == 0,
  na.rm = TRUE
)

summary(zero_within_by_id$zero_items)

# esempio: escludi solo chi ha >= 10 item costanti (decidi tu la soglia)
bad_ids <- zero_within_by_id$ID[zero_within_by_id$zero_items >= 10]

model_data3 <- subset(model_data2, !(ID %in% bad_ids))

length(unique(model_data2$ID))
length(unique(model_data3$ID))

# ==============================================================================
# 2. Compute reliability
# ==============================================================================

mod_syntax <- '
  level: 1
    # Unico fattore within
    f_w =~ pid_1 + pid_2 + pid_3 + pid_4 + pid_5 + pid_6 + pid_7 + 
           pid_8 + pid_9 + pid_10 + pid_11 + pid_12 + pid_13 + 
           pid_14 + pid_15
  
  level: 2
    # Unico fattore between
    f_b =~ pid_1 + pid_2 + pid_3 + pid_4 + pid_5 + pid_6 + pid_7 + 
           pid_8 + pid_9 + pid_10 + pid_11 + pid_12 + pid_13 + 
           pid_14 + pid_15
'

items <- paste0("y", 1:15)

rel <- multilevel_alpha(
  data = model_data3[items], # item numerici
  id = model_data3$ID,
  nsim = 5000,
  conf_level = .95,
  se = "robust.huber.white"
)

# Risultati chiave
rel$alpha
rel$omega
rel$alpha_ci
rel$omega_ci

# > # Risultati chiave
# > rel$alpha
# alpha2l    alphab    alphaw
# 0.8256431 0.8663292 0.7273300
# > rel$omega
# omega2l    omegab    omegaw
# 0.8138011 0.8521924 0.7237887
# > rel$alpha_ci
# 2.5%     97.5%
# alpha2l 0.8018422 0.8450367
# alphab  0.8353568 0.8889806
# alphaw  0.6867836 0.7576999

#' Multilevel composite reliability estimates indicated adequate
#' within-person reliability and high between-person reliability.
#' Following Lai (2021), reliability was therefore evaluated at
#' the level of observed scores rather than latent factor fit.
