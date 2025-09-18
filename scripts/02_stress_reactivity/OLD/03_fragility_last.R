# ============================================================
# Analisi "fragilità psicologica" — versione STEP-BY-STEP
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(lubridate)
  library(rio)
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(bayesplot)
  library(conflicted)
})

# ------------------------------------------------------------
# 0) Opzioni, preferenze di conflitto, utilità
# ------------------------------------------------------------
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("sd", "stats")
conflict_prefer("var", "stats")
conflict_prefer("mad", "stats")
conflict_prefer("match", "base")
conflict_prefer("lag", "dplyr")

# Z-score "sicuro": evita divisione per 0 restituendo dati centrati
z_ <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - m)
  (x - m) / s
}

prepare_fragility_data_enhanced <- function(d, min_obs_per_period = 1) {
  # =============================
  # Utility interne
  # =============================
  # Item 0–100 / altro -> 1..7 (per ordered_logistic)
  as_item_to_1_7 <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      x <- trimws(x)
      x[x == ""] <- NA_character_
      suppressWarnings(
        x <- readr::parse_number(
          x,
          locale = readr::locale(decimal_mark = ".", grouping_mark = ",")
        )
      )
    }
    x <- as.numeric(x)
    if (all(is.na(x))) return(as.integer(x))

    # già 1..7?
    if (all(x[is.finite(x)] %in% 1:7)) return(as.integer(x))

    xmin <- suppressWarnings(min(x, na.rm = TRUE))
    xmax <- suppressWarnings(max(x, na.rm = TRUE))

    # chiaramente 0..100 -> binning uniforme in 7 classi
    if (is.finite(xmin) && is.finite(xmax) && xmin >= 0 && xmax <= 100) {
      brk <- seq(0, 100, length.out = 8) # 8 breakpoints => 7 classi
      y <- findInterval(x, brk, rightmost.closed = TRUE, all.inside = TRUE)
      return(as.integer(y)) # 1..7
    }

    # fallback: riscalamento lineare del range osservato a 1..7
    y <- 1 + 6 * (x - xmin) / (xmax - xmin)
    y <- pmax(1, pmin(7, y))
    as.integer(round(y))
  }

  z_ <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(x - m)
    (x - m) / s
  }

  find_first <- function(df, candidates) {
    nm <- names(df)
    hit <- nm[tolower(nm) %in% tolower(candidates)]
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  # =============================
  # 1) ID soggetto (-> 1..I)
  # =============================
  subj_col <- d %>%
    dplyr::select(dplyr::any_of(c(
      "user_id",
      "id",
      "subject_id",
      "participant_id"
    ))) %>%
    names() %>%
    .[1]
  if (is.na(subj_col))
    stop("Nessuna colonna soggetto trovata (es. 'user_id'/'id').")
  d$.__subj__ <- as.integer(factor(d[[subj_col]]))

  # =============================
  # 2) Periodo (1=baseline, 2=pre, 3=post)
  # =============================
  per_col <- d %>%
    dplyr::select(dplyr::any_of(c("exam_period", "period", "phase"))) %>%
    names() %>%
    .[1]
  if (!is.na(per_col)) {
    per_raw <- d[[per_col]]
    if (is.numeric(per_raw)) {
      per <- as.integer(per_raw)
    } else {
      lv <- tolower(as.character(per_raw))
      per <- dplyr::case_when(
        lv %in% c("baseline", "base", "t0", "bl") ~ 1L,
        lv %in% c("pre", "pre_exam", "pre-exam", "preexam") ~ 2L,
        lv %in% c("post", "post_exam", "post-exam", "postexam") ~ 3L,
        TRUE ~ 1L
      )
    }
  } else {
    per <- rep(1L, nrow(d))
  }
  d$.__per__ <- as.integer(per)

  # =============================
  # 3) Sesso -> female (1=femmina, 0=maschio) a livello soggetto
  # =============================
  sex_col <- d %>%
    dplyr::select(dplyr::any_of(c("female", "sex", "gender", "sesso"))) %>%
    names() %>%
    .[1]
  if (!is.na(sex_col)) {
    sx <- d[[sex_col]]
    if (is.numeric(sx)) {
      female_by_row <- as.integer(sx != 0)
    } else {
      lv <- tolower(as.character(sx))
      female_by_row <- as.integer(
        lv %in% c("f", "female", "donna", "femmina", "woman", "1", "true")
      )
    }
  } else {
    warning(
      "Colonna sesso assente: imposto female=1 per tutti come placeholder."
    )
    female_by_row <- rep(1L, nrow(d))
  }
  # media per soggetto (≥.5 -> femmina)
  fem_mean <- tapply(
    female_by_row,
    d$.__subj__,
    function(v) mean(v, na.rm = TRUE)
  )
  fem_mean[is.na(fem_mean)] <- 1
  idx_subj <- as.integer(names(fem_mean))
  vals_bin <- as.integer(fem_mean >= 0.5)

  # =============================
  # 4) Items (happy, sad, satisfied, angry) -> 1..7
  # =============================
  happy_col <- find_first(d, c("happy", "happiness", "felice", "contento"))
  sad_col <- find_first(d, c("sad", "sadness", "triste"))
  satisfied_col <- find_first(
    d,
    c("satisfied", "satisfaction", "soddisfatto", "soddisfazione")
  )
  angry_col <- find_first(d, c("angry", "anger", "arrabbiato", "rabbia"))
  if (
    any(sapply(list(happy_col, sad_col, satisfied_col, angry_col), is.null))
  ) {
    stop("Non trovo tutte e 4 le colonne item (happy, sad, satisfied, angry).")
  }
  d$.__happy__ <- as_item_to_1_7(d[[happy_col]])
  d$.__sad__ <- as_item_to_1_7(d[[sad_col]])
  d$.__satisfied__ <- as_item_to_1_7(d[[satisfied_col]])
  d$.__angry__ <- as_item_to_1_7(d[[angry_col]])

  # =============================
  # 5) EMA per-beep (dimensioni 2..5) — **nuovi nomi**
  #    dim=2 detachment, 3 antagonism, 4 disinhibition, 5 psychoticism
  # =============================
  ema_cols <- list(
    det = find_first(
      d,
      c("pid5_detachment", "ema_detachment", "detachment", "z_det", "det")
    ),
    ant = find_first(
      d,
      c("pid5_antagonism", "ema_antagonism", "antagonism", "z_ant", "ant")
    ),
    dis = find_first(
      d,
      c(
        "pid5_disinhibition",
        "ema_disinhibition",
        "disinhibition",
        "z_dis",
        "dis"
      )
    ),
    psy = find_first(
      d,
      c("pid5_psychoticism", "ema_psychoticism", "psychoticism", "z_psy", "psy")
    )
  )
  has_ema <- !all(sapply(ema_cols, is.null))

  # =============================
  # 6) Baseline PID-5 (5 dimensioni) a livello soggetto
  #    (si assume che tu abbia già rinominato *_baseline a monte)
  # =============================
  base_cols <- c(
    "pid5_negative_affect_baseline",
    "pid5_detachment_baseline",
    "pid5_antagonism_baseline",
    "pid5_disinhibition_baseline",
    "pid5_psychoticism_baseline"
  )
  if (!all(base_cols %in% names(d))) {
    stop(
      "Mancano baseline PID-5 attese: ",
      paste(setdiff(base_cols, names(d)), collapse = ", ")
    )
  }

  # =============================
  # 7) Tieni solo le osservazioni con i 4 item presenti
  # =============================
  keep_obs <- with(
    d,
    !is.na(.__happy__) &
      !is.na(.__sad__) &
      !is.na(.__satisfied__) &
      !is.na(.__angry__)
  )
  d2 <- d[keep_obs, , drop = FALSE]
  if (nrow(d2) == 0)
    stop(
      "Nessuna osservazione con i 4 item completi (dopo conversione 0–100→1–7)."
    )

  # (opzionale) almeno min_obs_per_period per soggetto in (almeno) un periodo
  if (min_obs_per_period > 0) {
    cnt <- d2 %>%
      dplyr::count(.__subj__, .__per__) %>%
      tidyr::pivot_wider(
        names_from = .__per__,
        values_from = n,
        values_fill = 0
      )
    for (lvl in c("1", "2", "3")) if (!lvl %in% names(cnt)) cnt[[lvl]] <- 0L
    valid_subj <- cnt$`2` >= min_obs_per_period |
      cnt$`3` >= min_obs_per_period |
      cnt$`1` >= min_obs_per_period
    keep_ids <- cnt$`.__subj__`[valid_subj]
    d2 <- d2[d2$.__subj__ %in% keep_ids, , drop = FALSE]
  }

  # =============================
  # 8) Indici finali
  # =============================
  N_obs <- nrow(d2)
  obs_id <- seq_len(N_obs)
  subject <- as.integer(d2$.__subj__)
  period <- as.integer(d2$.__per__)
  I <- max(subject)

  # =============================
  # 9) Long vectors per gli item (N_items=4*N_obs)
  # =============================
  y_wide <- cbind(d2$.__happy__, d2$.__sad__, d2$.__satisfied__, d2$.__angry__)
  colnames(y_wide) <- c("happy", "sad", "satisfied", "angry")
  N_items <- 4L * N_obs
  y_item <- as.integer(c(t(y_wide)))
  item_id <- as.integer(rep(1:4, times = N_obs))
  obs_id_long <- as.integer(rep(obs_id, each = 4L))

  # =============================
  # 10) Baseline X_base (I x 5), z-score per soggetto
  # =============================
  # =============================
  # 10) Baseline X_base (I x 5), z-score per soggetto
  #     - Se più righe per soggetto: prendi la prima non-NA
  #     - Z-score per colonna
  #     - NA residui (soggetti senza baseline) -> imputa a 0 (media)
  # =============================
  base_by_subj <- d %>%
    dplyr::select(.__subj__, dplyr::all_of(base_cols)) %>%
    dplyr::group_by(.__subj__) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(base_cols),
        ~ {
          tmp <- .
          idx <- which(!is.na(tmp))[1]
          if (is.na(idx)) NA_real_ else as.numeric(tmp[idx])
        }
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.__subj__)

  Xb <- as.matrix(base_by_subj[, base_cols, drop = FALSE])
  Xb <- apply(Xb, 2, z_)
  colnames(Xb) <- c("z_naff_b", "z_det_b", "z_ant_b", "z_dis_b", "z_psy_b")

  if (nrow(Xb) < I) {
    Xb_full <- matrix(NA_real_, nrow = I, ncol = 5)
    Xb_full[base_by_subj$.__subj__, ] <- Xb
    Xb <- apply(Xb_full, 2, z_)
    colnames(Xb) <- c("z_naff_b", "z_det_b", "z_ant_b", "z_dis_b", "z_psy_b")
  }

  # >>> PATCH: imputa a 0 gli NA residui (soggetti privi di baseline)
  n_na_before <- sum(!is.finite(Xb))
  if (n_na_before > 0) {
    Xb[!is.finite(Xb)] <- 0
    message(sprintf(
      "[prepare_fragility_data_enhanced] Imputati a 0 (post z-score) %d valori baseline mancanti.",
      n_na_before
    ))
  }

  # =============================
  # 11) female (I) con indici preservati
  # =============================
  female_vec <- integer(I)
  female_vec[idx_subj] <- vals_bin

  # =============================
  # 12) EMA per-beep (dim 2..5) — vettori sparsificati
  # =============================
  M_ema <- 0L
  ema_val <- numeric(0)
  ema_dim <- integer(0)
  ema_obs <- integer(0)

  if (has_ema) {
    EMA_df <- tibble::tibble(
      det = if (!is.null(ema_cols$det)) as.numeric(d2[[ema_cols$det]]) else
        NA_real_,
      ant = if (!is.null(ema_cols$ant)) as.numeric(d2[[ema_cols$ant]]) else
        NA_real_,
      dis = if (!is.null(ema_cols$dis)) as.numeric(d2[[ema_cols$dis]]) else
        NA_real_,
      psy = if (!is.null(ema_cols$psy)) as.numeric(d2[[ema_cols$psy]]) else
        NA_real_
    )
    # z-score su non-NA (per mettere su scala comparabile)
    EMA_df <- EMA_df %>%
      dplyr::mutate(
        det = z_(det),
        ant = z_(ant),
        dis = z_(dis),
        psy = z_(psy)
      )
    append_dim <- function(vals, dim_index) {
      ok <- which(is.finite(vals))
      if (length(ok) > 0) {
        ema_val <<- c(ema_val, vals[ok])
        ema_dim <<- c(ema_dim, rep(dim_index, length(ok)))
        ema_obs <<- c(ema_obs, obs_id[ok])
      }
    }
    append_dim(EMA_df$det, 2L)
    append_dim(EMA_df$ant, 3L)
    append_dim(EMA_df$dis, 4L)
    append_dim(EMA_df$psy, 5L)
    M_ema <- length(ema_val)
  }

  # =============================
  # 13) stan_data
  # =============================
  stan_data <- list(
    I = I,
    N_obs = N_obs,
    K = 4L,
    P = 3L,
    D = 5L,

    N_items = N_items,
    y_item = y_item,
    item_id = item_id,
    obs_id = obs_id_long,

    subject = subject,
    period = period,

    M_ema = M_ema,
    ema_val = if (M_ema > 0) as.vector(ema_val) else numeric(),
    ema_dim = if (M_ema > 0) as.integer(ema_dim) else integer(),
    ema_obs = if (M_ema > 0) as.integer(ema_obs) else integer(),

    X_base = unname(Xb),
    female = as.integer(female_vec),

    N_interact = as.integer(choose(5L, 2L)) # D*(D-1)/2
    # use_ema lo imposti tu (0/1) prima del fit
  )

  # =============================
  # 14) Meta-info
  # =============================
  meta <- list(
    item_names = c("happy", "sad", "satisfied", "angry"),
    ema_names = c("detachment", "antagonism", "disinhibition", "psychoticism"),
    baseline_names = colnames(Xb),
    female_table = table(female_vec),
    period_counts = table(factor(
      period,
      levels = 1:3,
      labels = c("baseline", "pre", "post")
    ))
  )

  list(stan_data = stan_data, meta = meta)
}

set.seed(20250917)

# ------------------------------------------------------------
# 1) Caricamento dati grezzi + rinomina colonne baseline
#    WHY: normalizziamo i nomi per rendere robusto lo step
#         di preparazione (dataset con header leggermente diversi)
# ------------------------------------------------------------
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
stan_file_path <- here::here(
  "scripts",
  "02_stress_reactivity",
  "fragility_v3_rhs.stan"
)

d <- rio::import(data_path) %>%
  dplyr::rename(
    pid5_negative_affect_baseline = any_of(c(
      "domain_negative_affect_baseline",
      "pid5_negative_affect_baseline"
    )),
    pid5_detachment_baseline = any_of(c(
      "domain_detachment_baseline",
      "pid5_detachment_baseline"
    )),
    pid5_antagonism_baseline = any_of(c(
      "domain_antagonism_baseline",
      "pid5_antagonism_baseline"
    )),
    pid5_disinhibition_baseline = any_of(c(
      "domain_disinhibition_baseline",
      "pid5_disinhibition_baseline"
    )),
    pid5_psychoticism_baseline = any_of(c(
      "domain_psychoticism_baseline",
      "pid5_psychoticism_baseline"
    ))
  )

# ------------------------------------------------------------
# 2) Preparazione dati per Stan
#    WHY: costruisce gli indici long (item, osservazione, soggetto, periodo),
#         z-score, matrici baseline e mappa EMA (dimensioni 2..5), ecc.
#         (Assumo che tu abbia già definito prepare_fragility_data_enhanced())
# ------------------------------------------------------------
frag_data_enh <- prepare_fragility_data_enhanced(d, min_obs_per_period = 1)
sd <- frag_data_enh$stan_data

cat("=== Preparazione Dati Completata ===\n")
cat(sprintf(
  "Soggetti: %d | Osservazioni: %d | Item: %d | EMA: %d\n",
  sd$I,
  sd$N_obs,
  sd$N_items,
  sd$M_ema
))
cat(sprintf("Interazioni baseline×EMA: %d\n", sd$N_interact))
cat(sprintf(
  "Genere - Femmine: %d | Maschi: %d\n",
  sum(sd$female),
  sum(1 - sd$female)
))
cat(sprintf(
  "Periodo: baseline=%d | pre=%d | post=%d\n",
  sum(sd$period == 1),
  sum(sd$period == 2),
  sum(sd$period == 3)
))

# ------------------------------------------------------------
# 3) Compilazione modello Stan
#    WHY: compiliamo una sola volta e riutilizziamo per A e B
# ------------------------------------------------------------
cat("\n2) Compilazione modello Stan…\n")
model_frag <- cmdstan_model(stan_file_path)

# ------------------------------------------------------------
# 4) Fit Modello A (baseline-only): use_ema = 0
#    WHY: modello base di riferimento (senza predittori EMA dinamici)
# ------------------------------------------------------------
cat("\n3) Fit Modello A (baseline-only)…\n")
data_A <- sd
data_A$use_ema <- 0.0

fit_A <- model_frag$variational(
  data = data_A,
  seed = 20250917,
  algorithm = "meanfield",
  elbo_samples = 100,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 0
)

# ------------------------------------------------------------
# 5) Fit Modello B (baseline + EMA): use_ema = 1
#    WHY: modello "esteso" che include theta e interazioni baseline×EMA
# ------------------------------------------------------------
cat("4) Fit Modello B (baseline + EMA)…\n")
data_B <- sd
data_B$use_ema <- 1.0

fit_B <- model_frag$variational(
  data = data_B,
  seed = 20250917,
  algorithm = "meanfield",
  elbo_samples = 100,
  adapt_engaged = TRUE,
  tol_rel_obj = 0.001,
  eval_elbo = 100,
  output_samples = 2000,
  refresh = 0
)

# ------------------------------------------------------------
# 6) Confronto R^2 bayesiano (in-sample)
#    WHY: misura "interna" della frazione di varianza spiegata del target
#         (meno importante di LOO, ma utile come indicatore descrittivo)
# ------------------------------------------------------------
cat("\n5) Confronto R² Bayesiano…\n")
R2_A <- as_draws_df(fit_A$draws("R2_frag"))$R2_frag
R2_B <- as_draws_df(fit_B$draws("R2_frag"))$R2_frag

cat(sprintf("R² baseline-only:  %.3f (SD=%.3f)\n", mean(R2_A), sd(R2_A)))
cat(sprintf("R² baseline+EMA:  %.3f (SD=%.3f)\n", mean(R2_B), sd(R2_B)))
delta_R2 <- R2_B - R2_A
cat(sprintf(
  "ΔR² (B-A):        %.3f [95%%: %.3f, %.3f];  P(Δ>0)=%.3f\n",
  mean(delta_R2),
  quantile(delta_R2, .025),
  quantile(delta_R2, .975),
  mean(delta_R2 > 0)
))

# ------------------------------------------------------------
# 7) Coefficienti chiave (Modello B)
#    WHY: interpretiamo effetti principali/di interazione di maggiore interesse
# ------------------------------------------------------------
cat("\n6) Coefficienti chiave (Modello B)…\n")
b_female <- as_draws_df(fit_B$draws("b_female"))$b_female
c_ema_draws <- as_draws_df(fit_B$draws("c_ema"))
b_var_ema <- as_draws_df(fit_B$draws("b_var_ema"))$b_var_ema
b_var_female <- as_draws_df(fit_B$draws("b_var_female"))$b_var_female

cat(sprintf(
  "Effetto genere:  %.3f [%.3f, %.3f] | P(>0)=%.3f\n",
  mean(b_female),
  quantile(b_female, .025),
  quantile(b_female, .975),
  mean(b_female > 0)
))

dim_names <- c(
  "NA_latent",
  "Detachment",
  "Antagonism",
  "Disinhibition",
  "Psychoticism"
)
cat("Coefficienti EMA:\n")
for (i in 1:5) {
  v <- c_ema_draws[[paste0("c_ema[", i, "]")]]
  cat(sprintf(
    "  %-13s: %.3f [%.3f, %.3f] | P(>0)=%.3f\n",
    dim_names[i],
    mean(v),
    quantile(v, .025),
    quantile(v, .975),
    mean(v > 0)
  ))
}
cat(sprintf(
  "Effetto var_EMA: %.3f [%.3f, %.3f] | P(>0)=%.3f\n",
  mean(b_var_ema),
  quantile(b_var_ema, .025),
  quantile(b_var_ema, .975),
  mean(b_var_ema > 0)
))
cat(sprintf(
  "Genere×Var_EMA: %.3f [%.3f, %.3f] | P(>0)=%.3f\n",
  mean(b_var_female),
  quantile(b_var_female, .025),
  quantile(b_var_female, .975),
  mean(b_var_female > 0)
))

# ------------------------------------------------------------
# 8) PSIS-LOO a livello osservazione
#    WHY: misura la qualità predittiva "out-of-sample" (prioritaria su R²)
#         moment_match migliora la stabilità con VI
# ------------------------------------------------------------
cat("\n7) LOO (osservazione)…\n")
ll_A <- fit_A$draws("log_lik_obs", format = "matrix") # draws x N_obs
ll_B <- fit_B$draws("log_lik_obs", format = "matrix")

loo_A <- loo::loo(ll_A, moment_match = TRUE)
loo_B <- loo::loo(ll_B, moment_match = TRUE)
cmp_obs <- loo::loo_compare(list(B = loo_B, A = loo_A))
print(cmp_obs)

delta_elpd_obs <- as.numeric(cmp_obs["A", "elpd_diff"]) * -1
se_diff_obs <- as.numeric(cmp_obs["A", "se_diff"])
cat(sprintf(
  "\nΔELPD (B-A): %.1f ± %.1f  | decisivo se |Δ|>2SE → %s\n",
  delta_elpd_obs,
  se_diff_obs,
  ifelse(abs(delta_elpd_obs) > 2 * se_diff_obs, "SÌ", "NO")
))

# ------------------------------------------------------------
# 9) PSIS-LOO a livello soggetto (aggregato)
#    WHY: sommiamo le log-lik delle osservazioni appartenenti allo stesso
#         soggetto per valutare la qualità predittiva a "grain" soggetto
#         (evitiamo %*% sparse; usiamo rowsum su transposte)
# ------------------------------------------------------------
cat("\n8) LOO (soggetto)…\n")
subj_idx <- sd$subject # vettore di lunghezza N_obs con l'ID soggetto per osservazione
I <- sd$I
# ll_* sono (draws x N_obs). Trasponiamo, sommiamo per soggetto, ritrasponiamo:
ll_A_subj <- t(rowsum(t(ll_A), group = subj_idx, reorder = FALSE)) # draws x I
ll_B_subj <- t(rowsum(t(ll_B), group = subj_idx, reorder = FALSE))

loo_A_subj <- loo::loo(ll_A_subj, moment_match = TRUE)
loo_B_subj <- loo::loo(ll_B_subj, moment_match = TRUE)
cmp_subj <- loo::loo_compare(list(B = loo_B_subj, A = loo_A_subj))
print(cmp_subj)

delta_elpd_subj <- as.numeric(cmp_subj["A", "elpd_diff"]) * -1
se_diff_subj <- as.numeric(cmp_subj["A", "se_diff"])
cat(sprintf(
  "\nΔELPD_sog (B-A): %.1f ± %.1f  | decisivo: %s\n",
  delta_elpd_subj,
  se_diff_subj,
  ifelse(abs(delta_elpd_subj) > 2 * se_diff_subj, "SÌ", "NO")
))

# ------------------------------------------------------------
# 10) Verdetto descrittivo
#     WHY: riassunto leggibile dei tre blocchi informativi (R², LOO obs, LOO subj)
# ------------------------------------------------------------
cat("\n=== VERDETTO FINALE ===\n")
winner_obs <- rownames(cmp_obs)[1]
winner_subj <- rownames(cmp_subj)[1]
cat(sprintf(
  "R²:         vince %s (Δ=%.3f, P(Δ>0)=%.2f)\n",
  ifelse(mean(delta_R2) > 0, "B", "A"),
  abs(mean(delta_R2)),
  ifelse(mean(delta_R2) > 0, mean(delta_R2 > 0), mean(delta_R2 <= 0))
))
cat(sprintf(
  "LOO obs:    vince %s (ΔELPD=%.1f)\n",
  winner_obs,
  abs(delta_elpd_obs)
))
cat(sprintf(
  "LOO sogg.:  vince %s (ΔELPD=%.1f)\n",
  winner_subj,
  abs(delta_elpd_subj)
))

# ------------------------------------------------------------
# 11) PPC sugli item ordinali
#     WHY: verifica calibrazione del blocco di misurazione ordinale
# ------------------------------------------------------------
cat("\n9) PPC — item ordinali…\n")
y_item_obs <- sd$y_item
y_item_rep_mat <- posterior::as_draws_matrix(
  fit_B$draws("y_item_rep") # richiede che in Stan tu abbia generato y_item_rep
)

# distribuzione per categoria
ppc_bars(y_item_obs, y_item_rep_mat[1:200, ])
# rootogram (discreti)
ppc_rootogram(y_item_obs, y_item_rep_mat[1:200, ])

# ------------------------------------------------------------
# 12) PPC sulla "fragilità" (diff pre–post) costruita dagli item (mele vs mele)
#     WHY: confrontiamo la stessa statistica sui dati osservati e sui replicati;
#          evitiamo C stack usando rowsum(), non matrici sparse
# ------------------------------------------------------------
cat("\n10) PPC — fragilità pre–post per soggetto…\n")

# 12.1) Dalla GQ: draws degli item replicati
Yrep <- posterior::as_draws_matrix(fit_B$draws("y_item_rep")) # ndraws x N_items

# 12.2) Allinea item verso "negative affect": inverti happy(1) e satisfied(3)
invert <- sd$item_id %in% c(1L, 3L)
Yrep_na <- Yrep
Yrep_na[, invert] <- 8 - Yrep_na[, invert]

# 12.3) Media per osservazione (beep) = media dei 4 item
#       t(Yrep_na) è (N_items x ndraws); rowsum somma per obs_id; poi /4
NA_obs_t <- rowsum(t(Yrep_na), group = sd$obs_id, reorder = FALSE) / 4 # (N_obs x ndraws)

# 12.4) Medie per soggetto in PRE (period==2) e POST (period==3)
I <- sd$I
is_pre <- sd$period == 2L
is_post <- sd$period == 3L
ndraws <- ncol(NA_obs_t)

avg_by_subj <- function(rows_flag) {
  out <- matrix(NA_real_, nrow = I, ncol = ndraws)
  idx <- sd$subject[rows_flag]
  sums <- rowsum(
    NA_obs_t[rows_flag, , drop = FALSE],
    group = idx,
    reorder = FALSE
  )
  cnt <- tabulate(idx, nbins = I)
  present <- which(cnt > 0)
  out[present, ] <- sums[present, , drop = FALSE] / cnt[present]
  out # (I x ndraws)
}

PRE_t <- avg_by_subj(is_pre) # (I x ndraws)
POST_t <- avg_by_subj(is_post) # (I x ndraws)

# 12.5) Differenza pre - post per soggetto e draw → (ndraws x I)
DIFF_t <- PRE_t - POST_t # (I x ndraws)
diff_rep_mat <- t(DIFF_t) # (ndraws x I)

# 12.6) Costruisci la proxy osservata y_obs dagli item reali (stesso criterio)
df_items <- tibble(
  n = sd$obs_id,
  subj = sd$subject[sd$obs_id],
  per = sd$period[sd$obs_id],
  item = sd$item_id,
  y = sd$y_item
) %>%
  mutate(y_na = if_else(item %in% c(1, 3), 8L - y, y))

m_obs <- df_items %>%
  group_by(n, subj, per) %>%
  summarise(na_mean = mean(y_na), .groups = "drop")

y_obs <- m_obs %>%
  group_by(subj) %>%
  summarise(
    pre = if (any(per == 2)) mean(na_mean[per == 2]) else NA_real_,
    post = if (any(per == 3)) mean(na_mean[per == 3]) else NA_real_
  ) %>%
  mutate(diff = pre - post) %>%
  arrange(subj) %>%
  pull(diff)

# 12.7) Limita a soggetti con pre & post; allinea le colonne di diff_rep_mat
keep <- is.finite(y_obs)
y_obs <- y_obs[keep]
diff_rep_mat <- diff_rep_mat[, keep, drop = FALSE]

# 12.8) PPC: scatter e intervalli predittivi
ppc_scatter_avg(y = y_obs, yrep = diff_rep_mat)
ppc_intervals(y = y_obs, yrep = diff_rep_mat, prob = 0.5, prob_outer = 0.9)

# ------------------------------------------------------------
# 13) (Facoltativo) Oggetti riassuntivi per report
#     WHY: utile se vuoi serializzare risultati chiave o costruire tabelle
# ------------------------------------------------------------
results_summary <- list(
  data_summary = list(
    n_subjects = sd$I,
    n_observations = sd$N_obs,
    n_items = sd$N_items,
    n_female = sum(sd$female),
    n_male = sum(1 - sd$female)
  ),
  bayes_r2 = list(
    baseline_only = c(mean = mean(R2_A), sd = sd(R2_A)),
    baseline_ema = c(mean = mean(R2_B), sd = sd(R2_B)),
    delta = c(
      mean = mean(delta_R2),
      l = quantile(delta_R2, .025),
      u = quantile(delta_R2, .975),
      p_pos = mean(delta_R2 > 0)
    )
  ),
  loo = list(
    obs = list(
      delta_elpd = delta_elpd_obs,
      se = se_diff_obs,
      winner = rownames(cmp_obs)[1]
    ),
    subj = list(
      delta_elpd = delta_elpd_subj,
      se = se_diff_subj,
      winner = rownames(cmp_subj)[1]
    )
  ),
  key_coeffs = list(
    gender = c(
      mean = mean(b_female),
      l = quantile(b_female, .025),
      u = quantile(b_female, .975),
      p_pos = mean(b_female > 0)
    ),
    var_ema = c(
      mean = mean(b_var_ema),
      l = quantile(b_var_ema, .025),
      u = quantile(b_var_ema, .975),
      p_pos = mean(b_var_ema > 0)
    ),
    var_gxsex = c(
      mean = mean(b_var_female),
      l = quantile(b_var_female, .025),
      u = quantile(b_var_female, .975),
      p_pos = mean(b_var_female > 0)
    )
  )
)

# (se vuoi: saveRDS(results_summary, "results_summary.rds"))

# ------------------------------------------------------------
# 14) DIAGNOSTICA SOGGETTO-LEVEL: chi “fittta” peggio e perché
#     - Costruiamo y_obs (diff pre-post) dagli item osservati
#     - Costruiamo diff_rep_mat (stessa statistica dalle repliche)
#     - Calcoliamo misure di misfit per soggetto: mu_pred, sd_pred, z, ppp
#     - Tabella ordinata e grafici rapidi
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(posterior)
library(ggplot2)


# ----- 14.1) y_obs: proxy osservata della differenza pre–post per soggetto -----
df_obs <- tibble(
  n = sd$obs_id,
  subj = sd$subject[sd$obs_id],
  per = sd$period[sd$obs_id],
  item = sd$item_id,
  y = sd$y_item
) |>
  mutate(
    # allinea tutti verso "negative affect"
    y_na = if_else(item %in% c(1L, 3L), 8L - y, y)
  )

# media NA per beep (4 item)
m_obs <- df_obs |>
  group_by(n, subj, per) |>
  summarise(na_mean = mean(y_na), .groups = "drop")

# differenza pre - post (per soggetto)
y_obs_df <- m_obs |>
  group_by(subj) |>
  summarise(
    pre = if (any(per == 2L)) mean(na_mean[per == 2L]) else NA_real_,
    post = if (any(per == 3L)) mean(na_mean[per == 3L]) else NA_real_
  ) |>
  mutate(diff = pre - post) |>
  arrange(subj)

keep <- is.finite(y_obs_df$diff)
y_obs <- y_obs_df$diff[keep]
subj_keep <- y_obs_df$subj[keep]

# ----- 14.2) diff_rep_mat: stessa statistica dai dati replicati -----
# Draws degli item replicati: (ndraws x N_items)
Yrep <- posterior::as_draws_matrix(fit_B$draws("y_item_rep"))

# Allinea item verso NA (inverti happy=1 e satisfied=3)
invert <- sd$item_id %in% c(1L, 3L)
Yrep_na <- Yrep
Yrep_na[, invert] <- 8 - Yrep_na[, invert]

# Media per osservazione (beep): t(Yrep_na) è (N_items x ndraws)
# Sommiamo per obs_id con rowsum, poi dividiamo per 4 item
NA_obs_t <- rowsum(t(Yrep_na), group = sd$obs_id, reorder = FALSE) / 4 # (N_obs x ndraws)

# Medie per soggetto in PRE e POST (no loop)
I <- sd$I
is_pre <- sd$period == 2L
is_post <- sd$period == 3L
ndraws <- ncol(NA_obs_t)

avg_by_subj <- function(rows_flag) {
  sums <- rowsum(
    NA_obs_t[rows_flag, , drop = FALSE],
    group = sd$subject[rows_flag],
    reorder = FALSE
  ) # (S x ndraws), S≈I
  cnt <- tabulate(sd$subject[rows_flag], nbins = I) # beep per soggetto
  out <- sweep(sums, 1, ifelse(cnt == 0, NA_real_, cnt), "/")
  # garantisci dimensione I x ndraws
  if (nrow(out) < I) {
    tmp <- matrix(NA_real_, nrow = I, ncol = ndraws)
    tmp[as.integer(rownames(sums)), ] <- out
    out <- tmp
  }
  out
}

PRE_t <- avg_by_subj(is_pre) # (I x ndraws)
POST_t <- avg_by_subj(is_post) # (I x ndraws)

# Differenza pre - post per soggetto e draw → (ndraws x I)
diff_rep_mat <- t(PRE_t - POST_t)
diff_rep_mat <- diff_rep_mat[, keep, drop = FALSE] # allinea a y_obs

# ----- 14.3) Metriche di misfit per soggetto -----
# z-residuals (con protezione per sd = 0)
mu_pred <- apply(diff_rep_mat, 2, mean) # length = I_keep
sd_pred <- apply(diff_rep_mat, 2, sd) # length = I_keep
sd_pred[sd_pred == 0 | !is.finite(sd_pred)] <- NA_real_
z_resid <- (y_obs - mu_pred) / sd_pred

# Posterior predictive p-value bilaterale (correzione orientamento)
y_mat <- matrix(
  y_obs,
  nrow = nrow(diff_rep_mat),
  ncol = ncol(diff_rep_mat),
  byrow = TRUE
)
p_ge <- colMeans(diff_rep_mat >= y_mat) # length = I_keep
p_le <- colMeans(diff_rep_mat <= y_mat) # length = I_keep
ppp <- 2 * pmin(p_ge, p_le)
ppp[ppp > 1] <- 1

# Metadati utili per lettura
n_pre <- tabulate(sd$subject[sd$period == 2L], nbins = I)
n_post <- tabulate(sd$subject[sd$period == 3L], nbins = I)

diag_tbl <- tibble::tibble(
  subj = subj_keep,
  female = sd$female[subj_keep],
  n_pre = n_pre[subj_keep],
  n_post = n_post[subj_keep],
  y_obs = y_obs,
  mu_pred = mu_pred,
  sd_pred = sd_pred,
  z = z_resid,
  ppp = ppp
)

# --- Lookup subj -> user_id (coerente con prepare_* ) ---
lookup_uid <- d %>%
  transmute(
    subj = as.integer(factor(user_id)), # stesso mapping usato in prepare_*
    user_id = user_id
  ) %>%
  distinct(subj, user_id) %>%
  arrange(subj)

# Arricchisci la diagnostica con user_id
diag_tbl <- diag_tbl %>%
  left_join(lookup_uid, by = "subj") %>%
  relocate(user_id, .after = subj)

# Stampa le top list con user_id
worst_by_absz <- diag_tbl %>% arrange(desc(abs(z)))
worst_by_ppp <- diag_tbl %>% arrange(ppp)

cat("\n--- Top 10 peggiori per |z| (con user_id) ---\n")
print(head(worst_by_absz, 10))

cat("\n--- Top 10 peggiori per PPP (con user_id) ---\n")
print(head(worst_by_ppp, 10))

cat("\nSoggetti peggiori (|z|):\n")
print(worst_by_absz$user_id[1:10])

cat("\nSoggetti peggiori (PPP):\n")
print(worst_by_ppp$user_id[1:10])


# ----- 14.4) Grafici rapidi -----
# Caterpillar degli z-residuals
ggplot(
  diag_tbl,
  aes(x = reorder(factor(subj), z), y = z, color = factor(female))
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  coord_flip() +
  labs(
    x = "Soggetto",
    y = "z-residual (y_obs - mu_pred) / sd_pred",
    color = "Female",
    title = "Posterior Predictive z-residuals per soggetto"
  ) +
  theme_minimal()

# PPP vs |z|
ggplot(diag_tbl, aes(x = abs(z), y = ppp)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  labs(
    x = "|z|",
    y = "Posterior predictive p-value (bilaterale)",
    title = "Misfit per soggetto: |z| vs PPP"
  ) +
  theme_minimal()

# (Opzionale) Esporta tabella diagnostica
# write.csv(diag_tbl, file = "diagnostica_soggetti_ppc.csv", row.names = FALSE)

# -------------------------------
# 1) Setup: metriche e etichette
# -------------------------------
comp_vars <- c(
  "n_total",
  "n_per_2",
  "n_per_3",
  "sd_within_1_7",
  "pct_same_consec_1_7",
  "max_run_len_1_7",
  "pct_extremes_1_7",
  "entropy_1_7",
  "pct_mult10",
  "diff_pre_post"
)

metric_labels <- tibble::tibble(
  metric = comp_vars,
  label = c(
    "Beep totali",
    "Beep PRE",
    "Beep POST",
    "SD(NA 1–7)",
    "% risposte identiche consecutive",
    "Run max identiche",
    "% estremi (1 o 7)",
    "Entropia (1–7)",
    "% multipli di 10 (0–100)",
    "Reattività (PRE–POST)"
  ),
  higher_is = c(
    "neutro",
    "meglio",
    "meglio", # più dati è meglio
    "neutro",
    "peggio",
    "peggio", # tanta ripetizione = peggio
    "peggio",
    "meglio",
    "peggio", # tanti estremi/mult10 = peggio; più entropia = meglio
    "neutro" # dipende dal contesto
  )
)

# -------------------------------
# 2) Riassunti per gruppo
# -------------------------------
summarize_var <- function(x)
  c(
    n = sum(is.finite(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q25 = quantile(x, .25, na.rm = TRUE),
    q75 = quantile(x, .75, na.rm = TRUE)
  )

diag_flag2 <- diag_flag %>% mutate(group = "FLAGGED")
diag_ok2 <- diag_ok %>% mutate(group = "OTHERS")

long_comp <- bind_rows(diag_flag2, diag_ok2) %>%
  select(group, all_of(comp_vars)) %>%
  pivot_longer(-group, names_to = "metric", values_to = "value")

desc_comp <- long_comp %>%
  group_by(group, metric) %>%
  summarise(
    n = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, .25, na.rm = TRUE),
    q75 = quantile(value, .75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(metric_labels, by = "metric") %>%
  select(metric, label, group, n, mean, sd, median, q25, q75)

# -------------------------------
# 3) Effect size + p-value per metrica
# -------------------------------
hedges_g <- function(x, y) {
  nx <- sum(is.finite(x))
  ny <- sum(is.finite(y))
  mx <- mean(x, na.rm = TRUE)
  my <- mean(y, na.rm = TRUE)
  sx <- sd(x, na.rm = TRUE)
  sy <- sd(y, na.rm = TRUE)
  sp <- sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / max(1, (nx + ny - 2)))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  J <- 1 - (3 / (4 * (nx + ny) - 9))
  J * (mx - my) / sp
}

es_pval <- map_dfr(comp_vars, function(v) {
  x <- diag_flag[[v]]
  y <- diag_ok[[v]]
  p <- tryCatch(
    wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  tibble(
    metric = v,
    hedges_g = hedges_g(x, y),
    p_wilcox = p,
    mean_flagged = mean(x, na.rm = TRUE),
    mean_others = mean(y, na.rm = TRUE),
    sd_flagged = sd(x, na.rm = TRUE),
    sd_others = sd(y, na.rm = TRUE)
  )
})

nice_tbl <- es_pval %>%
  left_join(metric_labels, by = "metric") %>%
  mutate(
    g_dir = case_when(
      hedges_g > 0.5 ~ "↑ (moderato+)",
      hedges_g < -0.5 ~ "↓ (moderato+)",
      abs(hedges_g) <= 0.5 & abs(hedges_g) > 0.2 ~
        if_else(hedges_g >= 0, "↑ (piccolo)", "↓ (piccolo)"),
      TRUE ~ "≈ (trascurabile)"
    ),
    p_sig = case_when(
      is.na(p_wilcox) ~ "NA",
      p_wilcox < .001 ~ "<.001",
      p_wilcox < .01 ~ "<.01",
      p_wilcox < .05 ~ "<.05",
      TRUE ~ sprintf("%.3f", p_wilcox)
    )
  ) %>%
  select(
    Metrica = label,
    metric,
    `Media FLAGGED` = mean_flagged,
    `SD FLAGGED` = sd_flagged,
    `Media OTHERS` = mean_others,
    `SD OTHERS` = sd_others,
    `Hedges g` = hedges_g,
    `Direzione` = g_dir,
    `p (Wilcoxon)` = p_sig
  ) %>%
  arrange(match(metric, comp_vars)) %>%
  select(-metric)

cat("\n=== Confronto FLAGGED vs OTHERS (riassunto per metrica) ===\n")
data.frame(nice_tbl, n = nrow(nice_tbl), digits = 3)

# -------------------------------
# 4) “Red flags” per soggetto
#    (soglie ragionevoli: adatta se vuoi)
# -------------------------------
rules <- list(
  few_post = ~ n_per_3 < 2, # pochi beep post
  low_sd = ~ sd_within_1_7 < 0.5, # variabilità molto bassa
  high_mult10 = ~ pct_mult10 > 0.60, # moltissimi multipli di 10
  extreme_rx = ~ abs(diff_pre_post) > 3 # reattività molto grande (scala 1–7)
)

apply_rule <- function(df, fn) rlang::as_function(fn)(df)

# -------------------------------
# 4) “Red flags” per soggetto — versione semplice e robusta
# -------------------------------
# -------------------------------
# 4) “Red flags” per soggetto — versione semplice e robusta (FIX)
# -------------------------------
flagged_redflags <- diag_flag %>%
  mutate(
    few_post = coalesce(n_per_3, 0L) < 2, # pochi beep post
    low_sd = sd_within_1_7 < 0.5, # variabilità molto bassa
    high_mult10 = pct_mult10 > 0.60, # moltissimi multipli di 10
    extreme_rx = abs(diff_pre_post) > 3, # reattività molto grande (scala 1–7)
    n_flags = few_post + low_sd + high_mult10 + extreme_rx
  ) %>%
  arrange(
    desc(n_flags),
    desc(abs(diff_pre_post)),
    desc(pct_mult10),
    sd_within_1_7
  ) %>%
  select(
    subj,
    user_id,
    n_total,
    n_per_2,
    n_per_3,
    sd_within_1_7,
    pct_mult10,
    diff_pre_post,
    few_post,
    low_sd,
    high_mult10,
    extreme_rx,
    n_flags
  )

print(flagged_redflags, n = nrow(flagged_redflags))

cat("\n=== Red flags per soggetto (FLAGGED) ===\n")
data.frame(flagged_redflags, n = nrow(flagged_redflags), digits = 3)

# -------------------------------
# 5) Mini-grafici (confronto gruppi)
# -------------------------------
plot_vars <- c(
  "n_per_2",
  "n_per_3",
  "sd_within_1_7",
  "pct_mult10",
  "diff_pre_post"
)
plot_labels <- metric_labels %>% filter(metric %in% plot_vars)

plot_df <- bind_rows(diag_flag2, diag_ok2) %>%
  select(group, all_of(plot_vars)) %>%
  pivot_longer(-group, names_to = "metric", values_to = "value") %>%
  left_join(plot_labels, by = "metric")

ggplot(plot_df, aes(x = group, y = value)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = .15) +
  facet_wrap(~label, scales = "free_y") +
  labs(
    x = NULL,
    y = "Media (IC ~95%)",
    title = "FLAGGED vs OTHERS: confronto sintetico per metrica"
  ) +
  theme_minimal(base_size = 12)
