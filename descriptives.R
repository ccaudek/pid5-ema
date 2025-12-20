df_pid5 <- d %>%
  dplyr::select(starts_with("pid5_") & ends_with("baseline"))


summary(df_pid5)


library(dplyr)
library(purrr)

# 1) Numero di item previsti per ciascuna facet dopo esclusioni
facet_items_clean <- lapply(
  facet_items,
  function(pos) setdiff(pos, exclude_qpos)
)
facet_n_expected <- vapply(facet_items_clean, length, integer(1))

# 2) Trasforma le facet-sum in facet-mean (0–3) con regola >=50%
for (f in names(facet_items)) {
  sum_col <- paste0("facet_", f) # già nel tuo df_scores
  nvalid_col <- paste0("facet_", f, "_nvalid") # già nel tuo df_scores
  mean_col <- paste0("facet_", f, "_mean")

  need <- facet_n_expected[[f]]
  thresh <- ceiling(0.5 * need) # soglia 50%

  df_scores[[mean_col]] <- with(
    df_scores,
    ifelse(
      get(nvalid_col) >= thresh,
      get(sum_col) / get(nvalid_col), # media item (0–3)
      NA_real_
    )
  )
}

# 3) Definisci i domini come media delle 3 facet-mean (richiedi ≥2 facet valide)
domain_to_facets <- domain_facets # già nel tuo script

domain_mean_cols <- character(0)
for (dom in names(domain_to_facets)) {
  facet_means <- paste0("facet_", domain_to_facets[[dom]], "_mean")
  out_col <- paste0("domain_", dom, "_mean") # scala 0–3

  df_scores[[out_col]] <- apply(
    df_scores[, facet_means, drop = FALSE],
    1,
    function(x) if (sum(!is.na(x)) >= 2) mean(x, na.rm = TRUE) else NA_real_
  )
  domain_mean_cols <- c(domain_mean_cols, out_col)
}

# 4) Output "pulito" per baseline (rinomina come vuoi)
df_out_mean <- df_scores %>%
  select(user_id, all_of(domain_mean_cols)) %>%
  rename(
    pid5_negative_affect_baseline = domain_negative_affect_mean,
    pid5_detachment_baseline = domain_detachment_mean,
    pid5_antagonism_baseline = domain_antagonism_mean,
    pid5_disinhibition_baseline = domain_disinhibition_mean,
    pid5_psychoticism_baseline = domain_psychoticism_mean
  )

# 5) Medie di campione (0–3) da riportare
sample_means <- df_out_mean %>%
  summarise(across(starts_with("pid5_"), ~ mean(.x, na.rm = TRUE)))

print(sample_means)


"happy"
"sad"
"satisfied"
"angry"

"day"
"hour"
"calendar_day"
"bysubj_day"

"context_threat"
"tripm_boldness"
"tripm_meanness"

"boldness_baseline"
"meanness_baseline"


"dass_stress_baseline"
"dass_anxiety_baseline"
"dass_depression_baseline"

"sex"
"exam_period"

"user_id"


cols <- c(
  "happy",
  "sad",
  "satisfied",
  "angry",
  "day",
  "hour",
  "calendar_day",
  "bysubj_day",
  "context_threat",
  "tripm_boldness",
  "tripm_meanness",
  "boldness_baseline",
  "meanness_baseline",
  "dass_stress_baseline",
  "dass_anxiety_baseline",
  "dass_depression_baseline",
  "sex",
  "exam_period",
  "user_id",
  "dass_stress",
  "dass_depression",
  "dass_anxiety"
)

# Avviso per eventuali colonne mancanti
missing <- setdiff(cols, names(d))
if (length(missing)) {
  warning("Colonne non trovate in d: ", paste(missing, collapse = ", "))
}

# Selezione (non fallisce se alcune mancano)
d_sel <- d %>% dplyr::select(dplyr::any_of(cols))
