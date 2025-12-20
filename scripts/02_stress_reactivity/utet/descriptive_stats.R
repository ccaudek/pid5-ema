# ------------------------------------------------------------
# Descrittive: Negative Affect per periodo e per genere (dataset: d2)
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

# --- helper per portare gli item su scala 1..7 se servisse ---
as_item_1_7 <- function(x) {
  x <- suppressWarnings(readr::parse_number(as.character(x)))
  if (all(is.na(x))) return(as.integer(x))
  if (all(x[is.finite(x)] %in% 1:7)) return(as.integer(x))
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  y <- 1 + 6 * (x - xmin) / (xmax - xmin)
  as.integer(round(pmax(1, pmin(7, y))))
}

# --- 1) Selezione/ricostruzione degli item 1..7 ---
has_dot_items <- all(
  c(".__happy__", ".__satisfied__", ".__sad__", ".__angry__") %in% names(d2)
)
if (has_dot_items) {
  happy_1_7 <- d2$.__happy__
  satis_1_7 <- d2$.__satisfied__
  sad_1_7 <- d2$.__sad__
  angry_1_7 <- d2$.__angry__
} else {
  # partiamo dai campi "happy/satisfied/sad/angry" (0-100) e li riportiamo a 1..7
  happy_1_7 <- as_item_1_7(d2$happy)
  satis_1_7 <- as_item_1_7(d2$satisfied)
  sad_1_7 <- as_item_1_7(d2$sad)
  angry_1_7 <- as_item_1_7(d2$angry)
}

# --- 2) Periodo (Baseline/Pre/Post) e Soggetto ---
# Se ci sono le variabili ".__per__" (1=baseline, 2=pre, 3=post) e ".__subj__", usiamole.
if ("__.__per__" %in% names(d2)) {
  per_fac <- factor(
    d2$.__per__,
    levels = c(1, 2, 3),
    labels = c("Baseline", "Pre-exam", "Post-exam")
  )
} else {
  per_raw <- tolower(as.character(d2$exam_period))
  per_fac <- dplyr::case_when(
    per_raw %in% c("baseline", "base", "t0", "bl") ~ "Baseline",
    per_raw %in% c("pre", "pre_exam", "pre-exam", "preexam") ~ "Pre-exam",
    per_raw %in% c("post", "post_exam", "post-exam", "postexam") ~ "Post-exam",
    TRUE ~ NA_character_
  ) |>
    factor(levels = c("Baseline", "Pre-exam", "Post-exam"))
}

if ("__.__subj__" %in% names(d2)) {
  subj <- as.integer(d2$.__subj__)
} else {
  subj <- as.integer(factor(d2$user_id))
}

# --- 3) Genere binario (Male/Female) a partire da 'sex' ---
sex_raw <- tolower(as.character(d2$sex))
female <- as.integer(
  sex_raw %in% c("f", "female", "donna", "femmina", "woman", "1", "true")
)
gender <- factor(
  ifelse(female == 1, "Female", "Male"),
  levels = c("Male", "Female")
)

# --- 4) Costruzione del punteggio di Negative Affect (composito 1..7) ---
#    NA_score = (7 - happy) + (7 - satisfied) + sad + angry
NA_score <- (7 - happy_1_7) + (7 - satis_1_7) + sad_1_7 + angry_1_7

df <- tibble::tibble(
  subj = subj,
  period = per_fac,
  gender = gender,
  NA_score = as.numeric(NA_score)
) |>
  drop_na(period, gender, NA_score)

# --- 5) Media per soggetto e periodo (per non far pesare i soggetti con più EMA) ---
sp_means <- df |>
  group_by(subj, period, gender) |>
  summarise(NA_mean = mean(NA_score), .groups = "drop")

# --- 6) Descrittive globali per periodo (collassando i soggetti) ---
per_summ <- sp_means |>
  group_by(period) |>
  summarise(
    n_subj = n(),
    mean = mean(NA_mean),
    sd = sd(NA_mean),
    se = sd / sqrt(n_subj),
    l95 = mean - 1.96 * se,
    u95 = mean + 1.96 * se,
    .groups = "drop"
  )

# --- 7) Descrittive per periodo × genere ---
per_gender_summ <- sp_means |>
  group_by(period, gender) |>
  summarise(
    n_subj = n(),
    mean = mean(NA_mean),
    sd = sd(NA_mean),
    se = sd / sqrt(n_subj),
    l95 = mean - 1.96 * se,
    u95 = mean + 1.96 * se,
    .groups = "drop"
  )

# --- 8) Figure (medie con CI 95%) ---
p_period <- ggplot(per_summ, aes(x = period, y = mean)) +
  geom_pointrange(
    aes(ymin = l95, ymax = u95),
    position = position_dodge(width = 0.3)
  ) +
  labs(
    title = "Negative affect per periodo (media su soggetto)",
    x = NULL,
    y = "Negative Affect (composito 1–7)"
  ) +
  theme_minimal(base_size = 12)

p_period_gender <- ggplot(
  per_gender_summ,
  aes(x = period, y = mean, color = gender, group = gender)
) +
  geom_pointrange(
    aes(ymin = l95, ymax = u95),
    position = position_dodge(width = 0.35)
  ) +
  labs(
    title = "Negative affect per periodo × genere (media su soggetto)",
    x = NULL,
    y = "Negative Affect (composito 1–7)",
    color = "Genere"
  ) +
  theme_minimal(base_size = 12)

# Stampa a video
print(per_summ)
print(per_gender_summ)
print(p_period)
print(p_period_gender)

# (opzionale) salvataggi
# ggsave("figs/na_by_period.png", p_period, width = 6.5, height = 4.5, dpi = 300)
# ggsave("figs/na_by_period_gender.png", p_period_gender, width = 7, height = 4.5, dpi = 300)
