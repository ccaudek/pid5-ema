# Overview ----------------------------------------------------------------
# Associated project: PID-5 and EMA
# Script purpose: Fit a continuous-time model.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-07-26
# Last update:
# Status: In progress
# Notes:
# The scripts are provided here: https://osf.io/tk3pf

# Load necessary libraries ------------------------------------------------

# devtools::install_github("ryanoisin/ctnet")

if (!requireNamespace("pacman")) install.packages("pacman")

pacman::p_load(
  ctsem,
  tibble,
  qgraph,
  RColorBrewer,
  ctnet,
  dplyr,
  igraph,
  ggraph,
  tidygraph,
  ggplot2
)


# Import data --------------------------------------------------------------
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


# Data wrangling -----------------------------------------------------------
# Format time portion of data for continuous-time models
ema_raw$day <- as.POSIXct(ema_raw$day, format = "%m/%d/%Y %H:%M")

filtered_data <- ema_raw %>%
  group_by(user_id) %>%
  mutate(time_hours = as.numeric(difftime(day, min(day), units = "hours"))) %>%
  ungroup() %>%
  mutate(
    user_id = as.character(user_id),
    time_hours = as.numeric(time_hours)
  )

filtered_data_ctsem <- filtered_data %>%
  dplyr::select(
    user_id,
    time_hours,
    pid5_negative_affectivity,
    pid5_detachment,
    pid5_antagonism,
    pid5_disinhibition,
    pid5_psychoticism
  )

filtered_data_ctsem <- filtered_data_ctsem %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    user_id = as.character(user_id),
    time_hours = as.numeric(time_hours)
  ) %>%
  as.data.frame()

filtered_data_ctsem_std <- filtered_data_ctsem %>%
  mutate(across(
    .cols = c(
      pid5_negative_affectivity,
      pid5_detachment,
      pid5_antagonism,
      pid5_disinhibition,
      pid5_psychoticism
    ),
    .fns = ~ as.numeric(scale(.x))
  ))


# Continuous-time model estimation -----------------------------------------
set.seed(123)
# make variable labels
lab <- c("NA", "DE", "AN", "DI", "PS")
p <- 5

#Estimate Continuous-time model
ctmodel <- ctModel(
  type = 'ct',
  manifestNames = c(
    "pid5_negative_affectivity",
    "pid5_detachment",
    "pid5_antagonism",
    "pid5_disinhibition",
    "pid5_psychoticism"
  ),
  latentNames = paste0(c("NA", "DE", "AN", "DI", "PS"), "_l"),
  LAMBDA = diag(5),
  DRIFT = "auto",
  MANIFESTMEANS = matrix(0, nrow = 5, ncol = 1),
  MANIFESTVAR = diag(0, 5),
  CINT = "auto",
  DIFFUSION = "auto"
)

ctmodel$subjectIDname <- "user_id"
ctmodel$timeName <- "time_hours"


ctfit_vi <- ctStanFit(
  datalong = filtered_data_ctsem_std,
  ctstanmodel = ctmodel,
  vb = TRUE, # <- qui specifichi l'uso della variational Bayes
  cores = 8,
  iter = 20000 # aumenta il numero di iterazioni per una stima più stabile
)

# saveRDS(ctfit_vi, file = "ctfit_vi.rds")
# ctfit_vi <- readRDS(here::here("scripts", "ctfit_vi.rds"))

summary_vi <- summary(ctfit_vi)
# Estrai le righe relative alla DRIFT matrix
drift_df <- subset(summary_vi$parmatrices, matrix == "DRIFT")
head(drift_df)

make_drift_table_from_parmatrices <- function(
  drift_df,
  digits = 2,
  var_names = NULL
) {
  p <- length(unique(drift_df$row)) # dimensione della matrice
  if (is.null(var_names)) var_names <- paste0("V", 1:p)

  drift_table <- matrix("", nrow = p, ncol = p)
  rownames(drift_table) <- var_names
  colnames(drift_table) <- var_names

  for (i in 1:nrow(drift_df)) {
    r <- drift_df$row[i]
    c <- drift_df$col[i]
    m <- round(drift_df$Mean[i], digits)
    l <- round(drift_df$`2.5%`[i], digits)
    u <- round(drift_df$`97.5%`[i], digits)
    drift_table[r, c] <- paste0(m, " (", l, ", ", u, ")")
  }

  return(as.data.frame(drift_table))
}

drift_table <- make_drift_table_from_parmatrices(
  drift_df = drift_df,
  digits = 2,
  var_names = c("NA", "DE", "AN", "DI", "PS")
)

# * NA = Negative Affectivity
# * DE = Detachment
# * AN = Antagonism
# * DI = Disinhibition
# * PS = Psychoticism

print(drift_table)
# NA                      DE                   AN                   DI                   PS
# NA -7.93 (-7.93, -7.93)       1.07 (1.07, 1.07) -1.56 (-1.56, -1.56) -0.08 (-0.08, -0.08)    1.74 (1.74, 1.74)
# DE -0.73 (-0.73, -0.73) -10.23 (-10.23, -10.23)    0.85 (0.85, 0.85) -0.28 (-0.28, -0.28)    0.18 (0.18, 0.18)
# AN -0.99 (-0.99, -0.99)       -0.9 (-0.9, -0.9) -6.93 (-6.93, -6.93) -1.07 (-1.07, -1.07) -0.97 (-0.97, -0.97)
# DI    -0.8 (-0.8, -0.8)    -1.05 (-1.05, -1.05) -0.08 (-0.08, -0.08) -7.67 (-7.67, -7.67) -1.87 (-1.87, -1.87)
# PS    1.11 (1.11, 1.11)       1.26 (1.26, 1.26) -1.08 (-1.08, -1.08) -1.95 (-1.95, -1.95) -6.97 (-6.97, -6.97)

# Drift matrix and local dependency interpretation

#' The drift matrix parameters are summarized in Table 1, and represent a local
#' dependency network among the five PID-5 domains (Negative Affectivity,
#' Detachment, Antagonism, Disinhibition, and Psychoticism). All parameters were
#' estimated using continuous-time dynamic modeling with variational inference.
#'
#' We first examined the **auto-effects** (i.e., the diagonal elements of the
#' drift matrix), which indicate how strongly each variable returns to its own
#' equilibrium after a perturbation. All auto-effects were negative, consistent
#' with a *dynamically stable system*. However, the magnitude of these
#' self-regulation tendencies varied considerably.
#'
#' Interpretazione clinica: la regolazione rapida può indicare bassa cronicità
#' nei tratti maladattivi — un segnale positivo. Suggerisce che, pur con alcune
#' vulnerabilità, questi studenti non sembrano mostrare pattern persistenti di
#' disregolazione tipici della psicopatologia conclamata.
#'
#' Specifically, **Negative Affectivity (−7.93)**, **Detachment (−10.23)**,
#' **Antagonism (−6.93)**, **Disinhibition (−7.67)**, and **Psychoticism (−6.97)**
#' all showed strong negative auto-effects, suggesting that fluctuations in
#' these traits are quickly regulated and tend to return to baseline levels
#' over time. Among these, **Detachment** exhibited the most rapid return to
#' equilibrium, while **Antagonism** and **Psychoticism** showed comparatively
#' slower decay rates.

### Cross-effects: dynamics of mutual influence

#' We then explored the **cross-effects**, which capture how the level of one
#' trait influences the *instantaneous rate of change* of another. These effects
#' characterize the dynamic interactions between traits within individuals over
#' time.
#'
#' Several notable effects emerged:
#' * **Negative Affectivity** was *positively influenced* by **Detachment (1.07)**
#' and **Psychoticism (1.74)**, suggesting that increases in these two traits
#' tend to *accelerate the growth* of negative affectivity. In contrast,
#' **Antagonism (−1.56)** had a suppressive effect on Negative Affectivity,
#' acting as a regulatory force.
#'
#' * **Detachment** was increased by **Antagonism (0.85)** but decreased by
#' **Negative Affectivity (−0.73)** and **Disinhibition (−0.28)**. These results
#' indicate that antagonistic traits might promote social withdrawal, while
#' emotional distress and impulsivity may buffer against it, perhaps by
#' activating interpersonal engagement.
#'
#' * **Antagonism** itself was inhibited by all other traits, especially
#' **Detachment (−0.90)** and **Disinhibition (−1.07)**. This suggests that
#' individuals experiencing greater withdrawal or disinhibition may show less
#' antagonistic tendencies over time.
#'
#' * **Disinhibition** was negatively affected by all other traits, including
#' **Psychoticism (−1.87)** and **Negative Affectivity (−0.80)**. These results
#' imply that increased emotional distress and perceptual disturbances might
#' contribute to a reduction in disinhibited behavior over time.
#'
#' * Interestingly, **Psychoticism** was *enhanced* by **Negative Affectivity (1.11)**
#' and **Detachment (1.26)**, and *reduced* by **Antagonism (−1.08)** and
#' **Disinhibition (−1.95)**. This profile suggests that internalizing traits
#' may fuel perceptual and cognitive disturbances, while more outward-directed
#' behaviors serve as a protective buffer.
#'
#' Psychoticism aumenta con NA (+1.11) e Detachment (+1.26)
#' → Gli studenti che vivono affetto negativo e ritiro sociale hanno maggiori
#' probabilità di mostrare esperienze percettive o cognitive atipiche (p. es.,
#' sospettosità, pensieri strani) nel breve termine.
#' → Questo è compatibile con un profilo di vulnerabilità psicotica non clinica,
#' come riscontrato in letteratura sugli studenti ad alto rischio (es.
#' schizotipia).
#'
#' Antagonism diminuisce con Detachment (−0.90) e Disinhibition (−1.07)
#' → Più una persona è ritirata o impulsiva, meno mostra atteggiamenti ostili o
#' provocatori.
#' → Potrebbe suggerire che gli stili di disregolazione nei contesti sociali
#' variano: alcuni studenti possono reagire al disagio ritirandosi, altri con
#' impulsività, ma non necessariamente con aggressività manifesta.
#'
#' Negative Affectivity aumenta con Detachment (+1.07) e Psychoticism (+1.74),
#' ma diminuisce con Antagonism (−1.56)
#' → Chi si isola o ha esperienze percettive disturbate tende a mostrare
#' maggiore distress affettivo, ma l’antagonismo sembra agire da “contenitore” o
#' da forza compensatoria (forse perché canalizza l’affetto negativo verso
#' l’esterno?).
#' Interpretazione clinica:
#' Questi pattern indicano forme diverse di vulnerabilità psicologica all'interno
#' del campione:
#' - Alcuni studenti sembrano mostrare profili di ritiro + affetto negativo, con
#' rischio per forme internalizzanti (ansia, depressione, tratti schizotipici).
#' - Altri manifestano impulsività e antagonismo, che potrebbero essere
#' espressione di strategie meno adattive ma transitorie, forse influenzate dal
#' contesto accademico/stress universitario.

### Summary of the dynamic structure

#' Together, the pattern of results describes a highly interconnected network of
#' psychopathological traits. Negative Affectivity appears particularly sensitive
#' to external influences, acting as a node where multiple traits converge.
#' Psychoticism is dynamically linked to both internal distress and disinhibition,
#' potentially making it a useful marker of general dysregulation.
#'
#' All traits demonstrated negative auto-effects, confirming the overall
#' *stationarity* and *regulatory tendencies* of the psychological system. Yet,
#' the strength of cross-links between traits suggests that deviations from
#' equilibrium in one domain may rapidly propagate to others, reinforcing the
#' importance of modeling these processes dynamically and idiographically.
#'
#' Conclusione integrativa
#' Il modello suggerisce che la psicopatologia nella vita quotidiana degli
#' studenti universitari non è statica né rigida, ma:
#' - si manifesta come fluttuazioni regolabili di tratti disfunzionali,
#' - mostra pattern dinamici che possono rivelare vulnerabilità latenti,
#' - non mostra però caratteristiche di disregolazione cronica, tipiche dei
#' disturbi clinici conclamati.
#'
#' Implicazioni applicative:
#' Questi dati supportano l’utilità di un monitoraggio dinamico (EMA) nella
#' prevenzione della psicopatologia in popolazioni non cliniche.
#' Le connessioni dinamiche suggeriscono che interventi mirati su NA o Detachment
#' potrebbero avere effetti a cascata positivi su altri domini.

# Estrai solo i parametri DRIFT
drift_df <- summary(ctfit_vi)$parmatrices
drift_only <- drift_df[drift_df$matrix == "DRIFT", ]

# Numero di variabili
p <- length(unique(drift_only$row))

# Crea matrice DRIFT
drift_matrix <- matrix(drift_only$Mean, nrow = p, byrow = TRUE)

# Nomi delle variabili
var_names <- c("NA", "DE", "AN", "DI", "PS")
colnames(drift_matrix) <- var_names
rownames(drift_matrix) <- var_names

# Costruisci la edge list
edge_list <- as.data.frame(as.table(drift_matrix))
colnames(edge_list) <- c("from", "to", "weight")

# Filtra per pesi "non nulli"
edge_list <- edge_list[edge_list$weight != 0, ]

# Aggiungi informazioni su direzione del segno
edge_list$sign <- ifelse(edge_list$weight > 0, "positive", "negative")
edge_list$abs_weight <- abs(edge_list$weight)

# Crea oggetto grafo
graph <- graph_from_data_frame(edge_list, directed = TRUE)

# Visualizzazione con ggraph
ggraph(graph, layout = "circle") +
  geom_edge_link(
    aes(width = abs_weight, color = sign),
    arrow = arrow(length = unit(3, 'mm')),
    end_cap = circle(4, 'mm')
  ) +
  geom_node_label(aes(label = name), size = 6) +
  scale_edge_color_manual(values = c("positive" = "red", "negative" = "blue")) +
  scale_edge_width(range = c(0.5, 2)) +
  theme_void() +
  labs(title = "Dynamic Network of Latent Traits (DRIFT Matrix)")


# Crea layout Fruchterman-Reingold
layout_fr <- layout_with_fr(graph, weights = abs(E(graph)$weight))

# Crea un data.frame con coordinate x, y e nomi dei nodi
layout_df <- as.data.frame(layout_fr)
colnames(layout_df) <- c("x", "y")
layout_df$name <- V(graph)$name

# Crea tidygraph (richiesto da ggraph)
graph_tbl <- as_tbl_graph(graph)

# GGRAPH con layout manuale
ggraph(graph_tbl, layout = "manual", x = layout_df$x, y = layout_df$y) +
  geom_edge_link(
    aes(width = abs_weight, color = sign),
    arrow = arrow(length = unit(3, 'mm')),
    end_cap = circle(4, 'mm')
  ) +
  geom_node_label(aes(label = name), size = 6) +
  scale_edge_color_manual(values = c("positive" = "red", "negative" = "blue")) +
  scale_edge_width(range = c(0.5, 2)) +
  theme_void() +
  labs(title = "Dynamic Network (Fruchterman-Reingold layout)")


netplot <- function(
  mat,
  greyscale = FALSE,
  maximum = .5,
  asize = 4,
  edge.labels = TRUE,
  edge.label.cex = 1.1,
  fade = FALSE,
  shape = "circle",
  labels = NULL,
  vsize = 18,
  esize = 10,
  layout = "circle"
) {
  # Numero di variabili
  p <- nrow(mat)

  # Etichette automatiche se non fornite
  if (is.null(labels)) labels <- rownames(mat)
  if (is.null(labels)) labels <- paste0("X", 1:p)

  # Line types
  m_lty <- matrix(1, p, p)
  m_lty[mat < 0] <- 2

  # Colori
  m_col <- matrix("blue", p, p)
  m_col[mat > 0] <- "firebrick2"

  # Plot
  qgraph::qgraph(
    t(mat),
    edge.color = t(m_col),
    layout = layout,
    directed = TRUE,
    edge.labels = edge.labels,
    edge.label.cex = edge.label.cex,
    lty = t(m_lty),
    vsize = vsize,
    esize = esize,
    asize = asize,
    mar = c(8, 8, 8, 8),
    fade = fade,
    shape = shape,
    maximum = maximum,
    labels = labels
  )
}


drift <- matrix(
  c(
    -7.93,
    1.07,
    -1.56,
    -0.08,
    1.74,
    -0.73,
    -10.23,
    0.85,
    -0.28,
    0.18,
    -0.99,
    -0.90,
    -6.93,
    -1.07,
    -0.97,
    -0.80,
    -1.05,
    -0.08,
    -7.67,
    -1.87,
    1.11,
    1.26,
    -1.08,
    -1.95,
    -6.97
  ),
  nrow = 5,
  byrow = TRUE
)

# Nomi variabili
lab <- c("NA", "DE", "AN", "DI", "PS")

# Visualizza il grafo
netplot(drift, labels = lab)

pdf("CTnet_emp.pdf", width = 5.5, height = 5.5)
netplot(drift, labels = lab)
dev.off()

# MCMC sampling ------------------------------------------------------------
# Richiede almento 3 giorni!

# ctfit_mcmc <- ctStanFit(
#   datalong = filtered_data_ctsem_std,
#   ctstanmodel = ctmodel,
#   optimize = FALSE,
#   cores = 8,
#   iter = 2000
# )
