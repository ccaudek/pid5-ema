# In R, fai questo:
library(tidyverse)
results <- read_csv("results/exports/02_moderation_results.csv")

# Filtra solo significativi
sig_results <- results %>%
  filter(significant == TRUE) %>%
  arrange(type, outcome, vowel)

# Conta per tipo
sig_results %>%
  count(type)

# Vedi interazioni significative
sig_results %>%
  filter(type == "Interaction") %>%
  select(parameter, outcome, vowel, estimate, ci_lower, ci_upper)
