# Crea variabili piecewise
df_piecewise <- df_analysis %>%
  mutate(
    # Phase 1: Stress (baseline=0, pre=1, post=1)
    stress_phase = ifelse(timepoint %in% c("pre", "post"), 1, 0),
    # Phase 2: Recovery (baseline=0, pre=0, post=1)
    recovery_phase = ifelse(timepoint == "post", 1, 0)
  )

formula_piecewise <- bf(
  f0_mean_u ~
    stress_phase *
      (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c) +
      recovery_phase * (pid5_negative_affectivity_bl_c + pid5_detachment_bl_c) +
      (1 + stress_phase + recovery_phase | ID)
)

m_piecewise <- brm(
  formula_piecewise,
  data = df_piecewise,
  family = gaussian(),
  prior = priors_moderation,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4
)
pp_check(m_piecewise)
summary(m_piecewise)
