# ==============================================================================
# DIAGNOSTIC_check_ema_data.R
# Quick diagnostic to verify EMA data structure for between-within analysis
# ==============================================================================

library(tidyverse)
library(rio)

cat("\n", rep("=", 70), "\n", sep = "")
cat("DIAGNOSTIC CHECK: EMA Data for Between-Within Analysis\n")
cat(rep("=", 70), "\n\n")

# Load data
d_ema <- rio::import("data/processed/ema_plus_scales_cleaned.csv")

# Rename ID column
if ("user_id" %in% names(d_ema)) {
  d_ema$ID <- d_ema$user_id
}

pid5_vars <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

# ==============================================================================
# CHECK 1: exam_period Coverage
# ==============================================================================

cat("CHECK 1: exam_period Distribution\n")
cat(rep("-", 70), "\n")

length(unique(d_ema$ID))


if ("exam_period" %in% names(d_ema)) {
  period_dist <- d_ema %>%
    filter(!is.na(exam_period)) %>%
    count(exam_period) %>%
    arrange(exam_period)

  print(period_dist)

  # By person
  period_by_person <- d_ema %>%
    filter(!is.na(exam_period), !is.na(pid5_negative_affectivity)) %>%
    group_by(ID, exam_period) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    group_by(ID) %>%
    summarise(
      n_periods = n(),
      has_baseline = "baseline" %in% exam_period,
      has_pre = "pre_exam" %in% exam_period,
      has_post = "post_exam" %in% exam_period,
      .groups = "drop"
    )

  cat("\nSubjects by period coverage:\n")
  cat("  All 3 periods:", sum(period_by_person$n_periods == 3), "subjects\n")
  cat("  Only 2 periods:", sum(period_by_person$n_periods == 2), "subjects\n")
  cat("  Only 1 period:", sum(period_by_person$n_periods == 1), "subjects\n")

  missing_periods <- period_by_person %>%
    summarise(
      missing_baseline = sum(!has_baseline),
      missing_pre = sum(!has_pre),
      missing_post = sum(!has_post)
    )

  cat("\nMissing periods:\n")
  cat("  No baseline:", missing_periods$missing_baseline, "subjects\n")
  cat("  No pre_exam:", missing_periods$missing_pre, "subjects\n")
  cat("  No post_exam:", missing_periods$missing_post, "subjects\n")
} else {
  cat("⚠ WARNING: No 'exam_period' variable found!\n")
  cat("Cannot proceed with between-within analysis.\n")
}

# ==============================================================================
# CHECK 2: Within-Person Variance in PID-5
# ==============================================================================

cat("\n", rep("-", 70), "\n")
cat("CHECK 2: Within-Person Variance in PID-5\n")
cat(rep("-", 70), "\n")

if ("exam_period" %in% names(d_ema)) {
  # Aggregate by person and period - FILTER NA FIRST
  pid5_by_period <- d_ema %>%
    dplyr::filter(
      !is.na(exam_period),
      exam_period %in% c("baseline", "pre_exam", "post_exam"),
      !is.na(pid5_negative_affectivity) # CRITICAL: Remove NA
    ) %>%
    group_by(ID, exam_period) %>%
    summarise(
      across(all_of(pid5_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # Check how many subjects have 2+ periods
  subjects_by_n_periods <- pid5_by_period %>%
    group_by(ID) %>%
    summarise(n_periods = n(), .groups = "drop") %>%
    count(n_periods)

  cat("\nSubjects by number of periods with PID-5 data:\n")
  print(subjects_by_n_periods)

  n_with_2plus <- subjects_by_n_periods %>%
    dplyr::filter(n_periods >= 2) %>%
    dplyr::summarise(n_with_2plus = sum(n), .groups = "drop") %>%
    dplyr::pull(n_with_2plus)

  if (length(n_with_2plus) == 0) n_with_2plus <- 0

  if (n_with_2plus == 0) {
    cat("\n✗ CRITICAL: NO subjects have 2+ timepoints!\n")
    cat("Between-within analysis is NOT POSSIBLE.\n")
    cat("All subjects have data in only 1 period.\n\n")
  } else {
    # For subjects with 2+ periods, calculate range
    within_person_range <- pid5_by_period %>%
      group_by(ID) %>%
      dplyr::filter(n() >= 2) %>% # Only those with 2+ periods
      summarise(
        across(
          all_of(pid5_vars),
          ~ max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE),
          .names = "range_{.col}"
        ),
        n_periods = n(),
        .groups = "drop"
      )

    cat(
      "\nWithin-person range statistics (N =",
      nrow(within_person_range),
      "subjects with 2+ periods):\n"
    )
    cat("(How much do individuals vary across timepoints?)\n\n")

    for (var in pid5_vars) {
      range_var <- paste0("range_", var)
      avg_range <- mean(within_person_range[[range_var]], na.rm = TRUE)
      max_range <- max(within_person_range[[range_var]], na.rm = TRUE)
      prop_zero <- mean(within_person_range[[range_var]] == 0, na.rm = TRUE)

      var_short <- str_replace(var, "pid5_", "")

      cat(sprintf("  %s:\n", var_short))
      cat(sprintf("    Mean range: %.3f\n", avg_range))
      cat(sprintf("    Max range: %.3f\n", max_range))
      cat(sprintf("    %% with zero variance: %.1f%%\n", prop_zero * 100))
    }

    # Overall assessment
    avg_ranges <- within_person_range %>%
      dplyr::select(starts_with("range_")) %>%
      summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "var", values_to = "avg_range")

    overall_range <- mean(avg_ranges$avg_range)

    cat("\nOVERALL ASSESSMENT:\n")
    if (overall_range > 1.0) {
      cat(
        "✓ GOOD: Substantial within-person variance (mean range = ",
        round(overall_range, 2),
        ")\n",
        sep = ""
      )
      cat("  Between-within analysis is feasible and informative.\n")
    } else if (overall_range > 0.5) {
      cat(
        "⚠ MODERATE: Some within-person variance (mean range = ",
        round(overall_range, 2),
        ")\n",
        sep = ""
      )
      cat("  Between-within analysis may show modest effects.\n")
    } else {
      cat(
        "✗ LOW: Minimal within-person variance (mean range = ",
        round(overall_range, 2),
        ")\n",
        sep = ""
      )
      cat("  Traits are very stable - trait-only model may be sufficient.\n")
    }
  }
}

# ==============================================================================
# CHECK 3: Sample 10 Subjects
# ==============================================================================

cat("\n", rep("-", 70), "\n")
cat("CHECK 3: Sample of 10 Subjects\n")
cat(rep("-", 70), "\n\n")

if ("exam_period" %in% names(d_ema)) {
  sample_ids <- pid5_by_period %>%
    group_by(ID) %>%
    filter(n() >= 2) %>% # Only those with 2+ periods
    ungroup() %>%
    distinct(ID) %>%
    slice_sample(n = min(10, n())) %>%
    dplyr::pull(ID)

  sample_data <- pid5_by_period %>%
    filter(ID %in% sample_ids) %>%
    select(ID, exam_period, pid5_negative_affectivity, pid5_detachment) %>%
    pivot_wider(
      names_from = exam_period,
      values_from = c(pid5_negative_affectivity, pid5_detachment),
      names_glue = "{.value}_{exam_period}"
    )

  cat("Negative Affectivity and Detachment for 10 random subjects:\n\n")
  print(sample_data, n = Inf)

  cat("\nLook for variation across columns (baseline/pre/post).\n")
  cat("If values are very similar → low within-person variance.\n")
  cat("If values differ substantially → good within-person variance.\n")
}

# ==============================================================================
# RECOMMENDATION
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("RECOMMENDATION\n")
cat(rep("=", 70), "\n\n")

if (!exists("overall_range")) {
  cat("⚠ Cannot make recommendation - missing exam_period data.\n\n")
} else if (overall_range > 1.0) {
  cat("✓ PROCEED with between-within analysis.\n")
  cat("  Your data shows sufficient within-person variance.\n")
  cat("  Run: source('11_prepare_between_within_data.R')\n\n")
} else if (overall_range > 0.5) {
  cat("⚠ PROCEED WITH CAUTION.\n")
  cat("  Within-person variance is modest.\n")
  cat("  Results may show limited state effects.\n")
  cat("  Worth trying, but trait-only model may be sufficient.\n\n")
} else {
  cat("✗ RECONSIDER between-within analysis.\n")
  cat("  Within-person variance is very low (traits very stable).\n")
  cat("  Trait-only model (already done) is likely sufficient.\n")
  cat("  Consider reporting high trait stability as finding.\n\n")
}

cat("Next step: Review output above and decide whether to proceed.\n\n")
