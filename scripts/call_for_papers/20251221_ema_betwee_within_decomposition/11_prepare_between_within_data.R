# ==============================================================================
# 11_prepare_between_within_data.R
# Prepare data with between-person and within-person components
# ==============================================================================
# PURPOSE:
#   Decompose EMA PID-5 scores into:
#   - Between-person (trait): Person's average across all timepoints
#   - Within-person (state): Deviation from person mean at each timepoint
#
# OUTPUT:
#   Dataset with both components for testing:
#   Voice[t] = Trait_between + State_within + Timepoint + Interactions
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(rio)
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("BETWEEN-WITHIN DECOMPOSITION: EMA PID-5\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 1. LOAD EMA DATA (with timepoint-specific PID-5 scores)
# ==============================================================================

cat("Loading EMA data...\n")

# This should be the original EMA file with PID-5 at each assessment
data_path <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
if (!file.exists(data_path)) {
  stop(
    "File not found: ",
    data_path,
    "\nRun 01_clean_after_merge_FINAL.R first."
  )
}

d_ema <- rio::import(data_path)

# Rename for consistency
d_ema <- d_ema %>%
  rename(
    ID = any_of(c("user_id", "UserID", "id_anon"))
  )

# PID-5 EMA variables
pid5_vars <- c(
  "pid5_negative_affectivity",
  "pid5_detachment",
  "pid5_antagonism",
  "pid5_disinhibition",
  "pid5_psychoticism"
)

# Check variables exist
missing_vars <- setdiff(pid5_vars, names(d_ema))
if (length(missing_vars) > 0) {
  stop("Missing PID-5 variables: ", paste(missing_vars, collapse = ", "))
}

cat(
  "✓ EMA data loaded: N =",
  nrow(d_ema),
  "observations,",
  n_distinct(d_ema$ID),
  "subjects\n"
)

# ==============================================================================
# 2. IDENTIFY VOICE RECORDING TIMEPOINTS
# ==============================================================================

cat("\nIdentifying voice recording timepoints...\n")

# Load voice data to get IDs and dates
voice_baseline <- readxl::read_excel(
  here::here("data", "raw", "acustic_features", "datiacustici", "AUDIO.xlsx"),
  sheet = "BASELINE"
)
voice_pre <- readxl::read_excel(
  here::here("data", "raw", "acustic_features", "datiacustici", "AUDIO.xlsx"),
  sheet = "PRE"
)
voice_post <- readxl::read_excel(
  here::here("data", "raw", "acustic_features", "datiacustici", "AUDIO.xlsx"),
  sheet = "POST"
)

# Correct ID names (same corrections as in script 02)
correct_ids <- function(df) {
  df %>%
    mutate(
      ID = case_when(
        ID == "am_bo_1988_08_24_166" ~ "an_bo_1988_08_24_166",
        ID == "as_li_2005_04_26_447" ~ "as_si_2005_04_26_447",
        ID == "cl_bo_1987_10_16_628" ~ "ca_bo_1987_10_16_628",
        ID == "hi_na_2005_03_08_339" ~ "gi_na_2005_03_08_339",
        ID == "ma_si_2003_10_31_940" ~ "si_ma_2003_10_31_940",
        TRUE ~ ID
      )
    )
}

# Function to correct date years
correct_date_years <- function(date_vector) {
  # Convert to Date if not already
  dates <- as.Date(date_vector)

  # Get year
  years <- lubridate::year(dates)

  # Check for problematic years (< 2020 or > 2026)
  problematic <- years < 2020 | years > 2026

  if (any(problematic, na.rm = TRUE)) {
    cat(sprintf(
      "  ⚠ Found %d dates with incorrect years\n",
      sum(problematic, na.rm = TRUE)
    ))

    # For dates with year < 2020, assume it should be 2024 or 2025
    # We'll change year to 2025 (current study year)
    dates[problematic] <- dates[problematic] +
      lubridate::years(2025 - years[problematic])

    cat("  ✓ Corrected years to 2025\n")
  }

  return(dates)
}

voice_baseline <- correct_ids(voice_baseline) %>%
  select(ID, Data) %>%
  mutate(
    timepoint = "baseline",
    date = correct_date_years(Data)
  )

voice_pre <- correct_ids(voice_pre) %>%
  select(ID, Data) %>%
  mutate(
    timepoint = "pre_exam",
    date = correct_date_years(Data)
  )

voice_post <- correct_ids(voice_post) %>%
  select(ID, Data) %>%
  mutate(
    timepoint = "post_exam",
    date = correct_date_years(Data)
  )

voice_dates <- bind_rows(voice_baseline, voice_pre, voice_post) %>%
  select(ID, timepoint, voice_date = date)

# Verify date ranges
date_range <- range(voice_dates$voice_date, na.rm = TRUE)
cat(
  "✓ Voice recording dates identified for",
  n_distinct(voice_dates$ID),
  "subjects\n"
)
cat(sprintf(
  "  Date range: %s to %s\n",
  format(date_range[1], "%Y-%m-%d"),
  format(date_range[2], "%Y-%m-%d")
))

# Check if any dates still look suspicious
suspicious_dates <- voice_dates %>%
  mutate(year = lubridate::year(voice_date)) %>%
  filter(year < 2024 | year > 2026)

if (nrow(suspicious_dates) > 0) {
  cat("\n⚠ WARNING: Still found suspicious dates:\n")
  print(suspicious_dates)
  cat("\nPlease verify these dates manually\n\n")
} else {
  cat("  ✓ All dates are in expected range (2024-2026)\n\n")
}

# ==============================================================================
# 3. MATCH EMA PID-5 TO VOICE RECORDING DATES
# ==============================================================================

cat("\nMatching EMA PID-5 to voice recording dates...\n")

# Prepare EMA data with dates
d_ema_dated <- d_ema %>%
  mutate(
    ema_date = as.Date(day),
    exam_period = tolower(as.character(exam_period))
  ) %>%
  filter(
    exam_period %in% c("baseline", "pre_exam", "post_exam"),
    !is.na(pid5_negative_affectivity) # At least one PID-5 value must be non-NA
  ) %>%
  select(ID, ema_date, exam_period, all_of(pid5_vars))

cat("  Total EMA observations available:", nrow(d_ema_dated), "\n")
cat("  Distribution by exam_period:\n")
print(table(d_ema_dated$exam_period))

# Convert timepoint to exam_period for matching
voice_dates_clean <- voice_dates %>%
  mutate(
    exam_period = case_when(
      timepoint == "baseline" ~ "baseline",
      timepoint == "pre_exam" ~ "pre_exam",
      timepoint == "post_exam" ~ "post_exam"
    )
  )

cat("\n✓ Voice recordings to match:", nrow(voice_dates_clean), "\n")
cat("  Expected: 141 × 3 = 423 observations\n")
cat("  Distribution by timepoint:\n")
print(table(voice_dates_clean$timepoint))

# Function to find closest EMA date for each voice recording
find_closest_ema <- function(voice_id, voice_period, voice_date, ema_data) {
  # Filter EMA for same ID and period
  ema_subset <- ema_data %>%
    filter(ID == voice_id, exam_period == voice_period)

  if (nrow(ema_subset) == 0) {
    # No EMA data for this ID/period
    return(list(
      matched = FALSE,
      date_diff = NA_real_,
      ema_date = as.Date(NA),
      pid5_data = rep(NA_real_, length(pid5_vars))
    ))
  }

  # Calculate date differences
  ema_subset <- ema_subset %>%
    mutate(
      date_diff = abs(as.numeric(difftime(
        ema_date,
        voice_date,
        units = "days"
      )))
    )

  # Find closest
  closest <- ema_subset %>%
    filter(date_diff == min(date_diff)) %>%
    slice(1) # If tie, take first

  # Extract PID-5 values
  pid5_values <- closest %>%
    select(all_of(pid5_vars)) %>%
    as.list() %>%
    unlist()

  return(list(
    matched = TRUE,
    date_diff = closest$date_diff,
    ema_date = closest$ema_date,
    pid5_data = pid5_values
  ))
}

# Match each voice recording to closest EMA
cat("\nMatching voice recordings to EMA assessments...\n")

matched_data <- voice_dates_clean %>%
  rowwise() %>%
  mutate(
    match_result = list(find_closest_ema(
      ID,
      exam_period,
      voice_date,
      d_ema_dated
    ))
  ) %>%
  ungroup() %>%
  mutate(
    matched = map_lgl(match_result, ~ .x$matched),
    date_diff_days = map_dbl(match_result, ~ .x$date_diff),
    ema_date = map(match_result, ~ .x$ema_date) %>%
      map(as.Date, origin = "1970-01-01") %>%
      unlist() %>%
      as.Date(origin = "1970-01-01")
  )

# Extract PID-5 values
for (i in seq_along(pid5_vars)) {
  var_name <- pid5_vars[i]
  matched_data[[var_name]] <- map_dbl(
    matched_data$match_result,
    ~ .x$pid5_data[i]
  )
}

# Remove temporary column
matched_data <- matched_data %>%
  select(-match_result)

# Report matching statistics
cat("\n", rep("=", 70), "\n", sep = "")
cat("MATCHING RESULTS\n")
cat(rep("=", 70), "\n")

n_total <- nrow(matched_data)
n_matched <- sum(matched_data$matched)
n_unmatched <- n_total - n_matched

cat(sprintf("Total voice recordings: %d\n", n_total))
cat(sprintf(
  "Successfully matched: %d (%.1f%%)\n",
  n_matched,
  100 * n_matched / n_total
))
cat(sprintf(
  "No EMA data available: %d (%.1f%%)\n",
  n_unmatched,
  100 * n_unmatched / n_total
))

# Among matched, report date differences
matched_subset <- matched_data %>% filter(matched)

if (n_matched > 0) {
  cat("\nDate matching quality (among matched observations):\n")

  n_perfect <- sum(matched_subset$date_diff_days == 0, na.rm = TRUE)
  n_within_1 <- sum(matched_subset$date_diff_days <= 1, na.rm = TRUE)
  n_within_3 <- sum(matched_subset$date_diff_days <= 3, na.rm = TRUE)
  n_within_7 <- sum(matched_subset$date_diff_days <= 7, na.rm = TRUE)

  cat(sprintf(
    "  Perfect match (same day): %d (%.1f%%)\n",
    n_perfect,
    100 * n_perfect / n_matched
  ))
  cat(sprintf(
    "  Within 1 day: %d (%.1f%%)\n",
    n_within_1,
    100 * n_within_1 / n_matched
  ))
  cat(sprintf(
    "  Within 3 days: %d (%.1f%%)\n",
    n_within_3,
    100 * n_within_3 / n_matched
  ))
  cat(sprintf(
    "  Within 7 days: %d (%.1f%%)\n",
    n_within_7,
    100 * n_within_7 / n_matched
  ))
  cat(sprintf(
    "  Mean difference: %.1f days (SD = %.1f)\n",
    mean(matched_subset$date_diff_days, na.rm = TRUE),
    sd(matched_subset$date_diff_days, na.rm = TRUE)
  ))
  cat(sprintf(
    "  Median difference: %.1f days\n",
    median(matched_subset$date_diff_days, na.rm = TRUE)
  ))
  cat(sprintf(
    "  Max difference: %.0f days\n",
    max(matched_subset$date_diff_days, na.rm = TRUE)
  ))
}

# Report by timepoint
cat("\nMatching by timepoint:\n")
matching_by_tp <- matched_data %>%
  group_by(timepoint) %>%
  summarise(
    n_total = n(),
    n_matched = sum(matched),
    pct_matched = 100 * n_matched / n_total,
    mean_date_diff = mean(date_diff_days, na.rm = TRUE),
    n_perfect_match = sum(date_diff_days == 0, na.rm = TRUE),
    .groups = "drop"
  )

print(matching_by_tp)

cat("\n")

# ==============================================================================
# 4. QUALITY CHECK: Temporal distance of matches
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("QUALITY CHECK: Temporal distance of matches\n")
cat(rep("=", 70), "\n\n")

# Define maximum acceptable temporal distance (in days)
max_diff_baseline <- 21 # 3 weeks for baseline
max_diff_pre_post <- 7 # 1 week for pre/post exam

# Check baseline matches
baseline_matches <- matched_data %>%
  filter(matched, timepoint == "baseline")

if (nrow(baseline_matches) > 0) {
  cat("BASELINE matches:\n")
  cat(sprintf("  Total matched: %d\n", nrow(baseline_matches)))

  # Distribution of date differences
  baseline_diffs <- baseline_matches$date_diff_days
  cat(sprintf(
    "  Min difference: %.0f days\n",
    min(baseline_diffs, na.rm = TRUE)
  ))
  cat(sprintf(
    "  Q1: %.0f days\n",
    quantile(baseline_diffs, 0.25, na.rm = TRUE)
  ))
  cat(sprintf("  Median: %.0f days\n", median(baseline_diffs, na.rm = TRUE)))
  cat(sprintf(
    "  Q3: %.0f days\n",
    quantile(baseline_diffs, 0.75, na.rm = TRUE)
  ))
  cat(sprintf(
    "  Max difference: %.0f days\n",
    max(baseline_diffs, na.rm = TRUE)
  ))
  cat(sprintf(
    "  Mean: %.1f days (SD = %.1f)\n\n",
    mean(baseline_diffs, na.rm = TRUE),
    sd(baseline_diffs, na.rm = TRUE)
  ))

  # Count problematic matches
  n_problematic <- sum(baseline_diffs > max_diff_baseline, na.rm = TRUE)

  if (n_problematic > 0) {
    cat(sprintf(
      "⚠ WARNING: %d baseline matches exceed %d days threshold:\n\n",
      n_problematic,
      max_diff_baseline
    ))

    problematic <- baseline_matches %>%
      filter(date_diff_days > max_diff_baseline) %>%
      arrange(desc(date_diff_days)) %>%
      select(ID, voice_date, ema_date, date_diff_days)

    print(problematic, n = min(20, nrow(problematic)))

    if (nrow(problematic) > 20) {
      cat(sprintf("\n... and %d more\n", nrow(problematic) - 20))
    }
    cat("\n")
  }
}

# Check pre/post matches
pre_post_matches <- matched_data %>%
  filter(matched, timepoint %in% c("pre_exam", "post_exam"))

if (nrow(pre_post_matches) > 0) {
  n_problematic_pp <- sum(
    pre_post_matches$date_diff_days > max_diff_pre_post,
    na.rm = TRUE
  )

  if (n_problematic_pp > 0) {
    cat(sprintf(
      "⚠ WARNING: %d pre/post matches exceed %d days threshold:\n\n",
      n_problematic_pp,
      max_diff_pre_post
    ))

    problematic_pp <- pre_post_matches %>%
      filter(date_diff_days > max_diff_pre_post) %>%
      arrange(desc(date_diff_days)) %>%
      select(ID, timepoint, voice_date, ema_date, date_diff_days)

    print(problematic_pp, n = min(20, nrow(problematic_pp)))
    cat("\n")
  }
}

# ==============================================================================
# 5. APPLY TEMPORAL THRESHOLD
# ==============================================================================

cat(rep("=", 70), "\n")
cat("APPLYING TEMPORAL THRESHOLDS\n")
cat(rep("=", 70), "\n\n")

cat(sprintf(
  "Excluding matches with temporal distance > %d days (baseline)\n",
  max_diff_baseline
))
cat(sprintf(
  "Excluding matches with temporal distance > %d days (pre/post)\n\n",
  max_diff_pre_post
))

# Flag valid matches
matched_data_filtered <- matched_data %>%
  mutate(
    valid_match = case_when(
      !matched ~ FALSE,
      timepoint == "baseline" ~ date_diff_days <= max_diff_baseline,
      timepoint %in% c("pre_exam", "post_exam") ~
        date_diff_days <= max_diff_pre_post,
      TRUE ~ FALSE
    )
  )

# Report exclusions
n_excluded_baseline <- sum(
  matched_data_filtered$timepoint == "baseline" &
    matched_data_filtered$matched &
    !matched_data_filtered$valid_match,
  na.rm = TRUE
)
n_excluded_pre <- sum(
  matched_data_filtered$timepoint == "pre_exam" &
    matched_data_filtered$matched &
    !matched_data_filtered$valid_match,
  na.rm = TRUE
)
n_excluded_post <- sum(
  matched_data_filtered$timepoint == "post_exam" &
    matched_data_filtered$matched &
    !matched_data_filtered$valid_match,
  na.rm = TRUE
)

cat("Excluded due to temporal threshold:\n")
cat(sprintf("  Baseline: %d observations\n", n_excluded_baseline))
cat(sprintf("  Pre-exam: %d observations\n", n_excluded_pre))
cat(sprintf("  Post-exam: %d observations\n\n", n_excluded_post))

# Summary after filtering
summary_filtered <- matched_data_filtered %>%
  group_by(timepoint) %>%
  summarise(
    n_total = n(),
    n_valid = sum(valid_match),
    pct_valid = 100 * n_valid / n_total,
    mean_date_diff = mean(date_diff_days[valid_match], na.rm = TRUE),
    median_date_diff = median(date_diff_days[valid_match], na.rm = TRUE),
    max_date_diff = max(date_diff_days[valid_match], na.rm = TRUE),
    .groups = "drop"
  )

cat("Valid matches after temporal filtering:\n")
print(summary_filtered)
cat("\n")

# Prepare dataset with valid matches only
ema_matched <- matched_data_filtered %>%
  filter(valid_match) %>%
  select(ID, timepoint, voice_date, ema_date, date_diff_days, all_of(pid5_vars))

cat("✓ Matched dataset (after temporal filtering):\n")
cat("  N observations:", nrow(ema_matched), "\n")
cat("  N unique subjects:", n_distinct(ema_matched$ID), "\n")
cat("\n")

# ==============================================================================
# 6. CHECK TIMEPOINT COVERAGE
# ==============================================================================

cat(rep("=", 70), "\n")
cat("TIMEPOINT COVERAGE CHECK\n")
cat(rep("=", 70), "\n\n")

# Check coverage per subject
coverage_check <- ema_matched %>%
  group_by(ID) %>%
  summarise(
    n_timepoints = n(),
    has_baseline = "baseline" %in% timepoint,
    has_pre = "pre_exam" %in% timepoint,
    has_post = "post_exam" %in% timepoint,
    .groups = "drop"
  )

n_with_3 <- sum(coverage_check$n_timepoints == 3)
n_with_2 <- sum(coverage_check$n_timepoints == 2)
n_with_1 <- sum(coverage_check$n_timepoints == 1)

cat(sprintf(
  "Subjects with 3 timepoints: %d (%.1f%%)\n",
  n_with_3,
  100 * n_with_3 / nrow(coverage_check)
))
cat(sprintf(
  "Subjects with 2 timepoints: %d (%.1f%%)\n",
  n_with_2,
  100 * n_with_2 / nrow(coverage_check)
))
cat(sprintf(
  "Subjects with 1 timepoint: %d (%.1f%%) ← Cannot compute within-person variance\n",
  n_with_1,
  100 * n_with_1 / nrow(coverage_check)
))

# Show missing pattern for subjects with 2 timepoints
if (n_with_2 > 0) {
  cat("\nMissing pattern for subjects with 2 timepoints:\n")
  missing_pattern <- coverage_check %>%
    filter(n_timepoints == 2) %>%
    summarise(
      missing_baseline = sum(!has_baseline),
      missing_pre = sum(!has_pre),
      missing_post = sum(!has_post)
    )

  cat(sprintf("  Missing baseline: %d\n", missing_pattern$missing_baseline))
  cat(sprintf("  Missing pre-exam: %d\n", missing_pattern$missing_pre))
  cat(sprintf("  Missing post-exam: %d\n\n", missing_pattern$missing_post))
}

# DECISION: Keep only subjects with 2+ timepoints for between-within decomposition
if (n_with_1 > 0) {
  cat("⚠ DECISION: Excluding", n_with_1, "subjects with only 1 timepoint\n")
  cat("  (within-person decomposition requires 2+ observations)\n\n")

  ids_to_keep <- coverage_check %>%
    filter(n_timepoints >= 2) %>%
    pull(ID)

  ema_timepoints <- ema_matched %>%
    filter(ID %in% ids_to_keep)

  cat("✓ After exclusion:\n")
  cat("  N subjects:", n_distinct(ema_timepoints$ID), "\n")
  cat("  N observations:", nrow(ema_timepoints), "\n\n")
} else {
  ema_timepoints <- ema_matched
  cat("✓ All subjects have 2+ timepoints\n\n")
}

cat("✓ Dataset ready for between-within decomposition\n")
cat(
  "  Mean timepoints per subject:",
  round(nrow(ema_timepoints) / n_distinct(ema_timepoints$ID), 2),
  "\n\n"
)

# ==============================================================================
# 7. BETWEEN-WITHIN DECOMPOSITION
# ==============================================================================

cat(rep("=", 70), "\n")
cat("BETWEEN-WITHIN DECOMPOSITION\n")
cat(rep("=", 70), "\n\n")

# Compute person means (between-person component)
# This is the average across the available timepoints for each person
person_means <- ema_timepoints %>%
  group_by(ID) %>%
  summarise(
    across(
      all_of(pid5_vars),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_between"
    ),
    n_timepoints = n(),
    .groups = "drop"
  )

cat("Between-person components computed for", nrow(person_means), "subjects\n")
cat(
  "  Mean timepoints per person:",
  round(mean(person_means$n_timepoints), 1),
  "\n\n"
)

# Add person means to data and compute within-person deviations
ema_decomposed <- ema_timepoints %>%
  left_join(person_means, by = "ID") %>%
  mutate(
    # Within-person deviations (state component)
    # This is how much the person deviates from their own average at each timepoint
    across(
      all_of(pid5_vars),
      ~ .x - get(paste0(cur_column(), "_between")),
      .names = "{.col}_within"
    )
  )

# Verify decomposition
cat("Verifying decomposition:\n")
cat("(Within-person deviations should average to ~0 for each person)\n\n")

for (var in pid5_vars) {
  var_between <- paste0(var, "_between")
  var_within <- paste0(var, "_within")

  # Within-person deviations should average to ~0 for each person
  check <- ema_decomposed %>%
    group_by(ID) %>%
    summarise(
      mean_within = mean(get(var_within), na.rm = TRUE),
      .groups = "drop"
    )

  max_deviation <- max(abs(check$mean_within), na.rm = TRUE)

  var_name_short <- str_replace(var, "pid5_", "")

  if (max_deviation < 0.001) {
    cat(sprintf(
      "  ✓ %s: Decomposition correct (max dev = %.6f)\n",
      var_name_short,
      max_deviation
    ))
  } else {
    cat(sprintf(
      "  ⚠ %s: Large deviation (%.6f)\n",
      var_name_short,
      max_deviation
    ))
  }
}

cat("\n")

# ==============================================================================
# 8. ADD VOICE DATA
# ==============================================================================

cat("Merging with voice data...\n")

# Load voice data (from script 02)
if (!file.exists("results/df_analysis.rds")) {
  stop(
    "Voice data not found: results/df_analysis.rds\n",
    "Please run 02_voice_personality_analysis_FINAL.R first."
  )
}

df_voice <- readRDS("results/df_analysis.rds") %>%
  select(
    ID,
    timepoint,
    date,
    starts_with("f0_"),
    starts_with("f2_"),
    starts_with("jitter_"),
    starts_with("nne_")
  )

cat("✓ Voice data loaded from results/df_analysis.rds\n")
cat("  N observations:", nrow(df_voice), "\n")
cat("  N subjects:", n_distinct(df_voice$ID), "\n\n")

# CRITICAL: Standardize timepoint names
cat("Standardizing timepoint names...\n")

# Convert to character first
df_voice <- df_voice %>%
  mutate(timepoint = as.character(timepoint))

ema_decomposed <- ema_decomposed %>%
  mutate(timepoint = as.character(timepoint))

# Standardize voice timepoint names to match EMA
df_voice <- df_voice %>%
  mutate(
    timepoint = case_when(
      timepoint == "baseline" ~ "baseline",
      timepoint == "pre" ~ "pre_exam",
      timepoint == "post" ~ "post_exam",
      TRUE ~ timepoint
    )
  )

cat("✓ Timepoint names standardized\n")
cat("Voice data timepoints:\n")
print(table(df_voice$timepoint))
cat("\nEMA data timepoints:\n")
print(table(ema_decomposed$timepoint))
cat("\n")

# Check overlap
common_ids <- intersect(ema_decomposed$ID, df_voice$ID)
cat("Common IDs between EMA and voice:", length(common_ids), "\n")

# Count potential matches
matches <- ema_decomposed %>%
  semi_join(df_voice, by = c("ID", "timepoint"))

cat("Potential matches:", nrow(matches), "out of", nrow(ema_decomposed), "\n\n")

# CRITICAL: Start from ema_decomposed (matched + filtered EMA)
# and add voice features
df_between_within <- ema_decomposed %>%
  left_join(
    df_voice,
    by = c("ID", "timepoint")
  ) %>%
  # Also keep the voice_date for reference
  rename(date_voice = date) %>%
  select(
    ID,
    timepoint,
    voice_date,
    ema_date,
    date_diff_days,
    date_voice,
    all_of(pid5_vars),
    ends_with("_between"),
    ends_with("_within"),
    starts_with("f0_"),
    starts_with("f2_"),
    starts_with("jitter_"),
    starts_with("nne_")
  )

cat("✓ Merged dataset: N =", nrow(df_between_within), "observations\n")
cat("  N subjects:", n_distinct(df_between_within$ID), "\n\n")

# Check for missing voice data
n_missing_voice_a <- sum(is.na(df_between_within$f0_mean_a))
n_missing_voice_i <- sum(is.na(df_between_within$f0_mean_i))
n_missing_voice_u <- sum(is.na(df_between_within$f0_mean_u))

if (n_missing_voice_a > 0 | n_missing_voice_i > 0 | n_missing_voice_u > 0) {
  cat("⚠ WARNING: Missing voice data:\n")
  cat(sprintf("  Vowel /a/: %d observations\n", n_missing_voice_a))
  cat(sprintf("  Vowel /i/: %d observations\n", n_missing_voice_i))
  cat(sprintf("  Vowel /u/: %d observations\n", n_missing_voice_u))
  cat(
    "  These observations will be excluded from respective voice analyses\n\n"
  )
} else {
  cat("✓ No missing voice data\n\n")
}

# ==============================================================================
# 9. CENTER PREDICTORS
# ==============================================================================

cat("Centering predictors...\n")

# Grand-mean center between-person components
for (var in pid5_vars) {
  var_between <- paste0(var, "_between")
  var_between_c <- paste0(var, "_between_c")

  df_between_within[[var_between_c]] <- scale(
    df_between_within[[var_between]],
    center = TRUE,
    scale = FALSE
  )[, 1]
}

# Within-person components are already deviations (person-mean centered)
# Just rename for clarity
for (var in pid5_vars) {
  var_within <- paste0(var, "_within")
  var_within_c <- paste0(var, "_within_c")

  df_between_within[[var_within_c]] <- df_between_within[[var_within]]
}

cat("✓ Predictors centered\n\n")

# ==============================================================================
# 10. ADD CONTRAST CODES
# ==============================================================================

cat("Adding contrast codes...\n")

df_between_within <- df_between_within %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("baseline", "pre_exam", "post_exam")
    ),
    c1_stress = case_when(
      timepoint == "baseline" ~ -0.5,
      timepoint == "pre_exam" ~ 0.5,
      timepoint == "post_exam" ~ 0
    ),
    c2_recovery = case_when(
      timepoint == "baseline" ~ 0,
      timepoint == "pre_exam" ~ -0.5,
      timepoint == "post_exam" ~ 0.5
    )
  )

cat("✓ Contrast codes added\n\n")

# ==============================================================================
# 11. SAVE OUTPUT
# ==============================================================================

cat("Saving output...\n")

dir.create("results/between_within", showWarnings = FALSE, recursive = TRUE)

saveRDS(df_between_within, "results/between_within/df_between_within.rds")
rio::export(df_between_within, "results/between_within/df_between_within.csv")

cat("✓ Saved: results/between_within/df_between_within.rds\n")
cat("✓ Saved: results/between_within/df_between_within.csv\n\n")

# ==============================================================================
# 12. SUMMARY STATISTICS
# ==============================================================================

cat(rep("=", 70), "\n")
cat("SUMMARY STATISTICS\n")
cat(rep("=", 70), "\n\n")

# Dataset overview
cat("Final dataset:\n")
cat("  N observations:", nrow(df_between_within), "\n")
cat("  N subjects:", n_distinct(df_between_within$ID), "\n")
cat(
  "  Observations per subject:",
  round(nrow(df_between_within) / n_distinct(df_between_within$ID), 1),
  "\n\n"
)

# Timepoint distribution
cat("Observations by timepoint:\n")
print(table(df_between_within$timepoint))
cat("\n")

# Between-person variance
cat("Between-person variance (trait):\n")
for (var in pid5_vars) {
  var_between <- paste0(var, "_between")
  m <- mean(df_between_within[[var_between]], na.rm = TRUE)
  sd <- sd(df_between_within[[var_between]], na.rm = TRUE)
  var_name_short <- str_replace(var, "pid5_", "")
  cat(sprintf("  %s: M = %.3f, SD = %.3f\n", var_name_short, m, sd))
}

cat("\nWithin-person variance (state):\n")
for (var in pid5_vars) {
  var_within <- paste0(var, "_within")

  # Compute SD of within-person deviations
  # This tells us how much individuals fluctuate around their own mean
  within_sd <- sd(df_between_within[[var_within]], na.rm = TRUE)

  var_name_short <- str_replace(var, "pid5_", "")
  cat(sprintf(
    "  %s: SD of within-person deviations = %.3f\n",
    var_name_short,
    within_sd
  ))
}

# ICC: Proportion of variance that is between-person
cat("\nIntraclass Correlation (ICC) - Proportion between-person:\n")
for (var in pid5_vars) {
  var_between <- paste0(var, "_between")
  var_within <- paste0(var, "_within")

  # Variance of person means (between-person variance)
  var_between_total <- var(df_between_within[[var_between]], na.rm = TRUE)

  # Variance of within-person deviations (within-person variance)
  var_within_total <- var(df_between_within[[var_within]], na.rm = TRUE)

  # ICC = between-variance / (between-variance + within-variance)
  icc <- var_between_total / (var_between_total + var_within_total)

  var_name_short <- str_replace(var, "pid5_", "")
  cat(sprintf("  %s: ICC = %.3f ", var_name_short, icc))

  # Interpretation
  if (icc > 0.90) {
    cat("(trait-dominant)\n")
  } else if (icc > 0.70) {
    cat("(mostly trait)\n")
  } else if (icc > 0.50) {
    cat("(balanced)\n")
  } else if (icc > 0.30) {
    cat("(mostly state)\n")
  } else {
    cat("(state-dominant)\n")
  }
}

cat("\nInterpretation:\n")
cat("- ICC close to 1: Mostly between-person differences (stable trait)\n")
cat("- ICC close to 0: Mostly within-person fluctuation (dynamic state)\n")
cat("- ICC ~0.5: Equal contribution of trait and state\n\n")

cat(rep("=", 70), "\n")
cat("DATA PREPARATION COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("Next steps:\n")
cat("1. Run 12_fit_between_within_models.R to estimate models\n")
cat("2. Run 13_compare_between_within.R to compare with trait-only models\n\n")
