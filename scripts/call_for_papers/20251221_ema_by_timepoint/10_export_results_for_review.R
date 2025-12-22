# ==============================================================================
# 10_export_results_SAFE.R
# Export results with automatic structure detection
# ==============================================================================

library(tidyverse)
library(rio)

cat("\n", rep("=", 70), "\n", sep = "")
cat("INSPECTING AND EXPORTING RESULTS\n")
cat(rep("=", 70), "\n\n")

# Create output directory
dir.create("results/exports", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# HELPER FUNCTION: Safe export
# ==============================================================================

safe_export <- function(rds_file, output_name, description) {
  cat("\n", rep("-", 70), "\n", sep = "")
  cat("Processing:", description, "\n")
  cat("Source:", rds_file, "\n")

  if (!file.exists(rds_file)) {
    cat("  ⚠ File not found, skipping...\n")
    return(NULL)
  }

  # Read file
  obj <- readRDS(rds_file)

  # Inspect structure
  cat("\nStructure:\n")
  cat("  Class:", class(obj), "\n")

  if (is.data.frame(obj)) {
    cat("  Dimensions:", nrow(obj), "rows ×", ncol(obj), "columns\n")
    cat(
      "  Column names:",
      paste(names(obj)[1:min(10, ncol(obj))], collapse = ", ")
    )
    if (ncol(obj) > 10) cat(", ... (", ncol(obj) - 10, " more)", sep = "")
    cat("\n")

    # Export as CSV
    csv_file <- paste0("results/exports/", output_name, ".csv")
    write_csv(obj, csv_file)
    cat("  ✓ Saved CSV:", csv_file, "\n")

    # Export summary as TXT
    txt_file <- paste0("results/exports/", output_name, "_SUMMARY.txt")
    sink(txt_file)
    cat(rep("=", 70), "\n", sep = "")
    cat(toupper(description), "\n")
    cat(rep("=", 70), "\n\n", sep = "")
    cat("Dimensions:", nrow(obj), "rows ×", ncol(obj), "columns\n\n")
    cat("Column names:\n")
    print(names(obj))
    cat("\n\nFirst 20 rows:\n")
    cat(rep("-", 70), "\n", sep = "")
    print(head(obj, 20))
    cat("\n\nSummary statistics:\n")
    cat(rep("-", 70), "\n", sep = "")
    print(summary(obj))
    sink()
    cat("  ✓ Saved summary:", txt_file, "\n")
  } else if (is.list(obj)) {
    cat("  List with", length(obj), "elements\n")
    cat("  Element names:", paste(names(obj), collapse = ", "), "\n")

    # Export each element
    for (i in seq_along(obj)) {
      elem_name <- names(obj)[i]
      if (is.null(elem_name) || elem_name == "")
        elem_name <- paste0("element_", i)

      elem <- obj[[i]]

      if (is.data.frame(elem)) {
        csv_file <- paste0(
          "results/exports/",
          output_name,
          "_",
          elem_name,
          ".csv"
        )
        write_csv(elem, csv_file)
        cat("  ✓ Saved element '", elem_name, "': ", csv_file, "\n", sep = "")
      }
    }

    # Export full structure as TXT
    txt_file <- paste0("results/exports/", output_name, "_STRUCTURE.txt")
    sink(txt_file)
    cat(rep("=", 70), "\n", sep = "")
    cat(toupper(description), "\n")
    cat(rep("=", 70), "\n\n", sep = "")
    cat("List structure:\n\n")
    print(str(obj, max.level = 2))
    sink()
    cat("  ✓ Saved structure:", txt_file, "\n")
  } else {
    cat("  ⚠ Unrecognized object type, saving as RDS only\n")
  }

  cat("  → Done\n")
  return(obj)
}

# ==============================================================================
# PROCESS ALL FILES
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PROCESSING FILES\n")
cat(rep("=", 70), "\n")

# 1. EMA Validation
validation <- safe_export(
  "results/ema_validation_correlations.rds",
  "01_ema_validation",
  "EMA Validation Correlations"
)

# 2. Moderation Results
moderation <- safe_export(
  "results/results_table.rds",
  "02_moderation_results",
  "Moderation Analysis Results"
)

# 3. Temporal Covariation
temporal <- safe_export(
  "results/temporal_covariation_results.rds",
  "03_temporal_covariation",
  "Temporal Covariation Results"
)

# 4. EMA vs Baseline Comparison
comparison <- safe_export(
  "results/ema_vs_baseline_comparison.rds",
  "04_ema_vs_baseline",
  "EMA vs Baseline Comparison"
)

# 5. Dataset
dataset <- safe_export(
  "results/df_analysis.rds",
  "05_dataset_analysis",
  "Analysis Dataset"
)

# ==============================================================================
# CREATE MASTER SUMMARY
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("CREATING MASTER SUMMARY\n")
cat(rep("=", 70), "\n\n")

sink("results/exports/00_MASTER_SUMMARY.txt")

cat(rep("=", 70), "\n", sep = "")
cat("VOICE & PERSONALITY ANALYSIS - MASTER SUMMARY\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat(rep("=", 70), "\n", sep = "")
cat("FILES EXPORTED\n")
cat(rep("=", 70), "\n\n", sep = "")

exported_files <- list.files("results/exports", pattern = "\\.csv$|\\.txt$")
for (f in sort(exported_files)) {
  cat("  -", f, "\n")
}

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("QUICK RESULTS OVERVIEW\n")
cat(rep("=", 70), "\n\n", sep = "")

# EMA Validation
if (!is.null(validation)) {
  cat("1. EMA VALIDATION\n")
  cat(rep("-", 70), "\n", sep = "")
  if ("r" %in% names(validation)) {
    cat(
      "  Mean correlation:",
      sprintf("%.3f", mean(validation$r, na.rm = TRUE)),
      "\n"
    )
    cat(
      "  Range:",
      sprintf(
        "%.3f - %.3f",
        min(validation$r, na.rm = TRUE),
        max(validation$r, na.rm = TRUE)
      ),
      "\n"
    )
    strong <- sum(validation$r > 0.7, na.rm = TRUE)
    moderate <- sum(validation$r >= 0.5 & validation$r <= 0.7, na.rm = TRUE)
    weak <- sum(validation$r < 0.5, na.rm = TRUE)
    cat("  Strong (r > 0.70):", strong, "domains\n")
    cat("  Moderate (r 0.50-0.70):", moderate, "domains\n")
    cat("  Weak (r < 0.50):", weak, "domains\n\n")
  }
}

# Moderation Results
if (!is.null(moderation)) {
  cat("2. MODERATION ANALYSIS\n")
  cat(rep("-", 70), "\n", sep = "")
  cat("  Total effects tested:", nrow(moderation), "\n")

  # Try to identify significant effects
  if (all(c("ci_lower", "ci_upper") %in% names(moderation))) {
    sig <- sum(
      (moderation$ci_lower > 0) | (moderation$ci_upper < 0),
      na.rm = TRUE
    )
    cat("  Significant effects (95% CI excludes 0):", sig, "\n\n")
  } else if (all(c("Q2.5", "Q97.5") %in% names(moderation))) {
    sig <- sum((moderation$Q2.5 > 0) | (moderation$Q97.5 < 0), na.rm = TRUE)
    cat("  Significant effects (95% CI excludes 0):", sig, "\n\n")
  }
}

# Temporal Covariation
if (!is.null(temporal)) {
  cat("3. TEMPORAL COVARIATION\n")
  cat(rep("-", 70), "\n", sep = "")
  if ("sig_within" %in% names(temporal)) {
    sig_within <- sum(temporal$sig_within, na.rm = TRUE)
    sig_between <- sum(temporal$sig_between, na.rm = TRUE)
    cat("  Significant within-person effects:", sig_within, "\n")
    cat("  Significant between-person effects:", sig_between, "\n\n")
  }
}

# EMA vs Baseline
if (!is.null(comparison) && is.list(comparison)) {
  cat("4. EMA VS BASELINE COMPARISON\n")
  cat(rep("-", 70), "\n", sep = "")
  if ("r2_comparison" %in% names(comparison)) {
    r2_comp <- comparison$r2_comparison
    if ("R2_Diff" %in% names(r2_comp)) {
      mean_diff <- mean(r2_comp$R2_Diff, na.rm = TRUE)
      cat("  Mean ΔR²:", sprintf("%.3f", mean_diff), "\n")
      if (abs(mean_diff) < 0.02) {
        cat("  → Models are EQUIVALENT\n\n")
      } else if (mean_diff > 0) {
        cat("  → EMA shows better fit\n\n")
      } else {
        cat("  → Baseline shows better fit\n\n")
      }
    }
  }
}

# Dataset info
if (!is.null(dataset)) {
  cat("5. DATASET\n")
  cat(rep("-", 70), "\n", sep = "")
  cat("  N observations:", nrow(dataset), "\n")
  if ("ID" %in% names(dataset)) {
    cat("  N subjects:", n_distinct(dataset$ID), "\n")
  }
  if ("timepoint" %in% names(dataset)) {
    cat(
      "  Timepoints:",
      paste(unique(dataset$timepoint), collapse = ", "),
      "\n"
    )
  }
  cat("\n")
}

cat(rep("=", 70), "\n", sep = "")
cat("NEXT STEPS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("1. Review the CSV files in results/exports/\n")
cat("2. Key file to examine: 02_moderation_results.csv\n")
cat("3. Share relevant CSV files with Claude for interpretation\n")
cat("4. Use summaries for manuscript writing\n\n")

sink()

cat("✓ Master summary saved: results/exports/00_MASTER_SUMMARY.txt\n\n")

# ==============================================================================
# FINAL MESSAGE
# ==============================================================================

cat(rep("=", 70), "\n", sep = "")
cat("EXPORT COMPLETE\n")
cat(rep("=", 70), "\n\n")

cat("All results exported to: results/exports/\n\n")

cat("Key files to review:\n")
cat("  → 00_MASTER_SUMMARY.txt - Start here!\n")
cat("  → 02_moderation_results.csv - Main findings\n")
cat("  → 01_ema_validation.csv - Validation results\n\n")

cat("Total files created:", length(list.files("results/exports")), "\n\n")

cat(rep("=", 70), "\n\n", sep = "")
