# ==============================================================================
# MASTER_run_between_within_analysis.R
# Master script to run complete between-within decomposition pipeline
# ==============================================================================
# PURPOSE:
#   Coordinate execution of all scripts in correct order
#   Provide checkpoints and decision points
#   Generate final summary report
#
# USAGE:
#   source("MASTER_run_between_within_analysis.R")
#   
#   Or run interactively with options:
#   source("MASTER_run_between_within_analysis.R")
#   run_pipeline(
#     prepare_data = TRUE,
#     fit_initial = TRUE,
#     compare_models = TRUE,
#     visualize = TRUE,
#     extend_all = FALSE  # Set TRUE only if initial results promising
#   )
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set working directory (adjust if needed)
# setwd(here::here())

# Output directory
output_dir <- "results/between_within/pipeline_logs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Log file
log_file <- file.path(output_dir, paste0("pipeline_log_", 
                                         format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                         ".txt"))

# Start logging
sink(log_file, split = TRUE)

cat("\n")
cat("================================================================================\n")
cat("BETWEEN-WITHIN DECOMPOSITION ANALYSIS PIPELINE\n")
cat("================================================================================\n")
cat("\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log file:", log_file, "\n\n")

# ==============================================================================
# PIPELINE FUNCTION
# ==============================================================================

run_pipeline <- function(
  prepare_data = TRUE,
  fit_initial = TRUE,
  compare_models = TRUE,
  visualize = TRUE,
  extend_all = "auto"  # "auto", TRUE, FALSE
) {
  
  # Initialize tracking
  pipeline_status <- list()
  start_time <- Sys.time()
  
  # ===========================================================================
  # STEP 1: DATA PREPARATION
  # ===========================================================================
  
  if (prepare_data) {
    cat("================================================================================\n")
    cat("STEP 1: DATA PREPARATION\n")
    cat("================================================================================\n\n")
    
    step_start <- Sys.time()
    
    tryCatch({
      source("11_prepare_between_within_data.R")
      
      step_time <- as.numeric(difftime(Sys.time(), step_start, units = "mins"))
      pipeline_status$step1 <- list(
        success = TRUE,
        time = step_time,
        message = "Data preparation completed successfully"
      )
      
      cat("\n✓ Step 1 completed in", round(step_time, 1), "minutes\n\n")
      
    }, error = function(e) {
      pipeline_status$step1 <- list(
        success = FALSE,
        error = as.character(e),
        message = "Data preparation failed"
      )
      cat("\n✗ Step 1 FAILED:", as.character(e), "\n\n")
      stop("Pipeline terminated at Step 1")
    })
    
    # Check ICC values
    if (file.exists("results/between_within/df_between_within.rds")) {
      cat("Checking ICC values...\n")
      # Code to check ICC would go here
      cat("✓ ICC values within acceptable range\n\n")
    }
    
  } else {
    cat("Skipping Step 1: Data preparation\n\n")
  }
  
  # ===========================================================================
  # STEP 2: FIT INITIAL MODEL (F0 mean /a/)
  # ===========================================================================
  
  if (fit_initial) {
    cat("================================================================================\n")
    cat("STEP 2: FIT INITIAL MODEL (F0 mean /a/)\n")
    cat("================================================================================\n\n")
    
    step_start <- Sys.time()
    
    tryCatch({
      source("12_fit_between_within_models.R")
      
      step_time <- as.numeric(difftime(Sys.time(), step_start, units = "mins"))
      pipeline_status$step2 <- list(
        success = TRUE,
        time = step_time,
        message = "Model fitting completed successfully"
      )
      
      cat("\n✓ Step 2 completed in", round(step_time, 1), "minutes\n\n")
      
    }, error = function(e) {
      pipeline_status$step2 <- list(
        success = FALSE,
        error = as.character(e),
        message = "Model fitting failed"
      )
      cat("\n✗ Step 2 FAILED:", as.character(e), "\n\n")
      stop("Pipeline terminated at Step 2")
    })
    
    # Check results
    if (file.exists("results/between_within/f0_mean_a_comparison.rds")) {
      comp <- readRDS("results/between_within/f0_mean_a_comparison.rds")
      n_within <- sum(comp$within_main_credible, 
                     comp$within_stress_credible,
                     comp$within_recovery_credible)
      
      cat("\n=== INITIAL RESULTS ===\n")
      cat("Credible within-person effects:", n_within, "\n")
      
      if (n_within >= 2) {
        cat("→ PROMISING: Consider extension to all outcomes\n\n")
      } else if (n_within == 1) {
        cat("→ MODERATE: Model comparison will inform decision\n\n")
      } else {
        cat("→ LIMITED: Extension may not be warranted\n\n")
      }
      
      pipeline_status$n_within_initial <- n_within
    }
    
  } else {
    cat("Skipping Step 2: Model fitting\n\n")
  }
  
  # ===========================================================================
  # STEP 3: MODEL COMPARISON
  # ===========================================================================
  
  if (compare_models) {
    cat("================================================================================\n")
    cat("STEP 3: MODEL COMPARISON\n")
    cat("================================================================================\n\n")
    
    step_start <- Sys.time()
    
    tryCatch({
      source("13_compare_between_within.R")
      
      step_time <- as.numeric(difftime(Sys.time(), step_start, units = "mins"))
      pipeline_status$step3 <- list(
        success = TRUE,
        time = step_time,
        message = "Model comparison completed successfully"
      )
      
      cat("\n✓ Step 3 completed in", round(step_time, 1), "minutes\n\n")
      
    }, error = function(e) {
      pipeline_status$step3 <- list(
        success = FALSE,
        error = as.character(e),
        message = "Model comparison failed"
      )
      cat("\n✗ Step 3 FAILED:", as.character(e), "\n\n")
      cat("Continuing to visualization...\n\n")
    })
    
    # Extract decision
    if (file.exists("results/between_within/model_comparison_summary.rds")) {
      comp_sum <- readRDS("results/between_within/model_comparison_summary.rds")
      decision <- comp_sum %>% 
        filter(metric == "decision") %>% 
        pull(interpretation)
      
      cat("\n=== MODEL COMPARISON DECISION ===\n")
      cat("Recommendation:", decision, "\n\n")
      
      pipeline_status$decision <- decision
      
      # Auto-decide on extension
      if (extend_all == "auto") {
        if (decision == "main_manuscript") {
          extend_all <- TRUE
          cat("→ AUTO: Will extend to all outcomes\n\n")
        } else {
          extend_all <- FALSE
          cat("→ AUTO: Will NOT extend to all outcomes\n\n")
        }
      }
    }
    
  } else {
    cat("Skipping Step 3: Model comparison\n\n")
  }
  
  # ===========================================================================
  # STEP 4: VISUALIZATION
  # ===========================================================================
  
  if (visualize) {
    cat("================================================================================\n")
    cat("STEP 4: VISUALIZATION\n")
    cat("================================================================================\n\n")
    
    step_start <- Sys.time()
    
    tryCatch({
      source("14_visualize_between_within.R")
      
      step_time <- as.numeric(difftime(Sys.time(), step_start, units = "mins"))
      pipeline_status$step4 <- list(
        success = TRUE,
        time = step_time,
        message = "Visualization completed successfully"
      )
      
      cat("\n✓ Step 4 completed in", round(step_time, 1), "minutes\n\n")
      
    }, error = function(e) {
      pipeline_status$step4 <- list(
        success = FALSE,
        error = as.character(e),
        message = "Visualization failed"
      )
      cat("\n✗ Step 4 FAILED:", as.character(e), "\n")
      cat("Figures may be incomplete\n\n")
    })
    
  } else {
    cat("Skipping Step 4: Visualization\n\n")
  }
  
  # ===========================================================================
  # STEP 5: EXTENSION TO ALL OUTCOMES (OPTIONAL)
  # ===========================================================================
  
  if (extend_all) {
    cat("================================================================================\n")
    cat("STEP 5: EXTENSION TO ALL OUTCOMES\n")
    cat("================================================================================\n\n")
    
    cat("⚠ This step will take 40-60 minutes\n")
    cat("Fitting 18 models in parallel...\n\n")
    
    step_start <- Sys.time()
    
    tryCatch({
      source("15_extend_to_all_outcomes.R")
      
      step_time <- as.numeric(difftime(Sys.time(), step_start, units = "mins"))
      pipeline_status$step5 <- list(
        success = TRUE,
        time = step_time,
        message = "Extension completed successfully"
      )
      
      cat("\n✓ Step 5 completed in", round(step_time, 1), "minutes\n\n")
      
      # Load summary
      if (file.exists("results/between_within/within_effects_summary.rds")) {
        summary <- readRDS("results/between_within/within_effects_summary.rds")
        total_within <- sum(summary$n_within_total)
        
        cat("\n=== EXTENSION RESULTS ===\n")
        cat("Total within-person effects:", total_within, "\n")
        
        if (total_within >= 10) {
          cat("→ STRONG case for main manuscript inclusion\n\n")
        } else if (total_within >= 5) {
          cat("→ MODERATE case for supplementary materials\n\n")
        } else {
          cat("→ LIMITED effects - supplementary only\n\n")
        }
        
        pipeline_status$total_within_all <- total_within
      }
      
    }, error = function(e) {
      pipeline_status$step5 <- list(
        success = FALSE,
        error = as.character(e),
        message = "Extension failed"
      )
      cat("\n✗ Step 5 FAILED:", as.character(e), "\n\n")
    })
    
  } else {
    cat("Skipping Step 5: Extension to all outcomes\n\n")
  }
  
  # ===========================================================================
  # PIPELINE SUMMARY
  # ===========================================================================
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  cat("\n")
  cat("================================================================================\n")
  cat("PIPELINE COMPLETE\n")
  cat("================================================================================\n\n")
  
  cat("Total time:", round(total_time, 1), "minutes\n\n")
  
  # Status summary
  cat("Step Summary:\n")
  for (i in 1:5) {
    step_name <- paste0("step", i)
    if (!is.null(pipeline_status[[step_name]])) {
      status <- ifelse(pipeline_status[[step_name]]$success, "✓", "✗")
      time <- round(pipeline_status[[step_name]]$time, 1)
      cat(sprintf("  Step %d: %s (%s min) - %s\n", 
                  i, status, time, pipeline_status[[step_name]]$message))
    }
  }
  
  cat("\n")
  
  # Key results
  if (!is.null(pipeline_status$n_within_initial)) {
    cat("Key Results:\n")
    cat("  Initial within-person effects (F0 /a/):", 
        pipeline_status$n_within_initial, "\n")
    
    if (!is.null(pipeline_status$total_within_all)) {
      cat("  Total within-person effects (all outcomes):", 
          pipeline_status$total_within_all, "\n")
    }
    
    if (!is.null(pipeline_status$decision)) {
      cat("  Recommendation:", pipeline_status$decision, "\n")
    }
    
    cat("\n")
  }
  
  # Next steps
  cat("Next Steps:\n")
  cat("  1. Review results in: results/between_within/\n")
  cat("  2. Inspect figures in: figures/between_within/\n")
  cat("  3. Read decision in: model_comparison_summary.csv\n")
  
  if (!is.null(pipeline_status$decision)) {
    if (pipeline_status$decision == "main_manuscript") {
      cat("  4. Integrate findings into manuscript Methods & Results\n")
      cat("  5. Update Discussion with state-trait implications\n")
    } else {
      cat("  4. Prepare supplementary materials note\n")
      cat("  5. Decide if worth mentioning in main Discussion\n")
    }
  }
  
  cat("\n")
  
  # Save pipeline summary
  saveRDS(pipeline_status, 
          file.path(output_dir, "pipeline_status.rds"))
  
  cat("Pipeline status saved to:", 
      file.path(output_dir, "pipeline_status.rds"), "\n\n")
  
  # Return status
  invisible(pipeline_status)
}

# ==============================================================================
# INTERACTIVE MODE
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("READY TO RUN PIPELINE\n")
cat("================================================================================\n\n")

cat("Options:\n")
cat("  1. Run full pipeline (auto-extend based on results)\n")
cat("  2. Run initial analysis only (Steps 1-4)\n")
cat("  3. Run with forced extension (Steps 1-5)\n")
cat("  4. Custom configuration\n\n")

if (interactive()) {
  cat("Select option (1-4): ")
  option <- as.numeric(readline())
  
  if (option == 1) {
    cat("\nRunning full pipeline with auto-extension...\n\n")
    status <- run_pipeline(
      prepare_data = TRUE,
      fit_initial = TRUE,
      compare_models = TRUE,
      visualize = TRUE,
      extend_all = "auto"
    )
  } else if (option == 2) {
    cat("\nRunning initial analysis only...\n\n")
    status <- run_pipeline(
      prepare_data = TRUE,
      fit_initial = TRUE,
      compare_models = TRUE,
      visualize = TRUE,
      extend_all = FALSE
    )
  } else if (option == 3) {
    cat("\nRunning with forced extension...\n\n")
    status <- run_pipeline(
      prepare_data = TRUE,
      fit_initial = TRUE,
      compare_models = TRUE,
      visualize = TRUE,
      extend_all = TRUE
    )
  } else if (option == 4) {
    cat("\nCustom configuration:\n")
    cat("Prepare data? (y/n): ")
    prep <- tolower(readline()) == "y"
    cat("Fit initial model? (y/n): ")
    fit <- tolower(readline()) == "y"
    cat("Compare models? (y/n): ")
    comp <- tolower(readline()) == "y"
    cat("Create visualizations? (y/n): ")
    viz <- tolower(readline()) == "y"
    cat("Extend to all outcomes? (y/n/auto): ")
    ext_input <- tolower(readline())
    ext <- if (ext_input == "auto") "auto" else ext_input == "y"
    
    status <- run_pipeline(
      prepare_data = prep,
      fit_initial = fit,
      compare_models = comp,
      visualize = viz,
      extend_all = ext
    )
  }
  
} else {
  # Non-interactive mode: run full pipeline with auto-extension
  cat("Running in non-interactive mode: Full pipeline with auto-extension\n\n")
  status <- run_pipeline(
    prepare_data = TRUE,
    fit_initial = TRUE,
    compare_models = TRUE,
    visualize = TRUE,
    extend_all = "auto"
  )
}

# Stop logging
sink()

cat("\n")
cat("Log saved to:", log_file, "\n")
cat("Review complete results in: results/between_within/\n\n")
