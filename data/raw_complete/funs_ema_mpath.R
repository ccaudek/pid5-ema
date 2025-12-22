delete_first_row_xlsx <- function(FOLDER_NAME) {
  library(readxl)
  library(writexl)
  library(fs)

  # Directory path where Excel files are located
  directory_path <- here("data", "raw", FOLDER_NAME)

  # List all Excel files in the directory
  excel_files <- dir_ls(directory_path, regexp = "\\.xlsx$")

  # Loop through each Excel file
  for (file in excel_files) {
    # Read the Excel file
    data <- read_excel(file, skip = 1)

    # Get the original column names
    # original_column_names <- readxl::col_names(readxl::read_excel(file, n_max = 1))

    # Set the original column names to the modified data
    # colnames(data) <- original_column_names

    # Write the modified data back to the Excel file
    write_xlsx(data, path = file)

    # Print the processed file name
    cat("Deleted first row in:", basename(file), "\n")
  }
}


# Change the columns' names of the file with the raw data.
change_cols_names <- function(d) {
  # col_names <- c(
  #   "date_time", "context",
  #   "scs_pos_1", "happy",
  #   "dec_2", "dec_4",
  #   "scs_neg_5", "scs_pos_3",
  #   "satisfied", "scs_neg_8",
  #   "dec_3", "scs_neg_4",
  #   "dec_1", "scs_neg_2",
  #   "nervous", "scs_pos_6",
  #   "upset", "scs_pos_7",
  #   "present_emotion", "exam",
  #   "after_exam_emotion", "exam_preparation",
  #   "grade_exp", "before_exam_emotion",
  #   "before_exam_emotion_intensity", "real_grade",
  #   "after_exam_emotion_intensity", "user_id"
  # )

  col_names <- c(
    "context", "scs_pos_1",
    "happy", "dec_2",
    "dec_4", "scs_neg_5",
    "scs_pos_3", "satisfied",
    "scs_neg_8", "dec_3",
    "scs_neg_4", "dec_1",
    "scs_neg_2", "nervous",
    "scs_pos_6", "upset",
    "scs_pos_7", "emotion_now",
    "exam", "emotion_after_exam",
    "exam_preparation", "grade_exp",
    "emotion_before_exam", "emotion_intensity_before_exam",
    "real_grade", "emotion_intensity_after_exam",
    "user_id", "day",
    "hour", "minute",
    "time_window", "calendar_day",
    "bysubj_day"
  )

  names(d) <- col_names

  return(d)
}


center3L <- function(dataname, varname, idname, dayname) {
  within <- dataname %>%
    group_by({{ idname }}, {{ dayname }}) %>%
    mutate("{{varname}}_DddM" := mean({{ varname }}, na.rm = TRUE)) %>%
    mutate(DddM = mean({{ varname }}, na.rm = TRUE)) %>%
    mutate("{{varname}}_Moment" := {{ varname }} - mean({{ varname }}, na.rm = TRUE)) %>%
    ungroup()
  within2 <- dataname %>%
    group_by({{ idname }}, {{ dayname }}) %>%
    summarize(temp_dddm = mean({{ varname }}, na.rm = TRUE)) %>%
    summarize("{{varname}}_PppM" := mean(temp_dddm, na.rm = TRUE)) %>%
    ungroup()
  within3 <- dataname %>%
    group_by({{ idname }}, {{ dayname }}) %>%
    summarize(temp_dddm = mean({{ varname }}, na.rm = TRUE)) %>%
    summarize(PppM = mean(temp_dddm, na.rm = TRUE)) %>%
    ungroup()
  combinewithin <- merge(within2, within3)
  allwithin <- merge(within, combinewithin)
  allwithinc <- allwithin %>% 
    mutate("{{varname}}_Day" := DddM - PppM)
  between <- dataname %>%
    group_by({{ idname }}, {{ dayname }}) %>%
    summarize(temp_dddm = mean({{ varname }}, na.rm = TRUE)) %>%
    summarize(temp_pppm = mean(temp_dddm, na.rm = TRUE)) %>%
    summarize("{{varname}}_SssM" := mean(temp_pppm, na.rm = TRUE)) %>%
    ungroup()
  between2 <- dataname %>%
    group_by({{ idname }}, {{ dayname }}) %>%
    summarize(temp_dddm = mean({{ varname }}, na.rm = TRUE)) %>%
    summarize(temp_pppm = mean(temp_dddm, na.rm = TRUE)) %>%
    summarize(SssM = mean(temp_pppm, na.rm = TRUE)) %>%
    ungroup()
  between3 <- full_join(between, between2, by = character())
  combined <- full_join(allwithinc, between3, by = character())
  output <- combined %>% 
    mutate("{{varname}}_Person" := PppM - SssM)
  out <- select(output, -c(DddM, PppM, SssM))
  return(out)
}



#--------------------------------------------------------------------

# https://osf.io/fdnm7/files/osfstorage

reportMLM <- function(model, decimal = 3, intercept = FALSE, format = "apa", partR2 = FALSE, showTable = TRUE) {
  # Last modified by Hause Lin 31-10-17 20:06
  # APA format
  # cat("Requires p-values from lmerTest if model has been fitted using lme4\nFixed effects for MLM with effect sizes r, d, and semi-partial R (only works with lme and lmer fitted models)\n")
  # cat("r: .10 (small), .30 (medium), .50 (large) (Cohen, 1992)\n")
  # cat("d: 0.20 (small), 0.50 (medium), .80 (large) (Cohen, 1992)\n")
  # cat("R2: .02 (small), .13 (medium), .26 (large) (Cohen, 1992)\n")

  estimates <- data.frame(coef(summary(model)))
  effectNames <- rownames(estimates)
  colnames(estimates) <- tolower(colnames(estimates))
  colnames(estimates) <- c("estimate", "se", "df", "statistic", "p")

  # effect sizes (Kashdan & Steger, 2006)
  estimates$es.r <- sqrt((estimates$statistic^2 / (estimates$statistic^2 + estimates$df))) # r
  estimates$es.d <- (2 * estimates$statistic) / sqrt(estimates$df) # d

  # if partR2 == TRUE, compute partial R square (DOESN'T WORK IF MODEL HAS CATEGORICAL VARAIBLES WITH > 2 LEVELS!)
  if (partR2) {
    # effect size semi partial R (Edwards et al., 2008)
    anovaModel <- anova(model)
    Fs <- anovaModel$F # F-values for each effect (marginal = type 3 SS with Satterthwaite (requires lmerTest package))
    numDF <- anovaModel$NumDF # numerator DFs
    denDF <- anovaModel$DenDF # denominator DFs
    semiPartialREffect <- (numDF / denDF * Fs) / (1 + (numDF / denDF * Fs)) # effect sizes
    estimates$es.partR2 <- c(NA, semiPartialREffect)
  }

  # print formatted results
  if (intercept) {
    startingRow <- 1 # if intercept == TRUE, then report intercept too
  } else {
    startingRow <- 2 # if intercept == FALSE, then don't report intercept
  }

  for (i in startingRow:nrow(estimates)) { # for each fixed effect, format output accordingly (print line by line)

    # i <- i
    estimateV <- estimates[i, "estimate"]
    estimateV <- ifelse(abs(estimateV) < 0.01, round(estimateV, 3), round(estimateV, 2))
    estimateSign <- sign(estimateV) # sign (negative or positive)
    estimateV <- ifelse(abs(estimateV) < 0.01,
      sprintf("%.3f", abs(estimateV)),
      sprintf("%.2f", abs(estimateV))
    ) # character (absolute value)

    if (.Platform$OS.type == "unix") { # if linux/mac, ensure negative sign is minus, not hyphens
      estimateV <- ifelse(estimateSign == -1, paste0("−", estimateV), estimateV) # if negative, add minus sign to character vector
    } else {
      estimateV <- ifelse(estimateSign == -1, paste0("-", estimateV), estimateV) # hyphen
    }

    seV <- estimates[i, "se"]
    seV <- ifelse(abs(seV) < 0.01, round(seV, 3), round(seV, 2))
    seV <- ifelse(abs(seV) < 0.01,
      sprintf("%.3f", abs(seV)),
      sprintf("%.2f", abs(seV))
    ) # character
    dfV <- round(estimates[i, "df"])
    statisticV <- estimates[i, "statistic"]

    statisticSign <- sign(statisticV)
    statisticV <- ifelse(abs(statisticV) < 0.01,
      sprintf("%.3f", abs(statisticV)),
      sprintf("%.2f", abs(statisticV))
    ) # character

    if (.Platform$OS.type == "unix") { # if linux/mac, ensure negative signs is minus, not hyphens
      statisticV <- ifelse(statisticSign == -1, paste0("−", statisticV), statisticV) # minus
    } else {
      statisticV <- ifelse(statisticSign == -1, paste0("-", statisticV), statisticV) # hyphen
    }

    pV <- estimates[i, "p"]
    rV <- estimates[i, "es.r"]

    if (pV >= 0.001) { # if p value is >= .001, always print exact p value

      if (format == "nn") { # if nature neuro format

        message(sprintf(
          "%.0f: b = %s, SE = %s, t%.0f = %s, P = %.3f, r = %.2f",
          i, estimateV, seV, dfV, statisticV, pV, rV
        ))
      } else { # else always return APA format

        message(sprintf(
          "%.0f: b = %s, SE = %s, t(%.0f) = %s, p = %s, r = %s",
          i, estimateV, seV, dfV, statisticV, substring(sprintf("%.3f", pV), 2), substring(sprintf("%.2f", rV), 2)
        ))
      }
    } else if (pV < 0.001) { # if p value < .001, then always print p < .001

      pV <- 0.001

      if (format == "nn") {
        message(sprintf(
          "%.0f: b = %s, SE = %s, t%.0f = %s, P < %.3f, r = %.2f",
          i, estimateV, seV, dfV, statisticV, pV, rV
        ))
      } else {
        message(sprintf(
          "%.0f: b = %s, SE = %s, t(%.0f) = %s, p < %s, r = %s",
          i, estimateV, seV, dfV, statisticV, substring(sprintf("%.3f", pV), 2), substring(sprintf("%.2f", rV), 2)
        ))
      }
    }
  }

  rownames(estimates) <- paste(1:nrow(estimates), rownames(estimates)) # add row number to summary table

  if (showTable) {
    return(round(estimates, decimal))
  }
}
