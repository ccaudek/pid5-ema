# Script purpose: Read individual EMA data files and save an RDS file.
# @author: Ilaria Colpizzi <ilaria.colpizzi@unifi.it>

suppressPackageStartupMessages({
  library(psych)
  library(ggplot2)
  library(sjPlot)
  library(dplyr)
  library(rio)
  library(here)
  library(purrr)
  library(lubridate)
  library(readxl)
  library(janitor)
  library(stringr)
})

source(here::here("data", "raw_complete", "funs_ema_mpath.R"))

# All individual EMA files are stored in this directory.
# dir <- here::here("data", "raw", "mpath2023")

### This is a test with only correct files!!!!
### It must be replaced with the folder name containing all the data files,
### when all files will be corrected.
# FOLDER_NAME <- "ema_data"

# The first row of each xlsx file is empyt. It must be deleted. But
# this operation must only be performed once. To insure that, I use
# the following trick.
# if (file.exists(here::here("data", "raw", FOLDER_NAME, "boop.txt"))) {
#   # Delete the first row of each xlsx file.
#   delete_first_row_xlsx(FOLDER_NAME)
#   file.remove(here::here("data", "raw", FOLDER_NAME, "boop.txt"))
# }

# Folder with all raw data files.
dir_path <- here::here("data", "raw_complete", "ema_quest")

# prendo i nomi file e i path completi in modo consistente
file_names <- list.files(
  path = dir_path,
  pattern = "\\.xlsx$",
  full.names = FALSE
)
file_paths <- file.path(dir_path, file_names)

n_files <- length(file_names)
n_files

d_list <- vector("list", n_files)
errors <- tibble::tibble(
  file = character(),
  step = character(),
  message = character()
)

safe_import <- function(path) {
  # 1) tentativo con rio
  x <- tryCatch(
    rio::import(path, skip = 1),
    error = function(e) e
  )
  if (!inherits(x, "error")) return(x)

  # 2) fallback con readxl (a volte riesce dove rio no)
  y <- tryCatch(
    readxl::read_excel(path, skip = 1),
    error = function(e) e
  )
  if (!inherits(y, "error")) return(as.data.frame(y))

  # se entrambi falliscono, restituisco l'errore di rio (più informativo per .xlsx)
  x
}

for (index_file in seq_along(file_paths)) {
  fp <- file_paths[index_file]
  fn <- file_names[index_file]

  d <- safe_import(fp)
  if (inherits(d, "error")) {
    errors <- dplyr::add_row(
      errors,
      file = fn,
      step = "import",
      message = conditionMessage(d)
    )
    next
  }

  # subject code (toglie ".xlsx" alla fine)
  d$subj_code <- sub("\\.xlsx$", "", fn)

  # parsing datetime: più robusto di strptime
  if (!("Date and time" %in% names(d))) {
    errors <- dplyr::add_row(
      errors,
      file = fn,
      step = "datetime",
      message = "Missing column `Date and time`"
    )
    next
  }

  d$datetime <- suppressWarnings(
    lubridate::parse_date_time(
      d$`Date and time`,
      orders = c("a, d b Y H:M:S", "d b Y H:M:S", "Y-m-d H:M:S"),
      tz = "UTC"
    )
  )

  # se non parse-a nulla, segna e salta
  if (all(is.na(d$datetime))) {
    errors <- dplyr::add_row(
      errors,
      file = fn,
      step = "datetime",
      message = "All datetime values are NA after parsing"
    )
    next
  }

  d$day <- as.Date(d$datetime)
  d$hour <- lubridate::hour(d$datetime)
  d$minute <- lubridate::minute(d$datetime)

  # calendar_day e bysubj_day
  # (calendar_day: indice di giorno nell'intero dataset del file; bysubj_day: dense_rank dentro al soggetto)
  d$calendar_day <- as.integer(factor(d$day))

  d <- d %>%
    mutate(bysubj_day = dense_rank(calendar_day)) %>%
    select(-`Date and time`, -datetime)

  d_list[[index_file]] <- d
}

# tieni solo gli elementi non-null (file letti)
d_list <- purrr::compact(d_list)

# verifica: quanti file letti vs falliti
length(d_list)
nrow(errors)

# se vuoi vedere quali file hanno problemi
errors %>% arrange(step, file)


#
#
# d_list <- list()
#
# for (index_file in 1:n_files) {
#
#   d <-
#     rio::import(
#       here(dir, file_names[index_file]), skip = 1
#       #header = FALSE
#     )
#
#   # Get subject code
#   d$subj_code <-
#     substr(file_names[index_file], 1, nchar(file_names[index_file]) - 5)
#
#   # Change columns' names.
#   # d1 <- change_cols_names(d)
#   # d <- d1
#
#   # Get time window
#   d$datetime <- strptime(d$`Date and time`, "%a, %d %b %Y %H:%M:%S")
#   d$day <- ymd(substring(d$datetime, 1, 10))
#   d$hour <- hour(d$datetime)
#   d$minute <- minute(d$datetime)
#
#   # d$time_window <- character(length(d$datetime))
#   #
#   # for (i in seq_along(d$datetime)) {
#   #   if (d$hour[i] >= 0 & d$hour[i] < 12) {
#   #     d$time_window[i] <- 1 # "10-11"
#   #   } else if (d$hour[i] >= 12 & d$hour[i] < 16) {
#   #     d$time_window[i] <- 2 # "14-15"
#   #   } else if (d$hour[i] >= 16 & d$hour[i] < 19) {
#   #     d$time_window[i] <- 3 # "17-18"
#   #   } else if (d$hour[i] >= 19 & d$hour[i] < 21) {
#   #     d$time_window[i] <- 4 # "19-20"
#   #   } else if (d$hour[i] >= 21 & d$hour[i] < 24) {
#   #     d$time_window[i] <- 5 # "21-23"
#   #   } else {
#   #     d$time_window[i] <- "NA"
#   #   }
#   # }
#   # table(d$user_id, d$time_window)
#
#   foo <- d$day |>
#     as.numeric()
#   d$calendar_day <- as.numeric(factor(foo))
#
#   d1 <- d %>%
#     mutate(bysubj_day = dense_rank(calendar_day))
#
#   # data.frame(
#   #   d1$user_id, d1$calendar_day, d1$day, d1$bysubj_day
#   # )[1:200, ]
#
#   d2 <- d1 |>
#     dplyr::select(
#       -c(`Date and time`, datetime)
#     )
#
#   # d3 <- d2 %>%
#   #   janitor::clean_names()
#   #
#   # if ("momentary_emotion_no_multi_smiley" %in% names(d3)) {
#   #   d3$momentary_emotion_no_multi_smiley <-
#   #     as.character(d3$momentary_emotion_no_multi_smiley)
#   # }
#   #
#   # if ("momentary_emotion_yes_multi_smiley" %in% names(d3)) {
#   #   d3$momentary_emotion_yes_multi_smiley <-
#   #     as.character(d3$momentary_emotion_yes_multi_smiley)
#   # }
#
#   d_list[[index_file]] <- d2
# }

# Get a vector of unique column names across all data frames
all_columns <- unique(unlist(lapply(d_list, names)))

# Create a new data frame with consistent columns and fill missing columns with NA
combined_df <- do.call(
  bind_rows,
  lapply(d_list, function(df) {
    missing_columns <- setdiff(all_columns, names(df))
    df[, missing_columns] <- NA
    return(df)
  })
)

combined_df1 <- combined_df %>%
  dplyr::select(
    -all_of(c(
      "tripm_3 (1) (multipleChoice)",
      "labelExample (basic)",
      "reminder_voice_recording (basic)",
      "system_devspecs (devicespecs)",
      "text (basic)",
      "text3 (basic)",
      "text_2 (basic)",
      "umore (smiley)",
      "minute"
    ))
  ) %>%
  mutate(
    subj_code = recode(
      subj_code,
      "alro20051109715" = "al_ro_2005_11_09_715",
      "clpe20000217166" = "cl_pe_2000_02_17_166",
      "fla_me_2004_02_17_358" = "fl_me_2004_02_17_358",
      "ga_ga_2002_11_26_429 copia" = "ga_ga_2002_11_26_429",
      "gi_st_2005_07_311" = "gi_st_2005_07_22_311",
      "luri20051012907" = "lu_ri_2005_10_12_907",
      "ma_qu_09_15_711" = "ma_qu_2003_09_15_711",
      "maie20000731515" = "ma_ie_2000_07_31_515",
      "masc040120409" = "ma_sc_04_01_20_409",
      "masi20050211538" = "ma_si_2005_02_11_538",
      "reve20041021036" = "re_ve_2004_10_21_036",
      "so_fa_1997_04_21_013 copia" = "so_fa_1997_04_21_013",
      "st_ma_04_21_426" = "st_ma_2004_04_21_426",
      "vi_da_1999_11_499" = "vi_da_1999_02_11_499",
      "tommaso vaselli" = "to_va_2005_11_12_377",
      "laura.mal" = "la_ma_2004_05_16_992",
      "jennifer" = "je_ma_2005_01_13_962",
      "Gianluca" = "gi_be_2001_01_18_537",
      "Elisa lupetti del re" = "el_lu_2005_12_21_023",
      "CostanzaCiappi" = "co_ci_2003_07_29_927",
      "an_gr_2003_02_23" = "an_gr_2003_02_23_266",
      "elto" = "el_to_2005_09_02_232",
      "NiccoloVenturi" = "ni_ve_2003_12_26_353"
    )
  )

combined_df1$subj_code <- tolower(combined_df1$subj_code)

# "noemi lorena"
# "Alessandro"
# "Alice"
# "alino"
# "Alessia"
# "arianna"
# "Francesco"

# Save complete raw data
# rio::export(
#   combined_df1,
#   here::here("_data", "08.07.2025", "voice_data", "lookup_files", "ema_data_controls.xlsx")
# )

scs_to_convert <- c(
  "scs1_pos (multipleChoice)",
  "scs2_neg (multipleChoice)",
  "scs3_pos (multipleChoice)",
  "scs4_neg (multipleChoice)",
  "scs5_neg (multipleChoice)",
  "scs6_pos (multipleChoice)",
  "scs7_pos (multipleChoice)",
  "scs8_neg (multipleChoice)"
)

combined_df2 <- combined_df1 %>%
  mutate(across(
    all_of(scs_to_convert),
    ~ case_when(
      . == "Totalmente falso" ~ -3,
      . == "Moderatamente falso" ~ -2,
      . == "Leggermente falso" ~ -1,
      . == "Leggermente vero" ~ 1,
      . == "Moderatamente vero" ~ 2,
      . == "Totalmente vero" ~ 3
    )
  ))

pid5_cols <- c(
  "pid5_13 (multipleChoice)",
  "pid5_15 (multipleChoice)",
  "pid5_11 (multipleChoice)",
  "pid5_3 (multipleChoice)",
  "pid5_2 (multipleChoice)",
  "pid5_7 (multipleChoice)",
  "pid5_14 (multipleChoice)",
  "pid5_6 (multipleChoice)",
  "pid5_4 (multipleChoice)",
  "pid5_12 (multipleChoice)",
  "pid5_1 (multipleChoice)",
  "pid5_9 (multipleChoice)",
  "pid5_5 (multipleChoice)",
  "pid5_8 (multipleChoice)",
  "pid5_10 (multipleChoice)"
)

combined_df3 <- combined_df2 %>%
  mutate(across(
    all_of(pid5_cols),
    ~ case_when(
      str_trim(.) == "Questa affermazione è totalmente falsa" ~ 0,
      str_trim(.) == "Questa affermazione è abbastanza falsa" ~ 1,
      str_trim(.) == "Questa affermazione è abbastanza vera" ~ 2,
      str_trim(.) == "Questa affermazione è totalmente vera" ~ 3
    )
  ))

tripm_cols <- c(
  "tripm_1 (multipleChoice)",
  "tripm_3 (multipleChoice)",
  "tripm_2 (multipleChoice)",
  "tripm_4 (multipleChoice)"
)

combined_df4 <- combined_df3 %>%
  mutate(across(
    all_of(tripm_cols),
    ~ case_when(
      str_trim(.) == "Questa affermazione è falsa" ~ 1,
      str_trim(.) == "Questa affermazione è abbastanza falsa" ~ 2,
      str_trim(.) == "Questa affermazione è abbastanza vera" ~ 3,
      str_trim(.) == "Questa affermazione è vera" ~ 4
    )
  ))

dass_cols <- c(
  "dass21_2 (multipleChoice)",
  "dass21_5 (multipleChoice)",
  "dass21_3 (multipleChoice)",
  "dass21_4 (multipleChoice)",
  "dass21_6 (multipleChoice)",
  "dass21_1 (multipleChoice)"
)

combined_df5 <- combined_df4 %>%
  mutate(across(
    all_of(dass_cols),
    ~ case_when(
      str_trim(.) == "Questa affermazione è falsa" ~ 0,
      str_trim(.) == "Questa affermazione è abbastanza falsa" ~ 1,
      str_trim(.) == "Questa affermazione è abbastanza vera" ~ 2,
      str_trim(.) == "Questa affermazione è vera" ~ 3
    )
  ))

vq_cols <- c(
  "vq_2 (multipleChoice)",
  "vq_1 (multipleChoice)",
  "vq_4 (multipleChoice)",
  "vq_3 (multipleChoice)"
)

combined_df6 <- combined_df5 %>%
  mutate(across(
    all_of(vq_cols),
    ~ case_when(
      str_trim(.) == "Questa affermazione è falsa" ~ 0,
      str_trim(.) == "Questa affermazione è abbastanza falsa" ~ 1,
      str_trim(.) == "Questa affermazione è abbastanza vera" ~ 2,
      str_trim(.) == "Questa affermazione è vera" ~ 3
    )
  ))

cope_cols <- c(
  "cope_nvi_2 (multipleChoice)",
  "cope_nvi_7 (multipleChoice)",
  "cope_nvi_5 (multipleChoice)",
  "cope_nvi_8 (multipleChoice)",
  "cope_nvi_9 (multipleChoice)",
  "cope_nvi_10 (multipleChoice)",
  "cope_nvi_1 (multipleChoice)",
  "cope_nvi_3 (multipleChoice)",
  "cope_nvi_4 (multipleChoice)",
  "cope_nvi_6 (multipleChoice)"
)

combined_df7 <- combined_df6 %>%
  mutate(across(
    all_of(cope_cols),
    ~ case_when(
      str_trim(.) == "Di solito non lo faccio" ~ 1,
      str_trim(.) == "Lo faccio qualche volta" ~ 2,
      str_trim(.) == "Lo faccio con una certa frequenza" ~ 3,
      str_trim(.) == "Lo faccio quasi sempre" ~ 4
    )
  ))

combined_df7 <- combined_df7 %>%
  mutate(
    context_quality = recode(
      `context1 (multipleChoice)`,
      "Per niente favorevole" = -2,
      "Poco favorevole" = -1,
      "Né favorevole né sfavorevole" = 0,
      "Moderatamente favorevole" = 1,
      "Molto favorevole" = 2
    )
  ) %>%
  mutate(
    context_control = recode(
      `context2 (multipleChoice)`,
      "Per niente" = 0,
      "Poco" = 1,
      "Parzialmente" = 2,
      "Abbastanza" = 3,
      "Molto" = 4,
      "Completamente" = 5
    )
  ) %>%
  mutate(
    context_support = recode(
      `context3 (multipleChoice)`,
      "Per niente" = 0,
      "Poco" = 1,
      "Moderatamente" = 2,
      "Abbastanza" = 3,
      "Molto" = 4,
      "Completamente" = 5,
      "Non applicabile (la situazione non richiede supporto)" = 999
    )
  ) %>%
  mutate(
    context_threat = recode(
      `context4 (multipleChoice)`,
      "Per niente" = 0,
      "Poco" = 1,
      "Moderatamente" = 2,
      "Abbastanza" = 3,
      "Molto" = 4,
      "Estremamente" = 5
    )
  )


combined_df7$`context1 (multipleChoice)` <- NULL
combined_df7$`context2 (multipleChoice)` <- NULL
combined_df7$`context3 (multipleChoice)` <- NULL
combined_df7$`context4 (multipleChoice)` <- NULL


pattern <- " ?\\([^()]*\\)"

combined_df7 <- combined_df7 %>%
  rename_with(~ str_remove_all(.x, pattern), matches("\\("))

combined_df7 <- combined_df7 %>%
  mutate(context_support = na_if(context_support, 999))

# Save complete raw data
saveRDS(
  combined_df7,
  here::here("data", "raw_complete", "data_interim", "ema_data_clean.RDS")
)

# eof ---
