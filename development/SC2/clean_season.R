# =============================================================================
# SC2 current-season data cleaning (season 2025/26)
# Input:  SC2_25_26_raw_merged and SC2_20_25
# Output: SC2_25_26_clean, SC2_25_26_prefilter, SC2_25_26,
#         SC2db_prefilter, SC2db
# =============================================================================

resolve_script_dir <- function() {
  frame_paths <- character()
  frame_list <- sys.frames()
  if (length(frame_list) > 0) {
    for (idx in rev(seq_along(frame_list))) {
      ofile <- tryCatch(frame_list[[idx]]$ofile, error = function(e) "")
      if (!is.null(ofile) && nzchar(ofile)) {
        frame_paths <- c(
          frame_paths,
          normalizePath(ofile, winslash = "/", mustWork = FALSE)
        )
      }
    }
  }
  frame_paths <- unique(frame_paths[nzchar(frame_paths)])
  if (length(frame_paths) > 0) {
    return(dirname(frame_paths[[1]]))
  }

  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
  }

  this_file <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) ""
  )
  if (nzchar(this_file)) {
    return(dirname(this_file))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

bundle_scripts_dir <- resolve_script_dir()

if (!exists("normalize_geography_columns")) {
  source(file.path(bundle_scripts_dir, "..", "Source_files", "common_report_utils.R"))
}

if (!exists("SC2_25_26_raw_merged")) {
  stop("Object 'SC2_25_26_raw_merged' is missing. Source SC2_SQLquery_25-26.R first.")
}
if (!exists("SC2_20_25")) {
  stop("Object 'SC2_20_25' is missing. Source SC2_DataCleaning_BNCOVID19.R first.")
}

# Source-level text normalization is centralized in common_report_utils.R.
normalize_sc2_text <- normalize_norwegian_text

exclude_prove_kategori_pattern <- stringr::regex("^\\s*(?:3|P3(?:_.*)?)\\s*$", ignore_case = TRUE)
exclude_reference_pattern <- stringr::regex("ref", ignore_case = TRUE)

# -----------------------------------------------------------------------------
# Clean and harmonize the current-season extract before report filtering.
# -----------------------------------------------------------------------------
SC2_25_26_clean <- SC2_25_26_raw_merged %>%
  select(-any_of(c("levelid", "prkey", "endtcreated", "endtmodif", "objactionid", "objlck", "objowner", "objshared"))) %>%
  filter(is.na(ngs_report) | trimws(ngs_report) == "") %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), exclude_prove_kategori_pattern)) %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), exclude_reference_pattern)) %>%
  mutate(
    mut_s_1 = na_if(mut_s_1, "NA"),
    mut_s_2 = na_if(mut_s_2, "NA"),
    mut_s_3 = na_if(mut_s_3, "NA"),
    mut_s_4 = na_if(mut_s_4, "NA")
  ) %>%
  tidyr::unite("spike_mut", mut_s_1, mut_s_2, mut_s_3, mut_s_4, sep = ";", na.rm = TRUE) %>%
  mutate(
    week = lubridate::week(as.Date(prove_tatt)),
    year = lubridate::year(as.Date(prove_tatt)),
    wy = tsibble::yearweek(as.Date(prove_tatt)),
    my = tsibble::yearmonth(as.Date(prove_tatt))
  ) %>%
  mutate(across(where(is.character), normalize_sc2_text)) %>%
  filter(prove_tatt != "") %>%
  normalize_geography_columns()

SC2_25_26_prefilter <- SC2_25_26_clean

# -----------------------------------------------------------------------------
# Apply report inclusion criteria for the current season.
# -----------------------------------------------------------------------------
SC2_25_26 <- SC2_25_26_clean %>%
  filter((nc_coverage == "NA" | nc_coverage >= 0.7) & spike_mut != "")

# -----------------------------------------------------------------------------
# Append historical rows that are not already present in the current-season data.
# -----------------------------------------------------------------------------
keys_in_sc2_25_26 <- SC2_25_26$key

SC2_20_25_filtered <- SC2_20_25 %>%
  filter(!key %in% keys_in_sc2_25_26) %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), exclude_prove_kategori_pattern)) %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), exclude_reference_pattern))

SC2_25_26 <- bind_rows(SC2_25_26, SC2_20_25_filtered) %>%
  mutate(across(where(is.character), normalize_sc2_text)) %>%
  normalize_geography_columns()

# Standard report objects used downstream.
SC2db_prefilter <- SC2_25_26_prefilter
SC2db <- SC2_25_26

if (exists("close_sql_connections")) close_sql_connections()

rm(keys_in_sc2_25_26)