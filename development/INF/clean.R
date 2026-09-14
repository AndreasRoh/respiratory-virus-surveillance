# =============================================================================
# WORK IN PROGRESS — reviewed, behaviour-preserving copy
#
# Purpose: This file is a maintainability-focused review copy of INF_DataCleaning_25-26.R.
# Changes in the WIP series are limited to structure, documentation, and WIP
# dependency isolation. The calculations, filters, object names, and exported
# outputs are retained so results can be compared directly with production.
# =============================================================================

# INF 25-26 data cleaning
# Input: INF_25_26_raw_merged
# Output: INF_25_26_clean, INF_25_26_sequences, fludb

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

if (!exists("INF_25_26_raw_merged")) {
  stop("Object 'INF_25_26_raw_merged' is missing. Source INF_SQLquery_25-26.R first.")
}

# The temporary Sesongrapport includes only GISAID-submitted sequences, while the
# normal report retains its existing ngs_report screening rule.
temporary_sesongrapport <- tolower(Sys.getenv("INF_TEMPORARY_SESONGRAPPORT", unset = "false")) %in% c("1", "true", "yes")
gisaid_submission_cols <- intersect(
  c("gisaid_isolate_id", "gisaid_ha_id", "gisaid_na_id", "gisaid_m_id", "gisaid_ns_id", "gisaid_np_id", "gisaid_pa_id", "gisaid_pb1_id", "gisaid_pb2_id"),
  names(INF_25_26_raw_merged)
)
if (temporary_sesongrapport && length(gisaid_submission_cols) == 0) {
  stop("Temporary Sesongrapport requires at least one GISAID identifier column.")
}
gisaid_submitted <- rep(FALSE, nrow(INF_25_26_raw_merged))
if (length(gisaid_submission_cols) > 0) {
  gisaid_submitted <- apply(INF_25_26_raw_merged[, gisaid_submission_cols, drop = FALSE], 1, function(values) {
    any(!is.na(values) & nzchar(trimws(as.character(values))))
  })
}

INF_25_26_clean <- INF_25_26_raw_merged %>%
  mutate(
    prove_tatt = as.Date(prove_tatt, format = "%Y-%m-%d"),
    pasient_alder = suppressWarnings(as.numeric(trimws(as.character(pasient_alder)))),
    .sesongrapport_gisaid_submitted = gisaid_submitted,
    .sesongrapport_include = if (temporary_sesongrapport) .sesongrapport_gisaid_submitted else (is.na(ngs_report) | trimws(ngs_report) == "")
  ) %>%
  filter(.sesongrapport_include) %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), stringr::regex("^\\s*(?:3|P3(?:_.*)?)\\s*$", ignore_case = TRUE))) %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), stringr::regex("ref", ignore_case = TRUE))) %>%
  filter(trimws(coalesce(tessy_reportable_variable, "")) != "") %>%
  filter(!stringr::str_detect(dplyr::coalesce(tessy_reportable_variable, ""), stringr::regex("ref", ignore_case = TRUE))) %>%
  select(-.sesongrapport_gisaid_submitted, -.sesongrapport_include) %>%
  as.data.frame() %>%

  normalize_geography_columns()
seq_data_raw <- tbl(conFLU2526, "SEQUENCEDATA") %>%
  collect() %>%
  janitor::clean_names()

seq_filtered <- seq_data_raw %>%
  inner_join(INF_25_26_raw_merged %>% select(key), by = "key") %>%
  filter(grepl("01-HA|02-NA|03-M|04-PB1|05-PB2|07-PA|06-NP|08-NS", experiment)) %>%
  mutate(
    experiment = stringr::str_remove(experiment, "\\d+-"),
    experiment = ifelse(experiment == "M", "MP", experiment)
  ) %>%
  filter(type == "SEQUENCE") %>%
  select(key, experiment, data)

process_entry <- function(data_gz_base64) {
  tryCatch({
    data_gz_raw <- base64enc::base64decode(data_gz_base64)
    data_decompressed <- memDecompress(data_gz_raw, type = "gzip")
    data_decompressed_no_nulls <- data_decompressed[data_decompressed != as.raw(0)]
    data_text_no_nulls <- rawToChar(data_decompressed_no_nulls)
    stringr::str_remove(data_text_no_nulls, "B0t")
  }, error = function(e) {
    message("Error processing entry: ", e$message)
    NA
  })
}

INF_25_26_sequences <- seq_filtered %>%
  mutate(sequence = purrr::map_chr(data, process_entry)) %>%
  mutate(sequence = stringr::str_sub(sequence, start = 3)) %>%
  select(key, experiment, sequence)

# Final pathogen DB object name
fludb <- INF_25_26_clean

if (exists("close_sql_connections")) close_sql_connections()

rm(process_entry, seq_data_raw, seq_filtered)
