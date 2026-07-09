# RSV data cleaning
# Input: RSV_23_24_raw_merged
# Output: rsvdb (cleaned pathogen DB object)

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

if (!exists("RSV_23_24_raw_merged")) {
  stop("Object 'RSV_23_24_raw_merged' is missing. Source RSV_SQLquery.R first.")
}

safe_utf8 <- function(x) {
  x <- as.character(x)
  y <- iconv(x, from = "", to = "UTF-8", sub = NA)
  bad <- is.na(y)
  if (any(bad)) y[bad] <- iconv(x[bad], from = "latin1", to = "UTF-8", sub = "")
  y[is.na(y)] <- ""
  y
}

safe_parse_date <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  parsed <- suppressWarnings(lubridate::parse_date_time(x, orders = c("Y-m-d", "d.m.Y", "d/m/Y", "Y/m/d")))
  as.Date(parsed)
}

RSV_23_24_clean <- RSV_23_24_raw_merged %>%
  filter(is.na(ngs_report) | trimws(ngs_report) == "") %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), stringr::regex("^\\s*(?:3|P3(?:_.*)?)\\s*$", ignore_case = TRUE))) %>%
  filter(!stringr::str_detect(dplyr::coalesce(prove_kategori, ""), stringr::regex("ref", ignore_case = TRUE))) %>%
  select(-any_of(c("levelid", "prkey", "endtcreated", "endtmodif", "objactionid", "objlck", "objowner", "objshared"))) %>%
  mutate(across(where(is.character), safe_utf8)) %>%
  normalize_geography_columns()

if ("prove_tatt" %in% names(RSV_23_24_clean)) {
  RSV_23_24_clean <- RSV_23_24_clean %>%
    mutate(
      prove_tatt = na_if(prove_tatt, ""),
      prove_tatt = safe_parse_date(prove_tatt),
      week = lubridate::week(prove_tatt),
      year = lubridate::year(prove_tatt),
      wy = sprintf("%d-W%02d", lubridate::isoyear(prove_tatt), lubridate::isoweek(prove_tatt)),
      my = format(prove_tatt, "%Y-%m")
    )
}

# Final pathogen DB object name kept as requested
rsvdb <- RSV_23_24_clean

if (exists("close_sql_connections")) close_sql_connections()
