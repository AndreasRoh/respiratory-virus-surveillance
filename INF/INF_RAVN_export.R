resolve_script_dir <- function() {
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

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(stringr)
})

source(file.path(bundle_scripts_dir, "..", "Source_files", "common_report_utils.R"))
source(file.path(bundle_scripts_dir, "INF_SQLquery_25-26.R"))
source(file.path(bundle_scripts_dir, "INF_DataCleaning_25-26.R"))

if (!exists("fludb")) {
  stop("Object 'fludb' is missing. Source INF_DataCleaning_25-26.R first.")
}

output_dir <- "N:/Virologi/Influensa/2526/WGS_Analyse/RAVN"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

init_locale <- function() {
  locale_candidates <- c("nb_NO.UTF-8", "Norwegian (Bokmal)_Norway.utf8", "Norwegian")
  for (loc in locale_candidates) {
    ok <- tryCatch(Sys.setlocale(category = "LC_ALL", locale = loc), error = function(e) NA_character_)
    if (!is.na(ok) && nzchar(ok)) {
      return(invisible(ok))
    }
  }
  warning("Could not set Norwegian UTF-8 locale; continuing with system default locale.")
  invisible(NULL)
}

init_locale()

write_csv2_utf8_bom <- function(df, path, na = "") {
  tmp_path <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_path), add = TRUE)

  utils::write.table(
    df,
    file = tmp_path,
    sep = ";",
    dec = ",",
    row.names = FALSE,
    col.names = TRUE,
    qmethod = "double",
    na = na,
    fileEncoding = "UTF-8"
  )

  raw_data <- readBin(tmp_path, what = "raw", n = file.info(tmp_path)$size)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(raw_data, con)
}

write_utf8_text <- function(lines, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(lines), con = con, useBytes = FALSE)
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- str_replace_all(x, "[[:cntrl:]]", "")
  str_squish(x)
}

blank_to_na <- function(x) {
  x <- clean_text(x)
  x[x == ""] <- NA_character_
  x
}

normalize_sequence_result <- function(x) {
  x <- toupper(clean_text(x))
  case_when(
    x %in% c("H1", "A/H1", "A/H1N1", "AH1N1", "AH1N1PDM09") ~ "A/H1N1",
    x %in% c("H3", "A/H3", "A/H3N2", "AH3N2") ~ "A/H3N2",
    x %in% c("B/VICTORIA", "BVICTORIA", "BVIC", "B/VIC") ~ "B/Victoria",
    TRUE ~ ifelse(x == "", NA_character_, x)
  )
}

normalize_sex <- function(x) {
  x <- toupper(clean_text(x))
  case_when(
    x %in% c("M", "MALE", "MANN") ~ "M",
    x %in% c("F", "FEMALE", "KVINNE") ~ "F",
    TRUE ~ "UNK"
  )
}

normalize_material <- function(x) {
  x <- toupper(clean_text(x))
  case_when(
    x == "SE" ~ "Sekret",
    x == "" ~ NA_character_,
    TRUE ~ x
  )
}

normalize_resistance <- function(x) {
  x <- blank_to_na(x)
  ifelse(is.na(x), "NA", x)
}

normalize_mutation <- function(x) {
  x <- blank_to_na(x)
  ifelse(is.na(x), "NA", x)
}

derive_na_subtype <- function(x) {
  seq_result <- normalize_sequence_result(x)
  case_when(
    seq_result == "A/H1N1" ~ "N1",
    seq_result == "A/H3N2" ~ "N2",
    seq_result == "B/Victoria" ~ "B",
    TRUE ~ NA_character_
  )
}

derive_birth_year <- function(sample_date, age_years) {
  sample_date <- as.Date(sample_date)
  age_num <- suppressWarnings(as.numeric(age_years))
  birth_year <- ifelse(is.na(sample_date) | is.na(age_num), NA_real_, year(sample_date) - age_num)
  ifelse(is.na(birth_year), NA_character_, format(as.integer(birth_year), scientific = FALSE, trim = TRUE))
}

format_export_date <- function(x) {
  x <- as.Date(x)
  ifelse(is.na(x), NA_character_, format(x, "%d.%m.%Y"))
}

current_ravn_season_start_year <- function(reference_date = Sys.Date()) {
  reference_date <- as.Date(reference_date)
  iso_week <- isoweek(reference_date)
  iso_year <- isoyear(reference_date)
  if (iso_week >= 35) iso_year else iso_year - 1
}

filter_ravn_season <- function(df, start_year = current_ravn_season_start_year()) {
  df %>%
    mutate(
      prove_tatt = as.Date(prove_tatt),
      iso_week = isoweek(prove_tatt),
      iso_year = isoyear(prove_tatt)
    ) %>%
    filter(
      !is.na(prove_tatt),
      (iso_year == start_year & iso_week >= 35) |
        (iso_year == start_year + 1 & iso_week <= 20)
    )
}

translation_sheet <- tibble::tribble(
  ~old_export_column, ~new_export_column, ~source_column, ~transformation, ~note,
  "key", "key", "key", "Direct copy", "Unchanged sample ID.",
  "fylkenavn", "pasient_fylke_name", "pasient_fylke_name", "Text cleanup", "Uses the current county column.",
  "inf_res", "ngs_sekvens_resultat", "ngs_sekvens_resultat", "Normalizes H1/H3 to full subtype", "Current influenza type/subtype nomenclature.",
  "sted", "prove_innsender_navn", "prove_innsender_navn", "Text cleanup", "Current full sender name replaces old short code.",
  "ar", "pasient_birth_year", "prove_tatt + pasient_alder", "Derived as sample year minus age", "Exact birth year is not a standalone field in current fludb.",
  "", "pasient_alder", "pasient_alder", "Numeric cleanup", "New explicit field in the current export.",
  "k", "pasient_kjnn", "pasient_kjnn", "Standardized to M/F/UNK", "Current sex column.",
  "lokalisasjon", "prove_material", "prove_material", "SE -> Sekret", "Current material field is coarser than the old fludb field.",
  "na_subt", "na_subtype", "ngs_sekvens_resultat", "Derived from subtype", "N1 for A/H1N1, N2 for A/H3N2, B for B/Victoria.",
  "inf_strain_name", "gisaid_isolate_name", "gisaid_isolate_name", "Direct copy", "Current isolate name.",
  "interpretation_oseltamivir_resistance_testing", "dr_res_oseltamivir", "dr_res_oseltamivir", "Blank -> NA", "Current resistance column.",
  "interpretation_zanamivir_resistance_testing", "dr_res_zanamivir", "dr_res_zanamivir", "Blank -> NA", "Current resistance column.",
  "interpretation_m2blocker_resistance_testing", "dr_res_adamantine", "dr_res_adamantine", "Blank -> NA", "Current resistance column.",
  "resistens_mutasjon", "dr_na_mut", "dr_na_mut", "Blank -> NA", "Current NA mutation column.",
  "resistens_kommentar", "dr_kommentar", "dr_kommentar", "'.' kept, blank -> NA", "Comment from resistance assessment.",
  "interpretation_adamantane_resistanse_testing", "dr_res_adamantine", "dr_res_adamantine", "Blank -> NA", "Historical duplicate replaced by the current column.",
  "interpretation_balaxovir_resistance_testing", "dr_res_baloxavir", "dr_res_baloxavir", "Blank -> NA", "Old typo replaced by the current column name.",
  "interpretation_peramivir_testing", "dr_res_peramivir", "dr_res_peramivir", "Blank -> NA", "Current resistance column.",
  "materiale", "prove_material", "prove_material", "SE -> Sekret", "Same mapping as above.",
  "na_sek_acc", "gisaid_na_id", "gisaid_na_id", "Direct copy", "Current accession ID for NA.",
  "resistens", "dr_res_laninamivir", "dr_res_laninamivir", "Blank -> NA", "Old ambiguous summary column replaced with explicit laninamivir result.",
  "resistens_mut_m2", "dr_m2_mut", "dr_m2_mut", "Blank -> NA", "Current M2 mutation column.",
  "resistens_mut_pa", "dr_pa_mut", "dr_pa_mut", "Blank -> NA", "Current PA mutation column.",
  "sekvens_resultat", "tessy_reportable_variable", "tessy_reportable_variable", "Direct copy", "Current reportable nomenclature.",
  "seq_id_pa", "gisaid_pa_id", "gisaid_pa_id", "Direct copy", "Current accession ID for PA.",
  "", "nc_ha_subclade", "nc_ha_subclade", "Direct copy", "New explicit HA subclade field in the current export.",
  "", "nc_na_clade", "nc_na_clade", "Direct copy", "New explicit NA clade field in the current export.",
  "", "gisaid_isolate_id", "gisaid_isolate_id", "Direct copy", "New explicit isolate ID field in the current export.",
  "prove_tatt", "prove_tatt", "prove_tatt", "Formatted as dd.mm.yyyy", "Same field with standardized date format."
)

ravn_df <- fludb %>%
  filter_ravn_season() %>%
  transmute(
    key = clean_text(key),
    pasient_fylke_name = blank_to_na(pasient_fylke_name),
    ngs_sekvens_resultat = normalize_sequence_result(ngs_sekvens_resultat),
    prove_innsender_navn = blank_to_na(prove_innsender_navn),
    pasient_birth_year = derive_birth_year(prove_tatt, pasient_alder),
    pasient_alder = suppressWarnings(as.numeric(clean_text(pasient_alder))),
    pasient_kjnn = normalize_sex(pasient_kjnn),
    prove_material = normalize_material(prove_material),
    na_subtype = derive_na_subtype(ngs_sekvens_resultat),
    gisaid_isolate_name = blank_to_na(gisaid_isolate_name),
    dr_res_oseltamivir = normalize_resistance(dr_res_oseltamivir),
    dr_res_zanamivir = normalize_resistance(dr_res_zanamivir),
    dr_res_adamantine = normalize_resistance(dr_res_adamantine),
    dr_na_mut = normalize_mutation(dr_na_mut),
    dr_kommentar = blank_to_na(dr_kommentar),
    dr_res_baloxavir = normalize_resistance(dr_res_baloxavir),
    dr_res_peramivir = normalize_resistance(dr_res_peramivir),
    gisaid_na_id = blank_to_na(gisaid_na_id),
    dr_res_laninamivir = normalize_resistance(dr_res_laninamivir),
    dr_m2_mut = normalize_mutation(dr_m2_mut),
    dr_pa_mut = normalize_mutation(dr_pa_mut),
    tessy_reportable_variable = blank_to_na(tessy_reportable_variable),
    gisaid_pa_id = blank_to_na(gisaid_pa_id),
    nc_ha_subclade = blank_to_na(nc_ha_subclade),
    nc_na_clade = blank_to_na(nc_na_clade),
    gisaid_isolate_id = blank_to_na(gisaid_isolate_id),
    prove_tatt = format_export_date(prove_tatt)
  ) %>%
  mutate(
    across(
      where(is.character),
      normalize_norwegian_text
    )
  ) %>%
  arrange(
    dmy(prove_tatt),
    pasient_fylke_name,
    ngs_sekvens_resultat,
    key
  )

timestamp_label <- tolower(format(Sys.time(), "%d-%b-%Y %H.%M"))
export_path <- file.path(output_dir, paste0(timestamp_label, "_RAVN_INFcsv2.csv"))
translation_csv_path <- file.path(output_dir, paste0(timestamp_label, "_RAVN_INF_translation_sheet.csv"))
translation_txt_path <- file.path(output_dir, paste0(timestamp_label, "_RAVN_INF_translation_sheet.txt"))

write_csv2_utf8_bom(ravn_df, export_path, na = "")
translation_sheet <- translation_sheet %>%
  mutate(
    across(
      where(is.character),
      normalize_norwegian_text
    )
  )
write_csv2_utf8_bom(translation_sheet, translation_csv_path, na = "")

translation_email_lines <- c(
  "Translation old RAVN influenza export -> new export with current fludb nomenclature:",
  paste0(
    "- ",
    ifelse(translation_sheet$old_export_column == "", "[ny]", translation_sheet$old_export_column),
    " -> ",
    translation_sheet$new_export_column,
    " (kilde: ",
    translation_sheet$source_column,
    ifelse(
      translation_sheet$note == "",
      ")",
      paste0("; note: ", translation_sheet$note, ")")
    )
  )
)
write_utf8_text(translation_email_lines, translation_txt_path)

message("Wrote RAVN export: ", export_path)
message("Wrote translation sheet CSV: ", translation_csv_path)
message("Wrote translation sheet TXT: ", translation_txt_path)
message("Rows exported: ", nrow(ravn_df))
