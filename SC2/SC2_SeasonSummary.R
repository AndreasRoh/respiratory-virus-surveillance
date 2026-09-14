# Generate an introduction table for the normal SARS-CoV-2 reporting season.

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(stringr)
})
source(file.path("Source_files", "common_report_utils.R"))

# Override with SC2_REPORTING_SEASON_START_YEAR when preparing a later season.
reporting_season_start_year <- suppressWarnings(as.integer(
  Sys.getenv("SC2_REPORTING_SEASON_START_YEAR", unset = "2025")
))
if (is.na(reporting_season_start_year)) {
  stop("SC2_REPORTING_SEASON_START_YEAR must be a four-digit year.")
}

iso_week_start <- function(iso_year, iso_week) {
  jan_fourth <- as.Date(sprintf("%d-01-04", iso_year))
  week_one_monday <- jan_fourth - (as.integer(format(jan_fourth, "%u")) - 1L)
  week_one_monday + (iso_week - 1L) * 7L
}

season_start <- iso_week_start(reporting_season_start_year, 35L)
season_end <- iso_week_start(reporting_season_start_year + 1L, 34L) + 6L
season_label <- sprintf("%d/%02d", reporting_season_start_year, (reporting_season_start_year + 1L) %% 100L)

source(file.path("Source_files", "common_report_utils.R"))
source(file.path("SC2", "SC2_SQLquery_25-26.R"))

# The selected season is wholly in the current 2025/26 database. Supplying an
# empty historical extract lets the established current-report cleaning script
# apply its normal inclusion rules without bringing in legacy observations.
SC2_20_25 <- SC2_25_26_raw_merged[0, ]
source(file.path("SC2", "SC2_DataCleaning.R"))

season_data <- SC2db %>%
  mutate(
    sample_date = as.Date(prove_tatt),
    pangolin = trimws(as.character(nc_pangolin_short)),
    gisaid_isolate = trimws(as.character(gisaid_isolate_name)),
    overvaking = case_when(
      str_detect(trimws(as.character(prove_kategori)), "^(P1\\b|P1_|1\\b)") ~ "Sentinel",
      TRUE ~ "Generell overv\u00e5king"
    ),
    aldergruppe_raw = as.character(
      age_to_group_standard(suppressWarnings(as.numeric(trimws(as.character(pasient_alder)))))
    ),
    aldergruppe = case_when(
      is.na(aldergruppe_raw) | trimws(aldergruppe_raw) == "" ~ "Ukjent/ikke oppgitt",
      TRUE ~ aldergruppe_raw
    ),
    pasientstatus = case_when(
      str_detect(toupper(trimws(as.character(pasient_status))), "^INNELIGGENDE") ~ "Innlagt (hospitalised)",
      str_detect(toupper(trimws(as.character(pasient_status))), "^POLIKLINISK") ~ "Poliklinisk",
      TRUE ~ "Ukjent/ikke oppgitt"
    ),
    fylke = case_when(
      is.na(pasient_fylke_name) | trimws(as.character(pasient_fylke_name)) == "" ~ "Ukjent/ikke oppgitt",
      TRUE ~ as.character(pasient_fylke_name)
    ),
    has_pangolin = !is.na(pangolin) & nzchar(pangolin),
    has_gisaid_isolate = !is.na(gisaid_isolate) & nzchar(gisaid_isolate),
    maaned_start = as.Date(format(sample_date, "%Y-%m-01"))
  ) %>%
  filter(sample_date >= season_start, sample_date <= season_end)

season_data <- normalize_sex_column(season_data, candidate_cols = c("pasient_kjnn", "pasient_kjonn")) %>%
  mutate(
    kjonn = recode(
      pasient_kjonn_std,
      Female = "Kvinne",
      Male = "Mann",
      Ukjent = "Ukjent/ikke oppgitt",
      .default = "Ukjent/ikke oppgitt"
    )
  )

if (anyDuplicated(season_data$key)) {
  stop("The normal-season extract contains duplicate virus keys.")
}

total_sequences <- nrow(season_data)
if (total_sequences == 0L) {
  stop("No records found for the selected reporting season.")
}

distribution_table <- function(data, variable, levels = NULL) {
  if (!is.null(levels)) {
    out <- tibble(
      Kategori = levels,
      Antall = as.integer(table(factor(as.character(data[[variable]]), levels = levels)))
    )
  } else {
    out <- data %>%
      count(Kategori = .data[[variable]], name = "Antall") %>%
      arrange(desc(Antall), Kategori)
  }

  out %>% mutate(Andel = Antall / total_sequences)
}

overview_table <- tibble(
  Del = "Sekvensering",
  Kategori = c(
    "Totalt antall sekvenserte SARS-CoV-2-virus",
    "Pangolin-klassifisert",
    "GISAID-isolatnavn"
  ),
  Antall = c(
    total_sequences,
    sum(season_data$has_pangolin),
    sum(season_data$has_gisaid_isolate)
  ),
  Andel = Antall / total_sequences
)

overvaking_table <- distribution_table(season_data, "overvaking", c("Sentinel", "Generell overv\u00e5king"))
kjonn_table <- distribution_table(season_data, "kjonn", c("Kvinne", "Mann", "Ukjent/ikke oppgitt"))
alder_table <- distribution_table(season_data, "aldergruppe", c("0-4", "5-14", "15-24", "25-59", "60+", "Ukjent/ikke oppgitt"))
pasientstatus_table <- distribution_table(
  season_data,
  "pasientstatus",
  c("Innlagt (hospitalised)", "Poliklinisk", "Ukjent/ikke oppgitt")
)
fylke_table <- distribution_table(season_data, "fylke")

introduction_table <- bind_rows(
  overview_table,
  mutate(overvaking_table, Del = "Overv\u00e5king") %>% select(Del, Kategori, Antall, Andel),
  mutate(kjonn_table, Del = "Kj\u00f8nn") %>% select(Del, Kategori, Antall, Andel),
  mutate(alder_table, Del = "Aldersfordeling") %>% select(Del, Kategori, Antall, Andel),
  mutate(pasientstatus_table, Del = "Pasientstatus") %>% select(Del, Kategori, Antall, Andel)
)

season_month_starts <- seq(as.Date(format(season_start, "%Y-%m-01")), as.Date(format(season_end, "%Y-%m-01")), by = "month")
monthly_distribution <- function(data, variable, categories) {
  monthly_counts <- data %>%
    count(Kategori = .data[[variable]], maaned_start, name = "Antall")

  tidyr::expand_grid(Kategori = categories, maaned_start = season_month_starts) %>%
    left_join(monthly_counts, by = c("Kategori", "maaned_start")) %>%
    mutate(Antall = coalesce(Antall, 0L))
}

introduction_data <- introduction_table
monthly_for_overview <- function(variable, section) {
  season_data %>%
    count(Kategori = .data[[variable]], maaned_start, name = "Antall") %>%
    mutate(Del = section) %>%
    select(Del, Kategori, maaned_start, Antall)
}

introduction_monthly_counts <- bind_rows(
  season_data %>%
    count(maaned_start, name = "Antall") %>%
    transmute(Del = "Sekvensering", Kategori = "Totalt antall sekvenserte SARS-CoV-2-virus", maaned_start, Antall),
  season_data %>%
    filter(has_pangolin) %>%
    count(maaned_start, name = "Antall") %>%
    transmute(Del = "Sekvensering", Kategori = "Pangolin-klassifisert", maaned_start, Antall),
  season_data %>%
    filter(has_gisaid_isolate) %>%
    count(maaned_start, name = "Antall") %>%
    transmute(Del = "Sekvensering", Kategori = "GISAID-isolatnavn", maaned_start, Antall),
  monthly_for_overview("overvaking", "Overv\u00e5king"),
  monthly_for_overview("kjonn", "Kj\u00f8nn"),
  monthly_for_overview("aldergruppe", "Aldersfordeling"),
  monthly_for_overview("pasientstatus", "Pasientstatus")
)
introduction_monthly <- introduction_data %>%
  distinct(Del, Kategori) %>%
  tidyr::crossing(maaned_start = season_month_starts) %>%
  left_join(introduction_monthly_counts, by = c("Del", "Kategori", "maaned_start")) %>%
  mutate(Antall = coalesce(Antall, 0L), trend_key = paste(Del, Kategori, sep = "||"))

section_values <- unique(introduction_data$Del)
introduction_table <- bind_rows(lapply(section_values, function(section_value) {
  bind_rows(
    tibble(
      Kategori = section_value,
      Antall = NA_integer_,
      Andel = NA_real_,
      is_section = TRUE,
      trend_key = NA_character_
    ),
    introduction_data %>%
      filter(Del == section_value) %>%
      transmute(Kategori, Antall, Andel, is_section = FALSE, trend_key = paste(Del, Kategori, sep = "||"))
  )
}))
metadata_table <- tibble(
  Felt = c("Sesong", "Periode", "Datapopulasjon", "Sentinel-definisjon", "Pasientstatus"),
  Verdi = c(
    sprintf("%s (uke 35 til og med uke 34)", season_label),
    sprintf("%s til %s", format(season_start, "%d.%m.%Y"), format(season_end, "%d.%m.%Y")),
    "Normal SARS-CoV-2-rapport: normal ngs_report-regel, ikke Proveprosjekt 3, ikke referanseprove, NC-coverage >= 0,7 eller NA, og ikke-tom spike-mutasjon.",
    "Prove-kategori som starter med P1, P1_ eller 1.",
    "Innlagt = INNELIGGENDE; poliklinisk = POLIKLINISK; resterende verdier = ukjent/ikke oppgitt."
  )
)

workbook <- createWorkbook()
fhi_light_blue <- fhi_colour("blue_grid")
fhi_sparkline_blue <- fhi_colour("text_muted")
table_border_colour <- fhi_colour("blue_border")
table_border_style <- "dotted"
header_style <- createStyle(
  fontName = "Calibri",
  fontSize = 9,
  fontColour = fhi_colour("ink"),
  fgFill = fhi_light_blue,
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = c("Top", "Bottom", "Left", "Right"),
  borderColour = table_border_colour,
  borderStyle = table_border_style
)
section_style <- createStyle(
  fontName = "Calibri",
  fontSize = 9,
  fontColour = fhi_colour("ink"),
  fgFill = fhi_light_blue,
  textDecoration = "bold",
  valign = "center",
  border = c("Top", "Bottom", "Left", "Right"),
  borderColour = table_border_colour,
  borderStyle = table_border_style
)
body_style <- createStyle(
  fontName = "Calibri",
  fontSize = 9,
  fontColour = fhi_colour("ink"),
  valign = "center",
  border = c("Top", "Bottom", "Left", "Right"),
  borderColour = table_border_colour,
  borderStyle = table_border_style
)
title_style <- createStyle(fontName = "Calibri", fontSize = 9, fontColour = fhi_colour("ink"), textDecoration = "bold")
count_style <- createStyle(numFmt = "0", halign = "right")
percent_style <- createStyle(numFmt = "0.0%", halign = "right")

with_percent_header <- function(table) {
  names(table)[names(table) == "Andel"] <- "Andel (%)"
  table
}

insert_monthly_sparklines <- function(sheet_name, table, monthly_data, key_column = "Kategori") {
  trend_image_dir <- file.path(tempdir(), "sc2_season_summary_sparklines")
  dir.create(trend_image_dir, recursive = TRUE, showWarnings = FALSE)

  for (idx in seq_len(nrow(table))) {
    trend_key <- table[[key_column]][idx]
    if (is.na(trend_key)) {
      next
    }

    trend_data <- monthly_data %>%
      filter(.data[[key_column]] == .env$trend_key) %>%
      arrange(maaned_start)
    image_file <- tempfile(pattern = "trend_", tmpdir = trend_image_dir, fileext = ".png")

    grDevices::png(image_file, width = 390, height = 42, res = 150, bg = "transparent")
    graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
    max_count <- max(c(trend_data$Antall, 1L))
    graphics::plot(
      seq_len(nrow(trend_data)),
      trend_data$Antall,
      type = "n",
      axes = FALSE,
      xlab = "",
      ylab = "",
      ylim = c(0, max_count * 1.1)
    )
    graphics::abline(h = 0, col = fhi_light_blue, lwd = 1)
    graphics::lines(seq_len(nrow(trend_data)), trend_data$Antall, col = fhi_sparkline_blue, lwd = 1.4)
    grDevices::dev.off()

    insertImage(
      workbook,
      sheet_name,
      image_file,
      startRow = idx + 3,
      startCol = 4,
      width = 2.6,
      height = 0.2,
      units = "in",
      dpi = 150
    )
  }
}

write_distribution_sheet <- function(sheet_name, title, table, monthly_data = NULL) {
  export_table <- with_percent_header(table)
  has_sparklines <- !is.null(monthly_data)
  if (has_sparklines) {
    export_table[["Mnd. utvikling (antall)"]] <- ""
  }

  addWorksheet(workbook, sheet_name)
  writeData(workbook, sheet_name, title, startRow = 1, startCol = 1)
  addStyle(workbook, sheet_name, title_style, rows = 1, cols = 1)
  writeData(workbook, sheet_name, export_table, startRow = 3, startCol = 1, headerStyle = header_style)
  addStyle(workbook, sheet_name, body_style, rows = 4:(nrow(table) + 3), cols = seq_len(ncol(export_table)), gridExpand = TRUE)
  addStyle(workbook, sheet_name, count_style, rows = 4:(nrow(table) + 3), cols = 2, gridExpand = TRUE, stack = TRUE)
  addStyle(workbook, sheet_name, percent_style, rows = 4:(nrow(table) + 3), cols = 3, gridExpand = TRUE, stack = TRUE)
  setColWidths(
    workbook,
    sheet_name,
    cols = seq_len(ncol(export_table)),
    widths = if (has_sparklines) c(42, 14, 14, 31) else c(42, 14, 14)
  )
  if (has_sparklines) {
    setRowHeights(workbook, sheet_name, rows = 4:(nrow(table) + 3), heights = 20)
    insert_monthly_sparklines(sheet_name, table, monthly_data)
  }
  freezePane(workbook, sheet_name, firstActiveRow = 4)
}

write_overview_sheet <- function() {
  sheet_name <- "Innledning"
  export_table <- introduction_table %>%
    select(-is_section, -trend_key) %>%
    with_percent_header()
  export_table[["Mnd. utvikling (antall)"]] <- ""

  addWorksheet(workbook, sheet_name)
  writeData(workbook, sheet_name, sprintf("SARS-CoV-2 - sesongoversikt %s", season_label), startRow = 1, startCol = 1)
  addStyle(workbook, sheet_name, title_style, rows = 1, cols = 1)
  writeData(workbook, sheet_name, export_table, startRow = 3, startCol = 1, headerStyle = header_style)
  addStyle(workbook, sheet_name, body_style, rows = 4:(nrow(introduction_table) + 3), cols = 1:4, gridExpand = TRUE)
  section_rows <- which(introduction_table$is_section) + 3
  value_rows <- which(!introduction_table$is_section) + 3
  addStyle(workbook, sheet_name, section_style, rows = section_rows, cols = 1:4, gridExpand = TRUE)
  addStyle(workbook, sheet_name, count_style, rows = value_rows, cols = 2, gridExpand = TRUE, stack = TRUE)
  addStyle(workbook, sheet_name, percent_style, rows = value_rows, cols = 3, gridExpand = TRUE, stack = TRUE)
  setColWidths(workbook, sheet_name, cols = 1:4, widths = c(46, 14, 14, 31))
  setRowHeights(workbook, sheet_name, rows = 4:(nrow(introduction_table) + 3), heights = 20)
  insert_monthly_sparklines(sheet_name, introduction_table, introduction_monthly, key_column = "trend_key")
  freezePane(workbook, sheet_name, firstActiveRow = 4)
}
write_overview_sheet()
write_distribution_sheet("Fylkesfordeling", "Fylkesfordeling", fylke_table, monthly_distribution(season_data, "fylke", fylke_table$Kategori))
write_distribution_sheet("Kj\u00f8nn", "Kj\u00f8nnsfordeling", kjonn_table, monthly_distribution(season_data, "kjonn", kjonn_table$Kategori))
write_distribution_sheet("Alder", "Aldersfordeling", alder_table, monthly_distribution(season_data, "aldergruppe", alder_table$Kategori))
write_distribution_sheet("Pasientstatus", "Pasientstatus", pasientstatus_table, monthly_distribution(season_data, "pasientstatus", pasientstatus_table$Kategori))

addWorksheet(workbook, "Metadata")
writeData(workbook, "Metadata", metadata_table, startRow = 1, startCol = 1, headerStyle = header_style)
setColWidths(workbook, "Metadata", cols = 1:2, widths = c(24, 110))
addStyle(workbook, "Metadata", body_style, rows = 2:(nrow(metadata_table) + 1), cols = 1:2, gridExpand = TRUE)
addStyle(workbook, "Metadata", createStyle(wrapText = TRUE, valign = "top"), rows = 2:(nrow(metadata_table) + 1), cols = 1:2, gridExpand = TRUE, stack = TRUE)

results_dir <- file.path(dirname(dirname(dirname(normalizePath(getwd())))), "Results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
output_filename <- Sys.getenv("SC2_SEASON_SUMMARY_OUTPUT_FILE", unset = sprintf("SARSCoV2_%s_sesongoversikt.xlsx", gsub("/", "-", season_label)))
output_path <- if (report_audit_mode()) {
  file.path(report_output_dir(results_dir, "season_summaries"), output_filename)
} else {
  file.path(results_dir, output_filename)
}
saveWorkbook(workbook, output_path, overwrite = TRUE)
message("Wrote ", output_path)
