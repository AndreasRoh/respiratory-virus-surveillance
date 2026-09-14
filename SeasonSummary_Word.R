# Create a Word version of the influenza and SARS-CoV-2 season overview tables.

suppressPackageStartupMessages({
  library(flextable)
  library(officer)
  library(openxlsx)
})
source(file.path("Source_files", "common_report_utils.R"))

results_dir <- file.path(dirname(dirname(dirname(normalizePath(getwd())))), "Results")
flu_source <- Sys.getenv(
  "INF_SEASON_SUMMARY_WORD_SOURCE",
  unset = file.path(results_dir, "Influenza_2025-26_sesongoversikt.xlsx")
)
sc2_source <- Sys.getenv(
  "SC2_SEASON_SUMMARY_WORD_SOURCE",
  unset = file.path(results_dir, "SARSCoV2_2025-26_sesongoversikt.xlsx")
)
rsv_source <- Sys.getenv(
  "RSV_SEASON_SUMMARY_WORD_SOURCE",
  unset = file.path(results_dir, "RSV_2025-26_sesongoversikt.xlsx")
)
output_path <- Sys.getenv(
  "SEASON_SUMMARY_WORD_OUTPUT_FILE",
  unset = file.path(results_dir, "Sesongoversikt_2025-26_influensa_sc2_og_rsv.docx")
)
if (report_audit_mode()) {
  output_path <- file.path(report_output_dir(dirname(output_path), "season_summary_word"), basename(output_path))
}

extract_overview_sparklines <- function(workbook_path, n_values, label) {
  archive_entries <- unzip(workbook_path, list = TRUE)$Name
  image_entries <- archive_entries[grepl("^xl/media/image[0-9]+\\.png$", archive_entries)]
  image_entries <- image_entries[order(as.integer(sub(".*image([0-9]+)\\.png$", "\\1", image_entries)))]

  if (length(image_entries) < n_values) {
    stop("The ", label, " workbook does not contain all overview sparklines.")
  }

  image_dir <- file.path(tempdir(), paste0("season_summary_word_", label))
  dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)
  unzip(workbook_path, files = image_entries[seq_len(n_values)], exdir = image_dir)
  file.path(image_dir, image_entries[seq_len(n_values)])
}

format_overview_table <- function(overview) {
  overview <- overview[, c("Kategori", "Mnd. utvikling (antall)", "Antall", "Andel (%)")]
  overview[["Antall"]] <- ifelse(
    is.na(overview[["Antall"]]),
    "",
    formatC(as.numeric(overview[["Antall"]]), format = "d", big.mark = " ")
  )
  overview[["Andel (%)"]] <- ifelse(
    is.na(overview[["Andel (%)"]]),
    "",
    sprintf("%.1f%%", 100 * as.numeric(overview[["Andel (%)"]]))
  )
  overview[["Mnd. utvikling (antall)"]] <- ""
  overview
}

make_overview_table <- function(workbook_path, label) {
  overview <- read.xlsx(workbook_path, sheet = "Innledning", startRow = 3, check.names = FALSE)
  expected_columns <- c("Kategori", "Antall", "Andel (%)", "Mnd. utvikling (antall)")
  imported_columns <- c("Kategori", "Antall", "Andel.(%)", "Mnd..utvikling.(antall)")
  if (identical(names(overview), imported_columns)) {
    names(overview) <- expected_columns
  } else if (!identical(names(overview), expected_columns)) {
    stop("Unexpected overview-table columns in ", label, ".")
  }

  section_rows <- which(is.na(overview[["Antall"]]))
  value_rows <- which(!is.na(overview[["Antall"]]))
  sparkline_files <- extract_overview_sparklines(workbook_path, length(value_rows), label)
  display_table <- format_overview_table(overview)

  table_border <- fp_border(color = fhi_colour("blue_border"), width = 0.5, style = "dotted")
  ft <- flextable(display_table)
  ft <- font(ft, fontname = "Calibri", part = "all")
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- color(ft, color = fhi_colour("ink"), part = "all")
  ft <- bold(ft, part = "header")
  ft <- bg(ft, bg = fhi_colour("blue_grid"), part = "header")
  ft <- align(ft, align = "center", part = "header")
  ft <- align(ft, j = "Kategori", align = "left", part = "body")
  ft <- align(ft, j = c("Antall", "Andel (%)"), align = "right", part = "body")
  ft <- align(ft, j = "Mnd. utvikling (antall)", align = "center", part = "body")
  ft <- border_outer(ft, border = table_border)
  ft <- border_inner_h(ft, border = table_border)
  ft <- border_inner_v(ft, border = table_border)
  ft <- padding(ft, padding = 2, part = "all")
  ft <- height_all(ft, height = 0.25, part = "body")
  ft <- height(ft, i = 1, height = 0.25, part = "header")
  ft <- width(ft, j = "Kategori", width = 2.8)
  ft <- width(ft, j = "Antall", width = 0.75)
  ft <- width(ft, j = "Andel (%)", width = 0.85)
  ft <- width(ft, j = "Mnd. utvikling (antall)", width = 1.75)
  ft <- set_table_properties(ft, layout = "fixed")

  if (length(section_rows) > 0L) {
    ft <- bg(ft, i = section_rows, bg = fhi_colour("blue_grid"), part = "body")
    ft <- bold(ft, i = section_rows, bold = TRUE, part = "body")
  }

  for (idx in seq_along(value_rows)) {
    ft <- compose(
      ft,
      i = value_rows[idx],
      j = "Mnd. utvikling (antall)",
      value = as_paragraph(as_image(src = sparkline_files[idx], width = 1.65, height = 0.18))
    )
  }

  ft
}

flu_table <- make_overview_table(flu_source, "influenza")
sc2_table <- make_overview_table(sc2_source, "sc2")
rsv_table <- make_overview_table(rsv_source, "rsv")

heading_format <- fp_text(font.family = "Calibri", font.size = 12, bold = TRUE, color = fhi_colour("ink"))
title_format <- fp_text(font.family = "Calibri", font.size = 14, bold = TRUE, color = fhi_colour("ink"))
figure_format <- fp_text(font.family = "Calibri", font.size = 9, italic = TRUE, color = fhi_colour("ink"))

doc <- read_docx()
doc <- body_add_fpar(doc, fpar(ftext("Sesongoversikt 2025/26", prop = title_format)))
doc <- body_add_fpar(doc, fpar(ftext("Figurtekst: Trendlinjene viser antall sekvenserte virus per mnd. i sesongen 2025/26 (uke 35 til uke 34).", prop = figure_format)))
doc <- body_add_par(doc, "")
doc <- body_add_fpar(doc, fpar(ftext("Influensa", prop = heading_format)))
doc <- body_add_flextable(doc, value = flu_table)
doc <- body_add_par(doc, "")
doc <- body_add_fpar(doc, fpar(ftext("SARS-CoV-2", prop = heading_format)))
doc <- body_add_flextable(doc, value = sc2_table)
doc <- body_add_par(doc, "")
doc <- body_add_fpar(doc, fpar(ftext("RSV", prop = heading_format)))
doc <- body_add_flextable(doc, value = rsv_table)

print(doc, target = output_path)
message("Wrote ", output_path)
