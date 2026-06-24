audit_ppt_lines <- c()

audit_log <- function(...) {
  line <- paste0(...)
  audit_ppt_lines <<- c(audit_ppt_lines, line)
  message(line)
}

extract_slide_text <- function(ppt_path, slide_name) {
  con <- unz(ppt_path, slide_name, open = "rb")
  on.exit(close(con), add = TRUE)
  xml <- readLines(con, warn = FALSE, encoding = "UTF-8")
  text <- paste(xml, collapse = "\n")
  hits <- stringr::str_match_all(text, "<a:t>(.*?)</a:t>")[[1]]
  if (nrow(hits) == 0) character(0) else hits[, 2]
}

extract_slide_title <- function(texts) {
  texts <- texts[nzchar(trimws(texts))]
  if (length(texts) == 0) return(NA_character_)
  texts[[length(texts)]]
}

slide_num_from_name <- function(x) {
  as.integer(sub("^ppt/slides/slide([0-9]+)\\.xml$", "\\1", x))
}

ppt_path_audit <- file.path(Sys.getenv("INF_RESULTS_DIR"), file_name_result)
slide_entries <- utils::unzip(ppt_path_audit, list = TRUE)
slide_entries <- slide_entries[grepl("^ppt/slides/slide[0-9]+\\.xml$", slide_entries$Name), , drop = FALSE]
slide_entries <- slide_entries[order(slide_num_from_name(slide_entries$Name)), , drop = FALSE]

slide_audit <- lapply(slide_entries$Name, function(slide_name) {
  texts <- extract_slide_text(ppt_path_audit, slide_name)
  title <- extract_slide_title(texts)
  bad_texts <- unique(texts[grepl("Ã|Â|�", texts)])
  pct_values <- stringr::str_extract_all(texts, "[0-9]+(?:[.,][0-9]+)?%")[[1]]
  pct_numeric <- suppressWarnings(as.numeric(sub("%$", "", gsub(",", ".", pct_values))))
  pct_numeric <- pct_numeric[!is.na(pct_numeric)]
  data.frame(
    slide = slide_num_from_name(slide_name),
    title = ifelse(is.na(title), "", title),
    has_mojibake = length(bad_texts) > 0,
    mojibake_text = if (length(bad_texts) > 0) paste(bad_texts, collapse = " | ") else "",
    max_pct_label = if (length(pct_numeric) > 0) max(pct_numeric) else NA_real_,
    stringsAsFactors = FALSE
  )
})
slide_audit <- dplyr::bind_rows(slide_audit)

mojibake_slides <- slide_audit %>% dplyr::filter(has_mojibake)
suspicious_pct_slides <- slide_audit %>% dplyr::filter(!is.na(max_pct_label), max_pct_label > 100)

audit_log("[AUDIT] PPT path: ", ppt_path_audit)
audit_log("[AUDIT] Slide count: ", nrow(slide_audit))
audit_log("[AUDIT] Slides with mojibake: ", nrow(mojibake_slides))
if (nrow(mojibake_slides) > 0) {
  apply(mojibake_slides, 1, function(row) {
    audit_log(
      "[AUDIT][MOJIBAKE] slide ", row[["slide"]],
      " title='", row[["title"]], "' text='", row[["mojibake_text"]], "'"
    )
  })
}

audit_log("[AUDIT] Slides with percent labels above 100: ", nrow(suspicious_pct_slides))
if (nrow(suspicious_pct_slides) > 0) {
  apply(suspicious_pct_slides, 1, function(row) {
    audit_log(
      "[AUDIT][PCT>100] slide ", row[["slide"]],
      " title='", row[["title"]], "' max_pct_label=", row[["max_pct_label"]]
    )
  })
}

subclade_color_col_audit <- intersect(c("nc_ha_subclade", "NC_HA_Subclade"), names(fludb))[1]
virus_col_audit <- if ("ngs_sekvens_resultat" %in% names(fludb)) "ngs_sekvens_resultat" else NULL
virus_map_audit <- c("A/H1N1" = "H1N1", "A/H3N2" = "H3N2", "B/Victoria" = "BVIC")

if (!is.na(subclade_color_col_audit) && !is.null(virus_col_audit) && "prove_tatt" %in% names(fludb)) {
  for (virus_key_audit in names(virus_map_audit)) {
    virus_label_audit <- virus_map_audit[[virus_key_audit]]
    flu_v_audit <- fludb %>%
      dplyr::filter(
        season == current_season_label,
        .data[[virus_col_audit]] == virus_key_audit,
        !is.na(.data[[subclade_color_col_audit]]),
        trimws(as.character(.data[[subclade_color_col_audit]])) != ""
      )
    landsdel_col_audit <- if ("pasient_landsdel" %in% names(flu_v_audit)) {
      "pasient_landsdel"
    } else if ("pasient_landsdel_from_fylke" %in% names(flu_v_audit)) {
      "pasient_landsdel_from_fylke"
    } else {
      NULL
    }
    if (is.null(landsdel_col_audit) || nrow(flu_v_audit) == 0) next

    plots_audit <- build_subclade_landsdel_month_bars(
      flu_v_audit,
      subtype_name = virus_key_audit,
      date_col = "prove_tatt",
      subtype_col = virus_col_audit,
      landsdel_col = landsdel_col_audit,
      subclade_col = subclade_color_col_audit,
      title_prefix = "Pasient landsdel per måned og subklade",
      palette_base = kvalitativ_comb,
      top_n_subclades = 3,
      facet_ncol = 2
    )

    pct_sums <- plots_audit$percent_plot$data %>%
      dplyr::group_by(month_date, landsdel_label) %>%
      dplyr::summarise(total_percent = sum(percent), .groups = "drop")
    max_total_percent <- max(pct_sums$total_percent, na.rm = TRUE)
    over_100_n <- sum(pct_sums$total_percent > 100.0001, na.rm = TRUE)

    audit_log(
      "[AUDIT][SUBTYPE] ", virus_label_audit,
      " max_total_percent=", sprintf("%.4f", max_total_percent),
      " groups_over_100=", over_100_n
    )
  }
}

audit_file <- file.path(
  Sys.getenv("INF_RESULTS_DIR"),
  sub("\\.pptx$", "_ppt_audit.txt", file_name_result)
)
writeLines(audit_ppt_lines, audit_file, useBytes = TRUE)
message("[AUDIT] Summary written: ", audit_file)
