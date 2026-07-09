# Shared PowerPoint export helpers used by pathogen-specific runners.

normalize_slide_title <- function(title_value, fallback_title) {
  if (is.null(title_value) || length(title_value) == 0) {
    return(fallback_title)
  }

  title_chr <- trimws(paste(as.character(title_value), collapse = " "))
  if (!nzchar(title_chr) || identical(title_chr, "NULL")) {
    return(fallback_title)
  }

  title_chr
}

build_slide_subtitle <- function(data, date_col = NULL, obs_n = NULL, details = NULL) {
  parts <- character()
  if (!is.null(date_col) && date_col %in% names(data)) {
    date_vals <- data[[date_col]]
    date_vals <- suppressWarnings(as.Date(date_vals))
    if (all(is.na(date_vals)) && is.character(data[[date_col]])) {
      parsed_month <- tryCatch(parse_month_key_nb(data[[date_col]]), error = function(e) rep(as.Date(NA), length(data[[date_col]])))
      if (!all(is.na(parsed_month))) {
        date_vals <- parsed_month
      }
    }
    date_vals <- date_vals[!is.na(date_vals)]
    if (length(date_vals) > 0) {
      parts <- c(parts, paste0("Periode: ", format(min(date_vals), "%b %Y"), "-", format(max(date_vals), "%b %Y")))
    }
  }
  if (is.null(obs_n)) {
    obs_n <- nrow(data)
  }
  if (!is.na(obs_n) && length(obs_n) == 1) {
    parts <- c(parts, paste0("Observasjoner: n=", format(obs_n, big.mark = ",")))
  }
  if (!is.null(details) && nzchar(details)) {
    parts <- c(parts, details)
  }
  paste(parts, collapse = " | ")
}

write_ppt_safe <- function(ppt_obj, target_path) {
  tryCatch({
    print(ppt_obj, target = target_path)
    target_path
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("is open", msg, fixed = TRUE)) {
      alt_path <- file.path(
        dirname(target_path),
        paste0(tools::file_path_sans_ext(basename(target_path)), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pptx")
      )
      warning("Target PPTX is open; writing to fallback file: ", alt_path)
      print(ppt_obj, target = alt_path)
      alt_path
    } else {
      stop(e)
    }
  })
}

add_table_output <- function(presentation, table_obj, slide_title = NULL, title = NULL) {
  if (is.null(table_obj)) {
    return(presentation)
  }

  resolved_title <- normalize_slide_title(
    if (!is.null(slide_title)) slide_title else title,
    "Tabell uten tittel"
  )
  save_table_to_ppt(presentation, table_obj, resolved_title)
}

add_plot_output <- function(presentation, plot_obj, slide_title = NULL, title = NULL) {
  if (is.null(plot_obj)) {
    return(presentation)
  }

  resolved_title <- slide_title
  if (is.null(resolved_title) || length(resolved_title) == 0) {
    resolved_title <- title
  }
  if (is.null(resolved_title) || length(resolved_title) == 0) {
    plot_label_title <- tryCatch(plot_obj$labels$title, error = function(e) NULL)
    resolved_title <- plot_label_title
  }
  resolved_title <- normalize_slide_title(resolved_title, "Figur uten tittel")

  save_plot_to_ppt(presentation, plot_obj, title = resolved_title)
}

add_section_output <- function(presentation, section_title, section_subtitle = NULL) {
  add_section_slide(presentation, section_title, section_subtitle)
}

extract_section_label <- function(section_title) {
  section_chr <- trimws(as.character(section_title))
  sub("^Seksjon:[[:space:]]*", "", section_chr)
}

add_report_title_slide <- function(presentation, report_title, report_subtitle = NULL, slide_title = NULL) {
  report_title <- normalize_norwegian_text(report_title)
  report_subtitle <- if (is.null(report_subtitle)) NULL else normalize_norwegian_text(report_subtitle)
  accent_color <- "#4f6b72"
  title_plot <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::scale_x_continuous(NULL, breaks = NULL) +
    ggplot2::scale_y_continuous(NULL, breaks = NULL) +
    ggplot2::annotate("segment", x = 0.24, xend = 0.76, y = 0.30, yend = 0.30, linewidth = 1.4, colour = accent_color) +
    ggplot2::annotate("text", x = 0.5, y = 0.62, label = report_title, size = 13, fontface = "bold", family = "sans") +
    ggplot2::annotate("text", x = 0.5, y = 0.44, label = ifelse(is.null(report_subtitle), "", report_subtitle), size = 8, family = "sans", colour = "#4a4a4a") +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  save_plot_to_ppt(
    presentation = presentation,
    plot = title_plot,
    title = if (is.null(slide_title)) report_title else normalize_norwegian_text(slide_title)
  )
}

add_report_index_slide <- function(presentation, section_titles, slide_title = "Rapportindeks", slide_subtitle = "Kort oversikt over hovedseksjonene i rapporten") {
  section_labels <- unique(vapply(section_titles, extract_section_label, character(1)))
  section_labels <- normalize_norwegian_text(section_labels[nzchar(section_labels)])
  if (length(section_labels) == 0) return(presentation)

  index_plot <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::scale_x_continuous(NULL, breaks = NULL) +
    ggplot2::scale_y_continuous(NULL, breaks = NULL) +
    ggplot2::annotate("text", x = 0.10, y = 0.92, label = normalize_norwegian_text(slide_subtitle), hjust = 0, vjust = 1, size = 4.2, family = "sans", colour = "#4a4a4a") +
    ggplot2::annotate("segment", x = 0.10, xend = 0.90, y = 0.84, yend = 0.84, linewidth = 0.9, colour = "#4f6b72") +
    ggplot2::annotate(
      "text",
      x = 0.10, y = 0.78,
      label = paste(paste0(seq_along(section_labels), ". ", section_labels), collapse = "\n"),
      hjust = 0, vjust = 1, size = 5.0, family = "sans"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  save_plot_to_ppt(
    presentation = presentation,
    plot = index_plot,
    title = normalize_norwegian_text(slide_title)
  )
}

