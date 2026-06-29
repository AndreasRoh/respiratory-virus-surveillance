# Combined single-file SC2 analysis generated from SC2_Analysis.R and inlined SC2 module scripts.
# Keeps external sources: common_report_utils.R, SC2 SQL query scripts, and SC2_Classification.R.

# ============================================================================
# SARS-COV2 Surveillance
# Author: AR
# Date: 29.04.2026
# ============================================================================

resolve_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
  }
  this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE), error = function(e) "")
  if (nzchar(this_file)) {
    return(dirname(this_file))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
bundle_scripts_dir <- resolve_script_dir()
default_results_root <- "N:/Virologi/Influensa/2526/WGS_Analyse/Results"
results_root <- Sys.getenv("SC2_RESULTS_DIR", unset = default_results_root)
results_stats_dir <- Sys.getenv("SC2_RESULTS_STATS_DIR", unset = file.path(results_root, "Statistikk"))
results_share_dir <- Sys.getenv(
  "SC2_RESULTS_SHARE_DIR",
  unset = "C:/Users/aroh/OneDrive - Folkehelseinstituttet/Sesong 2025_26"
)

# ============================================================================
# SETUP - Install and load required packages
# ============================================================================

# Function to check, update, and install packages
check_install_update_packages <- function(packages) {
  installed_pkgs <- rownames(installed.packages())
  missing_packages <- setdiff(packages, installed_pkgs)
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }

  # Avoid in-session package updates: they can fail on Windows due to locked DLLs
  # and can trigger namespace unload conflicts (e.g., DBI imported by odbc/RSQLite).
  if (isTRUE(getOption("sc2_update_packages", FALSE))) {
    outdated_packages <- old.packages()
    if (!is.null(outdated_packages)) {
      outdated_required <- intersect(rownames(outdated_packages), packages)
      if (length(outdated_required) > 0) {
        message("Updating outdated required packages: ", paste(outdated_required, collapse = ", "))
        update.packages(oldPkgs = outdated_required, ask = FALSE, checkBuilt = TRUE)
      }
    }
  }
}

# List of required packages
required_packages <- c(
  "odbc", "RSQLite", "DBI", "tidyverse", "lubridate", "ggrepel", "scales",
  "openxlsx", "RColorBrewer", "officer", "tsibble", "patchwork",
  "flextable", "cowplot", "zoo", "reshape2", "janitor", "rvg", "treemapify",
  "data.table", "tools", "knitr"
)

# Use CRAN mirror in Germany for package installs/updates.
options(repos = c(CRAN = "https://cran.uni-muenster.de/"))

# Show warnings immediately and keep error messages visible in terminal.
options(warn = 1, show.error.messages = TRUE)

# Locale guard: set UTF-8 Norwegian locale as early as possible for stable
# rendering of Norwegian characters in plots and slide titles.
init_locale <- function() {
  locale_candidates <- c("nb_NO.UTF-8", "Norwegian (Bokmal)_Norway.utf8", "Norwegian")
  for (loc in locale_candidates) {
    ok <- tryCatch(Sys.setlocale(category = "LC_ALL", locale = loc), error = function(e) NA_character_)
    if (!is.na(ok) && nzchar(ok)) {
      message("Locale set to: ", ok)
      return(invisible(ok))
    }
  }
  warning("Could not set Norwegian UTF-8 locale; continuing with system default locale.")
  invisible(NA_character_)
}
init_locale()

# Restore terminal output if previous runs left active sinks.
while (sink.number() > 0) sink()
while (sink.number(type = "message") != 2) sink(type = "message")

analysis_started_at <- Sys.time()
log_timed_message <- function(...) {
  message(sprintf("[%s]", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), " ", paste0(..., collapse = ""))
  flush.console()
}
timed_step <- function(step_name, expr) {
  step_started_at <- Sys.time()
  log_timed_message("START: ", step_name)
  result <- force(expr)
  step_elapsed <- as.numeric(difftime(Sys.time(), step_started_at, units = "secs"))
  log_timed_message("DONE: ", step_name, " (", sprintf("%.2f", step_elapsed), "s)")
  result
}

load_required_libraries <- function(packages) {
  lapply(packages, function(pkg) {
    withCallingHandlers(
      library(pkg, character.only = TRUE),
      warning = function(w) {
        if (grepl("was built under R version", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  })
}

# Call the function to install and update packages
invisible(timed_step("Package install/update", check_install_update_packages(required_packages)))

# Load libraries
invisible(timed_step("Load libraries", load_required_libraries(required_packages)))
invisible(utils::globalVariables(c(
  ".", ".data", "Tessy_plot", "plot_date", "age_value", "age_group_plot",
  "age_group_raw", "prove_kategori_group", "group_plot", "n", "n_raw",
  "tessy_n", "percent", "kvalitativ_a", "sc2_palette"
)))

# Set data.table week option to legacy mode to maintain current behavior
options(datatable.week = "legacy")

# Execute shared report utilities (includes FHI palettes).
invisible(timed_step("Source common report utilities", source("Source_files/common_report_utils.R")))
if (!exists("fhi_discrete_palette", mode = "function")) {
  fhi_discrete_palette <- function(n, palette_name = NULL) {
    base_palette <- if (exists("kvalitativ_comb", inherits = TRUE)) {
      kvalitativ_comb
    } else if (exists("kvalitativ_a", inherits = TRUE)) {
      kvalitativ_a
    } else {
      c("#ec7c73", "#40436d", "#61d2b2", "#a93c38", "#f9dc8c", "#7176c9")
    }
    rep_len(base_palette, n)
  }
  warning("fhi_discrete_palette() not found after sourcing Color palettes.R; using kvalitativ palette fallback.")
}

# ============================================================================
# VARIABLES - Define variables and helper functions
# ============================================================================

# Define the mutations of interest in single or combination heatmap (Spike only)
mutations <- c("Q493E", "F456L", "V1104L", "R346T", "S31-")

start_date <- yearweek(Sys.Date() - weeks(8)) # 4 week period start
end_date <- yearweek(Sys.Date()) # 4 week period end

set_flextable_defaults(font.size = 6)

# Dynamic season boundaries (week 35 -> week 34) based on today's date.
season_info <- current_and_previous_seasons(Sys.Date())
current_season_label <- season_info$current_label
previous_season_label <- season_info$previous_label
current_season_bounds <- season_window_bounds(season_info$current_start_year)
data_window_start <- min(current_season_bounds$start, Sys.Date() %m-% months(6))

Seqlim <- 10 # How many sequences need to be valid per week to include in the analysis
export_graph <- read_pptx() # power point placeholder for the results

slide_plot_theme <- function() {
  theme(
    plot.margin = ggplot2::margin(3, 5, 3, 3),
    axis.text = element_text(size = 8.5),
    axis.title = element_text(size = 9.5),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.8),
    legend.key.height = grid::unit(0.14, "in"),
    legend.key.width = grid::unit(0.16, "in"),
    legend.spacing.y = grid::unit(0.01, "in"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.just = "center",
    panel.spacing = grid::unit(0.10, "in")
  )
}

apply_slide_plot_theme <- function(content) {
  if (inherits(content, "patchwork")) {
    return(content & slide_plot_theme())
  }
  if (inherits(content, "ggplot")) {
    return(content + slide_plot_theme())
  }
  content
}

collapse_minor_categories <- function(data, category_col, value_col, other_label = "Andre", min_peak_percent = 5, top_n = 8, preserve_categories = NULL) {
  category_sym <- rlang::sym(category_col)
  value_sym <- rlang::sym(value_col)

  keep_categories <- data %>%
    filter(!is.na(!!category_sym), trimws(as.character(!!category_sym)) != "") %>%
    group_by(!!category_sym) %>%
    summarise(
      max_value = suppressWarnings(max(!!value_sym, na.rm = TRUE)),
      total_value = sum(!!value_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      max_value = ifelse(is.finite(max_value), max_value, 0),
      total_value = ifelse(is.finite(total_value), total_value, 0)
    ) %>%
    arrange(desc(max_value), desc(total_value), as.character(!!category_sym)) %>%
    mutate(rank_order = row_number()) %>%
    filter(
      max_value > min_peak_percent |
        rank_order <= top_n |
        as.character(!!category_sym) %in% preserve_categories
    ) %>%
    pull(!!category_sym) %>%
    as.character() %>%
    unique()

  data %>%
    mutate(
      !!paste0(category_col, "_plot") := ifelse(
        as.character(!!category_sym) %in% keep_categories,
        as.character(!!category_sym),
        other_label
      )
    )
}
resolve_layout_placeholder <- function(presentation, layout, master, placeholder_type) {
  layout_df <- officer::layout_properties(presentation, layout = layout, master = master)
  match_row <- layout_df %>%
    filter(type == placeholder_type) %>%
    slice(1)

  if (nrow(match_row) == 0) {
    return(NULL)
  }

  match_row
}

placeholder_to_location <- function(placeholder_row, pad_left = 0, pad_top = 0, pad_right = 0, pad_bottom = 0) {
  officer::ph_location(
    left = placeholder_row$offx + pad_left,
    top = placeholder_row$offy + pad_top,
    width = max(0.2, placeholder_row$cx - pad_left - pad_right),
    height = max(0.2, placeholder_row$cy - pad_top - pad_bottom)
  )
}

# Unified PowerPoint export helper for both ggplot objects and tables.
export_to_ppt <- function(presentation, content, slide_title, slide_subtitle = NULL, layout = "Title and Content", master = "Office Theme", preserve_plot_labels = FALSE) {
  sanitize_xml_text <- function(x) {
    x <- as.character(x)
    if (exists("normalize_norwegian_text", mode = "function")) {
      x <- normalize_norwegian_text(x)
    }
    gsub("[\\x00-\\x08\\x0B\\x0C\\x0E-\\x1F]", "", x, perl = TRUE)
  }

  slide_title <- sanitize_xml_text(slide_title)
  slide_subtitle <- if (is.null(slide_subtitle)) NULL else sanitize_xml_text(slide_subtitle)
  if (length(layout) != 1) layout <- layout[1]

  title_location <- officer::ph_location_type(type = "title")
  body_left <- 0.50
  body_top <- 1.75
  body_width <- 9.00
  body_height <- 4.95
  body_location <- officer::ph_location(left = body_left, top = body_top, width = body_width, height = body_height)


  title_value <- if (is.null(slide_subtitle)) {
    officer::fpar(
      officer::ftext(
        slide_title,
        officer::fp_text(font.size = 22, bold = TRUE, font.family = "Arial", color = "#000000")
      )
    )
  } else {
    officer::block_list(
      officer::fpar(
        officer::ftext(
          slide_title,
          officer::fp_text(font.size = 22, bold = TRUE, font.family = "Arial", color = "#000000")
        )
      ),
      officer::fpar(
        officer::ftext(
          slide_subtitle,
          officer::fp_text(font.size = 10.5, bold = FALSE, font.family = "Arial", color = "#444444")
        )
      )
    )
  }

  presentation <- officer::add_slide(presentation, layout = layout, master = master)
  presentation <- officer::ph_with(presentation, value = title_value, location = title_location)

  if (inherits(content, "ggplot") || inherits(content, "patchwork")) {
    if (exists("normalize_object_text", mode = "function")) {
      content <- normalize_object_text(content)
    }
    content <- apply_slide_plot_theme(content)
    if (!isTRUE(preserve_plot_labels)) {
      content <- content + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
    }
    presentation <- officer::ph_with(
      presentation,
      value = rvg::dml(code = print(content)),
      location = body_location
    )
  } else {
    if (inherits(content, "flextable")) {
      ft <- content
    } else if (is.data.frame(content)) {
      if (exists("normalize_object_text", mode = "function")) {
        content <- normalize_object_text(content)
      }
      char_cols <- vapply(content, is.character, logical(1))
      content[char_cols] <- lapply(content[char_cols], sanitize_xml_text)
      ft <- flextable::flextable(as.data.frame(content)) |> flextable::autofit()
    } else if (is.character(content)) {
      if (exists("normalize_norwegian_text", mode = "function")) {
        content <- normalize_norwegian_text(content)
      }
      content <- sanitize_xml_text(content)
      ft <- flextable::flextable(as.data.frame(content)) |> flextable::autofit()
    } else {
      ft <- flextable::flextable(as.data.frame(content)) |> flextable::autofit()
    }
    presentation <- officer::ph_with(presentation, value = ft, location = body_location)
  }

  presentation
}

# Backward-compatible wrappers that route through the single export function.
save_plot <- function(plot, slide_title, export_graph) {
  export_to_ppt(export_graph, plot, slide_title)
}

build_slide_subtitle <- function(data, date_col = NULL, obs_n = NULL, details = NULL) {
  parts <- character()
  if (!is.null(date_col) && date_col %in% names(data)) {
    date_vals <- as.Date(data[[date_col]])
    date_vals <- date_vals[!is.na(date_vals)]
    if (length(date_vals) > 0) {
      parts <- c(parts, paste0("Periode: ", format_month_year_nb(min(date_vals)), " - ", format_month_year_nb(max(date_vals))))
    }
  }
  if (is.null(obs_n)) {
    obs_n <- nrow(data)
  }
  if (!is.na(obs_n) && length(obs_n) == 1) {
    parts <- c(parts, paste0("Observasjoner: n=", scales::label_number(big.mark = " ", decimal.mark = ",", accuracy = 1)(obs_n)))
  }
  if (!is.null(details) && nzchar(details)) {
    parts <- c(parts, details)
  }
  paste(parts, collapse = " | ")
}

make_monthly_heatmap <- function(
  data,
  date_col,
  category_col,
  value_col,
  fill_label,
  x_breaks = "1 month",
  bottom_categories = c("Andre undervarianter", "Andre SARS CoV 2", "Andre", "Ukjent"),
  base_y_text_size = 8
) {
  plot_data <- data %>%
    transmute(
      date_value = as.Date(.data[[date_col]]),
      category_value = as.character(.data[[category_col]]),
      fill_value = as.numeric(.data[[value_col]])
    ) %>%
    filter(!is.na(date_value), !is.na(category_value), nzchar(category_value), is.finite(fill_value)) %>%
    group_by(date_value, category_value) %>%
    summarise(fill_value = sum(fill_value, na.rm = TRUE), .groups = "drop")

  if (nrow(plot_data) == 0) {
    return(ggplot() + theme_void())
  }

  month_values <- sort(unique(plot_data$date_value))
  category_values <- plot_data %>%
    group_by(category_value) %>%
    summarise(total_value = sum(fill_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_value), category_value) %>%
    pull(category_value) %>%
    move_named_levels_to_end(bottom_categories)

  y_text_size <- heatmap_axis_text_size(
    length(category_values),
    base_size = base_y_text_size,
    min_size = max(5, base_y_text_size - 3)
  )

  plot_grid <- tidyr::expand_grid(
    date_value = month_values,
    category_value = category_values
  ) %>%
    left_join(plot_data, by = c("date_value", "category_value")) %>%
    mutate(
      fill_value = dplyr::coalesce(fill_value, 0),
      category_value = factor(category_value, levels = rev(category_values))
    )

  ggplot(plot_grid, aes(x = date_value, y = category_value, fill = fill_value)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_x_date(labels = format_month_label, date_breaks = x_breaks) +
    scale_fill_gradient(low = "#F7FBFF", high = "#0B5FA5") +
    labs(x = "M\u00e5ned", y = NULL, fill = fill_label) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = y_text_size),
      legend.position = "right"
    )
}

normalize_percent_values <- function(x) {
  x_num <- as.numeric(x)
  finite_x <- x_num[is.finite(x_num)]
  if (length(finite_x) == 0) {
    return(x_num)
  }
  if (max(finite_x, na.rm = TRUE) <= 1.0001) {
    x_num * 100
  } else {
    x_num
  }
}

make_percent_count_plot <- function(
  data,
  x_col,
  y_percent_col,
  fill_col,
  count_data,
  count_y_col,
  fill_values,
  fill_label,
  x_breaks = "1 month",
  legend_nrow = 2,
  heatmap_bottom_categories = c("Andre undervarianter", "Andre SARS CoV 2", "Andre", "Ukjent"),
  heatmap_y_text_size = 8
) {
  top_plot <- make_monthly_heatmap(
    data = data,
    date_col = x_col,
    category_col = fill_col,
    value_col = y_percent_col,
    fill_label = "Andel (%)",
    x_breaks = x_breaks,
    bottom_categories = heatmap_bottom_categories,
    base_y_text_size = heatmap_y_text_size
  ) +
    theme(axis.text.x = element_blank())

  bottom_plot <- ggplot(count_data, aes(x = .data[[x_col]], y = .data[[count_y_col]])) +
    geom_col(fill = kvalitativ_b[2]) +
    labs(x = "M\u00e5ned", y = "Antall (n)") +
    scale_x_date(labels = format_month_label, date_breaks = x_breaks) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "none"
    )

  top_plot / bottom_plot + patchwork::plot_layout(heights = c(3, 1))
}

make_stacked_percent_count_plot <- function(
  data,
  date_col,
  percent_col,
  count_col,
  category_col,
  fill_values,
  fill_label,
  x_breaks = "1 month",
  legend_ncol = 2,
  legend_position = "right"
) {
  plot_data <- data %>%
    mutate(
      percent_plot = normalize_percent_values(.data[[percent_col]]),
      category_plot = as.character(.data[[category_col]])
    ) %>%
    filter(!is.na(.data[[date_col]]), !is.na(category_plot), trimws(category_plot) != "")

  count_df <- plot_data %>%
    transmute(date_value = as.Date(.data[[date_col]]), count_value = as.numeric(.data[[count_col]])) %>%
    distinct() %>%
    arrange(date_value)

  top_plot <- ggplot(plot_data, aes(x = .data[[date_col]], y = percent_plot, fill = category_plot)) +
    geom_col() +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    labs(x = NULL, y = "Andel (%)", fill = fill_label) +
    scale_x_date(labels = format_month_label, date_breaks = x_breaks) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      legend.position = legend_position,
      legend.box = "vertical"
    ) +
    guides(fill = guide_legend(ncol = min(legend_ncol, 2), byrow = TRUE))

  bottom_plot <- ggplot(count_df, aes(x = date_value, y = count_value)) +
    geom_col(fill = "grey75", color = "black", alpha = 0.5) +
    labs(x = "M\u00e5ned", y = "Antall (n)") +
    scale_x_date(labels = format_month_label, date_breaks = x_breaks) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "none"
    )

  (top_plot / bottom_plot + patchwork::plot_layout(guides = "collect", heights = c(3, 1.5))) &
    theme(
      legend.position = legend_position,
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "center"
    )
}

make_percent_trend_plot <- function(
  data,
  date_col,
  percent_col,
  count_col,
  category_col,
  color_values,
  color_label,
  x_breaks = "1 month",
  legend_nrow = 2,
  legend_position = "right",
  legend_ncol = 2
) {
  plot_data <- data %>%
    mutate(
      percent_plot = normalize_percent_values(.data[[percent_col]]),
      category_plot = as.character(.data[[category_col]])
    ) %>%
    filter(!is.na(.data[[date_col]]), !is.na(category_plot), trimws(category_plot) != "")

  count_df <- plot_data %>%
    transmute(date_value = as.Date(.data[[date_col]]), count_value = as.numeric(.data[[count_col]])) %>%
    distinct() %>%
    arrange(date_value)

  max_count <- suppressWarnings(max(count_df$count_value, na.rm = TRUE))
  if (!is.finite(max_count) || max_count <= 0) {
    max_count <- 1
  }

  top_plot <- ggplot(plot_data, aes(x = .data[[date_col]], y = percent_plot, fill = category_plot)) +
    geom_col() +
    scale_fill_manual(values = color_values, drop = FALSE) +
    labs(x = NULL, y = "Andel (%)", fill = color_label) +
    scale_x_date(labels = format_month_label, date_breaks = x_breaks) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      legend.position = legend_position,
      legend.box = "vertical"
    ) +
    guides(fill = guide_legend(ncol = min(legend_ncol, 2), byrow = TRUE))

  bottom_plot <- ggplot() +
    geom_col(
      data = count_df,
      aes(x = date_value, y = count_value / max_count * 100, fill = "Antall (n)"),
      alpha = 0.35,
      color = "black",
      position = "identity"
    ) +
    geom_line(
      data = plot_data,
      aes(x = .data[[date_col]], y = percent_plot, color = category_plot, group = category_plot),
      linewidth = 1.2
    ) +
    scale_y_continuous(
      name = "Andel (%)",
      sec.axis = sec_axis(~ . * max_count / 100, name = "Antall (n)")
    ) +
    scale_color_manual(values = color_values, name = color_label, drop = FALSE) +
    scale_fill_manual(values = c("Antall (n)" = "grey75"), name = NULL) +
    scale_x_date(labels = format_month_label, date_breaks = x_breaks) +
    labs(x = "M\u00e5ned") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = legend_position,
      legend.box = "vertical"
    ) +
    guides(
      color = guide_legend(ncol = min(legend_ncol, 2), byrow = TRUE),
      fill = guide_legend(ncol = 1, byrow = FALSE)
    )

  (top_plot / bottom_plot + patchwork::plot_layout(guides = "collect", heights = c(3, 2))) &
    theme(
      legend.position = legend_position,
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "center"
    )
}

add_section_slide <- function(presentation, section_title, section_subtitle = NULL, main_topics = NULL) {
  topic_text <- if (!is.null(main_topics) && length(main_topics) > 0) {
    paste(paste0("- ", main_topics), collapse = "\n")
  } else {
    ""
  }
  section_plot <- ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    scale_x_continuous(NULL, breaks = NULL) +
    scale_y_continuous(NULL, breaks = NULL) +
    annotate("text", x = 0.5, y = 0.68, label = section_title, size = 11, fontface = "bold", family = "sans") +
    annotate("text", x = 0.5, y = 0.50, label = ifelse(is.null(section_subtitle), "", section_subtitle), size = 5, family = "sans") +
    annotate("text", x = 0.18, y = 0.27, label = topic_text, hjust = 0, vjust = 1, size = 4.4, family = "sans") +
    theme_void() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  export_to_ppt(presentation, section_plot, section_title, section_subtitle)
}

current_week_title <- week(Sys.Date())
current_year_title <- year(Sys.Date())
title_plot <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  scale_x_continuous(NULL, breaks = NULL) +
  scale_y_continuous(NULL, breaks = NULL) +
  annotate("text", x = 0.5, y = 0.60, label = "SARS-CoV-2-overv\u00E5king", size = 13, fontface = "bold", family = "sans") +
  annotate("text", x = 0.5, y = 0.42, label = paste0("Uke ", current_week_title, " - ", current_year_title), size = 8, family = "sans") +
  theme_void() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

export_graph <- export_to_ppt(export_graph, title_plot, paste0("SARS-CoV-2 Uke ", current_week_title))
index_topics <- c(
  "Datakvalitet og runkvalitet",
  "Sekvensering og dekning",
  "Pangolin og Tessy",
  "Mutasjoner og indeler",
  "PCR CT og pasientdata",
  "Fylker og resistens"
)
index_plot <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  scale_x_continuous(NULL, breaks = NULL) +
  scale_y_continuous(NULL, breaks = NULL) +
  annotate("text", x = 0.10, y = 0.88, label = paste(paste0(seq_along(index_topics), ". ", index_topics), collapse = "\n"), hjust = 0, vjust = 1, size = 5.2, family = "sans") +
  theme_void() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )
export_graph <- export_to_ppt(
  export_graph,
  index_plot,
  "Rapportindeks",
  "Kort oversikt over hovedseksjonene i ukesrapporten"
)
export_graph <- add_section_slide(
  export_graph,
  "Datakvalitet",
  "Dataintegritet og kompletthet i hele datasettet",
  c("Kvalitetsm\u00e5linger", "Avvik og outliers")
)


# ============================================================================
# DATA LOAD - SQL (BN)
# ============================================================================

invisible(timed_step("Source SC2_SQLquery_BNCOVID19.R", source(file.path(bundle_scripts_dir, "SC2_SQLquery_BNCOVID19.R"))))

invisible(timed_step("Source SC2_DataCleaning_BNCOVID19.R", source(file.path(bundle_scripts_dir, "SC2_DataCleaning_BNCOVID19.R"))))
invisible(timed_step("Source SC2_SQLquery_25-26.R", source(file.path(bundle_scripts_dir, "SC2_SQLquery_25-26.R"))))

invisible(timed_step("Source SC2_DataCleaning_25-26.R", source(file.path(bundle_scripts_dir, "SC2_DataCleaning_25-26.R"))))

# Root-cause fix for DBI warning:
# SQL address scripts can leave live connection objects in the global environment.
# If they are garbage-collected later, DBI warns: "call dbDisconnect() when finished working with a connection".
# We proactively disconnect all known SC2 connection handles before classification.
disconnect_if_valid <- function(obj_name, env = .GlobalEnv) {
  if (!exists(obj_name, envir = env, inherits = FALSE)) {
    return(invisible(NULL))
  }
  obj <- get(obj_name, envir = env, inherits = FALSE)
  if (inherits(obj, "DBIConnection")) {
    tryCatch(
      {
        if (DBI::dbIsValid(obj)) DBI::dbDisconnect(obj)
      },
      error = function(e) NULL
    )
  }
  rm(list = obj_name, envir = env)
  invisible(NULL)
}
invisible(lapply(c("conBNCOVID", "conSC22526", "conSC2", "con"), disconnect_if_valid))
closeAllConnections()

invisible(timed_step("Source SC2_Classification.R", source(file.path(bundle_scripts_dir, "SC2_Classification.R"))))

# Harmonized sample-category classification:
# P1/P1_* or 1* -> Sentinel; everything else -> Non-Sentinel.
classify_prove_kategori_group <- function(x) {
  x_chr <- trimws(as.character(x))
  dplyr::case_when(
    grepl("^(P1\\b|P1_|1\\b)", x_chr) ~ "Sentinel",
    TRUE ~ "Non-Sentinel"
  )
}
clean_project_code <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_chr <- gsub("\\s+", "", x_chr)
  x_chr <- ifelse(grepl("^[0-9]+", x_chr), sub("^([0-9]+).*$", "P\\1", x_chr), x_chr)
  x_chr <- ifelse(grepl("^P[0-9]+", x_chr), sub("^(P[0-9]+).*$", "\\1", x_chr), NA_character_)
  x_chr
}

# -----------------------------------------------------------------------------
# Post-cleaning harmonization for report use.
# Text cleanup belongs in the sourced cleaning files; keep only report-specific
# date and project-group enrichment here.
# -----------------------------------------------------------------------------
if (exists("SC2db")) {
  SC2db <- SC2db %>%
    mutate(
      prove_tatt = as.Date(prove_tatt),
      season = season_label_from_date(prove_tatt),
      prove_kategori_group = classify_prove_kategori_group(prove_kategori),
      prove_project_clean = ifelse(
        prove_kategori_group == "Non-Sentinel",
        clean_project_code(prove_kategori),
        NA_character_
      )
    )
}

if (exists("SC2db_v")) {
  SC2db_v <- SC2db_v %>%
    mutate(
      prove_tatt = as.Date(prove_tatt),
      season = season_label_from_date(prove_tatt),
      prove_kategori_group = classify_prove_kategori_group(prove_kategori),
      prove_project_clean = ifelse(
        prove_kategori_group == "Non-Sentinel",
        clean_project_code(prove_kategori),
        NA_character_
      )
    )
}



# ============================================================================
# SC2 DATA QUALITY AND EARLY EXPLORATION
# Output in this section:
#   - QC summary tables
#   - Missingness profile
#   - Numeric outlier scan
# ============================================================================

eda_source_df <- if (exists("SC2db")) SC2db else SC2db
results_dir_stats <- results_stats_dir
current_week_eda <- week(Sys.Date())
current_year_eda <- year(Sys.Date())

eda_tessy_col <- intersect(c("Tessy", "tessy"), names(eda_source_df))[1]
eda_date_col <- intersect(c("prove_tatt", "PROVE_TATT", "sample_date", "Sampledate"), names(eda_source_df))[1]
eda_cov_col <- intersect(c("nc_coverage", "coverage_breadth_artic", "coverage_breadth_swift", "coverage_breadth_eksterne", "coverage_breadth_nano"), names(eda_source_df))[1]
eda_key_col <- intersect(c("key", "KEY"), names(eda_source_df))[1]
eda_origin_col <- intersect(c("Origin", "origin"), names(eda_source_df))[1]

if (!is.na(eda_tessy_col)) {
  has_pk_group_col <- "prove_kategori_group" %in% names(eda_source_df)
  has_prove_kat_col <- "prove_kategori" %in% names(eda_source_df)
  has_pas_status_col <- "pasient_status" %in% names(eda_source_df)
  has_pas_vaks_col <- "pasient_vaks" %in% names(eda_source_df)
  has_pas_vaks_2u_col <- "pasient_vaks_2uipt" %in% names(eda_source_df)

  eda_df <- eda_source_df %>%
    mutate(
      plot_date_window = if (!is.na(eda_date_col)) as.Date(.data[[eda_date_col]]) else as.Date(NA),
      Tessy_plot = as.character(.data[[eda_tessy_col]]),
      prove_kategori_raw = if (has_prove_kat_col) as.character(prove_kategori) else NA_character_,
      prove_kategori_group = ifelse(
        if (has_pk_group_col) is.na(prove_kategori_group) else TRUE,
        classify_prove_kategori_group(prove_kategori_raw),
        if (has_pk_group_col) as.character(prove_kategori_group) else "Non-Sentinel"
      ),
      pasient_status_plot = if (has_pas_status_col) ifelse(is.na(pasient_status) | pasient_status == "" | pasient_status == "IKKE_SATT", "Ukjent", as.character(pasient_status)) else "Ukjent",
      pasient_vaks_plot = if (has_pas_vaks_col) ifelse(is.na(pasient_vaks) | trimws(as.character(pasient_vaks)) == "", "Ukjent", as.character(pasient_vaks)) else "Ukjent",
      pasient_vaks_2uipt_plot = if (has_pas_vaks_2u_col) ifelse(is.na(pasient_vaks_2uipt) | trimws(as.character(pasient_vaks_2uipt)) == "", "Ukjent", as.character(pasient_vaks_2uipt)) else "Ukjent"
    ) %>%
    filter(!is.na(plot_date_window), plot_date_window >= data_window_start)

  column_profile <- data.frame(
    column_name = names(eda_df),
    class = vapply(eda_df, function(x) paste(class(x), collapse = ","), character(1)),
    non_missing_n = vapply(eda_df, function(x) sum(!is.na(x) & trimws(as.character(x)) != ""), numeric(1)),
    unique_n = vapply(eda_df, function(x) dplyr::n_distinct(x, na.rm = TRUE), numeric(1))
  ) %>%
    mutate(
      missing_n = nrow(eda_df) - non_missing_n,
      missing_pct = round((missing_n / nrow(eda_df)) * 100, 2)
    ) %>%
    arrange(desc(missing_pct))

  # Data quality checks and outliers.
  qc_issues <- list()

  if (!is.na(eda_key_col)) {
    dup_key_df <- eda_df %>%
      filter(!is.na(.data[[eda_key_col]]), trimws(as.character(.data[[eda_key_col]])) != "") %>%
      count(.data[[eda_key_col]], name = "n") %>%
      filter(n > 1) %>%
      arrange(desc(n))
    qc_issues[["duplicate_keys"]] <- nrow(dup_key_df)
    if (nrow(dup_key_df) > 0) {
      export_graph <- export_to_ppt(export_graph, head(dup_key_df, 30), "SC2 duplicate key candidates")
    }
  }

  if (!is.na(eda_date_col)) {
    date_qc_df <- eda_df %>%
      mutate(plot_date = as.Date(.data[[eda_date_col]])) %>%
      summarise(
        missing_date_n = sum(is.na(plot_date)),
        future_date_n = sum(!is.na(plot_date) & plot_date > Sys.Date()),
        pre_2020_n = sum(!is.na(plot_date) & plot_date < as.Date("2020-01-01"))
      )
    qc_issues[["missing_date_n"]] <- date_qc_df$missing_date_n
    qc_issues[["future_date_n"]] <- date_qc_df$future_date_n
    qc_issues[["pre_2020_n"]] <- date_qc_df$pre_2020_n
  }

  age_qc_col <- intersect(c("pasient_alder", "pasient_age", "age", "alder"), names(eda_df))[1]
  if (!is.na(age_qc_col)) {
    age_qc_df <- eda_df %>%
      mutate(age_value = suppressWarnings(as.numeric(as.character(.data[[age_qc_col]])))) %>%
      summarise(
        missing_age_n = sum(is.na(age_value)),
        age_lt_0_n = sum(!is.na(age_value) & age_value < 0),
        age_gt_110_n = sum(!is.na(age_value) & age_value > 110)
      )
    qc_issues[["missing_age_n"]] <- age_qc_df$missing_age_n
    qc_issues[["age_lt_0_n"]] <- age_qc_df$age_lt_0_n
    qc_issues[["age_gt_110_n"]] <- age_qc_df$age_gt_110_n
  }

  if (!is.na(eda_cov_col)) {
    cov_qc_df <- eda_df %>%
      mutate(coverage_value = suppressWarnings(as.numeric(as.character(.data[[eda_cov_col]])))) %>%
      summarise(
        missing_cov_n = sum(is.na(coverage_value)),
        cov_lt_0_n = sum(!is.na(coverage_value) & coverage_value < 0),
        cov_gt_100_n = sum(!is.na(coverage_value) & coverage_value > 100),
        cov_between_1_100_n = sum(!is.na(coverage_value) & coverage_value > 1 & coverage_value <= 100),
        cov_between_0_1_n = sum(!is.na(coverage_value) & coverage_value >= 0 & coverage_value <= 1)
      )
    qc_issues[["missing_cov_n"]] <- cov_qc_df$missing_cov_n
    qc_issues[["cov_lt_0_n"]] <- cov_qc_df$cov_lt_0_n
    qc_issues[["cov_gt_100_n"]] <- cov_qc_df$cov_gt_100_n
    qc_issues[["cov_between_1_100_n"]] <- cov_qc_df$cov_between_1_100_n
    qc_issues[["cov_between_0_1_n"]] <- cov_qc_df$cov_between_0_1_n
  }

  if (!is.na(eda_tessy_col)) {
    tessy_vec <- as.character(eda_df[[eda_tessy_col]])
    qc_issues[["missing_tessy_n"]] <- sum(is.na(tessy_vec) | trimws(tessy_vec) == "")
  }

  qc_summary <- data.frame(
    metric = names(qc_issues),
    value = as.numeric(unlist(qc_issues))
  ) %>%
    arrange(desc(value))
  export_graph <- export_to_ppt(export_graph, qc_summary, "SC2 data quality summary")
  export_graph <- export_to_ppt(export_graph, head(column_profile, 30), "SC2 SC2db: highest missingness columns")

  # Numeric outlier scan (IQR-based):
  # include numeric-like columns (comma/dot aware), exclude ID/code fields.
  parse_numeric_locale <- function(v) {
    x <- as.character(v)
    x <- trimws(x)
    x[x == "" | toupper(x) %in% c("NA", "NAN")] <- NA_character_
    x <- gsub(",", ".", x, fixed = TRUE)
    suppressWarnings(as.numeric(x))
  }
  is_id_or_code_col <- function(nm) {
    nm_l <- tolower(nm)
    grepl("(key|lwid|lw_id|sekv.*id|sample.*id|run.*id|(^|_)id($|_))", nm_l) ||
      nm_l %in% c("week", "year", "isoweek", "isoyear", "prove_uke", "prove_sesong", "prove_kategori", "pasient_fylke_nr", "pasient_no")
  }

  all_cols <- names(eda_df)
  candidate_cols <- all_cols[!vapply(all_cols, is_id_or_code_col, logical(1))]
  if (length(candidate_cols) > 0) {
    parsed_map <- lapply(candidate_cols, function(col_name) {
      raw <- eda_df[[col_name]]
      if (is.numeric(raw)) {
        x <- as.numeric(raw)
      } else {
        x <- parse_numeric_locale(raw)
      }
      non_missing_raw <- sum(!is.na(raw) & trimws(as.character(raw)) != "")
      parse_rate <- if (non_missing_raw > 0) sum(!is.na(x)) / non_missing_raw else 0
      list(col = col_name, x = x, parse_rate = parse_rate)
    })
    names(parsed_map) <- candidate_cols
    numeric_cols <- names(parsed_map)[vapply(parsed_map, function(z) z$parse_rate >= 0.8, logical(1))]

    outlier_scan <- lapply(numeric_cols, function(col_name) {
      x <- parsed_map[[col_name]]$x
      x <- x[!is.na(x) & is.finite(x)]
      if (length(x) < 10) {
        return(NULL)
      }
      q1 <- as.numeric(quantile(x, 0.25, na.rm = TRUE))
      q3 <- as.numeric(quantile(x, 0.75, na.rm = TRUE))
      iqr <- q3 - q1
      lower <- q1 - 1.5 * iqr
      upper <- q3 + 1.5 * iqr
      out_n <- sum(x < lower | x > upper, na.rm = TRUE)
      data.frame(
        column_name = col_name,
        n = length(x),
        outlier_n = out_n,
        outlier_pct = round((out_n / length(x)) * 100, 2),
        min = min(x, na.rm = TRUE),
        p25 = q1,
        median = median(x, na.rm = TRUE),
        p75 = q3,
        max = max(x, na.rm = TRUE)
      )
    }) %>%
      bind_rows() %>%
      arrange(desc(outlier_pct), desc(outlier_n))

    if (nrow(outlier_scan) > 0) {
      export_graph <- export_to_ppt(export_graph, head(outlier_scan, 30), "SC2 numeric outlier scan (IQR)")
    } else {
      outlier_msg <- data.frame(
        note = "No numeric columns with >=10 non-missing values after exclusions (week/year removed).",
        numeric_columns_found = length(numeric_cols),
        stringsAsFactors = FALSE
      )
      export_graph <- export_to_ppt(export_graph, outlier_msg, "SC2 numeric outlier scan (IQR)")
    }
  } else {
    outlier_msg <- data.frame(
      note = "No numeric columns found for outlier scan after exclusions (week/year removed).",
      numeric_columns_found = 0L,
      stringsAsFactors = FALSE
    )
    export_graph <- export_to_ppt(export_graph, outlier_msg, "SC2 numeric outlier scan (IQR)")
  }
}

# ============================================================================
# CT DISTRIBUTION BY TESSY (LAST 6 MONTHS)
# Output in this section:
#   - PCR SC2 EXT CT per Tessy-kategori (siste 6 måneder)
# ============================================================================

ct_source_df <- if (exists("SC2db")) SC2db else SC2db
ct_col <- intersect(c("pcr_sc2_ext_ct", "PCR_SC2_EXT_CT"), names(ct_source_df))[1]
tessy_col <- intersect(c("Tessy", "tessy"), names(ct_source_df))[1]
pango_col <- intersect(c("nc_pangolin_short", "NC_Pangolin Short"), names(ct_source_df))[1]
date_col <- intersect(c("prove_tatt", "PROVE_TATT", "sample_date", "Sampledate"), names(ct_source_df))[1]

if (!is.na(ct_col) && !is.na(tessy_col) && !is.na(pango_col) && !is.na(date_col)) {
  ct_plot_df <- ct_source_df %>%
    mutate(
      plot_date = as.Date(.data[[date_col]]),
      ct_raw = as.character(.data[[ct_col]]),
      ct_value = suppressWarnings(as.numeric(.data[[ct_col]])),
      tessy_group = as.character(.data[[tessy_col]]),
      pangolin_short = as.character(.data[[pango_col]])
    ) %>%
    filter(
      !is.na(plot_date),
      plot_date >= (Sys.Date() %m-% months(6)),
      !is.na(ct_raw),
      trimws(ct_raw) != "",
      is.finite(ct_value),
      !is.na(tessy_group),
      trimws(tessy_group) != "",
      !is.na(pangolin_short),
      trimws(pangolin_short) != ""
    )

  if (nrow(ct_plot_df) > 0) {
    ct_x_min <- floor(max(0, quantile(ct_plot_df$ct_value, 0.02, na.rm = TRUE) - 1))
    ct_x_max <- ceiling(min(40, quantile(ct_plot_df$ct_value, 0.98, na.rm = TRUE) + 1))
    if (!is.finite(ct_x_min) || !is.finite(ct_x_max) || ct_x_min >= ct_x_max) {
      ct_x_min <- 0
      ct_x_max <- 40
    }
    ct_x_breaks <- pretty(c(ct_x_min, ct_x_max), n = 8)

    ct_summary_df <- ct_plot_df %>%
      group_by(tessy_group) %>%
      summarise(
        n = n(),
        mean_ct = mean(ct_value, na.rm = TRUE),
        median_ct = median(ct_value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(is.finite(mean_ct), is.finite(median_ct)) %>%
      mutate(tessy_group = forcats::fct_reorder(tessy_group, median_ct))

    top_pangolin_short <- ct_plot_df %>%
      count(pangolin_short, sort = TRUE) %>%
      slice_head(n = 12) %>%
      pull(pangolin_short)

    ct_plot_df <- ct_plot_df %>%
      mutate(
        tessy_group = factor(tessy_group, levels = levels(ct_summary_df$tessy_group)),
        pangolin_short_plot = ifelse(pangolin_short %in% top_pangolin_short, pangolin_short, "Andre")
      )

    ct_by_tessy_plot <- ggplot(ct_plot_df, aes(x = ct_value, y = tessy_group)) +
      geom_boxplot(outlier.shape = NA, width = 0.55, fill = kvantitativ_b2[3], alpha = 0.45, color = kvantitativ_b2[8]) +
      geom_jitter(aes(color = pangolin_short_plot), height = 0.15, width = 0, alpha = 0.45, size = 1.6) +
      geom_point(data = ct_summary_df, aes(x = mean_ct, y = tessy_group), color = kvalitativ_b[2], size = 2.8) +
      geom_point(data = ct_summary_df, aes(x = median_ct, y = tessy_group), color = kvalitativ_b[3], size = 3.0, shape = 18) +
      geom_text(
        data = ct_summary_df,
        aes(x = ct_x_max, y = tessy_group, label = paste0("n=", n)),
        hjust = 1.05,
        size = 2.9,
        color = "grey25"
      ) +
      scale_x_continuous(breaks = ct_x_breaks) +
      scale_color_manual(
        values = fhi_discrete_palette(n_distinct(ct_plot_df$pangolin_short_plot), sc2_palette),
        guide = guide_legend(ncol = 2, override.aes = list(alpha = 1, size = 3))
      ) +
      coord_cartesian(xlim = c(ct_x_min, ct_x_max)) +
      labs(
        title = "PCR SC2 EXT CT per Tessy-kategori (siste 6 måneder)",
        subtitle = "Prikker: prøve | Svart prikk: gjennomsnitt | Gul diamant: median",
        x = "PCR_SC2_EXT_CT",
        y = "Tessy-kategori",
        color = "Pangolin-kortnavn"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")

    export_graph <- export_to_ppt(export_graph, ct_by_tessy_plot, "PCR SC2 EXT CT per Tessy-kategori (siste 6 måneder)")
  } else {
    message("No rows available for Tessy CT plots in the last 6 months.")
  }
} else {
  message("Missing one or more required columns for Tessy CT plots: ct, tessy, pangolin, or date.")
}

# ============================================================================
# NGS COVERAGE PERFORMANCE BY RUN SETUP / RUN ID
# ============================================================================
export_graph <- add_section_slide(export_graph, "Runkvalitet", "Runoppsett, inkludering og operasjonelle avvik per run", c("Inkludering per run", "Dekning per run"))

ngs_qc_source <- if (exists("SC2db")) SC2db else if (exists("SC2db_prefilter")) SC2db_prefilter else if (exists("SC2db")) SC2db else NULL

if (!is.null(ngs_qc_source)) {
  current_week_ngs <- week(Sys.Date())
  current_year_ngs <- year(Sys.Date())
  ngs_date_col <- intersect(c("prove_tatt", "PROVE_TATT", "sample_date", "Sampledate"), names(ngs_qc_source))[1]
  ngs_cov_col <- intersect(c("nc_coverage", "coverage_breadth_artic", "coverage_breadth_swift", "coverage_breadth_eksterne", "coverage_breadth_nano"), names(ngs_qc_source))[1]
  ngs_run_col <- intersect(c("ngs_run_id", "sekv_oppsett_run_artic", "sekv_oppsett_nano", "sekv_oppsett_sanger", "sekv_oppsett_swift"), names(ngs_qc_source))[1]
  ngs_tessy_col <- intersect(c("Tessy", "tessy"), names(ngs_qc_source))[1]
  ngs_subclade_col <- intersect(c("nc_pangolin_short", "NC_Pangolin Short", "nc_clade", "nc_nextclade", "Tessy", "tessy"), names(ngs_qc_source))[1]

  if (!is.na(ngs_cov_col) && !is.na(ngs_run_col)) {
    ngs_qc_df <- ngs_qc_source %>%
      mutate(
        run_setup = as.character(.data[[ngs_run_col]]),
        plot_date = if (!is.na(ngs_date_col)) as.Date(.data[[ngs_date_col]]) else as.Date(NA),
        Tessy_plot = if (!is.na(ngs_tessy_col)) as.character(.data[[ngs_tessy_col]]) else NA_character_,
        cov_raw = as.character(.data[[ngs_cov_col]]),
        cov_num = suppressWarnings(as.numeric(cov_raw)),
        cov_norm = ifelse(!is.na(cov_num) & cov_num > 1.5, cov_num / 100, cov_num),
        spike_ok = if ("spike_mut" %in% names(.)) (!is.na(spike_mut) & trimws(as.character(spike_mut)) != "") else TRUE,
        # Match SC2db inclusion rule from SC2_SQLquery_25-26.
        include_by_coverage = (!is.na(cov_raw) & cov_raw == "NA") | (!is.na(cov_norm) & cov_norm >= 0.7),
        include_flag = include_by_coverage & spike_ok,
        qc_status = ifelse(include_flag, "Included_in_SC2db", "Failed_threshold_or_missing")
      ) %>%
      filter(!is.na(run_setup), trimws(run_setup) != "", run_setup != "Ukjent")

    ngs_perf <- ngs_qc_df %>%
      count(run_setup, qc_status, name = "n") %>%
      group_by(run_setup) %>%
      mutate(percent = (n / sum(n)) * 100) %>%
      ungroup()

    if (nrow(ngs_perf) > 0) {
      run_levels <- ngs_perf %>%
        distinct(run_setup) %>%
        mutate(
          run_num = suppressWarnings(as.numeric(stringr::str_extract(run_setup, "\\d+"))),
          run_suffix = toupper(stringr::str_extract(run_setup, "[A-Za-z]+$")),
          suffix_rank = dplyr::case_when(
            is.na(run_suffix) | run_suffix == "" ~ 0,
            run_suffix == "A" ~ 1,
            run_suffix == "B" ~ 2,
            TRUE ~ 3
          )
        ) %>%
        arrange(is.na(run_num), run_num, suffix_rank, run_setup) %>%
        pull(run_setup)

      # Keep only run setups represented by samples from the last 12 months.
      run_levels_last12m <- ngs_qc_df %>%
        filter(!is.na(plot_date), plot_date >= (Sys.Date() %m-% months(12))) %>%
        distinct(run_setup) %>%
        pull(run_setup) %>%
        as.character()
      run_levels_last12m <- run_levels[run_levels %in% run_levels_last12m]

      if (length(run_levels_last12m) == 0) {
        run_levels_last12m <- run_levels
      }
      ngs_qc_df <- ngs_qc_df %>%
        filter(run_setup %in% run_levels_last12m) %>%
        mutate(run_setup = factor(run_setup, levels = run_levels_last12m))
      ngs_perf <- ngs_perf %>%
        filter(run_setup %in% run_levels_last12m) %>%
        mutate(run_setup = factor(run_setup, levels = run_levels_last12m))

      p_ngs_perf_pct <- ggplot(ngs_perf, aes(x = run_setup, y = percent, fill = qc_status)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = fhi_discrete_palette(n_distinct(ngs_perf$qc_status), sc2_palette)) +
        scale_y_continuous(labels = scales::percent_format(scale = 1)) +
        coord_cartesian(ylim = c(0, 100)) +
        labs(
          title = "NGS-runoppsett: andel inkludert vs. ikke inkludert",
          x = "NGS-runoppsett / run-id",
          y = "Andel (%)",
          fill = "Status"
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      p_ngs_perf_n <- ggplot(ngs_perf, aes(x = run_setup, y = n, fill = qc_status)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = fhi_discrete_palette(n_distinct(ngs_perf$qc_status), sc2_palette)) +
        labs(
          title = "NGS-runoppsett: antall inkludert vs. ikke inkludert",
          x = "NGS-runoppsett / run-id",
          y = "Antall (n)",
          fill = "Status"
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p_ngs_perf_combined <- p_ngs_perf_pct / p_ngs_perf_n +
        patchwork::plot_layout(heights = c(1, 1), guides = "collect") &
        theme(legend.position = "right")
      export_graph <- export_to_ppt(export_graph, p_ngs_perf_combined, "NGS-runoppsett: inkludert vs. ikke inkludert")
    }

    # Coverage by run setup: box-and-whisker per run with per-sample points colored by Tessy.
    ngs_cov_run_df <- ngs_qc_df %>%
      filter(!is.na(cov_norm), !is.na(run_setup), trimws(run_setup) != "", run_setup != "Ukjent")

    if (nrow(ngs_cov_run_df) > 0) {
      cov_run_summary <- ngs_qc_df %>%
        group_by(run_setup) %>%
        summarise(
          n_samples = n(),
          n_cov_samples = sum(!is.na(cov_norm)),
          mean_cov = round(mean(cov_norm, na.rm = TRUE), 3),
          median_cov = round(median(cov_norm, na.rm = TRUE), 3),
          p10_cov = round(as.numeric(stats::quantile(cov_norm, probs = 0.10, na.rm = TRUE)), 3),
          p90_cov = round(as.numeric(stats::quantile(cov_norm, probs = 0.90, na.rm = TRUE)), 3),
          included_n = sum(qc_status == "Included_in_SC2db", na.rm = TRUE),
          failed_n = sum(qc_status == "Failed_threshold_or_missing", na.rm = TRUE),
          included_pct = round(100 * included_n / n_samples, 1),
          .groups = "drop"
        ) %>%
        arrange(desc(included_pct), desc(median_cov), desc(n_samples), run_setup)

      export_graph <- export_to_ppt(
        export_graph,
        cov_run_summary,
        "Deknings-QC per NGS-run (oppsummeringstabell)"
      )

      if ("Tessy_plot" %in% names(ngs_cov_run_df) && any(!is.na(ngs_cov_run_df$Tessy_plot) & trimws(ngs_cov_run_df$Tessy_plot) != "")) {
        top_tessy_run <- ngs_cov_run_df %>%
          filter(!is.na(Tessy_plot), trimws(Tessy_plot) != "") %>%
          count(Tessy_plot, sort = TRUE) %>%
          slice_head(n = 8) %>%
          pull(Tessy_plot)

        ngs_cov_run_df <- ngs_cov_run_df %>%
          mutate(
            Tessy_plot = ifelse(is.na(Tessy_plot) | trimws(Tessy_plot) == "", "Ukjent", Tessy_plot),
            Tessy_plot = ifelse(Tessy_plot %in% top_tessy_run, Tessy_plot, "Andre")
          )

        p_cov_run_box <- ggplot(ngs_cov_run_df, aes(x = run_setup, y = cov_norm)) +
          geom_boxplot(fill = kvantitativ_b2[2], color = kvantitativ_b2[8], outlier.shape = NA) +
          geom_jitter(aes(color = Tessy_plot), width = 0.2, height = 0, alpha = 0.7, size = 1.6) +
          scale_color_manual(values = fhi_discrete_palette(n_distinct(ngs_cov_run_df$Tessy_plot), sc2_palette)) +
          labs(
            title = "Dekning per NGS-runoppsett",
            subtitle = "Boksplott per run; hver pr\u00f8ve farget etter Tessy (topp 8 + Andre)",
            x = "NGS-runoppsett / run-id",
            y = "Normalisert dekningsgrad (0-1)",
            color = "Tessy"
          ) +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        p_cov_run_box <- ggplot(ngs_cov_run_df, aes(x = run_setup, y = cov_norm)) +
          geom_boxplot(fill = kvantitativ_b2[2], color = kvantitativ_b2[8], outlier.shape = NA) +
          geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 1.5, color = kvalitativ_b[2]) +
          labs(
            title = "Dekning per NGS-runoppsett",
            subtitle = "Boksplott per run; Tessy ikke tilgjengelig for punktfarging",
            x = "NGS-runoppsett / run-id",
            y = "Normalisert dekningsgrad (0-1)"
          ) +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }

      export_graph <- export_to_ppt(export_graph, p_cov_run_box, "Dekning per runoppsett (boksplott + pr\u00f8vepunkter)")
    }

    # Coverage by month colored by subclade/lineage.
    if (!is.na(ngs_subclade_col)) {
      ngs_cov_month_df <- ngs_qc_df %>%
        mutate(
          month_date = floor_date(plot_date, unit = "month"),
          subclade_plot = as.character(.data[[ngs_subclade_col]]),
          subclade_plot = ifelse(is.na(subclade_plot) | trimws(subclade_plot) == "", "Ukjent", subclade_plot)
        ) %>%
        filter(!is.na(cov_norm), !is.na(month_date))

      if (nrow(ngs_cov_month_df) > 0) {
        top_subclades <- ngs_cov_month_df %>%
          count(subclade_plot, sort = TRUE) %>%
          slice_head(n = 10) %>%
          pull(subclade_plot)

        ngs_cov_month_df <- ngs_cov_month_df %>%
          mutate(subclade_plot = ifelse(subclade_plot %in% top_subclades, subclade_plot, "Andre")) %>%
          group_by(month_date, subclade_plot) %>%
          summarise(
            mean_cov = mean(cov_norm, na.rm = TRUE),
            median_cov = median(cov_norm, na.rm = TRUE),
            n = n(),
            .groups = "drop"
          )

        p_cov_month_subclade <- ggplot(
          ngs_cov_month_df,
          aes(x = month_date, y = median_cov, color = subclade_plot, group = subclade_plot)
        ) +
          geom_line(linewidth = 1) +
          geom_point(aes(size = n), alpha = 0.85) +
          scale_color_manual(values = fhi_discrete_palette(n_distinct(ngs_cov_month_df$subclade_plot), sc2_palette)) +
          scale_size_continuous(range = c(1.5, 5)) +
          scale_x_date(date_breaks = "1 month", labels = format_month_label) +
          labs(
            title = "Dekning per m\u00e5ned farget etter subklade",
            subtitle = "Median dekningsgrad per m\u00e5ned; topp 10 subklader vist separat",
            x = "M\u00e5ned",
            y = "Normalisert dekningsgrad (0-1)",
            color = "Subklade",
            size = "n"
          ) +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        export_graph <- export_to_ppt(export_graph, p_cov_month_subclade, "Dekning per m\u00e5ned og subklade")
      }
    }
  }
}


# ============================================================================
# SEQUENCE SUMMARY - WEEK/MONTH
# ============================================================================

export_graph <- add_section_slide(
  export_graph,
  "Sekvenseringssted",
  "Hvor SC2-prøver ble sekvensert"
)
# ---- BEGIN INLINED: SC2/SC2_Seqs_per_month.R ----
###### Sekvenser per uke for prosentberegning: ######

if (!exists("SC2db")) {
  stop("Object 'SC2db' is missing. Run the SQL/classification scripts before sourcing SC2_Seqs_per_month.R.")
}

###### Sekvenser per m\u00e5ned for prosentberegning: ######

# Calculate Sequences per month for Spike protein sequence results
spm_spike <- SC2db %>%
  filter(nc_pangolin_short != "") %>%
  filter(spike_mut != "") %>%
  count(my, name = "TotalSeq") %>%
  ungroup() %>%
  mutate(
    Date = as.Date(my),
    YearMonth = format_month_key_nb(Date)
  ) %>%
  select(-Date) # Drop the temporary Date column if not needed

# Calculate Total Valid Sequences per month
v_seqs_per_month <- SC2db %>%
  count(my, name = "TotalSeq")

###### Sekvenser per opprinnelseslaboratorium #######

v_seqs_per_month_origin <- SC2db %>%
  group_by(my, Origin) %>%
  count(name = "TotalSeq") %>%
  ungroup() %>%
  mutate(my = parse_month_key_nb(my))


# Create a bar chart per month based on Origin
spmlabto <- ggplot(v_seqs_per_month_origin, aes(x = my, y = TotalSeq, fill = Origin)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Sekvenseringssted for SC2-prøver",
    x = "M\u00e5ned",
    y = "Antall (n)"
  ) +
  scale_x_date(
    breaks = "2 months", # Reduce crowding on long time series slides
    labels = format_month_year_nb # Norsk månedsnavn og år
  ) +
  scale_fill_manual(values = kvalitativ_a) + # Set custom colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) + # Center the title and increase its size
  theme(legend.position = "right") # Optional: Remove legend if not needed


export_graph <- export_to_ppt(export_graph, spmlabto, "Sekvenseringssted for SC2-prøver")


# Assuming 'my' is the column in the dataframe with date format "YYYY-MM-DD"
# Convert 'my' to Date type if not already
v_seqs_per_month_origin <- v_seqs_per_month_origin %>%
  mutate(my = as.Date(my))

# Get the system date
current_date <- Sys.Date()

# Calculate the date threshold for filtering (12 months prior to today)
threshold_date <- current_date %m-% months(12)

# Filter the dataframe for rows with 'my' from the last 12 months
v_seqs_per_month_origin12m <- v_seqs_per_month_origin %>%
  filter(my >= threshold_date)

# Create a bar chart per month based on Origin
spmlabto12m <- ggplot(v_seqs_per_month_origin12m, aes(x = my, y = TotalSeq, fill = Origin)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Sekvenseringssted for SC2-prøver",
    x = "M\u00e5ned",
    y = "Antall (n)"
  ) +
  scale_x_date(
    breaks = "1 month", # Show breaks every month
    labels = format_month_year_nb # Norsk månedsnavn og år
  ) +
  scale_fill_manual(values = kvalitativ_a) + # Set custom colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) + # Center the title and increase its size
  theme(legend.position = "right") # Optional: Remove legend if not needed



export_graph <- export_to_ppt(
  export_graph,
  spmlabto12m,
  "Sekvenseringssted for SC2-prøver",
  build_slide_subtitle(v_seqs_per_month_origin12m, "my", sum(v_seqs_per_month_origin12m$TotalSeq, na.rm = TRUE), "Siste 12 måneder")
)

# ---- END INLINED: SC2/SC2_Seqs_per_month.R ----


# ============================================================================
# PANGOLIN AND TESSY OVERVIEW
# Output in this section:
#   - Pangolin monthly overview
#   - Tessy monthly overview
#   - Top-line lineage tables for recent data windows
# ============================================================================

export_graph <- add_section_slide(
  export_graph,
  "Pangolin og Tessy",
  "Klassifiseringer per m\u00e5ned",
  c("Pangolin per m\u00e5ned", "Tessy siste 12 mnd", "Tessy siste 6 mnd")
)
# ---- BEGIN INLINED: SC2/SC2_Pangolin_p_m.R ----
sequencing_window_start <- if (exists("data_window_start")) as.Date(data_window_start) else (Sys.Date() %m-% months(6))

# Prepare data for the weekly sequence count
monthcount <- SC2db %>%
  dplyr::count(my, name = "TotalSeq")

# Calculate weekly counts for each Pangolin lineage
pangomtcount <- SC2db %>%
  group_by(my) %>%
  dplyr::count(Collapsed_pango, my, name = "count") %>%
  ungroup() %>%
  mutate(Sampledate = parse_month_key_nb(my)) %>%
  left_join(monthcount, by = "my") %>%
  mutate(Percent = count / TotalSeq)

# Count sequences per Pangolin lineage and Collapsed_pango variant for each month
fpangomtcount <- SC2db %>%
  group_by(my) %>%
  dplyr::count(nc_pangolin_short, Collapsed_pango, my, name = "count") %>%
  ungroup() %>%
  mutate(Sampledate = parse_month_key_nb(my)) %>%
  left_join(monthcount, by = "my") %>%
  mutate(Percent = count / TotalSeq)
# Generate and save whole-period Pangolin overview
collapsed_levels_full <- unique(as.character(pangomtcount$Collapsed_pango))
collapsed_colors_full <- fhi_discrete_palette(length(collapsed_levels_full), sc2_palette)
names(collapsed_colors_full) <- collapsed_levels_full

combined_plotwp <- make_stacked_percent_count_plot(
  data = pangomtcount,
  date_col = "Sampledate",
  percent_col = "Percent",
  count_col = "TotalSeq",
  category_col = "Collapsed_pango",
  fill_values = collapsed_colors_full,
  fill_label = "Pangolin",
  x_breaks = "2 months",
  legend_ncol = 2,
  legend_position = "right"
)

export_graph <- export_to_ppt(
  export_graph,
  combined_plotwp,
  "Pangolin per m\u00e5ned",
  build_slide_subtitle(pangomtcount, "Sampledate", nrow(SC2db), "Hele perioden | Stablede andeler øverst og volum nederst")
)


# Convert `my` to the desired format and correctly order the data
pangoxls <- pangomtcount %>%
  select(my, Collapsed_pango, count) %>%
  mutate(my = format_month_key_nb(as.Date(my))) %>%
  pivot_wider(
    names_from = my, # Column that defines the new columns
    values_from = count # Values for the new columns
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

# Convert back to the long format with correct ordering
pangostatistikk <- pangoxls %>%
  pivot_longer(
    cols = -Collapsed_pango, # Pivot all columns except 'Collapsed_pango'
    names_to = "my", # Name for the new key column
    values_to = "count" # Name for the new value column
  ) %>%
  mutate(flagg = 0) %>%
  mutate(my_ord = parse_month_key_nb(my)) %>%
  arrange(my_ord) %>%
  select(-my_ord) # Remove the ordering column

# Calculate total counts for each 'my'
total_counts <- pangostatistikk %>%
  group_by(my) %>%
  summarize(total_count = sum(count))

# Join total counts back to the data
pangostatistikk_with_totals <- pangostatistikk %>%
  left_join(total_counts, by = "my") %>%
  mutate(percent = round((count / total_count) * 100)) # Calculate and round percentage to whole number


# Select and arrange the final data
final_pangostatistikk <- pangostatistikk_with_totals %>%
  select(Collapsed_pango, my, count, percent, flagg)

# Ensure numeric values are formatted correctly for CSV export
final_pangostatistikk <- final_pangostatistikk %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ format(., scientific = FALSE) # Ensure numeric values are correctly formatted
  ))


recent_12m_start <- Sys.Date() %m-% months(12)
recent_6m_start <- Sys.Date() %m-% months(6)

pango_recent12 <- pangomtcount %>%
  filter(Sampledate >= recent_12m_start) %>%
  collapse_minor_categories(
    category_col = "Collapsed_pango",
    value_col = "Percent",
    other_label = "Andre SARS CoV 2",
    min_peak_percent = 5,
    top_n = 8,
    preserve_categories = c("Andre SARS CoV 2")
  ) %>%
  group_by(Sampledate, Collapsed_pango_plot) %>%
  summarise(count = sum(count), TotalSeq = first(TotalSeq), Percent = sum(Percent), .groups = "drop")

pango_recent6 <- pangomtcount %>%
  filter(Sampledate >= recent_6m_start) %>%
  collapse_minor_categories(
    category_col = "Collapsed_pango",
    value_col = "Percent",
    other_label = "Andre SARS CoV 2",
    min_peak_percent = 5,
    top_n = 8,
    preserve_categories = c("Andre SARS CoV 2")
  ) %>%
  group_by(Sampledate, Collapsed_pango_plot) %>%
  summarise(count = sum(count), TotalSeq = first(TotalSeq), Percent = sum(Percent), .groups = "drop")

pango_recent12_levels <- unique(as.character(pango_recent12$Collapsed_pango_plot))
pango_recent12_colors <- fhi_discrete_palette(length(pango_recent12_levels), sc2_palette)
names(pango_recent12_colors) <- pango_recent12_levels

pango_recent6_levels <- unique(as.character(pango_recent6$Collapsed_pango_plot))
pango_recent6_colors <- fhi_discrete_palette(length(pango_recent6_levels), sc2_palette)
names(pango_recent6_colors) <- pango_recent6_levels

combined_plot_12mo <- make_percent_trend_plot(
  data = pango_recent12,
  date_col = "Sampledate",
  percent_col = "Percent",
  count_col = "TotalSeq",
  category_col = "Collapsed_pango_plot",
  color_values = pango_recent12_colors,
  color_label = "Pangolin",
  x_breaks = "1 month",
  legend_nrow = 2
)

export_graph <- export_to_ppt(
  export_graph,
  combined_plot_12mo,
  "Pangolin siste 12 mnd",
  build_slide_subtitle(pango_recent12, "Sampledate", nrow(SC2db %>% filter(prove_tatt >= recent_12m_start)), "Andeler måned for måned | Stablede andeler øverst | Volum og trender nederst")
)

combined_plot_6mo <- make_percent_trend_plot(
  data = pango_recent6,
  date_col = "Sampledate",
  percent_col = "Percent",
  count_col = "TotalSeq",
  category_col = "Collapsed_pango_plot",
  color_values = pango_recent6_colors,
  color_label = "Pangolin",
  x_breaks = "1 month",
  legend_nrow = 2
)

export_graph <- export_to_ppt(
  export_graph,
  combined_plot_6mo,
  "Pangolin siste 6 mnd",
  build_slide_subtitle(pango_recent6, "Sampledate", nrow(SC2db %>% filter(prove_tatt >= recent_6m_start)), "Andeler måned for måned | Stablede andeler øverst | Volum og trender nederst")
)

# Create and save individual Pangolin variant plots for the last 12 months
subset_data4mofpango <- subset(fpangomtcount, Sampledate >= Sys.Date() %m-% months(4))
unique_collapsed_pangosrec <- subset_data4mofpango %>%
  pull(Collapsed_pango) %>%
  unique()


# Create a summary table of Pangolin variants for the shared sequencing window
subset_data6mopango <- subset(fpangomtcount, Sampledate >= sequencing_window_start)

# Loop through unique 'Collapsed_pango' values
for (collapsed_pango in unique_collapsed_pangosrec) {
  loop_started_at <- Sys.time()
  if (exists("log_timed_message", mode = "function")) {
    log_timed_message("Loop Pangolin START: ", collapsed_pango)
  }
  subset_data <- subset_data6mopango %>%
    filter(Collapsed_pango == collapsed_pango) %>%
    mutate(Sampledate = as.Date(Sampledate))

  subset_data_plot <- subset_data %>%
    mutate(Sampledate = as.Date(Sampledate)) %>%
    collapse_minor_categories(
      category_col = "nc_pangolin_short",
      value_col = "count",
      other_label = "Andre undervarianter",
      min_peak_percent = 3,
      top_n = 23,
      preserve_categories = c(
        "Andre undervarianter",
        subset_data %>%
          filter(Sampledate >= (Sys.Date() %m-% months(3))) %>%
          pull(nc_pangolin_short) %>%
          as.character() %>%
          unique()
      )
    ) %>%
    group_by(Sampledate, nc_pangolin_short_plot) %>%
    summarise(count = sum(count), Percent = sum(Percent), .groups = "drop")

  collapsed_pangosrecgr_percent <- make_monthly_heatmap(
    data = subset_data_plot,
    date_col = "Sampledate",
    category_col = "nc_pangolin_short_plot",
    value_col = "Percent",
    fill_label = "Andel (%)",
    x_breaks = "1 month",
    bottom_categories = c("Andre undervarianter"),
    base_y_text_size = 8
  ) +
    labs(title = paste("Andel av", collapsed_pango, "undervarianter de siste seks månedene"))

  export_graph <- export_to_ppt(
    export_graph,
    collapsed_pangosrecgr_percent,
    paste(collapsed_pango, "andel"),
    build_slide_subtitle(subset_data, "Sampledate", sum(subset_data$count, na.rm = TRUE), "Siste seks måneder | Heatmap av månedlige andeler per undervariant")
  )
  if (exists("log_timed_message", mode = "function")) {
    loop_elapsed <- as.numeric(difftime(Sys.time(), loop_started_at, units = "secs"))
    log_timed_message("Loop Pangolin DONE: ", collapsed_pango, " (", sprintf("%.2f", loop_elapsed), "s)")
  }
}

p6modata <- subset_data6mopango %>%
  mutate(n = paste0(count, " (", round(Percent, 2), "%", ")"))

pangolin_table_full <- p6modata %>%
  pivot_wider(
    id_cols = nc_pangolin_short,
    names_from = my,
    values_from = n,
  ) %>%
  rename("SARS-CoV2 Variants" = nc_pangolin_short)

pangolin_table_full <- pangolin_table_full %>%
  mutate(last_col_numeric = as.numeric(str_extract(.[[ncol(.)]], "^[0-9]+"))) %>%
  arrange(desc(last_col_numeric)) %>%
  select(-last_col_numeric)

ft <- flextable(pangolin_table_full %>% slice_head(n = 15))
ft <- autofit(ft)
export_graph <- export_to_ppt(
  export_graph,
  ft,
  "Pangolin topp 15",
  build_slide_subtitle(subset_data6mopango, "Sampledate", nrow(subset_data6mopango), "Siste seks måneder | Full tabell eksporteres til Excel")
)
#####

subset_data6mopango$count <- as.numeric(subset_data6mopango$count)

pangolin_count_table_full <- subset_data6mopango %>%
  pivot_wider(
    id_cols = nc_pangolin_short,
    names_from = my,
    values_from = count,
    values_fill = 0
  ) %>%
  rename("SARS-CoV2 Variants" = nc_pangolin_short) %>%
  arrange(desc(.[[ncol(.)]]))

# ---- END INLINED: SC2/SC2_Pangolin_p_m.R ----

# ---- BEGIN INLINED: SC2/SC2_Tessy_p_m.R ----
# Output in this section:
#   - Whole-period Tessy overview
#   - Per-Tessy collapsed-pango stacked chart + andel heatmap
#   - Tessy tables and recency trend slides

sequencing_window_start <- if (exists("data_window_start")) as.Date(data_window_start) else (Sys.Date() %m-% months(6))

# Count the number of sequences per month
monthcount <- SC2db %>%
  dplyr::count(my, name = "TotalSeq")

# Count the number of sequences for each tessy lineage per month
tessymtcount <- SC2db %>%
  group_by(my) %>%
  dplyr::count(Tessy, my, name = "count") %>%
  ungroup() %>%
  mutate(Sampledate = parse_month_key_nb(my))

# Merge the lineage counts with the total sequence counts
tessymtcount <- tessymtcount %>%
  left_join(monthcount, by = "my") %>%
  mutate(Percent = (count / TotalSeq) * 100)

# Add Sampledate to monthcount
monthcount <- monthcount %>%
  mutate(Sampledate = parse_month_key_nb(my))

# Subset data for shared sequencing window + short recency views
tessy12mo <- subset(SC2db, prove_tatt >= recent_12m_start & prove_tatt <= Sys.Date())
tessy6mo <- subset(SC2db, prove_tatt >= recent_6m_start & prove_tatt <= Sys.Date())
subset_data12mo <- subset(tessymtcount, Sampledate >= recent_12m_start)
subset_data2mo <- subset(tessymtcount, Sampledate >= Sys.Date() %m-% months(2))
subset_data4mo <- subset(tessymtcount, Sampledate >= Sys.Date() %m-% months(4))
subset_data6mo <- subset(tessymtcount, Sampledate >= recent_6m_start)
monthcount12mo <- subset(monthcount, Sampledate >= recent_12m_start)

subset_data_season <- subset(tessymtcount, Sampledate >= recent_6m_start)
subset_data_year <- subset(tessymtcount, Sampledate >= recent_12m_start)
subset_data_season_p <- subset_data_season %>%
  mutate(Percent = (count / TotalSeq) * 100)
subset_data_year_p <- subset_data_year %>%
  mutate(Percent = (count / TotalSeq) * 100)

tessy_levels_full <- unique(as.character(tessymtcount$Tessy))
tessy_colors_full <- fhi_discrete_palette(length(tessy_levels_full), sc2_palette)
names(tessy_colors_full) <- tessy_levels_full

tessy_stacked_plot <- make_stacked_percent_count_plot(
  data = tessymtcount,
  date_col = "Sampledate",
  percent_col = "Percent",
  count_col = "TotalSeq",
  category_col = "Tessy",
  fill_values = tessy_colors_full,
  fill_label = "Tessy-klassifisering",
  x_breaks = "2 months",
  legend_position = "right",
  legend_ncol = 2
)

export_graph <- export_to_ppt(
  export_graph,
  tessy_stacked_plot,
  "Tessy per måned - stablet",
  build_slide_subtitle(tessymtcount, "Sampledate", nrow(SC2db), "Hele perioden | Stablede andeler øverst og volum nederst")
)

tessy_heatmap_plot <- make_percent_count_plot(
  data = tessymtcount,
  x_col = "Sampledate",
  y_percent_col = "Percent",
  fill_col = "Tessy",
  count_data = monthcount,
  count_y_col = "TotalSeq",
  fill_values = tessy_colors_full,
  fill_label = "Tessy-klassifisering",
  x_breaks = "2 months",
  legend_nrow = 3,
  heatmap_bottom_categories = c("Andre SARS CoV 2", "Ukjent"),
  heatmap_y_text_size = 8
)

export_graph <- export_to_ppt(
  export_graph,
  tessy_heatmap_plot,
  "Tessy per måned - varmekart",
  build_slide_subtitle(tessymtcount, "Sampledate", nrow(SC2db), "Hele perioden | Heatmap øverst og volum nederst")
)

############################################### Tessy-kategori single plots

unique_tessyrec <- subset_data4mo %>%
  pull(Tessy) %>%
  unique() %>%
  as.character()

for (i in seq_along(unique_tessyrec)) {
  loop_started_at <- Sys.time()
  tessy <- unique_tessyrec[i]
  if (exists("log_timed_message", mode = "function")) {
    log_timed_message("Loop Tessy ", i, "/", length(unique_tessyrec), " START: ", tessy)
  }

  subset_data <- tessy6mo %>%
    filter(Tessy == tessy) %>%
    group_by(my) %>%
    dplyr::count(nc_pangolin_short, my, name = "Count") %>%
    ungroup() %>%
    left_join(SC2db %>% dplyr::count(my, name = "TotalSeq"), by = "my") %>%
    arrange(my) %>%
    mutate(Sampledate = parse_month_key_nb(my))

  if (nrow(subset_data) == 0 || sum(subset_data$Count, na.rm = TRUE) <= 0) {
    next
  }

  recent_sublineages <- subset_data %>%
    filter(Sampledate >= (Sys.Date() %m-% months(3))) %>%
    pull(nc_pangolin_short) %>%
    as.character() %>%
    unique()

  monthly_totals <- subset_data %>%
    group_by(Sampledate) %>%
    summarise(TotalCount = sum(Count), .groups = "drop")

  subset_data_percent_plot <- subset_data %>%
    group_by(Sampledate) %>%
    mutate(Percent = Count / sum(Count) * 100) %>%
    ungroup() %>%
    collapse_minor_categories(
      category_col = "nc_pangolin_short",
      value_col = "Count",
      other_label = "Andre undervarianter",
      min_peak_percent = 3,
      top_n = 23,
      preserve_categories = c("Andre undervarianter", recent_sublineages)
    ) %>%
    group_by(Sampledate, nc_pangolin_short_plot) %>%
    summarise(Percent = sum(Percent), Count = sum(Count), .groups = "drop") %>%
    left_join(monthly_totals, by = "Sampledate")

  tessy_levels_plot <- unique(as.character(subset_data_percent_plot$nc_pangolin_short_plot))
  tessy_fill_plot <- fhi_discrete_palette(length(tessy_levels_plot), sc2_palette)
  names(tessy_fill_plot) <- tessy_levels_plot

  collapsed_pangosrecgr <- make_stacked_percent_count_plot(
    data = subset_data_percent_plot,
    date_col = "Sampledate",
    percent_col = "Percent",
    count_col = "TotalCount",
    category_col = "nc_pangolin_short_plot",
    fill_values = tessy_fill_plot,
    fill_label = "Pangolin-nomenklatur",
    x_breaks = "1 month",
    legend_position = "right",
    legend_ncol = 2
  )

  export_graph <- export_to_ppt(
    export_graph,
    collapsed_pangosrecgr,
    paste(tessy, "pangolin per måned"),
    build_slide_subtitle(subset_data, "Sampledate", sum(subset_data$Count, na.rm = TRUE), "Siste seks måneder | Stablede pangolinandeler øverst og pangolinantall nederst")
  )

  collapsed_pangosrecgr_percent <- make_monthly_heatmap(
    data = subset_data_percent_plot,
    date_col = "Sampledate",
    category_col = "nc_pangolin_short_plot",
    value_col = "Percent",
    fill_label = "Andel (%)",
    x_breaks = "1 month",
    bottom_categories = c("Andre undervarianter"),
    base_y_text_size = 8
  ) +
    labs(title = paste("Pangolin-nomenklatur per Tessy-klassifisering", tessy))

  export_graph <- export_to_ppt(
    export_graph,
    collapsed_pangosrecgr_percent,
    paste(tessy, "andel"),
    build_slide_subtitle(subset_data, "Sampledate", sum(subset_data$Count, na.rm = TRUE), "Siste seks måneder | Heatmap av pangolinandeler per måned")
  )
  if (exists("log_timed_message", mode = "function")) {
    loop_elapsed <- as.numeric(difftime(Sys.time(), loop_started_at, units = "secs"))
    log_timed_message("Loop Tessy ", i, "/", length(unique_tessyrec), " DONE: ", tessy, " (", sprintf("%.2f", loop_elapsed), "s)")
  }
}

# ################Table last 4 months Tessy frequency ############################

p4modata <- subset_data4mo %>% mutate(n = paste0(count, " (", round(Percent, 2), "%", ")"))

tessy_table_full <- p4modata %>% pivot_wider(
  names_from = my, id_cols = Tessy, values_from = n,
  values_fill = list(count = 0, Percent = 0)
)

ft <- flextable(tessy_table_full %>% slice_head(n = 15))
export_graph <- export_to_ppt(
  export_graph,
  ft,
  "Tessy topp 15",
  build_slide_subtitle(subset_data4mo, "Sampledate", nrow(subset_data4mo), "Siste fire måneder | Full tabell eksporteres til Excel")
)

########################## Combined bar and line plot with tessy variables################

subset_data12mo_mod <- subset_data12mo %>%
  collapse_minor_categories(
    category_col = "Tessy",
    value_col = "Percent",
    other_label = "Andre SARS CoV 2",
    min_peak_percent = 5,
    top_n = 8,
    preserve_categories = c("Andre SARS CoV 2")
  ) %>%
  group_by(Sampledate, Tessy_plot) %>%
  summarise(
    count = sum(count),
    TotalSeq = first(TotalSeq),
    .groups = "drop"
  )

subset_data12mo_modp <- subset_data12mo_mod %>%
  mutate(Percent = (count / TotalSeq) * 100)

tessy_recent12_levels <- unique(as.character(subset_data12mo_modp$Tessy_plot))
tessy_recent12_colors <- fhi_discrete_palette(length(tessy_recent12_levels), sc2_palette)
names(tessy_recent12_colors) <- tessy_recent12_levels

combined_plot <- make_percent_trend_plot(
  data = subset_data12mo_modp,
  date_col = "Sampledate",
  percent_col = "Percent",
  count_col = "TotalSeq",
  category_col = "Tessy_plot",
  color_values = tessy_recent12_colors,
  color_label = "Tessy-klassifisering",
  x_breaks = "1 month",
  legend_position = "right",
  legend_ncol = 2
)

export_graph <- export_to_ppt(
  export_graph,
  combined_plot,
  "Tessy siste 12 måneder",
  build_slide_subtitle(subset_data12mo_modp, "Sampledate", nrow(tessy12mo), "Andeler måned for måned | Stablede andeler øverst | Volum og trender nederst")
)

subset_data_season_plot <- subset_data_season_p %>%
  collapse_minor_categories(
    category_col = "Tessy",
    value_col = "Percent",
    other_label = "Andre SARS CoV 2",
    min_peak_percent = 5,
    top_n = 8,
    preserve_categories = c("Andre SARS CoV 2")
  ) %>%
  group_by(Sampledate, Tessy_plot) %>%
  summarise(
    count = sum(count),
    TotalSeq = first(TotalSeq),
    Percent = sum(Percent),
    .groups = "drop"
  )

tessy_recent6_levels <- unique(as.character(subset_data_season_plot$Tessy_plot))
tessy_recent6_colors <- fhi_discrete_palette(length(tessy_recent6_levels), sc2_palette)
names(tessy_recent6_colors) <- tessy_recent6_levels

combined_plot <- make_percent_trend_plot(
  data = subset_data_season_plot,
  date_col = "Sampledate",
  percent_col = "Percent",
  count_col = "TotalSeq",
  category_col = "Tessy_plot",
  color_values = tessy_recent6_colors,
  color_label = "Tessy-klassifisering",
  x_breaks = "1 month",
  legend_position = "right",
  legend_ncol = 2
)

export_graph <- export_to_ppt(
  export_graph,
  combined_plot,
  "Tessy siste 6 måneder",
  build_slide_subtitle(subset_data_season_p, "Sampledate", nrow(tessy6mo), "Andeler måned for måned | Stablede andeler øverst | Volum og trender nederst")
)

# ---- END INLINED: SC2/SC2_Tessy_p_m.R ----


# ============================================================================
# MUTATION ANALYSIS - COMBINATIONS
# ============================================================================

export_graph <- add_section_slide(export_graph, "Mutasjoner", "Mutasjonskombinasjoner, frekvenser og proteindomener", c("Spike-mutasjoner", "Kombinasjoner", "Domenelollipop"))
# ---- BEGIN INLINED: SC2/SC2_spike_mut_of_interest.R ----
# Extract relevant mutation data and perform initial filtering and transformations
mutfr <- SC2db %>%
  select(prove_tatt, spike_mut, nc_pangolin_short, Collapsed_pango, Tessy) %>%
  filter(prove_tatt >= Sys.Date() - 365) %>%
  mutate(
    Sampledate = as.Date(prove_tatt), # Convert prove_tatt to Date format
    Substitution = gsub(";", ",", spike_mut), # Replace semicolons with commas
    month = format(Sampledate, "%Y-%m"), # Extract Year-Month
    Tessy_group = as.character(Tessy),
    Tessy_group = ifelse(is.na(Tessy_group) | trimws(Tessy_group) == "", "Ukjent", Tessy_group)
  ) %>%
  filter(Substitution != "") %>%
  group_by(month, Tessy_group, Substitution) %>%
  summarize(n = n(), .groups = "drop") %>%
  ungroup() %>%
  mutate(
    n = ifelse(is.na(n), 0, n),
    n_mut = str_count(Substitution, ","),
    date = as.Date(paste0(month, "-01")) # Set to first day of the month
  ) %>%
  arrange(date) # Sort by date

mutfr_monthly <- mutfr %>%
  group_by(date, Tessy_group) %>%
  summarise(count = sum(n), .groups = "drop") %>%
  group_by(date) %>%
  mutate(
    TotalDetections = sum(count, na.rm = TRUE),
    Percent = ifelse(TotalDetections > 0, (count / TotalDetections) * 100, 0)
  ) %>%
  ungroup()

mutfr_levels <- unique(as.character(mutfr_monthly$Tessy_group))
mutfr_colors <- fhi_discrete_palette(length(mutfr_levels), sc2_palette)
names(mutfr_colors) <- mutfr_levels

grnmutfr <- make_stacked_percent_count_plot(
  data = mutfr_monthly,
  date_col = "date",
  percent_col = "Percent",
  count_col = "TotalDetections",
  category_col = "Tessy_group",
  fill_values = mutfr_colors,
  fill_label = "Tessy-klassifisering",
  x_breaks = "1 month",
  legend_position = "right",
  legend_ncol = 2
)

export_graph <- export_to_ppt(
  export_graph,
  grnmutfr,
  "Spike-mutasjoner av interesse",
  build_slide_subtitle(mutfr_monthly, "date", sum(mutfr_monthly$count, na.rm = TRUE), "Siste 12 måneder | Stablede andeler øverst | Totale deteksjoner nederst")
)

# Create a heatmap for mutation combinations over the last 12 months
S_mut_data <- SC2db %>%
  filter(grepl(paste(mutations, collapse = "|"), spike_mut)) %>%
  mutate(
    Sampledate = as.Date(prove_tatt),
    Substitution = gsub(";", ",", spike_mut),
    YearMonth = format_month_key_nb(Sampledate)
  ) %>%
  select(Sampledate, spike_mut, nc_pangolin_short, Collapsed_pango, YearMonth) %>%
  filter(Sampledate >= Sys.Date() - months(12))

for (mutation in mutations) {
  loop_started_at <- Sys.time()
  if (exists("log_timed_message", mode = "function")) {
    log_timed_message("Loop Mutation column START: ", mutation)
  }
  S_mut_data[[mutation]] <- as.integer(str_detect(S_mut_data$spike_mut, mutation))
  if (exists("log_timed_message", mode = "function")) {
    loop_elapsed <- as.numeric(difftime(Sys.time(), loop_started_at, units = "secs"))
    log_timed_message("Loop Mutation column DONE: ", mutation, " (", sprintf("%.2f", loop_elapsed), "s)")
  }
}

S_mut_data <- S_mut_data %>%
  mutate(Combination = apply(S_mut_data[, mutations], 1, function(x) paste(names(x)[x == 1], collapse = ",")))

counts <- S_mut_data %>%
  group_by(YearMonth, Combination) %>%
  summarise(Count = n(), .groups = "drop")

countspango <- S_mut_data %>%
  group_by(YearMonth, nc_pangolin_short, Combination) %>%
  summarise(Count = n(), .groups = "drop")

counts <- counts %>%
  left_join(spm_spike, by = "YearMonth") %>%
  mutate(Percentage = (Count / TotalSeq) * 100)

countm <- counts %>%
  mutate(YearMonth = parse_month_key_nb(YearMonth))

hmapmut <- make_monthly_heatmap(
  data = countm,
  date_col = "YearMonth",
  category_col = "Combination",
  value_col = "Percentage",
  fill_label = "Andel (%)",
  x_breaks = "1 month",
  base_y_text_size = 8
) +
  labs(
    title = "Mutasjonskombinasjoner per måned",
    y = "Mutasjonskombinasjon"
  )

export_graph <- export_to_ppt(
  export_graph,
  hmapmut,
  "Mutasjonskombinasjoner per måned",
  build_slide_subtitle(countm, "YearMonth", sum(counts$Count, na.rm = TRUE), "Siste 12 måneder | Andel (%)")
)

countspango <- countspango %>%
  mutate(Sampledate = parse_month_key_nb(YearMonth))

unique_combinations <- unique(countspango$Combination)

for (combination in unique_combinations) {
  loop_started_at <- Sys.time()
  if (exists("log_timed_message", mode = "function")) {
    log_timed_message("Loop Combination START: ", combination)
  }

  subset_data <- countspango %>%
    filter(Combination == combination)

  collapsed_pangosrecgr <- make_monthly_heatmap(
    data = subset_data,
    date_col = "Sampledate",
    category_col = "nc_pangolin_short",
    value_col = "Count",
    fill_label = "Antall (n)",
    x_breaks = "1 month",
    base_y_text_size = 8
  ) +
    labs(title = paste("Pangolin per kombinasjon", combination))

  export_graph <- export_to_ppt(
    export_graph,
    collapsed_pangosrecgr,
    paste("Kombinasjon", combination, "antall"),
    build_slide_subtitle(subset_data, "Sampledate", sum(subset_data$Count, na.rm = TRUE), "Heatmap av pangolinfordeling per kombinasjon")
  )

  subset_data <- subset_data %>%
    group_by(Sampledate) %>%
    mutate(Percent = Count / sum(Count) * 100) %>%
    ungroup()

  collapsed_pangosrecgr_percent <- make_monthly_heatmap(
    data = subset_data,
    date_col = "Sampledate",
    category_col = "nc_pangolin_short",
    value_col = "Percent",
    fill_label = "Andel (%)",
    x_breaks = "1 month",
    base_y_text_size = 8
  ) +
    labs(title = paste("Andel av pangolin per kombinasjon", combination))

  export_graph <- export_to_ppt(
    export_graph,
    collapsed_pangosrecgr_percent,
    paste("Kombinasjon", combination, "andel"),
    build_slide_subtitle(subset_data, "Sampledate", sum(subset_data$Count, na.rm = TRUE), "Heatmap av pangolinandeler per kombinasjon")
  )

  if (exists("log_timed_message", mode = "function")) {
    loop_elapsed <- as.numeric(difftime(Sys.time(), loop_started_at, units = "secs"))
    log_timed_message("Loop Combination DONE: ", combination, " (", sprintf("%.2f", loop_elapsed), "s)")
  }
}

# ---- END INLINED: SC2/SC2_spike_mut_of_interest.R ----


# ============================================================================ 
# MUTATION ANALYSIS - PANGOLIN FOCUS
# ============================================================================

# ---- BEGIN INLINED: SC2/SC2_spike_mut_freq.R ----
# --- Filter and Prepare Linmut Data for Mutation Analysis ---

Linmut <- SC2db %>%
  # Select relevant columns for mutation analysis
  select(spike_mut, my, year, week, nc_pangolin_short, nc_pangolin_long, Collapsed_pango, Tessy, key, prove_tatt) %>%
  # Filter to include data within the last year
  filter(prove_tatt >= Sys.Date() - 365) %>%
  # Split mutation data into separate entries per substitution
  mutate(Substitution = str_split(spike_mut, ";|,", simplify = FALSE)) %>%
  unnest(Substitution) %>%
  # Remove redundant columns and empty substitutions
  select(-spike_mut) %>%
  mutate(Substitution = stringr::str_trim(as.character(Substitution))) %>%
  filter(Substitution != "")

# --- Count Weekly Mutations and Calculate Percentages ---

spikecount <- Linmut %>%
  count(Substitution, my, name = "count") %>%
  ungroup() %>%
  # Join to include total sequences per month
  left_join(v_seqs_per_month, by = "my") %>%
  # Calculate mutation percentage and create sample date for plotting
  mutate(
    Percent = count / TotalSeq,
    Sampledate = parse_month_key_nb(my)
  ) %>%
  filter(Sampledate >= Sys.Date() - 365)

# --- Count Mutations per Uke by Pangolin Variant ---

spikecountcpango <- Linmut %>%
  count(Substitution, my, Collapsed_pango, name = "count") %>%
  ungroup() %>%
  # Join to get total sequence counts by variant and time period
  left_join(
    SC2db %>%
      count(my, Collapsed_pango, name = "Total") %>%
      ungroup(),
    by = c("my", "Collapsed_pango")
  ) %>%
  # Calculate percentage and format date
  mutate(
    Percent = count / Total,
    Sampledate = parse_month_key_nb(my)
  ) %>%
  filter(Sampledate >= Sys.Date() - 180)

# --- Identify Unique Pangolin Variants of Interest ---

unique_collapsed_pangos <- spikecountcpango %>%
  # Filter variants with mutation percentages within specified range in the last 3 months
  filter(Percent > 0.10, Percent < 0.95, Sampledate >= Sys.Date() - 90) %>%
  distinct(Collapsed_pango)

# --- Generate and Save Mutation Line Plots by Pangolin Variant ---

for (collapsed_pango in unique_collapsed_pangos$Collapsed_pango) {
  loop_started_at <- Sys.time()
  if (exists("log_timed_message", mode = "function")) {
    log_timed_message("Loop SpikeFreq START: ", collapsed_pango)
  }
  # Filter data for the current variant
  plot_data <- spikecountcpango %>%
    filter(Collapsed_pango == collapsed_pango)

  # Generate line plot of mutation percentage over time
  plot <- ggplot(plot_data, aes(x = Sampledate, y = Percent, group = Substitution, colour = Substitution)) +
    geom_line(linewidth = 1) +
    scale_x_date(date_breaks = "1 month", labels = format_month_label, expand = c(0, 0)) +
    ylab("Andel (%)") +
    xlab("M\u00e5ned") +
    theme_classic() +
    theme(
      plot.title = element_text(color = "grey20", size = 20, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0.5, face = "plain"),
      axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0.5, face = "bold"),
      axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "bold")
    ) +
    scale_colour_manual(
      values = fhi_discrete_palette(n_distinct(plot_data$Substitution), sc2_palette),
      guide = "none"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    # Add labels for recent mutations of interest
    geom_text_repel(
      data = subset(plot_data, Sampledate == max(Sampledate) & Percent > 0.1 & Percent < 0.95),
      aes(label = sprintf("%s %s", scales::percent(Percent), toTitleCase(Substitution)), color = Substitution),
      direction = "y",
      force = 3,
      nudge_x = 70,
      na.rm = TRUE,
      segment.size = 0.2,
      segment.linetype = 2,
      segment.angle = 0,
      min.segment.length = 0,
      box.padding = 2,
      max.overlaps = Inf,
      show.legend = FALSE
    )

  # Use save_plot function to save the plot and insert into PowerPoint
  export_graph <- save_plot(plot, paste("Spike-mutasjoner -", collapsed_pango), export_graph)
  if (exists("log_timed_message", mode = "function")) {
    loop_elapsed <- as.numeric(difftime(Sys.time(), loop_started_at, units = "secs"))
    log_timed_message("Loop SpikeFreq DONE: ", collapsed_pango, " (", sprintf("%.2f", loop_elapsed), "s)")
  }
}

# --- Data Preparation for Latest Mutation Composition Analysis ---

# Convert 'my' column to date format
Linmut <- Linmut %>%
  mutate(my = parse_month_key_nb(my))

# Filter mutations with significant percentages in the latest time point
mut_int <- spikecount %>%
  filter(my == max(my), Percent > 0.05, Percent < 0.98) %>%
  select(Substitution)

# --- Create and Tree chart  for Each Significant Mutation ---


plots <- lapply(mut_int$Substitution, function(mut) {
  # Filter data for each mutation and count occurrences by Pangolin lineage
  x_i <- Linmut %>%
    filter(str_detect(Substitution, mut), my >= Sys.Date() - months(3)) %>%
    count(my, nc_pangolin_short, Substitution)

  # Create a tree map for the mutation composition
  tm <- ggplot(x_i, aes(area = n, fill = nc_pangolin_short, label = nc_pangolin_short)) +
    geom_treemap() +
    geom_treemap_text(color = "white", place = "centre", grow = TRUE) +
    ggtitle(mut) +
    theme(legend.position = "none") # Suppress the legend display

  tm
})

# Define the layout of the plots on each slide
num_plots <- length(plots)
plots_per_slide <- 2 * 2 # Number of plots per slide
num_slides <- ceiling(num_plots / plots_per_slide)

# Loop to Save Tree Map Grids Across Multiple Slides
for (slide_index in 1:num_slides) {
  loop_started_at <- Sys.time()
  if (exists("log_timed_message", mode = "function")) {
    log_timed_message("Loop Treemap slide START: ", slide_index, "/", num_slides)
  }
  # Identify the range of plots for this slide
  start_plot <- (slide_index - 1) * plots_per_slide + 1
  end_plot <- min(start_plot + plots_per_slide - 1, num_plots)

  # Subset the list of plots for the current slide
  plot_subset <- plots[start_plot:end_plot]

  # Arrange selected plots into a grid
  combined_plot <- plot_grid(plotlist = plot_subset, ncol = 2)

  # Use save_plot function to save combined plot and insert into PowerPoint
  export_graph <- save_plot(combined_plot, paste("Mutation Composition - Slide", slide_index), export_graph)
  if (exists("log_timed_message", mode = "function")) {
    loop_elapsed <- as.numeric(difftime(Sys.time(), loop_started_at, units = "secs"))
    log_timed_message("Loop Treemap slide DONE: ", slide_index, "/", num_slides, " (", sprintf("%.2f", loop_elapsed), "s)")
  }
}


# --- Spike lollipop map: mutation position + domain + Collapsed_pango ---

Linmut <- Linmut %>%
  mutate(Number = as.integer(gsub("\\D", "", Substitution))) %>%
  mutate(Sampledate = parse_month_key_nb(my)) %>%
  filter(!is.na(Number), Number >= 1, Number <= 1273) %>%
  mutate(Domain = case_when(
    Number < 13 ~ "SP",
    Number >= 13 & Number <= 205 ~ "NTD",
    Number >= 319 & Number <= 541 ~ "RBD",
    Number >= 788 & Number <= 806 ~ "FP",
    Number >= 912 & Number <= 984 ~ "HR1",
    Number >= 1163 & Number <= 1213 ~ "HR2",
    Number >= 1214 & Number <= 1237 ~ "TM",
    Number >= 1238 & Number <= 1273 ~ "CT",
    TRUE ~ "Andre"
  ))

domainmutcp <- Linmut %>%
  mutate(
    Tessy_group = as.character(Tessy),
    Tessy_group = ifelse(is.na(Tessy_group) | trimws(Tessy_group) == "", "Ukjent", Tessy_group)
  ) %>%
  group_by(Tessy_group, Substitution, Number, Domain) %>%
  summarise(n = n(), .groups = "drop")

# Keep only Tessy groups with observed mutations for this view.
domainmutcp <- domainmutcp %>%
  filter(!is.na(Tessy_group), trimws(as.character(Tessy_group)) != "")

if (nrow(domainmutcp) == 0) {
  message("No valid Spike mutation positions found for domain lollipop map.")
} else {
  domain_df <- tibble::tribble(
    ~Domain, ~xmin, ~xmax,
    "SP", 1, 12,
    "NTD", 13, 205,
    "RBD", 319, 541,
    "FP", 788, 806,
    "HR1", 912, 984,
    "HR2", 1163, 1213,
    "TM", 1214, 1237,
    "CT", 1238, 1273
  )

  domain_fill <- stats::setNames(kvalitativ_comb[seq_len(8)], c("SP", "NTD", "RBD", "FP", "HR1", "HR2", "TM", "CT"))

  tessy_levels <- domainmutcp %>%
    group_by(Tessy_group) %>%
    summarise(total_n = sum(n), .groups = "drop") %>%
    arrange(desc(total_n), Tessy_group) %>%
    pull(Tessy_group)
  domainmutcp$Tessy_group <- factor(domainmutcp$Tessy_group, levels = tessy_levels)

  for (tessy_name in tessy_levels) {
    tessy_df <- domainmutcp %>%
      filter(Tessy_group == tessy_name) %>%
      arrange(Number, Substitution)

    if (nrow(tessy_df) == 0) {
      next
    }

    # Plot in smaller AA windows to keep all labels readable and avoid crowding.
    label_df <- tessy_df %>%
      mutate(label_rank = row_number()) %>%
      mutate(
        label_band = ifelse(label_rank %% 2 == 0, "top", "bottom"),
        nudge_y = ifelse(
          label_band == "top",
          0.16 + (label_rank %% 20) * 0.014,
          -0.16 - (label_rank %% 20) * 0.014
        )
      )

    spike_lollipop_tessy <- ggplot(tessy_df, aes(x = Number, y = 1)) +
      geom_rect(
        data = domain_df,
        aes(xmin = xmin, xmax = xmax, ymin = 0.86, ymax = 1.14, fill = Domain),
        inherit.aes = FALSE,
        alpha = 0.25,
        color = NA
      ) +
      geom_segment(aes(xend = Number, yend = 0.86), linewidth = 0.3, alpha = 0.45, color = "grey35") +
      geom_point(aes(size = n, color = Domain), alpha = 0.9) +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = Number, y = 1, label = Substitution),
        color = "black",
        size = 2.7,
        direction = "both",
        nudge_y = label_df$nudge_y,
        min.segment.length = 0,
        seed = 2526,
        max.overlaps = Inf,
        force = 28,
        force_pull = 1.2,
        box.padding = 0.3,
        point.padding = 0.2
      ) +
      scale_x_continuous(
        limits = c(1, 1273),
        breaks = c(1, 13, 205, 319, 541, 788, 806, 912, 984, 1163, 1213, 1237, 1273)
      ) +
      scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL) +
      scale_fill_manual(values = domain_fill) +
      scale_color_manual(values = domain_fill) +
      scale_size_continuous(range = c(1.6, 6.5)) +
      guides(fill = "none") +
      labs(
        title = paste0("Spike mutation positions - Tessy: ", as.character(tessy_name)),
        subtitle = "Full Spike protein | all labels shown",
        x = "Spike amino-acid position",
        y = NULL,
        color = "Spike domain",
        fill = "Spike domain",
        size = "Antall (n)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.position = "right"
      )

    export_graph <- save_plot(
      spike_lollipop_tessy,
      paste0("Spike Mutation Domain Lollipop - Tessy ", as.character(tessy_name)),
      export_graph
    )
  }
}


# ---- END INLINED: SC2/SC2_spike_mut_freq.R ----


# ============================================================================
# FRAMESHIFT / INSERTION / DELETION ANALYSIS
# Output in this section:
#   - Per-gene indel heatmaps by Tessy
# ============================================================================

export_graph <- add_section_slide(
  export_graph,
  "Indeler",
  "Frameshift, insersjoner og delesjoner per gen og Tessy",
  c("Frameshift", "Delesjoner", "Insersjoner")
)

sc2_indel_source <- if (exists("SC2db")) SC2db else SC2db
sc2_indel_date_col <- intersect(c("prove_tatt", "PROVE_TATT", "sample_date", "Sampledate"), names(sc2_indel_source))[1]
sc2_indel_cols <- names(sc2_indel_source)[grepl("(frameshift|insertion|deletion)", names(sc2_indel_source), ignore.case = TRUE)]
sc2_tessy_col <- intersect(c("Tessy", "tessy"), names(sc2_indel_source))[1]

if (!is.na(sc2_indel_date_col) && length(sc2_indel_cols) > 0 && !is.na(sc2_tessy_col)) {
  sc2_indel_df <- sc2_indel_source %>%
    mutate(
      indel_plot_date = as.Date(.data[[sc2_indel_date_col]]),
      indel_month = floor_date(indel_plot_date, "month"),
      Tessy_group = as.character(.data[[sc2_tessy_col]])
    ) %>%
    filter(!is.na(indel_month), !is.na(Tessy_group), trimws(Tessy_group) != "", Tessy_group != "Ukjent")

  sc2_indel_df <- sc2_indel_df %>%
    mutate(Tessy_group = trimws(as.character(Tessy_group))) %>%
    filter(indel_plot_date >= (Sys.Date() %m-% months(6)))

  sc2_long <- sc2_indel_df %>%
    pivot_longer(cols = all_of(sc2_indel_cols), names_to = "mutation_col", values_to = "mutation_raw") %>%
    filter(!is.na(mutation_raw), trimws(as.character(mutation_raw)) != "") %>%
    separate_rows(mutation_raw, sep = ";|,") %>%
    mutate(
      mutation_raw = trimws(as.character(mutation_raw)),
      mutation_type = case_when(
        grepl("frameshift", mutation_col, ignore.case = TRUE) ~ "Frameshift",
        grepl("insertion", mutation_col, ignore.case = TRUE) ~ "Insertion",
        grepl("deletion", mutation_col, ignore.case = TRUE) ~ "Deletion",
        TRUE ~ "Andre"
      ),
      mutation_gene = sub("^(nc_[^_]+)_.*$", "\\1", mutation_col)
    ) %>%
    filter(
      mutation_raw != "",
      !tolower(mutation_raw) %in% c("na", "n/a", "none", "no mutations", "ikke_satt")
    )

  # Correct denominator: number of samples per month/Tessy in source data.
  # Use unique sample keys when available; otherwise fall back to row count.
  if ("key" %in% names(sc2_indel_df)) {
    sc2_month_totals <- sc2_indel_df %>%
      mutate(sample_key = as.character(key)) %>%
      mutate(sample_key = ifelse(is.na(sample_key) | trimws(sample_key) == "", NA_character_, sample_key)) %>%
      group_by(indel_month, Tessy_group) %>%
      summarise(
        total = ifelse(sum(!is.na(sample_key)) > 0, dplyr::n_distinct(sample_key, na.rm = TRUE), dplyr::n()),
        .groups = "drop"
      )
  } else {
    sc2_month_totals <- sc2_indel_df %>%
      count(indel_month, Tessy_group, name = "total")
  }

  mut_counts <- sc2_long %>%
    group_by(indel_month, Tessy_group, mutation_type, mutation_gene, mutation_raw) %>%
    summarise(n = n(), .groups = "drop") %>%
    left_join(sc2_month_totals, by = c("indel_month", "Tessy_group")) %>%
    mutate(percent = (n / total) * 100)

  if (nrow(mut_counts) > 0) {
    mutation_type_order <- c("Frameshift", "Deletion", "Insertion")
    for (m_type in mutation_type_order) {
      type_df <- mut_counts %>%
        filter(mutation_type == m_type)
      if (nrow(type_df) == 0) next

      gene_order <- type_df %>%
        count(mutation_gene, wt = n, name = "total_n") %>%
        arrange(desc(total_n), mutation_gene) %>%
        pull(mutation_gene)

      for (gene_name in gene_order) {
        gene_df <- type_df %>%
          filter(mutation_gene == gene_name) %>%
          mutate(
            mutation_raw = as.character(mutation_raw),
            mutation_raw = factor(mutation_raw, levels = rev(unique(mutation_raw[order(percent, na.last = TRUE)])))
          )
        if (nrow(gene_df) == 0) next

        y_text_size <- heatmap_axis_text_size(
          nlevels(gene_df$mutation_raw),
          base_size = 8,
          min_size = 5
        )

        indel_heatmap <- ggplot(gene_df, aes(x = indel_month, y = mutation_raw, fill = percent)) +
          geom_tile(color = "white") +
          facet_wrap(~Tessy_group, scales = "free_y", ncol = 3) +
          scale_fill_gradientn(colors = kvantitativ_b1, labels = scales::percent_format(scale = 1)) +
          scale_x_date(labels = format_month_label, date_breaks = "1 month") +
          labs(
            title = paste0("SC2 ", m_type, " - ", gene_name, " siste 6 måneder"),
            subtitle = "Fasettert etter Tessy",
            x = "",
            y = "Mutasjon",
            fill = "Andel (%)"
          ) +
          theme_minimal(base_size = 11) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = y_text_size)
          )

        export_graph <- export_to_ppt(
          export_graph,
          indel_heatmap,
          paste0("SC2 ", m_type, " - ", gene_name, " per Tessy")
        )
      }
    }
  } else {
    message("Skipping SC2 indel split plots: no plottable mutation rows.")
  }
} else {
  message("Skipping SC2 indel analysis: missing date/Tessy column or no frameshift/insertion/deletion columns found.")
}

# ============================================================================
# ORF ANALYSIS - ALL GENES
# ============================================================================


# ============================================================================
# DRUG RESISTANCE ANALYSIS
# Output in this section:
#   - 3CLpro resistance table
#   - RdRP resistance table
#   - Spike resistance table
# ============================================================================

export_graph <- add_section_slide(
  export_graph,
  "Resistens",
  "Resistensmutasjoner",
  c("3CLpro", "RdRP", "Spike")
)

# Paxlovid resistance
dr_pax <- SC2db %>%
  group_by(Tessy, dr_3c_lpro_mut, dr_res_paxlovid, dr_3c_lpro_fold) %>%
  filter(!is.na(dr_3c_lpro_mut) &
    dr_3c_lpro_mut != "NA" &
    dr_3c_lpro_mut != "No Mutations" &
    dr_3c_lpro_mut != "" &
    Tessy != "") %>%
  count() %>%
  ungroup() %>%
  as.data.frame()

# Remdesivir resistance
dr_remd <- SC2db %>%
  group_by(Tessy, dr_rd_rp_mut, dr_res_remdesevir, dr_rd_rp_fold) %>%
  filter(!is.na(dr_rd_rp_mut) &
    dr_rd_rp_mut != "NA" &
    dr_rd_rp_mut != "No Mutations" &
    dr_rd_rp_mut != "" &
    Tessy != "") %>%
  count() %>%
  ungroup() %>%
  as.data.frame()

# Antibody resistance
dr_ab <- SC2db %>%
  group_by(Tessy, dr_spike_m_abs, dr_spike_m_abs_fold) %>%
  filter(!is.na(dr_3c_lpro_mut) &
    dr_spike_m_abs != "NA" &
    dr_spike_m_abs != "No Mutations" &
    dr_spike_m_abs != "" &
    Tessy != "") %>%
  count() %>%
  ungroup() %>%
  as.data.frame()

table_data <- list(
  list(data = dr_pax, caption = "3CLpro-resistensmutasjoner"),
  list(data = dr_remd, caption = "RdRP-resistensmutasjoner"),
  list(data = dr_ab, caption = "Spike-relaterte resistensmutasjoner")
)

for (table_info in table_data) {
  export_graph <- export_to_ppt(export_graph, table_info$data, table_info$caption)
}

# ============================================================================
# AGE DISTRIBUTION BY TESSY (3M/6M)
# ============================================================================

age_source_df <- if (exists("SC2db")) SC2db else SC2db

age_date_col <- intersect(c("prove_tatt", "PROVE_TATT", "sample_date", "Sampledate"), names(age_source_df))[1]
age_tessy_col <- intersect(c("Tessy", "tessy"), names(age_source_df))[1]
age_col <- intersect(c("pasient_alder"), names(age_source_df))[1]

build_age_tessy_plot <- function(month_window, sentinel_only = FALSE) {
  if (is.na(age_date_col) || is.na(age_tessy_col)) {
    return(NULL)
  }
  if (is.na(age_col)) {
    return(NULL)
  }

  plot_df <- age_source_df %>%
    mutate(
      plot_date = as.Date(.data[[age_date_col]]),
      Tessy_plot = as.character(.data[[age_tessy_col]])
    ) %>%
    filter(
      !is.na(plot_date),
      plot_date >= (Sys.Date() %m-% months(month_window)),
      !is.na(Tessy_plot),
      trimws(Tessy_plot) != ""
    )

  if (sentinel_only) {
    if (!("prove_kategori_group" %in% names(plot_df))) {
      return(NULL)
    }
    plot_df <- plot_df %>% filter(prove_kategori_group == "Sentinel")
  }

  plot_df <- plot_df %>%
    mutate(
      pasient_alder_num = suppressWarnings(as.numeric(trimws(as.character(.data[[age_col]])))),
      age_group_raw = as.character(pasient_alder_num),
      pasient_aldersgruppe = as.character(age_to_group_standard(pasient_alder_num)),
      age_group_plot = case_when(
        pasient_aldersgruppe %in% c("0-4", "5-14", "15-24", "25-59", "60+") ~ pasient_aldersgruppe,
        TRUE ~ NA_character_
      )
    )

  age_levels <- c("0-4", "5-14", "15-24", "25-59", "60+")

  plot_df <- plot_df %>%
    mutate(age_group_plot = factor(age_group_plot, levels = age_levels)) %>%
    filter(!is.na(age_group_plot))

  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  age_tessy_df <- plot_df %>%
    count(Tessy_plot, age_group_plot, name = "n") %>%
    group_by(Tessy_plot) %>%
    mutate(
      tessy_n = sum(n),
      percent = (n / tessy_n) * 100
    ) %>%
    ungroup()

  tessy_levels <- age_tessy_df %>%
    distinct(Tessy_plot, tessy_n) %>%
    arrange(desc(tessy_n), Tessy_plot) %>%
    pull(Tessy_plot)

  tessy_labels <- age_tessy_df %>%
    distinct(Tessy_plot, tessy_n) %>%
    mutate(label = paste0(Tessy_plot, "\n(n=", tessy_n, ")")) %>%
    {
      setNames(.$label, .$Tessy_plot)
    }

  age_tessy_df <- age_tessy_df %>%
    mutate(Tessy_plot = factor(Tessy_plot, levels = tessy_levels))

  # Diagnostics for last 6 months: raw labels and normalized labels
  if (month_window == 6 && !sentinel_only) {
    current_week_age <- week(Sys.Date())
    current_year_age <- year(Sys.Date())
    results_dir_age <- results_stats_dir

    raw_age_counts <- plot_df %>%
      count(Tessy_plot, age_group_raw, name = "n_raw") %>%
      arrange(desc(n_raw))

    norm_age_counts <- plot_df %>%
      count(Tessy_plot, age_group_plot, name = "n_norm") %>%
      arrange(Tessy_plot, age_group_plot)
  }

  ggplot(age_tessy_df, aes(x = Tessy_plot, y = percent, fill = age_group_plot)) +
    geom_col(position = "stack") +
    scale_x_discrete(labels = tessy_labels) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    coord_cartesian(ylim = c(0, 100)) +
    scale_fill_manual(values = kvalitativ_a) +
    labs(
      title = paste0(
        "Aldersgruppefordeling per ECDC-variantklassifisering",
        ifelse(sentinel_only, " (Sentinel only)", ""),
        " (last ", month_window, " months)"
      ),
      x = "ECDC Variant Classification",
      y = "Andel (%)",
      fill = "Age group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.major.x = element_blank()
    )
}

## The age-by-Tessy plot builder above is kept for ad hoc diagnostics,
## but these exploratory slides are not part of the standard weekly report.

# ============================================================================
# PATIENT AND TESSY DISTRIBUTIONS
# Output in this section:
#   - Sex, age and vaccination summaries
#   - Current-season Tessy distributions by geography and patient group
# ============================================================================
export_graph <- add_section_slide(export_graph, "Pasientdata", "Demografi, geografi og kliniske fordelinger i overv\u00e5king", c("Kj\u00f8nn", "Alder", "Vaksinasjon"))

patient_source_df <- if (exists("SC2db")) SC2db else SC2db
patient_tessy_col <- intersect(c("Tessy", "tessy"), names(patient_source_df))[1]
patient_date_col <- intersect(c("prove_tatt", "PROVE_TATT", "sample_date", "Sampledate"), names(patient_source_df))[1]
patient_fylke_col <- intersect(c("pasient_fylke_name", "fylkenavn", "pasient_fylke"), names(patient_source_df))[1]
patient_landsdel_col <- intersect(c("pasient_landsdel", "landsdel"), names(patient_source_df))[1]
patient_age_col <- intersect(c("pasient_alder"), names(patient_source_df))[1]

# Build paired percentage/count slides for one current-season Tessy grouping.
build_tessy_group_plots <- function(df, x_col, x_label, export_graph_in) {
  if (is.na(patient_tessy_col) || is.na(patient_date_col) || is.na(x_col)) {
    return(export_graph_in)
  }

  plot_df <- df %>%
    mutate(
      plot_date = as.Date(.data[[patient_date_col]]),
      season_plot = season_label_from_date(plot_date),
      Tessy_plot = as.character(.data[[patient_tessy_col]]),
      group_plot = as.character(.data[[x_col]])
    ) %>%
    filter(
      !is.na(plot_date),
      season_plot == current_season_label,
      !is.na(Tessy_plot), trimws(Tessy_plot) != "",
      !is.na(group_plot), trimws(group_plot) != "",
      group_plot != "IKKE_SATT"
    )

  if (nrow(plot_df) == 0) {
    return(export_graph_in)
  }

  grouped_df <- plot_df %>%
    count(group_plot, Tessy_plot, name = "n") %>%
    group_by(group_plot) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    ungroup()

  x_labels_df <- grouped_df %>%
    group_by(group_plot) %>%
    summarise(group_n = sum(n), .groups = "drop") %>%
    arrange(desc(group_n), group_plot) %>%
    mutate(group_label = paste0(group_plot, " (n=", group_n, ")"))

  grouped_df <- grouped_df %>%
    left_join(x_labels_df, by = "group_plot") %>%
    mutate(group_label = factor(group_label, levels = x_labels_df$group_label))

  p_pct <- ggplot(grouped_df, aes(x = group_label, y = percent, fill = Tessy_plot)) +
    geom_col(position = "stack") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    coord_cartesian(ylim = c(0, 100)) +
    scale_fill_manual(values = fhi_discrete_palette(n_distinct(grouped_df$Tessy_plot), sc2_palette)) +
    labs(title = paste0("Tessy-fordeling per ", x_label, " (%) - gjeldende sesong"), x = x_label, y = "Andel (%)", fill = "Tessy") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_count <- ggplot(grouped_df, aes(x = group_label, y = n, fill = Tessy_plot)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = fhi_discrete_palette(n_distinct(grouped_df$Tessy_plot), sc2_palette)) +
    labs(title = paste0("Tessy-fordeling per ", x_label, " (antall) - gjeldende sesong"), x = x_label, y = "Antall (n)", fill = "Tessy") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  export_graph_out <- export_to_ppt(export_graph_in, p_pct, paste0("Tessy per ", x_label, " (%) - gjeldende sesong"))
  export_graph_out <- export_to_ppt(export_graph_out, p_count, paste0("Tessy per ", x_label, " (antall) - gjeldende sesong"))
  export_graph_out
}

if (!is.na(patient_tessy_col)) {
  patient_df <- patient_source_df
  if (!is.na(patient_age_col)) {
    patient_df <- patient_df %>%
      mutate(
        age_value_pd = suppressWarnings(as.numeric(trimws(as.character(.data[[patient_age_col]])))),
        pasient_aldersgruppe_sc2 = as.character(age_to_group_standard(age_value_pd))
      )
  } else {
    patient_df$pasient_aldersgruppe_sc2 <- NA_character_
  }
  patient_df <- normalize_sex_column(patient_df, candidate_cols = c("pasient_kjnn", "pasient_kjonn"))

  # Normalize patient-facing dimensions once, then reuse them across the patient section outputs.
  patient_df <- patient_df %>%
    mutate(
      pasient_status_plot = if ("pasient_status" %in% names(.)) ifelse(is.na(pasient_status) | pasient_status == "" | pasient_status == "IKKE_SATT", "Ukjent", as.character(pasient_status)) else "Ukjent",
      prove_kategori_group_plot = if ("prove_kategori_group" %in% names(.)) as.character(prove_kategori_group) else classify_prove_kategori_group(prove_kategori),
      pasient_vaks_plot = if ("pasient_vaks" %in% names(.)) ifelse(is.na(pasient_vaks) | trimws(as.character(pasient_vaks)) == "", "Ukjent", as.character(pasient_vaks)) else "Ukjent",
      pasient_vaks_2uipt_plot = if ("pasient_vaks_2uipt" %in% names(.)) ifelse(is.na(pasient_vaks_2uipt) | trimws(as.character(pasient_vaks_2uipt)) == "", "Ukjent", as.character(pasient_vaks_2uipt)) else "Ukjent"
    )

  norway_geojson_path <- resolve_norway_geojson_path()
  sc2_prev <- patient_df %>% dplyr::filter(season == previous_season_label)
  sc2_curr <- patient_df %>% dplyr::filter(season == current_season_label)
  if (!is.na(patient_fylke_col)) {
    p_fylke_prev <- build_fylke_map_plot_shared(
      sc2_prev,
      fylke_col = patient_fylke_col,
      shape_path = norway_geojson_path,
      fill_palette = kvantitativ_b2
    )
    p_fylke_curr <- build_fylke_map_plot_shared(
      sc2_curr,
      fylke_col = patient_fylke_col,
      shape_path = norway_geojson_path,
      fill_palette = kvantitativ_b2
    )
    if (!is.null(p_fylke_curr) && !is.null(p_fylke_prev)) {
      p_fylke_pair <- build_two_season_map_compare_shared(
        sc2_prev,
        sc2_curr,
        previous_label = previous_season_label,
        current_label = current_season_label,
        map_builder = build_fylke_map_plot_shared,
        fylke_col = patient_fylke_col,
        shape_path = norway_geojson_path,
        fill_palette = kvantitativ_b2
      )
      export_graph <- export_to_ppt(export_graph, p_fylke_pair, "Fylkesfordeling: forrige sesong vs. gjeldende sesong", preserve_plot_labels = TRUE)
    }
  }

  if (all(c("pasient_kjonn_std", "season") %in% names(patient_df))) {
    p_kjonn <- build_two_season_pie_compare(
      patient_df,
      season_col = "season",
      category_col = "pasient_kjonn_std",
      previous_label = previous_season_label,
      current_label = current_season_label,
      category_label = "Kj\u00f8nn",
      palette_base = sc2_palette
    )
    if (!is.null(p_kjonn)) {
      export_graph <- export_to_ppt(export_graph, p_kjonn, "Kj\u00f8nn: sesongsammenligning", preserve_plot_labels = TRUE)
    }
  }

  if ("pasient_aldersgruppe_sc2" %in% names(patient_df)) {
    p_alder <- build_two_season_pie_compare(
      patient_df,
      season_col = "season",
      category_col = "pasient_aldersgruppe_sc2",
      previous_label = previous_season_label,
      current_label = current_season_label,
      category_label = "Aldersgruppe",
      palette_base = sc2_palette
    )
    if (!is.null(p_alder)) {
      export_graph <- export_to_ppt(export_graph, p_alder, "Aldersgruppe: sesongsammenligning", preserve_plot_labels = TRUE)
    }
  }

  export_graph <- build_tessy_group_plots(patient_df, patient_fylke_col, "Fylke", export_graph)
  export_graph <- build_tessy_group_plots(patient_df, patient_landsdel_col, "Landsdel", export_graph)
  export_graph <- build_tessy_group_plots(patient_df, "pasient_aldersgruppe_sc2", "Aldersgruppe", export_graph)
  export_graph <- build_tessy_group_plots(patient_df, "pasient_status_plot", "Pasientstatus", export_graph)
  export_graph <- build_tessy_group_plots(patient_df, "prove_kategori_group_plot", "Pr\u00f8vekategorigruppe", export_graph)
  export_graph <- build_tessy_group_plots(patient_df, "pasient_vaks_plot", "Vaksinasjonsstatus", export_graph)
  export_graph <- build_tessy_group_plots(patient_df, "pasient_vaks_2uipt_plot", "Vaksinert <=2 uker f\u00f8r symptomstart", export_graph)
} else {
  message("Could not create patient Tessy plots (missing Tessy column).")
}

# ============================================================================
# PRINT AND SAVE FILES/GRAPHS
# ============================================================================

# Get the current week and year
current_week <- week(Sys.Date())
current_year <- year(Sys.Date())

# Create the file names with the current week and year
file_name_result <- paste0("SARSCOV2_", current_year, "_Week", current_week, "_result.pptx")
file_name_statistikk <- paste0("SARSCOV2_", current_year, "_Week", current_week, "_statistikk.csv")
file_name_tabeller <- paste0("SARSCOV2_", current_year, "_Week", current_week, "_tabeller.xlsx")

# Specify the full file paths
file_path_result <- file.path(results_root, file_name_result)
file_path_resultshare <- file.path(results_share_dir, file_name_result)
file_path_statistikk <- file.path(results_stats_dir, file_name_statistikk)

write_ppt_safe <- function(ppt_obj, target_path) {
  tryCatch(
    {
      print(ppt_obj, target = target_path)
      target_path
    },
    error = function(e) {
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
    }
  )
}

# Print outputs
result_written_path <- timed_step("Write PPTX to Results", write_ppt_safe(export_graph, file_path_result))
share_written_path <- timed_step("Write PPTX to OneDrive share", write_ppt_safe(export_graph, file_path_resultshare))
invisible(timed_step("Write statistikk CSV", write.csv2(final_pangostatistikk, file_path_statistikk, row.names = FALSE)))
file_path_tabeller <- file.path(results_root, file_name_tabeller)
invisible(timed_step("Write tabell-Excel", {
  wb <- openxlsx::createWorkbook()
  if (exists("pangolin_table_full")) {
    openxlsx::addWorksheet(wb, "Pangolin_prosent")
    openxlsx::writeData(wb, "Pangolin_prosent", pangolin_table_full)
  }
  if (exists("pangolin_count_table_full")) {
    openxlsx::addWorksheet(wb, "Pangolin_antall")
    openxlsx::writeData(wb, "Pangolin_antall", pangolin_count_table_full)
  }
  if (exists("tessy_table_full")) {
    openxlsx::addWorksheet(wb, "Tessy_prosent")
    openxlsx::writeData(wb, "Tessy_prosent", tessy_table_full)
  }
  openxlsx::saveWorkbook(wb, file_path_tabeller, overwrite = TRUE)
}))
log_timed_message("Result PPTX path: ", result_written_path)
log_timed_message("Share PPTX path: ", share_written_path)

total_elapsed_sec <- as.numeric(difftime(Sys.time(), analysis_started_at, units = "secs"))
log_timed_message("TOTAL RUNTIME: ", sprintf("%.2f", total_elapsed_sec), "s")

