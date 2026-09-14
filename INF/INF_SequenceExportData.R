# =============================================================================
# Influenza sequence-export data preparation
#
# Creates a sample-level data frame for a future Shiny filtering/FASTA-export
# application.  Unlike INF_DataCleaning.R, this script retains every
# stored sample in the base data.  The familiar operational exclusions are
# recorded as columns and applied only in `inf_sequence_default_view`.
#
# Objects created:
#   inf_sequence_export_data        One row per stored sample, including all
#                                   metadata, derived filter fields, segment
#                                   presence, and a representative sequence
#                                   for each influenza segment.
#   inf_sequence_records            One row per stored, decodable segment
#                                   sequence (the source for audits/exports).
#   inf_sequence_segment_summary    One longest decoded sequence per sample
#                                   and segment (the source for FASTA output).
#   inf_sequence_default_view       The usual operational view of the data.
#   inf_sequence_filter_fields      Suggested fields for Shiny filter widgets.
#   inf_sequence_filter_values      Available values for those filter fields.
#
# The script does not write files when sourced.  Call `write_inf_fasta()` with
# a filtered data frame and an explicitly chosen output directory to write one
# FASTA file per selected segment.
# =============================================================================

suppressPackageStartupMessages({
  library(base64enc)
  library(dplyr)
  library(janitor)
  library(purrr)
  library(stringr)
  library(tidyr)
})

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
  if (length(frame_paths) > 0) return(dirname(frame_paths[[1]]))

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

sequence_export_dir <- resolve_script_dir()

if (!exists("season_label_from_date", mode = "function")) {
  source(file.path(sequence_export_dir, "..", "Source_files", "common_report_utils.R"))
}

# Load the complete raw metadata before retrieving SEQUENCEDATA.  The existing
# query provides the mapped, wide metadata fields used by the analysis scripts.
if (!exists("INF_25_26_raw_merged")) {
  source(file.path(sequence_export_dir, "INF_SQLquery_25-26.R"))
}

if (!exists("conFLU2526")) {
  stop("The INF SQL connection `conFLU2526` is unavailable after loading INF_SQLquery_25-26.R.")
}

segment_levels <- c("HA", "NA", "M", "PB1", "PB2", "PA", "NP", "NS")

is_blank <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

parse_inf_date <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))

  x_chr <- trimws(as.character(x))
  parsed <- suppressWarnings(as.Date(x_chr))
  missing <- is.na(parsed) & !is.na(x_chr) & x_chr != ""
  if (any(missing)) {
    parsed[missing] <- suppressWarnings(as.Date(x_chr[missing], format = "%d.%m.%Y"))
  }
  parsed
}

clean_project_code <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_chr <- gsub("\\s+", "", x_chr)
  x_chr <- ifelse(grepl("^[0-9]+", x_chr), sub("^([0-9]+).*$", "P\\1", x_chr), x_chr)
  ifelse(grepl("^P[0-9]+", x_chr), sub("^(P[0-9]+).*$", "\\1", x_chr), NA_character_)
}

classify_prove_kategori_group <- function(x) {
  x_chr <- trimws(as.character(x))
  ifelse(grepl("^(P1\\b|P1_|1\\b)", x_chr), "Sentinel", "Non-Sentinel")
}

canonical_segment <- function(experiment) {
  experiment_chr <- toupper(trimws(as.character(experiment)))
  dplyr::case_when(
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*1[-_ ]*)?HA([^A-Z0-9]|$)") ~ "HA",
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*2[-_ ]*)?NA([^A-Z0-9]|$)") ~ "NA",
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*3[-_ ]*)?(M|MP|M1|M2)([^A-Z0-9]|$)") ~ "M",
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*4[-_ ]*)?PB1([^A-Z0-9]|$)") ~ "PB1",
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*5[-_ ]*)?PB2([^A-Z0-9]|$)") ~ "PB2",
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*7[-_ ]*)?PA([^A-Z0-9]|$)") ~ "PA",
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*6[-_ ]*)?NP([^A-Z0-9]|$)") ~ "NP",
    str_detect(experiment_chr, "(^|[^A-Z0-9])(?:0*8[-_ ]*)?NS([^A-Z0-9]|$)") ~ "NS",
    TRUE ~ NA_character_
  )
}

decode_sequence_data <- function(data_gz_base64) {
  if (is.na(data_gz_base64) || !nzchar(as.character(data_gz_base64))) {
    return(list(sequence = NA_character_, decode_error = "Missing sequence data"))
  }

  tryCatch({
    compressed_raw <- base64enc::base64decode(as.character(data_gz_base64))
    decompressed_raw <- memDecompress(compressed_raw, type = "gzip")
    raw_values <- as.integer(decompressed_raw)

    # The current database payloads start with binary technical bytes (for
    # example 165, 6, 0, 0) before the lower-case nucleotide sequence.  Decode
    # directly from raw bytes so a Windows text encoding never corrupts that
    # prefix or prevents the actual sequence from being recovered.
    if (length(raw_values) > 0 && raw_values[[1]] == utf8ToInt(">")) {
      header_end <- which(raw_values %in% c(10L, 13L))[[1]]
      if (!is.na(header_end)) raw_values <- raw_values[seq.int(header_end + 1L, length(raw_values))]
    }

    # Preserve valid IUPAC ambiguity codes and alignment gaps.  U is converted
    # to T so the result is a DNA-ready sequence for FASTA output.
    nucleotide_lookup <- c(
      `65` = "A", `67` = "C", `71` = "G", `84` = "T", `85` = "T",
      `82` = "R", `89` = "Y", `83` = "S", `87` = "W", `75` = "K",
      `77` = "M", `66` = "B", `68` = "D", `72` = "H", `86` = "V",
      `78` = "N", `45` = "-",
      `97` = "A", `99` = "C", `103` = "G", `116` = "T", `117` = "T",
      `114` = "R", `121` = "Y", `115` = "S", `119` = "W", `107` = "K",
      `109` = "M", `98` = "B", `100` = "D", `104` = "H", `118` = "V",
      `110` = "N"
    )
    sequence_chars <- unname(nucleotide_lookup[as.character(raw_values)])
    sequence <- paste0(sequence_chars[!is.na(sequence_chars)], collapse = "")
    if (!nzchar(sequence)) stop("No nucleotide characters remained after decoding")

    list(sequence = sequence, decode_error = NA_character_)
  }, error = function(e) {
    list(sequence = NA_character_, decode_error = conditionMessage(e))
  })
}

first_longest_sequence <- function(x) {
  valid <- x[!is.na(x) & nzchar(x)]
  if (length(valid) == 0) return(NA_character_)
  valid[[which.max(nchar(valid))]]
}

all_metadata <- INF_25_26_raw_merged %>%
  mutate(key = as.character(key)) %>%
  filter(!is.na(key), key != "") %>%
  distinct(key, .keep_all = TRUE)

# Only the fields needed for a later FASTA exporter are retrieved.  All stored
# SEQUENCE records are retained, including any with an unfamiliar experiment
# label, so that the latter can be audited instead of being silently discarded.
sequence_source <- tbl(conFLU2526, "SEQUENCEDATA") %>%
  select(KEY, EXPERIMENT, TYPE, DATA) %>%
  collect() %>%
  janitor::clean_names() %>%
  transmute(
    key = as.character(key),
    experiment = as.character(experiment),
    sequence_type = as.character(type),
    data = as.character(data)
  ) %>%
  filter(!is.na(key), key != "") %>%
  mutate(
    segment = canonical_segment(experiment),
    is_sequence_record = toupper(trimws(sequence_type)) == "SEQUENCE"
  )

decoded_sequence_records <- sequence_source %>%
  filter(is_sequence_record) %>%
  mutate(decoded = purrr::map(data, decode_sequence_data)) %>%
  transmute(
    key,
    experiment,
    sequence_type,
    segment,
    sequence = purrr::map_chr(decoded, "sequence"),
    decode_error = purrr::map_chr(decoded, "decode_error")
  )

# The sample universe is the union of metadata rows and SEQUENCEDATA rows.  It
# therefore includes a database sample even if either side is incomplete.
sample_keys <- tibble(
  key = sort(unique(c(all_metadata$key, sequence_source$key)))
) %>%
  filter(!is.na(key), key != "")

inf_sequence_records <- decoded_sequence_records %>%
  filter(!is.na(segment)) %>%
  mutate(segment = factor(segment, levels = segment_levels)) %>%
  arrange(key, segment, experiment)

inf_sequence_unrecognised_experiments <- decoded_sequence_records %>%
  filter(is.na(segment)) %>%
  arrange(key, experiment)

inf_sequence_decode_issues <- inf_sequence_records %>%
  filter(!is.na(decode_error)) %>%
  arrange(key, segment, experiment)

inf_sequence_segment_summary <- inf_sequence_records %>%
  group_by(key, segment) %>%
  summarise(
    n_stored_records = n(),
    n_decoded_records = sum(!is.na(sequence) & nzchar(sequence)),
    sequence = first_longest_sequence(sequence),
    .groups = "drop"
  ) %>%
  mutate(segment = as.character(segment))

segment_grid <- tidyr::crossing(sample_keys, segment = segment_levels) %>%
  left_join(inf_sequence_segment_summary, by = c("key", "segment")) %>%
  mutate(
    n_stored_records = tidyr::replace_na(n_stored_records, 0L),
    n_decoded_records = tidyr::replace_na(n_decoded_records, 0L),
    present = n_stored_records > 0L,
    decoded = n_decoded_records > 0L
  )

segment_columns <- segment_grid %>%
  select(key, segment, present, decoded, n_stored_records, n_decoded_records, sequence) %>%
  pivot_wider(
    names_from = segment,
    values_from = c(present, decoded, n_stored_records, n_decoded_records, sequence),
    names_glue = "{segment}_{.value}"
  )

presence_columns <- paste0(segment_levels, "_present")
decoded_columns <- paste0(segment_levels, "_decoded")

metadata_with_filters <- all_metadata %>%
  mutate(
    prove_tatt = if ("prove_tatt" %in% names(.)) parse_inf_date(prove_tatt) else as.Date(NA),
    season = season_label_from_date(prove_tatt),
    prove_kategori = if ("prove_kategori" %in% names(.)) as.character(prove_kategori) else NA_character_,
    prove_project = clean_project_code(prove_kategori),
    prove_kategori_group = classify_prove_kategori_group(prove_kategori),
    is_project_3 = !is.na(prove_project) & prove_project == "P3",
    ngs_report_blank = if ("ngs_report" %in% names(.)) is_blank(ngs_report) else TRUE,
    prove_kategori_reference = str_detect(coalesce(prove_kategori, ""), regex("ref", ignore_case = TRUE)),
    tessy_reportable_present = if ("tessy_reportable_variable" %in% names(.)) !is_blank(tessy_reportable_variable) else FALSE,
    tessy_reportable_reference = if ("tessy_reportable_variable" %in% names(.)) {
      str_detect(coalesce(tessy_reportable_variable, ""), regex("ref", ignore_case = TRUE))
    } else {
      FALSE
    }
  ) %>%
  mutate(
    default_visible = ngs_report_blank &
      !is_project_3 &
      !prove_kategori_reference &
      tessy_reportable_present &
      !tessy_reportable_reference,
    default_exclusion_reason = pmap_chr(
      list(ngs_report_blank, is_project_3, prove_kategori_reference, tessy_reportable_present, tessy_reportable_reference),
      function(report_blank, project_3, category_reference, tessy_present, tessy_reference) {
        reasons <- c(
          if (!report_blank) "ngs_report is populated",
          if (project_3) "prove_project is P3",
          if (category_reference) "prove_kategori is a reference",
          if (!tessy_present) "tessy_reportable_variable is blank",
          if (tessy_reference) "tessy_reportable_variable is a reference"
        )
        paste(reasons, collapse = "; ")
      }
    )
  )

inf_sequence_export_data <- sample_keys %>%
  left_join(metadata_with_filters, by = "key") %>%
  left_join(segment_columns, by = "key") %>%
  mutate(
    across(all_of(presence_columns), ~ tidyr::replace_na(.x, FALSE)),
    across(all_of(decoded_columns), ~ tidyr::replace_na(.x, FALSE)),
    sequence_segment_count = rowSums(across(all_of(presence_columns), as.integer)),
    decoded_segment_count = rowSums(across(all_of(decoded_columns), as.integer)),
    all_8_segments_present = sequence_segment_count == length(segment_levels)
  ) %>%
  rowwise() %>%
  mutate(
    segments_present = paste(segment_levels[as.logical(c_across(all_of(presence_columns)))], collapse = ";"),
    segments_decoded = paste(segment_levels[as.logical(c_across(all_of(decoded_columns)))], collapse = ";")
  ) %>%
  ungroup() %>%
  select(
    key,
    prove_tatt,
    season,
    ngs_sekvens_resultat,
    prove_kategori,
    prove_project,
    prove_kategori_group,
    default_visible,
    default_exclusion_reason,
    sequence_segment_count,
    decoded_segment_count,
    all_8_segments_present,
    segments_present,
    segments_decoded,
    everything()
  ) %>%
  arrange(desc(default_visible), desc(prove_tatt), key)

inf_sequence_default_view <- inf_sequence_export_data %>%
  filter(default_visible)

# These are intended for selectize/drop-down controls.  Sequence strings and
# internal flags are deliberately omitted; segment presence is still offered.
inf_sequence_filter_fields <- setdiff(
  names(inf_sequence_export_data),
  c(
    grep("_sequence$", names(inf_sequence_export_data), value = TRUE),
    "default_exclusion_reason"
  )
)

inf_sequence_filter_values <- lapply(
  inf_sequence_export_data[inf_sequence_filter_fields],
  function(x) {
    values <- unique(as.character(x))
    sort(values[!is.na(values) & values != ""])
  }
)

# Apply named Shiny-style filters without touching the underlying data.  For
# example: filter_inf_sequence_export(filters = list(season = "Season25_26",
# ngs_sekvens_resultat = c("A/H1N1", "A/H3N2"))).
filter_inf_sequence_export <- function(
  data = inf_sequence_export_data,
  filters = list(),
  include_default_excluded = FALSE
) {
  result <- as.data.frame(data)
  if (!include_default_excluded) {
    result <- result[!is.na(result$default_visible) & result$default_visible, , drop = FALSE]
  }

  for (field in names(filters)) {
    if (!field %in% names(result)) stop("Unknown filter field: ", field)
    selected <- as.character(filters[[field]])
    selected <- selected[!is.na(selected) & selected != ""]
    if (length(selected) > 0) {
      result <- result[as.character(result[[field]]) %in% selected, , drop = FALSE]
    }
  }
  result
}

fasta_header_value <- function(x) {
  value <- as.character(x)
  value[is.na(value) | value == ""] <- "NA"
  gsub("[[:space:]|]+", "_", value)
}

# Write the longest decoded sequence for each selected sample/segment.  Keeping
# this explicit makes the eventual download handler a thin, testable wrapper.
write_inf_fasta <- function(
  data,
  output_dir,
  segments = segment_levels,
  filename_prefix = "INF"
) {
  if (!is.data.frame(data)) stop("`data` must be a filtered sample data frame.")
  if (!"key" %in% names(data)) stop("`data` must contain a `key` column.")

  segments <- intersect(as.character(segments), segment_levels)
  if (length(segments) == 0) stop("Choose at least one valid influenza segment.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_metadata <- data %>%
    mutate(key = as.character(key)) %>%
    select(any_of(c("key", "season", "ngs_sekvens_resultat"))) %>%
    distinct(key, .keep_all = TRUE)

  output_paths <- character()
  for (segment_name in segments) {
    fasta_records <- inf_sequence_segment_summary %>%
      filter(segment == segment_name, !is.na(sequence), nzchar(sequence)) %>%
      inner_join(sample_metadata, by = "key") %>%
      arrange(key) %>%
      mutate(
        fasta_header = paste(
          fasta_header_value(key),
          segment_name,
          fasta_header_value(season),
          fasta_header_value(ngs_sekvens_resultat),
          sep = "|"
        )
      )

    if (nrow(fasta_records) == 0) next

    output_path <- file.path(output_dir, paste0(filename_prefix, "_", segment_name, ".fasta"))
    fasta_lines <- unlist(Map(function(header, sequence) c(paste0(">", header), sequence), fasta_records$fasta_header, fasta_records$sequence))
    writeLines(fasta_lines, con = output_path, useBytes = TRUE)
    output_paths <- c(output_paths, output_path)
  }
  unname(output_paths)
}

if (exists("close_sql_connections", mode = "function")) close_sql_connections()

# =============================================================================
# VCM FASTA export
#
# Exports the cleaned, GISAID-named influenza sequences for H1, H3 and BVIC.
# References are selected from prove_kategori and added separately because the
# regular cleaned view excludes reference samples.
# =============================================================================

inf_vcm_subtype_levels <- c("H1", "H3", "BVIC")

inf_vcm_fasta_header_fields <- c(
  "gisaid_isolate_name",
  "nc_ha_clade",
  "nc_ha_subclade",
  "tessy_reportable_variable",
  "prove_tatt"
)

inf_vcm_subtype <- function(ngs_sekvens_resultat, prove_kategori) {
  subtype_source <- paste(
    coalesce(as.character(ngs_sekvens_resultat), ""),
    coalesce(as.character(prove_kategori), "")
  )

  case_when(
    str_detect(subtype_source, regex("H1", ignore_case = TRUE)) ~ "H1",
    str_detect(subtype_source, regex("H3", ignore_case = TRUE)) ~ "H3",
    str_detect(subtype_source, regex("BVIC|B[/ _-]*VICTORIA", ignore_case = TRUE)) ~ "BVIC",
    TRUE ~ NA_character_
  )
}

inf_vcm_reference_flag <- function(prove_kategori) {
  str_detect(coalesce(as.character(prove_kategori), ""), regex("ref", ignore_case = TRUE))
}

inf_vcm_safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

prepare_inf_vcm_fasta_data <- function(
  data = inf_sequence_export_data,
  season = "Season25_26",
  subtypes = inf_vcm_subtype_levels
) {
  if (!is.data.frame(data)) stop("data must be a sample-level data frame.")

  required_columns <- c(
    "key",
    "season",
    "ngs_sekvens_resultat",
    "prove_kategori",
    "default_visible",
    "gisaid_isolate_name"
  )
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop("Missing required FASTA-export columns: ", paste(missing_columns, collapse = ", "))
  }

  header_missing_columns <- setdiff(inf_vcm_fasta_header_fields, names(data))
  for (field in header_missing_columns) data[[field]] <- NA_character_

  season_values <- as.character(season)
  season_values <- season_values[!is.na(season_values) & nzchar(season_values)]
  subtypes <- intersect(as.character(subtypes), inf_vcm_subtype_levels)
  if (length(subtypes) == 0) stop("Choose at least one of: ", paste(inf_vcm_subtype_levels, collapse = ", "))

  data %>%
    mutate(
      export_subtype = inf_vcm_subtype(ngs_sekvens_resultat, prove_kategori),
      reference_in_prove_kategori = inf_vcm_reference_flag(prove_kategori),
      gisaid_name_present = !is_blank(gisaid_isolate_name),
      selected_season = if (length(season_values) == 0) TRUE else .data$season %in% .env$season_values,
      cleaned_export = selected_season &
        default_visible &
        !reference_in_prove_kategori &
        gisaid_name_present,
      # Reference strains are intentionally retained across collection years.
      # They provide the phylogenetic context for the selected surveillance
      # season, and their historical `prove_tatt` dates must not exclude them.
      reference_export = reference_in_prove_kategori,
      export_selected = cleaned_export | reference_export,
      export_reason = case_when(
        reference_export ~ "reference",
        cleaned_export ~ "cleaned",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      export_selected,
      !is.na(export_subtype),
      export_subtype %in% .env$subtypes
    ) %>%
    distinct(key, .keep_all = TRUE)
}

write_inf_vcm_fasta_legacy <- function(
  data = inf_sequence_export_data,
  output_dir = Sys.getenv(
    "INF_VCM_FASTA_OUTPUT_DIR",
    unset = "N:/Virologi/Influensa/2526/WGS_Analyse/VCM 2026_27"
  ),
  season = Sys.getenv("INF_VCM_FASTA_SEASON", unset = "Season25_26"),
  subtypes = inf_vcm_subtype_levels,
  segments = segment_levels,
  header_fields = inf_vcm_fasta_header_fields
) {
  segments <- intersect(as.character(segments), segment_levels)
  subtypes <- intersect(as.character(subtypes), inf_vcm_subtype_levels)
  if (length(segments) == 0) stop("Choose at least one valid influenza segment.")
  if (length(subtypes) == 0) stop("Choose at least one valid influenza subtype.")

  unknown_header_fields <- setdiff(header_fields, names(data))
  if (length(unknown_header_fields) > 0) {
    stop("Unknown FASTA header fields: ", paste(unknown_header_fields, collapse = ", "))
  }
  if (!identical(header_fields[[1]], "gisaid_isolate_name")) {
    stop("The first FASTA header field must be gisaid_isolate_name.")
  }

  selected_data <- prepare_inf_vcm_fasta_data(
    data = data,
    season = season,
    subtypes = subtypes
  )
  if (nrow(selected_data) == 0) {
    stop("No cleaned, GISAID-named H1, H3 or BVIC samples matched the requested export.")
  }

  output_root <- file.path(
    output_dir,
    "FASTA",
    inf_vcm_safe_filename(if (length(season) == 1 && nzchar(season)) season else "all_seasons")
  )
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  sample_metadata <- selected_data %>%
    select(
      any_of(c(
        "key",
        "export_subtype",
        "export_reason",
        "reference_in_prove_kategori",
        "gisaid_name_present",
        "season",
        "prove_kategori",
        "ngs_sekvens_resultat",
        "decoded_segment_count",
        "segments_decoded",
        header_fields
      ))
    ) %>%
    distinct(key, .keep_all = TRUE)

  fasta_records <- inf_sequence_segment_summary %>%
    filter(segment %in% segments, !is.na(sequence), nzchar(sequence)) %>%
    inner_join(sample_metadata, by = "key")

  header_values <- lapply(
    header_fields,
    function(field) fasta_header_value(fasta_records[[field]])
  )
  fasta_records$fasta_header <- do.call(paste, c(header_values, sep = "|"))

  duplicate_headers <- fasta_records %>%
    count(export_subtype, segment, fasta_header, name = "n") %>%
    filter(n > 1)
  if (nrow(duplicate_headers) > 0) {
    stop(
      "Duplicate FASTA headers were found within a subtype/segment. ",
      "Resolve duplicate GISAID isolate names before export."
    )
  }

  output_paths <- character()
  file_manifest <- list()
  for (subtype in subtypes) {
    subtype_dir <- file.path(output_root, subtype)
    dir.create(subtype_dir, recursive = TRUE, showWarnings = FALSE)

    for (segment_name in segments) {
      segment_records <- fasta_records %>%
        filter(export_subtype == .env$subtype, segment == .env$segment_name) %>%
        arrange(desc(reference_in_prove_kategori), gisaid_isolate_name, key)

      output_path <- file.path(
        subtype_dir,
        paste0(
          "Influenza_",
          inf_vcm_safe_filename(if (length(season) == 1 && nzchar(season)) season else "all_seasons"),
          "_",
          subtype,
          "_",
          segment_name,
          ".fasta"
        )
      )
      fasta_lines <- unlist(Map(
        function(header, sequence) c(paste0(">", header), sequence),
        segment_records$fasta_header,
        segment_records$sequence
      ))
      writeLines(fasta_lines, con = output_path, useBytes = TRUE)

      output_paths <- c(output_paths, output_path)
      file_manifest[[length(file_manifest) + 1]] <- tibble(
        subtype = subtype,
        segment = segment_name,
        fasta_file = output_path,
        sequence_count = nrow(segment_records),
        reference_count = sum(segment_records$reference_in_prove_kategori)
      )
    }
  }

  season_values <- as.character(season)
  season_values <- season_values[!is.na(season_values) & nzchar(season_values)]
  reference_audit <- data %>%
    mutate(
      export_subtype = inf_vcm_subtype(ngs_sekvens_resultat, prove_kategori),
      reference_in_prove_kategori = inf_vcm_reference_flag(prove_kategori),
      gisaid_name_present = !is_blank(gisaid_isolate_name),
      selected_season = if (length(season_values) == 0) TRUE else .data$season %in% .env$season_values,
      included_in_export = selected_season &
        reference_in_prove_kategori &
        gisaid_name_present &
        export_subtype %in% .env$subtypes
    ) %>%
    filter(
      reference_in_prove_kategori,
      selected_season,
      !is.na(export_subtype),
      export_subtype %in% .env$subtypes
    ) %>%
    transmute(
      key,
      season,
      export_subtype,
      prove_kategori,
      gisaid_isolate_name,
      decoded_segment_count,
      segments_decoded,
      gisaid_name_present,
      included_in_export
    ) %>%
    arrange(export_subtype, gisaid_isolate_name)

  sample_manifest_path <- file.path(output_root, "Influenza_Fasta_sample_manifest.csv")
  segment_manifest_path <- file.path(output_root, "Influenza_Fasta_segment_manifest.csv")
  reference_audit_path <- file.path(output_root, "Influenza_Fasta_reference_audit.csv")
  write.csv(sample_metadata, sample_manifest_path, row.names = FALSE, na = "")
  write.csv(bind_rows(file_manifest), segment_manifest_path, row.names = FALSE, na = "")
  write.csv(reference_audit, reference_audit_path, row.names = FALSE, na = "")

  invisible(list(
    fasta_files = unname(output_paths),
    sample_manifest = sample_manifest_path,
    segment_manifest = segment_manifest_path,
    reference_audit = reference_audit_path
  ))
}


# =============================================================================
# Subtype-first VCM raw-sequence export
#
# GISAID isolate name is the header for surveillance samples. References use
# the BioNumerics key because their GISAID isolate name is intentionally blank.
# =============================================================================

inf_vcm_segment_order <- c("HA", "NA", "PB2", "PB1", "PA", "NP", "M", "NS")

inf_vcm_raw_filename <- function(subtype, segment) {
  paste0("raw_sequences_", subtype, "_", segment, ".fasta")
}

inf_vcm_aligned_filename <- function(subtype, segment) {
  paste0("aligned_sequences_", subtype, "_", segment, ".fasta")
}

inf_vcm_phylogeny_prefix <- function(subtype, segment) {
  paste0(subtype, "_", segment, "_phylogeny")
}

inf_vcm_stable_strain <- function(key) {
  strain <- trimws(as.character(key))
  if (any(is.na(strain) | !nzchar(strain))) {
    stop("A stable FASTA sample identifier is missing.")
  }
  fasta_header_value(strain)
}

write_inf_vcm_fasta <- function(
  data = inf_sequence_export_data,
  output_dir = Sys.getenv(
    "INF_VCM_FASTA_OUTPUT_DIR",
    unset = "N:/Virologi/Influensa/2526/WGS_Analyse/VCM 2026_27"
  ),
  season = Sys.getenv("INF_VCM_FASTA_SEASON", unset = "Season25_26"),
  subtypes = inf_vcm_subtype_levels,
  segments = inf_vcm_segment_order
) {
  subtypes <- intersect(as.character(subtypes), inf_vcm_subtype_levels)
  segments <- inf_vcm_segment_order[
    inf_vcm_segment_order %in% intersect(as.character(segments), segment_levels)
  ]
  if (length(subtypes) == 0) stop("Choose at least one valid influenza subtype.")
  if (length(segments) == 0) stop("Choose at least one valid influenza segment.")

  metadata_source_fields <- c(
    "pasient_landsdel",
    "pasient_fylke_name",
    "pasient_fylke_nr",
    "pasient_alder",
    "pasient_aldersgruppe",
    "pasient_kjonn",
    "pasient_vaks",
    "pasient_vaks_2uipt",
    "pasient_antiviralbehandling"
  )
  for (field in setdiff(metadata_source_fields, names(data))) {
    data[[field]] <- NA_character_
  }

  selected_data <- prepare_inf_vcm_fasta_data(
    data = data,
    season = season,
    subtypes = subtypes
  )
  if (nrow(selected_data) == 0) {
    stop("No cleaned H1, H3 or BVIC samples matched the requested export.")
  }

  sample_metadata <- selected_data %>%
    mutate(
      strain = if_else(
        reference_in_prove_kategori,
        inf_vcm_stable_strain(key),
        fasta_header_value(gisaid_isolate_name)
      ),
      header_source = if_else(
        reference_in_prove_kategori,
        "key",
        "gisaid_isolate_name"
      ),
      subtype = export_subtype,
      collection_date = if_else(
        is.na(prove_tatt),
        NA_character_,
        format(as.Date(prove_tatt), "%Y-%m-%d")
      ),
      collection_year_month = if_else(
        is.na(prove_tatt),
        NA_character_,
        format(as.Date(prove_tatt), "%Y-%m")
      )
    ) %>%
    transmute(
      key,
      strain,
      header_source,
      subtype,
      collection_date,
      collection_year_month,
      gisaid_isolate_name,
      clade = nc_ha_clade,
      subclade = nc_ha_subclade,
      tessy_category = tessy_reportable_variable,
      sample_category = prove_kategori,
      region = pasient_landsdel,
      county = pasient_fylke_name,
      county_number = pasient_fylke_nr,
      age = pasient_alder,
      age_group = pasient_aldersgruppe,
      sex = pasient_kjonn,
      vaccination_status = pasient_vaks,
      vaccination_status_2uipt = pasient_vaks_2uipt,
      antiviral_therapy = pasient_antiviralbehandling,
      export_reason,
      bn_key = key
    ) %>%
    distinct(key, .keep_all = TRUE)

  fasta_records <- inf_sequence_segment_summary %>%
    filter(segment %in% segments, !is.na(sequence), nzchar(sequence)) %>%
    inner_join(sample_metadata, by = "key")

  duplicate_strains <- fasta_records %>%
    count(subtype, segment, strain, name = "n") %>%
    filter(n > 1)
  if (nrow(duplicate_strains) > 0) {
    stop("Duplicate stable sample identifiers were found within a subtype/segment.")
  }

  output_root <- file.path(output_dir, "export")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  output_paths <- character()
  metadata_paths <- character()
  for (subtype_name in subtypes) {
    subtype_dir <- file.path(output_root, subtype_name)
    dir.create(subtype_dir, recursive = TRUE, showWarnings = FALSE)

    generated_fasta_paths <- file.path(
      subtype_dir,
      vapply(
        inf_vcm_segment_order,
        function(segment_name) inf_vcm_raw_filename(subtype_name, segment_name),
        character(1)
      )
    )
    unlink(generated_fasta_paths[file.exists(generated_fasta_paths)], force = TRUE)
    unlink(file.path(subtype_dir, "metadata.tsv"), force = TRUE)

    subtype_records <- fasta_records %>%
      filter(subtype == .env$subtype_name)

    for (segment_name in segments) {
      segment_records <- subtype_records %>%
        filter(segment == .env$segment_name) %>%
        arrange(desc(export_reason == "reference"), strain)

      if (nrow(segment_records) == 0) next

      output_path <- file.path(
        subtype_dir,
        inf_vcm_raw_filename(subtype_name, segment_name)
      )
      fasta_lines <- unlist(Map(
        function(strain, sequence) c(paste0(">", strain), sequence),
        segment_records$strain,
        segment_records$sequence
      ))
      writeLines(fasta_lines, con = output_path, useBytes = TRUE)
      output_paths <- c(output_paths, output_path)
    }

    # Metadata describes samples/strains, not segment records. Restrict to
    # strains which produced at least one FASTA record, then retain the single
    # sample-level row prepared above.
    subtype_metadata <- sample_metadata %>%
      filter(
        subtype == .env$subtype_name,
        key %in% unique(subtype_records$key)
      ) %>%
      transmute(
        strain,
        subtype,
        header_source,
        collection_date,
        collection_year_month,
        gisaid_isolate_name,
        clade,
        subclade,
        tessy_category,
        sample_category,
        region,
        county,
        county_number,
        age,
        age_group,
        sex,
        vaccination_status,
        vaccination_status_2uipt,
        antiviral_therapy,
        export_reason,
        bn_key
      ) %>%
      arrange(strain)

    if (nrow(subtype_metadata) == 0) next

    metadata_path <- file.path(subtype_dir, "metadata.tsv")
    write.table(
      subtype_metadata,
      file = metadata_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      na = ""
    )
    metadata_paths <- c(metadata_paths, metadata_path)
  }

  invisible(list(
    fasta_files = unname(output_paths),
    metadata_files = unname(metadata_paths)
  ))
}
