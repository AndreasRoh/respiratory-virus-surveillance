# FLuDB QC for 25-26
# Goal: flag consistency issues and quality outliers in fludb.

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

suppressPackageStartupMessages(suppressWarnings({
  library(data.table)
  library(dplyr)
  library(janitor)
  library(lubridate)
  library(stringr)
}))

if (!exists("fludb")) {
  source(file.path(bundle_scripts_dir, "INF_SQLquery_25-26.R"))
  source(file.path(bundle_scripts_dir, "INF_DataCleaning_25-26.R"))
}

if (!exists("fludb")) {
  stop("Missing input data: expected `fludb` after sourcing the INF data scripts.")
}

setDT(fludb)

output_dir <- file.path(bundle_scripts_dir, "..", "Results", "FLuDB_QC")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

row_any_nonblank <- function(dt, cols) {
  if (length(cols) == 0) {
    return(rep(FALSE, nrow(dt)))
  }
  data_list <- lapply(dt[, ..cols], normalize_chr)
  Reduce(`|`, lapply(data_list, function(x) x != ""))
}

row_x_count <- function(dt, cols) {
  if (length(cols) == 0) {
    return(rep.int(0L, nrow(dt)))
  }
  data_list <- lapply(dt[, ..cols], normalize_chr)
  counts <- lapply(data_list, function(x) stringr::str_count(x, "X"))
  as.integer(rowSums(as.data.frame(counts), na.rm = TRUE))
}

flag_map <- data.table(
  issue_code = c(
    "no_ngs_but_signal",
    "no_ngs_but_gisaid",
    "no_ngs_but_nextclade",
    "no_ngs_but_mutations",
    "x_in_mutations",
    "high_x_load",
    "adamantane_irrelevant_for_b",
    "missing_gisaid_isolate_name"
  ),
  issue_label = c(
    "NGS-Resultat is blank, but other sequence-linked fields are populated",
    "NGS-Resultat is blank, but at least one GISAID field is populated",
    "NGS-Resultat is blank, but at least one nc_ field is populated",
    "NGS-Resultat is blank, but mutation / resistance calls are populated",
    "Mutation columns contain X characters",
    "Mutation X-load is in the top 5 percent of non-zero rows",
    "Adamantane-related fields are populated for a B lineage sample",
    "GISAID isolate name is missing while sequence-linked fields are present"
  )
)

gisaid_cols <- intersect(
  c(
    "gisaid_isolate_name", "gisaid_isolate_id",
    "gisaid_ha_id", "gisaid_na_id", "gisaid_m_id", "gisaid_pa_id",
    "gisaid_pb1_id", "gisaid_pb2_id"
  ),
  names(fludb)
)

nc_cols <- grep("^nc_", names(fludb), value = TRUE)

mutation_cols <- intersect(
  c(
    "mut_ha1_1", "mut_ha1_vac", "mut_ha1_cluster_defining",
    "mut_ha1_without_subclade_signature", "dr_na_mut", "dr_pa_mut", "dr_m2_mut"
  ),
  names(fludb)
)

drug_cols <- intersect(
  c(
    "dr_res_adamantine", "dr_res_baloxavir", "dr_res_oseltamivir",
    "dr_res_peramivir", "dr_res_zanamivir", "dr_res_laninamivir"
  ),
  names(fludb)
)

sequence_fields <- unique(c(gisaid_cols, nc_cols, mutation_cols, drug_cols))

if (!"ngs_sekvens_resultat" %in% names(fludb)) {
  fludb[, ngs_sekvens_resultat := NA_character_]
}

fludb[, row_id := .I]
fludb[, ngs_blank := normalize_chr(ngs_sekvens_resultat) == ""]
fludb[, has_gisaid := row_any_nonblank(.SD, gisaid_cols), .SDcols = gisaid_cols]
fludb[, has_nc := row_any_nonblank(.SD, nc_cols), .SDcols = nc_cols]
fludb[, has_mutation := row_any_nonblank(.SD, mutation_cols), .SDcols = mutation_cols]
fludb[, has_drug_call := row_any_nonblank(.SD, drug_cols), .SDcols = drug_cols]
fludb[, x_count := row_x_count(.SD, mutation_cols), .SDcols = mutation_cols]
fludb[, x_quality_index := fifelse(x_count == 0L, 1, 1 / (1 + x_count))]
fludb[, b_lineage := !is.na(ngs_sekvens_resultat) & str_detect(ngs_sekvens_resultat, "^B/")]
fludb[, has_signal := row_any_nonblank(.SD, sequence_fields), .SDcols = sequence_fields]

if ("gisaid_isolate_name" %in% names(fludb)) {
  fludb[, gisaid_isolate_name := normalize_chr(gisaid_isolate_name)]
} else {
  fludb[, gisaid_isolate_name := ""]
}

if ("dr_res_adamantine" %in% names(fludb)) {
  fludb[, adamantine_nonblank := normalize_chr(dr_res_adamantine) != ""]
} else {
  fludb[, adamantine_nonblank := FALSE]
}

if ("dr_m2_mut" %in% names(fludb)) {
  fludb[, m2_mut_nonblank := normalize_chr(dr_m2_mut) != ""]
} else {
  fludb[, m2_mut_nonblank := FALSE]
}

fludb[, `:=`(
  flag_no_ngs_but_signal = ngs_blank & has_signal,
  flag_no_ngs_but_gisaid = ngs_blank & has_gisaid,
  flag_no_ngs_but_nextclade = ngs_blank & has_nc,
  flag_no_ngs_but_mutations = ngs_blank & (has_mutation | has_drug_call),
  flag_x_in_mutations = x_count > 0L,
  flag_adamantane_irrelevant_for_b = b_lineage & (adamantine_nonblank | m2_mut_nonblank),
  flag_missing_gisaid_isolate_name = has_signal & gisaid_isolate_name == ""
)]

non_zero_x <- fludb[x_count > 0L, x_count]
x_outlier_threshold <- if (length(non_zero_x) == 0L) Inf else as.integer(quantile(non_zero_x, 0.95, names = FALSE, type = 7))
fludb[, flag_high_x_load := x_count > 0L & x_count >= x_outlier_threshold]

flag_cols <- paste0("flag_", flag_map$issue_code)
flag_matrix <- as.matrix(fludb[, ..flag_cols])
qc_flagged_rows <- fludb[rowSums(flag_matrix, na.rm = TRUE) > 0L,
  c(
    "row_id",
    intersect(c("key", "gisaid_isolate_name", "prove_tatt", "season", "month_year", "ngs_sekvens_resultat"), names(fludb)),
    "x_count",
    "x_quality_index",
    flag_cols
  ),
  with = FALSE
]

if (nrow(qc_flagged_rows) > 0) {
  issue_long <- melt(
    qc_flagged_rows,
    id.vars = setdiff(names(qc_flagged_rows), flag_cols),
    measure.vars = flag_cols,
    variable.name = "issue_code",
    value.name = "flag"
  )[flag == TRUE]

  qc_summary <- merge(
    flag_map,
    issue_long[, .(flagged_n = .N), by = .(issue_code = sub("^flag_", "", issue_code))],
    by = "issue_code",
    all.x = TRUE
  )
  qc_summary[is.na(flagged_n), flagged_n := 0L]
  setorder(qc_summary, -flagged_n, issue_code)
} else {
  qc_summary <- flag_map[, .(issue_code, issue_label)][, flagged_n := 0L][]
  issue_long <- qc_flagged_rows[0]
}

fwrite(
  qc_flagged_rows,
  file.path(output_dir, "fludb_qc_flagged_rows_2526.csv2"),
  sep = ";",
  dec = ",",
  na = ""
)

fwrite(
  qc_summary,
  file.path(output_dir, "fludb_qc_summary_2526.csv2"),
  sep = ";",
  dec = ",",
  na = ""
)

qc_outliers <- fludb[flag_high_x_load == TRUE,
  .(
    row_id,
    key = if ("key" %in% names(fludb)) key else NA_character_,
    gisaid_isolate_name,
    ngs_sekvens_resultat,
    x_count,
    x_quality_index,
    season = if ("season" %in% names(fludb)) season else NA_character_
  )
]

fwrite(
  qc_outliers,
  file.path(output_dir, "fludb_qc_outliers_2526.csv2"),
  sep = ";",
  dec = ",",
  na = ""
)

message("FLuDB QC complete. Flagged rows: ", nrow(qc_flagged_rows), ", outliers: ", nrow(qc_outliers))

invisible(list(
  flagged_rows = qc_flagged_rows,
  summary = qc_summary,
  outliers = qc_outliers
))
