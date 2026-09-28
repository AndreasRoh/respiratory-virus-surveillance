# =============================================================================
# Influenza VCM FASTA export runner
#
# Writes subtype- and segment-specific FASTA files using the defaults in
# INF_SequenceExportData.R. Override output or season with environment
# variables INF_VCM_FASTA_OUTPUT_DIR and INF_VCM_FASTA_SEASON.
# =============================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE))
} else {
  project_inf_dir <- file.path(getwd(), "INF")
  if (file.exists(file.path(project_inf_dir, "INF_SequenceExportData.R"))) {
    project_inf_dir
  } else getwd()
}

source(file.path(script_dir, "INF_SequenceExportData.R"))
write_inf_vcm_fasta()
