# =============================================================================
# WORK IN PROGRESS — reviewed, behaviour-preserving copy
#
# Purpose: This file is a maintainability-focused review copy of INF_Analysis_Audit.R.
# Changes in the WIP series are limited to structure, documentation, and WIP
# dependency isolation. The calculations, filters, object names, and exported
# outputs are retained so results can be compared directly with production.
# =============================================================================

audit_results_dir <- file.path(getwd(), "Results_Audit")
dir.create(audit_results_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(INF_RESULTS_DIR = audit_results_dir)
Sys.setenv(INF_RESULTS_SHARE_DIR = audit_results_dir)

# Audit runner for influenza analysis output.
source(file.path(getwd(), "INF", "INF_Analysis.R"))
source(file.path(getwd(), "INF", "INF_PPT_Audit.R"))