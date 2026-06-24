audit_results_dir <- file.path(getwd(), "Results_Audit")
dir.create(audit_results_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(INF_RESULTS_DIR = audit_results_dir)
Sys.setenv(INF_RESULTS_SHARE_DIR = audit_results_dir)

# Audit runner for influenza analysis output.
source(file.path(getwd(), "INF", "INF_Analysis.R"))
source(file.path(getwd(), "INF", "INF_PPT_Audit.R"))