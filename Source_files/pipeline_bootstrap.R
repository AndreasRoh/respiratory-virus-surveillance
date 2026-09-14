# Shared pipeline bootstrap helpers used by pathogen-specific runners.

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
log_timed_message <- function(...) {
  message(sprintf("[%s]", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), " ", paste0(..., collapse = ""))
  flush.console()
}

timed_step <- function(step_name, expr) {
  step_started_at <- Sys.time()
  log_timed_message("START: ", step_name)
  result <- withCallingHandlers(
    force(expr),
    warning = function(w) {
      if (grepl("was built under R version", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  step_elapsed <- as.numeric(difftime(Sys.time(), step_started_at, units = "secs"))
  log_timed_message("DONE: ", step_name, " (", sprintf("%.2f", step_elapsed), "s)")
  result
}


load_required_libraries <- function(packages, suppress_startup = FALSE) {
  lapply(packages, function(pkg) {
    load_expr <- substitute(library(pkg_name, character.only = TRUE), list(pkg_name = pkg))
    withCallingHandlers(
      if (isTRUE(suppress_startup)) suppressPackageStartupMessages(eval(load_expr)) else eval(load_expr),
      warning = function(w) {
        if (grepl("was built under R version", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  })
}
check_install_update_packages <- function(packages, update_option_name = "sc2_update_packages") {
  installed_pkgs <- rownames(installed.packages())
  missing_packages <- setdiff(packages, installed_pkgs)
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }

  if (isTRUE(getOption(update_option_name, FALSE))) {
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

init_locale <- function(locale_candidates = c("nb_NO.UTF-8", "Norwegian (Bokmal)_Norway.utf8", "Norwegian")) {
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

reset_output_sinks <- function() {
  while (sink.number() > 0) sink()
  while (sink.number(type = "message") != 2) sink(type = "message")
  invisible(NULL)
}

init_report_pipeline <- function(
  required_packages,
  common_utils_path = "Source_files/common_report_utils.R",
  suppress_startup = FALSE,
  set_datatable_legacy = FALSE,
  global_vars = NULL,
  init_locale_first = FALSE,
  reset_sinks_first = FALSE
) {
  if (isTRUE(init_locale_first)) {
    init_locale()
  }
  if (isTRUE(reset_sinks_first)) {
    reset_output_sinks()
  }

  invisible(timed_step("Load libraries", load_required_libraries(required_packages, suppress_startup = suppress_startup)))

  if (!is.null(global_vars) && length(global_vars) > 0) {
    invisible(utils::globalVariables(global_vars))
  }
  if (isTRUE(set_datatable_legacy)) {
    options(datatable.week = "legacy")
  }

  invisible(timed_step("Source common report utilities", source(common_utils_path)))
  invisible(NULL)
}

