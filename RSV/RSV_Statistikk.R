resolve_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

bundle_scripts_dir <- resolve_script_dir()

options(repos = c(CRAN = "https://cran.uni-muenster.de/"))

check_install_update_packages <- function(packages) {
  installed_pkgs <- rownames(installed.packages())
  missing_packages <- setdiff(packages, installed_pkgs)
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }

  outdated_packages <- old.packages()
  if (!is.null(outdated_packages)) {
    outdated_required <- intersect(rownames(outdated_packages), packages)
    if (length(outdated_required) > 0) {
      message("Updating outdated required packages: ", paste(outdated_required, collapse = ", "))
      tryCatch(
        update.packages(oldPkgs = outdated_required, ask = FALSE, checkBuilt = TRUE),
        error = function(e) message("Package update skipped: ", conditionMessage(e))
      )
    }
  }
}

required_packages <- c("dplyr", "lubridate")
check_install_update_packages(required_packages)
suppressPackageStartupMessages(lapply(required_packages, library, character.only = TRUE))

source(file.path(bundle_scripts_dir, "RSV_SQLquery.R"))
source(file.path(bundle_scripts_dir, "RSV_DataCleaning_23-24.R"))

if (!exists("rsvdb")) {
  stop("Missing input data: expected `rsvdb` after sourcing RSV SQL query and data cleaning scripts.")
}

OUTPUT_DIR <- "N:/Virologi/Influensa/2526/WGS_Analyse/Results/Statistikk"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

season_label <- function(x) {
  d <- as.Date(x)
  y <- year(d)
  m <- month(d)
  start_year <- ifelse(m >= 9, y, y - 1)
  end_year_short <- sprintf("%02d", (start_year + 1) %% 100)
  paste0(start_year, "-", end_year_short)
}

write_stat_csv <- function(data, subtype_label, year_value, week_value, output_dir) {
  file_name <- paste0(subtype_label, "_", year_value, "_Week", week_value, "_statistikk.csv")
  file_path <- file.path(output_dir, file_name)
  write.table(
    data,
    file = file_path,
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    qmethod = "double"
  )
  message("Saved: ", file_path)
}

build_subtype_table <- function(data, subtype_regex, subtype_col_name, clade_col, current_season) {
  filtered <- data %>%
    filter(
      grepl(subtype_regex, toupper(as.character(ngs_sekvens_resultat))),
      !is.na(.data[[clade_col]]),
      trimws(as.character(.data[[clade_col]])) != "",
      !is.na(prove_tatt)
    ) %>%
    mutate(
      Sesong = season_label(prove_tatt),
      clade_value = trimws(as.character(.data[[clade_col]]))
    ) %>%
    filter(Sesong == current_season) %>%
    count(clade_value, Sesong, name = "Antall") %>%
    mutate(
      total_antall = sum(Antall),
      Prosent = round(if_else(total_antall > 0, (Antall / total_antall) * 100, 0), 2),
      FLAGG = 0L
    ) %>%
    select(-total_antall) %>%
    arrange(clade_value) %>%
    rename(!!subtype_col_name := clade_value)

  filtered
}

if (!("prove_tatt" %in% names(rsvdb))) {
  stop("Missing required column: `prove_tatt` in rsvdb.")
}
if (!("ngs_sekvens_resultat" %in% names(rsvdb))) {
  stop("Missing required column: `ngs_sekvens_resultat` in rsvdb.")
}

rsvdb <- rsvdb %>%
  mutate(prove_tatt = as.Date(prove_tatt))

clade_col <- if ("ngs_subclade" %in% names(rsvdb)) "ngs_subclade" else "ngs_clade"
if (!(clade_col %in% names(rsvdb))) {
  stop("Missing required clade column: expected `ngs_subclade` or `ngs_clade` in rsvdb.")
}

current_season <- season_label(Sys.Date())
rsv_current_season <- rsvdb %>%
  filter(season_label(prove_tatt) == current_season)

rsv_a <- build_subtype_table(
  data = rsv_current_season,
  subtype_regex = "A",
  subtype_col_name = "RSV A",
  clade_col = clade_col,
  current_season = current_season
)

rsv_b <- build_subtype_table(
  data = rsv_current_season,
  subtype_regex = "B",
  subtype_col_name = "RSV B",
  clade_col = clade_col,
  current_season = current_season
)

run_year <- format(Sys.Date(), "%Y")
run_week <- format(Sys.Date(), "%U")

write_stat_csv(rsv_a, "RSV_A", run_year, run_week, OUTPUT_DIR)
write_stat_csv(rsv_b, "RSV_B", run_year, run_week, OUTPUT_DIR)
