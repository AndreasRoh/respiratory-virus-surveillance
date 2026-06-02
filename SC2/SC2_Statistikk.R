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

required_packages <- c(
  "dplyr", "lubridate", "tidyr", "data.table",
  "stringr", "janitor", "dbplyr", "tsibble", "odbc", "DBI"
)
check_install_update_packages(required_packages)
suppressPackageStartupMessages(lapply(required_packages, library, character.only = TRUE))
like <- data.table::like
yearmonth <- function(x) as.Date(format(as.Date(x), "%Y-%m-01"))

source(file.path(bundle_scripts_dir, "SC2_SQLquery_BNCOVID19.R"))
source(file.path(bundle_scripts_dir, "SC2_SQLquery_25-26.R"))
source(file.path(bundle_scripts_dir, "SC2_DataCleaning_BNCOVID19.R"))
source(file.path(bundle_scripts_dir, "SC2_DataCleaning_25-26.R"))

if (!exists("SC2db")) {
  if (exists("SC2_20_25_raw_merged") && exists("SC2_25_26_raw_merged")) {
    SC2db <- dplyr::bind_rows(SC2_20_25_raw_merged, SC2_25_26_raw_merged)
  } else if (exists("SC2_25_26_raw_merged")) {
    SC2db <- SC2_25_26_raw_merged
  }
}
if (!exists("SC2db")) {
  stop("Missing input data: expected `SC2db` or `SC2_*_raw_merged` after sourcing SC2 SQL query scripts")
}
if (!("my" %in% names(SC2db)) && ("prove_tatt" %in% names(SC2db))) {
  SC2db$my <- as.Date(SC2db$prove_tatt)
}
if ("my" %in% names(SC2db)) {
  SC2db$my <- as.Date(paste0(format(as.Date(SC2db$my), "%Y-%m"), "-01"))
}

source(file.path(bundle_scripts_dir, "SC2_Classification.R"))
if (!exists("allvariants_v") && exists("SC2db")) {
  allvariants_v <- SC2db
}
if (!exists("allvariants_v")) {
  stop("Missing input data: expected `allvariants_v` (or SC2db fallback) after sourcing SC2_Classification.R")
}

Sys.setlocale("LC_TIME", "nb_NO.UTF-8")

OUTPUT_DIR <- "N:/Virologi/Influensa/2526/WGS_Analyse/Results/Statistikk"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

format_month_label <- function(x) {
  format(as.Date(x), "%Y %b") |> tolower()
}

month_seq <- function(from_date, to_date) {
  seq(from = floor_date(as.Date(from_date), "month"), to = floor_date(as.Date(to_date), "month"), by = "month") |>
    format_month_label()
}

build_monthly_variant_counts <- function(data, month_col, variant_col, months, variant_name) {
  month_sym <- rlang::sym(month_col)
  variant_sym <- rlang::sym(variant_col)

  variants <- data |>
    filter(!is.na(!!variant_sym), !!variant_sym != "") |>
    distinct(!!variant_sym) |>
    pull(!!variant_sym) |>
    sort()

  if (length(variants) == 0) {
    return(tibble(my = character(), variant = character(), Antall = integer(), flagg = integer()) |>
      rename(!!variant_name := variant))
  }

  data |>
    filter(!is.na(!!variant_sym), !!variant_sym != "") |>
    mutate(my = format_month_label(!!month_sym)) |>
    count(my, !!variant_sym, name = "Antall") |>
    rename(variant = !!variant_sym) |>
    complete(my = months, variant = variants, fill = list(Antall = 0)) |>
    mutate(flagg = 0L) |>
    arrange(my, variant) |>
    rename(!!variant_name := variant)
}

write_stat_csv <- function(data, suffix, year_value, week_value, output_dir) {
  file_name <- paste0("SARSCOV2_", year_value, "_Week", week_value, "_", suffix, ".csv")
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

sc2_base <- allvariants_v |>
  mutate(
    prove_tatt_date = as.Date(prove_tatt),
    my_date = as.Date(my)
  )

max_date <- max(sc2_base$prove_tatt_date, na.rm = TRUE)
start_6m <- max_date %m-% months(6)
months_6m <- month_seq(start_6m, max_date)

stat_monthly <- build_monthly_variant_counts(
  data = sc2_base |> filter(prove_tatt_date > start_6m),
  month_col = "prove_tatt_date",
  variant_col = "nc_pangolin_short",
  months = months_6m,
  variant_name = "nc_pangolin_short"
)

vum_monthly <- build_monthly_variant_counts(
  data = sc2_base |> filter(prove_tatt_date > start_6m),
  month_col = "prove_tatt_date",
  variant_col = "VUM",
  months = months_6m,
  variant_name = "VUM"
)

voi_monthly <- build_monthly_variant_counts(
  data = sc2_base |> filter(prove_tatt_date > start_6m),
  month_col = "prove_tatt_date",
  variant_col = "VOI",
  months = months_6m,
  variant_name = "VOI"
)

latest_my <- max(sc2_base$my_date, na.rm = TRUE)
start_2m <- floor_date(latest_my %m-% months(1), "month")
months_2m <- month_seq(start_2m, latest_my)

top10_variants <- sc2_base |>
  filter(my_date >= start_2m, my_date <= latest_my, !is.na(nc_pangolin_short), nc_pangolin_short != "") |>
  count(nc_pangolin_short, sort = TRUE) |>
  slice_head(n = 10) |>
  pull(nc_pangolin_short)

top10_monthly <- sc2_base |>
  filter(my_date >= start_2m, my_date <= latest_my, nc_pangolin_short %in% top10_variants) |>
  mutate(my = format_month_label(my_date)) |>
  count(my, nc_pangolin_short, name = "Antall") |>
  complete(my = months_2m, nc_pangolin_short = top10_variants, fill = list(Antall = 0)) |>
  mutate(flagg = 0L) |>
  arrange(my, nc_pangolin_short)

all_months <- month_seq(min(sc2_base$my_date, na.rm = TRUE), max(sc2_base$my_date, na.rm = TRUE))

stat_total <- sc2_base |>
  filter(!is.na(Collapsed_pango), Collapsed_pango != "") |>
  mutate(my = format_month_label(my_date)) |>
  count(my, Collapsed_pango, name = "Antall") |>
  complete(my = all_months, Collapsed_pango, fill = list(Antall = 0)) |>
  group_by(my) |>
  mutate(
    total_antall = sum(Antall),
    Prosent = round(if_else(total_antall > 0, (Antall / total_antall) * 100, 0), 0),
    flagg = 0L
  ) |>
  ungroup() |>
  select(-total_antall) |>
  arrange(my, Collapsed_pango)

run_year <- year(Sys.Date())
run_week <- week(Sys.Date())

write_stat_csv(stat_monthly, "Stat", run_year, run_week, OUTPUT_DIR)
write_stat_csv(vum_monthly, "VUM", run_year, run_week, OUTPUT_DIR)
write_stat_csv(voi_monthly, "VOI", run_year, run_week, OUTPUT_DIR)
write_stat_csv(top10_monthly, "top10", run_year, run_week, OUTPUT_DIR)
write_stat_csv(stat_total, "statistikk", run_year, run_week, OUTPUT_DIR)
