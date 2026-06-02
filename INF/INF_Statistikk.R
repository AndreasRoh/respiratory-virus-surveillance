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

required_packages <- c("dplyr", "lubridate", "tidyr")
check_install_update_packages(required_packages)
suppressPackageStartupMessages(lapply(required_packages, library, character.only = TRUE))

source("INF/INF_SQLquery_25-26.R")
source("INF/INF_DataCleaning_25-26.R")

if (!exists("fludb") && exists("INF_25_26_raw_merged")) {
  fludb <- INF_25_26_raw_merged
}
if (!exists("fludb")) {
  stop("Missing input data: expected `fludb` or `INF_25_26_raw_merged` after sourcing INF_SQLquery_25-26.R")
}

Sys.setlocale("LC_TIME", "nb_NO.utf8")

OUTPUT_DIR <- "N:/Virologi/Influensa/2526/WGS_Analyse/Results/Statistikk"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

format_month_label <- function(x) {
  format(as.Date(x), "%Y %b") |> tolower()
}

normalize_tessy <- function(x) {
  case_when(
    grepl("^genAH1", x) ~ sub("^genAH1/", "A/", x),
    grepl("^genAH3", x) ~ sub("^genAH3/", "A/", x),
    grepl("^genB", x) ~ sub("^genB/", "B/", x),
    TRUE ~ x
  )
}

build_reporting_category <- function(ngs_result, ha_clade, ha_subclade, tessy) {
  virus_type <- case_when(
    grepl("A/H1N1", ngs_result) ~ "A(H1)pdm09",
    grepl("A/H3N2", ngs_result) ~ "A(H3)",
    grepl("B/Victoria", ngs_result) ~ "B(Vic)",
    grepl("B/Yamagata", ngs_result) ~ "B(Yam)",
    TRUE ~ NA_character_
  )

  clade_part <- case_when(
    grepl("6B.1A.5a.2a.1", ha_clade) ~ "5a.2a.1",
    grepl("3C.2a1b.2a.2a.3a.1", ha_clade) ~ "2a.3a.1",
    TRUE ~ ha_clade
  )

  subclade_part <- if_else(!is.na(ha_subclade) & ha_subclade != "", paste0("(", ha_subclade, ")"), "")
  paste(virus_type, "klade", clade_part, subclade_part, "representativ", tessy)
}

build_subtype_table <- function(data, subtype) {
  subtype_data <- data |>
    filter(ngs_sekvens_resultat == subtype) |>
    mutate(my = format_month_label(prove_tatt))

  if (nrow(subtype_data) == 0) {
    return(tibble(my = character(), Rapporteringskategori = character(), Antall = integer(), Prosent = numeric(), flagg = integer()))
  }

  month_levels <- sort(unique(subtype_data$my))
  categories <- sort(unique(subtype_data$Rapporteringskategori))

  subtype_data |>
    count(my, Rapporteringskategori, name = "Antall") |>
    complete(my = month_levels, Rapporteringskategori = categories, fill = list(Antall = 0)) |>
    group_by(my) |>
    mutate(
      total_antall = sum(Antall),
      Prosent = round(if_else(total_antall > 0, (Antall / total_antall) * 100, 0), 2),
      flagg = 0L
    ) |>
    ungroup() |>
    select(-total_antall) |>
    arrange(my, Rapporteringskategori)
}

write_stat_csv <- function(data, subtype, year_value, week_value, output_dir) {
  file_name <- paste0(subtype, "_", year_value, "_Week", week_value, "_statistikk.csv")
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

fludb_clean <- fludb |>
  mutate(
    prove_tatt = as.Date(prove_tatt),
    week_year = isoweek(prove_tatt),
    year = year(prove_tatt)
  ) |>
  filter(!(year == 2025 & week_year < 35))

flu_filtered <- fludb_clean |>
  filter(
    ngs_sekvens_resultat %in% c("A/H1N1", "A/H3N2", "B/Victoria"),
    !is.na(tessy_reportable_variable),
    tessy_reportable_variable != "",
    tessy_reportable_variable != "NA",
    tessy_reportable_variable != "NULL"
  ) |>
  transmute(
    prove_tatt,
    ngs_sekvens_resultat,
    nc_ha_clade,
    nc_ha_subclade,
    tessy_reportable_variable = normalize_tessy(tessy_reportable_variable),
    Rapporteringskategori = build_reporting_category(
      ngs_sekvens_resultat,
      nc_ha_clade,
      nc_ha_subclade,
      normalize_tessy(tessy_reportable_variable)
    )
  )

h1n1_table <- build_subtype_table(flu_filtered, "A/H1N1")
h3n2_table <- build_subtype_table(flu_filtered, "A/H3N2")
bvic_table <- build_subtype_table(flu_filtered, "B/Victoria")

run_year <- format(Sys.Date(), "%Y")
run_week <- format(Sys.Date(), "%U")

write_stat_csv(h1n1_table, "A_H1N1", run_year, run_week, OUTPUT_DIR)
write_stat_csv(h3n2_table, "A_H3N2", run_year, run_week, OUTPUT_DIR)
write_stat_csv(bvic_table, "B_Victoria", run_year, run_week, OUTPUT_DIR)
