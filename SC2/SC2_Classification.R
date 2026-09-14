

######################   Collapsing Pangolins  ##################################
# https://www.who.int/activities/tracking-SARS-CoV-2-variants go here to update the list 
# Add a Collapsed pango column to the dataset based on the long Pangolin lineage 

SC2db <- SC2db %>%
  mutate(Collapsed_pango = case_when(
    like(nc_pangolin_long, "B.1.1.529.2.86.1.1.11.1.1.1.3.8.1") ~ "LP.8.1",
    like(nc_pangolin_long, "XDV.1.5.1.1.8.1") ~ "NB.1.8.1",
    like(nc_pangolin_long, "NB.1.8.1") ~ "NB.1.8.1",
    like(nc_pangolin_long, "KP.2.3") ~ "Andre SARS CoV 2",
    like(nc_pangolin_long, "XEC.") ~ "XEC",
    like(nc_pangolin_long, "KS.1.1KP.3.3") ~ "XEC",
    like(nc_pangolin_long, "B.1.1.529.2.86.1.1.11.1.2") ~ "KP.2",
    like(nc_pangolin_long, "B.1.1.529.2.86.1.1.11.1.3") ~ "KP.3",
    like(nc_pangolin_long, "BA.2.86.1.1.11.1.2") ~ "KP.2",
    like(nc_pangolin_long, "BA.2.86.1.1.11.1.3") ~ "KP.3",
    like(nc_pangolin_long, "B.1.1.529.2.86") ~ "BA.2.86",
    like(nc_pangolin_long, "BA.2.86") ~ "BA.2.86",
    like(nc_pangolin_long, "XBB.1.9.2.5.1") ~ "EG.5.1",
    like(nc_pangolin_long, "B.1.1.529.2.75.3.4.1.1.1.1") ~ "CH.1.1",
    like(nc_pangolin_long, "XBB.") ~ "XBB",
    like(nc_pangolin_long, "BJ.1BM.1.1.1") ~ "XBB",
    like(nc_pangolin_long, "B.1.1.529.5.3.1.1.1.1.1") ~ "BQ.1",
    like(nc_pangolin_long, "B.1.1.529.2.75") ~ "BA.2.75",
    like(nc_pangolin_long, "B.1.1.529.5") ~ "BA.5",
    like(nc_pangolin_long, "B.1.1.529.4") ~ "BA.4",
    like(nc_pangolin_long, "B.1.1.529.3") ~ "BA.3",
    like(nc_pangolin_long, "BA.3.2") ~ "BA.3.2",
    like(nc_pangolin_long, "B.1.1.529.2") ~ "BA.2",
    like(nc_pangolin_long, "B.1.1.529.1") ~ "BA.1",
    like(nc_pangolin_long, "B.1.1.529") ~ "Andre Omicron",
    like(nc_pangolin_long, "B.1.617") ~ "Delta",
    grepl("^B\\.1\\.1\\.7$", nc_pangolin_long) ~ "Alfa",  # Exact match for Alpha variant
    like(nc_pangolin_long, "B.1.351") ~ "Beta",
    like(nc_pangolin_long, "B.1.1.28.1") ~ "Gamma",
    like(nc_pangolin_long, "B.1.427") | like(nc_pangolin_long, "B.1.429") ~ "Epsilon",
    like(nc_pangolin_long, "B.1.1.28.2") ~ "Zeta",
    like(nc_pangolin_long, "B.1.525") ~ "Eta",
    like(nc_pangolin_long, "B.1.1.28.3") ~ "Theta",
    like(nc_pangolin_long, "B.1.526") ~ "Iota",
    like(nc_pangolin_long, "B.1.617.1") ~ "Kappa",
    like(nc_pangolin_long, "B.1.1.1.37.1") ~ "Lambda",
    like(nc_pangolin_long, "B.1.621") ~ "Mu",
    like(nc_pangolin_long, "LF.7LP.8.1.2LF.7") ~ "XFG",
    like(nc_pangolin_long, "XFG") ~ "XFG",
    TRUE ~ "Andre SARS CoV 2"
  ))


# FHI palette-driven variant colors (automatic, no manual per-variant hardcoding).
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

common_report_utils_path <- file.path(bundle_scripts_dir, "..", "Source_files", "common_report_utils.R")
if (file.exists(common_report_utils_path)) {
  source(common_report_utils_path)
}

variant_levels <- SC2db %>%
  dplyr::pull(Collapsed_pango) %>%
  as.character() %>%
  unique() %>%
  sort(na.last = TRUE)
variant_levels <- variant_levels[!is.na(variant_levels) & trimws(variant_levels) != ""]

base_palette <- if (exists("palette_all", inherits = TRUE)) {
  unname(palette_all)
} else if (exists("kvalitativ_comb", inherits = TRUE)) {
  kvalitativ_comb
} else {
  fhi_discrete_palette(6L, kvalitativ_comb)
}

variant_color <- stats::setNames(
  grDevices::colorRampPalette(base_palette)(length(variant_levels)),
  variant_levels
)

# Backward compatibility for existing references in downstream scripts.
custom_colors <- variant_color



# Classification by Origin of Sequences
SC2db <- SC2db %>%
  mutate(Origin = case_when(
    startsWith(key, "SUS-") ~ "Stavanger",
    startsWith(key, "HUS_n") ~ "Haukeland",
    startsWith(key, "OUS") ~ "OUS",
    startsWith(key, "STO") ~ "St. Olav",
    grepl("hCoV-19/Norway/Ahus", key) ~ "Ahus",
    grepl("hCoV_19_Norway_SIHF", key) ~ "SIHF",
    TRUE ~ "FHI"
  ))

# FHI palette-driven Origin colors (automatic, no manual hardcoding).
origin_levels <- SC2db %>%
  dplyr::pull(Origin) %>%
  as.character() %>%
  unique() %>%
  sort(na.last = TRUE)
origin_levels <- origin_levels[!is.na(origin_levels) & trimws(origin_levels) != ""]

origin_base_palette <- if (exists("palette_all", inherits = TRUE)) {
  unname(palette_all)
} else if (exists("kvalitativ_comb", inherits = TRUE)) {
  kvalitativ_comb
} else {
  fhi_discrete_palette(6L, kvalitativ_comb)
}

origin_color <- stats::setNames(
  grDevices::colorRampPalette(origin_base_palette)(length(origin_levels)),
  origin_levels
)





############################# TEssy variant mapping
# Load the first CSV file
mapping_cache_dir <- Sys.getenv(
  "SC2_TESSY_MAPPING_CACHE_DIR",
  unset = file.path(Sys.getenv("LOCALAPPDATA", unset = tempdir()), "respiratory-virus-surveillance", "ecdc")
)
mapping_urls <- c(
  "https://www.ecdc.europa.eu/sites/default/files/documents/PathogenVariant_public_mappings.csv",
  "https://www.ecdc.europa.eu/sites/default/files/documents/PathogenVariant_public_mappings_VUM.csv"
)

mapping_tables <- lapply(mapping_urls, function(url) {
  cache_path <- file.path(mapping_cache_dir, basename(url))
  downloaded <- tryCatch(utils::read.csv(url, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(downloaded)) {
    dir.create(mapping_cache_dir, recursive = TRUE, showWarnings = FALSE)
    try(utils::write.csv(downloaded, cache_path, row.names = FALSE), silent = TRUE)
    return(downloaded)
  }
  if (file.exists(cache_path)) {
    cached <- tryCatch(utils::read.csv(cache_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(cached)) return(cached)
  }
  warning("Unable to retrieve or read cached ECDC variant mapping: ", url)
  NULL
})

if (any(vapply(mapping_tables, is.null, logical(1)))) {
  is_audit_run <- exists("report_audit_mode", mode = "function") && report_audit_mode()
  if (!is_audit_run) {
    stop("ECDC variant mappings are unavailable. Retry later or provide cached CSV files in ", mapping_cache_dir)
  }
  warning("ECDC mappings are unavailable; audit mode is using Pangolin lineages as TESSy labels.")
  tessy_raw <- trimws(as.character(SC2db$nc_pangolin_short))
  SC2db$Tessy <- ifelse(is.na(tessy_raw) | tessy_raw == "", "Andre SARS CoV 2", tessy_raw)
} else {
  combined_variant_mappings <- do.call(rbind, Filter(Negate(is.null), mapping_tables))
  combined_variant_mappings$included.sub.lineages <- strsplit(combined_variant_mappings$included.sub.lineages, "\\|")
  combined_variant_mappings$included.sub.lineages <- lapply(combined_variant_mappings$included.sub.lineages, as.character)
  SC2db$Tessy <- vapply(SC2db$nc_pangolin_short, function(lineage) {
    matched <- combined_variant_mappings$VirusVariant[vapply(combined_variant_mappings$included.sub.lineages, function(values) lineage %in% values, logical(1))]
    if (length(matched) > 0) matched[[1]] else NA_character_
  }, character(1))
  SC2db$Tessy <- ifelse(is.na(SC2db$Tessy), "Andre SARS CoV 2", SC2db$Tessy)
}
# Mutation based VOI/VUM/VOC classification (not in use currently)


SC2db <- SC2db %>%
  mutate(
    VUM = case_when(
      grepl("XFG", Tessy)  ~ "XFG",
      grepl("NB.1.8.1", Tessy)  ~ "NB.1.8.1",
      grepl("BA.3.2", Tessy)  ~ "BA.3.2",
      TRUE ~ "" 
    )
  )

SC2db <- SC2db %>%
  mutate(
    VOI = case_when(
      grepl("BA.2.86", Tessy) ~ "BA.2.86",
      TRUE ~ "" 
    )
  )

############### Statistik.fhi.no - calcualte and autocomplete plus add the flagg for publications:

# Load necessary libraries
library(lubridate)
library(dplyr)
library(tidyr)

# Define the period for the last six months
last_six_months <- seq(
  from = floor_date(Sys.Date() %m-% months(6), unit = "month"),
  to = floor_date(Sys.Date(), unit = "month"),
  by = "month"
)

# Format these months consistently as lowercase and with leading zero in month (e.g., "2020 feb")
last_six_months <- format(last_six_months, "%Y %b")
last_six_months <- tolower(last_six_months)

# Process for VUM variant
## Extract unique values from VUM
unique_variants_vum <- na.omit(unique(SC2db$VUM))

## Generate all combinations of months and variants
complete_combined_vum <- expand.grid(
  my = last_six_months, 
  VUM = unique_variants_vum
)

## Process the actual data for VUM
processed_data_vum <- SC2db %>%
  filter(my > yearmonth(Sys.Date() %m-% months(6))) %>%
  filter(!is.na(VUM) & VUM != "") %>%
  mutate(my = format(as.Date(my), "%Y %b") %>% tolower()) %>%
  group_by(my, VUM) %>%
  count(name = "Antall") %>%
  mutate(flagg = 0)

## Join the complete set with actual VUM data
VUM <- complete_combined_vum %>%
  left_join(processed_data_vum, by = c("my", "VUM")) %>%
  replace_na(list(Antall = 0, flagg = 0)) %>%
  filter(VUM != "" & !is.na(VUM)) %>%
  mutate(month_ord = as.Date(paste0(my, " 01"), format = "%Y %b %d")) %>%
  arrange(month_ord, VUM) %>%
  select(-month_ord)

# Process for VOI variant
## Extract unique values from VOI
unique_variants_voi <- na.omit(unique(SC2db$VOI))

## Generate all combinations of months and variants
complete_combined_voi <- expand.grid(
  my = last_six_months, 
  VOI = unique_variants_voi
)

## Process the actual data for VOI
processed_data_voi <- SC2db %>%
  filter(my > yearmonth(Sys.Date() %m-% months(6))) %>%
  filter(!is.na(VOI) & VOI != "") %>%
  mutate(my = format(as.Date(my), "%Y %b") %>% tolower()) %>%
  group_by(my, VOI) %>%
  count(name = "Antall") %>%
  mutate(flagg = 0)

## Join the complete set with actual VOI data
VOI <- complete_combined_voi %>%
  left_join(processed_data_voi, by = c("my", "VOI")) %>%
  replace_na(list(Antall = 0, flagg = 0)) %>%
  filter(VOI != "" & !is.na(VOI)) %>%
  mutate(month_ord = as.Date(paste0(my, " 01"), format = "%Y %b %d")) %>%
  arrange(month_ord, VOI) %>%
  select(-month_ord)

# Build the complete Pangolin lineage table for the most recent six months.
last_six_months <- seq(
  from = floor_date(Sys.Date() %m-% months(6), unit = "month"),
  to = floor_date(Sys.Date(), unit = "month"),
  by = "month"
) |>
  format("%Y %b") |>
  tolower()

unique_variants <- SC2db %>%
  filter(my > yearmonth(Sys.Date() %m-% months(6))) %>%
  filter(!is.na(nc_pangolin_short), nc_pangolin_short != "") %>%
  distinct() %>%
  pull(nc_pangolin_short)

complete_combined <- expand.grid(
  my = last_six_months, 
  nc_pangolin_short = unique_variants
)

# Process the actual data with grouping
processed_data <- SC2db %>%
  filter(my > yearmonth(Sys.Date() %m-% months(6))) %>%
  filter(!is.na(nc_pangolin_short) & nc_pangolin_short != "") %>%
  mutate(my = format(as.Date(my), "%Y %b") %>% tolower()) %>%
  group_by(my, nc_pangolin_short) %>%
  count(name = "Antall") %>%
  mutate(flagg = 0)

# Join the complete set with actual data
Stat <- complete_combined %>%
  left_join(processed_data, by = c("my", "nc_pangolin_short")) %>%
  replace_na(list(Antall = 0, flagg = 0)) %>%
  filter(nc_pangolin_short != "" & !is.na(nc_pangolin_short)) %>%
  mutate(month_ord = as.Date(paste0(my, " 01"), format = "%Y %b %d")) %>%
  arrange(month_ord, nc_pangolin_short) %>%
  select(-month_ord)
