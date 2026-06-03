# INF Tessy cases 25-26

options(repos = c(CRAN = "https://cran.uni-muenster.de/"))

resolve_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path, winslash = "/", mustWork = FALSE)))
    }
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

check_install_update_packages <- function(packages) {
  installed_packages <- rownames(installed.packages())
  missing_packages <- setdiff(packages, installed_packages)
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

required_packages <- c("dplyr", "lubridate", "tidyr", "stringr", "magrittr")
check_install_update_packages(required_packages)
suppressPackageStartupMessages(lapply(required_packages, library, character.only = TRUE))

source(file.path(bundle_scripts_dir, "INF_SQLquery_25-26.R"))
source(file.path(bundle_scripts_dir, "INF_DataCleaning_25-26.R"))

if (exists("INF_25_26_clean")) {
  fludb <- INF_25_26_clean
}

if (!exists("fludb")) {
  stop("Object 'fludb' is missing. Source INF_DataCleaning_25-26.R first.")
}

season_start_date <- as.Date(
  sprintf("%d-10-01", ifelse(month(Sys.Date()) >= 10, year(Sys.Date()), year(Sys.Date()) - 1))
)

required_cols <- c(
  "ngs_sekvens_resultat", "pasient_kjnn", "prove_tatt", "pasient_alder",
  "tessy_reportable_variable", "gisaid_isolate_name", "pasient_status",
  "pasient_vaks", "prove_kategori", "mut_ha1_1", "gisaid_ha_id",
  "comment_gc_tessy", "gisaid_na_id", "gisaid_m_id", "gisaid_pa_id",
  "dr_res_adamantine", "dr_res_oseltamivir", "dr_res_baloxavir",
  "dr_res_zanamivir", "dr_pa_mut", "gisaid_isolate_id",
  "dr_oseltamivir_ic50", "dr_zanamivir_ic50", "dr_na_mut", "dr_m2_mut"
)

missing_cols <- setdiff(required_cols, names(fludb))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in fludb:", paste(missing_cols, collapse = ", ")))
}

normalize_tessy_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- stringr::str_replace_all(x, "[^[:print:]]", "")
  stringr::str_to_upper(stringr::str_squish(x))
}

warn_unmapped_values <- function(data, column, allowed_values, label = column) {
  values <- unique(na.omit(normalize_tessy_text(data[[column]])))
  values <- values[values != ""]
  unknown_values <- setdiff(values, allowed_values)
  if (length(unknown_values) > 0) {
    message(
      "Unmapped ", label, " values: ",
      paste(head(sort(unknown_values), 20), collapse = ", ")
    )
  }
  invisible(unknown_values)
}

warn_unmapped_values(
  fludb,
  "pasient_status",
  c("INNELIGGENDE", "INN", "INNE", "POLIKLINISK", "POL", "IKKE_SATT", "IKKE SATT", "UKJENT", "U", "UNK", "")
)
warn_unmapped_values(
  fludb,
  "pasient_vaks",
  c("J", "JA", "YES", "Y", "N", "NEI", "NO", "X", "UKJENT", "UNK", "U", "")
)
warn_unmapped_values(
  fludb,
  "pasient_kjnn",
  c("M", "MALE", "F", "FEMALE", "N", "UKJENT", "UNKNOWN", "UNK", "")
)

tessy_source <- fludb %>%
  mutate(
    ngs_sekvens_resultat_norm = normalize_tessy_text(ngs_sekvens_resultat),
    pasient_status_norm = normalize_tessy_text(pasient_status),
    pasient_vaks_norm = normalize_tessy_text(pasient_vaks),
    pasient_kjnn_norm = normalize_tessy_text(pasient_kjnn),
    prove_kategori_norm = normalize_tessy_text(prove_kategori),
    Sampledate = as.Date(prove_tatt, format = "%Y-%m-%d"),
    age = pasient_alder
  ) %>%
  mutate(
    Subtype = case_when(
      str_detect(ngs_sekvens_resultat_norm, "^A/H1N1") ~ "AH1N1pdm09",
      str_detect(ngs_sekvens_resultat_norm, "^A/H3N2") ~ "AH3N2",
      str_detect(ngs_sekvens_resultat_norm, "N2") ~ "AH3N2",
      str_detect(ngs_sekvens_resultat_norm, "^B/VICTORIA") ~ "BVic",
      str_detect(ngs_sekvens_resultat_norm, "^B/YAMAGATA") ~ "BYam",
      str_detect(ngs_sekvens_resultat_norm, "N1") ~ "AH1N1pdm09",
      str_detect(ngs_sekvens_resultat_norm, "H1") ~ "AH1pdm09",
      str_detect(ngs_sekvens_resultat_norm, "H3") ~ "AH3",
      TRUE ~ NA_character_
    ),
    INFType = case_when(
      str_detect(Subtype, "^A") ~ "A",
      str_detect(Subtype, "^B") ~ "B",
      TRUE ~ NA_character_
    ),
    Host_Gender = case_when(
      pasient_kjnn_norm %in% c("M", "MALE") ~ "M",
      pasient_kjnn_norm %in% c("F", "FEMALE") ~ "F",
      pasient_kjnn_norm %in% c("N", "U", "UKJENT", "UNKNOWN", "UNK", "") ~ "UNK",
      TRUE ~ "UNK"
    ),
    hospitalisation = case_when(
      pasient_status_norm %in% c("INNELIGGENDE", "INN", "INNE") ~ "Y",
      pasient_status_norm %in% c("POLIKLINISK", "POL") ~ "N",
      pasient_status_norm %in% c("IKKE_SATT", "IKKE SATT", "UKJENT", "U", "UNK", "") | is.na(pasient_status) ~ "UNK",
      TRUE ~ "UNK"
    ),
    vacstatus = case_when(
      pasient_vaks_norm %in% c("J", "JA", "YES", "Y") ~ "VACCINFULL",
      pasient_vaks_norm %in% c("N", "NEI", "NO") ~ "NOTVACC",
      pasient_vaks_norm %in% c("X", "UKJENT", "UNK", "U", "") | is.na(pasient_vaks) ~ "UNK",
      TRUE ~ "UNK"
    ),
    virussource = case_when(
      is.na(prove_kategori_norm) | prove_kategori_norm == "" ~ "UNK",
      str_detect(prove_kategori_norm, "^P1(?:$|_|[[:alpha:]].*|\\W.*)") ~ "S",
      TRUE ~ "N"
    ),
    vicat = case_when(
      virussource == "S" ~ "",
      virussource == "N" & hospitalisation == "Y" ~ "H",
      virussource == "N" & hospitalisation == "N" ~ "C",
      TRUE ~ "UNK"
    ),
    drop_reason = case_when(
      is.na(ngs_sekvens_resultat_norm) | ngs_sekvens_resultat_norm == "" ~ "blank_ngs_sekvens_resultat",
      is.na(Sampledate) ~ "missing_sampledate",
      Sampledate < season_start_date ~ "before_season_start",
      is.na(tessy_reportable_variable) | stringr::str_trim(tessy_reportable_variable) == "" ~ "missing_tessy_reportable_variable",
      is.na(gisaid_isolate_name) | stringr::str_trim(gisaid_isolate_name) == "" ~ "missing_gisaid_isolate_name",
      is.na(Subtype) ~ "unmapped_subtype",
      TRUE ~ NA_character_
    )
  )

message(
  "Tessy QC counts - total: ", nrow(tessy_source),
  "; kept after validation: ", sum(is.na(tessy_source$drop_reason))
)
if (any(!is.na(tessy_source$drop_reason))) {
  drop_summary <- tessy_source %>%
    filter(!is.na(drop_reason)) %>%
    count(drop_reason, sort = TRUE)
  message("Dropped rows by reason:")
  for (i in seq_len(nrow(drop_summary))) {
    message("  ", drop_summary$drop_reason[i], ": ", drop_summary$n[i])
  }
}

dropped_rows_table <- tessy_source %>%
  filter(!is.na(drop_reason)) %>%
  transmute(
    RecordId = gisaid_isolate_name,
    Key = key,
    SampleDate = Sampledate,
    ProveTatt = prove_tatt,
    Subtype,
    INFType,
    VirusSource = virussource,
    Hospitalisation = hospitalisation,
    VaccStatus = vacstatus,
    Gender = Host_Gender,
    NgsResult = ngs_sekvens_resultat,
    ProveKategori = prove_kategori,
    TessyReportableVariable = tessy_reportable_variable,
    DropReason = drop_reason
  )

is_blank_value <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  stringr::str_trim(x) == "" | x == "NA"
}

extract_gene_qc_token <- function(qc_sum, gene_label) {
  qc_sum <- as.character(qc_sum)
  qc_sum[is.na(qc_sum)] <- ""
  token <- stringr::str_match(qc_sum, paste0("(?:^|\\|)", gene_label, ":([^|]+)"))[, 2]
  token[is.na(token)] <- ""
  stringr::str_trim(token)
}

build_non_submission_comment <- function(gene, qc_token) {
  mapply(
    FUN = function(one_gene, one_token) {
      one_gene <- as.character(one_gene)
      one_token <- as.character(one_token)
      one_token[is.na(one_token)] <- ""
      one_token <- stringr::str_trim(one_token)
      if (one_token == "") {
        return("")
      }

      reason_parts <- character(0)
      if (stringr::str_detect(one_token, "LC")) {
        reason_parts <- c(reason_parts, "low coverage")
      }
      if (stringr::str_detect(one_token, "MS")) {
        reason_parts <- c(reason_parts, "mixed sites")
      }
      if (stringr::str_detect(one_token, "FS")) {
        reason_parts <- c(reason_parts, "frame shifts")
      }
      if (length(reason_parts) == 0) {
        reason_parts <- c(reason_parts, paste0("QC issue ", one_token))
      }

      paste0(
        one_gene,
        " - not submitted to GISAID due to ",
        paste(unique(reason_parts), collapse = " / "),
        " (",
        one_token,
        ")"
      )
    },
    gene,
    qc_token,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
}

gene_qc_rows <- tibble::tribble(
  ~Gene, ~EpiCol, ~MutCol, ~ResCol, ~QcGene,
  "HA", "gisaid_ha_id", "mut_ha1_1", NA_character_, "HA",
  "NA", "gisaid_na_id", "dr_na_mut", "dr_res_oseltamivir", "NA",
  "PA", "gisaid_pa_id", "dr_pa_mut", "dr_res_baloxavir", "PA",
  "PB1", "gisaid_pb1_id", "mut_pb1", NA_character_, "PB1",
  "PB2", "gisaid_pb2_id", "mut_pb2", NA_character_, "PB2",
  "M", "gisaid_m_id", "dr_m2_mut", "dr_res_adamantine", "M"
)

gene_qc_audit_detail <- purrr::pmap_dfr(
  gene_qc_rows,
  function(Gene, EpiCol, MutCol, ResCol, QcGene) {
    epi_missing <- is_blank_value(tessy_source[[EpiCol]])
    mut_present <- !is_blank_value(tessy_source[[MutCol]])
    res_present <- if (is.na(ResCol)) {
      rep(NA, nrow(tessy_source))
    } else {
      !is_blank_value(tessy_source[[ResCol]])
    }
    gene_token <- extract_gene_qc_token(tessy_source$ngs_qc_sum, QcGene)
    gene_token_present <- gene_token != ""
    gene_has_lc <- stringr::str_detect(gene_token, "LC")
    gene_has_ms <- stringr::str_detect(gene_token, "MS")
    gene_has_fs <- stringr::str_detect(gene_token, "FS")
    keep_after_tessy <- is.na(tessy_source$drop_reason)
    needs_gc_comment <- epi_missing & keep_after_tessy & gene_token_present & (gene_has_lc | gene_has_ms | gene_has_fs)

    tibble::tibble(
      Gene = Gene,
      RecordId = tessy_source$gisaid_isolate_name,
      Key = tessy_source$key,
      GisaidIsolateName = tessy_source$gisaid_isolate_name,
      SampleDate = tessy_source$Sampledate,
      ProveTatt = tessy_source$prove_tatt,
      NgsResult = tessy_source$ngs_sekvens_resultat,
      ProveKategori = tessy_source$prove_kategori,
      TessyDropReason = tessy_source$drop_reason,
      GeneEpiMissing = epi_missing,
      GeneEpiPresent = !epi_missing,
      GeneEpiStatus = ifelse(epi_missing, "missing", "present"),
      GeneMutationPresent = mut_present,
      GeneResistancePresent = res_present,
      GeneQcToken = gene_token,
      GeneQcReason = ifelse(gene_token == "", NA_character_, gene_token),
      GeneQcHasLC = gene_has_lc,
      GeneQcHasMS = gene_has_ms,
      GeneQcHasFS = gene_has_fs,
      SurvivesTessyFilters = keep_after_tessy,
      NeedsGcComment = needs_gc_comment
    )
  }
)

gene_qc_audit_summary <- gene_qc_audit_detail %>%
  group_by(Gene) %>%
  summarise(
    MissingEpiTotal = sum(GeneEpiMissing, na.rm = TRUE),
    SurvivesTessy = sum(GeneEpiMissing & SurvivesTessyFilters, na.rm = TRUE),
    DroppedByDateOrCleaning = sum(GeneEpiMissing & !SurvivesTessyFilters, na.rm = TRUE),
    DroppedBlankNgs = sum(GeneEpiMissing & TessyDropReason == "blank_ngs_sekvens_resultat", na.rm = TRUE),
    DroppedBeforeSeason = sum(GeneEpiMissing & TessyDropReason == "before_season_start", na.rm = TRUE),
    DroppedMissingId = sum(GeneEpiMissing & TessyDropReason == "missing_gisaid_isolate_name", na.rm = TRUE),
    HasGeneQcToken = sum(GeneEpiMissing & GeneQcToken != "", na.rm = TRUE),
    HasGeneLcMsFs = sum(GeneEpiMissing & GeneQcToken != "" & (GeneQcHasLC | GeneQcHasMS | GeneQcHasFS), na.rm = TRUE),
    NeedsGcComment = sum(NeedsGcComment, na.rm = TRUE),
    LeftWithoutGeneToken = sum(GeneEpiMissing & GeneQcToken == "" & SurvivesTessyFilters, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Gene)

gene_gc_comment_map <- gene_qc_audit_detail %>%
  filter(
    GeneEpiMissing,
    (GeneMutationPresent | GeneResistancePresent),
    GeneQcToken != "",
    GeneQcHasLC | GeneQcHasMS | GeneQcHasFS
  ) %>%
  mutate(GeneGcComment = build_non_submission_comment(Gene, GeneQcToken)) %>%
  group_by(RecordId) %>%
  summarise(
    GeneGcComment = paste(unique(GeneGcComment[GeneGcComment != ""]), collapse = " | "),
    .groups = "drop"
  )

fludb <- tessy_source %>%
  filter(is.na(drop_reason)) %>%
  select(
    -ngs_sekvens_resultat_norm, -pasient_status_norm, -pasient_vaks_norm,
    -pasient_kjnn_norm, -prove_kategori_norm, -drop_reason
  )

submission <- fludb %>%
  transmute(
    RecordId = gisaid_isolate_name,
    DateUsedForStatistics = Sampledate,
    Age = age,
    Subtype,
    GeneticClade = tessy_reportable_variable,
    Gender = Host_Gender,
    HAAAMutations = ifelse(is.na(mut_ha1_1), "", mut_ha1_1),
    HAISD = gisaid_ha_id,
    Hospitalisation = hospitalisation,
    CommentGC = comment_gc_tessy,
    VaccStatus = vacstatus,
    VirusSource = virussource,
    VirusCategoryIfNonSentinel = vicat,
    NAISD = gisaid_na_id,
    M2ISD = gisaid_m_id,
    PAISD = gisaid_pa_id,
    InterpretationM2BlockerResistanceTesting = dr_res_adamantine,
    InterpretationOseltamivirResistanceTesting = dr_res_oseltamivir,
    InterpretationPABlockerResistanceTesting = dr_res_baloxavir,
    InterpretationZanamivirResistanceTesting = dr_res_zanamivir,
    NAAAMutations = dr_na_mut,
    PAAAMutations = dr_pa_mut,
    M2AAMutations = dr_m2_mut,
    SequenceId = gisaid_isolate_id,
    OseltamivirMUNANA = dr_oseltamivir_ic50,
    ZanamivirMUNANA = dr_zanamivir_ic50,
    RecordType = "INFLANTIVIR",
    RecordTypeVersion = "9",
    Subject = "INFLANTIVIR",
    Status = "NEW/UPDATE",
    DataSource = "NO-FLUVIR",
    ReportingCountry = "NO",
    AgeMonth = "NA",
    Amantadine = "",
    AntigenicGroup = "",
    Comment = "",
    CommentAG = "",
    DateOfOnset = "UNK",
    ExposureDrug2weeksHouse = "UNK",
    ExposureDrug2weeksHouseType = "UNK",
    ExposureDrug2weeksPatient = "UNK",
    ExposureDrug2weeksPatientType = "UNK",
    ImmunoCompromised = "UNK",
    IMOVE = "",
    OseltamivirNAStar = "",
    Outcome = "UNK",
    ProbableCountryOfInfection = "",
    Progress4weeks = "",
    Rimantadine = "",
    ZanamivirNAStar = ""
  )

submission <- submission %>%
  left_join(gene_gc_comment_map, by = "RecordId") %>%
  mutate(
    CommentGC = case_when(
      !is.na(GeneGcComment) & GeneGcComment != "" & !is.na(CommentGC) & CommentGC != "" ~ paste(CommentGC, GeneGcComment, sep = " | "),
      !is.na(GeneGcComment) & GeneGcComment != "" ~ GeneGcComment,
      TRUE ~ CommentGC
    )
  ) %>%
  select(-GeneGcComment)

NAresH1 <- "S110;I117;E119;Q136;G147;R152;D199;I223;S247;H275;R293;N295;Q313;I427;I436;P458"
NAresH3 <- "E119;Q136;T148;D151;I222;R224;N245;N245-;A246-;T247-;G248-;K249-;A250-;K249;E276;R292;N294;N329;S331;S334;R371;Q391"
NAresBVIC <- "H101;G104;E105;T106;G108;E117;H134;Q138;P139;G140;Y142;N144;G145;T146;R150;N151;K152;P165;N169;K186;D197;A200;I221;A245;S246;G247;I262;H273;R292;N294;K360;I361;R374;A395;L396;G407;D432;H439;M464"

PAresH1 <- "E23; A36; I38; E119; E198; K34; A37; E199"
PAresH3 <- "L28; E23; K34; A36; A37; I38; E119; E198; E199"
PAresBVIC <- "E23; M34; F36; I38; E120; G199"

MresH1 <- "L26; V27; A30; S31; G34; L38"
MresH3 <- "L26; V27; A30; S31; G34; L38"
MresBVIC <- ""

submission <- submission %>%
  mutate(
    NAAAMutations = case_when(
      Subtype == "AH1N1pdm09" & str_trim(NAAAMutations) == "No Mutations" ~ NAresH1,
      Subtype == "AH3N2" & str_trim(NAAAMutations) == "No Mutations" ~ NAresH3,
      Subtype == "BVic" & str_trim(NAAAMutations) == "No Mutations" ~ NAresBVIC,
      TRUE ~ NAAAMutations
    ),
    PAAAMutations = case_when(
      Subtype == "AH1N1pdm09" & str_trim(PAAAMutations) == "No Mutations" ~ PAresH1,
      Subtype == "AH3N2" & str_trim(PAAAMutations) == "No Mutations" ~ PAresH3,
      Subtype == "BVic" & str_trim(PAAAMutations) == "No Mutations" ~ PAresBVIC,
      TRUE ~ PAAAMutations
    ),
    M2AAMutations = case_when(
      Subtype == "AH1N1pdm09" & str_trim(M2AAMutations) == "No Mutations" ~ MresH1,
      Subtype == "AH3N2" & str_trim(M2AAMutations) == "No Mutations" ~ MresH3,
      Subtype == "BVic" & str_trim(M2AAMutations) == "No Mutations" ~ MresBVIC,
      TRUE ~ M2AAMutations
    )
  )

columns_to_fix <- c("CommentGC", "OseltamivirMUNANA", "ZanamivirMUNANA", "PAISD", "M2ISD", "HAISD", "NAISD", "SequenceId")

submission <- submission %>%
  mutate(across(all_of(columns_to_fix), ~ replace(., is.na(.) | . == "NA", "")))

submission <- submission %>%
  mutate(
    HAAAMutations = case_when(
      HAAAMutations %in% c("No Sequence", "No Sequence/Partial Sequence", "NA") ~ "",
      TRUE ~ HAAAMutations
    ),
    NAAAMutations = case_when(
      NAAAMutations %in% c("No Sequence", "No Sequence/Partial Sequence", "NA") ~ "",
      TRUE ~ NAAAMutations
    ),
    PAAAMutations = case_when(
      PAAAMutations %in% c("No Sequence", "No Sequence/Partial Sequence", "NA") ~ "",
      TRUE ~ PAAAMutations
    ),
    M2AAMutations = case_when(
      M2AAMutations %in% c("No Sequence", "No Sequence/Partial Sequence", "NA") ~ "",
      TRUE ~ M2AAMutations
    )
  )

submission <- submission %>%
  mutate(
    M2ISD = if_else(
      !is.na(M2ISD) & str_trim(as.character(M2ISD)) != "" &
        (is.na(M2AAMutations) | str_trim(as.character(M2AAMutations)) == ""),
      "",
      as.character(M2ISD)
    )
  )

if (any(submission$Subtype == "BYam", na.rm = TRUE)) {
  message("BYam cases detected. Mutation default substitutions are not defined for B/Yamagata, so those rows are retained with raw mutation values.")
}

UKNcols <- c("Gender", "Hospitalisation", "VirusCategoryIfNonSentinel")
NAcols <- c(
  "InterpretationM2BlockerResistanceTesting",
  "InterpretationPABlockerResistanceTesting",
  "InterpretationOseltamivirResistanceTesting",
  "InterpretationM2BlockerResistanceTesting",
  "InterpretationZanamivirResistanceTesting"
)
Commacols <- c("OseltamivirMUNANA", "ZanamivirMUNANA")

for (col in UKNcols) {
  submission[[col]] <- ifelse(submission[[col]] == "" | submission[[col]] == "NA" | is.na(submission[[col]]), "UNK", submission[[col]])
}

for (col in NAcols) {
  submission[[col]] <- ifelse(submission[[col]] == "", "NA", submission[[col]])
}

for (col in Commacols) {
  submission[[col]] <- gsub(",", ".", as.character(submission[[col]]), fixed = TRUE)
  submission[[col]] <- gsub("^'+", "", submission[[col]])
}

output_dir <- "N:/Virologi/Influensa/2526/WGS_Analyse/ECDC_rapportering/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_filename <- paste0("Resist2TESSy1- ", format(Sys.Date(), "%U-%Y"), ".csv")
output_path <- file.path(output_dir, output_filename)
drop_audit_filename <- paste0("Resist2TESSy1- ", format(Sys.Date(), "%U-%Y"), "_dropped_rows.csv")
drop_audit_path <- file.path(output_dir, drop_audit_filename)
gene_qc_audit_summary_path <- file.path(
  output_dir,
  paste0("Resist2TESSy1- ", format(Sys.Date(), "%U-%Y"), "_gene_qc_summary.csv")
)
gene_qc_audit_detail_path <- file.path(
  output_dir,
  paste0("Resist2TESSy1- ", format(Sys.Date(), "%U-%Y"), "_gene_qc_detail.csv")
)

if (nrow(dropped_rows_table) > 0) {
  write.csv2(dropped_rows_table, drop_audit_path, row.names = FALSE)
  message("Wrote dropped-row audit table: ", drop_audit_path)
}

write.csv2(gene_qc_audit_summary, gene_qc_audit_summary_path, row.names = FALSE)
write.csv2(gene_qc_audit_detail, gene_qc_audit_detail_path, row.names = FALSE)
message("Wrote gene QC summary: ", gene_qc_audit_summary_path)
message("Wrote gene QC detail: ", gene_qc_audit_detail_path)

tryCatch(
  write.csv2(submission, output_path, row.names = FALSE),
  error = function(e) {
    fallback_filename <- paste0(
      "Resist2TESSy1- ",
      format(Sys.Date(), "%U-%Y"),
      "_WIP_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".csv"
    )
    fallback_path <- file.path(output_dir, fallback_filename)
    write.csv2(submission, fallback_path, row.names = FALSE)
    message("Primary output file is locked. Wrote fallback file: ", fallback_path)
  }
)
