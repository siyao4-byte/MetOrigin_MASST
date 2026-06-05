####################################################################################################
# MAPS + MetOriginDB matching AND unmatched-MS2 MGF export
#
# One-piece R workflow:
#   1) Read MAPS final annotation table
#   2) Read MetOriginDB table
#   3) Match MAPS features to MetOriginDB by:
#        PubChem CID -> SMILES -> InChIKey -> HMDB -> KEGG -> compound name
#   4) Preserve one row per MAPS feature/annotation
#   5) Aggregate multiple MetOriginDB rows per matching key before joining
#   6) Export full / compact / preview / unmatched MetOrigin outputs
#   7) Parse original MGF
#   8) Remove features already matched by MetOriginDB
#   9) Export unmatched features with MS2 as a clean MASST/domainMASST-ready MGF
#
# This script does NOT run MASST.
# You can submit the exported MGF separately in another program.
####################################################################################################


# ==================================================================================================
# 0. Libraries
# ==================================================================================================

required_packages <- c(
  "tidyverse",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "yaml"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(yaml)


# ==================================================================================================
# 1. User settings
# ==================================================================================================
# Adjust these paths for your project.

project_dir <- "C:/Users/siyao4/OneDrive - The University of Melbourne/Documents/MAPS_MetOrigin_Matching_Project/MAPS_MetOrigin_Matching_Project"

maps_file <- file.path(project_dir, "data/maps/final-annotation-df.csv")
metorigin_file <- file.path(project_dir, "data/metorigin/MetOriginDB.csv")
mgf_file <- file.path(project_dir, "data/mgf/data_iimn_gnps.mgf")

out_match_dir <- file.path(project_dir, "outputs/metorigin_matches")
out_qc_dir <- file.path(project_dir, "outputs/qc_reports")
out_mgf_dir <- file.path(project_dir, "outputs/domainmasst_unmatched_ms2")

dir.create(out_match_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_mgf_dir, showWarnings = FALSE, recursive = TRUE)

# Matching behaviour
case_insensitive_name_match <- TRUE

# TRUE values for origin flags and metorigin_matched
yes_values <- c("yes", "y", "true", "t", "1")

# Output file names
main_output_file <- file.path(out_match_dir, "MAPS_with_MetOriginDB_origin.csv")
summary_file <- file.path(out_qc_dir, "metorigin_matching_summary.csv")
unmatched_file <- file.path(out_qc_dir, "metorigin_unmatched_features.csv")
column_report_file <- file.path(out_qc_dir, "metorigin_column_detection_report.txt")

# MGF output files
unmatched_ms2_mgf_file <- file.path(out_mgf_dir, "MetOrigin_unmatched_features_with_MS2.mgf")
unmatched_ms2_manifest_file <- file.path(out_mgf_dir, "MetOrigin_unmatched_features_with_MS2_manifest.csv")
unmatched_ms2_feature_table_file <- file.path(out_mgf_dir, "MetOrigin_unmatched_features_with_MS2.csv")


# ==================================================================================================
# 2. General helper functions
# ==================================================================================================

clean_colname <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_remove("^_") %>%
    str_remove("_$")
}

make_unique_names <- function(x) {
  out <- clean_colname(x)
  make.unique(out, sep = "_")
}

clean_dataframe_columns <- function(df) {
  names(df) <- make_unique_names(names(df))
  df
}

clean_id <- function(x) {
  out <- as.character(x)
  out <- str_trim(out)
  out <- str_replace(out, "\\.0$", "")
  out[out %in% c("", "NA", "N/A", "na", "n/a", "nan", "NaN", "None", "none", "NULL", "null")] <- NA_character_
  out
}

clean_name <- function(x, case_insensitive = TRUE) {
  out <- as.character(x)
  out <- str_replace_all(out, "\\s+", " ")
  out <- str_trim(out)
  out[out %in% c("", "NA", "N/A", "na", "n/a", "nan", "NaN", "None", "none", "NULL", "null")] <- NA_character_
  if (case_insensitive) out <- str_to_lower(out)
  out
}

pick_col <- function(df, candidates) {
  candidates_clean <- clean_colname(candidates)
  hit <- candidates_clean[candidates_clean %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

is_yes <- function(x) {
  str_to_lower(str_trim(as.character(x))) %in% yes_values
}

collapse_yes_unknown <- function(values) {
  vals <- values %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower()
  vals <- vals[!is.na(vals) & vals != ""]
  if (any(vals == "yes")) return("Yes")
  if (any(vals == "unknown")) return("Unknown")
  if (any(vals == "no")) return("No")
  ""
}

collapse_unique_text <- function(values, sep = " | ") {
  vals <- values %>%
    as.character() %>%
    str_trim()
  vals <- vals[!is.na(vals) & vals != "" & !vals %in% c("NA", "NaN", "NULL", "None")]
  if (length(vals) == 0) return("")
  paste(sort(unique(vals)), collapse = sep)
}


# ==================================================================================================
# 3. Column detection
# ==================================================================================================

detect_maps_columns <- function(maps) {
  cols <- list(
    feature_id = pick_col(maps, c("feature_id", "featureid", "id")),
    usi = pick_col(maps, c("feature_usi", "usi")),
    name = pick_col(maps, c("compound_name", "compound.name", "name", "title")),
    smiles = pick_col(maps, c("smiles", "SMILES", "canonical_smiles")),
    inchikey = pick_col(maps, c("inchikey", "inchi_key", "inchi.key")),
    cid = pick_col(maps, c("CID", "cid", "pubchem_cid", "pubchem_id", "pubchem_compound_id")),
    hmdb = pick_col(maps, c("hmdb_id", "primary_hmdb_id", "hmdb")),
    kegg = pick_col(maps, c("kegg_id", "kegg")),
    confidence = pick_col(maps, c("confidence_level", "confidence.level", "annotation_confidence_level")),
    annotation_type = pick_col(maps, c("annotation_type", "annotation.type")),
    mz = pick_col(maps, c("mz")),
    rt = pick_col(maps, c("rt")),
    formula = pick_col(maps, c("formula")),
    iupac = pick_col(maps, c("iupac")),
    mono_mass = pick_col(maps, c("monoisotopic_mass", "monoisotopic.mass", "mass"))
  )

  if ("cid" %in% names(maps)) cols$cid <- "cid"
  if ("smiles" %in% names(maps)) cols$smiles <- "smiles"

  cols
}

detect_metorigin_columns <- function(metorigin) {
  cols <- list(
    name = pick_col(metorigin, c("compound_name", "compound_name_id", "compound.name", "name", "metabolite_name", "title", "chemical_name")),
    smiles = pick_col(metorigin, c("SMILES_ID", "smiles_id", "smiles", "canonical_smiles", "isomeric_smiles", "pubchem_smiles")),
    inchikey = pick_col(metorigin, c("inchikey", "inchi_key", "inchi.key", "standard_inchikey")),
    cid = pick_col(metorigin, c("PUBCHEM_COMPOUND_ID", "pubchem_compound_id", "cid", "pubchem_cid", "pubchem_id", "pubchem")),
    hmdb = pick_col(metorigin, c("hmdb_id", "primary_hmdb_id", "hmdb", "hmdbid")),
    kegg = pick_col(metorigin, c("kegg_id", "kegg", "keggid"))
  )

  if ("pubchem_compound_id" %in% names(metorigin)) cols$cid <- "pubchem_compound_id"
  if ("smiles_id" %in% names(metorigin)) cols$smiles <- "smiles_id"

  cols
}


# ==================================================================================================
# 4. MetOriginDB matching helpers
# ==================================================================================================

origin_flag_columns <- c(
  "from_human",
  "from_bacteria",
  "from_plant",
  "from_animal",
  "from_environment",
  "from_drug",
  "from_food"
)

origin_detail_columns <- c(
  "from_which_part",
  "from_which_bacteria",
  "bacteria_ncbi_id",
  "bacteria_phylum",
  "bacteria_class",
  "bacteria_order",
  "bacteria_family",
  "bacteria_genus",
  "bacteria_species",
  "from_which_plant",
  "from_which_animal",
  "from_which_environment",
  "from_which_drug",
  "from_which_food"
)

origin_columns_to_keep <- c(origin_flag_columns, origin_detail_columns)

metorigin_identity_columns <- c(
  "pubchem_compound_id",
  "smiles_id",
  "compound_name",
  "compound_name_id",
  "hmdb_id",
  "kegg_id",
  "inchikey",
  "mimedb_id",
  "lotus_id"
)

add_match_key <- function(df, source_col, key_col, kind, case_insensitive_names = TRUE) {
  if (is.na(source_col) || !source_col %in% names(df)) return(df)

  if (kind == "name") {
    df[[key_col]] <- clean_name(df[[source_col]], case_insensitive = case_insensitive_names)
  } else {
    df[[key_col]] <- clean_id(df[[source_col]])
  }

  df
}

prepare_match_keys <- function(maps, metorigin, maps_cols, met_cols, case_insensitive_names = TRUE) {
  key_specs <- tribble(
    ~logical,   ~key_col,         ~kind,
    "cid",      "match_cid",      "id",
    "smiles",   "match_smiles",   "id",
    "inchikey", "match_inchikey", "id",
    "hmdb",     "match_hmdb",     "id",
    "kegg",     "match_kegg",     "id",
    "name",     "match_name",     "name"
  )

  for (i in seq_len(nrow(key_specs))) {
    logical <- key_specs$logical[i]
    key_col <- key_specs$key_col[i]
    kind <- key_specs$kind[i]

    maps <- add_match_key(maps, maps_cols[[logical]], key_col, kind, case_insensitive_names)
    metorigin <- add_match_key(metorigin, met_cols[[logical]], key_col, kind, case_insensitive_names)
  }

  list(maps = maps, metorigin = metorigin)
}

build_metorigin_subset <- function(met_db, key) {
  keep_cols <- key

  for (c in origin_columns_to_keep) {
    if (c %in% names(met_db) && !c %in% keep_cols) keep_cols <- c(keep_cols, c)
  }

  for (c in metorigin_identity_columns) {
    if (c %in% names(met_db) && !c %in% keep_cols) keep_cols <- c(keep_cols, c)
  }

  met_tmp <- met_db %>%
    select(any_of(keep_cols)) %>%
    filter(!is.na(.data[[key]]), .data[[key]] != "")

  if (nrow(met_tmp) == 0) return(met_tmp)

  group_cols <- key

  flag_cols_present <- intersect(origin_flag_columns, names(met_tmp))
  detail_cols_present <- intersect(origin_detail_columns, names(met_tmp))
  identity_cols_present <- setdiff(intersect(metorigin_identity_columns, names(met_tmp)), key)

  met_small <- met_tmp %>%
    group_by(.data[[key]]) %>%
    summarise(
      across(all_of(flag_cols_present), collapse_yes_unknown),
      across(all_of(detail_cols_present), collapse_unique_text),
      across(all_of(identity_cols_present), collapse_unique_text),
      metorigin_n_database_records = n(),
      .groups = "drop"
    )

  names(met_small)[names(met_small) == key] <- key

  rename_map <- identity_cols_present
  for (c in rename_map) {
    if (c %in% names(met_small)) {
      names(met_small)[names(met_small) == c] <- paste0("metorigin_", c)
    }
  }

  met_small
}

join_by_key <- function(unmatched, met_db, key, match_method) {
  if (!key %in% names(unmatched) || !key %in% names(met_db)) {
    message("Skipping ", match_method, " because key column is missing.")
    return(list(matched = tibble(), unmatched = unmatched))
  }

  u <- unmatched %>% filter(!is.na(.data[[key]]), .data[[key]] != "")
  no_key <- unmatched %>% filter(is.na(.data[[key]]) | .data[[key]] == "")

  if (nrow(u) == 0) {
    message("Skipping ", match_method, " because no MAPS rows have this key.")
    return(list(matched = tibble(), unmatched = unmatched))
  }

  met_small <- build_metorigin_subset(met_db, key)

  if (nrow(met_small) == 0) {
    message("Skipping ", match_method, " because no MetOriginDB rows have this key.")
    return(list(matched = tibble(), unmatched = unmatched))
  }

  joined <- u %>%
    left_join(met_small, by = key) %>%
    mutate(.matched_now = !is.na(metorigin_n_database_records))

  matched <- joined %>%
    filter(.matched_now) %>%
    select(-.matched_now) %>%
    mutate(
      metorigin_matched = TRUE,
      metorigin_match_method = match_method
    )

  still_unmatched <- joined %>%
    filter(!.matched_now) %>%
    select(any_of(names(unmatched)))

  next_unmatched <- bind_rows(still_unmatched, no_key)

  message("Matched by ", match_method, ": ", nrow(matched))

  list(matched = matched, unmatched = next_unmatched)
}

run_matching <- function(maps, metorigin, maps_cols, met_cols, case_insensitive_names = TRUE) {
  keyed <- prepare_match_keys(maps, metorigin, maps_cols, met_cols, case_insensitive_names)
  maps2 <- keyed$maps
  met2 <- keyed$metorigin

  match_plan <- tribble(
    ~key,              ~method,
    "match_cid",       "PubChem CID",
    "match_smiles",    "SMILES exact",
    "match_inchikey",  "InChIKey",
    "match_hmdb",      "HMDB ID",
    "match_kegg",      "KEGG ID",
    "match_name",      "compound name exact"
  )

  unmatched <- maps2
  matched_parts <- list()

  for (i in seq_len(nrow(match_plan))) {
    res <- join_by_key(unmatched, met2, match_plan$key[i], match_plan$method[i])
    if (nrow(res$matched) > 0) {
      matched_parts[[length(matched_parts) + 1]] <- res$matched
    }
    unmatched <- res$unmatched
  }

  unmatched <- unmatched %>%
    mutate(
      metorigin_matched = FALSE,
      metorigin_match_method = "unmatched"
    )

  final <- bind_rows(matched_parts, unmatched)

  summary <- final %>%
    count(metorigin_matched, metorigin_match_method, name = "n_features") %>%
    arrange(desc(metorigin_matched), desc(n_features))

  list(final = final, summary = summary)
}

add_origin_summary_columns <- function(df) {
  category_map <- c(
    from_human = "human",
    from_bacteria = "bacteria",
    from_plant = "plant",
    from_animal = "animal",
    from_environment = "environment",
    from_drug = "drug",
    from_food = "food"
  )

  for (col in origin_flag_columns) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }

  df %>%
    rowwise() %>%
    mutate(
      metorigin_origin_categories = {
        cats <- c()
        for (col in names(category_map)) {
          if (is_yes(cur_data()[[col]])) cats <- c(cats, category_map[[col]])
        }
        if (length(cats) == 0) "" else paste(cats, collapse = ";")
      },
      metorigin_n_origin_categories = if_else(
        metorigin_origin_categories == "",
        0L,
        length(str_split(metorigin_origin_categories, ";")[[1]])
      )
    ) %>%
    ungroup()
}

drop_internal_match_columns <- function(df) {
  internal_cols <- c("match_cid", "match_smiles", "match_inchikey", "match_hmdb", "match_kegg", "match_name")
  df %>% select(-any_of(internal_cols))
}

make_compact_output <- function(df) {
  df <- add_origin_summary_columns(df)
  df <- drop_internal_match_columns(df)

  preferred_maps_cols <- c(
    "feature_id",
    "feature_usi",
    "usi",
    "mz",
    "rt",
    "compound_name",
    "smiles",
    "formula",
    "iupac",
    "monoisotopic_mass",
    "cid",
    "hmdb_id",
    "annotation_type",
    "confidence_level",
    "confidence_score"
  )

  match_cols <- c(
    "metorigin_matched",
    "metorigin_match_method",
    "metorigin_n_database_records",
    "metorigin_origin_categories",
    "metorigin_n_origin_categories"
  )

  metorigin_identity_cols <- c(
    "metorigin_pubchem_compound_id",
    "metorigin_smiles_id",
    "metorigin_compound_name",
    "metorigin_compound_name_id",
    "metorigin_hmdb_id",
    "metorigin_kegg_id",
    "metorigin_inchikey",
    "metorigin_mimedb_id",
    "metorigin_lotus_id"
  )

  compact_cols <- c(
    preferred_maps_cols,
    match_cols,
    origin_flag_columns,
    metorigin_identity_cols,
    origin_detail_columns
  )

  compact_cols <- compact_cols[compact_cols %in% names(df)]
  remaining_cols <- setdiff(names(df), compact_cols)

  df %>% select(all_of(compact_cols), all_of(remaining_cols))
}

make_excel_friendly_preview <- function(df, max_detail_chars = 300) {
  out <- df

  for (c in origin_detail_columns) {
    if (c %in% names(out)) {
      out[[c]] <- ifelse(
        is.na(out[[c]]) | nchar(as.character(out[[c]])) <= max_detail_chars,
        as.character(out[[c]]),
        paste0(substr(as.character(out[[c]]), 1, max_detail_chars), "...")
      )
    }
  }

  out
}

write_column_detection_report <- function(out_file, maps, metorigin, maps_cols, met_cols) {
  lines <- c(
    "MAPS -> MetOriginDB column detection report",
    strrep("=", 80),
    "",
    "Detected MAPS logical columns",
    strrep("-", 80),
    paste(names(maps_cols), unlist(maps_cols), sep = ": "),
    "",
    "Detected MetOriginDB logical columns",
    strrep("-", 80),
    paste(names(met_cols), unlist(met_cols), sep = ": "),
    "",
    "All MAPS columns after cleaning",
    strrep("-", 80),
    names(maps),
    "",
    "All MetOriginDB columns after cleaning",
    strrep("-", 80),
    names(metorigin)
  )

  writeLines(lines, out_file)
}


# ==================================================================================================
# 5. MGF parsing and unmatched MGF export
# ==================================================================================================

extract_feature_id_from_text <- function(text) {
  if (is.na(text) || text == "") return(NA_character_)

  patterns <- c(
    "feature_id[=:| _-]*([A-Za-z0-9_.-]+)",
    "feature\\.id[=:| _-]*([A-Za-z0-9_.-]+)",
    "featureID[=:| _-]*([A-Za-z0-9_.-]+)",
    "Feature[_ -]?([0-9]+)"
  )

  for (pat in patterns) {
    m <- str_match(text, regex(pat, ignore_case = TRUE))
    if (!is.na(m[1, 2])) return(as.character(m[1, 2]))
  }

  NA_character_
}

get_mgf_field <- function(block_lines, key) {
  hit <- block_lines[str_detect(str_to_upper(block_lines), paste0("^", str_to_upper(key), "="))]
  if (length(hit) == 0) return(NA_character_)
  str_replace(hit[[1]], paste0("^", key, "="), "")
}

set_mgf_field <- function(block_lines, key, value) {
  key_pat <- paste0("^", str_to_upper(key), "=")
  idx <- which(str_detect(str_to_upper(block_lines), key_pat))

  if (length(idx) > 0) {
    block_lines[idx[[1]]] <- paste0(key, "=", value)
  } else {
    insert_at <- ifelse(length(block_lines) > 0 && str_to_upper(str_trim(block_lines[[1]])) == "BEGIN IONS", 2, 1)
    block_lines <- append(block_lines, paste0(key, "=", value), after = insert_at - 1)
  }

  block_lines
}

parse_mgf_blocks <- function(mgf_path) {
  lines <- readLines(mgf_path, warn = FALSE)

  begin_idx <- which(str_to_upper(str_trim(lines)) == "BEGIN IONS")
  end_idx <- which(str_to_upper(str_trim(lines)) == "END IONS")

  if (length(begin_idx) != length(end_idx)) {
    stop("MGF parsing error: number of BEGIN IONS and END IONS lines does not match.")
  }

  map2_dfr(begin_idx, end_idx, function(i, j) {
    block <- lines[i:j]
    title <- get_mgf_field(block, "TITLE")
    scans <- get_mgf_field(block, "SCANS")

    feature_id <- extract_feature_id_from_text(title)
    if (is.na(feature_id)) feature_id <- extract_feature_id_from_text(scans)
    if (is.na(feature_id) && !is.na(scans) && scans != "") feature_id <- scans

    tibble(
      feature_id = as.character(feature_id),
      original_title = title,
      original_scans = scans,
      pepmass = get_mgf_field(block, "PEPMASS"),
      rtinseconds = get_mgf_field(block, "RTINSECONDS"),
      charge = get_mgf_field(block, "CHARGE"),
      block = list(block)
    )
  })
}

export_unmatched_ms2_mgf <- function(final_compact, mgf_path, out_mgf, out_manifest, out_unmatched_table) {
  if (!"feature_id" %in% names(final_compact)) {
    stop("final_compact does not contain feature_id.")
  }

  if (!"metorigin_matched" %in% names(final_compact)) {
    stop("final_compact does not contain metorigin_matched.")
  }

  unmatched_features <- final_compact %>%
    mutate(feature_id = as.character(feature_id)) %>%
    filter(!is_yes(metorigin_matched)) %>%
    distinct(feature_id, .keep_all = TRUE)

  unmatched_ids <- unmatched_features$feature_id

  message("MetOrigin-unmatched feature IDs retained for MGF check: ", length(unmatched_ids))

  mgf_blocks <- parse_mgf_blocks(mgf_path) %>%
    filter(!is.na(feature_id), feature_id %in% unmatched_ids)

  message("Unmatched MS2 spectra exported: ", nrow(mgf_blocks))
  message("Unmatched features with MS2: ", n_distinct(mgf_blocks$feature_id))

  if (nrow(mgf_blocks) == 0) {
    warning("No unmatched MS2 spectra found in MGF. Check feature_id format in MGF TITLE/SCANS.")
  }

  rewritten_blocks <- pmap(
    list(mgf_blocks$feature_id, mgf_blocks$original_title, mgf_blocks$block),
    function(feature_id, original_title, block) {
      if (is.na(original_title)) original_title <- ""
      new_title <- paste0("feature_id=", feature_id, "|original_title=", original_title)
      block <- set_mgf_field(block, "TITLE", new_title)
      block <- set_mgf_field(block, "SCANS", feature_id)
      block
    }
  )

  mgf_out_lines <- unlist(map(rewritten_blocks, function(x) c(x, "")), use.names = FALSE)
  writeLines(mgf_out_lines, out_mgf, useBytes = TRUE)

  manifest <- mgf_blocks %>%
    mutate(
      new_title = paste0("feature_id=", feature_id, "|original_title=", if_else(is.na(original_title), "", original_title)),
      scans = feature_id
    ) %>%
    select(
      feature_id,
      original_title,
      new_title,
      scans,
      original_scans,
      pepmass,
      rtinseconds,
      charge
    )

  write_csv(manifest, out_manifest)

  unmatched_with_ms2 <- unmatched_features %>%
    inner_join(
      manifest %>%
        group_by(feature_id) %>%
        summarise(
          n_ms2_spectra = n(),
          masst_query_titles = paste(unique(new_title), collapse = " || "),
          precursor_mz_mgf = paste(unique(na.omit(pepmass)), collapse = " || "),
          rt_seconds_mgf = paste(unique(na.omit(rtinseconds)), collapse = " || "),
          .groups = "drop"
        ),
      by = "feature_id"
    )

  write_csv(unmatched_with_ms2, out_unmatched_table)

  list(
    manifest = manifest,
    unmatched_with_ms2 = unmatched_with_ms2,
    out_mgf = out_mgf
  )
}


# ==================================================================================================
# 6. Run workflow
# ==================================================================================================

message("Starting MAPS -> MetOriginDB matching and unmatched-MS2 MGF export.")
message("MAPS file: ", maps_file)
message("MetOriginDB file: ", metorigin_file)
message("MGF file: ", mgf_file)

if (!file.exists(maps_file)) stop("MAPS file not found: ", maps_file)
if (!file.exists(metorigin_file)) stop("MetOriginDB file not found: ", metorigin_file)
if (!file.exists(mgf_file)) stop("MGF file not found: ", mgf_file)

maps <- read_csv(maps_file, show_col_types = FALSE, guess_max = 100000) %>%
  clean_dataframe_columns()

metorigin <- read_csv(metorigin_file, show_col_types = FALSE, guess_max = 100000) %>%
  clean_dataframe_columns()

message("MAPS rows: ", nrow(maps), " | columns: ", ncol(maps))
message("MetOriginDB rows: ", nrow(metorigin), " | columns: ", ncol(metorigin))

maps_cols <- detect_maps_columns(maps)
met_cols <- detect_metorigin_columns(metorigin)

write_column_detection_report(column_report_file, maps, metorigin, maps_cols, met_cols)

match_res <- run_matching(
  maps = maps,
  metorigin = metorigin,
  maps_cols = maps_cols,
  met_cols = met_cols,
  case_insensitive_names = case_insensitive_name_match
)

final_full <- add_origin_summary_columns(match_res$final)
final_full_no_helpers <- drop_internal_match_columns(final_full)
final_compact <- make_compact_output(final_full)
final_preview <- make_excel_friendly_preview(final_compact)

full_path <- file.path(out_match_dir, "MAPS_with_MetOriginDB_origin_FULL.csv")
compact_path <- file.path(out_match_dir, "MAPS_with_MetOriginDB_origin_COMPACT.csv")
preview_path <- file.path(out_match_dir, "MAPS_with_MetOriginDB_origin_EXCEL_PREVIEW_truncated_details.csv")

write_csv(final_compact, main_output_file)
write_csv(final_compact, compact_path)
write_csv(final_full_no_helpers, full_path)
write_csv(final_preview, preview_path)
write_csv(match_res$summary, summary_file)

unmatched <- final_compact %>%
  filter(!is_yes(metorigin_matched))

write_csv(unmatched, unmatched_file)

message("Saved main compact output: ", main_output_file)
message("Saved compact output: ", compact_path)
message("Saved full output: ", full_path)
message("Saved Excel-friendly preview: ", preview_path)
message("Saved summary: ", summary_file)
message("Saved unmatched features: ", unmatched_file)
message("Saved column detection report: ", column_report_file)

message("Matching summary:")
print(match_res$summary)

mgf_res <- export_unmatched_ms2_mgf(
  final_compact = final_compact,
  mgf_path = mgf_file,
  out_mgf = unmatched_ms2_mgf_file,
  out_manifest = unmatched_ms2_manifest_file,
  out_unmatched_table = unmatched_ms2_feature_table_file
)

message("Saved unmatched-MS2 MGF: ", unmatched_ms2_mgf_file)
message("Saved unmatched-MS2 manifest: ", unmatched_ms2_manifest_file)
message("Saved unmatched-MS2 feature table: ", unmatched_ms2_feature_table_file)

message("Done.")
