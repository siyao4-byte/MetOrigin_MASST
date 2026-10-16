# R compatibility entry point; one authoritative matching implementation.
config_path <- Sys.getenv("METORIGIN_CONFIG", "config/config.yaml")
python <- Sys.getenv("METORIGIN_PYTHON", ".venv310/Scripts/python.exe")
status <- system2(python, c("scripts/metorigin_match.py", "--config", shQuote(config_path)))
if (status != 0) stop("Indexed matching failed; no legacy or SID fallback is permitted.")
cfg <- yaml::read_yaml(config_path)
if (identical(Sys.getenv("METORIGIN_EXPORT_MGF", "false"), "true")) {
  full <- readr::read_csv(cfg$output_file, show_col_types = FALSE)
  if (!all(full$feature_id == full$original_feature_id)) {
    stop("MGF export requires unique original feature IDs; total-list identities cannot be mapped by scan number alone.")
  }
  status <- system2(python, c("scripts/prepare_unmatched_ms2_domainmasst_merge.py",
    "--metorigin_csv", shQuote(cfg$output_file), "--mgf", shQuote(cfg$domainmasst$input_mgf),
    "--out_dir", shQuote(cfg$domainmasst$output_dir)))
  if (status != 0) stop("Optional MGF preparation failed.")
}
message("Recorded source associations exported; domainMASST search was not submitted.")
