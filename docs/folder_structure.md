# Folder Structure

```text
MAPS_MetOrigin_Matching_Project/
  config/
    config.yaml
  data/
    maps/final-annotation-df.csv
    metorigin/MetOriginDB.csv
    mgf/data_iimn_gnps.mgf
  outputs/
    metorigin_matches/
    qc_reports/
    logs/
    domainmasst_unmatched_ms2/
  Graph1_MetOrigin_R_outputs/
  Graph2_domainMASST_R_outputs/
  Graph3_integrated_origin_R_outputs/
  runs/
    HGMD_1047/
      data/
      outputs/
      Graph1_MetOrigin_R_outputs/
      Graph2_domainMASST_R_outputs/
      Graph3_integrated_origin_R_outputs/
  scripts/
    setup_and_run_workflow.py
    metorigin_match.py
    run_domainmasst_unmatched.ps1
    Graph_All_compact_workflow.R
```

## Input Modes

### HGMD Mode

Use `scripts/setup_and_run_workflow.py --hgmd HGMD_xxxx` to fetch inputs from:

```text
Y:\MA_BPA_Microbiome\LCMS data\HGMD_xxxx*
```

The setup script copies:

```text
HGMD folder/final-annotation-df.csv       -> data/maps/final-annotation-df.csv
HGMD folder/mzmine/data_iimn_gnps.mgf     -> data/mgf/data_iimn_gnps.mgf
```

Existing files in `data/maps/` and `data/mgf/` are backed up under `backups/input_files_before_setup_*`.

HGMD analysis runs are saved under `runs/HGMD_xxxx/`. If that folder already exists, the old run folder is moved to `runs/HGMD_xxxx_previous_YYYYMMDD_HHMMSS/` before the fresh bundle is written. Setup-only commands do not save an output bundle.

### User-File Mode

If no `--hgmd` value is supplied, the workflow uses the files already present here:

```text
data/maps/final-annotation-df.csv
data/mgf/data_iimn_gnps.mgf
```

The MetOrigin database should remain here:

```text
data/metorigin/MetOriginDB.csv
```
