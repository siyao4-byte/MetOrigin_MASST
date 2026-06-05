# MAPS MetOrigin and domainMASST Workflow

This project matches MAPS annotations to MetOriginDB, prepares MetOrigin-unmatched MS2 spectra for domainMASST, summarizes domainMASST JSON/TSV evidence, and produces the current Graph1/Graph2/Graph3 plots.

## Setup

Create the Python 3.10 environment:

```powershell
py -3.10 -m venv .venv310
.\.venv310\Scripts\python.exe -m pip install --upgrade pip
.\.venv310\Scripts\python.exe -m pip install -r requirements.txt
```

For full domainMASST batch runs, keep the upstream code here:

```text
external/microbe_masst/
```

The MetOrigin database should be:

```text
data/metorigin/MetOriginDB.csv
```

This database is not included in this GitHub repository because the CSV is very large. Download the MetOriginDB MS1 compound database separately from the TidyMass databases page:

```text
https://www.tidymass.org/databases/
```

Choose the MetOriginDB CSV download, then place/rename it as:

```text
data/metorigin/MetOriginDB.csv
```

## Recommended Front Door

Use the setup workflow:

```powershell
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py
```

This has two input modes.

## Mode 1: HGMD Path

The base path is:

```text
Y:\MA_BPA_Microbiome\LCMS data
```

The script looks for folders whose names start with the requested HGMD ID, for example:

```text
HGMD_1047
HGMD_1047_pro1-434
```

Run:

```powershell
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047
```

It finds and confirms:

```text
HGMD folder/final-annotation-df.csv
HGMD folder/mzmine/data_iimn_gnps.mgf
```

Then it copies them into the standard project inputs:

```text
data/maps/final-annotation-df.csv
data/mgf/data_iimn_gnps.mgf
```

Existing input files are backed up under:

```text
backups/input_files_before_setup_*
```

At the end of an HGMD analysis run, the workflow saves a run bundle:

```text
runs/HGMD_1047/
  data/
  outputs/
  Graph1_MetOrigin_R_outputs/
  Graph2_domainMASST_R_outputs/
  Graph3_integrated_origin_R_outputs/
```

If `runs/HGMD_1047/` already exists, the old folder is moved to:

```text
runs/HGMD_1047_previous_YYYYMMDD_HHMMSS/
```

and a fresh `runs/HGMD_1047/` folder is written.

Setup-only commands do not save an output bundle, because no new analysis outputs were made.

## Mode 2: User Input Files

If no `--hgmd` is supplied, the script uses the existing files in `data/`:

```text
data/maps/final-annotation-df.csv
data/mgf/data_iimn_gnps.mgf
```

Run:

```powershell
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --interactive_run_prompts
```

The script shows the current default inputs and asks whether to continue.

## Running Analysis Steps

You can run setup only, or add flags:

```powershell
# HGMD setup only
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047

# HGMD setup + MetOrigin matching
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047 --run_metorigin

# HGMD setup + MetOrigin + prepare domainMASST input only
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047 --run_metorigin --prepare_domainmasst

# HGMD setup + MetOrigin + full domainMASST batch
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047 --run_metorigin --run_domainmasst

# HGMD setup + MetOrigin + full domainMASST batch + all plots, saved under runs/HGMD_1047/
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047 --run_metorigin --run_domainmasst --run_graphs

# Run graphs after existing outputs are available
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --run_graphs
```

For a fully prompted run:

```powershell
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047 --interactive_run_prompts
```

Use `--yes` only when you already know the detected HGMD folder and files are correct.

Run folder options:

```powershell
# Use a custom runs/ folder name
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047 --run_label HGMD_1047_test1 --run_metorigin --run_graphs

# Save default user-file mode outputs under runs/<label>/
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --run_label USER_TEST --save_run_outputs --run_graphs

# Disable run-folder saving
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --hgmd HGMD_1047 --no_save_run
```

## Graph Workflow

The compact graph workflow runs Graph1, Graph2, and Graph3:

```powershell
"C:\Program Files\R\R-4.5.1\bin\Rscript.exe" scripts\Graph_All_compact_workflow.R
```

Current plot outputs are PNG-only.

## Main Outputs

```text
outputs/metorigin_matches/MAPS_with_MetOriginDB_origin.csv
outputs/domainmasst_unmatched_ms2/
Graph1_MetOrigin_R_outputs/
Graph2_domainMASST_R_outputs/
Graph3_integrated_origin_R_outputs/
```

## Manual Fallback Commands

```powershell
.\.venv310\Scripts\python.exe scripts\metorigin_match.py --config config\config.yaml
powershell -ExecutionPolicy Bypass -File scripts\run_domainmasst_unmatched.ps1
powershell -ExecutionPolicy Bypass -File scripts\run_domainmasst_unmatched.ps1 -RunMasst
"C:\Program Files\R\R-4.5.1\bin\Rscript.exe" scripts\Graph_All_compact_workflow.R
```
