# MetOrigin and domainMASST Workflow (a combined workflow for origin inference with interactive visualization ouput)

This project matches annotated metabolites to MetOriginDB, prepares MetOrigin-unmatched MS2 spectra for domainMASST, summarizes domainMASST JSON/TSV evidence, and produces the current Graph1/Graph2/Graph3 plots. The main resources used are MetOriginDB https://github.com/tidymass/metorigindb , and domainMASST https://github.com/robinschmid/microbe_masst

It will work on both annotated features, or features with ms2 spectra (with provided .mgf file)

## Setup

Clone the repository and create the Python 3.10 environment:

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

This database is not included in this GitHub repository because the CSV is very large (at least 5 GB). Download the MetOriginDB MS1 compound database separately from the TidyMass databases page:

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

## Matching Logic and Origin Interpretation

This workflow combines two different kinds of origin evidence:

```text
MetOrigin  = compound identity -> curated compound-origin database
domainMASST = MS/MS spectrum -> public dataset occurrence by source/domain
```

### Workflow Flow Chart

```text
Input MAPS annotation table
        |
        v
Clean column names and detect key fields
        |
        v
Use available identifiers:
CID, SMILES, InChIKey, HMDB, KEGG, compound name
        |
        v
Try MetOriginDB matching in priority order
        |
        +--> 1. PubChem CID exact match
        |
        +--> 2. SMILES exact match
        |
        +--> 3. InChIKey exact match
        |
        +--> 4. HMDB ID exact match
        |
        +--> 5. KEGG ID exact match
        |
        +--> 6. Compound name exact match
        |
        v
If matched:
aggregate all MetOriginDB origin records for that compound
        |
        v
Export MetOrigin origin evidence:
food, bacteria, plant, animal, human, drug, environment
        |
        v
For unmatched features with MS2:
prepare MGF queries for domainMASST
        |
        v
Run domainMASST spectral search
        |
        v
Summarize public-dataset spectral occurrence by domain:
microbe, plant, animal, human, food/environment, etc.
        |
        v
Merge MetOrigin + domainMASST evidence
        |
        v
Generate static plots and interactive HTML report
```

### MetOrigin Matching

MetOrigin asks:

```text
Does this annotated compound identity exist in MetOriginDB?
If yes, what curated origins are associated with it?
```

Matching is sequential. Once a feature matches at one step, it is removed from the unmatched pool and is not counted again at later steps. For example, a PubChem CID match will not also be counted as a SMILES or name match.

MetOriginDB can contain multiple records for the same compound. The workflow aggregates all matching MetOriginDB records before joining them back to MAPS, so one MAPS feature can be linked to multiple origins and many source details.

### domainMASST Evidence

domainMASST asks:

```text
Does this MS/MS spectrum occur in public datasets assigned to certain source domains?
```

The domainMASST `occurrence fraction` is a relative occurrence measure:

```text
matched spectra or records / total spectra or records available in that source group
```

It should be interpreted as public-dataset spectral prevalence, not as a direct MetOrigin-style confidence score.

### Limitations

Origin annotations from this workflow are evidence-supported source associations, not definitive causal origins.

MetOrigin limitations:

- It requires a matching identifier or exact compound name.
- Exact SMILES/name matching can miss synonyms, salts, stereoisomers, tautomers, and different canonical SMILES.
- MetOriginDB is not complete; no hit does not mean no biological origin.
- One compound can occur in many domains, so Food, Bacteria, Plant, Animal, Human, Drug, or Environment labels are not exclusive.
- Bacterial taxonomy is structured, but food/plant/animal fields are often long text lists rather than full taxonomic trees.
- Aggregating MetOriginDB rows can collapse many species/genus records into long lists.

domainMASST limitations:

- Spectral occurrence is not the same as chemical identity.
- Public datasets are biased toward well-studied organisms and sample types.
- A spectrum found in a domain may reflect contamination, shared chemistry, diet, host-microbe exchange, or annotation ambiguity.
- Occurrence fraction depends on database composition, not only biology.
- One feature can appear in multiple domains, so domainMASST evidence is not mutually exclusive.
- MASST does not always distinguish structural isomers or close analogues.

Integrated interpretation:

- MetOrigin confidence and domainMASST evidence are not directly equivalent.
- Agreement between MetOrigin and domainMASST is stronger evidence, but disagreement is not necessarily an error.
- Results should be reported as possible or supported source associations.

Suggested wording:

```text
The workflow reports origin evidence from curated compound-origin annotations
and public spectral occurrence. Results indicate supported source associations
rather than exclusive or causal origins.
```

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
