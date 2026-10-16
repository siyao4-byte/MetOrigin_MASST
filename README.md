# MetOrigin and domainMASST Recorded Source Associations

Match metabolite annotations against a local MetOriginDB snapshot, retain all recorded source associations and bacterial lineages, and generate circular networks, UpSet plots, association diagrams and an interactive report. Optional domainMASST tools support a separate MS/MS spectral-occurrence workflow.

**These are recorded source associations, not definitive biological origins.** A human-only association means only human evidence is recorded under the matching criteria; it does not establish exclusive human production. Blood, urine and faeces are occurrence sites, not necessarily sites of synthesis.

## Quick Start: Windows PowerShell

Run commands from the repository root. Python 3.10 is recommended for compatibility with the optional domainMASST tools. Calling the virtual-environment executable directly avoids activation and execution-policy issues.

```powershell
git clone https://github.com/siyao4-byte/MetOrigin_MASST.git
Set-Location MetOrigin_MASST
py -3.10 -m venv .venv310
.\.venv310\Scripts\python.exe -m pip install --upgrade pip
.\.venv310\Scripts\python.exe -m pip install -r requirements.txt
```

For plots, install R and set its executable path to the version installed on your computer:

```powershell
$Rscript = 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe'
& $Rscript -e "install.packages(c('tidyverse', 'patchwork', 'yaml'), repos='https://cloud.r-project.org')"
```

The `&` is required in PowerShell when executing a quoted path or a path stored in a variable. You can omit R entirely for matching and table export.

### 1. Provide the Inputs

Download the MetOriginDB CSV through the [TidyMass database page](https://www.tidymass.org/databases/) and retain the unmodified file at:

```text
data/metorigin/MetOriginDB.csv
```

Place your annotation table at `data/maps/Total-List-Of-Annotations.csv`, or specify its actual path in a local configuration. Typical input columns are `compound.name`, `smiles`, `CID`, `InChIKey`, `HMDB_ID` and `KEGG_ID`; not all are required. At least one usable matching key is needed. A `feature.ID` is retained when globally unique; otherwise run-local row IDs preserve individual records.

MGF spectra are **not required** for MetOrigin matching. They are needed only for the optional domainMASST route.

### 2. Configure a Run

```powershell
Copy-Item config\config_total_annotations_metorigin.yaml config\total.local.yaml
notepad config\total.local.yaml
```

Keep local YAML files directly in `config/`, because paths are resolved relative to its parent directory, the project root. Use forward slashes inside YAML paths. Set:

- `maps_annotation_file`: your annotation CSV.
- `metorigin_file`: the original database CSV.
- `metorigin_index`: a reusable SQLite index location for this database snapshot.
- `output_file`: a **new** run path, such as `runs/my-run/outputs/metorigin_matches/MAPS_with_MetOriginDB_origin.csv`.
- Keep `input_pubchem_namespace: cid` and `allow_sid_fallback: false`.

Local `*.local.yaml` configurations are ignored by Git. Do not overwrite an existing local configuration when repeating these steps.

### 3. Match and Plot

```powershell
.\.venv310\Scripts\python.exe scripts\run_indexed_workflow.py --config config\total.local.yaml --plots --rscript "$Rscript"
```

The runner builds a missing database index, checks the snapshot, validates input CIDs through PubChem, matches annotations and exports tables and PNG/PDF plots. Subsequent runs reuse the index and CID validation cache. The database and index are large; allow several tens of GB of free disk space for indexing and complete evidence exports.

Matching only:

```powershell
.\.venv310\Scripts\python.exe scripts\run_indexed_workflow.py --config config\total.local.yaml
```

Add `--html` for the interactive report. Reports can be large for complete total lists. Add `--offline` to reuse cached CID validation without contacting PubChem; unvalidated CIDs are not automatically accepted, and SID fallback remains forbidden. No confidence-level or MSNovelist exclusion is applied automatically: prepare the intended input scope before running.

Existing matching outputs are protected against overwrite. For a new analysis, change `output_file` to a new run directory. To redraw existing results, do **not** rerun matching:

```powershell
$env:METORIGIN_RUN_DIR = 'runs/my-run'
$env:METORIGIN_INCLUDE_DOMAINMASST = 'false'
& $Rscript scripts\Graph_All_compact_workflow.R
.\.venv310\Scripts\python.exe scripts\generate_interactive_origin_report.py --run runs/my-run
```

## Matching and Indexed Evidence

The canonical Python matcher checks CID, SMILES, InChIKey, HMDB, KEGG and exact normalized name. Structural conflicts and unresolved ambiguity are retained in the audit rather than silently accepted. Compatible lower-priority records can contribute additional evidence when they agree with the anchored structure.

- **Input PubChem identifiers are always CID.** CID and SID are separate namespaces; a failed CID never triggers a SID lookup.
- Full structure keys preserve stereochemical distinctions; no connectivity-only relaxation or salt stripping is performed.
- Name-only evidence is explicitly labelled as weaker support.
- The original CSV is unchanged and recorded by SHA-256 and a snapshot manifest.
- `lineage_definitions` stores each complete bacterial lineage once; `bacterial_lineages` links database records to lineage IDs. IDs are snapshot-local, not permanent taxonomy identifiers.
- Organisms and parallel taxonomy fields are expanded together before deduplication. Missing ranks stay missing; incompatible list lengths are flagged rather than guessed.
- Other sources are retained in `source_associations`, including individual details, status and record-level provenance. These details do not yet use shared source-definition IDs.
- Input-to-record links are preserved in `match_evidence`. All associations remain available; top-N and report limits affect display only.

Counts distinguish input rows, normalized compound-name identities with SMILES fallback, and overlapping source associations. SMILES-only entries are not automatically confirmed metabolites. See [method details](docs/indexed-workflow.md).

## Outputs

```text
runs/my-run/
  config/                         Run configuration snapshot
  outputs/metorigin_matches/
    MAPS_with_MetOriginDB_origin.csv
    source_associations.csv       Complete identity-record-source evidence
    bacterial_lineages.csv        Complete synchronized lineage evidence
    match_evidence.csv
    match_evidence.sqlite         Indexed candidate and selected matches
  outputs/qc_reports/             Counts, unmatched rows and audit summaries
  Graph1_MetOrigin_R_outputs/      PNG/PDF plots and complete plot-source tables
  interactive_origin_report.html  Optional interactive report
```

`_FULL` and `_COMPACT` summary CSVs are compatibility aliases. Detailed source and taxonomy evidence is in the normalized tables, not repeated flattened strings. The compatibility `feature_id` in association tables is the identity ID; use the input summary and match evidence to trace individual input rows.

## Optional domainMASST Workflow

domainMASST tests public **spectral occurrence**, which is different from MetOriginDB identity matching. It is not run or submitted by the indexed runner. Keep upstream tools under `external/microbe_masst/` and follow their installation requirements from [microbe_masst](https://github.com/robinschmid/microbe_masst).

The existing `setup_and_run_workflow.py`, `run_domainmasst_unmatched.ps1`, and Graph 2/3 scripts remain for single-dataset workflows. Run `--help` on the Python preparation script before selecting inputs:

```powershell
.\.venv310\Scripts\python.exe scripts\prepare_unmatched_ms2_domainmasst_merge.py --help
.\.venv310\Scripts\python.exe scripts\setup_and_run_workflow.py --help
```

For HGMD input discovery, supply `--hgmd_base 'D:/LCMS-data'` or set `$env:HGMD_BASE`; no institutional network drive is assumed. Use unique dataset-local feature IDs that correspond to the MGF scans. A multi-dataset total list cannot be mapped to one MGF by repeated scan IDs alone. Review query inputs before using `--run_masst` or `-RunMasst`, which submits spectral queries to the external service.

The indexed MetOrigin workflow has regression tests. The optional external domainMASST search was not rerun as part of this upgrade.

## Validation and Troubleshooting

```powershell
.\.venv310\Scripts\python.exe -m unittest discover -s tests -v
.\.venv310\Scripts\python.exe scripts\audit_metorigin_index.py --help
.\.venv310\Scripts\python.exe scripts\compare_indexed_run.py --help
```

- **Python not found:** install Python 3.10 with the Windows launcher, then recreate `.venv310`.
- **Rscript not found:** correct `$Rscript` and pass `--rscript`, or put Rscript on `PATH`.
- **Output already exists:** choose a fresh run path or use the plot-only commands above.
- **Snapshot mismatch:** preserve the previous index and build the changed database into a new index path.
- **Few matches:** check the named versus SMILES-only denominators, identifier availability and conflict audit. No match is not proof of biological absence.

Source code, configuration templates, documentation and tests belong in Git. Local inputs, database indexes, runs, backups, environments and credentials do not. See [folder structure](docs/folder_structure.md).
