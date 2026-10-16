# Indexed recorded-source workflow

Historical runs and the downloaded CSV remain unchanged. See the [README](../README.md) for the recommended PowerShell setup and run commands.

## Identifier contract

Every PubChem number in the MAPS total list is a **CID**. CID failure NEVER invokes a SID endpoint, never queries the SID index, and never reinterprets the input as a SID. Ordinary structure/HMDB/KEGG/name matching remains available and is labelled separately.

Database CID and SID fields remain separate. Multiple values in the database's CID field are stored as unverified CID candidates, not guessed CID/SID pairs. PubChem compound-property lookup and exact structure agreement determine whether a CID candidate can support a match. Unavailable validation is recorded, not silently treated as success. Only public CID numbers are submitted to PubChem, not spectra or participant metadata.

Full InChIKeys derived from RDKit preserve stereochemical distinctions. No salt stripping, stereochemistry removal or connectivity-only relaxation is performed. Name-only matches remain explicitly weaker evidence. Conflicting structural candidates are logged; ambiguous multiple-structure matches without a resolving input structure are not assigned.

## Database snapshot

`data/metorigin/MetOriginDB.csv` is never modified. The indexed directory contains SQLite and a snapshot manifest with SHA-256, byte count, original path, modification time, schema and builder-code checksum. The index refuses to overwrite an existing path. Source checksum is checked before matching. The upstream commit is recorded as unknown when it is absent from the original download, not invented.

Tables:

- `database_records`: stable snapshot-local row ID and lossless compressed original record.
- `identifiers`: typed identifier, original field and explicit validation status; indexed by namespace/value.
- `source_associations`: all recorded source states and details, including Unknown and No.
- `lineage_definitions` and `bacterial_lineages`: normalized complete lineage tuples and their record associations.
- `validation_issues`: malformed identifiers and taxonomy-length mismatches.
- `pubchem_validation`: cached CID-only validation results and retrieval time.

Bacterial names and ranks are split together in their original order, retaining repeated ranks and missing positions. Known unequal-length lists are quarantined, not zipped or broadcast. A wholly missing rank column remains missing for every organism. Deduplication happens after complete tuples have been assembled.

## Execution

```powershell
.\.venv310\Scripts\python.exe scripts\metorigin_index.py --source data\metorigin\MetOriginDB.csv --index data\metorigin\indexed-20261016\metorigin.sqlite
.\.venv310\Scripts\python.exe scripts\metorigin_match.py --config config\total.local.yaml
$env:METORIGIN_RUN_DIR = "runs/my-run"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" scripts\Graph_All_compact_workflow.R
.\.venv310\Scripts\python.exe scripts\generate_interactive_origin_report.py --help
```

Create `config/total.local.yaml` from the tracked template and set its output to `runs/my-run/outputs/metorigin_matches/MAPS_with_MetOriginDB_origin.csv` before using the commands above. Alternatively, `scripts/run_indexed_workflow.py --config config/total.local.yaml --plots --html --rscript PATH_TO_RSCRIPT` builds a missing index, matches, and renders the new run. Replace `PATH_TO_RSCRIPT` with a quoted executable path; see the README. The HTML report keeps legacy runs unchanged when using `--all-runs`.

Skip the build command when that snapshot already exists. Use a new output run directory for subsequent matches. `--offline` reuses cached CID checks; it does not relax validation. The tracked total-list configuration uses an example relative input path; local overrides are ignored by Git. No confidence-level or MSNovelist filtering is introduced by this infrastructure change; all input rows are preserved. Study-specific inclusion filters remain separate.

## Matching and outputs

The Python indexed matcher is authoritative; the R matching entry point delegates to it. Evidence from weaker keys is added only when it agrees with the anchored exact structure. General flags are summarized with positive evidence taking precedence, but the complete table preserves conflicting and unknown records.

Each run exports a feature/annotation summary, complete source associations, complete bacterial lineages, all candidate-match evidence (CSV and indexed SQLite), and a machine-readable QC summary. No analytical table is truncated. Association and lineage exports are normalized by identity and database record, rather than repeating thousands of taxa for every assay row. Their feature_id compatibility column contains the identity ID; input_row_id is intentionally blank. Every input-to-record link remains in match_evidence, enabling reconstruction without duplication. Original feature IDs in the input summary are retained when present and unique; otherwise run-local row IDs are used and the original ID remains available. Compound identities use normalized name with SMILES fallback, matching the project's count convention; these IDs are not claims of chemical equivalence across unresolved isomers.

HTML limits and top-N plot grouping are display policies only. HTML links to full downloadable evidence tables. Graph 3 reads ungrouped source details rather than the grouped plot-label columns. DomainMASST is opt-in and remains distinct spectral-occurrence evidence; no domainMASST queries are submitted by the indexed rebuild.

## Interpretation and counting

Use **recorded source associations**, not definitive biological origins. Human-only recorded association means human is the only positive source category under the matching criteria; it does not demonstrate exclusive human production. Blood, urine and faeces are recorded occurrence sites, not necessarily biosynthetic sites. Unknown is not evidence of absence.

Keep counts of input rows, unique identities, identity-source links and identity-specific-source links separate. A metabolite associated with both human and bacteria contributes once to their UpSet intersection and once to each source set. General-to-specific association widths need not conserve unique-identity totals; figures show overlapping memberships, not material flow.

The QC summary separates named identities, SMILES-only entries and unidentified rows. SMILES-only entries are searchable structures, not automatically confirmed metabolites. Database matches without positive source categories remain in the results but do not enter the positive-source plots. Historical comparisons must use shared identity keys and account for changes in input scope.

## Validation

```powershell
.\.venv310\Scripts\python.exe -m unittest discover -s tests -v
```

Tests cover CID/SID separation, rejected CID/structure conflicts, repeated phyla, missing placeholders, malformed list lengths, multiple sites, compatible lower-priority evidence, snapshot provenance, complete exports and overwrite protection. Historical-run comparisons must distinguish changed inputs from changed matching; compare only shared identities and report the input checksums.

Bar/matrix x limits are explicitly identical and long association labels have line-aware spacing. Large interactive reports can load slowly; PNG/PDF plots and normalized CSV/SQLite tables are available separately. Review rendered figures before publication; automatic PDF collision checks can flag noncontiguous Cairo text spans and intentional network edges behind labels.
