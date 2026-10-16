# Folder Structure

```text
MetOrigin_MASST/
  config/                 Tracked templates; ignored *.local.yaml overrides
  scripts/                Indexed matcher, reporting and optional domainMASST tools
  tests/                  Synthetic regression tests (no study data)
  docs/                   Usage and method documentation
  data/                   Local inputs and database snapshots; ignored
    maps/
    metorigin/
      MetOriginDB.csv
      indexed-*/          SQLite index and snapshot manifest
    mgf/                  Optional single-dataset spectra
  runs/                   Per-run configuration, evidence, audits and plots; ignored
  backups/                Local recovery copies; ignored
  external/               Optional upstream domainMASST installation; ignored
  .venv310/               Python environment; ignored
```

Run PowerShell commands from the repository root. The recommended entry point is `scripts/run_indexed_workflow.py`; see the [README](../README.md) for setup, input configuration and plot-only reruns.

The legacy single-dataset setup helper can also use root-level `outputs/` and Graph 1/2/3 output folders, all ignored by Git. It can discover HGMD folders using `--hgmd_base` or `HGMD_BASE`. Original input files are backed up before replacement. These compatibility helpers are retained because other scripts still invoke them.

Do not delete local `data/`, `runs/` or `backups/` as part of repository cleanup. Generated files are excluded from version control, not automatically removed from disk.
