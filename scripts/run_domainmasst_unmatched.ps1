param(
    [switch]$RunMasst,
    [int]$ParallelQueries = 5,
    [double]$MinCos = 0.7,
    [double]$MzTol = 0.02,
    [double]$PrecursorMzTol = 0.02,
    [int]$MinMatchedSignals = 3
)

$ErrorActionPreference = "Stop"
$ProjectDir = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$Python = Join-Path $ProjectDir ".venv310\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python 3.10 virtual environment not found at $Python. Create it with: py -3.10 -m venv .venv310"
}

Set-Location $ProjectDir

$argsList = @(
    "scripts\prepare_unmatched_ms2_domainmasst_merge.py",
    "--metorigin_csv", "outputs\metorigin_matches\MAPS_with_MetOriginDB_origin.csv",
    "--mgf", "data\mgf\data_iimn_gnps.mgf",
    "--out_dir", "outputs\domainmasst_unmatched_ms2",
    "--min_cos", "$MinCos",
    "--mz_tol", "$MzTol",
    "--precursor_mz_tol", "$PrecursorMzTol",
    "--min_matched_signals", "$MinMatchedSignals",
    "--parallel_queries", "$ParallelQueries"
)

if ($RunMasst) {
    $argsList += "--run_masst"
}

& $Python @argsList
