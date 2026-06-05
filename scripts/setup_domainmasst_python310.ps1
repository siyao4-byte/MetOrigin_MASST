$ErrorActionPreference = "Stop"
$ProjectDir = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $ProjectDir

if (-not (Get-Command py -ErrorAction SilentlyContinue)) {
    throw "The Python launcher 'py' is not available. Install Python 3.10 first."
}

py -3.10 --version
py -3.10 -m venv .venv310
.\.venv310\Scripts\python.exe -m pip install --upgrade pip
.\.venv310\Scripts\python.exe -m pip install -r requirements.txt

if (-not (Test-Path "external\microbe_masst")) {
    New-Item -ItemType Directory -Force external | Out-Null
    git clone https://github.com/robinschmid/microbe_masst.git external\microbe_masst
}

.\.venv310\Scripts\python.exe -c "import pandas, yaml, pyteomics, tqdm, requests; print('Python/domainMASST basics OK')"
