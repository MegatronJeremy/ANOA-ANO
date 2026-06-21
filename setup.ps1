# setup.ps1 -- one-shot environment setup for the nanoplastic scRNA-seq pipeline (Windows).
# Finds Python 3.10+, creates .venv, installs requirements, and checks the raw data.
# Idempotent: safe to re-run.
$ErrorActionPreference = "Stop"
$root = $PSScriptRoot

Write-Host "=== Setup: nanoplastic scRNA-seq pipeline ===" -ForegroundColor Cyan

# --- 1. Find a Python 3.10+ interpreter -----------------------------------
function Find-Python {
    foreach ($cmd in @("py -3.12", "py -3.11", "py -3.10", "py -3", "python")) {
        $parts = $cmd.Split(" ")
        $exe = $parts[0]
        if (Get-Command $exe -ErrorAction SilentlyContinue) {
            try {
                $v = & $exe $parts[1..($parts.Length-1)] -c "import sys; print('%d.%d' % sys.version_info[:2])" 2>$null
                if ($v -match "^3\.(1[0-9]|[2-9][0-9])$") { return ,@($exe, $parts[1..($parts.Length-1)]) }
            } catch {}
        }
    }
    return $null
}

$py = Find-Python
if ($null -eq $py) {
    Write-Host "No Python 3.10+ found. Install from https://www.python.org/downloads/ and re-run." -ForegroundColor Red
    exit 1
}
Write-Host "Using Python: $($py[0]) $($py[1])" -ForegroundColor Green

# --- 2. Create venv -------------------------------------------------------
$venv = Join-Path $root ".venv"
$vpy  = Join-Path $venv "Scripts\python.exe"
if (-not (Test-Path $vpy)) {
    Write-Host "Creating .venv ..." -ForegroundColor Cyan
    & $py[0] $py[1] -m venv $venv
} else {
    Write-Host ".venv already exists -- reusing." -ForegroundColor DarkGray
}

# --- 3. Install dependencies ---------------------------------------------
Write-Host "Installing requirements.txt ..." -ForegroundColor Cyan
& $vpy -m pip install --upgrade pip | Out-Null
& $vpy -m pip install -r (Join-Path $root "requirements.txt")

# --- 4. Check raw data ----------------------------------------------------
$raw = Join-Path $root "data\raw"
$files = @(
    "filtered_feature_bc_matrix.h5ad",
    "filtered_feature_bc_matrix_Sample2.h5ad",
    "filtered_feature_bc_matrix_Sample3.h5ad",
    "filtered_feature_bc_matrix_Sample4.h5ad"
)
$missing = @()
foreach ($f in $files) { if (-not (Test-Path (Join-Path $raw $f))) { $missing += $f } }
if ($missing.Count -gt 0) {
    Write-Host "Missing raw data in data\raw\:" -ForegroundColor Yellow
    $missing | ForEach-Object { Write-Host "  - $_" -ForegroundColor Yellow }
    Write-Host "Download from Zenodo DOI 10.5281/zenodo.15866724 (see legacy\src\download_data.py)." -ForegroundColor Yellow
} else {
    Write-Host "All 4 raw .h5ad files present." -ForegroundColor Green
}

Write-Host ""
Write-Host "Setup done. Verify with:  .\run.ps1 check" -ForegroundColor Green
Write-Host "Then run a stage, e.g.:    .\run.ps1 qc -Smoke -Debug" -ForegroundColor Green
