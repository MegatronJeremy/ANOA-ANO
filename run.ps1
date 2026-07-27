# run.ps1 -- thin wrapper around run_pipeline.py (NOT a second code path).
# Forwards everything to the headless driver; new stages appear automatically
# once they are registered in STAGE_REGISTRY.
#
#   .\run.ps1                      # interactive menu (run_pipeline.py with no args)
#   .\run.ps1 check                # environment doctor (--check)
#   .\run.ps1 qc                   # run a stage
#   .\run.ps1 qc -Smoke -Debug     # ... on the smoke subsample, verbose
#   .\run.ps1 all                  # run every registered stage in order
#   .\run.ps1 data                 # download the raw dataset from Zenodo into data/raw/
#   .\run.ps1 test                 # pytest
#   .\run.ps1 menu                 # force the interactive menu
#
# run.bat / setup.bat are thin cmd.exe shims around these scripts (double-click
# or plain cmd), forwarding all arguments. If .venv is missing, run creates it
# via setup automatically on first use.
param(
    [string]$Command = "menu",
    [switch]$Smoke,
    [switch]$Debug
)
$ErrorActionPreference = "Stop"
$root = $PSScriptRoot

# Corporate SSL-inspecting proxies re-sign HTTPS with an internal root CA that a
# fresh venv's pip -- and tools that download over HTTPS, like celltypist's model
# fetch -- don't trust, producing CERTIFICATE_VERIFY_FAILED. If a local CA bundle
# has been exported (corp_ca_bundle.pem: machine-specific, git-ignored), point
# pip and requests-based downloads at it. Set before the venv check so setup.ps1
# inherits it too. No-op for anyone without the file (e.g. a grader on a normal
# network) -- the default system cert handling is used unchanged.
$caBundle = Join-Path $root "corp_ca_bundle.pem"
if (Test-Path $caBundle) {
    $env:PIP_CERT           = $caBundle
    $env:REQUESTS_CA_BUNDLE = $caBundle
    $env:SSL_CERT_FILE      = $caBundle
    Write-Host "Using local corporate CA bundle for pip/requests SSL." -ForegroundColor DarkGray
}

# Prefer the venv interpreter. It must (a) exist and (b) actually have the
# dependencies installed -- a half-built .venv (created but `pip install` never
# finished) would otherwise fall through to a confusing ModuleNotFoundError at
# runtime. So probe for a core dependency (numpy, the first thing the driver
# imports) and run setup when the venv is missing OR incomplete.
$py = Join-Path $root ".venv\Scripts\python.exe"

function Test-VenvReady {
    if (-not (Test-Path $py)) { return $false }
    # Probe for a core dependency. Local SilentlyContinue so the broken venv's
    # traceback on stderr doesn't get escalated into a terminating error by the
    # script-level $ErrorActionPreference = "Stop"; we judge purely by exit code.
    $ErrorActionPreference = "SilentlyContinue"
    & $py -c "import numpy" *> $null
    return ($LASTEXITCODE -eq 0)
}

if (-not (Test-VenvReady)) {
    if (Test-Path $py) {
        Write-Host ".venv exists but dependencies are missing -- running setup to finish it..." -ForegroundColor Yellow
    } else {
        Write-Host ".venv not found -- running setup first (one-time)..." -ForegroundColor Yellow
    }
    & (Join-Path $root "setup.ps1")
    if (-not (Test-VenvReady)) {
        Write-Host "Setup did not produce a working .venv (numpy still not importable) -- see the output above." -ForegroundColor Red
        exit 1
    }
}

# Ordered list of registered stages. Keep in sync with STAGE_REGISTRY as
# stages 3-6 are added (qc, integration, annotation, composition, de, size).
$Stages = @("qc", "integration", "annotation", "composition", "de", "size", "bonus")

$flags = @()
if ($Smoke) { $flags += "--smoke-test" }
if ($Debug) { $flags += "--debug" }

function Invoke-Stage($stage) {
    Write-Host ">> stage: $stage $flags" -ForegroundColor Cyan
    & $py (Join-Path $root "run_pipeline.py") --stage $stage @flags
}

switch ($Command.ToLower()) {
    "menu"  { & $py (Join-Path $root "run_pipeline.py") }
    "check" { & $py (Join-Path $root "run_pipeline.py") --check }
    "test"   { & $py -m pytest $root }
    "data"   { & $py -m src.download_data }
    "slides" { & $py -m src.make_slides }
    "all"    { foreach ($s in $Stages) { Invoke-Stage $s } }
    default {
        # Forward any other token straight to --stage. run_pipeline.py's argparse
        # is the source of truth for valid stage names, so newly registered stages
        # work here with no edit to this script (only `all` needs $Stages kept in sync).
        Invoke-Stage $Command.ToLower()
    }
}
