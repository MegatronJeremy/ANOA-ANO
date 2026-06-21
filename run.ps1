# run.ps1 -- thin wrapper around run_pipeline.py (NOT a second code path).
# Forwards everything to the headless driver; new stages appear automatically
# once they are registered in STAGE_REGISTRY.
#
#   .\run.ps1                      # interactive menu (run_pipeline.py with no args)
#   .\run.ps1 check                # environment doctor (--check)
#   .\run.ps1 qc                   # run a stage
#   .\run.ps1 qc -Smoke -Debug     # ... on the smoke subsample, verbose
#   .\run.ps1 all                  # run every registered stage in order
#   .\run.ps1 test                 # pytest
#   .\run.ps1 menu                 # force the interactive menu
param(
    [string]$Command = "menu",
    [switch]$Smoke,
    [switch]$Debug
)
$ErrorActionPreference = "Stop"
$root = $PSScriptRoot

# Prefer the venv interpreter; fall back to whatever `python` resolves to.
$py = Join-Path $root ".venv\Scripts\python.exe"
if (-not (Test-Path $py)) {
    Write-Host ".venv not found -- run .\setup.ps1 first (falling back to system python)." -ForegroundColor Yellow
    $py = "python"
}

# Ordered list of registered stages. Keep in sync with STAGE_REGISTRY as
# stages 3-6 are added (qc, integration, annotation, composition, de, size).
$Stages = @("qc", "integration")

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
    "test"  { & $py -m pytest $root }
    "all"   { foreach ($s in $Stages) { Invoke-Stage $s } }
    default {
        # Forward any other token straight to --stage. run_pipeline.py's argparse
        # is the source of truth for valid stage names, so newly registered stages
        # work here with no edit to this script (only `all` needs $Stages kept in sync).
        Invoke-Stage $Command.ToLower()
    }
}
