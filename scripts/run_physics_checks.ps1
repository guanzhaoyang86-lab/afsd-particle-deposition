param([string]$RunId = (Get-Date -Format 'yyyyMMdd_HHmmss'))
$ErrorActionPreference = 'Stop'
$project = Split-Path $PSScriptRoot -Parent
Set-Location -LiteralPath $project
$pythonExe = Join-Path $project '.envs\pysph\python.exe'
$logDir = Join-Path $project "logs\physics\$RunId"
$outputDir = Join-Path $project "outputs\debug\$RunId"
if ((Test-Path -LiteralPath $logDir) -or (Test-Path -LiteralPath $outputDir)) { throw 'RunId already exists; no overwrite.' }
New-Item -ItemType Directory -Path $logDir, $outputDir | Out-Null
$results = [System.Collections.Generic.List[object]]::new()
function Invoke-Check([string]$Name, [string[]]$Arguments, [int]$Expected = 0) {
    $ErrorActionPreference = 'Continue'
    & $pythonExe @Arguments 2>&1 | Tee-Object -FilePath (Join-Path $logDir "$Name.log") | Out-Host
    $resultCode = $LASTEXITCODE
    $results.Add([pscustomobject]@{name=$Name; executable=$pythonExe; argv=$Arguments; expected_exit=$Expected; actual_exit=$resultCode})
    $results | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $logDir 'commands_and_status.json') -Encoding utf8
    if ($resultCode -ne $Expected) { throw "$Name exit $resultCode; expected $Expected" }
}
Invoke-Check 'scaffold' @('-B', 'tests/run_tests.py')
Invoke-Check 'physics_gates' @('-B', 'scripts/run_physics_tests.py', '--gate', 'all')
foreach ($case in @('starved', 'copper-core', 'overfed')) {
    Invoke-Check $case @('-B', 'scripts/run_afsd_debug.py', '--case', "configs/debug/$case.json", '--output', (Join-Path $outputDir $case), '--steps', '200')
}
Invoke-Check 'dt_half' @('-B', 'scripts/run_afsd_debug.py', '--case', 'configs/debug/starved_dt_half.json', '--output', (Join-Path $outputDir 'dt_half'), '--steps', '400')
Invoke-Check 'restart' @('-B', 'scripts/run_afsd_debug.py', '--case', 'configs/debug/starved.json', '--output', (Join-Path $outputDir 'restart'), '--restart', (Join-Path $outputDir 'starved/checkpoint.npz'), '--steps', '20')
Invoke-Check 'wall_cap' @('-B', 'scripts/run_afsd_debug.py', '--case', 'configs/debug/starved.json', '--output', (Join-Path $outputDir 'wall_cap'), '--steps', '200', '--wall-seconds', '0.05') 2
# Explicit failure-injection artifact, not a calibrated physics configuration.
$failureConfig = Get-Content -Raw -LiteralPath 'configs/debug/starved.json' | ConvertFrom-Json
$failureConfig.extends = Join-Path $project 'configs/base.json'
$failureConfig.parameters.debug.density_limit.value = 1e-12
$failureConfig.parameters.debug.density_limit.source = 'Explicit diagnostic failure injection; not physical calibration'
$failurePath = Join-Path $logDir 'density_failure_config.json'
$failureConfig | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $failurePath -Encoding utf8NoBOM
Invoke-Check 'density_alarm' @('-B', 'scripts/run_afsd_debug.py', '--case', $failurePath, '--output', (Join-Path $outputDir 'density_alarm'), '--steps', '10') 1
$folders = @('starved', 'copper-core', 'overfed', 'dt_half', 'restart', 'wall_cap', 'density_alarm') | ForEach-Object { Join-Path $outputDir $_ }
Invoke-Check 'readback' (@('-B', 'scripts/validate_physics_outputs.py') + $folders)
Write-Output "CHECKS_PASS logs=$logDir outputs=$outputDir"
