param([ValidateRange(1,1000)][int]$Stages = 1)
$ErrorActionPreference = 'Stop'
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$pythonPath = Join-Path $projectRoot '.envs/pysph/python.exe'
$gpuRoot = Join-Path $projectRoot 'outputs/gpu_acceleration_20261007'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:PYTHONPATH = ''
Push-Location $projectRoot
try {
    for ($stageIndex = 0; $stageIndex -lt $Stages; $stageIndex++) {
        $candidates = @(Get-ChildItem -LiteralPath $gpuRoot -Directory | ForEach-Object {
            $manifestPath = Join-Path $_.FullName 'manifest.json'
            if (Test-Path -LiteralPath $manifestPath) {
                $meta = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
                if ($meta.status -eq 'PASS_GPU_SCREENED_DIAGNOSTIC_CONTINUATION') {
                    [PSCustomObject]@{ Name = $_.Name; Time = [double]$meta.screened_end_time_s }
                }
            }
        })
        $last = $candidates | Sort-Object -Property Time -Descending | Select-Object -First 1
        $newName = 'continue_' + (Get-Date -Format 'yyyyMMdd_HHmmss_fff')
        $arguments = @('-B', 'scripts/run_afsd_gpu_20261007.py', '--mode', 'continue',
            '--output', ('outputs/gpu_acceleration_20261007/' + $newName),
            '--cost', 'outputs/gpu_acceleration_20261007/benchmark_v1',
            '--validation', 'outputs/gpu_acceleration_20261007/replay_v1')
        if ($null -ne $last) {
            $arguments += @('--parent', ('outputs/gpu_acceleration_20261007/' + $last.Name))
        }
        & $pythonPath @arguments
        if ($LASTEXITCODE -ne 0) {
            throw 'GPU continuation stopped: inspect the new output receipt; no failed state will be continued.'
        }
    }
}
finally {
    Pop-Location
}
