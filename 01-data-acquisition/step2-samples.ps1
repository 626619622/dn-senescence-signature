# Step 2a: fetch sample/group metadata and supplementary file inventory for the 9 chosen datasets

$ErrorActionPreference = 'Continue'
$geo = Join-Path $env:USERPROFILE ".codex\skills\geo-search\scripts\geo.ps1"
if (-not (Test-Path -LiteralPath $geo)) { throw ("geo skill script not found: " + $geo) }
$accs = @('GSE30528','GSE30529','GSE30122','GSE96804','GSE104954','GSE131882','GSE195460','GSE142025','GSE183276')

$logDir = "codex_output\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

foreach ($acc in $accs) {
    $outCsv = "codex_output\geo\$acc\$acc-samples.csv"
    $log = Join-Path $logDir "$acc-samples.log"
    Start-Process -FilePath "powershell" -NoNewWindow -Wait `
        -ArgumentList @("-NoProfile","-ExecutionPolicy","Bypass","-File",$geo,"-Action","samples","-Accession",$acc,"-Out",$outCsv) `
        -RedirectStandardOutput $log | Out-Null
    $tail = (Get-Content $log -Tail 3) -join ' ; '
    Write-Output ("$acc | " + $tail)
}

Write-Output ""
Write-Output "=== supplementary file inventory ==="
foreach ($acc in $accs) {
    $log = Join-Path $logDir "$acc-suppl.log"
    Start-Process -FilePath "powershell" -NoNewWindow -Wait `
        -ArgumentList @("-NoProfile","-ExecutionPolicy","Bypass","-File",$geo,"-Action","suppl","-Accession",$acc) `
        -RedirectStandardOutput $log | Out-Null
    $lines = Get-Content $log | Where-Object { $_ -match '^\s{2}\S' }
    Write-Output ("--- $acc : " + ($lines.Count) + " files")
    $lines | Select-Object -First 8 | ForEach-Object { Write-Output ("     " + $_.Trim()) }
    if ($lines.Count -gt 8) { Write-Output ("     ... (+" + ($lines.Count - 8) + " more)") }
}
