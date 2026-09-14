# Step 2d: fetch + audit group info for promising external-validation candidates

$geo = Join-Path $env:USERPROFILE ".codex\skills\geo-search\scripts\geo.ps1"
$accs = @('GSE99325','GSE47183','GSE175759','GSE145747','GSE189005','GSE294519','GSE316126')
$logDir = "codex_output\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

foreach ($acc in $accs) {
    $outCsv = "codex_output\geo\$acc\$acc-samples.csv"
    $log = Join-Path $logDir "$acc-samples.log"
    Start-Process -FilePath "powershell" -NoNewWindow -Wait `
        -ArgumentList @("-NoProfile","-ExecutionPolicy","Bypass","-File",$geo,"-Action","samples","-Accession",$acc,"-Out",$outCsv) `
        -RedirectStandardOutput $log | Out-Null

    if (-not (Test-Path -LiteralPath $outCsv)) { Write-Output ("$acc | samples fetch failed"); continue }
    $rows = Import-Csv $outCsv
    Write-Output ("=== $acc  (n=" + $rows.Count + ")")
    Write-Output ("   source_name:")
    $rows | Group-Object source | Sort-Object Count -Descending | Select-Object -First 6 |
        ForEach-Object { Write-Output ("      " + $_.Count.ToString().PadLeft(4) + "  " + $_.Name) }
    $kv = @{}
    foreach ($r in $rows) {
        foreach ($pair in ($r.characteristics -split '\s*\|\s*')) {
            if ($pair -match '^\s*([^:]+):\s*(.+?)\s*$') {
                $k = $Matches[1].Trim().ToLower(); $v = $Matches[2].Trim()
                if (-not $kv.ContainsKey($k)) { $kv[$k] = @{} }
                if (-not $kv[$k].ContainsKey($v)) { $kv[$k][$v] = 0 }
                $kv[$k][$v]++
            }
        }
    }
    foreach ($k in ($kv.Keys | Where-Object { $_ -match 'disease|diagnosis|group|condition|state|type|stage' })) {
        $vals = ($kv[$k].GetEnumerator() | Sort-Object Value -Descending | Select-Object -First 8 | ForEach-Object { "$($_.Key)($($_.Value))" }) -join ' | '
        Write-Output ("   $k : $vals")
    }
    Write-Output ""
}
