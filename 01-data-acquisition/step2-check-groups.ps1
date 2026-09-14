# Step 2b: group-information audit (ASCII-only script; Chinese report assembled separately)

$accs = @('GSE30528','GSE30529','GSE30122','GSE96804','GSE104954','GSE131882','GSE195460','GSE142025','GSE183276')
$outDir = "codex_output\geo"
$groupKeys = 'disease|diagnosis|group|condition|state|type|subject|stage|class'

$result = @()
foreach ($acc in $accs) {
    $csv = Join-Path $outDir "$acc\$acc-samples.csv"
    if (-not (Test-Path -LiteralPath $csv)) {
        $result += [pscustomobject]@{ accession=$acc; n_samples='MISSING'; platform=''; group_key=''; group_values=''; other_fields='' }
        continue
    }
    $rows = Import-Csv $csv
    $n = $rows.Count
    $platform = ($rows | Group-Object platform | Sort-Object Count -Descending | Select-Object -First 1).Name

    $kv = @{}
    foreach ($r in $rows) {
        $src = $r.source
        if ($src) {
            if (-not $kv.ContainsKey('source_name')) { $kv['source_name'] = @{} }
            if (-not $kv['source_name'].ContainsKey($src)) { $kv['source_name'][$src] = 0 }
            $kv['source_name'][$src]++
        }
        foreach ($pair in ($r.characteristics -split '\s*\|\s*')) {
            if ($pair -match '^\s*([^:]+):\s*(.+?)\s*$') {
                $k = $Matches[1].Trim().ToLower()
                $v = $Matches[2].Trim()
                if (-not $kv.ContainsKey($k)) { $kv[$k] = @{} }
                if (-not $kv[$k].ContainsKey($v)) { $kv[$k][$v] = 0 }
                $kv[$k][$v]++
            }
        }
    }

    $candidate = $kv.Keys | Where-Object { $_ -match $groupKeys } |
        Sort-Object @{Expression = { $kv[$_].Count }; Descending = $true }
    $gKey = $candidate | Select-Object -First 1

    $gText = ''
    if ($gKey) {
        $gText = (($kv[$gKey].GetEnumerator() | Sort-Object Value -Descending | ForEach-Object { "$($_.Key)($($_.Value))" }) -join ' | ')
    }
    $othersText = (($kv.Keys | Where-Object { $_ -ne $gKey } | ForEach-Object { "$_=$($kv[$_].Count)" }) -join ', ')

    $result += [pscustomobject]@{
        accession=$acc; n_samples=$n; platform=$platform
        group_key=$gKey; group_values=$gText; other_fields=$othersText
    }
}

$outFile = Join-Path $outDir "group-audit.csv"
$result | Export-Csv -LiteralPath $outFile -NoTypeInformation -Encoding UTF8
Write-Output ("written: " + $outFile)
$result | Format-List | Out-String -Width 400 | Write-Output
