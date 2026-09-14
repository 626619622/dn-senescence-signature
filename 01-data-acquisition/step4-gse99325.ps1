# Inspect GSE99325 supplementary inventory and sizes (no series matrix available)

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$H = @{ 'User-Agent' = 'codex-geo-skill' }
$base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99325/suppl/"

$html = (Invoke-WebRequest $base -Headers $H -UseBasicParsing -TimeoutSec 60).Content
$files = [regex]::Matches($html, 'href="([^"/?][^"]*)"') | ForEach-Object { $_.Groups[1].Value } |
         Where-Object { $_ -notmatch '^[a-zA-Z]+:' -and $_ -ne '..' } | Select-Object -Unique | Sort-Object

Write-Output ("supplementary files: " + $files.Count)
foreach ($f in $files) {
    try {
        $r = Invoke-WebRequest -Uri ($base + $f) -Method Head -Headers $H -UseBasicParsing -TimeoutSec 40
        $len = $r.Headers['Content-Length']
        if ($len -is [array]) { $len = $len[0] }
        if ($len) { Write-Output ("{0,12:N1} MB  {1}" -f ([int64]$len / 1MB), $f) }
        else { Write-Output ("       ?      MB  " + $f) }
    } catch {
        Write-Output ("   HEAD-FAIL      " + $f + " -> " + $_.Exception.Message)
    }
    Start-Sleep -Milliseconds 300
}

Write-Output ""
Write-Output "=== also check GSE99340 (SuperSeries) matrix dir ==="
$mbase = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99340/matrix/"
try {
    $html2 = (Invoke-WebRequest $mbase -Headers $H -UseBasicParsing -TimeoutSec 40).Content
    $mf = [regex]::Matches($html2, 'href="([^"/?][^"]*)"') | ForEach-Object { $_.Groups[1].Value } | Select-Object -Unique
    $mf | ForEach-Object { Write-Output ("   " + $_) }
} catch { Write-Output ("   matrix dir error: " + $_.Exception.Message) }
