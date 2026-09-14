# Download GSE99340 SuperSeries matrices (contains GSE99325 samples) and inspect

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$H = @{ 'User-Agent' = 'codex-geo-skill' }
$base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99340/matrix/"
$dir = "codex_output\geo\GSE99325"
New-Item -ItemType Directory -Force -Path $dir | Out-Null

foreach ($f in @('GSE99340-GPL19109_series_matrix.txt.gz','GSE99340-GPL19184_series_matrix.txt.gz')) {
    $r = Invoke-WebRequest -Uri ($base + $f) -Method Head -Headers $H -UseBasicParsing -TimeoutSec 40
    $len = $r.Headers['Content-Length']
    if ($len -is [array]) { $len = $len[0] }
    Write-Output ("remote size: {0:N1} MB  {1}" -f ([int64]$len / 1MB), $f)

    $out = Join-Path $dir $f
    if (-not (Test-Path -LiteralPath $out)) {
        Invoke-WebRequest -Uri ($base + $f) -Headers $H -OutFile $out -UseBasicParsing -TimeoutSec 1800
        Write-Output ("downloaded: $f (" + [math]::Round((Get-Item -LiteralPath $out).Length/1MB,1) + " MB)")
    }
    Start-Sleep -Milliseconds 300
}

Write-Output ""
Write-Output "=== head of GPL19109 matrix ==="
Add-Type -AssemblyName System.IO.Compression.FileSystem
$p = Join-Path $dir 'GSE99340-GPL19109_series_matrix.txt.gz'
$fs = [System.IO.File]::OpenRead($p)
$gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
$sr = New-Object System.IO.StreamReader($gz)
$i = 0
while ($i -lt 200 -and -not $sr.EndOfStream) {
    $l = $sr.ReadLine()
    if ($l -match 'Series_title|Sample_title|Sample_geo_accession|Sample_source_name|Sample_characteristics|!series_matrix_table_begin|ID_REF') {
        $s = if ($l.Length -gt 200) { $l.Substring(0,200) + '...' } else { $l }
        Write-Output ("  " + $s)
    }
    $i++
}
$sr.Close(); $gz.Close(); $fs.Close()
