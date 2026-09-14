# Step 3c: inspect what is actually inside the downloaded archives/matrices

Add-Type -AssemblyName System.IO.Compression.FileSystem

function Show-GzHead([string]$path, [int]$lines) {
    $fs = [System.IO.File]::OpenRead($path)
    $gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
    $sr = New-Object System.IO.StreamReader($gz)
    Write-Output ("--- " + (Split-Path $path -Leaf))
    for ($i = 0; $i -lt $lines; $i++) {
        $l = $sr.ReadLine()
        if ($null -eq $l) { break }
        if ($l.Length -gt 160) { $l = $l.Substring(0, 160) + '...' }
        Write-Output ("    " + $l)
    }
    $sr.Close(); $gz.Close(); $fs.Close()
}

foreach ($p in @('codex_output\geo\GSE294519\GSE294519_series_matrix.txt.gz',
                 'codex_output\geo\GSE175759\GSE175759_series_matrix.txt.gz',
                 'codex_output\geo\GSE189005\GSE189005_series_matrix.txt.gz')) {
    Show-GzHead $p 6
    Write-Output ""
}

Write-Output "=== GSE142025_RAW.tar contents ==="
$tar = [System.IO.Compression.ZipFile]::OpenRead((Resolve-Path 'codex_output\geo\GSE142025\GSE142025_RAW.tar'))
Write-Output ("not a zip (expected for tar). entries via tar list:")
$tar.Dispose()

Write-Output ""
Write-Output "=== sizes of everything downloaded so far ==="
$root = (Get-Location).Path
Get-ChildItem "codex_output\geo" -Recurse -File |
    Where-Object { $_.Length -gt 1000 } |
    Sort-Object Length -Descending |
    ForEach-Object { Write-Output ("{0,10:N1} MB  {1}" -f ($_.Length/1MB), $_.FullName.Replace($root + '\','')) }
Write-Output ("TOTAL: " + [math]::Round(((Get-ChildItem "codex_output\geo" -Recurse -File | Measure-Object Length -Sum).Sum)/1MB,1) + " MB")
