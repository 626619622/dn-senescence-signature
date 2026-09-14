# Verify that GSE99325 biopsy samples are present in the downloaded GSE99340 matrices

Add-Type -AssemblyName System.IO.Compression.FileSystem

$targets = (Import-Csv "codex_output\geo\GSE99325\GSE99325-samples.csv" | ForEach-Object { $_.gsm })
Write-Output ("GSE99325 samples in metadata: " + $targets.Count)

foreach ($f in @('GSE99340-GPL19109_series_matrix.txt.gz','GSE99340-GPL19184_series_matrix.txt.gz')) {
    $p = "codex_output\geo\GSE99325\$f"
    $fs = [System.IO.File]::OpenRead($p)
    $gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
    $sr = New-Object System.IO.StreamReader($gz)

    $gsms = @()
    $titles = @()
    while (-not $sr.EndOfStream) {
        $l = $sr.ReadLine()
        if ($l -match '^!Sample_geo_accession') { $gsms += ([regex]::Matches($l, '(GSM\d+)') | ForEach-Object { $_.Groups[1].Value }) }
        elseif ($l -match '^!Sample_title') { $titles += ([regex]::Matches($l, '"([^"]*)"') | ForEach-Object { $_.Groups[1].Value }) }
        elseif ($l -match '^!series_matrix_table_begin') { break }
    }
    $sr.Close(); $gz.Close(); $fs.Close()

    $shared = @($gsms | Where-Object { $targets -contains $_ })
    Write-Output ("--- " + $f)
    Write-Output ("    samples in matrix: " + $gsms.Count)
    Write-Output ("    GSE99325 samples present: " + $shared.Count)
    if ($shared.Count -gt 0) {
        $idx = 0..($gsms.Count-1) | Where-Object { $targets -contains $gsms[$_] }
        Write-Output ("    sample titles (first 12):")
        $idx | Select-Object -First 12 | ForEach-Object { Write-Output ("       " + $gsms[$_] + "  " + $titles[$_]) }
    }
}
