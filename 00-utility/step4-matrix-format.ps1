# Inspect ID_REF format and value ranges of each downloaded matrix (determines preprocessing approach)

Add-Type -AssemblyName System.IO.Compression.FileSystem

function Get-MatrixInfo([string]$path) {
    $fs = [System.IO.File]::OpenRead($path)
    $gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
    $sr = New-Object System.IO.StreamReader($gz)

    $nSamples = 0
    $inTable = $false
    $firstIds = @()
    $firstVals = @()
    $rows = 0
    while (-not $sr.EndOfStream) {
        $l = $sr.ReadLine()
        if ($l -match '^!Sample_geo_accession') { $nSamples = ([regex]::Matches($l, 'GSM\d+')).Count; continue }
        if ($l -match '^!series_matrix_table_begin') { $inTable = $true; continue }
        if ($l -match '^!series_matrix_table_end') { break }
        if ($inTable) {
            $rows++
            if ($rows -eq 1) { continue }   # header row: "ID_REF" \t GSM...
            if ($rows -le 4) {
                $parts = $l -split "`t"
                $firstIds += $parts[0].Trim('"')
                $firstVals += $parts[1].Trim('"')
            }
            if ($rows -gt 200) { break }
        }
    }
    $sr.Close(); $gz.Close(); $fs.Close()

    Write-Output ("--- " + (Split-Path $path -Leaf))
    Write-Output ("    samples: $nSamples")
    Write-Output ("    probe IDs (first 3): " + ($firstIds -join ', '))
    Write-Output ("    values (first 3): " + ($firstVals -join ', '))
}

Get-MatrixInfo "codex_output\geo\GSE96804\GSE96804_series_matrix.txt.gz"
Get-MatrixInfo "codex_output\geo\GSE30528\GSE30528_series_matrix.txt.gz"
Get-MatrixInfo "codex_output\geo\GSE30529\GSE30529_series_matrix.txt.gz"
Get-MatrixInfo "codex_output\geo\GSE99325\GSE99340-GPL19184_series_matrix.txt.gz"
