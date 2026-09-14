# Counts + journal distribution for public-data / ML diabetes research

function Get-Esearch([string]$q, [int]$retmax) {
    $u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=$retmax&sort=relevance&term=" + [uri]::EscapeDataString($q)
    return Invoke-RestMethod $u -TimeoutSec 60
}

$queries = [ordered]@{
    "DKD/ DN + ML, 2024-2026"        = '("diabetic kidney disease"[tiab] OR "diabetic nephropathy"[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR "LASSO"[tiab] OR "XGBoost"[tiab]) AND 2024:2026[dp]'
    "Diabetes + bioinformatics, 2024-2026" = '(diabetes[tiab] OR diabetic[tiab]) AND (bioinformatics[tiab] OR "gene expression omnibus"[tiab] OR GEO[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR LASSO[tiab] OR "deep learning"[tiab]) AND 2024:2026[dp]'
    "Diabetes + MR, 2024-2026"        = '(diabetes[tiab] OR diabetic[tiab]) AND ("Mendelian randomization"[tiab]) AND 2024:2026[dp]'
    "Diabetes + single cell, 2024-2026" = '(diabetes[tiab] OR diabetic[tiab]) AND ("single-cell"[tiab] OR "single cell"[tiab]) AND 2024:2026[dp]'
    "Diabetes + proteomics/Olink, 2023-2026" = '(diabetes[tiab] OR diabetic[tiab]) AND (proteomic*[tiab] OR Olink[tiab] OR "proteome-wide"[tiab]) AND 2023:2026[dp]'
    "Diabetes + machine learning + diagnosis model" = '(diabetes[tiab] OR diabetic[tiab]) AND ("diagnostic model"[tiab] OR "predictive model"[tiab] OR nomogram[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR LASSO[tiab]) AND 2024:2026[dp]'
}

foreach ($k in $queries.Keys) {
    $r = Get-Esearch $queries[$k] 0
    Write-Output ("COUNT | {0} | {1}" -f $k, $r.esearchresult.count)
    Start-Sleep -Milliseconds 400
}

Write-Output ""
Write-Output "=== Journal distribution: diabetes + ML/bioinformatics 2024-2026 (top 200 by relevance) ==="
$q = '(diabetes[tiab] OR diabetic[tiab]) AND (bioinformatics[tiab] OR "gene expression omnibus"[tiab] OR "machine learning"[tiab] OR "random forest"[tiab] OR LASSO[tiab]) AND 2024:2026[dp]'
$r = Get-Esearch $q 200
$ids = $r.esearchresult.idlist
if ($ids.Count -gt 0) {
    $j = Invoke-RestMethod ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&tool=codex&id=" + ($ids -join ',')) -TimeoutSec 60
    $rows = foreach ($id in $ids) {
        $d = $j.result.$id
        if ($d) {
            [pscustomobject]@{
                journal = $d.fulljournalname
                year    = ($d.pubdate -split ' ')[0]
                title   = $d.title
                pmid    = $id
            }
        }
    }
    $rows | Group-Object journal | Sort-Object Count -Descending | Select-Object -First 30 |
        ForEach-Object { Write-Output ("{0,3}  {1}" -f $_.Count, $_.Name) }

    $out = "codex_output\recon\pubmed-diabetes-ml-2024-2026.csv"
    New-Item -ItemType Directory -Force -Path (Split-Path $out) | Out-Null
    $rows | Export-Csv -Path $out -NoTypeInformation -Encoding UTF8
    Write-Output ("saved: " + $out)
}
