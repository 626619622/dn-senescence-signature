# Inspect PubMed query translation to confirm counts are trustworthy

$qs = @(
    '("diabetic kidney disease"[tiab] OR "diabetic nephropathy"[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR "LASSO"[tiab] OR "XGBoost"[tiab]) AND 2024:2026[dp]',
    '(diabetes[tiab] OR diabetic[tiab]) AND (bioinformatics[tiab] OR "gene expression omnibus"[tiab] OR GEO[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR LASSO[tiab] OR "deep learning"[tiab]) AND 2024:2026[dp]',
    '(diabetes[tiab] OR diabetic[tiab]) AND ("Mendelian randomization"[tiab]) AND 2024:2026[dp]'
)

foreach ($q in $qs) {
    $u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=3&term=" + [uri]::EscapeDataString($q)
    $r = Invoke-RestMethod $u -TimeoutSec 60
    Write-Output ("count = " + $r.esearchresult.count)
    Write-Output ("translation = " + $r.esearchresult.querytranslation)
    Write-Output ("ids = " + ($r.esearchresult.idlist -join ','))
    Write-Output "----"
    Start-Sleep -Milliseconds 500
}
