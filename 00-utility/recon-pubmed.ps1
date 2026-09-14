# Reconnaissance: which journals publish public-data + machine learning diabetes work
# ASCII-only script (Windows PowerShell 5.1 safe)

function Get-Count([string]$q) {
    $u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=0&term=" + [uri]::EscapeDataString($q)
    try {
        $r = Invoke-RestMethod $u -TimeoutSec 40
        return $r.esearchresult.count
    } catch {
        return "NETERR"
    }
}

$queries = [ordered]@{
    "diabetic kidney disease + machine learning 2024-2026" = "('diabetic kidney disease'[tiab] OR 'diabetic nephropathy'[tiab]) AND ('machine learning'[tiab] OR 'random forest'[tiab] OR 'LASSO'[tiab]) AND ('2024'[dp] : '2026'[dp])"
    "diabetes + GEO bioinformatics + machine learning"     = "('diabetes mellitus'[tiab] OR 'diabetic'[tiab]) AND ('bioinformatics'[tiab] OR 'GEO'[tiab] OR 'gene expression omnibus'[tiab]) AND ('machine learning'[tiab]) AND ('2024'[dp] : '2026'[dp])"
    "diabetes + Mendelian randomization 2024-2026"         = "('diabetes'[tiab] OR 'diabetic'[tiab]) AND ('Mendelian randomization'[tiab]) AND ('2024'[dp] : '2026'[dp])"
    "diabetic + single cell + machine learning"            = "('diabetic'[tiab]) AND ('single-cell'[tiab] OR 'single cell'[tiab]) AND ('machine learning'[tiab] OR 'deep learning'[tiab]) AND ('2024'[dp] : '2026'[dp])"
}

foreach ($k in $queries.Keys) {
    $c = Get-Count $queries[$k]
    Write-Output ("{0} => {1}" -f $k, $c)
    Start-Sleep -Milliseconds 400
}
