# Step 1 novelty scan: which angles in DKD/ML literature are saturated vs open (2024-2026)
# ASCII-only

$ErrorActionPreference = 'Stop'
$outDir = "codex_output\recon"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$q = '("diabetic kidney disease"[tiab] OR "diabetic nephropathy"[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR "LASSO"[tiab] OR "XGBoost"[tiab] OR "bioinformatics"[tiab] OR "gene expression omnibus"[tiab] OR "single-cell"[tiab]) AND 2024:2026[dp]'
$u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=600&sort=pub_date&term=" + [uri]::EscapeDataString($q)
$r = Invoke-RestMethod $u -TimeoutSec 90
$ids = $r.esearchresult.idlist
Write-Output ("total hits: " + $r.esearchresult.count + " | fetched ids: " + $ids.Count)

$xmlPath = Join-Path $outDir "dkd-ml-abstracts-2024-2026.xml"
$chunkSize = 150
$allRecords = @()
for ($i = 0; $i -lt $ids.Count; $i += $chunkSize) {
    $end = [Math]::Min($i + $chunkSize - 1, $ids.Count - 1)
    $chunk = $ids[$i..$end]
    $u2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&rettype=abstract&tool=codex&id=" + ($chunk -join ',')
    $xmlText = (Invoke-WebRequest $u2 -UseBasicParsing -TimeoutSec 180).Content
    if ($i -eq 0) { [System.IO.File]::WriteAllText((Join-Path (Get-Location) $xmlPath), $xmlText, (New-Object System.Text.UTF8Encoding($false))) }
    $doc = [xml]$xmlText
    $allRecords += $doc.PubmedArticleSet.PubmedArticle
    Write-Output ("fetched chunk: " + $chunk.Count)
    Start-Sleep -Milliseconds 400
}

$angles = [ordered]@{
    "ferroptosis"        = 'ferroptosis'
    "cuproptosis"        = 'cuproptosis'
    "disulfidptosis"     = 'disulfidptosis'
    "pyroptosis"         = 'pyroptosis'
    "panoptosis"         = 'panoptosis'
    "m6A / methylation"  = 'm6a|methylation|epigenetic'
    "mitochondrial"      = 'mitochondri'
    "autophagy"          = 'autophagy|mitophagy'
    "ER stress"          = 'endoplasmic reticulum stress|er stress'
    "immune infiltration"= 'immune infiltrat|immune landscape|immune cell|immune microenvironment'
    "single-cell"        = 'single-cell|single cell|scrna|snrna'
    "Mendelian"          = 'mendelian'
    "drug prediction"    = 'drug|compound|molecular docking|cmap|therapeutic target'
    "diagnostic model"   = 'diagnos'
    "progression/prognosis" = 'progress|prognos|predict'
    "subtype / clustering" = 'subtype|consensus cluster|clustering|subcluster'
    "podocyte"           = 'podocyte'
    "tubular"            = 'tubul'
    "glomerular"         = 'glomerul'
    "biomarker"          = 'biomarker'
    "microbiome"         = 'microbiota|microbiome|gut'
    "exosome"            = 'exosom'
    "senescence"         = 'senescence'
    "lactylation"        = 'lactyl'
    "hypoxia"            = 'hypoxia'
    "macrophage"         = 'macrophage'
    "endothelial"        = 'endothelial'
    "multi-omics"        = 'multi-omics|multiomics|integrat.*omics'
    "deep learning"      = 'deep learning|neural network'
    "SHAP / explainable" = 'shap|explainab'
    "external validation"= 'external validation|validation cohort|independent cohort'
    "comorbidity"        = 'sepsis|covid|heart failure|retinopathy|sarcopenia|obesity'
}

$rows = @()
foreach ($a in $allRecords) {
    $pmid = $a.MedlineCitation.PMID.'#text'
    if (-not $pmid) { $pmid = [string]$a.MedlineCitation.PMID }
    $title = [string]$a.MedlineCitation.Article.ArticleTitle
    $journal = [string]$a.MedlineCitation.Article.Journal.Title
    $year = [string]$a.MedlineCitation.Article.Journal.JournalIssue.PubDate.Year
    $abs = ""
    $absNode = $a.MedlineCitation.Article.Abstract.AbstractText
    if ($absNode) { $abs = ($absNode | ForEach-Object { if ($_.'#text') { $_.'#text' } else { [string]$_ } }) -join ' ' }
    $text = ($title + " " + $abs).ToLower()
    $titleLower = $title.ToLower()
    $hit = @()
    foreach ($k in $angles.Keys) {
        if ($titleLower -match $angles[$k]) { $hit += $k }
    }
    $hitAll = @()
    foreach ($k in $angles.Keys) {
        if ($text -match $angles[$k]) { $hitAll += $k }
    }
    $rows += [pscustomobject]@{
        pmid = $pmid; year = $year; journal = $journal; title = $title
        angles_in_title = ($hit -join '; ')
        angles_in_text  = ($hitAll -join '; ')
    }
}

$csv = Join-Path $outDir "dkd-ml-records-2024-2026.csv"
$rows | Export-Csv -Path $csv -NoTypeInformation -Encoding UTF8
Write-Output ("records saved: " + $csv + " (n=" + $rows.Count + ")")

Write-Output ""
Write-Output "=== Angle frequency in TITLES (n = total records) ==="
$summary = foreach ($k in $angles.Keys) {
    $inTitle = ($rows | Where-Object { $_.angles_in_title -like "*$k*" }).Count
    $inText  = ($rows | Where-Object { $_.angles_in_text -like "*$k*" }).Count
    [pscustomobject]@{ angle = $k; in_title = $inTitle; in_title_or_abstract = $inText }
}
$summary | Sort-Object in_title -Descending | Format-Table -AutoSize | Out-String -Width 200 | Write-Output
$summary | Export-Csv -Path (Join-Path $outDir "dkd-angle-frequency.csv") -NoTypeInformation -Encoding UTF8
