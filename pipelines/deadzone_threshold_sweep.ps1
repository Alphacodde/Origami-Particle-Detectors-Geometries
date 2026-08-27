param(
    [string]$RootPath = ".",
    [string]$ExtractScriptPath = "extract_hits_for_deadzone.py",
    [string]$SweepScriptPath = "deadzone_threshold_sweep.py",
    [int]$NThresholds = 101,
    [double]$MarkThreshold = 1.0,
    [string]$OutCsv = "dead_zone_threshold_sweep.csv",
    [string]$OutPlot = "dead_zone_threshold_sweep.png",
    [switch]$Recurse,
    [switch]$SkipExtraction
)
if (-not (Test-Path -LiteralPath $RootPath)) {
    Write-Error "RootPath '$RootPath' does not exist."
    exit 1
}
if (-not (Test-Path -LiteralPath $SweepScriptPath)) {
    Write-Error "Sweep script '$SweepScriptPath' not found. Pass -SweepScriptPath if it's elsewhere."
    exit 1
}
if ($SkipExtraction) {
    Write-Host "Skipping extraction (-SkipExtraction set); using existing hits CSVs under '$RootPath'." -ForegroundColor Yellow
} else {
    if (-not (Test-Path -LiteralPath $ExtractScriptPath)) {
        Write-Error "Extraction script '$ExtractScriptPath' not found. Pass -ExtractScriptPath if it's elsewhere."
        exit 1
    }
    $rootFiles = Get-ChildItem -LiteralPath $RootPath -Filter "*.root" -File -Recurse:$Recurse |
        Select-Object -ExpandProperty FullName
    if (-not $rootFiles -or $rootFiles.Count -eq 0) {
        Write-Error "No .root files found under '$RootPath'$(if ($Recurse) { ' (recursively)' })."
        exit 1
    }
    Write-Host "Found $($rootFiles.Count) .root file(s) under '$RootPath'." -ForegroundColor Cyan
    Write-Host "Running: python `"$ExtractScriptPath`" <$($rootFiles.Count) file(s)>" -ForegroundColor Cyan
    python $ExtractScriptPath @rootFiles
    if ($LASTEXITCODE -ne 0) {
        Write-Error "extract_hits_for_deadzone.py reported failures (exit code $LASTEXITCODE). See BATCH SUMMARY above."
        exit $LASTEXITCODE
    }
}
$sweepOutCsv = Join-Path $RootPath $OutCsv
$sweepOutPlot = Join-Path $RootPath $OutPlot
Write-Host "`nRunning: python $SweepScriptPath `"$RootPath`" `"$RootPath`" --n-thresholds $NThresholds --mark-threshold $MarkThreshold --out-csv `"$sweepOutCsv`" --out-plot `"$sweepOutPlot`"" -ForegroundColor Cyan
python $SweepScriptPath $RootPath $RootPath --n-thresholds $NThresholds --mark-threshold $MarkThreshold --out-csv $sweepOutCsv --out-plot $sweepOutPlot
if ($LASTEXITCODE -ne 0) {
    Write-Error "deadzone_threshold_sweep.py failed (exit code $LASTEXITCODE)."
    exit $LASTEXITCODE
}
Write-Host "`nDone. Threshold-sweep outputs written to:" -ForegroundColor Cyan
Write-Host "  $sweepOutCsv"
Write-Host "  $sweepOutPlot"
