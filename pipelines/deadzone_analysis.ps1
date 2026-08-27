param(
    [string]$RootPath = ".",
    [string]$ExtractScriptPath = "extract_hits_for_deadzone.py",
    [string]$DeadZoneScriptPath = "deadzone_map.py",
    [double]$Threshold = 0.5,
    [switch]$Recurse,
    [switch]$SkipExtraction
)
if (-not (Test-Path -LiteralPath $RootPath)) {
    Write-Error "RootPath '$RootPath' does not exist."
    exit 1
}
if (-not (Test-Path -LiteralPath $DeadZoneScriptPath)) {
    Write-Error "Dead-zone script '$DeadZoneScriptPath' not found. Pass -DeadZoneScriptPath if it's elsewhere."
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
Write-Host "`nRunning: python $DeadZoneScriptPath original `"$RootPath`" `"$RootPath`" --threshold $Threshold" -ForegroundColor Cyan
python $DeadZoneScriptPath original $RootPath $RootPath --threshold $Threshold
if ($LASTEXITCODE -ne 0) {
    Write-Error "deadzone_map.py failed (exit code $LASTEXITCODE)."
    exit $LASTEXITCODE
}
Write-Host "`nDone. Summary tables written to:" -ForegroundColor Cyan
Write-Host "  $RootPath\original_four_dead_zone_summary.csv"
Write-Host "  $RootPath\original_four_dead_zone_summary_aggregate.csv"
