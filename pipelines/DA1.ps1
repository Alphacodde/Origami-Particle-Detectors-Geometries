param(
    [string]$RootPath = ".",
    [string]$ExtractScriptPath = "extract_hits_for_deadzone.py",
    [switch]$Recurse
)
if (-not (Test-Path -LiteralPath $RootPath)) {
    Write-Error "RootPath '$RootPath' does not exist."
    exit 1
}
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
$exitCode = $LASTEXITCODE
if ($exitCode -ne 0) {
    Write-Error "extract_hits_for_deadzone.py reported failures (exit code $exitCode). See BATCH SUMMARY above."
    exit $exitCode
}
Write-Host "`nDone. A <name>_hits.csv was written alongside each source .root file." -ForegroundColor Cyan
