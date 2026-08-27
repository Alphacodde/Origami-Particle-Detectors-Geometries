param(
    [ValidateSet("Release", "Debug")]
    [string]$BuildConfig = "Release",
    [ValidateSet("fixed", "scan", "both")]
    [string]$Mode = "both",
    [switch]$SkipGeometryGen,
    [switch]$SkipValidation,
    [switch]$Build,
    [double]$SiThicknessMM = 0.300,
    [double]$KaptonThicknessMM = 0.050,
    [double]$MomentumGeV = 5.0,
    [double]$GunZmm = 400,
    [double]$DiskRadiusMm = 150,
    [double]$VertexSmearMm = 3.0,
    [int]$NFixed = 1000,
    [int]$NScan = 500000
)
$env:PYTHONUTF8 = "1"
$ErrorActionPreference = "Stop"
$ProjectRoot = $PSScriptRoot
if (-not $ProjectRoot) { $ProjectRoot = Get-Location }
$BuildOutDir  = Join-Path $ProjectRoot "build\$BuildConfig"
$ExePath      = Join-Path $BuildOutDir "origamiDet.exe"
$ResultsDir   = Join-Path $ProjectRoot "results_diff_geom"
$GenMacroDir  = Join-Path $BuildOutDir "macros\_generated"
$GeomSrcDir   = Join-Path $ProjectRoot "geom_diff_mesh"
function Say($msg) { Write-Host $msg }
function SayStep($msg) { Write-Host "`n== $msg ==" -ForegroundColor Cyan }
function Invoke-Checked([string]$exe, [string[]]$cmdArgs, [string]$onFailMsg) {
    & $exe @cmdArgs
    if ($LASTEXITCODE -ne 0) {
        Write-Error "$onFailMsg (exit code $LASTEXITCODE)"
        exit 1
    }
}
function Get-UniqueDest([string]$dir, [string]$name) {
    $dest = Join-Path $dir $name
    if (-not (Test-Path $dest)) { return $dest }
    $stem = [System.IO.Path]::GetFileNameWithoutExtension($name)
    $ext  = [System.IO.Path]::GetExtension($name)
    $counter = 2
    while (Test-Path (Join-Path $dir "${stem}_${counter}${ext}")) { $counter++ }
    return Join-Path $dir "${stem}_${counter}${ext}"
}
function Move-LatestRoot([string]$tag, [string]$destPrefix) {
    $found = Get-ChildItem -Path $BuildOutDir -Filter "origami_${tag}_fold*_run*.root" -File -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending | Select-Object -First 1
    if (-not $found) {
        Write-Warning "  no origami_${tag}_fold*_run*.root produced - check the log above."
        return $null
    }
    $destName = "${destPrefix}_$($found.Name)"
    $dest = Get-UniqueDest $ResultsDir $destName
    Move-Item -Path $found.FullName -Destination $dest -Force
    Say "  -> moved to $(Split-Path $ResultsDir -Leaf)\$(Split-Path $dest -Leaf)"
    return $dest
}
Say "== origamiDet: run_everything (GDML) =="
Say "  ProjectRoot  = $ProjectRoot"
Say "  BuildOutDir  = $BuildOutDir"
Say "  Mode         = $Mode"
Say "  Si / Kapton  = $SiThicknessMM mm / $KaptonThicknessMM mm"
Say "  Momentum     = $MomentumGeV GeV"
Say "  Fixed mode gun (Mark 1 cone): gunZmm=$GunZmm mm, diskRadiusMm=$DiskRadiusMm mm"
Say "  Scan mode gun (Mark 2 isotropic4pi): vertexSmearMm=$VertexSmearMm mm"
Say "  N_Fixed / N_Scan = $NFixed / $NScan"
Say "  Shared geometry envelope: see diff_geom_macros\_geometry_config.py (SHARED_R_MM/SHARED_L_MM)"
Say "  Sweep now includes comparison_barrel.py's barrel_reference (R_f=1, K=0 baseline)"
Say "  Geometry format: GDML (indexed tessellated solids) - see this script's own header comment"
New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null
New-Item -ItemType Directory -Force -Path $GenMacroDir | Out-Null
New-Item -ItemType Directory -Force -Path $GeomSrcDir | Out-Null
if ($Build) {
    SayStep "Building ($BuildConfig)"
    Push-Location (Join-Path $ProjectRoot "build")
    cmake --build . --config $BuildConfig --parallel
    Pop-Location
}
if (-not (Test-Path $ExePath)) {
    Write-Error "origamiDet.exe not found at $ExePath - build it first (or pass -Build)."
    exit 1
}
if (-not $SkipGeometryGen) {
    SayStep "Step 1: generating geometry (GDML)"
    $env:ORIGAMIDET_GEOM_OUTPUT_DIR = $GeomSrcDir
    Say "  Generator output dir (ORIGAMIDET_GEOM_OUTPUT_DIR) = $GeomSrcDir"
    Push-Location $ProjectRoot
    try {
        Say "-- kresling.py --"
        Invoke-Checked "python" @("diff_geom_macros\kresling.py") "kresling.py failed"
        Say "-- yoshimura.py --"
        Invoke-Checked "python" @("diff_geom_macros\yoshimura.py") "yoshimura.py failed"
        Say "-- miura_cylindrical.py --"
        Invoke-Checked "python" @("diff_geom_macros\miura_cylindrical.py") "miura_cylindrical.py failed"
        Say "-- comparison_barrel.py --"
        Invoke-Checked "python" @("diff_geom_macros\comparison_barrel.py") "comparison_barrel.py failed"
    } finally {
        Pop-Location
        Remove-Item Env:\ORIGAMIDET_GEOM_OUTPUT_DIR -ErrorAction SilentlyContinue
    }
    $geomSrcFiles = Get-ChildItem -Path $GeomSrcDir -Filter "*_silicon.gdml" -File -ErrorAction SilentlyContinue
    if (-not $geomSrcFiles -or $geomSrcFiles.Count -eq 0) {
        Write-Error "Step 1 completed but $GeomSrcDir contains no *_silicon.gdml files - geometry generation did not write where this script expects, or is still writing .stl (check diff_geom_macros\_solid_export.py's export_solid() is the version that writes GDML). Check ORIGAMIDET_GEOM_OUTPUT_DIR support in the generator scripts before continuing."
        exit 1
    }
    Say "  Step 1 OK: $($geomSrcFiles.Count) silicon GDML(s) found in $GeomSrcDir"
} else {
    SayStep "Step 1: SKIPPED (-SkipGeometryGen)"
}
SayStep "Step 1b: mesh QA - SKIPPED (prep_stl.py is STL-specific, not yet updated for GDML; see this script's own comment above the console output already printed by export_solid() in Step 1)"
SayStep "Step 2: syncing geom_diff_mesh/macros into $BuildOutDir"
$BuildGeomDir = Join-Path $BuildOutDir "geometry"
New-Item -ItemType Directory -Force -Path $BuildGeomDir | Out-Null
$copySrcCount = (Get-ChildItem -Path $GeomSrcDir -File -Recurse -ErrorAction SilentlyContinue | Measure-Object).Count
if ($copySrcCount -eq 0) {
    Write-Error "$GeomSrcDir is empty - nothing to sync into $BuildGeomDir. Re-run Step 1 (do not pass -SkipGeometryGen) before continuing."
    exit 1
}
Copy-Item -Path (Join-Path $GeomSrcDir "*") -Destination $BuildGeomDir -Recurse -Force
Copy-Item -Path (Join-Path $ProjectRoot "macros\*.mac") -Destination (Join-Path $BuildOutDir "macros") -Force
$copiedSilicon = (Get-ChildItem -Path $BuildGeomDir -Filter "*_silicon.gdml" -File -ErrorAction SilentlyContinue | Measure-Object).Count
Say "  Step 2 OK: $copySrcCount file(s) synced from $(Split-Path $GeomSrcDir -Leaf); $copiedSilicon silicon GDML(s) now in $BuildGeomDir"
if (-not $SkipValidation) {
    SayStep "Step 3: Tier-1 validation (flat plate)"
    $plateShell = Join-Path $ProjectRoot "plate_shell.stl"
    Invoke-Checked "python" @((Join-Path $ProjectRoot "tools\make_validation_plate.py"), "--side", "50", "--out", $plateShell) `
        "make_validation_plate.py failed"
    Invoke-Checked "python" @((Join-Path $ProjectRoot "tools\thicken_stack.py"), $plateShell, `
        "--layer", $SiThicknessMM, (Join-Path $ProjectRoot "geometry\plate_silicon.gdml"), "silicon", `
        "--layer", $KaptonThicknessMM, (Join-Path $ProjectRoot "geometry\plate_kapton.gdml"), "kapton") `
        "thicken_stack.py failed on the validation plate"
    Remove-Item $plateShell -ErrorAction SilentlyContinue
    Copy-Item -Path (Join-Path $ProjectRoot "geometry\plate_silicon.gdml"), (Join-Path $ProjectRoot "geometry\plate_kapton.gdml") `
              -Destination (Join-Path $BuildOutDir "geometry") -Force
    Push-Location $BuildOutDir
    & ".\origamiDet.exe" "macros\validate_flat_plate.mac"
    Pop-Location
    $validationRoot = Get-ChildItem -Path $BuildOutDir -Filter "origami_plate_fold*_run*.root" -File -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending | Select-Object -First 1
    if (-not $validationRoot) {
        Write-Error "Validation run produced no origami_plate_*.root - check the log above."
        exit 1
    }
    python (Join-Path $ProjectRoot "tools\check_validation.py") $validationRoot.FullName --expected-mm $SiThicknessMM
    if ($LASTEXITCODE -ne 0) {
        Write-Error "ABORTING: Tier-1 flat-plate validation failed. Do not trust any real-geometry result until this passes - see README_FIXES.md."
        exit 1
    }
    $dest = Get-UniqueDest $ResultsDir $validationRoot.Name
    Move-Item -Path $validationRoot.FullName -Destination $dest -Force
    Say "Validation PASSED. Moved to $(Split-Path $ResultsDir -Leaf)\$(Split-Path $dest -Leaf)"
} else {
    SayStep "Step 3: Tier-1 validation SKIPPED (-SkipValidation, not recommended)"
}
SayStep "Step 4: sweep"
$siliconFiles = Get-ChildItem -Path $BuildGeomDir -Filter "*_silicon.gdml" -File |
    Where-Object { $_.BaseName -ne "plate_silicon" }
if ($siliconFiles.Count -eq 0) {
    Write-Warning "No *_silicon.gdml files found in geometry\ (besides the plate) - nothing to sweep."
}
$nRun = 0
foreach ($si in $siliconFiles) {
    $tag = $si.BaseName -replace "_silicon$", ""
    $kaptonName = "${tag}_kapton.gdml"
    $kaptonPath = Join-Path $si.DirectoryName $kaptonName
    if (-not (Test-Path $kaptonPath)) {
        Write-Warning "No matching $kaptonName for $($si.Name) - skipping."
        continue
    }
    Say "`n-- $tag --"
    $geomLines = @(
        "/origami/stlFile geometry/$($si.Name)"
        "/origami/kaptonStlFile geometry/$kaptonName"
        "/origami/sensorThickness $SiThicknessMM mm"
        "/origami/substrateThickness $KaptonThicknessMM mm"
        "/origami/foldState 1.0"
        ""
    )
    $tailLines = @(
        ""
        "/run/verbose 0"
        "/event/verbose 0"
        "/tracking/verbose 0"
        ""
        "/run/initialize"
    )
    if ($Mode -eq "fixed" -or $Mode -eq "both") {
        $macLines = $geomLines + @(
            "/origami/gun/mode mark1cone"
            "/origami/gun/momentumGeV $MomentumGeV"
            "/origami/gun/gunZmm $GunZmm"
            "/origami/gun/diskRadiusMm $DiskRadiusMm"
            "/origami/gun/fixedAngleDeg 0"
        ) + $tailLines + @(
            "/run/beamOn $NFixed"
        )
        $macPath = Join-Path $GenMacroDir "${tag}_fixed.mac"
        $macLines | Set-Content -Path $macPath
        Push-Location $BuildOutDir
        & ".\origamiDet.exe" "macros\_generated\${tag}_fixed.mac"
        Pop-Location
        Move-LatestRoot -tag $tag -destPrefix "fixed"
        $nRun++
    }
    if ($Mode -eq "scan" -or $Mode -eq "both") {
        $macLines = $geomLines + @(
            "/origami/gun/mode isotropic4pi"
            "/origami/gun/momentumGeV $MomentumGeV"
            "/origami/gun/vertexSmearMm $VertexSmearMm"
        ) + $tailLines + @(
            "/run/beamOn $NScan"
        )
        $macPath = Join-Path $GenMacroDir "${tag}_scan.mac"
        $macLines | Set-Content -Path $macPath
        Push-Location $BuildOutDir
        & ".\origamiDet.exe" "macros\_generated\${tag}_scan.mac"
        Pop-Location
        Move-LatestRoot -tag $tag -destPrefix "scan"
        $nRun++
    }
}
SayStep "Done: $nRun run(s) completed, results in $(Split-Path $ResultsDir -Leaf)\"
