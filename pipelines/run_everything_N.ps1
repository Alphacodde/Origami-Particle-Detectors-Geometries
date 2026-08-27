param(
    [ValidateSet("Release", "Debug")]
    [string]$BuildConfig = "Release",
    [ValidateSet("fixed", "scan", "both")]
    [string]$Mode = "both",
    [int[]]$NValues = @(6, 8, 12, 16, 24, 32),
    [ValidateSet("kresling", "yoshimura", "miura", "barrel")]
    [string[]]$Structures = @("kresling", "yoshimura", "miura", "barrel"),
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
$BuildOutDir   = Join-Path $ProjectRoot "build\$BuildConfig"
$ExePath       = Join-Path $BuildOutDir "origamiDet.exe"
$ResultsDir    = Join-Path $ProjectRoot "Results_N_Sweep"
$GenMacroDir   = Join-Path $BuildOutDir "macros\_generated"
$NSweepDir     = Join-Path $ProjectRoot "geometry_N_sweep"
$NSweepBuildDir = Join-Path $BuildOutDir "geometry_N_sweep"
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
    Say "  -> moved to Results_N_Sweep\$(Split-Path $dest -Leaf)"
    return $dest
}
Say "== origamiDet: run_everything_N (facet-count sweep) =="
Say "  ProjectRoot   = $ProjectRoot"
Say "  BuildOutDir   = $BuildOutDir"
Say "  NSweepDir     = $NSweepDir"
Say "  Mode          = $Mode"
Say "  N values      = $($NValues -join ', ')"
Say "  Structures    = $($Structures -join ', ')"
Say "  Si / Kapton   = $SiThicknessMM mm / $KaptonThicknessMM mm"
Say "  Momentum      = $MomentumGeV GeV"
Say "  Fixed mode gun (Mark 1 cone): gunZmm=$GunZmm mm, diskRadiusMm=$DiskRadiusMm mm"
Say "  Scan mode gun (Mark 2 isotropic4pi): vertexSmearMm=$VertexSmearMm mm"
Say "  N_Fixed / N_Scan = $NFixed / $NScan"
Say "  Shared geometry envelope: see diff_geom_macros\_geometry_config.py (SHARED_R_MM/SHARED_L_MM)"
Say "  Axial density (layers/m_rings/rows) held at each structure's CONFIG default throughout - only circumferential N varies (see sweep_geometry_by_N.py docstring)"
Say "  Geometry format: GDML (indexed tessellated solids) - see this script's own header comment"
New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null
New-Item -ItemType Directory -Force -Path $GenMacroDir | Out-Null
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
    SayStep "Step 1: generating N-sweep geometry (geometry_N_sweep\, GDML)"
    Push-Location $ProjectRoot
    try {
        Say "-- sweep_geometry_by_N.py --"
        Invoke-Checked "python" @("diff_geom_macros\sweep_geometry_by_N.py") "sweep_geometry_by_N.py failed"
    } finally {
        Pop-Location
    }
} else {
    SayStep "Step 1: SKIPPED (-SkipGeometryGen)"
}
if (-not (Test-Path $NSweepDir)) {
    Write-Error "$NSweepDir does not exist - run without -SkipGeometryGen at least once first."
    exit 1
}
$nSweepSiliconFiles = Get-ChildItem -Path $NSweepDir -Filter "*_silicon.gdml" -File -ErrorAction SilentlyContinue
if (-not $nSweepSiliconFiles -or $nSweepSiliconFiles.Count -eq 0) {
    Write-Error "$NSweepDir contains no *_silicon.gdml files - either geometry generation did not write where this script expects, or sweep_geometry_by_N.py is still writing .stl (see this script's own header ASSUMPTION note - it may need updating to import diff_geom_macros\_solid_export.py's GDML export path, the same way tools\thicken_stack.py was updated for Step 3 below)."
    exit 1
}
Say "  Step 1 OK: $($nSweepSiliconFiles.Count) silicon GDML(s) found in $NSweepDir"
SayStep "Step 1b: mesh QA - SKIPPED (prep_stl.py is STL-specific, not yet updated for GDML)"
SayStep "Step 2: syncing geometry_N_sweep/macros into $BuildOutDir"
$nSweepCopySrcCount = (Get-ChildItem -Path $NSweepDir -File -Recurse -ErrorAction SilentlyContinue | Measure-Object).Count
if ($nSweepCopySrcCount -eq 0) {
    Write-Error "$NSweepDir is empty - nothing to sync into $BuildOutDir. Re-run Step 1 (do not pass -SkipGeometryGen) before continuing."
    exit 1
}
Copy-Item -Path $NSweepDir -Destination $BuildOutDir -Recurse -Force
Copy-Item -Path (Join-Path $ProjectRoot "macros\*.mac") -Destination (Join-Path $BuildOutDir "macros") -Force
$nSweepCopiedSilicon = (Get-ChildItem -Path $NSweepBuildDir -Filter "*_silicon.gdml" -File -ErrorAction SilentlyContinue | Measure-Object).Count
Say "  Step 2 OK: $nSweepCopySrcCount file(s) synced from $(Split-Path $NSweepDir -Leaf); $nSweepCopiedSilicon silicon GDML(s) now in $NSweepBuildDir"
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
    New-Item -ItemType Directory -Force -Path (Join-Path $BuildOutDir "geometry") | Out-Null
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
    Say "Validation PASSED. Moved to Results_N_Sweep\$(Split-Path $dest -Leaf)"
} else {
    SayStep "Step 3: Tier-1 validation SKIPPED (-SkipValidation, not recommended)"
}
SayStep "Step 4: N-sweep"
$allSilicon = Get-ChildItem -Path $NSweepDir -Filter "*_silicon.gdml" -File
$siliconFiles = $allSilicon | Where-Object {
    $base = $_.BaseName -replace "_silicon$", ""
    if ($base -eq "barrel_reference") {
        return $Structures -contains "barrel"
    }
    if ($base -match "^(kresling|yoshimura|miura)_N(\d+)$") {
        $structureMatch = $Matches[1]
        $nMatch = [int]$Matches[2]
        return ($Structures -contains $structureMatch) -and ($NValues -contains $nMatch)
    }
    Write-Warning "Unrecognized N-sweep filename pattern, skipping: $($_.Name)"
    return $false
}
if ($siliconFiles.Count -eq 0) {
    Write-Warning "No *_silicon.gdml files matched -Structures/-NValues in $NSweepDir - nothing to sweep."
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
        "/origami/stlFile geometry_N_sweep/$($si.Name)"
        "/origami/kaptonStlFile geometry_N_sweep/$kaptonName"
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
SayStep "Done: $nRun run(s) completed, results in Results_N_Sweep\"
