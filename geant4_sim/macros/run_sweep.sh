#!/bin/bash
# run_sweep.sh - runs origamiDet across all 5 geometries automatically,
# generating a temp macro per geometry and renaming the ROOT output so
# results don't collide (see the naming-collision warning in
# template_run.mac).
#
# Usage:
#   ./run_sweep.sh
#
# Edit the GEOMETRIES array below to match your actual STL filenames once
# you've run prep_stl.py on them and confirmed they're valid.
#
# Run this from the GEANT4 BUILD directory (where origamiDet lives and
# where CMake copied the geometry/ and macros/ folders).

set -e  # stop on first error - don't silently continue with bad data if one
        # geometry's simulation fails

mkdir -p results

# --- Edit this list to match your actual 5 STL files ---
GEOMETRIES=(
  "miura_deployed"
  "kresling_deployed"
  "waterbomb_deployed"
  "kirigami_deployed"
  "fifth_geometry_deployed"
)

N_EVENTS=500000
MOMENTUM_GEV=4.0
MAX_ANGLE_DEG=60
SENSOR_THICKNESS_MM=0.15

for geom in "${GEOMETRIES[@]}"; do
  echo "=================================================="
  echo "Running: $geom"
  echo "=================================================="

  stl_path="geometry/${geom}.stl"
  if [ ! -f "$stl_path" ]; then
    echo "WARNING: $stl_path not found, skipping. Check the GEOMETRIES list"
    echo "at the top of this script matches your actual STL filenames in"
    echo "the geometry/ directory."
    continue
  fi

  macro_path="macros/_generated_${geom}.mac"
  cat > "$macro_path" <<EOF
/origami/stlFile ${stl_path}
/origami/sensorThickness ${SENSOR_THICKNESS_MM} mm
/origami/foldState 1.0

/origami/gun/randomAngleMaxDeg ${MAX_ANGLE_DEG}
/origami/gun/momentumGeV ${MOMENTUM_GEV}
/origami/gun/gunZmm 500

/run/verbose 0
/event/verbose 0
/tracking/verbose 0

/run/initialize
/run/beamOn ${N_EVENTS}
EOF

  ./origamiDet "$macro_path"

  # RunAction names output origami_run0.root every time (run ID resets per
  # process) - rename immediately so the next geometry doesn't overwrite it.
  if [ -f "origami_run0.root" ]; then
    mv origami_run0.root "results/${geom}.root"
    echo "-> results/${geom}.root"
  else
    echo "ERROR: expected output origami_run0.root not found after run for $geom"
    echo "Check GEANT4 stdout above for errors before trusting other results."
  fi
done

echo ""
echo "Sweep complete. Results in results/*.root"
echo "Next: run analysis/analyze_results.py to build the Pareto frontier."
