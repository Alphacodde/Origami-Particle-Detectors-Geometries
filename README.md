# Origami-Folded Sensor Geometries for Cylindrical Tracking Detectors: A Geant4 Comparison

[![Paper](https://img.shields.io/badge/JINST-Submitted-blue.svg)](paper/main.pdf)
[![Geant4](https://img.shields.io/badge/Geant4-10.7%2B-red.svg)](https://geant4.web.cern.ch/)
[![ROOT](https://img.shields.io/badge/CERN_ROOT-6.24%2B-orange.svg)](https://root.cern/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-green.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Official source code, parametric geometry generators, analysis pipeline, and manuscript source for the study:  
**"Origami-Folded Sensor Geometries for Cylindrical Tracking Detectors: A Geant4 Comparison"** (Kanishk Sharma, submitted to *JINST*).

---

## 📌 Highlights

- **Three Origami Topologies**: Implements **Yoshimura** (diamond frustum), **Miura-ori** (parallelogram cylinder), and **Kresling** (triangular chiral cylinder) against a standard flat-tile **Barrel reference**.
- **Material-Differentiated Meshes**: Bi-layer GDML models separating active silicon sensor facets ($300\ \mu\text{m}$ / $50\ \mu\text{m}$) and flexible Kapton substrate ($50\ \mu\text{m}$).
- **3-Objective Pareto Frontier**: Quantitative trade-off optimization between **Geometric Acceptance**, **Mean Radiation Length $\langle X/X_0 \rangle$**, and **Dead-Zone Fraction**.
- **Highland Multiple Scattering & Vertex Robustness**: Validated against analytical Highland predictions (Exp 2) and realistic 3D vertex smearing uncertainties (Exp 3).

---

## 📂 Raw Simulation Dataset (Google Drive)

Due to file size limits (~10.7 GB total), the raw Geant4 ROOT output files are hosted on Google Drive:

[![Google Drive Dataset](https://img.shields.io/badge/Google%20Drive-Download%20ROOT%20Ntuples%20(10.7%20GB)-4285F4?style=for-the-badge&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1ZbLm0H2FCLPQepytwLGCHj7_VS4sdF78?usp=sharing)

> **Direct Link**: [https://drive.google.com/drive/folders/1ZbLm0H2FCLPQepytwLGCHj7_VS4sdF78?usp=sharing](https://drive.google.com/drive/folders/1ZbLm0H2FCLPQepytwLGCHj7_VS4sdF78?usp=sharing)

### Folder Setup
Download and place the unzipped folders at the root of this repository:

```text
Paper Code/
├── results_diff_geom/      # Nominal 300 µm multi-run ground truth (Exp 1 + Dead-zone)
├── results_exp2/           # Exp 2: Multiple Coulomb scattering momentum scan
├── results_exp3/           # Exp 3: Vertex-smear Gaussian robustness scan
├── results_50um/           # 50 µm thin-silicon ITS3 comparison
├── Results_N_Sweep/        # Facet-count (N ∈ {6..20}) parameter sweep
└── Results_Theta_Sweep/    # Fold-angle (θ ∈ {45°..75°}) parameter sweep
```

---

## 🚀 Quickstart

### 1. Environment Setup
```bash
# Clone the repository
git clone https://github.com/<your-username>/origami-tracking-detectors.git
cd origami-tracking-detectors

# Create conda environment (recommended)
conda env create -f environment.yml
conda activate origami-tracker

# Or via pip
pip install -r requirements.txt
```

### 2. Reproduce All Paper Figures
Once the ROOT datasets are placed in the repository root:
```bash
python analysis/reproduce_all_figures.py --data-root .
```
All 27 publication-grade figures (300 DPI) will be regenerated in `paper/figures/`.

### 3. Compile Manuscript
```bash
cd paper
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```
Output: `paper/main.pdf` (32 pages, camera-ready).

---

## 🏗️ Repository Architecture

```text
Paper Code/
├── paper/                     # LaTeX manuscript source and compiled PDF
│   ├── main.tex               # JINST-format manuscript
│   ├── biblio.bib             # BibTeX bibliography (10 references)
│   ├── jinst-pub.sty          # JINST style package
│   ├── main.pdf               # Camera-ready compiled PDF (32 pages)
│   └── figures/               # 27 publication figures (300 DPI)
│
├── geant4_sim/                # Geant4 C++ simulation application
│   ├── CMakeLists.txt         # Build configuration (Geant4 10.7+, C++17)
│   ├── origamiDet.cc          # Application entry point
│   ├── include/               # 8 C++ header files
│   ├── src/                   # 7 C++ implementation files
│   ├── macros/                # Run macros (run_normal.mac, run_exp2.mac, etc.)
│   └── CADMesh/               # Single-header STL/OBJ mesh loader
│
├── geometry_generation/       # Parametric CAD generator library
│   ├── generate_all_geometries.py # Master one-click GDML batch builder
│   └── diff_geom_macros/      # Geometry tessellation modules
│       ├── _differentiated_mesh.py # Bi-layer Si/Kapton faceted mesh builder
│       ├── _gdml_export.py    # Geant4 GDML exporter
│       ├── _geometry_config.py# Geometry dimensions & envelope parameters
│       ├── comparison_barrel.py
│       ├── miura_cylindrical.py
│       ├── kresling.py
│       ├── yoshimura.py
│       ├── sweep_geometry_by_N.py
│       └── sweep_geometry_by_theta.py
│
├── geometries/                # Pre-generated GDML models (Ready-to-simulate)
│   ├── nominal_300um/         # Nominal 300 µm silicon GDML files (8 files)
│   ├── thin_50um/             # 50 µm thin-silicon GDML files (8 files)
│   └── validation/            # Tier-1 flat-plate validation geometry
│
├── pipelines/                 # Automated end-to-end PowerShell runners
│   ├── run_everything_diff_mesh.ps1 # Master pipeline for nominal 300 µm runs
│   ├── run_everything_50um.ps1      # 50 µm ITS3 comparison pipeline
│   ├── run_everything_N.ps1         # Facet-count (N) sweep pipeline
│   ├── run_everything_theta.ps1     # Fold-angle (θ) sweep pipeline
│   ├── run_exp2.ps1                 # Experiment 2 MCS momentum scan
│   ├── run_exp3.ps1                 # Experiment 3 vertex smearing sweep
│   ├── deadzone_analysis.ps1        # Dead-zone extraction pipeline
│   └── deadzone_threshold_sweep.ps1 # Threshold sensitivity sweep
│
├── analysis/                  # Python analysis and visualization suite
│   ├── reproduce_all_figures.py     # Master one-click plot generator
│   ├── pareto.py                    # 3-objective Pareto & trade-off analysis
│   ├── scatter_fit.py               # Highland MCS Rayleigh+Rutherford MLE fit
│   ├── deadzone_map.py              # 2D (θ, z) spatial dead-zone mapping
│   ├── deadzone_threshold_sweep.py  # Response surface & threshold robustness
│   ├── analyze_exp2.py              # Multiple scattering momentum analysis
│   ├── analyze_exp3.py              # Vertex smear Gaussian robustness
│   ├── make_its3_figure.py          # ALICE ITS3 thin-Si comparison plot
│   ├── analysis_m3_differentiated.py# Primary ntuple hit reader
│   └── extract_hits_for_deadzone.py # ROOT-to-CSV hit table extractor
│
└── tools/                     # Geometry validation and calibration tools
    ├── prep_stl.py                  # STL watertightness and normal validator
    ├── thicken_stack.py             # Volumetric solid mesh generator
    ├── make_validation_plate.py     # Flat-plate analytical benchmark builder
    └── check_validation.py         # Statistical validation gate vs. theory
```

---

## 🔬 Physics Specifications & Envelope

All simulated geometries adhere to identical global boundary envelopes:

| Parameter | Nominal Specification | Thin-Silicon Variant |
|:---|:---|:---|
| **Cylinder Envelope Radius ($R$)** | $50.0\text{ mm}$ | $50.0\text{ mm}$ |
| **Active Length ($L$)** | $200.0\text{ mm}$ | $200.0\text{ mm}$ |
| **Silicon Thickness ($t_{\text{Si}}$)** | $300\ \mu\text{m}$ | $50\ \mu\text{m}$ (ALICE ITS3 matched) |
| **Kapton Substrate Thickness** | $50\ \mu\text{m}$ | $50\ \mu\text{m}$ |
| **Primary Particle Beam** | $\pi^+$ pions ($500\text{ MeV}/c$) | $\pi^+$ pions ($500\text{ MeV}/c$) |
| **Nominal Vertex Smear** | 3D Gaussian ($\sigma = 1.0\text{ mm}$) | 3D Gaussian ($\sigma = 1.0\text{ mm}$) |
| **Events per Run** | $500{,}000$ primaries | $500{,}000$ primaries |

---

## ⚙️ Building & Running Geant4 Simulations

### 1. Compilation
```bash
mkdir build && cd build
cmake ../geant4_sim -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### 2. Standalone Simulation Run
Pre-generated GDML models in `geometries/` allow immediate execution without Python CAD generation:
```bash
./origamiDet ../geant4_sim/macros/run_normal.mac
```

### 3. Parametric Geometry Regeneration
To generate customized STL/GDML models with different fold parameters:
```bash
# Generate all baseline models
python geometry_generation/generate_all_geometries.py
```

---

## 📊 Summary of Experiments & Paper Figures

| Experiment / Study | Pipeline | Primary Output Figures | Paper Section |
|:---|:---|:---|:---|
| **Exp 1: Baseline Comparison** | `pipelines/run_everything_diff_mesh.ps1` | `fig6a`, `fig6b`, `fig7`, `fig8`, `fig9` | Sections 4.1, 4.2 |
| **Exp 2: Multiple Scattering (MCS)** | `pipelines/run_exp2.ps1` | `fig_scattering_mixture_*` (6 panels) | Section 4.5 |
| **Exp 3: Vertex Smear Robustness** | `pipelines/run_exp3.ps1` | `fig_vertex_smear_sigma`, `fig_vertex_smear_displacement` | Section 4.6 |
| **Facet Count ($N$) Sweep** | `pipelines/run_everything_N.ps1` | `fig4_n_sweep` | Section 4.3 |
| **Fold Angle ($\theta$) Sweep** | `pipelines/run_everything_theta.ps1` | `fig5_theta_sweep` | Section 4.4 |
| **Dead-Zone Mapping & Sweep** | `pipelines/deadzone_analysis.ps1` | `fig10_map`, `fig11_map`, `fig12_map`, `fig_deadzone_threshold_sweep` | Section 4.7 |
| **Thin-Silicon ($50\ \mu\text{m}$) ITS3** | `pipelines/run_everything_50um.ps1` | `fig_its3_comparison_v2` | Section 4.8 |

---

## 📖 Citation

If you use this simulation framework, geometry generator, or dataset, please cite:

```bibtex
@article{sharma2025origami,
  title   = {Origami-Folded Sensor Geometries for Cylindrical Tracking Detectors: A {Geant4} Comparison},
  author  = {Sharma, Kanishk},
  journal = {JINST},
  year    = {2025},
  note    = {submitted}
}
```

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
