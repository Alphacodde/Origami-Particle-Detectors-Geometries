# Origami-Folded Sensor Geometries for Cylindrical Tracking Detectors: A Geant4 Comparison

[![Geant4](https://img.shields.io/badge/Geant4-10.7%2B-red.svg)](https://geant4.web.cern.ch/)
[![ROOT](https://img.shields.io/badge/CERN_ROOT-6.24%2B-orange.svg)](https://root.cern/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-green.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Official source code, parametric geometry generators, analysis pipeline, and manuscript source for the study:  
**"Origami-Folded Sensor Geometries for Cylindrical Tracking Detectors: A Geant4 Comparison"**  
*Kanishk Sharma*

---

## 🌟 Highlights

- **Three Cylindrical Origami Topologies**: Implements **Cylindrical Miura-ori** (parallelogram tessellation), **Yoshimura** (diamond frustum), and **Kresling** (chiral triangular tower) against a matched fine-tiled **Barrel reference** ($N_{\text{cols}} = 141$).
- **Material-Differentiated Meshes**: Bi-layer GDML models separating active discrete silicon sensor facets ($300\ \mu\text{m}$ nominal / $50\ \mu\text{m}$ thin) and continuous flexible Kapton substrate ($50\ \mu\text{m}$) with normal-projected seam clearance.
- **Pareto Trade-Off Frontier**: Evaluates the multi-objective balance between **Geometric Hit Acceptance ($A$)**, **Traversed Radiation Length ($\langle X/X_0 \rangle$)**, and **Relative Under-Coverage Fraction ($\delta_{\text{under}}$)**.
- **Fiducial vs. Peripheral Decomposition**: True per-run primary polar angle decomposition ($|\eta| \le 0.80$ vs $|\eta| > 0.80$) revealing that acceptance gains are heavily driven by peripheral angular coverage (94.4%–98.9% peripheral share).
- **Multiple Scattering & Robustness Validation**: Core-plus-tail mixture model validation against the analytical Highland formula across $p \in [0.5, 10.0]\text{ GeV}/c$ (Exp 2) and 3D Gaussian vertex smearing stability up to $\sigma = 10\text{ mm}$ (Exp 3).

---

## 📊 Key Results & Headline Performance

Across five independent Monte Carlo replications ($500{,}000$ primary $5\text{ GeV}/c\ \pi^+$ events per run, 3D vertex smear $\sigma = 3.0\text{ mm}$), all candidate fold geometries demonstrate distinct Pareto-optimal trade-offs:

### 1. Headline Acceptance & Material Traversal (Table 4 & Table 5)
| Geometry | Hit Acceptance $A$ | $\Delta A$ vs. Barrel | Cohen's $h$ | Mean $X/X_0$ | $\Delta(X/X_0)$ | Under-Coverage ($R < 0.5$) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **Barrel (fine, $N=141$)** | $83.17 \pm 0.05\%$ | Baseline | — | $0.003977 \pm 10^{-6}$ | Baseline | $0.00\%$ (Reference) |
| **Miura-ori** | $83.65 \pm 0.05\%$ | $+0.49\text{ pp}$ ($p < 10^{-7}$) | $+0.013$ | $0.004051 \pm 10^{-6}$ | $+1.88\%$ | **$0.04\%$** (Near-uniform) |
| **Yoshimura** | $84.49 \pm 0.04\%$ | $+1.30\text{ pp}$ ($p < 10^{-9}$) | $+0.035$ | $0.004295 \pm 10^{-6}$ | $+8.01\%$ | **$1.86\%$** (Localized seams) |
| **Kresling** | $85.52 \pm 0.04\%$ | **$+2.36\text{ pp}$** ($p < 10^{-10}$) | $+0.065$ | $0.004203 \pm 10^{-6}$ | $+5.69\%$ | **$1.45\%$** (Helical seams) |

*Note: Flat origami plate geometry is permanently excluded from all Pareto and detector trade-off analyses.*

### 2. Angular Acceptance Decomposition (Table 11, Section 4.7)
Evaluating acceptance partitioned by generated particle polar angle into fiducial ($|\eta| \le 0.80$, $N_{\text{gen}} \approx 332{,}000$) and peripheral ($|\eta| > 0.80$, $N_{\text{gen}} \approx 168{,}000$) regions:

| Fold Geometry | $\Delta A_{\text{global}}$ | $\Delta A_{\text{fiducial}}$ | $z_{\text{fid}}$ | $\Delta A_{\text{peripheral}}$ | Peripheral Share of Gain |
|:---|:---:|:---:|:---:|:---:|:---:|
| **Miura-ori** | $+0.49\text{ pp}$ | $+0.04\text{ pp}$ | $+0.46$ (n.s.) | $+1.39\text{ pp}$ | **$94.4\%$** |
| **Kresling** | $+2.36\text{ pp}$ | $-0.09\text{ pp}$ | $-0.95$ (n.s.) | $+7.20\text{ pp}$ | **$97.5\%$** |
| **Yoshimura** | $+1.30\text{ pp}$ | $-0.04\text{ pp}$ | $-0.46$ (n.s.) | $+3.96\text{ pp}$ | **$98.9\%$** |

All three fold families retain near-identical fiducial acceptance to the flat barrel control ($\Delta A_{\text{fid}} \approx 0.00\text{ pp}$), with over 94% of their geometric gain originating in the peripheral pseudorapidity regions.

---

## 📦 Raw Simulation Dataset (Google Drive)

Due to file size limits (> 10 GB), the complete Geant4 ROOT output ntuples are hosted on Google Drive:

[![Google Drive Dataset](https://img.shields.io/badge/Google%20Drive-Download%20ROOT%20Ntuples%20(10.7%20GB)-4285F4?style=for-the-badge&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1ZbLm0H2FCLPQepytwLGCHj7_VS4sdF78?usp=sharing)

> **Direct Link**: [https://drive.google.com/drive/folders/1ZbLm0H2FCLPQepytwLGCHj7_VS4sdF78?usp=sharing](https://drive.google.com/drive/folders/1ZbLm0H2FCLPQepytwLGCHj7_VS4sdF78?usp=sharing)

### Folder Setup
Download and extract the dataset folders into the repository root:

```text
.
├── results_diff_geom/       # Nominal 300 µm multi-run ground truth (Exp 1 & dead-zone mapping)
├── results_exp2/            # Exp 2: Multiple Coulomb scattering momentum scan (0.5 - 10 GeV/c)
├── results_exp3/            # Exp 3: Vertex-smear Gaussian robustness scan (sigma 0 - 20 mm)
├── results_50um/            # 50 µm thin-silicon ALICE ITS3 comparison
├── Results_N_Sweep/         # Facet-count (N in {6..32}) parameter sweep
└── Results_Theta_Sweep/     # Fold-angle (theta) parameter sweep
```

---

## 🚀 Quickstart

### 1. Environment Setup
```bash
# Clone repository
git clone https://github.com/Alphacodde/Origami-Particle-Detectors-Geometries.git
cd Origami-Particle-Detectors-Geometries

# Create conda environment (recommended)
conda env create -f environment.yml
conda activate origami-tracker

# Or via pip
pip install -r requirements.txt
```

### 2. Reproduce Figures and Analysis
Once the ROOT datasets are placed in the repository root:
```bash
# Generate all publication figures (300 DPI)
python analysis/reproduce_all_figures.py --data-root .

# Run fiducial vs peripheral acceptance decomposition (Table 11)
python analysis/fiducial_acceptance.py

# Generate unrolled spatial dead-zone / hit-density ratio maps
python analysis/deadzone_map.py original results_diff_geom results_diff_geom
```

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

## 📁 Repository Architecture

```text
Origami-Particle-Detectors-Geometries/
├── paper/                          # LaTeX manuscript source and compiled camera-ready PDF
│   ├── main.tex                    # JINST-format camera-ready manuscript (32 pages)
│   ├── main.pdf                    # Compiled PDF document
│   ├── biblio.bib                  # BibTeX bibliography
│   ├── JHEP.bst                    # Bibliography style file
│   ├── jinstpub.sty                # JINST style package
│   └── figures/                    # Publication figures (300 DPI PNGs)
│
├── analysis/                       # Python analysis, statistical tests & plotting suite
│   ├── fiducial_acceptance.py      # Fiducial vs peripheral decomposition (Table 11)
│   ├── deadzone_map.py             # 2D (theta, z) spatial hit-density ratio mapping
│   ├── deadzone_threshold_sweep.py # Threshold sensitivity and response surfaces
│   ├── pareto.py                   # 3-objective Pareto optimization & trade-off analysis
│   ├── scatter_fit.py              # Highland MCS Rayleigh+PowerLaw mixture MLE fit
│   ├── analyze_exp2.py             # Multiple scattering momentum analysis (Exp 2)
│   ├── analyze_exp3.py             # Vertex smear Gaussian robustness analysis (Exp 3)
│   ├── make_its3_figure.py         # ALICE ITS3 thin-silicon comparison
│   ├── reproduce_all_figures.py    # Master figure generation script
│   ├── rerender_high_res.py        # 300-DPI figure re-renderer
│   ├── analysis_m3_differentiated.py # Differentiated mesh hit reader
│   └── extract_hits_for_deadzone.py  # ROOT-to-CSV hit table extractor
│
├── geant4_sim/                     # Geant4 C++ simulation application
│   ├── CMakeLists.txt              # Build configuration (Geant4 10.7+, C++17)
│   ├── origamiDet.cc               # Application entry point
│   ├── include/                    # C++ header files (DetectorConstruction, EventAction, etc.)
│   ├── src/                        # C++ implementation files
│   ├── macros/                     # Geant4 execution macros
│   └── CADMesh/                    # Tessellated mesh loader
│
├── geometry_generation/            # Parametric CAD generator library
│   ├── generate_all_geometries.py  # Master GDML generation script
│   └── diff_geom_macros/           # Tessellation and seam clearance generators
│       ├── _geometry_config.py     # Envelope dimensions (R=40mm, L=120mm)
│       ├── _differentiated_mesh.py # Bi-layer Si/Kapton faceted mesh builder
│       ├── _gdml_export.py         # GDML tessellated solid exporter
│       ├── comparison_barrel.py    # Fine-tiled cylindrical control generator
│       ├── kresling.py             # Kresling origami cylinder generator
│       ├── miura_cylindrical.py    # Cylindrical Miura-ori generator
│       ├── yoshimura.py            # Yoshimura diamond pattern generator
│       ├── sweep_geometry_by_N.py  # Facet count (N) sweep generator
│       └── sweep_geometry_by_theta.py # Fold angle (theta) sweep generator
│
├── geometries/                     # Pre-generated GDML models (ready-to-simulate)
│   ├── nominal_300um/              # Nominal 300 µm silicon GDML models
│   ├── thin_50um/                  # 50 µm thin-silicon GDML models
│   └── validation/                 # Tier-1 flat-plate validation geometry
│
├── pipelines/                      # Automated PowerShell simulation runners
│   ├── run_everything_diff_mesh.ps1# Nominal 300 µm multi-run replication pipeline
│   ├── run_everything_50um.ps1     # 50 µm thin-silicon ITS3 comparison pipeline
│   ├── run_everything_N.ps1        # Facet-count (N) parameter sweep pipeline
│   ├── run_everything_theta.ps1    # Fold-angle (theta) parameter sweep pipeline
│   ├── run_exp2.ps1                # Experiment 2 MCS momentum scan
│   ├── run_exp3.ps1                # Experiment 3 vertex smear sweep
│   └── deadzone_analysis.ps1       # Dead-zone extraction pipeline
│
└── tools/                          # Verification, QA, and packaging tools
    ├── prep_stl.py                 # STL watertightness and normal orientation QA
    ├── thicken_stack.py            # Solid extrusion and clearance validator
    ├── make_validation_plate.py    # Flat-plate analytical benchmark builder
    └── check_validation.py        # Statistical validation gate vs analytical theory
```

---

## 📐 Physics Specifications & Envelope

All simulated geometries adhere to identical global boundary envelopes matching the camera-ready manuscript:

| Parameter | Nominal Specification | Thin-Silicon Benchmark | Notes |
|:---|:---:|:---:|:---|
| **Cylinder Envelope Radius ($R$)** | $40.0\text{ mm}$ | $40.0\text{ mm}$ | Matched across all 4 geometries |
| **Active Length ($L$)** | $120.0\text{ mm}$ | $120.0\text{ mm}$ | Symmetric around origin ($z \in [-60, +60]\text{ mm}$) |
| **Silicon Sensor Thickness ($t_{\text{Si}}$)** | $300\ \mu\text{m}$ | $50\ \mu\text{m}$ | Active detection layer ($X_0 = 93.66\text{ mm}$) |
| **Kapton Substrate Thickness** | $50\ \mu\text{m}$ | $50\ \mu\text{m}$ | Passive flexible substrate ($X_0 = 285.75\text{ mm}$) |
| **Primary Particle Beam** | $\pi^+$ pions ($5.0\text{ GeV}/c$) | $\pi^+$ pions ($5.0\text{ GeV}/c$) | Isotropic $4\pi$ emission |
| **Nominal Vertex Smear** | 3D Gaussian ($\sigma = 3.0\text{ mm}$) | 3D Gaussian ($\sigma = 3.0\text{ mm}$) | Realistic IP beam-spot uncertainty |
| **Events per Replication** | $500{,}000$ primaries | $500{,}000$ primaries | Five independent runs per geometry |
| **Fiducial Cut** | $|\eta| \le 0.80$ | $|\eta| \le 0.80$ | $\theta \in [48.39^\circ, 131.61^\circ]$ |

---

## ⚙️ Building & Running Geant4 Simulations

### 1. Build Simulation
```bash
mkdir build && cd build
cmake ../geant4_sim -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release --parallel
```

### 2. Standalone Simulation Run
Pre-generated GDML models in `geometries/` allow immediate execution:
```bash
./origamiDet ../geant4_sim/macros/template_run.mac
```

### 3. Parametric Geometry Regeneration
To generate customized STL/GDML models with varied fold parameters:
```bash
python geometry_generation/generate_all_geometries.py
```

---

## 📋 Summary of Experiments & Paper Figures

| Experiment / Study | Pipeline | Primary Output Figures | Paper Section |
|:---|:---|:---|:---|
| **Exp 1: Baseline Comparison** | `pipelines/run_everything_diff_mesh.ps1` | `fig6a`, `fig6b`, `fig7`, `fig8`, `fig9` | Sections 4.1, 4.2 |
| **Exp 2: Multiple Scattering (MCS)** | `pipelines/run_exp2.ps1` | `fig_scattering_mixture_*` (6 panels), `fig_scattering_theta0` | Section 4.5 |
| **Exp 3: Vertex Smear Robustness** | `pipelines/run_exp3.ps1` | `fig_vertex_smear_sigma`, `fig_vertex_smear_displacement` | Section 4.6 |
| **Facet Count ($N$) Sweep** | `pipelines/run_everything_N.ps1` | `fig4_n_sweep`, `fig3_yoshimura_degenerate_n6`, `fig_kresling_degenerate` | Section 4.3 |
| **Fold Angle ($\theta$) Sweep** | `pipelines/run_everything_theta.ps1` | `fig5_theta_sweep` | Section 4.4 |
| **Spatial Hit-Density & Dead Zones** | `pipelines/deadzone_analysis.ps1` | `fig10_map_miura`, `fig11_map_yoshimura`, `fig12_map_kresling`, `fig_deadzone_threshold_sweep` | Section 4.7 |
| **Fiducial Decomposition (Table 11)** | `analysis/fiducial_acceptance.py` | Table 11 | Section 4.7 |
| **Thin-Silicon ($50\ \mu\text{m}$) ITS3** | `pipelines/run_everything_50um.ps1` | `fig_its3_comparison_v2` | Section 4.8 |

---

## 📄 Citation

```bibtex
@article{Sharma:2026origami,
  author        = {Sharma, Kanishk},
  title         = {Origami-Folded Sensor Geometries for Cylindrical Tracking Detectors: A {Geant4} Comparison},
  journal       = {Journal of Instrumentation},
  year          = {2026},
  eprint        = {arXiv:26xx.xxxxx},
  archivePrefix = {arXiv},
  primaryClass  = {physics.ins-det}
}
```

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
