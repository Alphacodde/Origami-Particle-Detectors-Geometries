#!/usr/bin/env python3
"""
item4_fiducial.py - Fiducial vs. angular-coverage acceptance decomposition.

Evaluates hit acceptance partitioned by the true primary emission polar angle theta
(labThetaDeg recorded at generation time in the Geant4 ROOT ntuples) across all
five independent simulation runs per geometry.

Fiducial definition:
  |eta| <= 0.80  <=>  theta in [48.392 deg, 131.608 deg]
Peripheral definition:
  |eta| > 0.80   <=>  theta outside [48.392 deg, 131.608 deg]

Key improvements over legacy hit-coordinate proxy:
  1. Bug A fix: True measured denominator (N_gen_fid and N_gen_periph) evaluated
     per run from all 500,000 generated primary tracks, eliminating theoretical
     denominator mismatch and guaranteeing A <= 1.0 by construction.
  2. Bug B fix: Primary emission angle taken directly from labThetaDeg (generated
     polar angle), avoiding distortions from facet displacement at cylinder seams.
"""

import glob
import math
import numpy as np
import uproot

ETA_MAX = 0.80
THETA_MIN_DEG = math.degrees(2.0 * math.atan(math.exp(-ETA_MAX)))  # ~48.392 deg
THETA_MAX_DEG = 180.0 - THETA_MIN_DEG                              # ~131.608 deg

GEOMETRIES = ["barrel_reference", "miura", "kresling", "yoshimura"]
TAG_MAP = {
    "barrel_reference": "scan_origami_barrel_reference_fold1p00_p5p00GeV_vtx3p0mm_run0*.root",
    "miura":            "scan_origami_miura_deployed_fold1p00_p5p00GeV_vtx3p0mm_run0*.root",
    "kresling":         "scan_origami_kresling_deployed_fold1p00_p5p00GeV_vtx3p0mm_run0*.root",
    "yoshimura":        "scan_origami_yoshimura_deployed_fold1p00_p5p00GeV_vtx3p0mm_run0*.root",
}


def analyze_runs(data_dir="results_diff_geom"):
    results = {}
    for g in GEOMETRIES:
        pattern = f"{data_dir}/{TAG_MAP[g]}"
        files = sorted(glob.glob(pattern))
        if not files:
            raise FileNotFoundError(f"No files found for pattern: {pattern}")

        runs_data = []
        for fpath in files:
            with uproot.open(fpath) as f:
                t = f["PionEvents"]
                theta = t["labThetaDeg"].array(library="np")
                hit = t["hitDetector"].array(library="np")

                fid_mask = (theta >= THETA_MIN_DEG) & (theta <= THETA_MAX_DEG)
                n_gen_fid = int(np.sum(fid_mask))
                n_gen_periph = int(len(theta) - n_gen_fid)

                n_hit_fid = int(np.sum((hit == 1) & fid_mask))
                n_hit_periph = int(np.sum((hit == 1) & (~fid_mask)))
                n_hit_global = int(np.sum(hit == 1))

                runs_data.append({
                    "n_total": len(theta),
                    "n_gen_fid": n_gen_fid,
                    "n_gen_periph": n_gen_periph,
                    "n_hit_global": n_hit_global,
                    "n_hit_fid": n_hit_fid,
                    "n_hit_periph": n_hit_periph,
                    "a_global": n_hit_global / len(theta),
                    "a_fid": n_hit_fid / n_gen_fid,
                    "a_periph": n_hit_periph / n_gen_periph,
                })
        results[g] = runs_data
    return results


def print_summary(results):
    print("=" * 95)
    print(f"FIDUCIAL vs. PERIPHERAL ACCEPTANCE DECOMPOSITION (|eta| <= {ETA_MAX:.2f}, theta in [{THETA_MIN_DEG:.2f} deg, {THETA_MAX_DEG:.2f} deg])")
    print("=" * 95)

    b_runs = results["barrel_reference"]
    b_glob = np.mean([r["a_global"] for r in b_runs])
    b_fid = np.mean([r["a_fid"] for r in b_runs])
    b_per = np.mean([r["a_periph"] for r in b_runs])
    b_hits = np.mean([r["n_hit_global"] for r in b_runs])
    b_hits_fid = np.mean([r["n_hit_fid"] for r in b_runs])
    b_hits_per = np.mean([r["n_hit_periph"] for r in b_runs])

    print(f"{'Geometry':<18} {'A_global':>9} {'dA_glob':>9} {'A_fid':>9} {'dA_fid':>9} {'A_periph':>9} {'dA_periph':>10} {'Periph Share':>14}")
    print("-" * 95)

    for g in GEOMETRIES:
        runs = results[g]
        mean_glob = np.mean([r["a_global"] for r in runs])
        mean_fid = np.mean([r["a_fid"] for r in runs])
        mean_per = np.mean([r["a_periph"] for r in runs])
        mean_hits = np.mean([r["n_hit_global"] for r in runs])
        mean_hits_fid = np.mean([r["n_hit_fid"] for r in runs])
        mean_gen_per = np.mean([r["n_gen_periph"] for r in runs])
        mean_tot = np.mean([r["n_total"] for r in runs])

        if g == "barrel_reference":
            print(f"{'Barrel (fine)':<18} {mean_glob:>9.4f} {'---':>9} {mean_fid:>9.4f} {'---':>9} {mean_per:>9.4f} {'---':>10} {'---':>14}")
        else:
            d_glob = (mean_glob - b_glob) * 100
            d_fid = (mean_fid - b_fid) * 100
            d_per = (mean_per - b_per) * 100
            f_periph = mean_gen_per / mean_tot
            periph_share = (d_per * f_periph / d_glob) * 100
            name = "Miura-ori" if g == "miura" else g.capitalize()
            print(f"{name:<18} {mean_glob:>9.4f} {d_glob:>+8.2f}p {mean_fid:>9.4f} {d_fid:>+8.2f}p {mean_per:>9.4f} {d_per:>+9.2f}p {periph_share:>13.1f}%")

    print("=" * 95)


if __name__ == "__main__":
    res = analyze_runs()
    print_summary(res)
