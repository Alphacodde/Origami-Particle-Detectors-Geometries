import sys
import glob
import os
import numpy as np
import uproot
import json
SILICON_X0_MM = 93.66
KAPTON_X0_MM = 285.75
NTUPLE_NAME = 'PionEvents'

def load_run(root_path):
    with uproot.open(root_path) as f:
        available = [k.split(';')[0] for k in f.keys()]
        name = NTUPLE_NAME if NTUPLE_NAME in available else available[0]
        if name != NTUPLE_NAME:
            print(f"  WARNING: '{NTUPLE_NAME}' not found in {root_path}; falling back to '{name}' (available: {available})")
        tree = f[name]
        data = tree.arrays(library='np')
    return data

def summarize_geometry(name, data):
    hit = data['hitDetector'].astype(bool)
    n_total = len(hit)
    n_hit = int(hit.sum())
    acceptance_fraction = n_hit / n_total if n_total > 0 else 0.0
    has_split = 'siliconPathLength_mm' in data and 'kaptonPathLength_mm' in data
    if has_split:
        si_mm = data['siliconPathLength_mm'][hit]
        kap_mm = data['kaptonPathLength_mm'][hit]
        path_X0 = si_mm / SILICON_X0_MM + kap_mm / KAPTON_X0_MM
        path_mm = si_mm + kap_mm
        si_zero_frac = float((si_mm == 0).mean()) if n_hit else None
        kap_zero_frac = float((kap_mm == 0).mean()) if n_hit else None
    else:
        path_mm = data['totalPathLength_mm'][hit]
        path_X0 = path_mm / SILICON_X0_MM
        si_zero_frac = None
        kap_zero_frac = None
    edep_MeV = data['totalEdep_MeV'][hit]
    lab_theta_all = data.get('labThetaDeg', None)
    local_inc_all = data.get('localIncidenceDeg', None)
    hit_lab_theta = lab_theta_all[hit] if lab_theta_all is not None else None
    if local_inc_all is not None:
        hit_local_inc = local_inc_all[hit]
        valid_local_mask = hit_local_inc >= 0.0
        hit_local_inc_valid = hit_local_inc[valid_local_mask]
        path_X0_valid_local = path_X0[valid_local_mask]
    else:
        hit_local_inc_valid = None
        path_X0_valid_local = None
    lab_theta_binned = _bin_by_lab_theta(lab_theta_all, hit, path_X0) if lab_theta_all is not None else None
    local_inc_binned = _bin_by_local_incidence(hit_local_inc_valid, path_X0_valid_local) if hit_local_inc_valid is not None and len(hit_local_inc_valid) > 0 else None
    if local_inc_binned is not None:
        legacy_angle_binned = local_inc_binned
    elif lab_theta_binned is not None:
        legacy_angle_binned = lab_theta_binned
    else:
        legacy_angle_binned = None
    result = {'geometry': name, 'n_events_total': n_total, 'n_events_hit': n_hit, 'acceptance_fraction': acceptance_fraction, 'differentiated_scoring_used': has_split, 'mean_path_length_mm': float(np.mean(path_mm)) if n_hit else None, 'mean_path_length_X0': float(np.mean(path_X0)) if n_hit else None, 'median_path_length_X0': float(np.median(path_X0)) if n_hit else None, 'std_path_length_X0': float(np.std(path_X0, ddof=1)) if n_hit > 1 else None, 'sem_path_length_X0': float(np.std(path_X0, ddof=1) / np.sqrt(n_hit)) if n_hit > 1 else None, 'p25_path_length_X0': float(np.percentile(path_X0, 25)) if n_hit else None, 'p75_path_length_X0': float(np.percentile(path_X0, 75)) if n_hit else None, 'mean_edep_MeV': float(np.mean(edep_MeV)) if n_hit else None, 'mean_lab_theta_deg': float(np.mean(hit_lab_theta)) if hit_lab_theta is not None and len(hit_lab_theta) > 0 else None, 'mean_local_incidence_deg': float(np.mean(hit_local_inc_valid)) if hit_local_inc_valid is not None and len(hit_local_inc_valid) > 0 else None, 'kapton_zero_path_fraction': kap_zero_frac, 'silicon_zero_path_fraction': si_zero_frac, 'path_X0_by_lab_theta': lab_theta_binned, 'path_X0_by_local_incidence': local_inc_binned, 'path_X0_by_angle_bin': legacy_angle_binned}
    return result

def _bin_by_lab_theta(lab_theta_all, hit_mask, path_X0_hit, n_bins=18, max_angle=180):
    bins = np.linspace(0, max_angle, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    means = []
    acceptance = []
    hit_theta = lab_theta_all[hit_mask]
    for lo, hi in zip(bins[:-1], bins[1:]):
        in_bin_all = (lab_theta_all >= lo) & (lab_theta_all < hi) if hi < max_angle else (lab_theta_all >= lo) & (lab_theta_all <= hi)
        n_in_bin_all = int(in_bin_all.sum())
        n_in_bin_hit = int((in_bin_all & hit_mask).sum())
        acc = n_in_bin_hit / n_in_bin_all if n_in_bin_all > 0 else 0.0
        acceptance.append(float(acc))
        in_bin_hit = (hit_theta >= lo) & (hit_theta < hi) if hi < max_angle else (hit_theta >= lo) & (hit_theta <= hi)
        if in_bin_hit.sum() > 0:
            means.append(float(np.mean(path_X0_hit[in_bin_hit])))
        else:
            means.append(None)
    return {'bin_centers_deg': bin_centers.tolist(), 'mean_path_X0': means, 'acceptance_fraction': acceptance}

def _bin_by_local_incidence(local_inc_valid, path_X0_valid, n_bins=9, max_angle=90):
    bins = np.linspace(0, max_angle, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    means = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (local_inc_valid >= lo) & (local_inc_valid < hi) if hi < max_angle else (local_inc_valid >= lo) & (local_inc_valid <= hi)
        means.append(float(np.mean(path_X0_valid[mask])) if mask.sum() > 0 else None)
    return {'bin_centers_deg': bin_centers.tolist(), 'mean_path_X0': means}
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python3 analysis_m3_differentiated.py results/*.root')
        sys.exit(1)
    root_files = []
    for pattern in sys.argv[1:]:
        matches = sorted(glob.glob(pattern))
        if matches:
            root_files.extend(matches)
        elif os.path.exists(pattern):
            root_files.append(pattern)
        else:
            print(f"WARNING: no files matched '{pattern}' - skipping")
    if not root_files:
        print('ERROR: no .root files found. Check the path/pattern you passed in.')
        sys.exit(1)
    print(f'Found {len(root_files)} ROOT file(s): {root_files}')
    all_summaries = []
    for path in root_files:
        name = os.path.basename(path).replace('.root', '')
        print(f"\nLoading {path} as geometry '{name}'...")
        data = load_run(path)
        summary = summarize_geometry(name, data)
        all_summaries.append(summary)
        print(json.dumps({k: v for k, v in summary.items() if not k.startswith('path_X0_by')}, indent=2))
    out_path = 'geant4_run_Theta.json'
    with open(out_path, 'w') as f:
        json.dump(all_summaries, f, indent=2)
    print(f'\nWrote {out_path}')