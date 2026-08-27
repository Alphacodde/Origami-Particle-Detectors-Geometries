import sys
import glob
import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
import json
SILICON_X0_MM = 93.7

def load_run(root_path):
    with uproot.open(root_path) as f:
        tree = f['MuonEvents']
        data = tree.arrays(library='np')
    return data

def summarize_geometry(name, data):
    hit = data['hitDetector'].astype(bool)
    n_total = len(hit)
    n_hit = int(hit.sum())
    acceptance_fraction = n_hit / n_total if n_total > 0 else 0.0
    path_mm = data['totalPathLength_mm'][hit]
    edep_MeV = data['totalEdep_MeV'][hit]
    path_X0 = path_mm / SILICON_X0_MM
    lab_theta_all = data.get('labThetaDeg', data.get('incidentAngleDeg', None))
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
    elif hit_lab_theta is not None:
        legacy_angle_binned = _bin_by_angle(hit_lab_theta, path_X0)
    else:
        legacy_angle_binned = None
    result = {'geometry': name, 'n_events_total': n_total, 'n_events_hit': n_hit, 'acceptance_fraction': acceptance_fraction, 'mean_path_length_mm': float(np.mean(path_mm)) if n_hit else None, 'mean_path_length_X0': float(np.mean(path_X0)) if n_hit else None, 'median_path_length_X0': float(np.median(path_X0)) if n_hit else None, 'std_path_length_X0': float(np.std(path_X0, ddof=1)) if n_hit > 1 else None, 'sem_path_length_X0': float(np.std(path_X0, ddof=1) / np.sqrt(n_hit)) if n_hit > 1 else None, 'p25_path_length_X0': float(np.percentile(path_X0, 25)) if n_hit else None, 'p75_path_length_X0': float(np.percentile(path_X0, 75)) if n_hit else None, 'mean_edep_MeV': float(np.mean(edep_MeV)) if n_hit else None, 'mean_lab_theta_deg': float(np.mean(hit_lab_theta)) if hit_lab_theta is not None and len(hit_lab_theta) > 0 else None, 'mean_local_incidence_deg': float(np.mean(hit_local_inc_valid)) if hit_local_inc_valid is not None and len(hit_local_inc_valid) > 0 else None, 'path_X0_by_lab_theta': lab_theta_binned, 'path_X0_by_local_incidence': local_inc_binned, 'path_X0_by_angle_bin': legacy_angle_binned}
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

def _bin_by_angle(angle_deg, path_X0, n_bins=6, max_angle=60):
    bins = np.linspace(0, max_angle, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    means = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (angle_deg >= lo) & (angle_deg < hi)
        means.append(float(np.mean(path_X0[mask])) if mask.sum() > 0 else None)
    return {'bin_centers_deg': bin_centers.tolist(), 'mean_path_X0': means}

def build_pareto_frontier(summaries, coverage_metrics):
    fig, ax = plt.subplots(figsize=(8, 6))
    for s in summaries:
        name = s['geometry']
        x0 = s['mean_path_length_X0']
        coverage = coverage_metrics.get(name)
        if x0 is None or coverage is None:
            print(f'WARNING: missing data for {name}, skipping in Pareto plot (x0={x0}, coverage={coverage})')
            continue
        ax.scatter(x0, coverage, s=120, label=name)
        ax.annotate(name, (x0, coverage), xytext=(5, 5), textcoords='offset points', fontsize=9)
    ax.set_xlabel('Mean material budget traversed (X0)')
    ax.set_ylabel('Coverage: area gain ratio (deployed area / stowed footprint)')
    ax.set_title('Coverage vs. Material Budget — Fold Family Comparison')
    ax.grid(True, alpha=0.3)
    ax.invert_xaxis()
    plt.tight_layout()
    plt.savefig('pareto_frontier.png', dpi=200)
    print('Saved pareto_frontier.png')
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python3 analyze_results.py results/*.root')
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
    with open('geant4_differentiated_mesh.json', 'w') as f:
        json.dump(all_summaries, f, indent=2)
    print('\nWrote geant4_differentiated_mesh.json')
    print('\nNOTE: to build the Pareto frontier plot, you also need coverage')
    print('metrics from measure_coverage.py for each geometry (run separately')
    print('on deployed+stowed STL pairs), then call build_pareto_frontier()')
    print('with both datasets - see combine_results.py for an automated')
    print('assembly helper.')