import sys
import glob
import os
import re
import json
import numpy as np
import uproot
import matplotlib.pyplot as plt
SILICON_X0_MM = 93.66
KAPTON_X0_MM = 285.75

def load_run(root_path):
    with uproot.open(root_path) as f:
        available = [k.split(';')[0] for k in f.keys()]
        if 'PionEvents' in available:
            name = 'PionEvents'
        elif 'MuonEvents' in available:
            name = 'MuonEvents'
        elif len(available) == 1:
            name = available[0]
            print(f"  NOTE: neither 'PionEvents' nor 'MuonEvents' found in {root_path}; using the only tree present, '{name}'.")
        else:
            raise KeyError(f"Can't find a ntuple tree in {root_path} - available keys: {available}")
        tree = f[name]
        data = tree.arrays(library='np')
    return data

def _parse_filename_tag(filename_stem, prefix, unit_suffix=''):
    pattern = re.escape(prefix) + '(-?\\d+p?\\d*)' + (re.escape(unit_suffix) if unit_suffix else '(?:_|$)')
    m = re.search(pattern, filename_stem)
    if m is None:
        return None
    return float(m.group(1).replace('p', '.'))

def resolve_geometry_name(path, data):
    filename_stem = os.path.basename(path).replace('.root', '')
    tag_branch = data.get('structureTag')
    if tag_branch is None or len(tag_branch) == 0:
        print(f"  NOTE: no structureTag branch found - using filename stem '{filename_stem}' as geometry name.")
        return filename_stem
    unique_tags = np.unique(tag_branch)
    if len(unique_tags) != 1:
        print(f'  WARNING: file contains {len(unique_tags)} distinct structureTag values {unique_tags.tolist()} - not a single clean run. Falling back to filename stem.')
        return filename_stem
    tag = str(unique_tags[0])
    if tag not in filename_stem:
        print(f"  WARNING: internal structureTag='{tag}' does not appear in filename '{filename_stem}'. Using the structureTag (trust the data over the label) - double check this file.")
    return tag

def validate_vertex_smear(path, data, tol_mm=1.0):
    filename_stem = os.path.basename(path).replace('.root', '')
    claimed_mm = _parse_filename_tag(filename_stem, '_vtx', 'mm')
    if claimed_mm is None:
        return None
    issues = []
    axis_stds = {}
    for axis in ('vertexX_mm', 'vertexY_mm', 'vertexZ_mm'):
        col = data.get(axis)
        if col is None or len(col) == 0:
            continue
        axis_stds[axis] = float(np.std(col))
    if not axis_stds:
        return None
    for axis, measured in axis_stds.items():
        if abs(measured - claimed_mm) > tol_mm:
            issues.append(f'{axis} std={measured:.2f}mm vs filename-claimed {claimed_mm:.1f}mm')
    if issues:
        print(f"  WARNING: vertex smear mismatch in '{filename_stem}':")
        for issue in issues:
            print(f'    - {issue}')
        return {'claimed_mm': claimed_mm, 'measured_std_mm': axis_stds, 'ok': False}
    return {'claimed_mm': claimed_mm, 'measured_std_mm': axis_stds, 'ok': True}

def summarize_robustness(name, sigma_mm, path, data):
    hit = data['hitDetector'].astype(bool)
    n_total = len(hit)
    n_hit = int(hit.sum())
    acceptance_fraction = n_hit / n_total if n_total > 0 else 0.0
    vx = data.get('vertexX_mm')
    vy = data.get('vertexY_mm')
    vz = data.get('vertexZ_mm')
    if vx is None or vy is None or vz is None:
        print(f'  WARNING: no vertexX/Y/Z_mm columns in this file - was it produced by a Mark 3 build? Skipping robustness summary.')
        return None
    r_mm = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    path_mm = data['totalPathLength_mm']
    local_inc = data.get('localIncidenceDeg')
    n_hits_primary = data.get('nHitsPrimary')
    has_split = 'siliconPathLength_mm' in data and 'kaptonPathLength_mm' in data
    if has_split:
        si_mm = data['siliconPathLength_mm']
        kap_mm = data['kaptonPathLength_mm']
        x_over_x0_all = si_mm / SILICON_X0_MM + kap_mm / KAPTON_X0_MM
    else:
        x_over_x0_all = path_mm / SILICON_X0_MM
    path_X0_hit = x_over_x0_all[hit]
    r_hit = r_mm[hit]
    local_inc_valid = None
    r_local_inc_valid = None
    if local_inc is not None:
        valid_local = hit & (local_inc >= 0.0)
        local_inc_valid = local_inc[valid_local]
        r_local_inc_valid = r_mm[valid_local]
    result = {'geometry': name, 'vertex_smear_sigma_mm': sigma_mm, 'n_events_total': n_total, 'n_events_hit': n_hit, 'acceptance_fraction': acceptance_fraction, 'mean_realized_displacement_mm': float(np.mean(r_mm)), 'std_realized_displacement_mm': float(np.std(r_mm)), 'mean_path_length_X0': float(np.mean(path_X0_hit)) if n_hit > 0 else None, 'std_path_length_X0': float(np.std(path_X0_hit, ddof=1)) if n_hit > 1 else None, 'mean_local_incidence_deg': float(np.mean(local_inc_valid)) if local_inc_valid is not None and len(local_inc_valid) > 0 else None, 'mean_nHitsPrimary': float(np.mean(n_hits_primary[hit])) if n_hits_primary is not None and n_hit > 0 else None, '_r_hit_mm': r_hit, '_path_X0_hit': path_X0_hit, '_r_local_inc_mm': r_local_inc_valid, '_local_inc_deg': local_inc_valid, '_r_all_mm': r_mm, '_hit_mask': hit}
    return result

def _bin_by_displacement(r_mm, values, n_bins=10, r_max=None, agg='mean'):
    r_mm = np.asarray(r_mm, dtype=float)
    if len(r_mm) == 0:
        return {'bin_centers_mm': [], 'values': [], 'n_in_bin': []}
    if r_max is None:
        r_max = float(np.percentile(r_mm, 99))
    bins = np.linspace(0, r_max, n_bins + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    out_values = []
    n_in_bin = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        in_bin = (r_mm >= lo) & (r_mm < hi) if hi < r_max else (r_mm >= lo) & (r_mm <= hi)
        n = int(in_bin.sum())
        n_in_bin.append(n)
        if n == 0:
            out_values.append(None)
            continue
        if agg == 'mean':
            out_values.append(float(np.mean(np.asarray(values)[in_bin])))
        elif agg == 'fraction':
            out_values.append(float(np.mean(np.asarray(values)[in_bin])))
        else:
            raise ValueError(f"unknown agg '{agg}'")
    return {'bin_centers_mm': centers.tolist(), 'values': out_values, 'n_in_bin': n_in_bin}

def plot_vs_sigma(summaries, out_path='exp3_vs_sigma.png'):
    by_geometry = {}
    for s in summaries:
        if s is None:
            continue
        by_geometry.setdefault(s['geometry'], []).append(s)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    colors = plt.cm.tab10.colors
    panels = [('acceptance_fraction', 'Acceptance (fraction of events hitting sensor)'), ('mean_path_length_X0', 'Mean path length ($X/X_0$)'), ('mean_local_incidence_deg', 'Mean local incidence angle (deg)'), ('mean_nHitsPrimary', 'Mean primary hits per event')]
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for i, (geom, points) in enumerate(sorted(by_geometry.items())):
            points = sorted(points, key=lambda s: s['vertex_smear_sigma_mm'])
            sigmas = [pt['vertex_smear_sigma_mm'] for pt in points]
            vals = [pt[key] for pt in points]
            color = colors[i % len(colors)]
            ax.plot(sigmas, vals, 'o-', color=color, label=geom)
        ax.set_xlabel('Vertex smear $\\sigma$ (mm)')
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle('Experiment 3: outcome vs. vertex-smear sigma (per structure)')
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f'Saved {out_path}')

def plot_vs_displacement(summaries, out_path='exp3_vs_displacement.png'):
    by_geometry = {}
    for s in summaries:
        if s is None:
            continue
        by_geometry.setdefault(s['geometry'], []).append(s)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    colors = plt.cm.tab10.colors
    for i, (geom, points) in enumerate(sorted(by_geometry.items())):
        r_all = np.concatenate([p['_r_all_mm'] for p in points])
        hit_all = np.concatenate([p['_hit_mask'] for p in points])
        r_hit = np.concatenate([p['_r_hit_mm'] for p in points])
        path_X0_hit = np.concatenate([p['_path_X0_hit'] for p in points])
        color = colors[i % len(colors)]
        acc_binned = _bin_by_displacement(r_all, hit_all, agg='fraction')
        ax1.plot(acc_binned['bin_centers_mm'], acc_binned['values'], 'o-', color=color, label=geom)
        path_binned = _bin_by_displacement(r_hit, path_X0_hit, agg='mean')
        ax2.plot(path_binned['bin_centers_mm'], path_binned['values'], 'o-', color=color, label=geom)
    ax1.set_xlabel('Realized vertex displacement r (mm)')
    ax1.set_ylabel('Acceptance (fraction hitting sensor)')
    ax1.set_title('Acceptance vs. per-event displacement')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)
    ax2.set_xlabel('Realized vertex displacement r (mm)')
    ax2.set_ylabel('Mean path length ($X/X_0$), hit events only')
    ax2.set_title('Path length vs. per-event displacement')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)
    fig.suptitle('Experiment 3: outcome vs. realized per-event vertex displacement (pooled across all sigma points)')
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f'Saved {out_path}')

def _json_safe(summary):
    return {k: v for k, v in summary.items() if not k.startswith('_')}
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python3 analyze_exp3.py results_exp3/*.root')
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
        print(f'\nLoading {path}...')
        data = load_run(path)
        name = resolve_geometry_name(path, data)
        vertex_check = validate_vertex_smear(path, data)
        if vertex_check is None:
            print(f"  WARNING: filename has no '_vtxNmm' vertex-smear tag - this file may not be Experiment 3 output. Skipping.")
            continue
        sigma_mm = vertex_check['claimed_mm']
        print(f"  -> geometry: '{name}'   vertex smear sigma: {sigma_mm} mm   [{('ok' if vertex_check['ok'] else 'SUSPECT - filename does not match file contents')}]")
        summary = summarize_robustness(name, sigma_mm, path, data)
        if summary is None:
            continue
        summary['vertex_smear_check'] = vertex_check
        summary['filename_matches_contents'] = vertex_check['ok']
        all_summaries.append(summary)
        print(f"  -> acceptance={summary['acceptance_fraction']:.4f}   mean_path_X0={summary['mean_path_length_X0']}   mean_local_inc_deg={summary['mean_local_incidence_deg']}   realized_disp(mean/std)={summary['mean_realized_displacement_mm']:.3f}/{summary['std_realized_displacement_mm']:.3f} mm")
    if not all_summaries:
        print('\nERROR: no valid Experiment 3 summaries produced - check warnings above (missing vertexX/Y/Z_mm columns, no vertex-smear tag, etc).')
        sys.exit(1)
    n_suspect = sum((1 for s in all_summaries if not s['filename_matches_contents']))
    if n_suspect:
        print(f"\n{'=' * 70}")
        print(f'WARNING: {n_suspect} of {len(all_summaries)} file(s) have a filename that does not match their actual vertex-smear contents - excluded from the plots below. Review by hand.')
        print(f"{'=' * 70}")
    clean_summaries = [s for s in all_summaries if s['filename_matches_contents']]
    with open('exp3_robustness_summary.json', 'w') as f:
        json.dump([_json_safe(s) for s in all_summaries], f, indent=2)
    print(f'\nWrote exp3_robustness_summary.json ({len(all_summaries)} run(s), all flags included)')
    plot_vs_sigma(clean_summaries)
    plot_vs_displacement(clean_summaries)
    print("\nNOTE: 'realized displacement' (per-event r = sqrt(vx^2+vy^2+vz^2)) will not exactly equal the nominal sigma a run was configured with - for a 3D isotropic Gaussian smear with per-axis sigma, r follows a Maxwell-Boltzmann-like distribution with mean ~1.6*sigma, not sigma itself. This is expected and is why plot_vs_displacement() bins on the actual per-event r rather than assuming r==sigma.")