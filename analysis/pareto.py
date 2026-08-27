import sys
import os
import re
import glob
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
EXCLUDED_N_SWEEP_POINTS = {'yoshimura_N6'}
EXCLUDED_THETA_SWEEP_POINTS = {'kresling_theta60'}
HITS_KEYWORDS = {'Reference Barrel (Control)': ['barrel'], 'Kresling': ['kresling'], 'Cylindrical Miura Ori': ['miura'], 'Yoshimura': ['yoshimura', 'yoshimua'], 'Origami Plate (Control)': ['plate']}
NO_DEADZONE_GEOMETRIES = {'Reference Barrel (Control)', 'Origami Plate (Control)'}

def normalize_geom_name(name):
    s = str(name).lower()
    if 'barrel' in s:
        return 'Reference Barrel (Control)'
    elif 'kresling' in s:
        return 'Kresling'
    elif 'miura' in s:
        return 'Cylindrical Miura Ori'
    elif 'yoshimura' in s or 'yoshimua' in s:
        return 'Yoshimura'
    elif 'plate' in s:
        return 'Origami Plate (Control)'
    return name

def _load_json_list(path):
    if not os.path.exists(path):
        if os.path.exists(os.path.join('results', path)):
            path = os.path.join('results', path)
    if not os.path.exists(path):
        return None
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)

def load_geant4_multi_run(multi_run_path):
    data = _load_json_list(multi_run_path)
    if not data:
        return {}
    groups = {}
    for item in data:
        raw_name = item.get('geometry', 'Unknown')
        base = re.sub('_run_\\d+$', '', raw_name)
        groups.setdefault(base, []).append(item)
    candidates = {}
    for base, items in groups.items():
        norm_name = normalize_geom_name(base)
        accs = np.array([it['acceptance_fraction'] for it in items if it.get('acceptance_fraction') is not None], dtype=float)
        x0s = np.array([it['mean_path_length_X0'] for it in items if it.get('mean_path_length_X0') is not None], dtype=float)
        n_runs = len(items)

        def _mean_sd_sem(arr):
            if len(arr) == 0:
                return (None, None, None)
            mean = float(np.mean(arr))
            if len(arr) > 1:
                sd = float(np.std(arr, ddof=1))
                sem = sd / np.sqrt(len(arr))
            else:
                sd, sem = (0.0, 0.0)
            return (mean, sd, sem)
        acc_mean, acc_sd, acc_sem = _mean_sd_sem(accs)
        x0_mean, x0_sd, x0_sem = _mean_sd_sem(x0s)
        n_thrown_total = sum((it.get('n_events_total', 0) or 0 for it in items))
        n_accepted_total = sum((it.get('n_events_hit', 0) or 0 for it in items))
        entry = {'name': norm_name, 'raw_name': base, 'acceptance': acc_mean if acc_mean is not None else 0.0, 'acceptance_sd': acc_sd, 'acceptance_sem': acc_sem, 'x0': x0_mean, 'x0_sd': x0_sd, 'x0_sem': x0_sem, 'n_thrown': n_thrown_total, 'n_accepted': n_accepted_total, 'n_runs': n_runs}
        candidates[norm_name] = entry
    return candidates

def load_geant4_summaries(geant4_path):
    if not os.path.exists(geant4_path):
        if os.path.exists(os.path.join('results', geant4_path)):
            geant4_path = os.path.join('results', geant4_path)
    if not os.path.exists(geant4_path):
        print(f'Error: {geant4_path} not found.')
        return {}
    with open(geant4_path, 'r', encoding='utf-8') as f:
        geant4_data = json.load(f)
    candidates = {}
    for item in geant4_data:
        raw_name = item.get('geometry', 'Unknown')
        norm_name = normalize_geom_name(raw_name)
        acc = item.get('acceptance_fraction', None)
        n_thrown = item.get('n_thrown')
        n_accepted = item.get('n_accepted')
        acc_sem = None
        if acc is not None and n_thrown:
            acc_sem = float(np.sqrt(max(acc * (1.0 - acc), 0.0) / n_thrown))
        entry = {'name': norm_name, 'raw_name': raw_name, 'acceptance': float(acc) if acc is not None else 0.0, 'acceptance_sem': acc_sem, 'x0': item.get('mean_path_length_X0'), 'x0_sem': item.get('mean_path_length_X0_sem'), 'n_thrown': n_thrown, 'n_accepted': n_accepted}
        if norm_name not in candidates or entry['acceptance'] > candidates[norm_name]['acceptance']:
            candidates[norm_name] = entry
    return candidates

def load_deadzone_csv(deadzone_path):
    deadzone_data = {}
    if not deadzone_path or not os.path.exists(deadzone_path):
        return deadzone_data
    df = pd.read_csv(deadzone_path)
    cols_lower = {c.strip().lower(): c for c in df.columns}
    struct_col = cols_lower.get('structure')
    if struct_col is None:
        print(f"Warning: could not find 'structure' column in {deadzone_path}; found columns: {list(df.columns)}")
        return deadzone_data
    if 'mean' in cols_lower and 'sem' in cols_lower:
        mean_col, sem_col = (cols_lower['mean'], cols_lower['sem'])
        for _, row in df.iterrows():
            norm_k = normalize_geom_name(row[struct_col])
            dz = row[mean_col]
            if pd.isna(dz):
                continue
            deadzone_data[norm_k] = {'deadzone': float(dz), 'deadzone_sem': float(row[sem_col]) if pd.notna(row[sem_col]) else None}
        return deadzone_data
    dz_col = next((c for c in df.columns if 'dead_zone_fraction' in c.strip().lower()), None)
    if dz_col is None:
        print(f"Warning: could not find 'dead_zone_fraction' column in {deadzone_path}; found columns: {list(df.columns)}")
        return deadzone_data
    by_struct = {}
    for _, row in df.iterrows():
        norm_k = normalize_geom_name(row[struct_col])
        dz = row[dz_col]
        if pd.isna(dz):
            continue
        by_struct.setdefault(norm_k, []).append(float(dz))
    for norm_k, vals in by_struct.items():
        arr = np.array(vals)
        mean = float(np.mean(arr))
        sem = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr) > 1 else 0.0
        deadzone_data[norm_k] = {'deadzone': mean, 'deadzone_sem': sem}
    return deadzone_data

def match_hits_file(norm_name, hits_files):
    keywords = HITS_KEYWORDS.get(norm_name, [])
    for fpath in hits_files:
        base = os.path.basename(fpath).lower()
        if any((kw in base for kw in keywords)):
            return fpath
    return None

def load_x0_from_hits_dir(hits_dir, candidates):
    if not hits_dir or not os.path.isdir(hits_dir):
        return candidates
    hits_files = glob.glob(os.path.join(hits_dir, '*.csv'))
    for norm_name, c in candidates.items():
        fpath = match_hits_file(norm_name, hits_files)
        if fpath is None:
            continue
        try:
            df = pd.read_csv(fpath)
        except Exception as e:
            print(f'Warning: failed to read {fpath}: {e}')
            continue
        col = next((col for col in df.columns if col.strip().lower() == 'path_x0'), None)
        if col is None:
            print(f'Warning: no path_X0 column found in {fpath}')
            continue
        x0_vals = df[col].dropna().to_numpy()
        if len(x0_vals) == 0:
            continue
        c['x0'] = float(np.mean(x0_vals))
        c['x0_sem'] = float(np.std(x0_vals, ddof=1) / np.sqrt(len(x0_vals)))
        c['x0_n'] = int(len(x0_vals))
        c['hits_file'] = fpath
    return candidates

def merge_datasets(geant4_candidates, deadzone_data):
    candidates = []
    for norm_name, c in geant4_candidates.items():
        c = dict(c)
        if norm_name in NO_DEADZONE_GEOMETRIES:
            c['deadzone'] = None
            c['deadzone_sem'] = None
        else:
            dz_entry = deadzone_data.get(norm_name)
            if dz_entry is None:
                c['deadzone'] = None
                c['deadzone_sem'] = None
            else:
                c['deadzone'] = dz_entry['deadzone']
                c['deadzone_sem'] = dz_entry.get('deadzone_sem')
        c['x0'] = c['x0'] if c['x0'] is not None else 0.0
        candidates.append(c)
    return candidates
SWEEP_FAMILIES = {'kresling': 'Kresling', 'miura': 'Cylindrical Miura Ori', 'yoshimura': 'Yoshimura'}
SWEEP_COLORS = {'kresling': '#ff7f0e', 'miura': '#2ca02c', 'yoshimura': '#d62728'}
SWEEP_SHORT_LABELS = {'kresling': 'Kresling', 'miura': 'Miura-ori', 'yoshimura': 'Yoshimura'}

def _load_sweep_points(sweep_path, param_name, excluded_points, stats_family='scan'):
    data = _load_json_list(sweep_path)
    if not data:
        return {}
    pattern = re.compile(f'^(fixed|scan)_origami_(kresling|miura|yoshimura)_{param_name}(\\d+)_fold')
    points = {fam: [] for fam in SWEEP_FAMILIES}
    for item in data:
        raw_name = item.get('geometry', '')
        m = pattern.match(raw_name)
        if not m:
            continue
        stats_kind, family, value_str = m.groups()
        if stats_kind != stats_family:
            continue
        exclusion_key = f'{family}_{param_name}{value_str}'
        if exclusion_key in excluded_points:
            continue
        acc = item.get('acceptance_fraction')
        x0 = item.get('mean_path_length_X0')
        x0_sem = item.get('sem_path_length_X0')
        if acc is None or x0 is None:
            continue
        points[family].append((int(value_str), float(acc), float(x0), float(x0_sem) if x0_sem is not None else None))
    for family in points:
        points[family].sort(key=lambda t: t[0])
    return points

def load_n_sweep(n_sweep_path):
    return _load_sweep_points(n_sweep_path, param_name='N', excluded_points=EXCLUDED_N_SWEEP_POINTS, stats_family='scan')

def load_theta_sweep(theta_sweep_path):
    return _load_sweep_points(theta_sweep_path, param_name='theta', excluded_points=EXCLUDED_THETA_SWEEP_POINTS, stats_family='scan')

def compute_pareto_dominance(candidates):
    n = len(candidates)
    for c in candidates:
        c['dominated_by'] = []
        c['is_pareto'] = True
    for i in range(n):
        c1 = candidates[i]
        for j in range(n):
            if i == j:
                continue
            c2 = candidates[j]
            both_have_dz = c1['deadzone'] is not None and c2['deadzone'] is not None
            no_worse = c2['acceptance'] >= c1['acceptance'] and c2['x0'] <= c1['x0'] and (not both_have_dz or c2['deadzone'] <= c1['deadzone'])
            strictly_better = c2['acceptance'] > c1['acceptance'] or c2['x0'] < c1['x0'] or (both_have_dz and c2['deadzone'] < c1['deadzone'])
            if no_worse and strictly_better:
                c1['is_pareto'] = False
                c1['dominated_by'].append(c2['name'])
    return candidates
COLORS = {'Reference Barrel (Control)': '#1f77b4', 'Kresling': '#ff7f0e', 'Cylindrical Miura Ori': '#2ca02c', 'Yoshimura': '#d62728', 'Origami Plate (Control)': '#9467bd'}
SHORT_LABELS = {'Reference Barrel (Control)': 'Barrel reference', 'Kresling': 'Kresling', 'Cylindrical Miura Ori': 'Miura-ori', 'Yoshimura': 'Yoshimura', 'Origami Plate (Control)': 'Origami Plate'}

def _short_label(name):
    return SHORT_LABELS.get(name, name)

def plot_acceptance_vs_x0(candidates, outdir, n_sweep=None, theta_sweep=None):
    fig, ax = plt.subplots(figsize=(9, 6.5))
    core_x0 = [c['x0'] for c in candidates]
    core_acc = [c['acceptance'] for c in candidates]
    x0_lo, x0_hi = (min(core_x0), max(core_x0))
    x0_span = max(x0_hi - x0_lo, 1e-06)
    clip_lo = x0_lo - 1.5 * x0_span
    clip_hi = x0_hi + 3.0 * x0_span

    def _plot_sweep_overlay(sweep, linestyle, overlay_label_suffix):
        if not sweep:
            return
        any_offscale = False
        for family, pts in sweep.items():
            if not pts:
                continue
            color = SWEEP_COLORS.get(family, '#333333')
            xs = np.array([p[2] for p in pts])
            ys = np.array([p[1] for p in pts])
            in_range = xs <= clip_hi
            if not in_range.all():
                any_offscale = True
            n_keep = int(in_range.sum())
            n_keep = min(n_keep + 1, len(xs))
            ax.plot(xs[:n_keep], ys[:n_keep], linestyle=linestyle, linewidth=1.3, color=color, alpha=0.55, zorder=2, marker='d', markersize=5, markeredgecolor='black', markeredgewidth=0.5, clip_on=True)
        if any_offscale:
            ax.annotate('', xy=(0.985, 0.5), xycoords='axes fraction', xytext=(0.94, 0.5), textcoords='axes fraction', arrowprops=dict(arrowstyle='-|>', color='#999999', lw=1.5))
        label = overlay_label_suffix + (' (some points off-scale ->)' if any_offscale else '')
        ax.plot([], [], linestyle=linestyle, linewidth=1.3, color='#555555', alpha=0.7, marker='d', markersize=5, markeredgecolor='black', markeredgewidth=0.5, label=label)
    _plot_sweep_overlay(n_sweep, linestyle=':', overlay_label_suffix='Facet-count (N) sweep trajectory')
    _plot_sweep_overlay(theta_sweep, linestyle='--', overlay_label_suffix='Fold-angle (theta) sweep trajectory')
    for c in candidates:
        name = c['name']
        color = COLORS.get(name, '#333333')
        marker = '*' if c['is_pareto'] else 'o'
        size = 260 if c['is_pareto'] else 180
        xerr = 3 * c['x0_sem'] if c.get('x0_sem') else None
        yerr = 3 * c['acceptance_sem'] if c.get('acceptance_sem') else None
        ax.errorbar(c['x0'], c['acceptance'], xerr=xerr, yerr=yerr, fmt=marker, color=color, markersize=14 if marker == '*' else 10, markeredgecolor='black', markeredgewidth=1.0, ecolor=color, elinewidth=1.6, capsize=4, capthick=1.6, zorder=4)
        ax.annotate(_short_label(name), (c['x0'], c['acceptance']), xytext=(10, 6), textcoords='offset points', fontsize=11)
    ax.set_xlabel('Mean radiation-length traversal (X$_0$)', fontsize=11)
    ax.set_ylabel('Acceptance fraction', fontsize=11)
    ax.set_title('Acceptance vs. Radiation-Length Traversal\n(error bars = 3xSE/SEM, for visibility)', fontsize=13)
    ax.grid(True, alpha=0.4)
    if n_sweep or theta_sweep:
        ax.set_xlim(clip_lo, clip_hi)
        ax.legend(loc='best', fontsize=9)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig_acceptance_vs_x0.png')
    plt.savefig(path, dpi=300)
    plt.close()
    print(f'Saved figure: {path}')

def plot_tradeoff_matrix_3panel(candidates, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fold_candidates = [c for c in candidates if c['deadzone'] is not None]
    ax = axes[0]
    for c in candidates:
        marker = '*' if c['is_pareto'] else 'o'
        color = COLORS.get(c['name'], '#333333')
        ax.scatter(c['x0'], c['acceptance'], s=220, color=color, marker=marker, edgecolor='black')
        ax.annotate(_short_label(c['name']), (c['x0'], c['acceptance']), xytext=(6, 6), textcoords='offset points', fontsize=9)
    ax.set_xlabel('Mean X$_0$')
    ax.set_ylabel('Acceptance fraction')
    ax.set_title('Acceptance vs. X$_0$')
    ax.grid(True, alpha=0.3)
    ax = axes[1]
    for c in fold_candidates:
        marker = '*' if c['is_pareto'] else 'o'
        color = COLORS.get(c['name'], '#333333')
        ax.scatter(c['deadzone'] * 100, c['acceptance'], s=220, color=color, marker=marker, edgecolor='black')
        ax.annotate(_short_label(c['name']), (c['deadzone'] * 100, c['acceptance']), xytext=(6, 6), textcoords='offset points', fontsize=9)
    ax.set_xlabel('Dead-zone fraction (%)')
    ax.set_ylabel('Acceptance fraction')
    ax.set_title('Acceptance vs. Dead-Zone')
    ax.grid(True, alpha=0.3)
    ax = axes[2]
    for c in fold_candidates:
        marker = '*' if c['is_pareto'] else 'o'
        color = COLORS.get(c['name'], '#333333')
        ax.scatter(c['deadzone'] * 100, c['x0'], s=220, color=color, marker=marker, edgecolor='black')
        ax.annotate(_short_label(c['name']), (c['deadzone'] * 100, c['x0']), xytext=(6, 6), textcoords='offset points', fontsize=9)
    ax.set_xlabel('Dead-zone fraction (%)')
    ax.set_ylabel('Mean X$_0$')
    ax.set_title('X$_0$ vs. Dead-Zone')
    ax.grid(True, alpha=0.3)
    plt.suptitle('Three-Objective Trade-off Matrix (Acceptance, X$_0$, Dead-Zone)', fontsize=15, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(outdir, 'fig_tradeoff_matrix_3panel.png')
    plt.savefig(path, dpi=300)
    plt.close()
    print(f'Saved figure: {path}')

def plot_parallel_coords_3obj(candidates, outdir):
    fig, ax = plt.subplots(figsize=(10, 6.5))
    objs = ['Acceptance\n(Max)', 'Mean X$_0$\n(Min)', 'Dead-Zone Fraction\n(Min, fold geometries only)']
    x_axis = np.arange(len(objs))
    acc_vals = [c['acceptance'] for c in candidates]
    x0_vals = [c['x0'] for c in candidates]
    dz_vals = [c['deadzone'] for c in candidates if c['deadzone'] is not None]
    acc_min, acc_max = (min(acc_vals), max(acc_vals))
    x0_min, x0_max = (min(x0_vals), max(x0_vals))
    dz_min, dz_max = (min(dz_vals), max(dz_vals)) if dz_vals else (0.0, 1.0)
    for c in candidates:
        n_acc = (c['acceptance'] - acc_min) / (acc_max - acc_min) if acc_max > acc_min else 1.0
        n_x0 = (x0_max - c['x0']) / (x0_max - x0_min) if x0_max > x0_min else 1.0
        color = COLORS.get(c['name'], '#333333')
        linestyle = '-' if c['is_pareto'] else '--'
        linewidth = 2.5 if c['is_pareto'] else 1.2
        alpha = 0.95 if c['is_pareto'] else 0.6
        if c['deadzone'] is not None:
            n_dz = (dz_max - c['deadzone']) / (dz_max - dz_min) if dz_max > dz_min else 1.0
            scores = [n_acc, n_x0, n_dz]
            markers = ['o', 'o', 'o']
            label = _short_label(c['name'])
            ax.plot(x_axis, scores, marker='o', linestyle=linestyle, linewidth=linewidth, color=color, alpha=alpha, label=label)
        else:
            scores = [n_acc, n_x0]
            ax.plot(x_axis[:2], scores, marker='o', linestyle=linestyle, linewidth=linewidth, color=color, alpha=alpha)
            ax.plot(x_axis[2], 1.0, marker='X', markersize=12, color=color, markeredgecolor='black', alpha=alpha)
            label = f"{_short_label(c['name'])} (no dead-zone value)"
            ax.plot([], [], marker='o', linestyle=linestyle, linewidth=linewidth, color=color, alpha=alpha, label=label)
    ax.set_xticks(x_axis)
    ax.set_xticklabels(objs, fontweight='bold')
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel('Normalized objective score (1.0 = best)', fontweight='bold')
    ax.set_title('Parallel Coordinates - Normalized 3-Objective Performance', fontweight='bold', pad=12)
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
    plt.tight_layout()
    path = os.path.join(outdir, 'fig_parallel_coords_3obj.png')
    plt.savefig(path, dpi=300)
    plt.close()
    print(f'Saved figure: {path}')

def plot_acceptance_vs_deadzone(candidates, outdir):
    fold_candidates = [c for c in candidates if c['deadzone'] is not None]
    if not fold_candidates:
        print('Skipping fig_acceptance_vs_deadzone.png: no geometries have a dead-zone value.')
        return
    fig, ax = plt.subplots(figsize=(9, 7))
    for c in fold_candidates:
        marker = '*' if c['is_pareto'] else 'o'
        color = COLORS.get(c['name'], '#333333')
        ax.scatter(c['deadzone'] * 100, c['acceptance'], s=260, color=color, marker=marker, edgecolor='black')
        ax.annotate(_short_label(c['name']), (c['deadzone'] * 100, c['acceptance']), xytext=(8, 8), textcoords='offset points', fontsize=11)
    ax.set_xlabel('Dead-zone fraction (%)')
    ax.set_ylabel('Acceptance fraction')
    ax.set_title('Acceptance vs. Dead-Zone Fraction\n(fold geometries only - barrel dead-zone fraction is 0 by definition)', fontsize=13)
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig_acceptance_vs_deadzone.png')
    plt.savefig(path, dpi=300)
    plt.close()
    print(f'Saved figure: {path}')

def plot_x0_vs_deadzone(candidates, outdir):
    fold_candidates = [c for c in candidates if c['deadzone'] is not None]
    if not fold_candidates:
        print('Skipping fig_x0_vs_deadzone.png: no geometries have a dead-zone value.')
        return
    fig, ax = plt.subplots(figsize=(9, 7))
    for c in fold_candidates:
        marker = '*' if c['is_pareto'] else 'o'
        color = COLORS.get(c['name'], '#333333')
        ax.scatter(c['deadzone'] * 100, c['x0'], s=260, color=color, marker=marker, edgecolor='black')
        ax.annotate(_short_label(c['name']), (c['deadzone'] * 100, c['x0']), xytext=(8, 8), textcoords='offset points', fontsize=11)
    ax.set_xlabel('Dead-zone fraction (%)')
    ax.set_ylabel('Mean radiation-length traversal (X$_0$)')
    ax.set_title('Radiation-Length Traversal vs. Dead-Zone Fraction\n(fold geometries only)', fontsize=13)
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    path = os.path.join(outdir, 'fig_x0_vs_deadzone.png')
    plt.savefig(path, dpi=300)
    plt.close()
    print(f'Saved figure: {path}')

def _plot_sweep_2panel(sweep, outdir, filename, param_label, title, x0_err_available=True):
    families_with_data = [fam for fam, pts in sweep.items() if pts]
    if not families_with_data:
        print(f'Skipping {filename}: no sweep data found.')
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax_acc, ax_x0 = axes
    for family in families_with_data:
        pts = sweep[family]
        xs = [p[0] for p in pts]
        accs = [p[1] for p in pts]
        x0s = [p[2] for p in pts]
        x0_sems = [p[3] for p in pts]
        color = SWEEP_COLORS.get(family, '#333333')
        label = SWEEP_SHORT_LABELS.get(family, family)
        ax_acc.plot(xs, accs, marker='o', color=color, linewidth=2, label=label)
        if x0_err_available and all((s is not None for s in x0_sems)):
            ax_x0.errorbar(xs, x0s, yerr=[3 * s for s in x0_sems], marker='o', color=color, linewidth=2, capsize=3, label=label)
        else:
            ax_x0.plot(xs, x0s, marker='o', color=color, linewidth=2, label=label)
    ax_acc.set_xlabel(param_label)
    ax_acc.set_ylabel('Acceptance fraction')
    ax_acc.set_title(f'Acceptance vs. {param_label}')
    ax_acc.grid(True, alpha=0.3)
    ax_acc.legend()
    ax_x0.set_xlabel(param_label)
    ax_x0.set_ylabel('Mean X$_0$')
    ax_x0.set_title(f'Radiation-Length Traversal vs. {param_label}')
    ax_x0.grid(True, alpha=0.3)
    ax_x0.legend()
    plt.suptitle(title, fontsize=15, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(outdir, filename)
    plt.savefig(path, dpi=300)
    plt.close()
    print(f'Saved figure: {path}')

def plot_n_sweep(n_sweep, outdir):
    _plot_sweep_2panel(n_sweep, outdir, 'fig_n_sweep.png', param_label='Fold count N', title='Fold-Count (N) Sweep - Acceptance and X$_0$ vs. N')

def plot_theta_sweep(theta_sweep, outdir):
    _plot_sweep_2panel(theta_sweep, outdir, 'fig_theta_sweep.png', param_label='Fold angle theta (deg)', title='Fold-Angle (theta) Sweep - Acceptance and X$_0$ vs. theta')

def plot_all(candidates, outdir, n_sweep=None, theta_sweep=None):
    os.makedirs(outdir, exist_ok=True)
    plot_acceptance_vs_x0(candidates, outdir, n_sweep=n_sweep, theta_sweep=theta_sweep)
    plot_tradeoff_matrix_3panel(candidates, outdir)
    plot_parallel_coords_3obj(candidates, outdir)
    plot_acceptance_vs_deadzone(candidates, outdir)
    plot_x0_vs_deadzone(candidates, outdir)
    if n_sweep:
        plot_n_sweep(n_sweep, outdir)
    if theta_sweep:
        plot_theta_sweep(theta_sweep, outdir)

def print_summary(candidates):
    has_run_stats = any((c.get('n_runs') for c in candidates))
    print('\n' + '=' * 110)
    print('3D PARETO FRONTIER ANALYSIS SUMMARY TABLE')
    if has_run_stats:
        print('(Acceptance/X0 are across-run means +/- SEM where multi-run data was available; dead-zone is the 5-run aggregate mean +/- SEM where available.)')
    print('=' * 110)
    header = f"{'Geometry Name':<30} | {'Acceptance':<16} | {'Mean X0':<16} | {'DeadZone Frac':<18} | {'Runs':<5} | {'Status':<15}"
    print(header)
    print('-' * 110)
    for c in candidates:
        status = 'Pareto Optimal' if c['is_pareto'] else 'Dominated'
        acc_str = f"{c['acceptance']:.4f}"
        if c.get('acceptance_sem'):
            acc_str += f" +/-{c['acceptance_sem']:.4f}"
        x0_str = f"{c['x0']:.6f}"
        if c.get('x0_sem'):
            x0_str += f" +/-{c['x0_sem']:.1e}"
        if c['deadzone'] is not None:
            dz_str = f"{c['deadzone']:.4f}"
            if c.get('deadzone_sem'):
                dz_str += f" +/-{c['deadzone_sem']:.4f}"
        else:
            dz_str = 'n/a'
        n_runs_str = str(c['n_runs']) if c.get('n_runs') else '1'
        row = f"{c['name']:<30} | {acc_str:<16} | {x0_str:<16} | {dz_str:<18} | {n_runs_str:<5} | {status:<15}"
        print(row)
    print('=' * 110)
    print('\nPARETO OPTIMAL SET (NON-DOMINATED FRONTIER):')
    for c in candidates:
        if c['is_pareto']:
            dz_str = f"{c['deadzone']:.2%}" if c['deadzone'] is not None else 'n/a'
            print(f"  * {c['name']} (Acc: {c['acceptance']:.2%}, X0: {c['x0']:.6f}, DZ: {dz_str})")
    dominated_set = [c for c in candidates if not c['is_pareto']]
    if dominated_set:
        print('\nDOMINATED DESIGN CANDIDATES:')
        for c in dominated_set:
            dom_by = ', '.join(c['dominated_by'])
            print(f"  * {c['name']} -> Dominated by [{dom_by}]")
    print()

def parse_args():
    p = argparse.ArgumentParser(description='3D Pareto Frontier Analysis for Origami Tracker Geometries.')
    p.add_argument('geant4_summaries', nargs='?', default='geant4_final_result_summaries_multi_run_.json', help='Path to the primary GEANT4 summary file.')
    p.add_argument('deadzone_csv', nargs='?', default='original_four_dead_zone_summary_aggregate.csv', help='Path to the dead-zone summary CSV.')
    p.add_argument('--n-sweep', default='geant4_N_Sweep_summaries.json', help='Path to the facet-count (N) sweep summary file.')
    p.add_argument('--theta-sweep', default='geant4_N_Sweep_results_.json', help='Path to the fold-angle (theta) sweep summary file.')
    p.add_argument('--hits-dir', default=None, help='Optional directory of raw per-hit CSVs.')
    p.add_argument('--outdir', default='.', help='Directory to save output figures into.')
    return p.parse_args()

def main():
    args = parse_args()
    geant4_candidates = load_geant4_multi_run(args.geant4_summaries)
    if geant4_candidates:
        print(f'Loaded multi-run GEANT4 summaries from {args.geant4_summaries} ({len(geant4_candidates)} geometries, across-run mean/sd/sem).')
    else:
        print(f"'{args.geant4_summaries}' not found or empty as a multi-run file; trying legacy single-run geant4_summaries.json ...")
        geant4_candidates = load_geant4_summaries('geant4_summaries.json')
    if not geant4_candidates:
        print('ERROR: No valid geometry candidates found in input JSON file(s).')
        sys.exit(1)
    deadzone_data = load_deadzone_csv(args.deadzone_csv)
    if not deadzone_data:
        print(f"'{args.deadzone_csv}' not found or empty; trying legacy deadzone_summary.csv ...")
        deadzone_data = load_deadzone_csv('deadzone_summary.csv')
    if args.hits_dir:
        geant4_candidates = load_x0_from_hits_dir(args.hits_dir, geant4_candidates)
    candidates = merge_datasets(geant4_candidates, deadzone_data)
    candidates = compute_pareto_dominance(candidates)
    print_summary(candidates)
    n_sweep = load_n_sweep(args.n_sweep) if args.n_sweep else {}
    theta_sweep = load_theta_sweep(args.theta_sweep) if args.theta_sweep else {}
    plot_all(candidates, args.outdir, n_sweep=n_sweep, theta_sweep=theta_sweep)
if __name__ == '__main__':
    main()