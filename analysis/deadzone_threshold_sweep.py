import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
STRUCTURE_COLORS = {'miura': '#0072B2', 'kresling': '#D55E00', 'yoshimura': '#009E73'}
FALLBACK_COLORS = ['#CC79A7', '#E69F00', '#56B4E9', '#000000']

def _color_for(structure: str, fallback_idx: int) -> str:
    key = str(structure).lower()
    for name, color in STRUCTURE_COLORS.items():
        if key == name or key.startswith(name + '_'):
            return color
    return FALLBACK_COLORS[fallback_idx % len(FALLBACK_COLORS)]

def _use_clean_style():
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'], 'axes.edgecolor': '#333333', 'axes.labelcolor': '#222222', 'text.color': '#222222', 'xtick.color': '#444444', 'ytick.color': '#444444', 'axes.titlesize': 13, 'axes.titleweight': 'bold', 'axes.labelsize': 11, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
from deadzone_map import _load_all_hits, is_structure, is_excluded_structure, build_maps, dead_zone_fraction

def compute_ratio_maps(hits_dir: str, barrel_hits_dir: str=None) -> dict:
    if barrel_hits_dir is None:
        barrel_hits_dir = hits_dir
    all_hits = _load_all_hits(hits_dir)
    barrel_pool = all_hits if barrel_hits_dir == hits_dir else _load_all_hits(barrel_hits_dir)
    is_barrel_mask = barrel_pool['structure'].apply(lambda s: is_structure(str(s), 'barrel'))
    barrel_rows = barrel_pool[is_barrel_mask]
    if barrel_rows.empty:
        seen = sorted(barrel_pool['structure'].unique()) if not barrel_pool.empty else []
        raise FileNotFoundError(f"No rows identifying as 'barrel' found in '{barrel_hits_dir}'. structure values actually present: {seen}. Can't build per-run baselines.")
    barrel_by_run = {run: df for run, df in barrel_rows.groupby('run')}
    print(f'Found barrel baselines for runs: {sorted(barrel_by_run)}')
    fold_hits_pool = all_hits[~is_barrel_mask] if barrel_hits_dir == hits_dir else all_hits[~all_hits['structure'].apply(lambda s: is_structure(str(s), 'barrel'))]
    ratio_maps = {}
    for (structure, run), fold_hits in fold_hits_pool.groupby(['structure', 'run']):
        if is_excluded_structure(str(structure)):
            print(f"skipping excluded structure '{structure}' (run {run})")
            continue
        if run not in barrel_by_run:
            print(f"skipping {structure} run {run}: no matching barrel run in '{barrel_hits_dir}' (available: {sorted(barrel_by_run)}) -- refusing to fall back to a different run's baseline.")
            continue
        barrel_hits = barrel_by_run[run]
        ratio, _, _, _ = build_maps(fold_hits, barrel_hits)
        ratio_maps[structure, run] = ratio
        print(f'built ratio map for {structure} run {run}')
    if not ratio_maps:
        raise RuntimeError(f"No (structure, run) ratio maps could be built from '{hits_dir}'. Check that hits CSVs have 'structure'/'run' columns and matching barrel runs exist.")
    return ratio_maps

def sweep(ratio_maps: dict, n_thresholds: int=101) -> tuple[pd.DataFrame, pd.DataFrame]:
    thresholds = np.linspace(0.0, 1.0, n_thresholds)
    rows = []
    for (structure, run), ratio in ratio_maps.items():
        for t in thresholds:
            dz = dead_zone_fraction(ratio, threshold=t)
            rows.append({'structure': structure, 'run': run, 'threshold': t, 'dead_zone_fraction': dz})
    per_run = pd.DataFrame(rows)
    agg = per_run.groupby(['structure', 'threshold'])['dead_zone_fraction'].agg(mean='mean', std='std', n='count').reset_index()
    agg['sem'] = agg['std'] / np.sqrt(agg['n'])
    return (per_run, agg)

def plot_sweep(agg: pd.DataFrame, out_path: str, mark_threshold: float=1.0, display_names: dict=None):
    _use_clean_style()
    structures = sorted(agg['structure'].unique())
    display_names = display_names or {}
    labels = [display_names.get(s, s) for s in structures]
    colors = [_color_for(s, i) for i, s in enumerate(structures)]
    fig = plt.figure(figsize=(16, 6.4))
    ax0 = fig.add_subplot(1, 2, 1, projection='3d')
    for i, structure in enumerate(structures):
        sub = agg[agg['structure'] == structure].sort_values('threshold')
        x = sub['threshold'].to_numpy()
        z = (sub['mean'] * 100).to_numpy()
        y = np.full(len(sub), i)
        verts = [list(zip(x, y, z)) + [(x[-1], y[-1], 0), (x[0], y[0], 0)]]
        wall = Poly3DCollection(verts, facecolor=colors[i], edgecolor='none', alpha=0.28)
        ax0.add_collection3d(wall)
        ax0.plot(x, y, z, color=colors[i], linewidth=2.4, zorder=5)
    ax0.set_xlabel('Dead-zone threshold (ratio cutoff)', labelpad=10)
    ax0.set_yticks(range(len(structures)))
    ax0.set_yticklabels(labels)
    ax0.set_zlabel('Dead-zone fraction (%)', labelpad=8)
    ax0.set_xlim(0, 1)
    ax0.set_ylim(-0.5, len(structures) - 0.5)
    ax0.set_zlim(bottom=0)
    ax0.set_title('Dead-zone fraction vs. threshold, per geometry', pad=14)
    ax0.view_init(elev=22, azim=-55)
    for pane in (ax0.xaxis.pane, ax0.yaxis.pane, ax0.zaxis.pane):
        pane.set_facecolor((1, 1, 1, 0.0))
        pane.set_edgecolor((0, 0, 0, 0.08))
    ax0.xaxis._axinfo['grid']['color'] = (0, 0, 0, 0.08)
    ax0.yaxis._axinfo['grid']['color'] = (0, 0, 0, 0.08)
    ax0.zaxis._axinfo['grid']['color'] = (0, 0, 0, 0.08)
    ax1 = fig.add_subplot(1, 2, 2)
    ax1.set_xlim(0, 1.05)
    for i, structure in enumerate(structures):
        sub = agg[agg['structure'] == structure].sort_values('threshold')
        x = sub['threshold'].to_numpy()
        mean = (sub['mean'] * 100).to_numpy()
        sem = (sub['sem'] * 100).to_numpy()
        ax1.plot(x, mean, color=colors[i], linewidth=2.2, solid_capstyle='round')
        ax1.fill_between(x, mean - sem, mean + sem, color=colors[i], alpha=0.16, linewidth=0)
    ax1.set_xlabel('Dead-zone threshold (ratio cutoff)')
    ax1.set_ylabel('Dead-zone fraction (%)')
    ax1.set_title('Same trend, 2D (mean ± SEM across runs)', pad=14)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(axis='y', color='#DDDDDD', linewidth=0.8, zorder=0)
    ax1.set_axisbelow(True)
    if mark_threshold is not None:
        ax1.axvline(mark_threshold, color='#999999', linestyle='--', linewidth=1, zorder=0)
        ax1.text(mark_threshold, 0.01, 'current fixed\nthreshold', ha='center', va='bottom', fontsize=8, color='#999999', style='italic', transform=ax1.get_xaxis_transform())
    y0, y1 = ax1.get_ylim()
    end_vals = []
    for i, structure in enumerate(structures):
        sub = agg[agg['structure'] == structure].sort_values('threshold')
        end_vals.append((sub['mean'].iloc[-1] * 100, labels[i], colors[i]))
    end_vals.sort(reverse=True)
    label_y_step = 0.1 * (y1 - y0)
    top_y = y1 * 0.94
    for rank, (val, label, color) in enumerate(end_vals):
        y_pos = top_y - rank * label_y_step
        ax1.text(1.02, y_pos, f'{label}\n{val:.1f}%', color=color, fontsize=9.5, fontweight='bold', ha='left', va='top', linespacing=1.35, transform=ax1.get_yaxis_transform(), clip_on=False)
    zoom_lo = max(0.0, (mark_threshold if mark_threshold is not None else 1.0) - 0.1)
    zoom_hi = min(1.0, (mark_threshold if mark_threshold is not None else 1.0) + 0.02)
    zoom_mask_any = (agg['threshold'] >= zoom_lo) & (agg['threshold'] <= zoom_hi)
    if zoom_mask_any.any():
        axins = inset_axes(ax1, width='100%', height='100%', bbox_to_anchor=(1.06, -0.02, 0.34, 0.5), bbox_transform=ax1.transAxes, loc='lower left', borderpad=0)
        for i, structure in enumerate(structures):
            sub = agg[agg['structure'] == structure].sort_values('threshold')
            zsub = sub[(sub['threshold'] >= zoom_lo) & (sub['threshold'] <= zoom_hi)]
            if zsub.empty:
                continue
            x = zsub['threshold'].to_numpy()
            mean = (zsub['mean'] * 100).to_numpy()
            sem = (zsub['sem'] * 100).to_numpy()
            axins.plot(x, mean, color=colors[i], linewidth=2.0, solid_capstyle='round')
            axins.fill_between(x, mean - sem, mean + sem, color=colors[i], alpha=0.18, linewidth=0)
        zoom_rows = agg[zoom_mask_any]
        y_pad = 0.06 * max(zoom_rows['mean'].max() * 100 - zoom_rows['mean'].min() * 100, 1.0)
        axins.set_xlim(zoom_lo, zoom_hi)
        axins.set_ylim(zoom_rows['mean'].min() * 100 - y_pad, zoom_rows['mean'].max() * 100 + y_pad)
        axins.set_title('Zoom', fontsize=8.5, color='#666666', pad=4)
        axins.tick_params(labelsize=7.5, length=2.5)
        axins.grid(axis='y', color='#DDDDDD', linewidth=0.6, zorder=0)
        axins.set_axisbelow(True)
        for spine in axins.spines.values():
            spine.set_color('#888888')
            spine.set_linewidth(0.9)
        axins.patch.set_facecolor('white')
        axins.patch.set_alpha(1.0)
        if mark_threshold is not None and zoom_lo <= mark_threshold <= zoom_hi:
            axins.axvline(mark_threshold, color='#999999', linestyle='--', linewidth=0.8, zorder=0)
        mark_inset(ax1, axins, loc1=1, loc2=2, fc='none', ec='#999999', linewidth=0.8, alpha=0.6)
    fig.subplots_adjust(left=0.05, right=0.8, bottom=0.11, top=0.9, wspace=0.38)
    fig.savefig(out_path, dpi=300, bbox_inches=None)
    plt.close(fig)
    print(f'Saved plot -> {out_path}')

def clean_label(structure: str) -> str:
    key = str(structure).lower()
    for suffix in ('_deployed', '_reference', '_flat'):
        if key.endswith(suffix):
            key = key[:-len(suffix)]
    special = {'miura': 'Miura-ori', 'kresling': 'Kresling', 'yoshimura': 'Yoshimura'}
    return special.get(key, key.replace('_', ' ').title())

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('hits_dir', help='Directory containing the *_hits.csv files.')
    parser.add_argument('barrel_hits_dir', nargs='?', default=None, help='Directory containing barrel hits CSVs (defaults to hits_dir).')
    parser.add_argument('--n-thresholds', type=int, default=101, help='Number of threshold points sampled from 0 to 1 inclusive (default 101, i.e. step 0.01).')
    parser.add_argument('--mark-threshold', type=float, default=1.0, help="Threshold value to mark with a vertical line in the 2D panel (default 1.0, matching the paper's Table 12/13 definition).")
    parser.add_argument('--out-csv', default='dead_zone_threshold_sweep.csv')
    parser.add_argument('--out-plot', default='dead_zone_threshold_sweep.png')
    args = parser.parse_args()
    ratio_maps = compute_ratio_maps(args.hits_dir, args.barrel_hits_dir)
    per_run, agg = sweep(ratio_maps, n_thresholds=args.n_thresholds)
    agg.to_csv(args.out_csv, index=False)
    print(f'\nSaved sweep table ({len(agg)} rows) -> {args.out_csv}')
    print(f'(full resolution: {args.n_thresholds} threshold points from 0 to 1, step {1.0 / (args.n_thresholds - 1):.4f} -- every value is in the CSV, this is just a coarser preview)')
    nearest_tenth = np.round(agg['threshold'] / 0.1) * 0.1
    preview = agg[np.isclose(agg['threshold'], nearest_tenth, atol=1e-09)]
    print(preview.to_string(index=False))
    display_names = {s: clean_label(s) for s in agg['structure'].unique()}
    plot_sweep(agg, args.out_plot, mark_threshold=args.mark_threshold, display_names=display_names)
if __name__ == '__main__':
    main()