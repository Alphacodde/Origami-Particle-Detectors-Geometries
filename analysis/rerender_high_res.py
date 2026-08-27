import os
import sys
import glob
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
OUT_DIR = os.path.abspath('Test Render')
os.makedirs(OUT_DIR, exist_ok=True)
DPI = 300
print(f"=== High-Resolution Figure Re-Render (300 DPI) -> '{OUT_DIR}' ===")
print('\n[1/5] Re-rendering Experiment 2 multiple-scattering plots (300 DPI)...')
import analyze_exp2
exp2_roots = sorted(glob.glob('results_exp2/*.root'))
if exp2_roots:
    summaries_exp2 = []
    for rpath in exp2_roots:
        data = analyze_exp2.load_run(rpath)
        name = analyze_exp2.resolve_geometry_name(rpath, data)
        mom = analyze_exp2.validate_momentum(rpath, data)
        if mom:
            s = analyze_exp2.summarize_scattering(name, mom['claimed_GeV'], rpath, data)
            if s:
                summaries_exp2.append(s)
else:
    import json
    with open('exp2_scattering_summary.json') as f:
        summaries_exp2 = json.load(f)
by_geom = {}
for s in summaries_exp2:
    if s:
        by_geom.setdefault(s['geometry'], []).append(s)
fig, ax = plt.subplots(figsize=(8, 6))
colors = plt.cm.tab10.colors
for i, (geom, points) in enumerate(sorted(by_geom.items())):
    points = sorted(points, key=lambda pt: pt['momentum_GeV'])
    p_vals = [pt['momentum_GeV'] for pt in points]
    sim_mrad = [pt['sim_theta0']['gaussian_fit_mrad'] if pt['sim_theta0']['gaussian_fit_mrad'] is not None else pt['sim_theta0']['rms_mrad'] for pt in points]
    highland_mrad = [pt['highland_theta0_mean_mrad'] for pt in points]
    color = colors[i % len(colors)]
    ax.plot(p_vals, sim_mrad, 'o-', color=color, label=f'{geom} (GEANT4)')
    ax.plot(p_vals, highland_mrad, '--', color=color, alpha=0.6, label=f'{geom} (Highland)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Momentum p (GeV/c)')
ax.set_ylabel('$\\theta_0$ (mrad, projected-angle equivalent)')
ax.set_title('Experiment 2: Multiple scattering $\\theta_0$ vs. momentum')
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
p1 = os.path.join(OUT_DIR, 'fig_scattering_theta0.png')
plt.savefig(p1, dpi=DPI)
plt.close(fig)
print(f'  -> Saved {p1}')
fig, ax = plt.subplots(figsize=(8, 5))
for i, (geom, points) in enumerate(sorted(by_geom.items())):
    points = sorted(points, key=lambda pt: pt['momentum_GeV'])
    p_vals = [pt['momentum_GeV'] for pt in points if pt.get('sim_over_highland_ratio') is not None]
    ratios = [pt['sim_over_highland_ratio'] for pt in points if pt.get('sim_over_highland_ratio') is not None]
    if ratios:
        ax.plot(p_vals, ratios, 'o-', color=colors[i % len(colors)], label=geom)
ax.axhline(1.0, color='black', linewidth=0.8, linestyle=':')
ax.set_xscale('log')
ax.set_xlabel('Momentum p (GeV/c)')
ax.set_ylabel('$\\theta_0^{sim} / \\theta_0^{Highland}$')
ax.set_title('Experiment 2: GEANT4 vs. Highland agreement')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
plt.tight_layout()
p2 = os.path.join(OUT_DIR, 'fig_scattering_ratio.png')
plt.savefig(p2, dpi=DPI)
plt.close(fig)
print(f'  -> Saved {p2}')
print('\n[2/5] Re-rendering Experiment 3 vertex smear plots (300 DPI)...')
import analyze_exp3
exp3_roots = sorted(glob.glob('results_exp3/*.root'))
if exp3_roots:
    summaries_exp3 = []
    for rpath in exp3_roots:
        data = analyze_exp3.load_run(rpath)
        name = analyze_exp3.resolve_geometry_name(rpath, data)
        sigma = analyze_exp3.validate_vertex_smear(rpath, data)
        if sigma:
            s = analyze_exp3.summarize_robustness(name, sigma['claimed_mm'], rpath, data)
            if s:
                summaries_exp3.append(s)
    p3 = os.path.join(OUT_DIR, 'fig_vertex_smear_sigma.png')
    by_geom3 = {}
    for s in summaries_exp3:
        if s:
            by_geom3.setdefault(s['geometry'], []).append(s)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    panels = [('acceptance_fraction', 'Acceptance (fraction of events hitting sensor)'), ('mean_path_length_X0', 'Mean path length ($X/X_0$)'), ('mean_local_incidence_deg', 'Mean local incidence angle (deg)'), ('mean_nHitsPrimary', 'Mean primary hits per event')]
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for i, (geom, points) in enumerate(sorted(by_geom3.items())):
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
    plt.savefig(p3, dpi=DPI)
    plt.close(fig)
    print(f'  -> Saved {p3}')
    p4 = os.path.join(OUT_DIR, 'fig_vertex_smear_displacement.png')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    for i, (geom, points) in enumerate(sorted(by_geom3.items())):
        r_all = np.concatenate([p['_r_all_mm'] for p in points])
        hit_all = np.concatenate([p['_hit_mask'] for p in points])
        r_hit = np.concatenate([p['_r_hit_mm'] for p in points])
        path_X0_hit = np.concatenate([p['_path_X0_hit'] for p in points])
        color = colors[i % len(colors)]
        acc_binned = analyze_exp3._bin_by_displacement(r_all, hit_all, agg='fraction')
        ax1.plot(acc_binned['bin_centers_mm'], acc_binned['values'], 'o-', color=color, label=geom)
        path_binned = analyze_exp3._bin_by_displacement(r_hit, path_X0_hit, agg='mean')
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
    plt.savefig(p4, dpi=DPI)
    plt.close(fig)
    print(f'  -> Saved {p4}')
print('\n[3/5] Re-rendering Pareto analysis & parameter sweep figures (300 DPI)...')
import pareto
json_multi = 'JSONs, Miscellenous/geant4_final_result_summaries_multi_run_.json'
csv_dz = 'results_diff_geom/original_four_dead_zone_summary_aggregate.csv'
n_sweep_json = 'JSONs, Miscellenous/geant4_N_Sweep_summaries.json'
theta_sweep_json = 'JSONs, Miscellenous/geant4_N_Sweep_results_.json'
geant4_cands = pareto.load_geant4_multi_run(json_multi)
deadzone_d = pareto.load_deadzone_csv(csv_dz)
candidates = pareto.merge_datasets(geant4_cands, deadzone_d)
candidates = pareto.compute_pareto_dominance(candidates)
n_sw = pareto.load_n_sweep(n_sweep_json)
theta_sw = pareto.load_theta_sweep(theta_sweep_json)
pareto.plot_all(candidates, OUT_DIR, n_sweep=n_sw, theta_sweep=theta_sw)
rename_map = {'fig_acceptance_vs_deadzone.png': 'fig6a_acceptance_vs_deadzone.png', 'fig_acceptance_vs_x0.png': 'fig6b_acceptance_vs_x0.png', 'fig_x0_vs_deadzone.png': 'fig7_x0_vs_deadzone.png', 'fig_tradeoff_matrix_3panel.png': 'fig8_tradeoff_matrix.png', 'fig_parallel_coords_3obj.png': 'fig9_parallel_coords.png', 'fig_n_sweep.png': 'fig4_n_sweep.png', 'fig_theta_sweep.png': 'fig5_theta_sweep.png'}
for src_name, dst_name in rename_map.items():
    src_f = os.path.join(OUT_DIR, src_name)
    dst_f = os.path.join(OUT_DIR, dst_name)
    if os.path.exists(src_f):
        if os.path.exists(dst_f):
            os.remove(dst_f)
        os.rename(src_f, dst_f)
        print(f'  -> Renamed {src_name} to {dst_name}')
print('\n[4/5] Re-rendering Dead-Zone 3D Response Surface (300 DPI)...')
import deadzone_threshold_sweep
dz_sweep_csv = 'results_diff_geom/dead_zone_threshold_sweep.csv'
if os.path.exists(dz_sweep_csv):
    agg_dz = pd.read_csv(dz_sweep_csv)
    display_names = {s: deadzone_threshold_sweep.clean_label(s) for s in agg_dz['structure'].unique()}
    p_dz = os.path.join(OUT_DIR, 'fig_deadzone_threshold_sweep.png')
    deadzone_threshold_sweep.plot_sweep(agg_dz, p_dz, mark_threshold=1.0, display_names=display_names)
    print(f'  -> Saved {p_dz}')
print('\n[5/5] Re-rendering 2D Spatial Hit-Density Maps (300 DPI)...')
import deadzone_map
hits_dir = 'results_diff_geom'
barrel_hits_dir = hits_dir
frames = [deadzone_map.load_hits(p) for p in sorted(glob.glob(f'{hits_dir}/*_hits.csv'))]
all_hits = pd.concat(frames, ignore_index=True)
fold_hits = all_hits[~all_hits['structure'].apply(lambda s: any((deadzone_map.is_structure(s, p) for p in deadzone_map.EXCLUDED_STRUCTURE_PREFIXES)))]
barrel_hits = all_hits[all_hits['structure'].apply(lambda s: deadzone_map.is_structure(s, 'barrel'))]
map_rename = {'kresling': 'fig12_map_kresling.png', 'miura': 'fig10_map_miura.png', 'yoshimura': 'fig11_map_yoshimura.png'}
for struct_name, group in fold_hits.groupby('structure'):
    run_group = group[group['run'] == '0_4']
    if run_group.empty:
        run_group = group[group['run'] == group['run'].iloc[0]]
    matched_barrel = barrel_hits[barrel_hits['run'] == run_group['run'].iloc[0]]
    if matched_barrel.empty:
        matched_barrel = barrel_hits
    ratio, mean_pathX0, theta_edges, z_edges = deadzone_map.build_maps(run_group, matched_barrel)
    base_key = struct_name.lower().split('_')[0]
    out_fname = map_rename.get(base_key, f'map_{struct_name}.png')
    out_fpath = os.path.join(OUT_DIR, out_fname)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    im0 = axes[0].imshow(ratio.T, origin='lower', aspect='auto', extent=[theta_edges[0], theta_edges[-1], z_edges[0], z_edges[-1]], cmap='RdBu_r', vmin=0, vmax=2)
    axes[0].set_title(f'{deadzone_threshold_sweep.clean_label(struct_name)}: Hit-density ratio vs. barrel\n(< 1 = dead zone, > 1 = overlap-ish)')
    axes[0].set_xlabel('$\\theta$ (deg)')
    axes[0].set_ylabel('z (mm)')
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(mean_pathX0.T, origin='lower', aspect='auto', extent=[theta_edges[0], theta_edges[-1], z_edges[0], z_edges[-1]], cmap='viridis')
    axes[1].set_title(f'{deadzone_threshold_sweep.clean_label(struct_name)}: Mean path length ($X/X_0$) per bin\n(elevated = overlap seam)')
    axes[1].set_xlabel('$\\theta$ (deg)')
    axes[1].set_ylabel('z (mm)')
    fig.colorbar(im1, ax=axes[1])
    plt.tight_layout()
    plt.savefig(out_fpath, dpi=DPI)
    plt.close(fig)
    print(f'  -> Saved {out_fpath}')
print(f"\n{'=' * 70}")
print(f'All high-resolution figures successfully rendered into: {OUT_DIR}')
for f in sorted(os.listdir(OUT_DIR)):
    p = os.path.join(OUT_DIR, f)
    print(f'  {f:42s}  {os.path.getsize(p):>10,} bytes')
print(f"{'=' * 70}")