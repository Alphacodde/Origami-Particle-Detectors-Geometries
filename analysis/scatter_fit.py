import argparse
import csv
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import uproot
from scipy.optimize import minimize
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
try:
    from tqdm import tqdm
    _HAVE_TQDM = True
except ImportError:
    _HAVE_TQDM = False
COLORS = dict(paper='#FFFFFF', ink='#1B2430', ink_soft='#5B6472', data='#2B3A55', core='#D97706', core_soft='#F3C57D', tail='#B91C1C', tail_wash='#FEE2E2', grid='#E5E7EB')
FONT_DISPLAY = 'Lora'
FONT_BODY = 'DejaVu Sans'

def _apply_house_style():
    plt.rcParams.update({'figure.facecolor': COLORS['paper'], 'axes.facecolor': COLORS['paper'], 'savefig.facecolor': COLORS['paper'], 'axes.edgecolor': COLORS['ink_soft'], 'axes.labelcolor': COLORS['ink'], 'axes.linewidth': 1.0, 'text.color': COLORS['ink'], 'xtick.color': COLORS['ink_soft'], 'ytick.color': COLORS['ink_soft'], 'font.family': FONT_BODY, 'font.size': 10.5, 'axes.titlesize': 13, 'axes.titleweight': 'bold', 'grid.color': COLORS['grid'], 'grid.linewidth': 0.7, 'legend.frameon': False, 'legend.fontsize': 9, 'mathtext.fontset': 'cm'})
PION_MASS_GEV = 0.13957
X0_SILICON_MM = 93.7
X0_KAPTON_MM = 285.8
FILENAME_RE = re.compile('\n    (?P<tag_prefix>.+?)_\n    fold(?P<fold>[0-9p]+)_\n    p(?P<momentum>[0-9p]+)GeV_\n    vtx(?P<vtx>[0-9p]+)mm_\n    tag(?P<tag>[A-Za-z0-9]+)_\n    run(?P<run>[0-9_]+)\n    ', re.VERBOSE)

def _p_to_float(token: str) -> float:
    return float(token.replace('p', '.'))

def parse_run_metadata(filepath: Path) -> dict:
    stem = filepath.stem
    m = FILENAME_RE.search(stem)
    if not m:
        print(f'  [warn] filename does not match expected pattern: {filepath.name}')
        return dict(momentum_gev=np.nan, fold=np.nan, vtx_mm=np.nan, tag='unknown', run='unknown', structure='unknown', label=stem)
    momentum = _p_to_float(m.group('momentum'))
    fold = _p_to_float(m.group('fold'))
    vtx = _p_to_float(m.group('vtx'))
    tag = m.group('tag')
    run = m.group('run')
    tag_prefix = m.group('tag_prefix')
    structure = tag_prefix
    for lead in (f'{tag}_origami_', f'{tag}_'):
        if structure.startswith(lead):
            structure = structure[len(lead):]
            break
    if structure.startswith('origami_'):
        structure = structure[len('origami_'):]
    structure = structure or 'unknown'
    label = f'{tag}_{structure}_run{run}_p{momentum:g}GeV'
    return dict(momentum_gev=momentum, fold=fold, vtx_mm=vtx, tag=tag, run=run, structure=structure, label=label)

def highland_theta0_rad(x_over_x0: np.ndarray, pc_beta_mev: float, z: float=1.0) -> np.ndarray:
    x = np.clip(x_over_x0, 1e-12, None)
    return 13.6 / pc_beta_mev * z * np.sqrt(x) * (1.0 + 0.038 * np.log(x * z ** 2))

def beta_of_momentum(p_gev: float, mass_gev: float=PION_MASS_GEV) -> float:
    e = np.sqrt(p_gev ** 2 + mass_gev ** 2)
    return p_gev / e

def rayleigh_pdf(theta, sigma):
    return theta / sigma ** 2 * np.exp(-theta ** 2 / (2 * sigma ** 2))

def rayleigh_cdf(theta, sigma):
    return 1.0 - np.exp(-theta ** 2 / (2 * sigma ** 2))

def powerlaw_pdf(theta, theta_min, alpha):
    return (alpha - 1) / theta_min * (theta / theta_min) ** (-alpha)

def powerlaw_cdf(theta, theta_min, alpha):
    theta = np.clip(theta, theta_min, None)
    return 1.0 - (theta / theta_min) ** (-(alpha - 1))

def mixture_pdf(theta, sigma, f, theta_min, alpha):
    return (1 - f) * rayleigh_pdf(theta, sigma) + f * np.where(theta > theta_min, powerlaw_pdf(theta, theta_min, alpha), 0.0)

def mixture_cdf(theta, sigma, f, theta_min, alpha):
    return np.where(theta < theta_min, (1 - f) * rayleigh_cdf(theta, sigma), (1 - f) * rayleigh_cdf(theta_min, sigma) + f * powerlaw_cdf(theta, theta_min, alpha))

def _neg_log_likelihood(params, data):
    sigma, f, alpha = params
    if sigma <= 0 or f < 0 or f >= 1 or (alpha <= 1):
        return 10000000000.0
    theta_min = 3 * sigma
    pdf = mixture_pdf(data, sigma, f, theta_min, alpha)
    pdf = np.clip(pdf, 1e-300, None)
    return -np.sum(np.log(pdf))

def fit_mixture(scatter_deg: np.ndarray, sigma_guess: float):
    data = scatter_deg[scatter_deg > 0]
    x0 = [sigma_guess, 0.02, 3.0]
    res = minimize(_neg_log_likelihood, x0, args=(data,), method='Nelder-Mead', options={'xatol': 1e-09, 'fatol': 1e-07, 'maxiter': 8000})
    sigma_fit, f_fit, alpha_fit = res.x
    theta_min = 3 * sigma_fit
    ks_stat, ks_p = stats.kstest(data, lambda x: mixture_cdf(x, sigma_fit, f_fit, theta_min, alpha_fit))
    core = data[data <= theta_min]
    tail = data[data > theta_min]
    ms_core = np.sum(core ** 2)
    ms_tail = np.sum(tail ** 2)
    ms_total = ms_core + ms_tail
    return dict(converged=bool(res.success), sigma_fit_deg=float(sigma_fit), tail_fraction=float(f_fit), tail_alpha=float(alpha_fit), theta_min_deg=float(theta_min), ks_stat=float(ks_stat), ks_pvalue=float(ks_p), n_events_fit=int(len(data)), n_core=int(len(core)), n_tail=int(len(tail)), tail_ms_fraction=float(ms_tail / ms_total) if ms_total > 0 else np.nan)

def process_file(filepath: Path, out_dir: Path, tree_name: str='PionEvents') -> dict:
    meta = parse_run_metadata(filepath)
    with uproot.open(filepath) as f:
        keys = f.keys()
        key = tree_name if any((k.startswith(tree_name) for k in keys)) else keys[0]
        tree = f[key]
        arrs = tree.arrays(library='np')
    hit = arrs['hitDetector'] == 1
    n_total = len(hit)
    n_hit = int(hit.sum())
    scat = arrs['scatterAngleDeg'][hit]
    si = arrs['siliconPathLength_mm'][hit]
    kap = arrs['kaptonPathLength_mm'][hit]
    p_gev = meta['momentum_gev']
    if np.isnan(p_gev):
        raise ValueError(f'Could not determine beam momentum from filename {filepath.name}; pass --momentum-gev to override for non-conforming filenames.')
    beta = beta_of_momentum(p_gev)
    pc_beta_mev = beta * p_gev * 1000.0
    x_over_x0 = si / X0_SILICON_MM + kap / X0_KAPTON_MM
    theta0_rad = highland_theta0_rad(x_over_x0, pc_beta_mev)
    theta0_deg = np.degrees(theta0_rad)
    highland_mean_deg = float(theta0_deg.mean())
    fit = fit_mixture(scat, sigma_guess=highland_mean_deg)
    hist = compute_binned_histogram(scat)
    chi2_red = make_plot(scat, highland_mean_deg, fit, meta, out_dir, hist=hist)
    summary = dict(file=filepath.name, **meta, beta=beta, n_total=n_total, n_hit=n_hit, hit_fraction=n_hit / n_total if n_total else np.nan, scat_mean_deg=float(scat.mean()), scat_median_deg=float(np.median(scat)), scat_rms_deg=float(np.sqrt(np.mean(scat ** 2))), highland_theta0_deg=highland_mean_deg, **fit, binned_chi2_per_dof=chi2_red)
    summary['sigma_fit_over_highland'] = summary['sigma_fit_deg'] / highland_mean_deg
    summary['_hist'] = hist
    return summary

def compute_binned_histogram(scat_deg: np.ndarray, n_bins: int=40) -> dict:
    data = scat_deg[scat_deg > 0]
    lo = max(np.percentile(data, 0.1), 1e-06)
    hi = data.max()
    bins = np.logspace(np.log10(lo), np.log10(hi), n_bins)
    raw_counts, edges = np.histogram(data, bins=bins)
    counts, edges = np.histogram(data, bins=bins, density=True)
    centers = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    keep = counts > 0
    n_total_hits = len(data)
    alpha_ci = 1 - 0.6827
    lower_count = np.where(raw_counts > 0, stats.chi2.ppf(alpha_ci / 2, 2 * raw_counts) / 2, 0.0)
    upper_count = stats.chi2.ppf(1 - alpha_ci / 2, 2 * (raw_counts + 1)) / 2
    err_lo = (raw_counts - lower_count) / (n_total_hits * widths)
    err_hi = (upper_count - raw_counts) / (n_total_hits * widths)
    return dict(centers=centers, counts=counts, err_lo=err_lo, err_hi=err_hi, raw_counts=raw_counts, keep=keep, lo=lo, hi=hi)

def make_plot(scat_deg, highland_theta0_deg, fit, meta, out_dir: Path, hist: dict=None):
    _apply_house_style()
    if hist is None:
        hist = compute_binned_histogram(scat_deg)
    centers = hist['centers']
    counts = hist['counts']
    err_lo = hist['err_lo']
    err_hi = hist['err_hi']
    keep = hist['keep']
    lo, hi = (hist['lo'], hist['hi'])
    ymin_data = counts[keep].min()
    ymax_data = counts[keep].max()
    sigma_fit = fit['sigma_fit_deg']
    f_fit = fit['tail_fraction']
    alpha_fit = fit['tail_alpha']
    theta_min = fit['theta_min_deg']
    model_vals = mixture_pdf(centers, sigma_fit, f_fit, theta_min, alpha_fit)
    highland_vals = rayleigh_pdf(centers, highland_theta0_deg)
    density_err = 0.5 * (err_lo + err_hi)
    resid = (counts[keep] - model_vals[keep]) / density_err[keep]
    n_free_params = 3
    dof = max(keep.sum() - n_free_params, 1)
    chi2_red = float(np.sum(resid ** 2) / dof)
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    ax.axvspan(theta_min, hi * 1.5, color=COLORS['tail_wash'], zorder=0)
    ax.plot(centers, highland_vals, '--', color=COLORS['tail'], lw=1.8, alpha=0.85, zorder=2, label=f'Highland only ($\\theta_0$ = {highland_theta0_deg:.3g}°)')
    ax.plot(centers, model_vals, '-', color=COLORS['core'], lw=2.2, zorder=3, label='Gaussian core + Rutherford tail')
    ax.errorbar(centers[keep], counts[keep], yerr=[err_lo[keep], err_hi[keep]], fmt='o', mfc=COLORS['paper'], mec=COLORS['data'], mew=1.4, ms=5, ecolor=COLORS['data'], elinewidth=1.1, capsize=2.5, capthick=1.1, alpha=0.95, zorder=4, label='Simulated data (Poisson $1\\sigma$)')
    ax.axvline(theta_min, color=COLORS['ink_soft'], ls=':', lw=1, zorder=1)
    ax.text(theta_min, ymax_data * 2.2, '  tail region', color=COLORS['tail'], fontsize=8.5, style='italic', va='bottom', family=FONT_BODY)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Scatter angle (degrees)', fontsize=10.5)
    ax.set_ylabel('Probability density', fontsize=10.5)
    ax.set_ylim(ymin_data / 5, ymax_data * 6)
    ax.set_xlim(lo * 0.7, hi * 2.2)
    ax.set_title(meta['label'], fontsize=14, family=FONT_DISPLAY, weight='bold', color=COLORS['ink'], pad=30)
    ax.text(0.0, 1.045, f"p = {meta['momentum_gev']:g} GeV/c   ·   barrel fold {meta['fold']:g}   ·   vtx {meta['vtx_mm']:g} mm", transform=ax.transAxes, fontsize=9.5, color=COLORS['ink_soft'], family=FONT_BODY)
    ax.grid(alpha=0.6, which='major', lw=0.7)
    ax.grid(alpha=0.25, which='minor', lw=0.4)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.legend(loc='lower left', fontsize=9, handlelength=2.2)
    stats_text = f'$\\sigma_{{fit}}$ = {sigma_fit:.3g}°   ($\\sigma_{{fit}}$/Highland = {sigma_fit / highland_theta0_deg:.2f})\ntail fraction f = {100 * f_fit:.2f}%\ntail index $\\alpha$ = {alpha_fit:.2f}  (Rutherford: 3)\nbinned $\\chi^2$/dof = {chi2_red:.2f}'
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes, ha='right', va='top', fontsize=8.7, family=FONT_BODY, color=COLORS['ink'], bbox=dict(boxstyle='round,pad=0.5', facecolor=COLORS['paper'], edgecolor=COLORS['grid'], linewidth=1))
    fig.tight_layout()
    out_path = out_dir / f"{meta['label']}_scatter_fit.png"
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    return chi2_red

def make_geometry_comparison_plot(rows: list, out_dir: Path):
    rows = [r for r in rows if not np.isnan(r['momentum_gev']) and r['structure'] != 'unknown']
    if not rows:
        print('No runs with valid momentum + structure; skipping geometry comparison plot.')
        return
    structures = sorted(set((r['structure'] for r in rows)))
    if len(structures) < 2:
        print('Only one geometry present; use make_trend_plots instead of the comparison plot.')
        return
    _apply_house_style()
    CANONICAL_GEO_COLORS = {'barrel_reference': '#1f77b4', 'kresling_deployed': '#ff7f0e', 'miura_deployed': '#2ca02c', 'yoshimura_deployed': '#d62728'}
    fallback_palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    geo_color = {s: CANONICAL_GEO_COLORS.get(s, fallback_palette[i % len(fallback_palette)]) for i, s in enumerate(structures)}
    marker_kw = dict(ms=6, mew=1.4, mfc=COLORS['paper'])
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5))
    ax = axes[0, 0]
    for s in structures:
        sub = sorted([r for r in rows if r['structure'] == s], key=lambda r: r['momentum_gev'])
        p = [r['momentum_gev'] for r in sub]
        sigma_fit = [r['sigma_fit_deg'] for r in sub]
        ax.plot(p, sigma_fit, 'o-', color=geo_color[s], label=s.replace('_', ' '), **marker_kw)
    ref = sorted(rows, key=lambda r: r['momentum_gev'])
    ax.plot([r['momentum_gev'] for r in ref], [r['highland_theta0_deg'] for r in ref], '--', color=COLORS['ink_soft'], lw=1.5, zorder=1, label='Highland $\\theta_0$ (reference)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('Fitted core $\\sigma$ (degrees)')
    ax.set_title('Core width by geometry', family=FONT_DISPLAY, fontsize=13)
    ax.legend(loc='best', fontsize=8.5)
    ax = axes[0, 1]
    for s in structures:
        sub = sorted([r for r in rows if r['structure'] == s], key=lambda r: r['momentum_gev'])
        p = [r['momentum_gev'] for r in sub]
        ratio = [r['sigma_fit_over_highland'] for r in sub]
        ax.plot(p, ratio, 'o-', color=geo_color[s], label=s.replace('_', ' '), **marker_kw)
    ax.axhline(1.0, color=COLORS['ink_soft'], ls=':', lw=1.2)
    ax.set_xscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('$\\sigma_{fit}$ / Highland $\\theta_0$')
    ax.set_title('Core-fit agreement with Highland', family=FONT_DISPLAY, fontsize=13)
    ax.legend(loc='best', fontsize=8.5)
    ax = axes[1, 0]
    for s in structures:
        sub = sorted([r for r in rows if r['structure'] == s], key=lambda r: r['momentum_gev'])
        p = [r['momentum_gev'] for r in sub]
        tail_ms = [100 * r['tail_ms_fraction'] for r in sub]
        ax.plot(p, tail_ms, 'o-', color=geo_color[s], label=s.replace('_', ' '), **marker_kw)
    ax.set_xscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('Tail share of $\\Sigma\\theta^2$ (%)')
    ax.set_title('Tail contribution by geometry', family=FONT_DISPLAY, fontsize=13)
    ax.legend(loc='best', fontsize=8.5)
    ax = axes[1, 1]
    ax.axhspan(2.85, 3.15, color=COLORS['tail_wash'], zorder=0)
    for s in structures:
        sub = sorted([r for r in rows if r['structure'] == s], key=lambda r: r['momentum_gev'])
        p = [r['momentum_gev'] for r in sub]
        alpha = [r['tail_alpha'] for r in sub]
        ax.plot(p, alpha, 'o-', color=geo_color[s], zorder=2, label=s.replace('_', ' '), **marker_kw)
    ax.axhline(3.0, color=COLORS['tail'], ls=':', lw=1.2, zorder=1, label='Rutherford ($\\alpha$=3)')
    ax.set_xscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('Tail power index $\\alpha$')
    ax.set_title('Tail slope by geometry', family=FONT_DISPLAY, fontsize=13)
    ax.legend(loc='best', fontsize=8.5)
    for ax in axes.flat:
        ax.grid(alpha=0.6, which='major', lw=0.7)
        ax.grid(alpha=0.25, which='minor', lw=0.4)
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
    fig.suptitle('Scattering behavior across origami fold geometries', fontsize=16.5, family=FONT_DISPLAY, weight='bold', color=COLORS['ink'], y=1.01)
    fig.tight_layout()
    out_path = out_dir / 'geometry_comparison.png'
    fig.savefig(out_path, dpi=170, bbox_inches='tight')
    plt.close(fig)
    print(f'saved geometry comparison plot -> {out_path}')
_MOMENTUM_PALETTE = ['#1B2430', '#2B5C8A', '#0D9488', '#D97706', '#B91C1C', '#7C3AED']

def make_overlay_plots(rows: list, out_dir: Path):
    rows = [r for r in rows if not np.isnan(r['momentum_gev']) and r['structure'] != 'unknown']
    have_hist = [r for r in rows if '_hist' in r]
    if len(have_hist) < len(rows):
        print(f'  [note] {len(rows) - len(have_hist)} run(s) missing binned histogram data; skipping those in the overlay plots.')
    rows = have_hist
    if not rows:
        print('No runs with histogram data available; skipping multi-momentum overlay plots.')
        return
    _apply_house_style()
    for structure in sorted(set((r['structure'] for r in rows))):
        sub = sorted([r for r in rows if r['structure'] == structure], key=lambda r: r['momentum_gev'])
        if len(sub) < 2:
            continue
        colors = [_MOMENTUM_PALETTE[i % len(_MOMENTUM_PALETTE)] for i in range(len(sub))]
        fig, ax = plt.subplots(figsize=(7.5, 5.8))
        for r, c in zip(sub, colors):
            h = r['_hist']
            keep = h['keep']
            ax.errorbar(h['centers'][keep], h['counts'][keep], yerr=[h['err_lo'][keep], h['err_hi'][keep]], fmt='o', ms=4.5, mfc=COLORS['paper'], mec=c, mew=1.3, ecolor=c, elinewidth=1.0, capsize=2, alpha=0.9, label=f"{r['momentum_gev']:g} GeV/c", zorder=3)
            model_vals = mixture_pdf(h['centers'], r['sigma_fit_deg'], r['tail_fraction'], r['theta_min_deg'], r['tail_alpha'])
            ax.plot(h['centers'], model_vals, '-', color=c, lw=1.6, alpha=0.75, zorder=2)
            highland_vals = rayleigh_pdf(h['centers'], r['highland_theta0_deg'])
            ymin_plot = h['counts'][keep].min() / 5
            highland_visible = highland_vals > ymin_plot
            ax.plot(h['centers'][highland_visible], highland_vals[highland_visible], ':', color=c, lw=1.8, alpha=0.9, zorder=2)
        style_handles = [plt.Line2D([0], [0], color=COLORS['ink_soft'], lw=1.8, ls='-', label='Core + tail fit'), plt.Line2D([0], [0], color=COLORS['ink_soft'], lw=1.8, ls=':', label='Highland only (theory)')]
        style_legend = ax.legend(handles=style_handles, loc='upper right', fontsize=8, frameon=True, facecolor=COLORS['paper'], edgecolor=COLORS['grid'])
        ax.add_artist(style_legend)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Scatter angle (degrees)')
        ax.set_ylabel('Probability density')
        ax.set_title(f"{structure.replace('_', ' ')} — all momenta", fontsize=14, family=FONT_DISPLAY, weight='bold', color=COLORS['ink'], pad=14)
        ax.legend(title='Momentum', loc='lower left', fontsize=8.5, title_fontsize=8.5)
        ax.add_artist(style_legend)
        ax.grid(alpha=0.6, which='major', lw=0.7)
        ax.grid(alpha=0.25, which='minor', lw=0.4)
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
        fig.tight_layout()
        out_path = out_dir / f'{structure}_overlay_raw.png'
        fig.savefig(out_path, dpi=170)
        plt.close(fig)
        print(f'saved raw-angle overlay -> {out_path}')
        fig, ax = plt.subplots(figsize=(7.5, 5.8))
        all_x, all_y = ([], [])
        for r, c in zip(sub, colors):
            h = r['_hist']
            keep = h['keep']
            theta0 = r['highland_theta0_deg']
            x = h['centers'][keep] / theta0
            y = h['counts'][keep] * theta0
            ax.errorbar(x, y, yerr=[h['err_lo'][keep] * theta0, h['err_hi'][keep] * theta0], fmt='o', ms=4.5, mfc=COLORS['paper'], mec=c, mew=1.3, ecolor=c, elinewidth=1.0, capsize=2, alpha=0.9, label=f"{r['momentum_gev']:g} GeV/c", zorder=3)
            all_x.append(x)
            all_y.append(y)
        all_x = np.concatenate(all_x)
        all_y = np.concatenate(all_y)
        x_ref = np.logspace(np.log10(all_x.min() * 0.7), np.log10(all_x.max() * 1.3), 200)
        ax.plot(x_ref, rayleigh_pdf(x_ref, 1.0), '--', color=COLORS['ink_soft'], lw=1.8, zorder=1, label='Rayleigh(1) reference')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Scatter angle / Highland $\\theta_0$')
        ax.set_ylabel('Probability density $\\times\\ \\theta_0$')
        ax.set_ylim(all_y.min() / 5, all_y.max() * 6)
        ax.set_xlim(all_x.min() * 0.7, all_x.max() * 1.4)
        ax.set_title(f"{structure.replace('_', ' ')} — Highland-normalized", fontsize=14, family=FONT_DISPLAY, weight='bold', color=COLORS['ink'], pad=14)
        ax.legend(title='Momentum', loc='lower left', fontsize=8.5, title_fontsize=8.5)
        ax.grid(alpha=0.6, which='major', lw=0.7)
        ax.grid(alpha=0.25, which='minor', lw=0.4)
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
        fig.tight_layout()
        out_path = out_dir / f'{structure}_overlay_scaled.png'
        fig.savefig(out_path, dpi=170)
        plt.close(fig)
        print(f'saved Highland-normalized overlay -> {out_path}')

def make_trend_plots(rows: list, out_dir: Path):
    rows_valid = [r for r in rows if not np.isnan(r['momentum_gev'])]
    structures = set((r.get('structure', 'unknown') for r in rows_valid))
    if len(structures - {'unknown'}) >= 2:
        make_geometry_comparison_plot(rows, out_dir)
        return
    rows = rows_valid
    if len(rows) < 2:
        print('Fewer than 2 runs with valid momentum; skipping trend plots.')
        return
    _apply_house_style()
    rows = sorted(rows, key=lambda r: r['momentum_gev'])
    p = np.array([r['momentum_gev'] for r in rows])
    ratio = np.array([r['sigma_fit_over_highland'] for r in rows])
    tail_f = np.array([r['tail_fraction'] for r in rows])
    alpha = np.array([r['tail_alpha'] for r in rows])
    tail_ms = np.array([r['tail_ms_fraction'] for r in rows])
    highland = np.array([r['highland_theta0_deg'] for r in rows])
    sigma_fit = np.array([r['sigma_fit_deg'] for r in rows])
    marker_kw = dict(ms=6, mew=1.4, mfc=COLORS['paper'])
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5))
    ax = axes[0, 0]
    ax.plot(p, highland, 'o--', color=COLORS['tail'], label='Highland $\\theta_0$', **marker_kw)
    ax.plot(p, sigma_fit, 'o-', color=COLORS['core'], label='Fitted core $\\sigma$', **marker_kw)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('Angle (degrees)')
    ax.set_title('Core width vs. momentum', family=FONT_DISPLAY, fontsize=12.5)
    ax.legend(loc='best')
    ax = axes[0, 1]
    ax.plot(p, ratio, 'o-', color=COLORS['data'], **marker_kw)
    ax.axhline(1.0, color=COLORS['ink_soft'], ls=':', lw=1.2)
    ax.set_xscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('$\\sigma_{fit}$ / Highland $\\theta_0$')
    ax.set_title('Core-fit agreement with Highland', family=FONT_DISPLAY, fontsize=12.5)
    ax = axes[1, 0]
    ax.plot(p, 100 * tail_f, 'o-', color=COLORS['core'], label='Tail fraction (events)', **marker_kw)
    ax.plot(p, 100 * tail_ms, 'o--', color=COLORS['tail'], label='Tail share of $\\Sigma\\theta^2$', **marker_kw)
    ax.set_xscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('Percent (%)')
    ax.set_title('Tail contribution vs. momentum', family=FONT_DISPLAY, fontsize=12.5)
    ax.legend(loc='best')
    ax = axes[1, 1]
    ax.axhspan(2.85, 3.15, color=COLORS['tail_wash'], zorder=0)
    ax.plot(p, alpha, 'o-', color=COLORS['data'], zorder=2, **marker_kw)
    ax.axhline(3.0, color=COLORS['tail'], ls=':', lw=1.2, zorder=1, label='Rutherford expectation ($\\alpha$ = 3)')
    ax.set_xscale('log')
    ax.set_xlabel('Momentum (GeV/c)')
    ax.set_ylabel('Tail power index $\\alpha$')
    ax.set_title('Tail slope vs. momentum', family=FONT_DISPLAY, fontsize=12.5)
    ax.legend(loc='best')
    for ax in axes.flat:
        ax.grid(alpha=0.6, which='major', lw=0.7)
        ax.grid(alpha=0.25, which='minor', lw=0.4)
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
    fig.suptitle('Highland core vs. Rutherford tail across the momentum scan', fontsize=15.5, family=FONT_DISPLAY, weight='bold', color=COLORS['ink'], y=1.01)
    fig.tight_layout()
    out_path = out_dir / 'cross_run_trends.png'
    fig.savefig(out_path, dpi=170, bbox_inches='tight')
    plt.close(fig)
    print(f'saved cross-run trend plot -> {out_path}')

def _process_file_worker(args_tuple):
    fp, out_dir, tree_name = args_tuple
    try:
        return ('ok', fp.name, process_file(fp, out_dir, tree_name=tree_name))
    except Exception as e:
        return ('error', fp.name, str(e))

def _run_batch(files, out_dir, tree_name, n_jobs):
    rows = []
    total = len(files)
    progress = tqdm(total=total, unit='file', desc='Processing runs') if _HAVE_TQDM else None

    def _report_plain(done, name, status):
        marker = 'ok' if status == 'ok' else 'FAILED'
        sys.stdout.write(f'\r[{done:>4}/{total}] {marker:<7} {name[:60]:<60}')
        sys.stdout.flush()
    if n_jobs <= 1:
        for i, fp in enumerate(files, start=1):
            status, name, result = _process_file_worker((fp, out_dir, tree_name))
            if status == 'ok':
                rows.append(result)
            else:
                (progress.write if progress else print)(f'  [error] failed on {name}: {result}')
            if progress:
                progress.update(1)
            else:
                _report_plain(i, name, status)
    else:
        work = [(fp, out_dir, tree_name) for fp in files]
        done = 0
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {pool.submit(_process_file_worker, w): w[0] for w in work}
            for future in as_completed(futures):
                status, name, result = future.result()
                done += 1
                if status == 'ok':
                    rows.append(result)
                else:
                    (progress.write if progress else print)(f'  [error] failed on {name}: {result}')
                if progress:
                    progress.update(1)
                else:
                    _report_plain(done, name, status)
    if progress:
        progress.close()
    elif not _HAVE_TQDM:
        sys.stdout.write('\n')
        sys.stdout.flush()
    return rows

def find_root_files(input_dir: str) -> list:
    return sorted(Path(input_dir).glob('*.root'))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-dir', type=str, default=None, help='Directory containing .root files to process (all *.root files are used).')
    ap.add_argument('--files', nargs='+', default=None, help='Explicit list of .root file paths (alternative to --input-dir).')
    ap.add_argument('--out-dir', type=str, default='scatter_fit_results', help='Directory to write plots and summary CSV into.')
    ap.add_argument('--tree', type=str, default='PionEvents', help='Tree name inside the ROOT files (default: PionEvents).')
    ap.add_argument('--jobs', '-j', type=int, default=1, help="Number of files to process in parallel (default: 1, sequential). Each file's read + fit + plot is CPU-bound and independent, so this scales well up to roughly your core count. Use -j0 to auto-detect.")
    args = ap.parse_args()
    if not args.input_dir and (not args.files):
        ap.error('Provide either --input-dir or --files')
    files = find_root_files(args.input_dir) if args.input_dir else [Path(f) for f in args.files]
    if not files:
        print('No .root files found.', file=sys.stderr)
        sys.exit(1)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n_jobs = args.jobs
    if n_jobs == 0:
        import os
        n_jobs = os.cpu_count() or 1
    n_jobs = max(1, min(n_jobs, len(files)))
    if not _HAVE_TQDM:
        print('(tip: `pip install tqdm` for a nicer progress bar -- using a plain fallback for now)')
    print(f"Processing {len(files)} file(s) with --jobs {n_jobs}{(' (parallel)' if n_jobs > 1 else ' (sequential)')}\n")
    rows = _run_batch(files, out_dir, args.tree, n_jobs)
    if not rows:
        print('No files processed successfully.', file=sys.stderr)
        sys.exit(1)
    csv_path = out_dir / 'summary.csv'
    fieldnames = [k for k in rows[0].keys() if not k.startswith('_')]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(f'\nsaved summary CSV -> {csv_path}')
    make_trend_plots(rows, out_dir)
    make_overlay_plots(rows, out_dir)
    print(f'\nDone. Processed {len(rows)}/{len(files)} files. Results in: {out_dir}/')
if __name__ == '__main__':
    main()