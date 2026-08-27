import argparse, os, re, sys, glob, pickle, warnings
import numpy as np
from scipy import optimize, integrate
warnings.filterwarnings('ignore', category=RuntimeWarning)
FNAME_PATTERN = re.compile('^(?:(?P<tag_prefix>exp\\d+)_)?(?P<structure>.+?)_fold(?P<fold>[^_]+)_p(?P<momentum>[^_]+)GeV_vtx(?P<vertex>[^_]+)mm(?:_tag(?P<tag>[^_]+))?_run(?P<run>\\d+)(?:_(?P<idx>\\d+))?\\.root$')

def decode_pfloat(s):
    if s is None:
        return None
    try:
        return float(s.replace('p', '.'))
    except ValueError:
        return None

def parse_filename(path):
    base = os.path.basename(path)
    m = FNAME_PATTERN.match(base)
    if not m:
        return dict(label=base.replace('.root', ''), momentum_gev=None, fold=None, vertex_mm=None, tag=None, structure=None, run=None, idx=None, parsed=False)
    g = m.groupdict()
    mom = decode_pfloat(g.get('momentum'))
    tag = g.get('tag') or g.get('tag_prefix') or 'run'
    idx_part = f"_{g['idx']}" if g.get('idx') is not None else ''
    run_part = f"_run{g['run']}" if g.get('run') is not None else ''
    structure = g.get('structure')
    label = f'{tag}_{structure}_p{mom:g}GeV{run_part}{idx_part}' if mom is not None else base.replace('.root', '')
    return dict(label=label, momentum_gev=mom, fold=decode_pfloat(g.get('fold')), vertex_mm=decode_pfloat(g.get('vertex')), tag=tag, structure=structure, run=g.get('run'), idx=g.get('idx'), parsed=True)

def rayleigh_pdf(theta, theta0):
    return theta / theta0 ** 2 * np.exp(-theta ** 2 / (2 * theta0 ** 2))

def rayleigh_cdf(theta, theta0):
    return 1 - np.exp(-theta ** 2 / (2 * theta0 ** 2))

def smooth_mixture_unnorm(theta, theta0, theta_min, n, s):
    logistic = 1.0 / (1.0 + np.exp(-(np.log(theta) - np.log(theta_min)) / s))
    core = rayleigh_pdf(theta, theta0)
    Rmin = rayleigh_pdf(theta_min, theta0)
    tail = Rmin * (theta / theta_min) ** (-n)
    return (1 - logistic) * core + logistic * tail

def smooth_mixture_norm(theta0, theta_min, n, s):
    part1, _ = integrate.quad(lambda t: smooth_mixture_unnorm(t, theta0, theta_min, n, s), 1e-09, theta_min * 50, limit=200)
    Rmin = rayleigh_pdf(theta_min, theta0)
    tail_start = theta_min * 50
    part2 = Rmin * theta_min ** n * tail_start ** (1 - n) / (n - 1)
    return part1 + part2

def smooth_mixture_pdf(theta, theta0, theta_min, n, s):
    norm = smooth_mixture_norm(theta0, theta_min, n, s)
    return smooth_mixture_unnorm(theta, theta0, theta_min, n, s) / norm

def binned_nll_smooth(params, centers, widths, counts):
    theta0, theta_min, n, s = params
    if theta0 <= 0 or theta_min <= 0 or n <= 1.01 or (s <= 0.02) or (s > 3) or (theta_min < 0.5 * theta0):
        return 10000000000.0
    norm = smooth_mixture_norm(theta0, theta_min, n, s)
    if norm <= 0 or not np.isfinite(norm):
        return 10000000000.0
    pdf_vals = np.clip(smooth_mixture_unnorm(centers, theta0, theta_min, n, s) / norm, 1e-300, None)
    return -np.sum(counts * np.log(pdf_vals * widths))

def nll_rayleigh_only(theta0, data):
    if theta0 <= 0:
        return 10000000000.0
    pdf = np.clip(rayleigh_pdf(data, theta0), 1e-300, None)
    return -np.sum(np.log(pdf))

def binned_chi2(data, pdf_func, n_bins=60):
    edges = np.logspace(np.log10(np.percentile(data, 0.02)), np.log10(data.max() * 1.001), n_bins)
    obs, _ = np.histogram(data, bins=edges)
    N = len(data)
    exp = np.array([N * integrate.quad(pdf_func, edges[i], edges[i + 1], limit=200)[0] for i in range(len(edges) - 1)])
    mask = exp > 5
    chi2 = np.sum((obs[mask] - exp[mask]) ** 2 / exp[mask])
    return (chi2, mask.sum())

def load_scatter_angles(root_path):
    import uproot
    arr = uproot.open(root_path)['PionEvents'].arrays(['scatterAngleDeg', 'nHitsPrimary'], library='np')
    mask = arr['nHitsPrimary'] == 1
    sa = arr['scatterAngleDeg'][mask]
    return sa[sa > 0]

def fit_one_run(data, theta0_seed=None, theta_min_seed=None, n_bins=300):
    if theta0_seed is None:
        theta0_seed = np.median(data) / 1.1774
    if theta_min_seed is None:
        theta_min_seed = np.percentile(data, 99.0)
    edges = np.logspace(np.log10(data.min()), np.log10(data.max() * 1.001), n_bins)
    counts, edges = np.histogram(data, bins=edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    keep = counts > 0
    res = optimize.minimize(binned_nll_smooth, [theta0_seed, theta_min_seed, 3.5, 0.25], args=(centers[keep], widths[keep], counts[keep]), method='Nelder-Mead', options={'maxiter': 30000, 'maxfev': 30000, 'xatol': 1e-09, 'fatol': 1e-07, 'adaptive': True})
    theta0, theta_min, n, s = res.x
    res_ray = optimize.minimize_scalar(nll_rayleigh_only, bracket=(theta0_seed * 0.5, theta0_seed, theta0_seed * 2), args=(data,))
    theta0_ray = res_ray.x
    chi2_mix, nb_mix = binned_chi2(data, lambda x: smooth_mixture_pdf(x, theta0, theta_min, n, s))
    chi2_ray, nb_ray = binned_chi2(data, lambda x: rayleigh_pdf(x, theta0_ray))
    return dict(theta0=theta0, theta_min=theta_min, n=n, s=s, theta0_rayleigh_only=theta0_ray, chi2_mixture=chi2_mix, ndof_mixture=nb_mix - 4, chi2_rayleigh=chi2_ray, ndof_rayleigh=nb_ray - 1, n_events=len(data), fit_success=bool(res.success))

def rayleigh_only_pdf_from_fit(theta0):
    return lambda x: rayleigh_pdf(x, theta0)

def make_per_run_plot(label, data, fit, out_path, style='paper'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    for f in ['/usr/share/texmf/fonts/opentype/public/lm/lmroman10-regular.otf', '/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf', '/usr/share/texmf/fonts/opentype/public/lm/lmroman10-italic.otf']:
        try:
            fm.fontManager.addfont(f)
        except Exception:
            pass
    SERIF = ['Latin Modern Roman', 'DejaVu Serif', 'Times New Roman', 'serif']
    if style == 'paper':
        plt.rcParams.update({'font.family': 'serif', 'font.serif': SERIF, 'font.size': 9, 'mathtext.fontset': 'cm', 'text.color': 'black', 'axes.edgecolor': 'black', 'axes.labelcolor': 'black', 'xtick.color': 'black', 'ytick.color': 'black', 'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True, 'xtick.minor.visible': True, 'ytick.minor.visible': True, 'axes.linewidth': 0.8, 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'legend.frameon': False, 'legend.fontsize': 8})
        C_FILL, C_EDGE = ('#CDE0F5', '#4C7EB3')
        C_MODEL, C_BASE, C_PTS, C_VLINE = ('#D55E00', '#666666', '#0B3D6E', '#999999')
    else:
        plt.rcParams.update({'font.family': 'serif', 'font.serif': SERIF, 'mathtext.fontset': 'cm', 'text.color': '#1C2B47', 'axes.edgecolor': '#1C2B47', 'axes.labelcolor': '#1C2B47', 'xtick.color': '#4A5A78', 'ytick.color': '#4A5A78', 'axes.linewidth': 1.1, 'figure.facecolor': '#FAF6EE', 'axes.facecolor': '#FFFFFF', 'savefig.facecolor': '#FAF6EE'})
        C_FILL, C_EDGE = ('#B9C7DE', 'none')
        C_MODEL, C_BASE, C_PTS, C_VLINE = ('#B5502D', '#7C8B5A', '#25406B', '#4A5A78')
    theta0, theta_min, n, s = (fit['theta0'], fit['theta_min'], fit['n'], fit['s'])
    theta0_ray = fit['theta0_rayleigh_only']
    xmax = np.percentile(data, 99.99) * 1.5
    xs = np.linspace(1e-06, xmax, 3000)
    model_vals = smooth_mixture_pdf(xs, theta0, theta_min, n, s)
    fig, (ax_top, ax_bot) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    ax_top.hist(data, bins=120, range=(0, np.percentile(data, 99.5)), density=True, histtype='stepfilled', facecolor=C_FILL, edgecolor=C_EDGE, linewidth=0.6, zorder=2, label='Simulation')
    ax_top.plot(xs, rayleigh_pdf(xs, theta0_ray), color=C_BASE, ls='--', lw=1.4, zorder=3, label='Gaussian (Highland)')
    ax_top.plot(xs, model_vals, color=C_MODEL, lw=1.6, zorder=4, label='Core + tail fit')
    ax_top.set_xlim(0, np.percentile(data, 99.5))
    ax_top.set_ylim(bottom=0)
    ax_top.set_xlabel('$\\theta$ (deg)')
    ax_top.set_ylabel('$f(\\theta)$')
    ax_top.legend(fontsize=7.5, loc='upper right')
    n_entries = len(data)
    x_lo, x_hi = (data.min() * 0.7, data.max() * 1.3)
    bins_log = np.logspace(np.log10(np.percentile(data, 0.01)), np.log10(data.max() * 1.05), 40)
    counts, edges = np.histogram(data, bins=bins_log)
    centers = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    dens = counts / (n_entries * widths)
    err = np.sqrt(np.clip(counts, 1, None)) / (n_entries * widths)
    mask = counts >= 3
    xs_full = np.logspace(np.log10(x_lo), np.log10(x_hi), 3000)
    ax_bot.plot(xs_full, rayleigh_pdf(xs_full, theta0_ray), color=C_BASE, ls='--', lw=1.3, zorder=3)
    ax_bot.plot(xs_full, smooth_mixture_pdf(xs_full, theta0, theta_min, n, s), color=C_MODEL, lw=1.6, zorder=4)
    ax_bot.errorbar(centers[mask], dens[mask], yerr=err[mask], fmt='o', ms=2.8, mfc=C_PTS, mec=C_PTS, ecolor=C_PTS, elinewidth=0.7, capsize=1.5, zorder=5)
    ax_bot.axvline(theta_min, color=C_VLINE, ls=':', lw=1.0, zorder=1)
    ax_bot.set_xscale('log')
    ax_bot.set_yscale('log')
    ax_bot.set_xlim(x_lo, x_hi)
    ax_bot.set_ylim(bottom=dens[mask].min() * 0.3, top=dens[mask].max() * 3)
    ax_bot.set_xlabel('$\\theta$ (deg)')
    ax_bot.set_ylabel('$f(\\theta)$')
    fig.suptitle(label.replace('_', ' '), fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)

def make_grid_overview(fit_cache, out_path, style='paper', ncols=4):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    for f in ['/usr/share/texmf/fonts/opentype/public/lm/lmroman10-regular.otf', '/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf']:
        try:
            fm.fontManager.addfont(f)
        except Exception:
            pass
    SERIF = ['Latin Modern Roman', 'DejaVu Serif', 'Times New Roman', 'serif']
    if style == 'paper':
        plt.rcParams.update({'font.family': 'serif', 'font.serif': SERIF, 'font.size': 8, 'mathtext.fontset': 'cm', 'text.color': 'black', 'axes.edgecolor': 'black', 'axes.labelcolor': 'black', 'xtick.color': 'black', 'ytick.color': 'black', 'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 0.7, 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white'})
        C_MODEL, C_BASE, C_PTS, C_VLINE = ('#D55E00', '#666666', '#0B3D6E', '#999999')
        bg = 'white'
    else:
        plt.rcParams.update({'font.family': 'serif', 'font.serif': SERIF, 'mathtext.fontset': 'cm', 'text.color': '#1C2B47', 'axes.edgecolor': '#1C2B47', 'axes.labelcolor': '#1C2B47', 'xtick.color': '#4A5A78', 'ytick.color': '#4A5A78', 'axes.linewidth': 0.9, 'figure.facecolor': '#FAF6EE', 'axes.facecolor': '#FFFFFF', 'savefig.facecolor': '#FAF6EE'})
        C_MODEL, C_BASE, C_PTS, C_VLINE = ('#B5502D', '#7C8B5A', '#25406B', '#4A5A78')
        bg = '#FAF6EE'
    items = sorted(fit_cache.items(), key=lambda kv: (kv[1]['meta']['momentum_gev'] is None, kv[1]['meta']['momentum_gev'] or 0, kv[0]))
    n = len(items)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.6 * ncols, 2.3 * nrows), squeeze=False)
    for idx, (label, entry) in enumerate(items):
        ax = axes[idx // ncols, idx % ncols]
        data, fit = (entry['data'], entry['fit'])
        theta0, theta_min, n_tail, s = (fit['theta0'], fit['theta_min'], fit['n'], fit['s'])
        theta0_ray = fit['theta0_rayleigh_only']
        n_entries = len(data)
        x_lo, x_hi = (data.min() * 0.7, data.max() * 1.3)
        bins_log = np.logspace(np.log10(np.percentile(data, 0.01)), np.log10(data.max() * 1.05), 30)
        counts, edges = np.histogram(data, bins=bins_log)
        centers = np.sqrt(edges[:-1] * edges[1:])
        widths = np.diff(edges)
        dens = counts / (n_entries * widths)
        mask = counts >= 3
        xs_full = np.logspace(np.log10(x_lo), np.log10(x_hi), 500)
        ax.plot(xs_full, rayleigh_pdf(xs_full, theta0_ray), color=C_BASE, ls='--', lw=0.9, zorder=3)
        ax.plot(xs_full, smooth_mixture_pdf(xs_full, theta0, theta_min, n_tail, s), color=C_MODEL, lw=1.1, zorder=4)
        ax.plot(centers[mask], dens[mask], 'o', ms=1.8, color=C_PTS, zorder=5)
        ax.axvline(theta_min, color=C_VLINE, ls=':', lw=0.7, zorder=1)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(bottom=dens[mask].min() * 0.3, top=dens[mask].max() * 3)
        mom = entry['meta']['momentum_gev']
        title = f'p={mom:g} GeV/c' if mom is not None else label
        ax.set_title(title, fontsize=7.5, pad=2)
        ax.tick_params(labelsize=5.5)
    for j in range(n, nrows * ncols):
        axes[j // ncols, j % ncols].axis('off')
    fig.suptitle('Core + Rutherford-tail fit — all runs (log-log, tail region)', fontsize=10, y=1.0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder', help='Folder containing .root files (searched recursively)')
    ap.add_argument('--out', default='./mcs_batch_output', help='Output directory')
    ap.add_argument('--pattern', default='*.root', help='Glob pattern for input files')
    ap.add_argument('--style', choices=['paper', 'pretty'], default='paper', help='Plot style for per-run figures')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(os.path.join(args.out, 'per_run_plots'), exist_ok=True)
    if any((c in args.folder for c in '*?[]')) or os.path.isfile(args.folder):
        files = sorted(glob.glob(args.folder))
    elif os.path.isdir(args.folder):
        files = sorted(glob.glob(os.path.join(args.folder, '**', args.pattern), recursive=True))
    else:
        files = sorted(glob.glob(os.path.join(args.folder, '**', args.pattern), recursive=True))
        if not files:
            files = sorted(glob.glob(args.folder))
    if not files:
        print(f'No files matching {args.pattern} found under {args.folder}', file=sys.stderr)
        sys.exit(1)
    print(f'Found {len(files)} file(s).')
    rows = []
    fit_cache = {}
    for i, path in enumerate(files):
        meta = parse_filename(path)
        label = meta['label']
        print(f'[{i + 1}/{len(files)}] {os.path.basename(path)}  ->  label={label}')
        try:
            data = load_scatter_angles(path)
            if len(data) < 200:
                print(f'    SKIP: only {len(data)} valid events after selection (need >=200)')
                continue
            fit = fit_one_run(data)
        except Exception as e:
            print(f'    FAILED: {e}')
            continue
        row = {**meta, **fit, 'path': path}
        rows.append(row)
        fit_cache[label] = dict(fit=fit, meta=meta, data=data)
        try:
            plot_path = os.path.join(args.out, 'per_run_plots', f'{label}.png')
            make_per_run_plot(label, data, fit, plot_path, style=args.style)
        except Exception as e:
            print(f'    (plot failed: {e})')
        chi2_ratio = fit['chi2_rayleigh'] / max(fit['ndof_rayleigh'], 1) / max(fit['chi2_mixture'] / max(fit['ndof_mixture'], 1), 1e-09)
        print(f"    theta0={fit['theta0']:.5f} deg   n_tail={fit['n']:.2f}   chi2/ndof: Rayleigh-only={fit['chi2_rayleigh'] / max(fit['ndof_rayleigh'], 1):.1f}  mixture={fit['chi2_mixture'] / max(fit['ndof_mixture'], 1):.2f}  (improvement {chi2_ratio:.0f}x)")
    if not rows:
        print('No successful fits -- nothing to write.', file=sys.stderr)
        sys.exit(1)
    import csv
    csv_path = os.path.join(args.out, 'summary.csv')
    fieldnames = ['label', 'momentum_gev', 'fold', 'vertex_mm', 'tag', 'structure', 'run', 'idx', 'n_events', 'theta0', 'theta_min', 'n', 's', 'theta0_rayleigh_only', 'chi2_mixture', 'ndof_mixture', 'chi2_rayleigh', 'ndof_rayleigh', 'fit_success', 'path']
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader()
        for r in sorted(rows, key=lambda r: r['momentum_gev'] if r['momentum_gev'] is not None else 1000000000.0):
            w.writerow(r)
    print(f'\nWrote summary table: {csv_path}')
    with open(os.path.join(args.out, 'fit_cache.pkl'), 'wb') as f:
        pickle.dump(fit_cache, f)
    print(f"Wrote fit cache: {os.path.join(args.out, 'fit_cache.pkl')}")
    moms = [(r['momentum_gev'], r['theta0']) for r in rows if r['momentum_gev']]
    if len(moms) >= 3:
        moms.sort()
        p_arr = np.array([m[0] for m in moms])
        th0_arr = np.array([m[1] for m in moms])
        logp, logth0 = (np.log(p_arr), np.log(th0_arr))
        k, logA = np.polyfit(logp, logth0, 1)
        print(f'\nMomentum scaling across {len(moms)} runs: theta0 ~ p^({k:.3f})  (Highland predicts ~ -1, i.e. 1/p)')
    print(f"\nDone. Per-run plots go in {os.path.join(args.out, 'per_run_plots')} (run with --style paper|pretty to choose the look).")
    if len(fit_cache) >= 2:
        grid_path = os.path.join(args.out, 'grid_overview.png')
        try:
            make_grid_overview(fit_cache, grid_path, style=args.style)
            print(f'Wrote grid overview: {grid_path}')
        except Exception as e:
            print(f'Grid overview failed: {e}')
if __name__ == '__main__':
    main()