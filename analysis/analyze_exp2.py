import sys
import glob
import os
import re
import json
import numpy as np
import uproot
import matplotlib.pyplot as plt
try:
    from scipy.optimize import curve_fit, minimize
    from scipy.integrate import quad
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
SILICON_X0_MM = 93.66
KAPTON_X0_MM = 285.75
PION_MASS_GEV = 0.13957
THETA_MIN_FACTOR = 3.0

def highland_theta0_rad(p_GeV, x_over_x0, mass_GeV=PION_MASS_GEV, z=1):
    p_GeV = np.asarray(p_GeV, dtype=float)
    x_over_x0 = np.asarray(x_over_x0, dtype=float)
    E_GeV = np.sqrt(p_GeV ** 2 + mass_GeV ** 2)
    beta = p_GeV / E_GeV
    theta0 = 0.0136 / (beta * p_GeV) * z * np.sqrt(x_over_x0) * (1.0 + 0.038 * np.log(x_over_x0))
    return theta0

def theta0_from_space_angle(space_angle_rad):
    return np.asarray(space_angle_rad, dtype=float) / np.sqrt(2.0)

def _synthetic_rayleigh_selftest(n=200000, true_theta0_rad=0.0005, seed=0):
    rng = np.random.default_rng(seed)
    theta_x = rng.normal(0.0, true_theta0_rad, n)
    theta_y = rng.normal(0.0, true_theta0_rad, n)
    space = np.sqrt(theta_x ** 2 + theta_y ** 2)
    result = fit_theta0(space, core_fraction=0.9)
    ok = True
    print(f'Self-test: injected theta_0 = {true_theta0_rad * 1000:.4f} mrad, n={n}')
    if result['rms_rad'] is None:
        print('  FAIL: rms_rad is None')
        ok = False
    else:
        rel_err = abs(result['rms_rad'] - true_theta0_rad) / true_theta0_rad
        status = 'PASS' if rel_err < 0.02 else 'FAIL'
        if status == 'FAIL':
            ok = False
        print(f"  RMS estimate:      {result['rms_mrad']:.4f} mrad (rel. err {rel_err * 100:.2f}%)  [{status}]")
    if result['gaussian_fit_rad'] is None:
        print('  NOTE: Rayleigh-fit estimate unavailable (scipy missing or fit did not converge) - only RMS estimate was checked.')
    else:
        rel_err = abs(result['gaussian_fit_rad'] - true_theta0_rad) / true_theta0_rad
        status = 'PASS' if rel_err < 0.02 else 'FAIL'
        if status == 'FAIL':
            ok = False
        print(f"  Rayleigh-fit est.: {result['gaussian_fit_mrad']:.4f} mrad (rel. err {rel_err * 100:.2f}%)  [{status}]")
    print(f"  Self-test overall: {('PASS' if ok else 'FAIL')}")
    return ok

def rayleigh_pdf(theta, sigma0):
    theta = np.asarray(theta, dtype=float)
    sigma0 = float(sigma0)
    return theta / sigma0 ** 2 * np.exp(-0.5 * (theta / sigma0) ** 2)

def powerlaw_tail_pdf(theta, theta_min, alpha):
    theta = np.asarray(theta, dtype=float)
    theta_min = float(theta_min)
    alpha = float(alpha)
    pdf = np.zeros_like(theta)
    mask = theta >= theta_min
    if np.any(mask):
        pdf[mask] = (alpha - 1.0) / theta_min * (theta[mask] / theta_min) ** (-alpha)
    return pdf

def mixture_pdf(theta, sigma0, f, alpha, theta_min_factor=THETA_MIN_FACTOR):
    theta = np.asarray(theta, dtype=float)
    theta_min = theta_min_factor * float(sigma0)
    core = rayleigh_pdf(theta, sigma0)
    tail = powerlaw_tail_pdf(theta, theta_min, alpha)
    return (1.0 - f) * core + f * tail

def _mixture_neg_log_likelihood(params, theta, theta_min_factor):
    sigma0, f, alpha = params
    if sigma0 <= 0.0 or not 0.0 <= f <= 1.0 or alpha <= 1.001:
        return np.inf
    pdf = mixture_pdf(theta, sigma0, f, alpha, theta_min_factor)
    if np.any(pdf <= 0.0) or np.any(~np.isfinite(pdf)):
        return np.inf
    return -np.sum(np.log(pdf))

def _fit_mixture_point_estimate(theta, x0, bounds, theta_min_factor):
    res = minimize(_mixture_neg_log_likelihood, x0, args=(theta, theta_min_factor), method='Nelder-Mead', bounds=bounds, options={'xatol': 1e-09, 'fatol': 1e-07, 'maxiter': 4000, 'maxfev': 4000})
    return res

def fit_gaussian_core_powerlaw_tail(space_angle_rad, theta_min_factor=THETA_MIN_FACTOR, n_bootstrap=200, seed=0, verbose=True):
    if not HAVE_SCIPY:
        if verbose:
            print('    NOTE: scipy not available - skipping mixture core+tail fit entirely (need scipy.optimize.minimize).')
        return None
    theta = np.asarray(space_angle_rad, dtype=float)
    theta = theta[np.isfinite(theta) & (theta > 0.0)]
    n = len(theta)
    if n < 200:
        if verbose:
            print(f'    NOTE: only {n} valid events - too few for a stable 3-parameter core+tail fit (need >=200). Skipping.')
        return None
    core_cut = np.percentile(theta, 80.0)
    core_events = theta[theta <= core_cut]
    if len(core_events) >= 10:
        sigma0_guess = float(np.sqrt(np.mean(core_events ** 2) / 2.0))
    else:
        sigma0_guess = float(np.sqrt(np.mean(theta ** 2) / 2.0))
    sigma0_guess = max(sigma0_guess, 1e-09)
    theta_min_guess = theta_min_factor * sigma0_guess
    f_guess = float(np.mean(theta > theta_min_guess))
    f_guess = min(max(f_guess, 0.001), 0.4)
    alpha_guess = 3.0
    x0 = np.array([sigma0_guess, f_guess, alpha_guess])
    bounds = [(sigma0_guess * 0.05, sigma0_guess * 20.0), (0.0001, 0.5), (1.05, 8.0)]
    res = _fit_mixture_point_estimate(theta, x0, bounds, theta_min_factor)
    if not res.success and verbose:
        print(f'    NOTE: mixture MLE fit reported non-convergence ({res.message}) - using best point found; treat uncertainties/chi2 below with extra caution for this run.')
    sigma0_fit, f_fit, alpha_fit = [float(v) for v in res.x]
    rng = np.random.default_rng(seed)
    boot_params = np.full((n_bootstrap, 3), np.nan)
    for i in range(n_bootstrap):
        sample = rng.choice(theta, size=n, replace=True)
        r = _fit_mixture_point_estimate(sample, res.x, bounds, theta_min_factor)
        if np.isfinite(r.fun):
            boot_params[i] = r.x
    valid_boot = boot_params[np.all(np.isfinite(boot_params), axis=1)]
    if len(valid_boot) >= 20:
        sigma0_err = float(np.std(valid_boot[:, 0], ddof=1))
        f_err = float(np.std(valid_boot[:, 1], ddof=1))
        alpha_err = float(np.std(valid_boot[:, 2], ddof=1))
    else:
        sigma0_err = f_err = alpha_err = None
        if verbose:
            print(f'    NOTE: only {len(valid_boot)}/{n_bootstrap} bootstrap replicates converged - parameter uncertainties and the plotted confidence band are unavailable for this run.')
    theta_min_fit = theta_min_factor * sigma0_fit
    chi2, ndof = (None, None)
    try:
        lo = max(float(np.min(theta)), 1e-09)
        hi = float(np.max(theta))
        edges = np.geomspace(lo, hi, 26)
        obs, _ = np.histogram(theta, bins=edges)
        exp = np.empty(len(obs))
        for bi, (e_lo, e_hi) in enumerate(zip(edges[:-1], edges[1:])):
            val, _ = quad(lambda t: mixture_pdf(np.array([t]), sigma0_fit, f_fit, alpha_fit, theta_min_factor)[0], e_lo, e_hi, limit=100)
            exp[bi] = val * n
        use = exp >= 5.0
        if use.sum() > 3:
            chi2 = float(np.sum((obs[use] - exp[use]) ** 2 / exp[use]))
            ndof = int(use.sum() - 3)
    except Exception as e:
        if verbose:
            print(f'    NOTE: binned chi2 diagnostic failed ({e}) - point estimates and errors above are unaffected.')
    return {'n_events': n, 'sigma0_rad': sigma0_fit, 'sigma0_mrad': sigma0_fit * 1000.0, 'sigma0_err_rad': sigma0_err, 'sigma0_err_mrad': sigma0_err * 1000.0 if sigma0_err is not None else None, 'f_tail': f_fit, 'f_tail_err': f_err, 'alpha': alpha_fit, 'alpha_err': alpha_err, 'theta_min_factor': theta_min_factor, 'theta_min_rad': theta_min_fit, 'theta_min_mrad': theta_min_fit * 1000.0, 'n_bootstrap_requested': n_bootstrap, 'n_bootstrap_successful': int(len(valid_boot)), 'chi2': chi2, 'ndof': ndof, 'chi2_per_ndof': chi2 / ndof if chi2 is not None and ndof else None, 'converged': bool(res.success), '_bootstrap_params': valid_boot}

def _mixture_curve_band(theta_grid_rad, boot_params, theta_min_factor=THETA_MIN_FACTOR, percentiles=(16.0, 84.0)):
    curves = np.array([mixture_pdf(theta_grid_rad, s, f, a, theta_min_factor) for s, f, a in boot_params])
    lo = np.percentile(curves, percentiles[0], axis=0)
    hi = np.percentile(curves, percentiles[1], axis=0)
    return (lo, hi)

def _synthetic_mixture_selftest(n=300000, true_sigma0_rad=0.0004, true_f=0.05, true_alpha=3.0, seed=0):
    rng = np.random.default_rng(seed)
    n_tail = rng.binomial(n, true_f)
    n_core = n - n_tail
    u = rng.uniform(1e-12, 1.0, n_core)
    core_sample = true_sigma0_rad * np.sqrt(-2.0 * np.log(u))
    theta_min = THETA_MIN_FACTOR * true_sigma0_rad
    u2 = rng.uniform(1e-12, 1.0, n_tail)
    tail_sample = theta_min * u2 ** (-1.0 / (true_alpha - 1.0))
    sample = np.concatenate([core_sample, tail_sample])
    rng.shuffle(sample)
    print(f'Self-test: injected sigma0={true_sigma0_rad * 1000:.4f} mrad, f={true_f:.3f}, alpha={true_alpha:.2f}, n={n} (~{n_tail} tail events)')
    fit = fit_gaussian_core_powerlaw_tail(sample, n_bootstrap=0, verbose=False)
    if fit is None:
        print('  FAIL: fit returned None')
        return False
    ok = True
    checks = [('sigma0', fit['sigma0_rad'], true_sigma0_rad, 0.05), ('f_tail', fit['f_tail'], true_f, 0.25), ('alpha', fit['alpha'], true_alpha, 0.25)]
    for label, got, truth, tol in checks:
        rel_err = abs(got - truth) / truth
        status = 'PASS' if rel_err < tol else 'FAIL'
        if status == 'FAIL':
            ok = False
        print(f'  {label:8s}: fit={got:.4f}  truth={truth:.4f}  rel.err={rel_err * 100:.1f}%  (tol {tol * 100:.0f}%)  [{status}]')
    print(f"  Self-test overall: {('PASS' if ok else 'FAIL')}")
    return ok

def _rayleigh_pdf_counts(x, sigma, amplitude, bin_width):
    x = np.asarray(x, dtype=float)
    return amplitude * bin_width * (x / sigma ** 2) * np.exp(-0.5 * (x / sigma) ** 2)

def fit_theta0(space_angle_rad, core_fraction=0.9):
    space = np.asarray(space_angle_rad, dtype=float)
    space = space[np.isfinite(space) & (space >= 0.0)]
    n = len(space)
    if n < 10:
        return {'n_events': n, 'rms_rad': None, 'gaussian_fit_rad': None, 'gaussian_fit_error_rad': None}
    rms_rad = float(theta0_from_space_angle(np.sqrt(np.mean(space ** 2))))
    gaussian_fit_rad = None
    gaussian_fit_error_rad = None
    fit_failed_reason = None
    if HAVE_SCIPY and n >= 50:
        hi_pct = core_fraction * 100.0
        hi_cut = np.percentile(space, hi_pct)
        core = space[space <= hi_cut]
        n_bins = max(20, min(80, int(np.sqrt(len(core)))))
        counts, edges = np.histogram(core, bins=n_bins)
        centers = 0.5 * (edges[:-1] + edges[1:])
        bin_width = edges[1] - edges[0]
        try:
            sigma0 = rms_rad if rms_rad else np.mean(core)
            p0 = [sigma0, len(core)]
            fit_func = lambda x, sigma, amplitude: _rayleigh_pdf_counts(x, sigma, amplitude, bin_width)
            popt, pcov = curve_fit(fit_func, centers, counts, p0=p0, bounds=([1e-12, 1.0], [np.inf, np.inf]), maxfev=5000)
            fitted_sigma = float(popt[0])
            if fitted_sigma <= 0.0 or not np.isfinite(fitted_sigma):
                fit_failed_reason = f'non-physical fitted sigma ({fitted_sigma})'
            else:
                gaussian_fit_rad = fitted_sigma
                if np.isfinite(pcov[0, 0]) and pcov[0, 0] >= 0:
                    gaussian_fit_error_rad = float(np.sqrt(pcov[0, 0]))
        except Exception as e:
            fit_failed_reason = str(e)
        if fit_failed_reason is not None:
            print(f'    NOTE: Rayleigh core fit did not converge to a physical result ({fit_failed_reason}) - falling back to RMS only for this point.')
    return {'n_events': n, 'rms_rad': rms_rad, 'rms_mrad': rms_rad * 1000.0, 'gaussian_fit_rad': gaussian_fit_rad, 'gaussian_fit_mrad': gaussian_fit_rad * 1000.0 if gaussian_fit_rad is not None else None, 'gaussian_fit_error_rad': gaussian_fit_error_rad}

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

def validate_momentum(path, data, tol_frac=0.02):
    filename_stem = os.path.basename(path).replace('.root', '')
    claimed_GeV = _parse_filename_tag(filename_stem, '_p', 'GeV')
    if claimed_GeV is None:
        return None
    if not 0.1 <= claimed_GeV <= 20.0:
        print(f"  WARNING: filename-claimed momentum {claimed_GeV} GeV/c in '{filename_stem}' is outside the expected 0.1-20 GeV/c range for this sweep - check for a hand-edited macro or unit mixup.")
        return {'claimed_GeV': claimed_GeV, 'ok': False}
    return {'claimed_GeV': claimed_GeV, 'ok': True}

def summarize_scattering(name, momentum_GeV, path, data, make_tail_plots=True, tail_plot_dir='tail_fit_plots', n_bootstrap=200):
    scatter_deg = data.get('scatterAngleDeg')
    if scatter_deg is None:
        print(f'  WARNING: no scatterAngleDeg column in this file - was it produced by a Mark 3 build? Skipping scattering summary.')
        return None
    hit = data['hitDetector'].astype(bool)
    has_split = 'siliconPathLength_mm' in data and 'kaptonPathLength_mm' in data
    if has_split:
        si_mm = data['siliconPathLength_mm']
        kap_mm = data['kaptonPathLength_mm']
        path_mm = si_mm + kap_mm
        x_over_x0_all = si_mm / SILICON_X0_MM + kap_mm / KAPTON_X0_MM
    else:
        path_mm = data['totalPathLength_mm']
        x_over_x0_all = path_mm / SILICON_X0_MM
    valid = hit & (scatter_deg >= 0.0) & (path_mm > 0.0)
    n_valid = int(valid.sum())
    if n_valid < 10:
        print(f'  WARNING: only {n_valid} valid scattering events (hitDetector=1, scatterAngleDeg>=0) - too few to fit theta_0 reliably at p={momentum_GeV} GeV/c.')
    n_hits_diag = data.get('nHitsPrimary')
    if n_hits_diag is not None:
        n_hits_diag = np.asarray(n_hits_diag)
        m_ref = valid & (n_hits_diag == 1)
        path_ref = float(np.mean(path_mm[m_ref])) if m_ref.sum() >= 30 else None
        for k in sorted(set(n_hits_diag[valid].tolist())):
            if k == 1:
                continue
            m = valid & (n_hits_diag == k)
            cnt = int(m.sum())
            if cnt < 30:
                continue
            med = float(np.median(scatter_deg[m]))
            mean = float(np.mean(scatter_deg[m]))
            mean_path = float(np.mean(path_mm[m]))
            path_note = ''
            if path_ref is not None and path_ref > 0:
                path_note = f", path {mean_path / path_ref:.2f}x nHits==1's"
            if mean > 0 and med > 0 and (mean / med > 2.0 or mean / med < 0.5):
                print(f'    NOTE: nHitsPrimary=={k} subgroup (n={cnt}) has mean/median scatterAngleDeg = {mean:.4f}/{med:.4f} deg (ratio {mean / med:.2f}{path_note}) - consistent with Geant4 stepping-granularity noise (see the FIX comment above this block), not extra material. Do not re-introduce an nHitsPrimary-based selection here without checking this first.')
    scatter_rad = np.deg2rad(scatter_deg[valid])
    x_over_x0 = x_over_x0_all[valid]
    x_over_x0_floored = np.maximum(x_over_x0, 0.001)
    n_floored = int(np.sum(x_over_x0 < 0.001))
    if n_floored > 0:
        print(f'    NOTE: {n_floored} event(s) had X/X0 < 1e-3 (very short/grazing path) - floored to 1e-3 for the Highland prediction only; theta_0 extraction from scatterAngleDeg itself is unaffected.')
    sim_fit = fit_theta0(scatter_rad)
    highland_theta0_per_event_rad = highland_theta0_rad(momentum_GeV, x_over_x0_floored)
    highland_theta0_mean_rad = float(np.mean(highland_theta0_per_event_rad)) if n_valid > 0 else None
    tail_fit = None
    if n_valid >= 200:
        tail_fit = fit_gaussian_core_powerlaw_tail(scatter_rad, n_bootstrap=n_bootstrap)
    tail_fit_json = None
    if tail_fit is not None:
        tail_fit_json = {k: v for k, v in tail_fit.items() if not k.startswith('_')}
        print(f"    -> tail fit: sigma0={tail_fit['sigma0_mrad']:.4f}" + (f"+/-{tail_fit['sigma0_err_mrad']:.4f}" if tail_fit['sigma0_err_mrad'] is not None else '') + f" mrad   f_tail={tail_fit['f_tail']:.4f}" + (f"+/-{tail_fit['f_tail_err']:.4f}" if tail_fit['f_tail_err'] is not None else '') + f"   alpha={tail_fit['alpha']:.3f}" + (f"+/-{tail_fit['alpha_err']:.3f}" if tail_fit['alpha_err'] is not None else '') + (f"   chi2/ndof={tail_fit['chi2_per_ndof']:.2f}" if tail_fit['chi2_per_ndof'] is not None else ''))
        if make_tail_plots:
            os.makedirs(tail_plot_dir, exist_ok=True)
            safe_name = re.sub('[^A-Za-z0-9_.-]+', '_', name)
            tag = f'{safe_name}_p{momentum_GeV:g}GeV'
            try:
                plot_tail_fit_linear(scatter_rad, tail_fit, name, momentum_GeV, os.path.join(tail_plot_dir, f'{tag}_tailfit_linear.png'))
                plot_tail_fit_loglog(scatter_rad, tail_fit, name, momentum_GeV, os.path.join(tail_plot_dir, f'{tag}_tailfit_loglog.png'))
            except Exception as e:
                print(f'    NOTE: tail-fit plotting failed for this run ({e}) - the fit result above is unaffected.')
    result = {'geometry': name, 'momentum_GeV': momentum_GeV, 'n_events_total': len(hit), 'n_events_hit': int(hit.sum()), 'n_events_valid_for_scattering': n_valid, 'differentiated_scoring_used': has_split, 'mean_x_over_x0': float(np.mean(x_over_x0)) if n_valid > 0 else None, 'sim_theta0': sim_fit, 'highland_theta0_mean_mrad': highland_theta0_mean_rad * 1000.0 if highland_theta0_mean_rad is not None else None, 'sim_over_highland_ratio': None, 'tail_fit': tail_fit_json}
    sim_best = sim_fit['gaussian_fit_rad'] if sim_fit['gaussian_fit_rad'] is not None else sim_fit['rms_rad']
    if sim_best is not None and highland_theta0_mean_rad not in (None, 0.0):
        result['sim_over_highland_ratio'] = float(sim_best / highland_theta0_mean_rad)
    return result

def plot_theta0_vs_momentum(summaries, out_path='theta0_vs_momentum.png'):
    by_geometry = {}
    for s in summaries:
        if s is None:
            continue
        by_geometry.setdefault(s['geometry'], []).append(s)
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.tab10.colors
    for i, (geom, points) in enumerate(sorted(by_geometry.items())):
        points = sorted(points, key=lambda s: s['momentum_GeV'])
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
    plt.savefig(out_path, dpi=200)
    print(f'Saved {out_path}')

def plot_sim_over_highland(summaries, out_path='theta0_ratio_vs_momentum.png'):
    by_geometry = {}
    for s in summaries:
        if s is None or s['sim_over_highland_ratio'] is None:
            continue
        by_geometry.setdefault(s['geometry'], []).append(s)
    if not by_geometry:
        print(f'  NOTE: no valid sim/Highland ratios to plot - skipping {out_path}.')
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.tab10.colors
    for i, (geom, points) in enumerate(sorted(by_geometry.items())):
        points = sorted(points, key=lambda s: s['momentum_GeV'])
        p_vals = [pt['momentum_GeV'] for pt in points]
        ratios = [pt['sim_over_highland_ratio'] for pt in points]
        ax.plot(p_vals, ratios, 'o-', color=colors[i % len(colors)], label=geom)
    ax.axhline(1.0, color='black', linewidth=0.8, linestyle=':')
    ax.set_xscale('log')
    ax.set_xlabel('Momentum p (GeV/c)')
    ax.set_ylabel('$\\theta_0^{sim} / \\theta_0^{Highland}$')
    ax.set_title('Experiment 2: GEANT4 vs. Highland agreement')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f'Saved {out_path}')

def plot_tail_fit_linear(theta_rad, fit, geometry, momentum_GeV, out_path):
    theta_mrad = np.asarray(theta_rad, dtype=float) * 1000.0
    sigma0_mrad = fit['sigma0_mrad']
    f_tail = fit['f_tail']
    alpha = fit['alpha']
    theta_min_mrad = fit['theta_min_mrad']
    boot = fit.get('_bootstrap_params')
    x_max = float(np.percentile(theta_mrad, 99.0)) * 1.3
    x_max = max(x_max, theta_min_mrad * 1.5)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.hist(theta_mrad[theta_mrad <= x_max], bins=110, range=(0.0, x_max), density=True, color='#6f9fd8', alpha=0.55, edgecolor='none', label='Simulated data (space angle)')
    grid_mrad = np.linspace(1e-06, x_max, 600)
    grid_rad = grid_mrad / 1000.0
    fit_curve = mixture_pdf(grid_rad, fit['sigma0_rad'], f_tail, alpha, fit['theta_min_factor']) / 1000.0
    core_curve = rayleigh_pdf(grid_rad, fit['sigma0_rad']) / 1000.0
    if boot is not None and len(boot) >= 20:
        lo, hi = _mixture_curve_band(grid_rad, boot, fit['theta_min_factor'])
        ax.fill_between(grid_mrad, lo / 1000.0, hi / 1000.0, color='crimson', alpha=0.15, label='Fit 68% band (bootstrap)')
    ax.plot(grid_mrad, fit_curve, color='crimson', lw=2.2, label=f'Core+tail fit ($\\sigma_0$={sigma0_mrad:.3f} mrad, f={f_tail:.3f}, $\\alpha$={alpha:.2f})')
    ax.plot(grid_mrad, core_curve, color='crimson', lw=1.3, ls='--', alpha=0.75, label='Rayleigh core only (Highland-equivalent)')
    ax.axvline(theta_min_mrad, color='dimgray', lw=1.0, ls=':', label=f'$\\theta_{{min}}=3\\sigma_0$={theta_min_mrad:.3f} mrad')
    ax.set_xlim(0.0, x_max)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel('Space scattering angle [mrad]')
    ax.set_ylabel('Probability density')
    ax.set_title(f'{geometry}, p={momentum_GeV:g} GeV/c: linear scale - core region')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f'    Saved {out_path}')

def plot_tail_fit_loglog(theta_rad, fit, geometry, momentum_GeV, out_path):
    theta_mrad = np.asarray(theta_rad, dtype=float) * 1000.0
    n = len(theta_mrad)
    theta_min_mrad = fit['theta_min_mrad']
    lo_edge = max(float(np.percentile(theta_mrad, 0.05)), 0.0001)
    hi_edge = float(np.percentile(theta_mrad, 99.9))
    edges = np.geomspace(lo_edge, hi_edge, 46)
    counts, _ = np.histogram(theta_mrad, bins=edges)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    density = counts / (n * widths)
    density_err = np.sqrt(counts) / (n * widths)
    keep = counts > 0
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.errorbar(centers[keep], density[keep], yerr=density_err[keep], fmt='o', ms=4, color='#2f6fab', ecolor='#2f6fab', elinewidth=1.0, capsize=2.5, label='Simulated data (Poisson errors)')
    grid_mrad = np.geomspace(lo_edge, hi_edge, 500)
    grid_rad = grid_mrad / 1000.0
    fit_curve = mixture_pdf(grid_rad, fit['sigma0_rad'], fit['f_tail'], fit['alpha'], fit['theta_min_factor']) / 1000.0
    core_curve = rayleigh_pdf(grid_rad, fit['sigma0_rad']) / 1000.0
    boot = fit.get('_bootstrap_params')
    if boot is not None and len(boot) >= 20:
        blo, bhi = _mixture_curve_band(grid_rad, boot, fit['theta_min_factor'])
        ax.fill_between(grid_mrad, blo / 1000.0, bhi / 1000.0, color='crimson', alpha=0.15)
    ax.plot(grid_mrad, fit_curve, color='crimson', lw=2.2, label='Core+tail fit')
    ax.plot(grid_mrad, core_curve, color='crimson', lw=1.3, ls='--', alpha=0.75, label='Highland/Rayleigh only')
    ax.axvline(theta_min_mrad, color='dimgray', lw=1.0, ls=':')
    chi2_note = ''
    if fit.get('chi2_per_ndof') is not None:
        chi2_note = f"   $\\chi^2$/ndof={fit['chi2_per_ndof']:.2f}"
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=float(np.min(density[keep])) * 0.3, top=float(np.max(density[keep])) * 3.0)
    ax.set_xlabel('Space scattering angle [mrad] (log)')
    ax.set_ylabel('Probability density (log)')
    ax.set_title(f'{geometry}, p={momentum_GeV:g} GeV/c: log-log - Rutherford tail vs. Highland-only{chi2_note}')
    ax.legend(fontsize=8, loc='lower left')
    ax.grid(True, which='both', alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f'    Saved {out_path}')
if __name__ == '__main__':
    args = sys.argv[1:]
    if '--selftest' in args:
        print('Running fit_theta0() regression self-test against a known-truth synthetic Rayleigh sample...\n')
        ok = True
        for true_mrad in (0.05, 0.34, 1.4):
            ok = _synthetic_rayleigh_selftest(true_theta0_rad=true_mrad / 1000.0) and ok
            print()
        print('Running fit_gaussian_core_powerlaw_tail() regression self-test against a known-truth synthetic mixture sample...\n')
        for true_sigma0_mrad, true_f, true_alpha in ((0.34, 0.03, 3.0), (0.34, 0.08, 2.5)):
            ok = _synthetic_mixture_selftest(true_sigma0_rad=true_sigma0_mrad / 1000.0, true_f=true_f, true_alpha=true_alpha) and ok
            print()
        sys.exit(0 if ok else 1)
    tail_plot_dir = 'tail_fit_plots'
    make_tail_plots = True
    n_bootstrap = 200
    file_args = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == '--tail-plots-dir':
            tail_plot_dir = args[i + 1]
            i += 2
        elif a == '--no-tail-plots':
            make_tail_plots = False
            i += 1
        elif a == '--n-bootstrap':
            n_bootstrap = int(args[i + 1])
            i += 2
        else:
            file_args.append(a)
            i += 1
    if not file_args:
        print('Usage: python3 analyze_exp2.py results_exp2/*.root [--tail-plots-dir DIR] [--no-tail-plots] [--n-bootstrap N]')
        print('       python3 analyze_exp2.py --selftest')
        sys.exit(1)
    root_files = []
    for pattern in file_args:
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
    if not HAVE_SCIPY:
        print('NOTE: scipy not found - falling back to RMS-only theta_0 estimates (no Gaussian-core fit), and the core+tail mixture fit/plots will be skipped entirely. pip install scipy --break-system-packages for both.')
    if make_tail_plots and HAVE_SCIPY:
        print(f"Tail-fit plots will be written to '{tail_plot_dir}/' ({n_bootstrap} bootstrap replicates per run).")
    all_summaries = []
    for path in root_files:
        print(f'\nLoading {path}...')
        data = load_run(path)
        name = resolve_geometry_name(path, data)
        mom_check = validate_momentum(path, data)
        if mom_check is None:
            print(f"  WARNING: filename has no '_pNGeV' momentum tag - this file may not be Experiment 2 output. Skipping.")
            continue
        momentum_GeV = mom_check['claimed_GeV']
        print(f"  -> geometry: '{name}'   momentum: {momentum_GeV} GeV/c")
        summary = summarize_scattering(name, momentum_GeV, path, data, make_tail_plots=make_tail_plots, tail_plot_dir=tail_plot_dir, n_bootstrap=n_bootstrap)
        if summary is None:
            continue
        all_summaries.append(summary)
        sim = summary['sim_theta0']
        best_label = 'gaussian-fit' if sim['gaussian_fit_mrad'] is not None else 'RMS'
        best_val = sim['gaussian_fit_mrad'] if sim['gaussian_fit_mrad'] is not None else sim['rms_mrad']
        print(f"  -> theta_0: sim({best_label})={best_val:.3f} mrad   Highland={summary['highland_theta0_mean_mrad']:.3f} mrad   ratio={summary['sim_over_highland_ratio']}")
    if not all_summaries:
        print('\nERROR: no valid Experiment 2 summaries produced - check warnings above (missing scatterAngleDeg column, no momentum tag, etc).')
        sys.exit(1)
    with open('exp2_scattering_summary.json', 'w') as f:
        json.dump(all_summaries, f, indent=2)
    print(f'\nWrote exp2_scattering_summary.json ({len(all_summaries)} run(s))')
    plot_theta0_vs_momentum(all_summaries)
    plot_sim_over_highland(all_summaries)
    print("\nNOTE: Highland's formula (and the space-angle -> projected-angle conversion this script applies, theta0_from_space_angle()) are SMALL-ANGLE approximations. At the low end of the sweep (p=0.5 GeV/c) scattering angles are largest and both approximations are least reliable - expect sim_over_highland_ratio to drift from 1.0 there; this is expected physics, not a simulation bug.")
    print("\nNOTE: this script previously restricted the theta_0 extraction to nHitsPrimary>=2 events. That filter was removed (see the FIX comment in summarize_scattering()) after tracing EventAction.cc directly: nHitsPrimary tracks GEANT4 STEPPING GRANULARITY, not distinct material crossings (mean silicon path length is essentially unchanged between nHitsPrimary==1 and ==2 events on the exp2 barrel p=2 GeV/c reference run). Events that happen to get an extra step boundary for the SAME physical path pick up an extra condensed-history MSC angle resample at that boundary, inflating scatterAngleDeg with simulation noise rather than real scattering - this, not a real physical departure from Highland, was the dominant cause of the large, run-to-run-unstable (and sign-flipping) sim/Highland ratio in earlier results. The proper fix is upstream (the step-size limiter or geometry-safety stepping in DetectorConstruction.cc / the physics list that's causing the extra step boundary in the first place) - this script's fix only avoids amplifying that artifact at analysis time. Run with --selftest to check fit_theta0() AND fit_gaussian_core_powerlaw_tail() against known-truth synthetic samples, and watch for the nHitsPrimary diagnostic NOTE above if this dependency reappears in future data.")
    print(f"\nNOTE: the core+tail mixture fit (sigma0, f_tail, alpha per run) is written into exp2_scattering_summary.json under each run's 'tail_fit' key, and per-run diagnostic plots (linear 'core' view + log-log 'tail' view, both with error bars/bands) are saved to '{tail_plot_dir}/'. theta_min=3*sigma0 is DERIVED per-run from that run's own fitted sigma0, not a fixed absolute angle - don't compare theta_min_mrad values directly across momentum points without accounting for that.")