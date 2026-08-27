import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot
GEOMETRY_FILES = {'Barrel': 'results_diff_geom/barrel_scan_run_1.root', 'Miura': 'results_diff_mesh/miura_scan_run_1.root', 'Yoshimura': 'results_diff_mesh/yoshimura_scan_run_1.root', 'Kresling': 'results_diff_mesh/kresling_scan_run_1.root'}
BRANCH_OVERRIDES = {}
TREE_NAME = None
BEAM_MOMENTUM_GEV = 5.0
PARTICLE_MASS_GEV = 0.13957
OUT_DIR = 'ms_analysis_out'
_AXIS = {'px': 'x', 'py': 'y', 'pz': 'z'}
_MOM_WORDS = {'p', 'mom', 'momentum', 'momenta'}
_PHASE_TOKENS = {'in': {'in', 'before', 'entry', 'enter', 'start', 'initial'}, 'out': {'out', 'after', 'exit', 'leave', 'end', 'final'}}
_X0_TOKENS = {'x0', 'radlen', 'radiationlength', 'xx0'}

def _tokenize(name: str):
    parts = re.split('[^a-zA-Z0-9]+', name)
    tokens = []
    for p in parts:
        if not p:
            continue
        sub = re.findall('[A-Z]?[a-z0-9]+|[A-Z]+(?=[A-Z]|$)', p)
        tokens.extend(sub if sub else [p])
    return [t.lower() for t in tokens]

def _guess(names, component, phase):
    axis = _AXIS[component]
    phase_tokens = _PHASE_TOKENS[phase]
    hits = []
    for n in names:
        toks = set(_tokenize(n))
        if not toks & phase_tokens:
            continue
        if component in toks:
            hits.append(n)
        elif axis in toks and toks & _MOM_WORDS:
            hits.append(n)
    return hits
_X0_SUBSTRINGS = ['x0', 'radlen', 'radiationlength', 'thicknessoverx0']

def _guess_x0(names):
    hits = []
    for n in names:
        toks = _tokenize(n)
        toks_set = set(toks)
        joined = ''.join(toks)
        if toks_set & _X0_TOKENS or ('x' in toks_set and '0' in toks_set) or any((s in joined for s in _X0_SUBSTRINGS)):
            hits.append(n)
    return hits

def detect_branches(file_path, override=None):
    override = override or {}
    with uproot.open(file_path) as f:
        tree_key = TREE_NAME
        if tree_key is None:
            trees = [k for k, cls in f.classnames().items() if cls.startswith('TTree') or 'RNTuple' in cls]
            if not trees:
                raise RuntimeError(f'No TTree/RNTuple found in {file_path}. Top-level keys: {list(f.keys())}')
            tree_key = trees[0]
        tree = f[tree_key]
        names = tree.keys()
    resolved = {}
    needed = ['px_in', 'py_in', 'pz_in', 'px_out', 'py_out', 'pz_out']
    for key in needed:
        if key in override:
            resolved[key] = override[key]
            continue
        component, phase = key.split('_')
        candidates = _guess(names, component, phase)
        if len(candidates) == 1:
            resolved[key] = candidates[0]
        elif len(candidates) > 1:
            candidates.sort(key=len)
            resolved[key] = candidates[0]
            warnings.warn(f"[{file_path}] multiple candidates for '{key}': {candidates} -> picked '{resolved[key]}'. Set BRANCH_OVERRIDES if wrong.")
        else:
            raise RuntimeError(f"Could not auto-detect branch for '{key}' in {file_path}.\nAvailable branches:\n  " + '\n  '.join(names) + f'\n\nAdd an entry to BRANCH_OVERRIDES for this geometry, e.g.\nBRANCH_OVERRIDES["<name>"] = {{"px_in": "...", "py_in": "...", "pz_in": "...", "px_out": "...", "py_out": "...", "pz_out": "...", "x_over_x0": "..."}}')
    if 'x_over_x0' in override:
        resolved['x_over_x0'] = override['x_over_x0']
    else:
        x0_candidates = _guess_x0(names)
        resolved['x_over_x0'] = x0_candidates[0] if x0_candidates else None
    resolved['_tree'] = tree_key
    return resolved

def scattering_angle_mrad(px_in, py_in, pz_in, px_out, py_out, pz_out):
    p_in = np.column_stack((px_in, py_in, pz_in))
    p_out = np.column_stack((px_out, py_out, pz_out))
    dot = np.sum(p_in * p_out, axis=1)
    mag_in = np.linalg.norm(p_in, axis=1)
    mag_out = np.linalg.norm(p_out, axis=1)
    valid = (mag_in > 0) & (mag_out > 0)
    cos_theta = np.full_like(dot, np.nan, dtype=float)
    cos_theta[valid] = np.clip(dot[valid] / (mag_in[valid] * mag_out[valid]), -1.0, 1.0)
    theta_rad = np.arccos(cos_theta)
    return theta_rad * 1000.0

def highland_theta0_mrad(p_gev, mass_gev, x_over_x0):
    x_over_x0 = np.asarray(x_over_x0, dtype=float)
    beta = p_gev / np.sqrt(p_gev ** 2 + mass_gev ** 2)
    p_mev = p_gev * 1000.0
    with np.errstate(divide='ignore', invalid='ignore'):
        log_term = np.where(x_over_x0 > 0, 1 + 0.038 * np.log(x_over_x0), np.nan)
        theta0_rad = 13.6 / (beta * p_mev) * np.sqrt(x_over_x0) * log_term
    return theta0_rad * 1000.0

def load_geometry(name, file_path):
    override = BRANCH_OVERRIDES.get(name)
    branches = detect_branches(file_path, override)
    read_list = [branches[k] for k in ['px_in', 'py_in', 'pz_in', 'px_out', 'py_out', 'pz_out']]
    if branches['x_over_x0']:
        read_list.append(branches['x_over_x0'])
    with uproot.open(file_path) as f:
        tree = f[branches['_tree']]
        arrays = tree.arrays(read_list, library='np')
    px_in = arrays[branches['px_in']]
    py_in = arrays[branches['py_in']]
    pz_in = arrays[branches['pz_in']]
    px_out = arrays[branches['px_out']]
    py_out = arrays[branches['py_out']]
    pz_out = arrays[branches['pz_out']]
    x_over_x0 = None
    if branches['x_over_x0']:
        x_over_x0 = arrays[branches['x_over_x0']]
    theta_mrad = scattering_angle_mrad(px_in, py_in, pz_in, px_out, py_out, pz_out)
    good = ~np.isnan(theta_mrad)
    theta_mrad = theta_mrad[good]
    if x_over_x0 is not None:
        x_over_x0 = x_over_x0[good]
    print(f"[{name}] tree='{branches['_tree']}'  n_particles={len(theta_mrad)}  branches used: { {k: v for k, v in branches.items() if not k.startswith('_')}}")
    return {'theta_mrad': theta_mrad, 'x_over_x0': x_over_x0, 'branches': branches}

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = {}
    for geom_name, path in GEOMETRY_FILES.items():
        if not os.path.exists(path):
            print(f"[WARN] {geom_name}: file not found at '{path}' -- skipping.")
            continue
        try:
            results[geom_name] = load_geometry(geom_name, path)
        except RuntimeError as e:
            print(f"\n[ERROR] Branch auto-detection failed for '{geom_name}':\n{e}\n")
    if not results:
        print('No geometries loaded successfully. Fix paths/branches and re-run.')
        return
    rows = []
    for geom_name, data in results.items():
        theta = data['theta_mrad']
        x_over_x0 = data['x_over_x0']
        rms = float(np.sqrt(np.mean(theta ** 2)))
        std = float(np.std(theta, ddof=1))
        p68 = float(np.percentile(theta, 68))
        p95 = float(np.percentile(theta, 95))
        mean_x_x0 = float(np.mean(x_over_x0)) if x_over_x0 is not None else np.nan
        if x_over_x0 is not None:
            theta0_highland = float(np.mean(highland_theta0_mrad(BEAM_MOMENTUM_GEV, PARTICLE_MASS_GEV, x_over_x0)))
            pct_diff = 100.0 * (rms - theta0_highland) / theta0_highland
        else:
            theta0_highland = np.nan
            pct_diff = np.nan
        rows.append({'Geometry': geom_name, 'N': len(theta), 'mean_X_X0': mean_x_x0, 'theta_RMS_mrad': rms, 'theta_std_mrad': std, 'theta_p68_mrad': p68, 'theta_p95_mrad': p95, 'Highland_theta0_mrad': theta0_highland, 'pct_diff_sim_vs_Highland': pct_diff})
    summary = pd.DataFrame(rows).sort_values('Geometry')
    csv_path = os.path.join(OUT_DIR, 'summary_table.csv')
    summary.to_csv(csv_path, index=False)
    print('\n=== Summary ===')
    print(summary.to_string(index=False))
    print(f'\nSaved: {csv_path}')
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(summary['Geometry'], summary['theta_RMS_mrad'])
    ax.set_ylabel('$\\theta_{\\rm RMS}$ [mrad]')
    ax.set_title('Multiple-scattering RMS angle by fold geometry')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'theta_rms_vs_geometry.png'), dpi=150)
    plt.close(fig)
    if summary['mean_X_X0'].notna().any():
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(summary['mean_X_X0'], summary['theta_RMS_mrad'])
        for _, r in summary.iterrows():
            ax.annotate(r['Geometry'], (r['mean_X_X0'], r['theta_RMS_mrad']), textcoords='offset points', xytext=(5, 5))
        ax.set_xlabel('$X/X_0$')
        ax.set_ylabel('$\\theta_{\\rm RMS}$ [mrad]')
        ax.set_title('$\\theta_{\\rm RMS}$ vs material budget')
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, 'theta_rms_vs_X_X0.png'), dpi=150)
        plt.close(fig)
    for geom_name, data in results.items():
        theta = data['theta_mrad']
        x_over_x0 = data['x_over_x0']
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(theta, bins=100, density=True, alpha=0.7, label='GEANT4 simulation')
        if x_over_x0 is not None:
            theta0 = float(np.mean(highland_theta0_mrad(BEAM_MOMENTUM_GEV, PARTICLE_MASS_GEV, x_over_x0)))
            xs = np.linspace(0, np.percentile(theta, 99.5), 300)
            gauss = 1.0 / (theta0 * np.sqrt(2 * np.pi)) * np.exp(-xs ** 2 / (2 * theta0 ** 2))
            ax.plot(xs, gauss * 2, label=f'Highland ($\\theta_0$={theta0:.2f} mrad)')
        ax.set_xlabel('Scattering angle [mrad]')
        ax.set_ylabel('Probability density')
        ax.set_title(f'Multiple Coulomb Scattering -- {geom_name}')
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f'distribution_{geom_name}.png'), dpi=150)
        plt.close(fig)
    print(f'\nAll plots saved under ./{OUT_DIR}/')
if __name__ == '__main__':
    main()