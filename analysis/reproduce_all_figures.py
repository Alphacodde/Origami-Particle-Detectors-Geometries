import argparse
import subprocess
import sys
from pathlib import Path
REPO_ROOT = Path(__file__).parent.parent

def run(cmd):
    print(f"  >> {' '.join((str(c) for c in cmd))}")
    result = subprocess.run(cmd, check=True)
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-root', default='.', help='Path to the directory containing results_diff_geom/, results_exp2/, etc.')
    args = parser.parse_args()
    data = Path(args.data_root).resolve()
    scripts = REPO_ROOT / 'analysis'
    figures = REPO_ROOT / 'paper' / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    print('\n=== Reproducing all publication figures ===\n')
    run([sys.executable, scripts / 'pareto.py', '--results', str(data / 'results_diff_geom'), '--n-sweep', str(data / 'results_N_sweep'), '--theta-sweep', str(data / 'results_theta_sweep'), '--out', str(figures)])
    run([sys.executable, scripts / 'scatter_fit.py', '--results', str(data / 'results_exp2'), '--out', str(figures)])
    run([sys.executable, scripts / 'deadzone_map.py', '--results', str(data / 'results_diff_geom'), '--out', str(figures)])
    run([sys.executable, scripts / 'deadzone_threshold_sweep.py', '--results', str(data / 'results_diff_geom'), '--out', str(figures)])
    run([sys.executable, scripts / 'analyze_exp3.py', '--results', str(data / 'results_exp3'), '--out', str(figures)])
    run([sys.executable, scripts / 'make_its3_figure.py', '--out', str(figures)])
    print(f'\n=== Done. Figures written to {figures} ===')
if __name__ == '__main__':
    main()