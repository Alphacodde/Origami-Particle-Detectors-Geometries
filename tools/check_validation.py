import argparse
import sys

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('root_file')
    p.add_argument('--expected-mm', type=float, required=True, help='Expected silicon thickness in mm (e.g. 0.300)')
    p.add_argument('--tolerance-frac', type=float, default=0.25, help='Allowed fractional deviation from expected (default 0.25)')
    p.add_argument('--tree', default='PionEvents', help="Ntuple/tree name (matches RunAction.cc's booked ntuple)")
    args = p.parse_args()
    try:
        import uproot
    except ImportError:
        sys.exit('uproot is required: pip install uproot --break-system-packages')
    try:
        f = uproot.open(args.root_file)
    except Exception as e:
        sys.exit(f"ERROR: could not open '{args.root_file}': {e}")
    tree_key = None
    for key in f.keys():
        if key.split(';')[0] == args.tree:
            tree_key = key
            break
    if tree_key is None:
        sys.exit(f"ERROR: tree '{args.tree}' not found in '{args.root_file}'. Available keys: {list(f.keys())}")
    tree = f[tree_key]
    hit = tree['hitDetector'].array(library='np')
    path = tree['totalPathLength_mm'].array(library='np')
    hit_mask = hit == 1
    n_hit = int(hit_mask.sum())
    if n_hit == 0:
        sys.exit("ERROR: zero events hit the plate - check gunZmm/diskRadiusMm against the plate's actual position/size. This is a geometry-miss, not a thickness-measurement failure.")
    mean_path = float(path[hit_mask].mean())
    expected = args.expected_mm
    tol = args.tolerance_frac
    lo, hi = (expected * (1 - tol), expected * (1 + tol))
    ok = lo <= mean_path <= hi
    print(f'[check_validation] {n_hit} events hit the plate')
    print(f'[check_validation] mean totalPathLength_mm = {mean_path:.5f}')
    print(f'[check_validation] expected ~ {expected:.5f} mm (pass range [{lo:.5f}, {hi:.5f}])')
    if not ok:
        print("[check_validation] FAIL - mean path length is outside the expected range. Near 0 => thickness applied along the wrong axis (the old BuildThickenedShell failure mode). Near the plate's 50mm side length => the offset went in-plane instead of along the surface normal.")
        return 1
    print('[check_validation] PASS')
    return 0
if __name__ == '__main__':
    sys.exit(main())