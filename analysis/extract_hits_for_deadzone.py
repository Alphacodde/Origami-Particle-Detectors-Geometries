import sys
import re
import argparse
import numpy as np
import pandas as pd
import uproot
SILICON_X0_MM = 93.7
TREE_NAME = 'PionEvents'
RUN_PATTERN = re.compile('run_?(\\d+(?:_\\d+)?)')
LOCAL_PATTERNS = [('^u.*local', '^v.*local', 'facet.*id'), ('^local.*u', '^local.*v', 'facet.*id')]
GLOBAL_PATTERNS = [('^entry.*x_mm$', '^entry.*y_mm$', '^entry.*z_mm$'), ('x_mm$', 'y_mm$', 'z_mm$'), ('^hit.?x', '^hit.?y', '^hit.?z'), ('^x$', '^y$', '^z$')]
NON_HIT_POSITION_PATTERNS = ['^vertex']

def find_matching_branches(branch_names, pattern_sets):
    candidates = [b for b in branch_names if not any((re.search(p, b, re.IGNORECASE) for p in NON_HIT_POSITION_PATTERNS))]
    for patterns in pattern_sets:
        matches = []
        for pat in patterns:
            hit = [b for b in candidates if re.search(pat, b, re.IGNORECASE)]
            if len(hit) == 1:
                matches.append(hit[0])
            else:
                matches = None
                break
        if matches:
            return matches
    return None

def parse_run(root_path, run_override=None):
    if run_override is not None:
        return run_override
    m = RUN_PATTERN.search(root_path)
    return m.group(1) if m else None

def process_file(root_path, run_override=None):
    try:
        run = parse_run(root_path, run_override)
        if run is None:
            print(f"ERROR: couldn't parse a run number out of '{root_path}' (expected e.g. '..._run_3...' or '..._run0_2...'). Pass --run to set it explicitly for this invocation.")
            return 'no_run'
        with uproot.open(root_path) as f:
            if TREE_NAME not in f:
                print(f"WARNING: '{TREE_NAME}' tree not found. Keys available: {f.keys()}")
                return 'no_tree'
            tree = f[TREE_NAME]
            branch_names = tree.keys()
            print(f'=== Branches in {root_path} : {TREE_NAME} ===')
            for b in branch_names:
                print(f'  {b}')
            print()
            local_match = find_matching_branches(branch_names, LOCAL_PATTERNS)
            global_match = find_matching_branches(branch_names, GLOBAL_PATTERNS)
            if local_match:
                print(f'FOUND local facet coordinates: {local_match}')
                print('  -> preferred: resolves hinge lines cleanly, no unrolling needed.')
                coord_type, cols = ('local', local_match)
            elif global_match:
                print(f'FOUND global position branches: {global_match}')
                print('  -> usable, deadzone_map.py will unroll (theta, z) from these.')
                coord_type, cols = ('global', global_match)
            else:
                print('NO position branches found under any expected naming pattern.')
                print()
                print("This means position isn't currently being written to the ntuple")
                print('at all -- not a read-side gap, a write-side one. The GEANT4')
                print('sensitive detector (referenced as SensorSD.cc in your project)')
                print('has the hit position available at the moment it registers a hit')
                print('(e.g. via step->GetPreStepPoint()->GetPosition()) -- it just')
                print("isn't currently being pushed into the ntuple as a branch.")
                print()
                print("Fix is on the C++ side: in SensorSD.cc's ProcessHits (or")
                print('wherever hitDetector/totalPathLength_mm/etc. get filled),')
                print('add three more analysis-manager columns (e.g. hitX_mm,')
                print('hitY_mm, hitZ_mm) filled from the hit position, matching')
                print('however totalPathLength_mm etc. are currently being filled.')
                print('Share SensorSD.cc if you want help wiring that in directly.')
                return 'no_position'
            has_structure_tag = 'structureTag' in branch_names
            if not has_structure_tag:
                print("WARNING: 'structureTag' branch not found -- falling back to 'unknown' for the structure column. Expected on the updated PionEvents schema; check if this file predates that change.")
            read_cols = cols + ['hitDetector', 'totalPathLength_mm']
            if has_structure_tag:
                read_cols.append('structureTag')
            data = tree.arrays(read_cols, library='np')
        hit = data['hitDetector'].astype(bool)
        path_X0 = data['totalPathLength_mm'][hit] / SILICON_X0_MM
        if has_structure_tag:
            structure_per_hit = data['structureTag'][hit]
            uniq = np.unique(structure_per_hit)
            if len(uniq) > 1:
                print(f"WARNING: multiple structureTag values in one file {list(uniq)} -- this shouldn't happen for a single scan run. Using the first hit's tag for the whole file; check your .root generation.")
            structure = str(uniq[0]) if len(uniq) else 'unknown'
        else:
            structure = 'unknown'
        base = {'structure': structure, 'run': run}
        if coord_type == 'global':
            xcol, ycol, zcol = cols
            df = pd.DataFrame({**base, 'x': data[xcol][hit], 'y': data[ycol][hit], 'z': data[zcol][hit], 'path_X0': path_X0})
        else:
            ucol, vcol, facetcol = cols
            df = pd.DataFrame({**base, 'u_local': data[ucol][hit], 'v_local': data[vcol][hit], 'facet_id': data[facetcol][hit], 'path_X0': path_X0})
        out_path = root_path.replace('.root', '_hits.csv')
        df.to_csv(out_path, index=False)
        print(f"\nWrote {len(df)} hit rows (structure='{structure}', run='{run}') -> {out_path}")
        print('Feed this into deadzone_map.py as the fold-geometry hits file')
        print('(and run this same script on your barrel_reference.root for the baseline).')
        return 'ok'
    except Exception as exc:
        print(f'ERROR processing {root_path}: {exc}')
        return 'error'

def main():
    parser = argparse.ArgumentParser(description='Extract per-hit CSVs from PionEvents ROOT ntuples for deadzone_map.py.')
    parser.add_argument('root_paths', nargs='+', help='One or more .root files to process.')
    parser.add_argument('--run', default=None, help='Override the run number for ALL files in this invocation, instead of parsing it from each filename. Use only when processing files one run at a time -- for a mixed batch, rename files or run separately per run.')
    args = parser.parse_args()
    root_paths = args.root_paths
    batch = len(root_paths) > 1
    results = {}
    for i, root_path in enumerate(root_paths, start=1):
        if batch:
            print(f"\n{'=' * 70}")
            print(f'[{i}/{len(root_paths)}] {root_path}')
            print('=' * 70)
        results[root_path] = process_file(root_path, run_override=args.run)
    if batch:
        print(f"\n{'=' * 70}")
        print('BATCH SUMMARY')
        print('=' * 70)
        for root_path, status in results.items():
            label = {'ok': 'OK - hits CSV written', 'no_tree': f'SKIPPED - no {TREE_NAME} tree', 'no_position': 'SKIPPED - no position branches', 'no_run': "SKIPPED - couldn't parse run number", 'error': 'FAILED - see error above'}[status]
            print(f'  {root_path}: {label}')
        n_ok = sum((1 for s in results.values() if s == 'ok'))
        n_failed = sum((1 for s in results.values() if s in ('error', 'no_tree', 'no_run')))
        print(f'\n{n_ok}/{len(root_paths)} file(s) produced a hits CSV.')
        if n_failed:
            sys.exit(1)
    elif results[root_paths[0]] != 'ok':
        sys.exit(1)
if __name__ == '__main__':
    main()