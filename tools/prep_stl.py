import sys
import trimesh
import numpy as np

def check_stl(path):
    print(f"\n{'=' * 70}\nChecking: {path}\n{'=' * 70}")
    mesh = trimesh.load(path, file_type='stl')
    bbox = mesh.bounds
    extent = mesh.extents
    print(f'Bounding box min: {bbox[0]}')
    print(f'Bounding box max: {bbox[1]}')
    print(f'Extents (size):   {extent}')
    print(f'Characteristic size (max extent): {extent.max():.4f} [unlabeled units]')
    print('\n--> ACTION REQUIRED: does the characteristic size above match')
    print('    what you expect in millimeters for this fold, roughly the')
    print('    physical size you designed in CAD? If not, STOP and fix units')
    print('    in your CAD tool (re-export with mm units) before continuing.')
    print(f'\nIs watertight (whole mesh): {mesh.is_watertight}  (informational only for multi-component solids - see below)')
    components = mesh.split(only_watertight=False)
    n_components = len(components)
    n_bad_components = 0
    if n_components > 1:
        for c in components:
            if not (c.is_watertight and c.is_winding_consistent):
                n_bad_components += 1
        print(f'Per-component check: {n_components} component(s), {n_components - n_bad_components} individually watertight + consistently wound, {n_bad_components} NOT')
        components_ok = n_bad_components == 0
    else:
        components_ok = mesh.is_watertight
    if not components_ok:
        print(f'  FAIL: {n_bad_components} component(s) are individually broken (not watertight or inconsistently wound) - this is a real defect, not the expected multi-prism pattern.')
        print('  Common fixes:')
        print("    - In your CAD tool, check for 'naked edges' / open boundary")
        print('      loops and stitch/heal them.')
        print('    - trimesh can attempt an automatic fix (see below) but')
        print('      ALWAYS re-inspect visually after auto-fix - it can hide')
        print('      real modeling errors rather than fix them, and does NOT')
        print('      understand the intentional multi-prism structure (it')
        print('      may try to weld components that are meant to stay')
        print('      independent) - prefer fixing the generator script over')
        print("      auto-fixing its output for this project's geometry.")
    print(f'Is winding consistent: {mesh.is_winding_consistent}')
    if not mesh.is_winding_consistent:
        print('  FAIL: inconsistent triangle winding -> inconsistent normals.')
        print('  This WILL cause GEANT4 to misidentify inside/outside for the')
        print('  solid, silently breaking the thickening step. Fix in CAD or')
        print('  let trimesh attempt mesh.fix_normals() (see below).')
    areas = mesh.area_faces
    n_degenerate = int(np.sum(areas < 1e-09))
    print(f'Degenerate (zero-area) triangles: {n_degenerate} / {len(areas)}')
    if n_degenerate > 0:
        print('  FAIL: degenerate triangles present - remove before export.')
    ok = components_ok and mesh.is_winding_consistent and (n_degenerate == 0)
    print(f"\n{('PASS' if ok else 'FAIL - fix before GEANT4 import')}")
    return (mesh, ok)

def attempt_autofix(mesh, out_path):
    mesh.fix_normals()
    trimesh.repair.fill_holes(mesh)
    mesh.remove_degenerate_faces()
    mesh.remove_duplicate_faces()
    mesh.export(out_path, file_type='stl_ascii')
    print(f'\nAuto-fix attempted. Wrote: {out_path}')
    print('RE-RUN this script on the output and visually re-inspect before use.')
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('stl_path')
    parser.add_argument('--non-interactive', action='store_true', help='Never prompt for auto-fix; exit non-zero on FAIL instead.')
    args = parser.parse_args()
    mesh, ok = check_stl(args.stl_path)
    if not ok:
        if args.non_interactive:
            print('\n--non-interactive: not offering auto-fix. Exiting non-zero.')
            sys.exit(1)
        fix_path = args.stl_path.replace('.stl', '_fixed.stl')
        resp = input(f'\nAttempt auto-fix and save to {fix_path}? [y/N]: ')
        if resp.lower() == 'y':
            attempt_autofix(mesh, fix_path)