import math
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _solid_export import thicken_and_export_stack
from _differentiated_export import thicken_and_export_differentiated_stack, weld_report
from _geometry_config import SHARED_R_MM, SHARED_L_MM, SI_THICKNESS_MM, KAPTON_THICKNESS_MM, center_verts_z, UNIT_CELL_TARGET_MM, derive_repeat_count, OUTPUT_DIR
try:
    import FreeCAD, Mesh
    IN_FREECAD = True
except ImportError:
    IN_FREECAD = False
_CIRCUMFERENCE_MM = 2.0 * math.pi * SHARED_R_MM
_ROWS = derive_repeat_count(SHARED_L_MM, UNIT_CELL_TARGET_MM)
_COLS_FEA = derive_repeat_count(_CIRCUMFERENCE_MM, UNIT_CELL_TARGET_MM)
_CIRCLE_SAGITTA_TOL_MM = 0.01

def _min_segments_for_sagitta(radius_mm: float, tol_mm: float) -> int:
    if tol_mm >= radius_mm:
        return 3
    max_half_angle = math.acos(1.0 - tol_mm / radius_mm)
    n = math.ceil(math.pi / max_half_angle)
    return max(n, 3)
_COLS_SMOOTH = _min_segments_for_sagitta(SHARED_R_MM, _CIRCLE_SAGITTA_TOL_MM)
_COLS = max(_COLS_FEA, _COLS_SMOOTH)
CONFIG = dict(R=SHARED_R_MM, L=SHARED_L_MM, rows=_ROWS, cols=_COLS, cols_fea=_COLS_FEA, cols_smooth=_COLS_SMOOTH, sagitta_tol_mm=_CIRCLE_SAGITTA_TOL_MM, si_thickness_mm=SI_THICKNESS_MM, kapton_thickness_mm=KAPTON_THICKNESS_MM, output_tag='reference', use_differentiated_kapton=False)

def build_barrel(R, L, rows, cols):
    verts = []
    idx = {}
    n = 0
    for r in range(rows + 1):
        z = L / rows * r
        for c in range(cols):
            theta = 2.0 * math.pi * c / cols
            x = R * math.cos(theta)
            y = R * math.sin(theta)
            verts.append([x, y, z])
            idx[r, c] = n
            n += 1
    faces = []
    for r in range(rows):
        for c in range(cols):
            cp = (c + 1) % cols
            a = idx[r, c]
            b = idx[r, cp]
            cc = idx[r + 1, c]
            d = idx[r + 1, cp]
            faces.append((a, b, d))
            faces.append((a, d, cc))
    return (np.array(verts, dtype=float), faces)

def build_freecad_mesh(verts, faces):
    tris = [(FreeCAD.Vector(*verts[f[0]]), FreeCAD.Vector(*verts[f[1]]), FreeCAD.Vector(*verts[f[2]])) for f in faces]
    mesh = Mesh.Mesh(tris)
    doc = FreeCAD.newDocument('ComparisonBarrel')
    obj = doc.addObject('Mesh::Feature', 'ComparisonBarrel')
    obj.Mesh = mesh
    doc.recompute()
    print(f'[OK] FreeCAD doc (visualization only, Z-centered): {doc.Name}')

def print_stats(verts, faces, cfg):
    print('\n══════════════════════════════════════════')
    print('  COMPARISON BARREL  —  Geometry Report')
    print('══════════════════════════════════════════')
    max_sagitta = cfg['R'] * (1.0 - math.cos(math.pi / cfg['cols']))
    print(f"  Radius R      : {cfg['R']} mm  (SHARED across all structures)")
    print(f"  Length L      : {cfg['L']} mm  (SHARED across all structures)")
    print(f"  Mesh rows/cols: {cfg['rows']} / {cfg['cols']}")
    print(f"    cols if FEA-panel-parity only : {cfg['cols_fea']}  (would be a visibly faceted {cfg['cols_fea']}-gon)")
    print(f"    cols required for smooth circle: {cfg['cols_smooth']}  (sagitta tol {cfg['sagitta_tol_mm']} mm)")
    print(f'    actual worst-case chordal deviation from true R: {max_sagitta:.5f} mm')
    print(f'  Fold ratio Rf : 1.0  (no folding - stowed == deployed)')
    print(f'  Gaussian K    : 0  (developable, like a plain cylinder)')
    print(f'  Vertices      : {len(verts)}')
    print(f'  Triangles     : {len(faces)}')
    print(f"  Si thickness  : {cfg['si_thickness_mm']} mm")
    print(f"  Kapton thick. : {cfg['kapton_thickness_mm']} mm")
    print(f"  Kapton mode   : {('DIFFERENTIATED (per-vertex, sanity-check only)' if cfg['use_differentiated_kapton'] else 'SAME-MESH (per-facet, default -- no seams to bridge)')}")
    if cfg['use_differentiated_kapton']:
        n_in, n_out, n_welded = weld_report(verts, faces)
        print(f'  Weld report   : {n_in} -> {n_out} vertices ({n_welded} welded before Kapton offset)')
        print(f'    Expected ~0 for the barrel (already shares indices correctly) --')
        print(f'    a large weld count here would be unexpected and worth investigating.')
    print(f'  Z range (post-centering): [{verts[:, 2].min():.2f}, {verts[:, 2].max():.2f}] mm')
    print('══════════════════════════════════════════\n')

def main():
    cfg = CONFIG.copy()
    verts, faces = build_barrel(cfg['R'], cfg['L'], cfg['rows'], cfg['cols'])
    verts = center_verts_z(verts)
    print_stats(verts, faces, cfg)
    tag = cfg['output_tag']
    si_path = f'{OUTPUT_DIR}/barrel_{tag}_silicon.gdml'
    kapton_path = f'{OUTPUT_DIR}/barrel_{tag}_kapton.gdml'
    if cfg['use_differentiated_kapton']:
        thicken_and_export_differentiated_stack(verts, faces, si_thickness_mm=cfg['si_thickness_mm'], si_path=si_path, si_label='silicon', kapton_thickness_mm=cfg['kapton_thickness_mm'], kapton_path=kapton_path, kapton_label='kapton')
    else:
        thicken_and_export_stack(verts, faces, [(cfg['si_thickness_mm'], si_path, 'silicon'), (cfg['kapton_thickness_mm'], kapton_path, 'kapton')])
    if IN_FREECAD:
        build_freecad_mesh(verts, faces)
    else:
        print('[Sample vertices — Z-centered, post-thickening-input]')
        for i, v in enumerate(verts[:6]):
            print(f'  V[{i}] = ({v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f})')
if __name__ == '__main__':
    main()