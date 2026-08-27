import math
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _solid_export import thicken_and_export_stack
from _differentiated_export import thicken_and_export_differentiated_stack, weld_report
from _geometry_config import SHARED_R_MM, SHARED_L_MM, SI_THICKNESS_MM, KAPTON_THICKNESS_MM, center_verts_z, UNIT_CELL_TARGET_MM, derive_repeat_count, OUTPUT_DIR
try:
    import FreeCAD, Part, Mesh
    IN_FREECAD = True
except ImportError:
    IN_FREECAD = False
_LAYERS = derive_repeat_count(SHARED_L_MM, UNIT_CELL_TARGET_MM)
CONFIG = dict(n=6, layers=_LAYERS, R=SHARED_R_MM, H=SHARED_L_MM / _LAYERS, twist_deg=30.0, si_thickness_mm=SI_THICKNESS_MM, kapton_thickness_mm=KAPTON_THICKNESS_MM, export_path='', output_tag='deployed', use_differentiated_kapton=True)

def build_kresling(n, layers, R, H, twist_rad):
    V = np.zeros((layers + 1, n, 3))
    for m in range(layers + 1):
        for j in range(n):
            theta = 2 * math.pi * j / n + m * twist_rad
            V[m, j] = [R * math.cos(theta), R * math.sin(theta), m * H]

    def idx(m, j):
        return m * n + j % n
    F = []
    for m in range(layers):
        for j in range(n):
            jp = (j + 1) % n
            F.append([idx(m, j), idx(m, jp), idx(m + 1, j)])
            F.append([idx(m, jp), idx(m + 1, jp), idx(m + 1, j)])
    return (V, F)

def export_obj(V, F, path):
    flat = V.reshape(-1, 3)
    with open(path, 'w') as fh:
        fh.write('# Kresling Pattern\n')
        for v in flat:
            fh.write(f'v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n')
        for f in F:
            fh.write(f'f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n')
    print(f'[OK] OBJ → {path}')

def build_freecad_mesh(V, F):
    flat = V.reshape(-1, 3)
    tris = [(FreeCAD.Vector(*flat[f[0]]), FreeCAD.Vector(*flat[f[1]]), FreeCAD.Vector(*flat[f[2]])) for f in F]
    mesh = Mesh.Mesh(tris)
    doc = FreeCAD.newDocument('Kresling')
    obj = doc.addObject('Mesh::Feature', 'Kresling_Surface')
    obj.Mesh = mesh
    doc.recompute()
    print(f'[OK] FreeCAD doc (visualization only, Z-centered): {doc.Name}')
    return doc

def print_stats(V, F, cfg):
    flat = V.reshape(-1, 3)
    print('\n══════════════════════════════════════════')
    print('  KRESLING  —  Geometry Report')
    print('══════════════════════════════════════════')
    print(f"  Sides / layer : {cfg['n']}")
    print(f"  Layers        : {cfg['layers']}  (DERIVED: round(SHARED_L_MM/UNIT_CELL_TARGET_MM))")
    print(f"  Radius R      : {cfg['R']} mm  (SHARED across all structures)")
    print(f"  Height / unit : {cfg['H']:.3f} mm  (ACTUAL unit-cell size after rounding - target was {UNIT_CELL_TARGET_MM} mm)")
    print(f"  Twist / layer : {cfg['twist_deg']:.1f}°")
    print(f"  Total height  : {cfg['layers'] * cfg['H']:.1f} mm  (SHARED across all structures)")
    print(f'  Vertices      : {len(flat)}')
    print(f'  Triangles     : {len(F)}')
    print(f"  Si thickness  : {cfg['si_thickness_mm']} mm")
    print(f"  Kapton thick. : {cfg['kapton_thickness_mm']} mm")
    print(f"  Kapton mode   : {('DIFFERENTIATED (per-vertex, continuous shell)' if cfg['use_differentiated_kapton'] else 'SAME-MESH (per-facet, legacy)')}")
    if cfg['use_differentiated_kapton']:
        n_in, n_out, n_welded = weld_report(flat, F)
        print(f'  Weld report   : {n_in} -> {n_out} vertices ({n_welded} welded before Kapton offset)')
    print(f'  Z range (post-centering): [{flat[:, 2].min():.2f}, {flat[:, 2].max():.2f}] mm')
    print('══════════════════════════════════════════\n')

def main():
    cfg = CONFIG.copy()
    twist_rad = math.radians(cfg['twist_deg'])
    print('[Kresling] Building geometry …')
    V, F = build_kresling(cfg['n'], cfg['layers'], cfg['R'], cfg['H'], twist_rad)
    flat = center_verts_z(V.reshape(-1, 3))
    print_stats(flat.reshape(V.shape), F, cfg)
    if cfg['export_path']:
        export_obj(flat.reshape(V.shape), F, cfg['export_path'])
    tag = cfg['output_tag']
    si_path = f'{OUTPUT_DIR}/kresling_{tag}_silicon.gdml'
    kapton_path = f'{OUTPUT_DIR}/kresling_{tag}_kapton.gdml'
    if cfg['use_differentiated_kapton']:
        thicken_and_export_differentiated_stack(flat, F, si_thickness_mm=cfg['si_thickness_mm'], si_path=si_path, si_label='silicon', kapton_thickness_mm=cfg['kapton_thickness_mm'], kapton_path=kapton_path, kapton_label='kapton')
    else:
        thicken_and_export_stack(flat, F, [(cfg['si_thickness_mm'], si_path, 'silicon'), (cfg['kapton_thickness_mm'], kapton_path, 'kapton')])
    if IN_FREECAD:
        build_freecad_mesh(flat.reshape(V.shape), F)
    else:
        print('[Sample vertices — Z-centered, post-thickening-input]')
        for i, v in enumerate(flat[:6]):
            print(f'  V[{i}] = ({v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f})')
if __name__ == '__main__':
    main()