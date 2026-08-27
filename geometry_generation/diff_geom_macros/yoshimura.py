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
_M_RINGS = derive_repeat_count(SHARED_L_MM, UNIT_CELL_TARGET_MM)
CONFIG = dict(n=8, m_rings=_M_RINGS, R=SHARED_R_MM, L=SHARED_L_MM, fold_angle=math.pi / 6, si_thickness_mm=SI_THICKNESS_MM, kapton_thickness_mm=KAPTON_THICKNESS_MM, export_path='', output_tag='deployed', use_differentiated_kapton=True)

def build_yoshimura(n, m_rings, R, L, fold_angle):
    dz = L / m_rings
    dr = R * (1.0 - math.cos(fold_angle))
    r_in = R - dr
    all_verts = []
    ring_idx = {}
    curr_idx = 0
    for m in range(m_rings + 1):
        for j in range(n):
            theta = 2 * math.pi * j / n
            x = R * math.cos(theta)
            y = R * math.sin(theta)
            z = m * dz
            all_verts.append([x, y, z])
            ring_idx['bot', m, j] = curr_idx
            curr_idx += 1
    for m in range(m_rings):
        for j in range(n):
            theta = 2 * math.pi * (j + 0.5) / n
            x = r_in * math.cos(theta)
            y = r_in * math.sin(theta)
            z = (m + 0.5) * dz
            all_verts.append([x, y, z])
            ring_idx['mid', m, j] = curr_idx
            curr_idx += 1
    faces = []
    for m in range(m_rings):
        for j in range(n):
            jp = (j + 1) % n
            bot_lo_l = ring_idx['bot', m, j]
            bot_lo_r = ring_idx['bot', m, jp]
            mid = ring_idx['mid', m, j]
            bot_hi_l = ring_idx['bot', m + 1, j]
            bot_hi_r = ring_idx['bot', m + 1, jp]
            faces.append([bot_lo_l, bot_lo_r, mid])
            faces.append([mid, bot_hi_l, bot_lo_l])
            faces.append([mid, bot_lo_r, bot_hi_r])
            faces.append([mid, bot_hi_r, bot_hi_l])
    return (np.array(all_verts, dtype=float), faces)

def export_obj(verts, faces, path):
    with open(path, 'w') as fh:
        fh.write('# Yoshimura Pattern\n')
        for v in verts:
            fh.write(f'v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n')
        for f in faces:
            fh.write(f'f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n')
    print(f'[OK] OBJ (raw, zero-thickness, visualization only) → {path}')

def build_freecad_mesh(verts, faces):
    tris = [(FreeCAD.Vector(*verts[f[0]]), FreeCAD.Vector(*verts[f[1]]), FreeCAD.Vector(*verts[f[2]])) for f in faces]
    mesh = Mesh.Mesh(tris)
    doc = FreeCAD.newDocument('Yoshimura')
    obj = doc.addObject('Mesh::Feature', 'Yoshimura')
    obj.Mesh = mesh
    doc.recompute()
    print(f'[OK] FreeCAD doc (visualization only, Z-centered): {doc.Name}')

def print_stats(verts, faces, cfg):
    print('\n══════════════════════════════════════════')
    print('  YOSHIMURA  —  Geometry Report')
    print('══════════════════════════════════════════')
    print(f"  Facets/ring   : {cfg['n']}")
    print(f"  Diamond rings : {cfg['m_rings']}  (DERIVED: round(SHARED_L_MM/UNIT_CELL_TARGET_MM))")
    print(f"  Radius R      : {cfg['R']} mm  (SHARED across all structures)")
    print(f"  Length L      : {cfg['L']} mm  (SHARED across all structures)")
    print(f"  Ring spacing dz: {cfg['L'] / cfg['m_rings']:.3f} mm  (ACTUAL unit-cell size after rounding - target was {UNIT_CELL_TARGET_MM} mm)")
    print(f"  Fold angle    : {math.degrees(cfg['fold_angle']):.1f}°")
    print(f"  Inset dr      : {cfg['R'] * (1 - math.cos(cfg['fold_angle'])):.2f} mm")
    print(f'  Vertices      : {len(verts)}')
    print(f'  Triangles     : {len(faces)}')
    print(f"  Si thickness  : {cfg['si_thickness_mm']} mm")
    print(f"  Kapton thick. : {cfg['kapton_thickness_mm']} mm")
    print(f"  Kapton mode   : {('DIFFERENTIATED (per-vertex, continuous shell)' if cfg['use_differentiated_kapton'] else 'SAME-MESH (per-facet, legacy)')}")
    if cfg['use_differentiated_kapton']:
        n_in, n_out, n_welded = weld_report(verts, faces)
        print(f'  Weld report   : {n_in} -> {n_out} vertices ({n_welded} welded before Kapton offset)')
    print(f'  Z range (post-centering): [{verts[:, 2].min():.2f}, {verts[:, 2].max():.2f}] mm')
    print('══════════════════════════════════════════\n')

def main():
    cfg = CONFIG.copy()
    verts, faces = build_yoshimura(cfg['n'], cfg['m_rings'], cfg['R'], cfg['L'], cfg['fold_angle'])
    verts = center_verts_z(verts)
    print_stats(verts, faces, cfg)
    if cfg['export_path']:
        export_obj(verts, faces, cfg['export_path'])
    tag = cfg['output_tag']
    si_path = f'{OUTPUT_DIR}/yoshimura_{tag}_silicon.gdml'
    kapton_path = f'{OUTPUT_DIR}/yoshimura_{tag}_kapton.gdml'
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