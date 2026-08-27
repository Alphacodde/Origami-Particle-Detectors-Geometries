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
_COLS = derive_repeat_count(_CIRCUMFERENCE_MM, UNIT_CELL_TARGET_MM)
CONFIG = dict(radius=SHARED_R_MM, height=SHARED_L_MM, rows=_ROWS, cols=_COLS, crease_angle=45.0, si_thickness_mm=SI_THICKNESS_MM, kapton_thickness_mm=KAPTON_THICKNESS_MM, output_tag='deployed', use_differentiated_kapton=True)

def build_miura(radius, height, rows, cols, alpha_deg):
    alpha = math.radians(alpha_deg)
    circ = 2.0 * math.pi * radius
    cell_w = circ / cols
    cell_h = height / rows
    offset = cell_w * math.cos(alpha)
    zig_h = cell_h

    def flat_vertex(col, row, even_row):
        x = col * cell_w + (0 if even_row else offset * 0.5)
        y = row * zig_h
        return (x, y)

    def to_3d(x, y):
        theta = x / circ * 2.0 * math.pi
        return (radius * math.cos(theta), radius * math.sin(theta), y)
    verts = []
    faces = []

    def add_tri(a, b, c):
        i0 = len(verts)
        verts.append(a)
        i1 = len(verts)
        verts.append(b)
        i2 = len(verts)
        verts.append(c)
        faces.append((i0, i1, i2))
    for r in range(rows):
        even_top = r % 2 == 0
        for c in range(cols):
            bl_x, bl_y = flat_vertex(c, r, even_top)
            br_x, br_y = flat_vertex(c + 1, r, even_top)
            tr_x, tr_y = flat_vertex(c + 1, r + 1, not even_top)
            tl_x, tl_y = flat_vertex(c, r + 1, not even_top)
            bl = to_3d(bl_x, bl_y)
            br = to_3d(br_x, br_y)
            tr = to_3d(tr_x, tr_y)
            tl = to_3d(tl_x, tl_y)
            add_tri(bl, br, tr)
            add_tri(bl, tr, tl)
    return (np.array(verts), faces)

def build_freecad_mesh(verts, faces):
    tris = [(FreeCAD.Vector(*verts[f[0]]), FreeCAD.Vector(*verts[f[1]]), FreeCAD.Vector(*verts[f[2]])) for f in faces]
    mesh = Mesh.Mesh(tris)
    doc = FreeCAD.newDocument('CylindricalMiuraOri')
    obj = doc.addObject('Mesh::Feature', 'CylindricalMiuraOri')
    obj.Mesh = mesh
    doc.recompute()
    print(f'[OK] FreeCAD doc (visualization only, Z-centered): {doc.Name}')

def print_stats(verts, faces, cfg):
    print('\n══════════════════════════════════════════')
    print('  CYLINDRICAL MIURA-ORI  —  Geometry Report')
    print('══════════════════════════════════════════')
    print(f"  Radius        : {cfg['radius']} mm  (SHARED across all structures)")
    print(f"  Height        : {cfg['height']} mm  (SHARED across all structures)")
    print(f"  Rows / Cols   : {cfg['rows']} / {cfg['cols']}  (DERIVED from UNIT_CELL_TARGET_MM, both axes)")
    print(f"  Actual cell h : {cfg['height'] / cfg['rows']:.3f} mm  (target was {UNIT_CELL_TARGET_MM} mm)")
    print(f"  Actual cell w : {2 * math.pi * cfg['radius'] / cfg['cols']:.3f} mm  (target was {UNIT_CELL_TARGET_MM} mm)")
    print(f"  Crease angle  : {cfg['crease_angle']:.1f}°")
    print(f'  Vertices      : {len(verts)}  (raw, unwelded -- silicon path uses this directly)')
    print(f'  Triangles     : {len(faces)}')
    print(f"  Si thickness  : {cfg['si_thickness_mm']} mm")
    print(f"  Kapton thick. : {cfg['kapton_thickness_mm']} mm")
    print(f"  Kapton mode   : {('DIFFERENTIATED (welded + per-vertex, continuous shell)' if cfg['use_differentiated_kapton'] else 'SAME-MESH (per-facet, legacy)')}")
    if cfg['use_differentiated_kapton']:
        n_in, n_out, n_welded = weld_report(verts, faces)
        print(f'  Weld report   : {n_in} -> {n_out} vertices ({n_welded} welded before Kapton offset)')
        print(f'    NOTE: this weld count is expected to be LARGE for Miura specifically')
        print(f'    (build_miura() does not share indices at the source -- see docstring')
        print(f"    above) -- a near-zero weld count here would indicate build_miura()'s")
        print(f'    panel-corner coordinates no longer coincide as expected and should be')
        print(f'    investigated, not assumed benign.')
    print(f'  Z range (post-centering): [{verts[:, 2].min():.2f}, {verts[:, 2].max():.2f}] mm')
    print('══════════════════════════════════════════\n')

def main():
    cfg = CONFIG.copy()
    verts, faces = build_miura(cfg['radius'], cfg['height'], cfg['rows'], cfg['cols'], cfg['crease_angle'])
    verts = center_verts_z(verts)
    print_stats(verts, faces, cfg)
    tag = cfg['output_tag']
    si_path = f'{OUTPUT_DIR}/miura_{tag}_silicon.gdml'
    kapton_path = f'{OUTPUT_DIR}/miura_{tag}_kapton.gdml'
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