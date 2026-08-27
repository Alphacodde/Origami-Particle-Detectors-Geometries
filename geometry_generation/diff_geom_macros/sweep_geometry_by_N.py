import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _solid_export import thicken_and_export_stack
from _differentiated_export import thicken_and_export_differentiated_stack, weld_report
from _geometry_config import SHARED_R_MM, SHARED_L_MM, SI_THICKNESS_MM, KAPTON_THICKNESS_MM, center_verts_z
import kresling, yoshimura, miura_cylindrical, comparison_barrel
N_VALUES = [6, 8, 12, 16, 24, 32]
OUT_DIR = 'geometry_N_sweep'
USE_DIFFERENTIATED_KAPTON = True

def gen_kresling(N):
    layers = kresling.CONFIG['layers']
    H = SHARED_L_MM / layers
    twist_rad = math.radians(kresling.CONFIG['twist_deg'])
    V, F = kresling.build_kresling(N, layers, SHARED_R_MM, H, twist_rad)
    return (center_verts_z(V.reshape(-1, 3)), F)

def gen_yoshimura(N):
    m_rings = yoshimura.CONFIG['m_rings']
    verts, faces = yoshimura.build_yoshimura(N, m_rings, SHARED_R_MM, SHARED_L_MM, yoshimura.CONFIG['fold_angle'])
    return (center_verts_z(verts), faces)

def gen_miura(N):
    rows = miura_cylindrical.CONFIG['rows']
    verts, faces = miura_cylindrical.build_miura(SHARED_R_MM, SHARED_L_MM, rows, N, miura_cylindrical.CONFIG['crease_angle'])
    return (center_verts_z(verts), faces)

def gen_barrel(N):
    rows = comparison_barrel.CONFIG['rows']
    verts, faces = comparison_barrel.build_barrel(SHARED_R_MM, SHARED_L_MM, rows, N)
    return (center_verts_z(verts), faces)
SWEPT_GENERATORS = {'kresling': gen_kresling, 'yoshimura': gen_yoshimura, 'miura': gen_miura}

def export_pair(structure, tag, verts, faces):
    si_path = f'{OUT_DIR}/{tag}_silicon.gdml'
    kapton_path = f'{OUT_DIR}/{tag}_kapton.gdml'
    if structure != 'barrel' and USE_DIFFERENTIATED_KAPTON:
        n_in, n_out, n_welded = weld_report(verts, faces)
        print(f'  [{tag}] weld report: {n_in} -> {n_out} verts ({n_welded} welded)')
        thicken_and_export_differentiated_stack(verts, faces, si_thickness_mm=SI_THICKNESS_MM, si_path=si_path, si_label=f'{tag} silicon', kapton_thickness_mm=KAPTON_THICKNESS_MM, kapton_path=kapton_path, kapton_label=f'{tag} kapton')
    else:
        thicken_and_export_stack(verts, faces, [(SI_THICKNESS_MM, si_path, f'{tag} silicon'), (KAPTON_THICKNESS_MM, kapton_path, f'{tag} kapton')])

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for structure, gen_fn in SWEPT_GENERATORS.items():
        for N in N_VALUES:
            verts, faces = gen_fn(N)
            tag = f'{structure}_N{N}'
            export_pair(structure, tag, verts, faces)
    verts, faces = gen_barrel(comparison_barrel.CONFIG['cols'])
    export_pair('barrel', 'barrel_reference', verts, faces)
if __name__ == '__main__':
    main()