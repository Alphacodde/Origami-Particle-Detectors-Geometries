import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _solid_export import thicken_and_export_stack
from _differentiated_export import thicken_and_export_differentiated_stack, weld_report
from _geometry_config import SHARED_R_MM, SHARED_L_MM, SI_THICKNESS_MM, KAPTON_THICKNESS_MM, center_verts_z
import kresling, yoshimura, miura_cylindrical
THETA_DEG_VALUES = [15, 30, 45, 60, 75]
OUT_DIR = 'geometry_theta_sweep'
USE_DIFFERENTIATED_KAPTON = True

def gen_kresling(theta_deg):
    n = kresling.CONFIG['n']
    layers = kresling.CONFIG['layers']
    H = SHARED_L_MM / layers
    twist_rad = math.radians(theta_deg)
    V, F = kresling.build_kresling(n, layers, SHARED_R_MM, H, twist_rad)
    return (center_verts_z(V.reshape(-1, 3)), F)

def gen_yoshimura(theta_deg):
    n = yoshimura.CONFIG['n']
    m_rings = yoshimura.CONFIG['m_rings']
    fold_angle_rad = math.radians(theta_deg)
    verts, faces = yoshimura.build_yoshimura(n, m_rings, SHARED_R_MM, SHARED_L_MM, fold_angle_rad)
    return (center_verts_z(verts), faces)

def gen_miura(theta_deg):
    cols = miura_cylindrical.CONFIG['cols']
    rows = miura_cylindrical.CONFIG['rows']
    verts, faces = miura_cylindrical.build_miura(SHARED_R_MM, SHARED_L_MM, rows, cols, theta_deg)
    return (center_verts_z(verts), faces)
SWEPT_GENERATORS = {'kresling': gen_kresling, 'yoshimura': gen_yoshimura, 'miura': gen_miura}

def print_theta_stats(structure, theta_deg):
    if structure == 'yoshimura':
        fold_angle_rad = math.radians(theta_deg)
        dr = SHARED_R_MM * (1.0 - math.cos(fold_angle_rad))
        r_in = SHARED_R_MM - dr
        print(f'  [yoshimura] theta={theta_deg:.1f}deg -> fold_angle={fold_angle_rad:.4f} rad, dr={dr:.2f} mm, r_in={r_in:.2f} mm (R={SHARED_R_MM} mm)')
    else:
        print(f'  [{structure}] theta={theta_deg:.1f}deg')

def export_pair(tag, verts, faces):
    si_path = f'{OUT_DIR}/{tag}_silicon.gdml'
    kapton_path = f'{OUT_DIR}/{tag}_kapton.gdml'
    if USE_DIFFERENTIATED_KAPTON:
        n_in, n_out, n_welded = weld_report(verts, faces)
        print(f'  [{tag}] weld report: {n_in} -> {n_out} verts ({n_welded} welded)')
        thicken_and_export_differentiated_stack(verts, faces, si_thickness_mm=SI_THICKNESS_MM, si_path=si_path, si_label=f'{tag} silicon', kapton_thickness_mm=KAPTON_THICKNESS_MM, kapton_path=kapton_path, kapton_label=f'{tag} kapton')
    else:
        thicken_and_export_stack(verts, faces, [(SI_THICKNESS_MM, si_path, f'{tag} silicon'), (KAPTON_THICKNESS_MM, kapton_path, f'{tag} kapton')])

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for structure, gen_fn in SWEPT_GENERATORS.items():
        for theta_deg in THETA_DEG_VALUES:
            print_theta_stats(structure, theta_deg)
            verts, faces = gen_fn(theta_deg)
            tag = f'{structure}_theta{theta_deg:g}'
            export_pair(tag, verts, faces)
if __name__ == '__main__':
    main()