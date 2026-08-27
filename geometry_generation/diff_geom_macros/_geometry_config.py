import os
import numpy as np
OUTPUT_DIR = os.environ.get('ORIGAMIDET_GEOM_OUTPUT_DIR', 'geometry')
os.makedirs(OUTPUT_DIR, exist_ok=True)
SHARED_R_MM = 40.0
SHARED_L_MM = 120.0
SI_THICKNESS_MM = float(os.environ.get('ORIGAMIDET_SI_THICKNESS_MM', '0.300'))
KAPTON_THICKNESS_MM = 0.05
UNIT_CELL_TARGET_MM = 20.0

def derive_repeat_count(total_length_mm: float, unit_cell_target_mm: float, minimum: int=1) -> int:
    n = round(total_length_mm / unit_cell_target_mm)
    return max(n, minimum)

def center_verts_z(verts: np.ndarray) -> np.ndarray:
    verts = np.asarray(verts, dtype=float)
    zmin = verts[:, 2].min()
    zmax = verts[:, 2].max()
    z_mid = 0.5 * (zmin + zmax)
    out = verts.copy()
    out[:, 2] -= z_mid
    return out