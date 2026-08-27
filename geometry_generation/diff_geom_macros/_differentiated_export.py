import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _solid_export import thicken_triangles, thicken_triangles_union, export_trimesh_ascii
from _differentiated_mesh import build_differentiated_kapton, weld_vertices

def thicken_and_export_differentiated_stack(verts, faces, si_thickness_mm, si_path, si_label, kapton_thickness_mm, kapton_path, kapton_label, weld_tol=1e-07):
    si_mesh_raw = thicken_triangles_union(verts, faces, 0.0, si_thickness_mm)
    si_mesh = export_trimesh_ascii(si_mesh_raw, None, si_path, label=f'{si_label} ({si_thickness_mm}mm, per-facet, unioned)')
    kap_v, kap_f = build_differentiated_kapton(verts, faces, kapton_thickness_mm=kapton_thickness_mm, silicon_thickness_mm=si_thickness_mm, tol=weld_tol)
    kapton_mesh = export_trimesh_ascii(kap_v, kap_f, kapton_path, label=f'{kapton_label} ({kapton_thickness_mm}mm, per-vertex/differentiated)')
    return (si_mesh, kapton_mesh)

def weld_report(verts, faces, tol=1e-07):
    wv, wf, n_welded = weld_vertices(verts, faces, tol=tol)
    return (len(verts), len(wv), n_welded)