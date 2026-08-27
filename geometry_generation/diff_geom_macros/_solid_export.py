import numpy as np
try:
    from _differentiated_mesh import weld_vertices as _weld_vertices
except ImportError:
    _weld_vertices = None

def thicken_triangles(verts, faces, z_start, z_end, weld_tol=1e-07):
    verts = np.asarray(verts, dtype=float)
    out_verts = []
    out_faces = []

    def add_tri(a, b, c):
        i0 = len(out_verts)
        out_verts.append(a)
        i1 = len(out_verts)
        out_verts.append(b)
        i2 = len(out_verts)
        out_verts.append(c)
        out_faces.append((i0, i1, i2))
    n_degenerate = 0
    for i, j, k in faces:
        p0, p1, p2 = (verts[i], verts[j], verts[k])
        n = np.cross(p1 - p0, p2 - p0)
        norm = np.linalg.norm(n)
        if norm < 1e-12:
            n_degenerate += 1
            continue
        n = n / norm
        t0, t1, t2 = (p0 + n * z_end, p1 + n * z_end, p2 + n * z_end)
        b0, b1, b2 = (p0 + n * z_start, p1 + n * z_start, p2 + n * z_start)
        add_tri(t0, t1, t2)
        add_tri(b0, b2, b1)
        add_tri(b0, b1, t1)
        add_tri(b0, t1, t0)
        add_tri(b1, b2, t2)
        add_tri(b1, t2, t1)
        add_tri(b2, b0, t0)
        add_tri(b2, t0, t2)
    if n_degenerate:
        print(f'[thicken_triangles] skipped {n_degenerate} degenerate (zero-area) input triangle(s)')
    out_verts = np.array(out_verts)
    if _weld_vertices is not None:
        out_verts, out_faces, n_welded = _weld_vertices(out_verts, out_faces, tol=weld_tol)
        if n_welded:
            print(f'[thicken_triangles] welded {n_welded} coincident vertex/vertices (mitered prism corners) before returning')
    return (out_verts, out_faces)

def thicken_triangles_union(verts, faces, z_start, z_end):
    import trimesh
    verts = np.asarray(verts, dtype=float)
    prisms = []
    n_degenerate = 0
    prism_faces = [(0, 1, 2), (3, 5, 4), (3, 4, 1), (3, 1, 0), (4, 5, 2), (4, 2, 1), (5, 3, 0), (5, 0, 2)]
    for i, j, k in faces:
        p0, p1, p2 = (verts[i], verts[j], verts[k])
        n = np.cross(p1 - p0, p2 - p0)
        norm = np.linalg.norm(n)
        if norm < 1e-12:
            n_degenerate += 1
            continue
        n = n / norm
        t0, t1, t2 = (p0 + n * z_end, p1 + n * z_end, p2 + n * z_end)
        b0, b1, b2 = (p0 + n * z_start, p1 + n * z_start, p2 + n * z_start)
        prism_verts = np.array([t0, t1, t2, b0, b1, b2])
        prisms.append(trimesh.Trimesh(vertices=prism_verts, faces=prism_faces, process=False))
    if n_degenerate:
        print(f'[thicken_triangles_union] skipped {n_degenerate} degenerate (zero-area) input triangle(s)')
    if not prisms:
        raise ValueError('thicken_triangles_union: no valid (non-degenerate) input triangles')
    unioned = trimesh.boolean.union(prisms)
    unioned.fix_normals()
    print(f'[thicken_triangles_union] unioned {len(prisms)} prisms -> {len(unioned.faces)} faces, watertight={unioned.is_watertight}, winding_consistent={unioned.is_winding_consistent}')
    if not unioned.is_watertight:
        print('  WARNING: union result is NOT watertight - this is unexpected (every individual prism was verified watertight going in). Inspect this specific geometry/parameter point before trusting it in Geant4 - do not assume this is benign the way the old per-facet-list non-watertight note was.')
    return unioned

def export_solid(verts_or_mesh, faces=None, path=None, label=''):
    try:
        import trimesh
    except ImportError:
        raise SystemExit('trimesh is required: pip install trimesh')
    from _gdml_export import export_gdml_tessellated
    import os
    if isinstance(verts_or_mesh, trimesh.Trimesh):
        if faces is not None:
            raise TypeError('export_solid: pass faces=None when verts_or_mesh is already a trimesh.Trimesh')
        mesh = verts_or_mesh
    else:
        mesh = trimesh.Trimesh(vertices=verts_or_mesh, faces=faces, process=True)
        mesh.fix_normals()
    root, ext = os.path.splitext(path)
    if ext.lower() != '.gdml':
        print(f"[export_solid] NOTE: '{path}' does not end in .gdml - writing GDML content to '{root}.gdml' instead (STL is no longer written - see this function's own docstring, Option A GDML fix).")
        path = root + '.gdml'
    solid_name = os.path.basename(root)
    gdml_material = 'G4_KAPTON' if 'kapton' in label.lower() else 'G4_Si'
    print(f'[{label}] pre-export in-memory watertight={mesh.is_watertight} winding_consistent={mesh.is_winding_consistent}')
    if not mesh.is_watertight:
        print(f"  WARNING: '{label}' is NOT watertight going into GDML export. Unlike the old STL path, GDML preserves whatever connectivity is here exactly (see module docstring) - it will NOT silently heal a bad mesh on export the way nothing ever healed STL either, but you should not expect GDML to fix a genuinely broken input mesh. Investigate the mesh itself before trusting this file in Geant4.")
    export_gdml_tessellated(mesh, path, solid_name, material=gdml_material)
    print(f'[OK] {label} -> {path}  bbox_extent={mesh.extents}')
    return mesh
export_trimesh_ascii = export_solid

def thicken_and_export(verts, faces, thickness, path, label='', z_start=None, union=True):
    if z_start is None:
        z_start = -thickness / 2.0
    if union:
        mesh = thicken_triangles_union(verts, faces, z_start, z_start + thickness)
        return export_trimesh_ascii(mesh, None, path, label=f'{label} ({thickness}mm)')
    else:
        tv, tf = thicken_triangles(verts, faces, z_start, z_start + thickness)
        return export_trimesh_ascii(tv, tf, path, label=f'{label} ({thickness}mm)')

def thicken_and_export_stack(verts, faces, layers, direction=1.0, union=True):
    z = 0.0
    results = []
    for thickness, path, label in layers:
        z_start = z
        z_end = z + direction * thickness
        lo, hi = (z_start, z_end) if z_start <= z_end else (z_end, z_start)
        if union:
            mesh = thicken_triangles_union(verts, faces, lo, hi)
            mesh = export_trimesh_ascii(mesh, None, path, label=f'{label} ({thickness}mm)')
        else:
            tv, tf = thicken_triangles(verts, faces, lo, hi)
            mesh = export_trimesh_ascii(tv, tf, path, label=f'{label} ({thickness}mm)')
        results.append(mesh)
        z = z_end
    return results