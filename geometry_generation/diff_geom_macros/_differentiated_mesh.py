import numpy as np
from collections import defaultdict

def weld_vertices(verts, faces, tol=1e-07):
    verts = np.asarray(verts, dtype=float).reshape(-1, 3)
    n = len(verts)
    scale = 1.0 / tol
    keys = np.round(verts * scale).astype(np.int64)
    bucket = {}
    remap = np.empty(n, dtype=np.int64)
    new_verts = []
    for i in range(n):
        key = (keys[i, 0], keys[i, 1], keys[i, 2])
        if key in bucket:
            remap[i] = bucket[key]
        else:
            new_idx = len(new_verts)
            bucket[key] = new_idx
            new_verts.append(verts[i])
            remap[i] = new_idx
    new_verts = np.array(new_verts, dtype=float)
    new_faces = [(int(remap[i]), int(remap[j]), int(remap[k])) for i, j, k in faces]
    n_welded = n - len(new_verts)
    return (new_verts, new_faces, n_welded)

def compute_vertex_normals(verts, faces):
    verts = np.asarray(verts, dtype=float)
    n = len(verts)
    accum = np.zeros((n, 3), dtype=float)
    n_degenerate = 0
    for i, j, k in faces:
        p0, p1, p2 = (verts[i], verts[j], verts[k])
        fn = np.cross(p1 - p0, p2 - p0)
        norm = np.linalg.norm(fn)
        if norm < 1e-12:
            n_degenerate += 1
            continue
        accum[i] += fn
        accum[j] += fn
        accum[k] += fn
    if n_degenerate:
        print(f'[compute_vertex_normals] skipped {n_degenerate} degenerate face(s) in normal accumulation')
    lengths = np.linalg.norm(accum, axis=1)
    zero_mask = lengths < 1e-12
    if zero_mask.any():
        print(f'[compute_vertex_normals] WARNING: {zero_mask.sum()} vertex/vertices have no valid incident face area -- check mesh connectivity')
        lengths[zero_mask] = 1.0
    normals = accum / lengths[:, None]
    normals[zero_mask] = 0.0
    return normals

def compute_seam_clearance(verts, faces, vertex_normals, facet_thickness_mm, min_alignment=0.05, max_pad_multiple=3.0):
    max_pad_mm = facet_thickness_mm * max_pad_multiple
    n = len(verts)
    required = np.full(n, facet_thickness_mm, dtype=float)
    n_flagged = 0
    n_ceiling = 0
    max_pad = 0.0
    n_padded = 0
    for i, j, k in faces:
        p0, p1, p2 = (verts[i], verts[j], verts[k])
        fn = np.cross(p1 - p0, p2 - p0)
        norm = np.linalg.norm(fn)
        if norm < 1e-12:
            continue
        fn = fn / norm
        for idx in (i, j, k):
            nv = vertex_normals[idx]
            denom = float(np.dot(fn, nv))
            if denom > min_alignment:
                cand = facet_thickness_mm / denom
            else:
                cand = max_pad_mm
                n_flagged += 1
            if cand > max_pad_mm:
                cand = max_pad_mm
                n_ceiling += 1
            if cand > required[idx]:
                if required[idx] == facet_thickness_mm and cand > facet_thickness_mm:
                    n_padded += 1
                pad = cand - facet_thickness_mm
                if pad > max_pad:
                    max_pad = pad
                required[idx] = cand
    print(f'[compute_seam_clearance] {n_padded} of {n} vertex/vertices needed seam-clearance padding beyond the flat {facet_thickness_mm}mm value (max pad {max_pad:.4f} mm, ceiling {max_pad_mm:.4f} mm)')
    if n_flagged:
        print(f'[compute_seam_clearance] {n_flagged} face/vertex incidence(s) had face-normal/vertex-normal alignment <= {min_alignment} -- used the {max_pad_multiple}x ceiling there.')
    if n_ceiling:
        print(f"[compute_seam_clearance] WARNING: {n_ceiling} face/vertex incidence(s) hit the {max_pad_multiple}x clearance CEILING ({max_pad_mm:.4f} mm) rather than getting the full projected clearance -- these vertices are NOT guaranteed to clear silicon anymore. This means the seam there is too acute for a single continuous per-vertex Kapton offset to represent safely. Options: (1) accept a small residual overlap risk at these specific vertices if they're few and localized, (2) exclude this parameter point the same way you already exclude Kresling theta=60deg / Yoshimura N=6, or (3) split Kapton at these vertices too (duplicate, like silicon does) instead of forcing continuity through an angle this sharp.")
    return required

def diagnose_offset_validity(orig_verts, faces, outer, inner, required=None, flip_report_limit=20):
    verts = np.asarray(orig_verts, dtype=float)
    n_flipped = 0
    n_warn = 0
    n_checked = 0
    flipped_faces = []
    reports = []
    for f_idx, (i, j, k) in enumerate(faces):
        p0, p1, p2 = (verts[i], verts[j], verts[k])
        src_n = np.cross(p1 - p0, p2 - p0)
        src_norm = np.linalg.norm(src_n)
        if src_norm < 1e-12:
            continue
        src_n = src_n / src_norm
        o0, o1, o2 = (outer[i], outer[j], outer[k])
        out_n = np.cross(o1 - o0, o2 - o0)
        out_norm = np.linalg.norm(out_n)
        n_checked += 1
        if out_norm < 1e-12:
            severity = 'FLIP'
            alignment = float('nan')
            n_flipped += 1
            flipped_faces.append(f_idx)
        else:
            out_n = out_n / out_norm
            alignment = float(np.dot(src_n, out_n))
            if alignment < 0:
                severity = 'FLIP'
                n_flipped += 1
                flipped_faces.append(f_idx)
            elif alignment < 0.3:
                severity = 'WARN'
                n_warn += 1
            else:
                continue
        pad_note = ''
        if required is not None:
            pads = [required[i], required[j], required[k]]
            pad_note = f', padding spread {min(pads):.4f}-{max(pads):.4f}mm (delta {max(pads) - min(pads):.4f}mm)'
        reports.append((severity, f_idx, i, j, k, alignment, pad_note))
    print(f'[diagnose_offset_validity] {n_checked} faces checked, {n_flipped} inverted, {n_warn} badly distorted')
    for severity, f_idx, i, j, k, alignment, pad_note in reports[:flip_report_limit]:
        print(f'  [{severity}] face {f_idx} (verts {i},{j},{k}): alignment={alignment:.3f}{pad_note}')
    if len(reports) > flip_report_limit:
        print(f'  ... {len(reports) - flip_report_limit} more flagged faces (raise flip_report_limit to see all)')
    return (n_flipped, flipped_faces)

def offset_shell(verts, faces, z_start, z_end, tol=1e-07, clear_facet_thickness_mm=None, diagnose=True):
    verts, faces, n_welded = weld_vertices(verts, faces, tol=tol)
    if n_welded:
        print(f'[offset_shell] welded {n_welded} duplicate vertex/vertices before offsetting')
    normals = compute_vertex_normals(verts, faces)
    n = len(verts)
    z_start_arr = np.full(n, float(z_start), dtype=float)
    thickness = float(z_end) - float(z_start)
    if clear_facet_thickness_mm is not None:
        required = compute_seam_clearance(verts, faces, normals, clear_facet_thickness_mm)
        z_start_arr = np.maximum(z_start_arr, required)
    z_end_arr = z_start_arr + thickness
    outer = verts + normals * z_end_arr[:, None]
    inner = verts + normals * z_start_arr[:, None]
    if diagnose:
        required_report = z_start_arr if clear_facet_thickness_mm is not None else None
        diagnose_offset_validity(verts, faces, outer, inner, required=required_report)
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
    for i, j, k in faces:
        add_tri(outer[i], outer[j], outer[k])
        add_tri(inner[i], inner[k], inner[j])
    edge_faces = defaultdict(list)
    for f_idx, (i, j, k) in enumerate(faces):
        for a, b in [(i, j), (j, k), (k, i)]:
            edge_faces[tuple(sorted((a, b)))].append((f_idx, a, b))
    n_boundary = 0
    for edge, owners in edge_faces.items():
        if len(owners) != 1:
            continue
        n_boundary += 1
        _, a, b = owners[0]
        add_tri(inner[a], inner[b], outer[b])
        add_tri(inner[a], outer[b], outer[a])
    print(f'[offset_shell] {len(faces)} faces, {n_boundary} true boundary edges walled (vs. {3 * len(faces)} edges an all-facets wall like thicken_triangles would use)')
    return (np.array(out_verts), out_faces)

def build_differentiated_kapton(verts, faces, kapton_thickness_mm, silicon_thickness_mm, tol=1e-07):
    z_start = silicon_thickness_mm
    z_end = silicon_thickness_mm + kapton_thickness_mm
    return offset_shell(verts, faces, z_start, z_end, tol=tol, clear_facet_thickness_mm=silicon_thickness_mm)