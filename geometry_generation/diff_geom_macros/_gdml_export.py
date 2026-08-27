import numpy as np
_GDML_HEADER = '<?xml version="1.0" encoding="UTF-8"?>\n<gdml xmlns:gdml="http://cern.ch/2001/Schemas/GDML"\n      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n      xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">\n'

def export_gdml_tessellated(mesh, path, solid_name, material='G4_Si', length_unit='mm'):
    verts = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)
    lines = [_GDML_HEADER, '  <define>']
    for vi, (x, y, z) in enumerate(verts):
        lines.append(f'    <position name="{solid_name}_v{vi}" x="{x:.9g}" y="{y:.9g}" z="{z:.9g}" unit="{length_unit}"/>')
    lines.append('  </define>')
    lines.append('  <solids>')
    lines.append(f'    <tessellated name="{solid_name}_solid">')
    for a, b, c in faces:
        lines.append(f'      <triangular vertex1="{solid_name}_v{a}" vertex2="{solid_name}_v{b}" vertex3="{solid_name}_v{c}"/>')
    lines.append('    </tessellated>')
    lines.append('  </solids>')
    lines.append('  <structure>')
    lines.append(f'    <volume name="{solid_name}_volume">')
    lines.append(f'      <materialref ref="{material}"/>')
    lines.append(f'      <solidref ref="{solid_name}_solid"/>')
    lines.append('    </volume>')
    lines.append('  </structure>')
    lines.append(f'  <setup name="Default" version="1.0">')
    lines.append(f'    <world ref="{solid_name}_volume"/>')
    lines.append('  </setup>')
    lines.append('</gdml>')
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f'[OK] GDML (indexed tessellated solid) -> {path}  verts={len(verts)}  faces={len(faces)}  watertight(pre-export, in-memory)={mesh.is_watertight}')
    return path

def _self_test_roundtrip(gdml_path):
    import re
    import trimesh
    with open(gdml_path) as fh:
        text = fh.read()
    name_to_idx = {}
    verts = []
    for m in re.finditer('<position name="([^"]+)" x="([^"]+)" y="([^"]+)" z="([^"]+)"', text):
        name, x, y, z = m.groups()
        name_to_idx[name] = len(verts)
        verts.append([float(x), float(y), float(z)])
    faces = []
    for m in re.finditer('<triangular vertex1="([^"]+)" vertex2="([^"]+)" vertex3="([^"]+)"', text):
        a, b, c = m.groups()
        faces.append([name_to_idx[a], name_to_idx[b], name_to_idx[c]])
    mesh = trimesh.Trimesh(vertices=np.array(verts), faces=np.array(faces), process=False)
    print(f'[self-test] re-parsed {gdml_path}: verts={len(verts)} faces={len(faces)} watertight={mesh.is_watertight} winding_consistent={mesh.is_winding_consistent}')
    return mesh
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        _self_test_roundtrip(sys.argv[1])
    else:
        print('Usage: python _gdml_export.py <path-to-written.gdml>  (re-parses and checks watertightness, as a standalone diagnostic)')