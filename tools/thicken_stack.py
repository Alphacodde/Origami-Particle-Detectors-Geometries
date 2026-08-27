import argparse
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'diff_geom_macros'))

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input', help='Input shell STL (zero-thickness surface)')
    p.add_argument('--layer', nargs=3, action='append', required=True, metavar=('THICKNESS_MM', 'OUT_PATH', 'LABEL'), help='One layer: thickness (mm), output path, label. Repeat for multiple layers, in stacking order.')
    args = p.parse_args()
    try:
        import trimesh
    except ImportError:
        sys.exit('trimesh is required: pip install trimesh')
    try:
        from _solid_export import thicken_and_export_stack
    except ImportError:
        sys.exit('Could not import _solid_export.py - expected at diff_geom_macros/_solid_export.py (sibling of tools/).')
    mesh = trimesh.load(args.input, force='mesh')
    layers = [(float(t), path, label) for t, path, label in args.layer]
    thicken_and_export_stack(mesh.vertices, mesh.faces, layers)
if __name__ == '__main__':
    main()