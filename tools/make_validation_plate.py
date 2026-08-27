import argparse
import numpy as np
import trimesh

def make_plate(side_mm: float, out_path: str) -> None:
    h = side_mm / 2.0
    vertices = np.array([[-h, -h, 0.0], [h, -h, 0.0], [h, h, 0.0], [-h, h, 0.0]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    mesh.export(out_path)
    print(f"[make_validation_plate] wrote '{out_path}', {side_mm}mm square, normal along +Z, Z-extent=0 before thickening.")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--side', type=float, default=50.0, help='Plate side length in mm')
    p.add_argument('--out', default='plate_shell.stl')
    args = p.parse_args()
    make_plate(args.side, args.out)
if __name__ == '__main__':
    main()