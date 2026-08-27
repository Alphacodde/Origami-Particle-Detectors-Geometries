import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from diff_geom_macros import comparison_barrel, kresling, miura_cylindrical, yoshimura, _geometry_config
GEOMETRIES = {'nominal_300um': 0.3, 'thin_50um': 0.05}

def run():
    root = Path(__file__).parent.parent / 'geometries'
    for folder, si_thick in GEOMETRIES.items():
        out = root / folder
        out.mkdir(parents=True, exist_ok=True)
        print(f'Generating {folder} (si={si_thick * 1000.0:.0f} um) -> {out}')
        _geometry_config.SI_THICKNESS_MM = si_thick
        for mod, name in [(comparison_barrel, 'barrel_reference'), (kresling, 'kresling_deployed'), (miura_cylindrical, 'miura_deployed'), (yoshimura, 'yoshimura_deployed')]:
            print(f'  {name} ...', end=' ', flush=True)
            mod.export_gdml(str(out / f'{name}_silicon.gdml'), str(out / f'{name}_kapton.gdml'))
            print('done')
    print('\nAll geometries written.')
if __name__ == '__main__':
    run()