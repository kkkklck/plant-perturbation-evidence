"""Reproduce figure tables or the existing main-figure renderers from derived inputs."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
from materialize_runtime import materialize
from verify_repository import verify

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step', choices=['verify','tables','figures','extended-data'])
    args=parser.parse_args()
    verify()
    if args.step == 'verify': return
    materialize()
    if args.step == 'tables':
        subprocess.run([sys.executable,str(ROOT/'scripts/159_build_paper_figure_tables_v2.py')],cwd=ROOT,check=True)
        for figure in range(1,5):
            src = next((ROOT/f'source_data/figure{figure}').glob('fig*.tsv'))
            actual = ROOT/'data/paper_figures_v2'/src.name
            if actual.read_bytes() != src.read_bytes():
                raise ValueError('Rebuilt figure table differs: '+src.name)
        print('PASS: all four rebuilt source tables are byte-identical to archived tables.')
        return
    env=os.environ.copy()
    env['MPLBACKEND']='Agg'
    env['PYTHONUTF8']='1'
    if args.step == 'figures':
        programs=[ROOT/f'scripts/figures/legacy/figures_v{v}_preview/build_v{v}.py' for v in [15,11,16]]
    else:
        programs=[ROOT/'scripts/figures'/name for name in ['166_plot_robertson_expression_filter_sensitivity.py','167_plot_heat_generic_stress_sensitivity.py']]
    for script in programs:
        subprocess.run([sys.executable,str(script)],cwd=ROOT,env=env,check=True)
    if args.step == 'figures':
        import shutil
        output=ROOT/'outputs/main_figures';output.mkdir(parents=True,exist_ok=True)
        for number, version in [(1,15),(2,11),(3,11),(4,16)]:
            for ext in ['png','pdf','svg']:
                src=ROOT/f'scripts/figures/legacy/figures_v{version}_preview/Figure{number}_v{version}_preview.{ext}'
                shutil.copyfile(src,output/f'Figure{number}.{ext}')
        print('Main figures exported to outputs/main_figures.')
    verify()


if __name__ == '__main__':
    main()
