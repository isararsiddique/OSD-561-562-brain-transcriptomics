#!/usr/bin/env python3
"""Run the exact analysis and validation sequentially with frozen input data."""
import argparse
from pathlib import Path
import os
import subprocess
import sys
import analyse_gene_omnibus as a

parser=argparse.ArgumentParser()
parser.add_argument('--base',type=Path,default=a.DEFAULT_BASE)
args=parser.parse_args()
env=os.environ.copy()
env.setdefault('OPENBLAS_NUM_THREADS','4')
env.setdefault('OMP_NUM_THREADS','4')
for script in ['analyse_gene_omnibus.py','verify_gene_omnibus.py','sensitivity_dominance95.py',
               'descriptive_leave_one_out.py','verify_sensitivities.py','build_figure_report.py']:
    subprocess.run([sys.executable,str(a.HERE/'scripts'/script),'--base',str(args.base.resolve())],
                   cwd=a.HERE,env=env,check=True)
