#!/usr/bin/env python3
"""Prove that the supplied validator rejects numeric and rendered-content corruption.
Only run in an unsealed review workspace; inputs are restored in finally blocks.
"""
import json
from pathlib import Path
import subprocess
import sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
if (ROOT/'PACKAGE_FILE_MANIFEST.json').exists():
    raise SystemExit('Run fault probes before sealing, so hash checks cannot mask numerical/content detection.')
results=[]
for kind,path in [('BH corruption',ROOT/'results/panel_results.csv'),('Unrendered fabricated claim',ROOT/'results/manuscript_content.json')]:
    original=path.read_bytes()
    try:
        if path.suffix=='.csv':
            d=pd.read_csv(path);d.loc[0,'q_BH_all_tests_in_family']=.000001;d.to_csv(path,index=False)
        else:
            d=json.loads(path.read_text());d.insert(5,{'type':'paragraph','text':'FAULT PROBE fabricated claim 987654321 is established.'});path.write_text(json.dumps(d))
        run=subprocess.run([sys.executable,str(ROOT/'scripts/validate_revision.py')],cwd=ROOT,capture_output=True,text=True)
        report=json.loads((ROOT/'results/validation_report.json').read_text())
        target='independent BH'if path.suffix=='.csv'else'rendered text block'
        caught=run.returncode!=0 and any(target in x for x in report['computational_errors'])
        results.append({'probe':kind,'rejected':caught,'errors':report['computational_errors']})
    finally:path.write_bytes(original)
report={'all_probes_rejected':all(x['rejected']for x in results),'probes':results,'files_restored':True}
(ROOT/'results/validation_fault_injection_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
if not report['all_probes_rejected']:raise SystemExit(1)
