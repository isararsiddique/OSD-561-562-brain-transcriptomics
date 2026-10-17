"""Check published artifact hashes, local document links and excluded assay inputs."""
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]


def main():
    manifest=json.loads((ROOT/'results/publication_manifest.json').read_text())
    for row in manifest['files']:
        p=(ROOT/row['path']).resolve()
        assert p.is_relative_to(ROOT),row['path']
        assert p.stat().st_size==row['bytes'],row['path']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
    code=json.loads((ROOT/'pipeline/code_manifest.json').read_text())
    for row in code['scripts']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    links=0
    for p in ROOT.rglob('*.md'):
        if '.git' in p.parts: continue
        for m in re.finditer(r'\]\(([^)]+)\)',p.read_text()):
            target=m.group(1).strip().strip('<>')
            if target.startswith(('https:','http:','mailto:','#')):continue
            target=unquote(target.split('#')[0])
            assert (p.parent/target).exists(),f'{p.relative_to(ROOT)}: {target}'
            links+=1
    inputs=json.loads((ROOT/'data/input_manifest.json').read_text())
    sources=json.loads((ROOT/'data/download_manifest.json').read_text())
    assert len(inputs['files'])==len(sources['files'])==303
    wanted={r['path']:(r['bytes'],r['sha256']) for r in inputs['files']}
    assert len(wanted)==303
    for r in sources['files']:
        assert wanted[r['path']]==(r['bytes'],r['sha256'])
        assert not Path(r['path']).is_absolute() and '..' not in Path(r['path']).parts
        if r['method'].startswith('repository_'):
            assert (ROOT/r['repository_path']).is_file(),r['path']
            assert not re.search(r'(STAR_|RSEM_|\.mtx|\.h5|\.RCC)',r['repository_path']),r['path']
        else: assert r['source_url'].startswith('https://'),r['path']
    for folder in ('pipeline_inputs','pipeline_results','.venv','.renv'):
        assert not (ROOT/folder).exists(),f'Excluded folder present: {folder}'
    for p in ROOT.rglob('*'):
        if not p.is_file() or '.git' in p.parts:continue
        assert p.stat().st_size < 100*1024*1024,f'GitHub file-size limit: {p}'
        assert p.suffix.lower() not in {'.h5','.mtx','.fastq','.bam','.rds','.pkl'},f'Assay/model object: {p}'
        if p.suffix in {'.py','.R','.md','.json','.csv','.tsv','.txt','.svg','.yml','.log'}:
            text=p.read_text()
            assert not re.search(r'/(?:Users|private/tmp|tmp/OSDR)(?:/|$)',text),f'Workstation path: {p}'
    print(json.dumps({'published_artifact_hashes':len(manifest['files']),
                      'script_hashes':len(code['scripts']),'local_document_links':links,
                      'input_source_records':len(sources['files']),
                      'raw_assay_inputs_excluded':True,'passed':True},indent=2))


if __name__=='__main__':main()
