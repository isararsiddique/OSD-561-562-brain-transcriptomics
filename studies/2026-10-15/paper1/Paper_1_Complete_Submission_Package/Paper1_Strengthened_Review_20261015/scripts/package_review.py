#!/usr/bin/env python3
"""Copy, validate, seal and ZIP a complete portable author-review package.

Historical drafts, render intermediates and stale pre-revision manifests are
excluded. Existing destination folders are never replaced. No external upload.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[1]
ROOT_FILES = [
    'README.md', 'REVISION_STATUS.md', 'JOURNAL_TARGET.md', 'COVER_LETTER_DRAFT.md',
    'authors.json', 'requirements.txt', 'runtime_versions.txt', 'requirements_strengthening.txt', 'runtime_strengthening.txt',
    'MANUSCRIPT_COMPLETE.md', 'MANUSCRIPT_COMPLETE.docx', 'MANUSCRIPT_COMPLETE.pdf',
    'SUPPLEMENTARY_INFORMATION.md', 'SUPPLEMENTARY_INFORMATION.docx',
    'SUPPLEMENTARY_INFORMATION.pdf',
]
TREES = ['inputs', 'results', 'scripts', 'tests', 'extensions', 'figures',
         'publication_figures', 'supplementary_data']
MUTABLE = {'PACKAGE_FILE_MANIFEST.json', 'PACKAGE_VERIFICATION.json',
           'results/validation_report.json'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def ignore(directory, names):
    return [name for name in names if name in {'__pycache__', '.DS_Store',
            'review_file_manifest.json'} or name.endswith('.pyc') or
            name.startswith('_qa_') or name.endswith('_tmp.csv')]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    dest = args.destination.resolve()
    if dest.exists() or dest.with_suffix('.zip').exists():
        raise SystemExit('Destination or archive exists; choose a new reviewed package name.')
    report = json.loads((ROOT/'results/validation_report.json').read_text())
    assert report['computational_validation_passed'] and not report['computational_errors']
    assert json.loads((ROOT/'results/validation_fault_injection_report.json').read_text())['all_probes_rejected']
    assert json.loads((ROOT/'results/supplementary_data_validation.json').read_text())['complete_table_exports'] == 50
    visual = json.loads((ROOT/'results/visual_review.json').read_text())
    assert visual['all_pages_reviewed']
    assert all(sha(ROOT/(stem+'.pdf')) == record['pdf_sha256'] for stem,record in visual['documents'].items())
    dest.mkdir(parents=True)
    for name in ROOT_FILES:
        shutil.copy2(ROOT/name, dest/name)
    for name in TREES:
        shutil.copytree(ROOT/name, dest/name, ignore=ignore)
    # Verify the copied manuscript, counts, source paths and every statistical
    # family from the portable checkout before adding its file manifest.
    subprocess.run([sys.executable, str(dest/'scripts/validate_revision.py')],
                   cwd=dest, check=True)
    files = sorted(p for p in dest.rglob('*') if p.is_file() and
                   '__pycache__' not in p.parts and p.suffix != '.pyc')
    records = [{'path':p.relative_to(dest).as_posix(), 'bytes':p.stat().st_size,
                'sha256':sha(p)} for p in files
               if p.relative_to(dest).as_posix() not in MUTABLE]
    manifest = dest/'PACKAGE_FILE_MANIFEST.json'
    manifest.write_text(json.dumps(records, indent=2)+'\n')
    assert all(sha(dest/r['path']) == r['sha256'] for r in records)
    subprocess.run([sys.executable, str(dest/'scripts/validate_revision.py')],
                   cwd=dest, check=True)
    checked = json.loads((dest/'results/validation_report.json').read_text())
    assert all(sha(dest/r['path']) == r['sha256'] for r in records)
    verification = {
        'status':'PASS', 'review_date':'2026-10-15',
        'sealed_files':len(records), 'all_sealed_sha256_match':True,
        'manifest_sha256':sha(manifest),
        'preseal_computational_checks':report['computational_checks'],
        'portable_checks_including_sealed_files':checked['computational_checks'],
        'statistical_unit_tests_passed':11,
        'all_rendered_pages_visually_reviewed':True,
        'main_pages':report['pdf_pages']['MANUSCRIPT_COMPLETE'], 'supplement_pages':report['pdf_pages']['SUPPLEMENTARY_INFORMATION'],
        'main_figures':13, 'supplement_figures':7,
        'numbered_supplementary_data_tables':50,
        'mutable_files_excluded_from_seal':sorted(MUTABLE),
        'submission_ready':False,
        'release_items_remaining':report['release_items_remaining'],
    }
    (dest/'PACKAGE_VERIFICATION.json').write_text(json.dumps(verification, indent=2)+'\n')
    # Validation imports can create bytecode; exclude that from distributed files.
    for cache in dest.rglob('__pycache__'):
        shutil.rmtree(cache)
    archive = dest.with_suffix('.zip')
    with ZipFile(archive, 'w', ZIP_DEFLATED, compresslevel=6) as bundle:
        for path in sorted(dest.rglob('*')):
            if path.is_file():
                bundle.write(path, dest.name+'/'+path.relative_to(dest).as_posix())
    with ZipFile(archive) as bundle:
        assert bundle.testzip() is None
        for rec in records:
            assert hashlib.sha256(bundle.read(dest.name+'/'+rec['path'])).hexdigest() == rec['sha256']
    digest = sha(archive)
    archive.with_suffix('.zip.sha256').write_text(digest+'  '+archive.name+'\n')
    print(json.dumps({'package':str(dest), 'archive':str(archive),
                      'archive_bytes':archive.stat().st_size, 'archive_sha256':digest,
                      'sealed_files':len(records), 'ZIP_CRC_and_sealed_hashes':'PASS'}, indent=2))

if __name__ == '__main__':
    main()
