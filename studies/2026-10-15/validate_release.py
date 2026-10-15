#!/usr/bin/env python3
"""Verify distributed analysis sources and dataset exclusion without refitting."""
import ast
import csv
import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
checked = python_sources = 0
with (HERE / 'SOURCE_FILE_MANIFEST.tsv').open() as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['distribution'] != 'included':
            continue
        path = HERE / row['study'] / row['source_relative_path']
        data = path.read_bytes()
        assert len(data) == int(row['bytes']), path
        assert hashlib.sha256(data).hexdigest() == row['sha256'], path
        if path.suffix == '.py':
            ast.parse(data, filename=str(path))
            python_sources += 1
        checked += 1
for paper in ('paper1', 'paper2'):
    for path in (HERE / paper).rglob('*'):
        assert not path.is_symlink(), path
        if path.is_file():
            relative = path.relative_to(HERE / paper)
            assert not any(part in {'inputs', 'data', 'datasets', 'cache', 'sources', 'markers', 'briefs'} for part in relative.parts), path
            assert path.suffix not in {'.npz', '.gz', '.zip', '.tab', '.gtf', '.gff', '.bam', '.fastq'}, path
print(f'PASS: {checked} archive files match; {python_sources} Python sources parse; input datasets and binary dataset formats are excluded.')
