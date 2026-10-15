#!/usr/bin/env python3
"""Verify and restore a complete study archive locally; never commit its data."""
import argparse
import hashlib
import json
import shutil
import stat
import sys
import urllib.error
import urllib.request
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def file_hash(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(metadata, target):
    request = urllib.request.Request(
        f"https://zenodo.org/api/records/{metadata['record_id']}",
        headers={'User-Agent': 'spaceflight-brain-reproduction/2026-10-15'})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            record = json.load(response)
    except urllib.error.HTTPError as error:
        if error.code in (403, 404):
            raise ValueError('This record is not publicly downloadable yet. Use --archive with the verified local ZIP, or wait for Zenodo publication.') from error
        raise
    files = record.get('files', [])
    if isinstance(files, dict):
        entries = files.get('entries', {})
        files = list(entries.values()) if isinstance(entries, dict) else entries
    match = next((f for f in files if f.get('key') == metadata['filename']), None)
    if match is None:
        raise ValueError('Expected release ZIP not found in the public record.')
    links = match.get('links', {})
    url = links.get('content') or links.get('self')
    if not url or not url.startswith('https://zenodo.org/'):
        raise ValueError('Expected a Zenodo-hosted archive download URL.')
    partial = target.with_suffix(target.suffix + '.partial')
    with urllib.request.urlopen(url, timeout=120) as source, partial.open('wb') as output:
        shutil.copyfileobj(source, output, length=1024 * 1024)
    if file_hash(partial) != metadata['sha256']:
        raise ValueError(f'Archive checksum mismatch; untrusted download retained at {partial}')
    partial.replace(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paper', required=True, choices=('paper1', 'paper2'))
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--archive', type=Path)
    source.add_argument('--download', action='store_true')
    parser.add_argument('--check-only', action='store_true', help='Verify ZIP identity without extracting')
    args = parser.parse_args()
    metadata = json.loads((HERE / 'ARCHIVES.json').read_text())[args.paper]
    local = ROOT / '.local-reproduction'
    local.mkdir(exist_ok=True)
    archive = args.archive.resolve() if args.archive else local / metadata['filename']
    if args.download and not archive.exists():
        download(metadata, archive)
    if not archive.is_file() or file_hash(archive) != metadata['sha256']:
        raise ValueError('Missing archive or SHA-256 mismatch. Use the exact 15 October 2026 ZIP named in ARCHIVES.json.')
    print(f'{args.paper}: archive SHA-256 PASS')
    if args.check_only:
        return
    destination = local / args.paper
    destination.mkdir(exist_ok=True)
    with ZipFile(archive) as z:
        entries = z.infolist()
        prefix = entries[0].filename.split('/')[0] + '/'
        targets = []
        for entry in entries:
            path = PurePosixPath(entry.filename)
            if path.is_absolute() or '..' in path.parts or '\\' in entry.filename or not entry.filename.startswith(prefix):
                raise ValueError(f'Unsafe archive path: {entry.filename}')
            if stat.S_ISLNK(entry.external_attr >> 16):
                raise ValueError(f'Archive symlink is not allowed: {entry.filename}')
            relative = entry.filename[len(prefix):]
            if not relative or entry.is_dir():
                continue
            target = destination / relative
            if not target.resolve().is_relative_to(destination.resolve()):
                raise ValueError(f'Unsafe destination: {relative}')
            if target.exists():
                with z.open(entry) as member:
                    h = hashlib.sha256()
                    for block in iter(lambda: member.read(1024 * 1024), b''):
                        h.update(block)
                if not target.is_file() or file_hash(target) != h.hexdigest():
                    raise ValueError(f'Reproduction file has changed: {target}. Move the existing local reproduction folder aside before restoring again.')
            targets.append((entry, target))
        for entry, target in targets:
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with z.open(entry) as source, target.open('wb') as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
    print(f'Restored {len(targets)} files to {destination}')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, urllib.error.URLError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
