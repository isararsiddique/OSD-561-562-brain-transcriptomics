"""Acquire the recorded study inputs from primary repositories with SHA256 checks."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import time
import urllib.request
import zlib

REPOSITORY = Path(__file__).resolve().parents[1]


def within(root, name):
    path = (root / name).resolve()
    if Path(name).is_absolute() or not path.is_relative_to(root) or path == root:
        raise ValueError(f"Unsafe manifest path: {name}")
    return path


def valid(path, record):
    if not path.is_file() or path.stat().st_size != record['bytes']:
        return False
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest() == record['sha256']


def fetch(url, expected_bytes, target=None, byte_range=None):
    if not url.startswith('https://'):
        raise ValueError('Only recorded HTTPS source URLs are allowed')
    headers = {'User-Agent': 'OSDR-CNS-reproducibility/1.0', 'Accept-Encoding': 'identity'}
    if byte_range is not None:
        headers['Range'] = f'bytes={byte_range[0]}-{byte_range[1]}'
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=120) as response:
                if byte_range is not None:
                    content_range = response.headers.get('Content-Range', '')
                    wanted = f'bytes {byte_range[0]}-{byte_range[1]}/'
                    if response.status != 206 or not content_range.startswith(wanted):
                        raise ValueError('Source did not honor the exact ZIP byte range; no full archive was downloaded')
                size = 0
                data = bytearray()
                handle = target.open('wb') if target is not None else None
                try:
                    while chunk := response.read(1024 * 1024):
                        size += len(chunk)
                        if size > expected_bytes:
                            raise ValueError('Source exceeds the recorded length')
                        if handle is not None: handle.write(chunk)
                        else: data.extend(chunk)
                finally:
                    if handle is not None: handle.close()
                if size != expected_bytes:
                    raise ValueError(f'Source length differs: expected {expected_bytes}, received {size}')
                return bytes(data)
        except (OSError, ValueError):
            if attempt == 2: raise
            time.sleep(2 ** attempt)


def zip_bytes(record):
    if record['method'] == 'deflate_range':
        start, end = record['range']
        compressed = fetch(record['source_url'], end-start+1, byte_range=(start,end))
        raw = zlib.decompress(compressed, -15)
        crc = record['crc32']
    else:
        member = record['zip']
        offset = member['offset']
        header = fetch(record['source_url'], 30, byte_range=(offset, offset+29))
        signature, version, flags, compression, _, _, _, _, _, nlen, elen = struct.unpack('<IHHHHHIIIHH', header)
        if signature != 0x04034b50 or flags & 1 or compression != member['compression']:
            raise ValueError('ZIP local header disagrees with the recorded inventory')
        name = fetch(record['source_url'], nlen, byte_range=(offset+30,offset+29+nlen)).decode('utf-8')
        if name != record['archive_member']:
            raise ValueError('ZIP member name disagrees with the recorded inventory')
        start = offset + 30 + nlen + elen
        compressed = fetch(record['source_url'], member['compressed_size'], byte_range=(start,start+member['compressed_size']-1))
        raw = zlib.decompress(compressed,-15) if compression == 8 else compressed if compression == 0 else None
        if raw is None: raise ValueError(f'Unsupported ZIP compression: {compression}')
        crc = member['CRC']
    if len(raw) != record['bytes'] or zlib.crc32(raw) & 0xffffffff != crc:
        raise ValueError('ZIP member length or CRC32 differs from the recorded source')
    return raw


def acquire(destination, cache=None, offline=False, only=None):
    root = Path(destination).resolve()
    cache = Path(cache).resolve() if cache else None
    if root == REPOSITORY or REPOSITORY.is_relative_to(root):
        raise ValueError('The destination cannot contain or replace the code checkout')
    manifest = json.loads((REPOSITORY/'data/download_manifest.json').read_text())
    definition = json.loads((REPOSITORY/'data/input_manifest.json').read_text())
    if manifest.get('schema_version') != 1 or definition.get('schema_version') != 1:
        raise ValueError('Unsupported acquisition manifest')
    records = manifest['files']
    expected = {r['path']:(r['bytes'],r['sha256']) for r in definition['files']}
    if len(expected) != len(records) or len({r['path'] for r in records}) != len(records):
        raise ValueError('Duplicate or inconsistent source records')
    for record in records:
        within(root, record['path'])
        if expected.get(record['path']) != (record['bytes'], record['sha256']):
            raise ValueError('Acquisition and input manifests disagree')
    if only:
        selected = [r for r in records if any(x in r['path'] for x in only)]
        if not selected: raise ValueError('No inputs match --only')
    else: selected = records
    root.mkdir(parents=True,exist_ok=True)
    lock = root/'.acquisition.lock'
    fd = os.open(lock, os.O_CREAT|os.O_EXCL|os.O_WRONLY, 0o600)
    os.close(fd)
    try:
        for i, record in enumerate(selected,1):
            target = within(root,record['path'])
            if target.exists():
                if not valid(target,record): raise ValueError(f'Existing input differs: {target}; preserve it and use a fresh destination')
                print(f'[{i}/{len(selected)}] verified {record["path"]}',flush=True)
                continue
            target.parent.mkdir(parents=True,exist_ok=True)
            temp = target.with_name(target.name+'.partial')
            try:
                cached = within(cache,record['path']) if cache else None
                method = record['method']
                if cached and valid(cached,record): shutil.copyfile(cached,temp)
                elif method in ('repository_file','repository_gzip'):
                    source = within(REPOSITORY,record['repository_path'])
                    if method == 'repository_gzip':
                        with gzip.open(source,'rb') as f,temp.open('wb') as out: shutil.copyfileobj(f,out)
                    else: shutil.copyfile(source,temp)
                elif method == 'copy' and valid(within(root,record['copy_from']),record):
                    shutil.copyfile(within(root,record['copy_from']),temp)
                elif offline: raise ValueError(f'Offline input unavailable: {record["path"]}')
                elif method in ('zip_member','deflate_range'): temp.write_bytes(zip_bytes(record))
                elif method in ('download','copy'): fetch(record['source_url'],record['bytes'],target=temp)
                else: raise ValueError(f'Unknown source method: {method}')
                if not valid(temp,record): raise ValueError(f'SHA256 mismatch: {record["path"]}. The source version changed; no input was accepted')
                temp.replace(target)
            finally:
                if temp.exists(): temp.unlink()
            print(f'[{i}/{len(selected)}] acquired {record["path"]}',flush=True)
        if all(valid(within(root,r['path']),r) for r in records):
            (root/'input_manifest.json').write_text(json.dumps(definition,indent=2)+'\n')
            print(f'Complete input bundle: {len(records)} verified files in {root}',flush=True)
        else: print('Partial acquisition; no complete input manifest was written.',flush=True)
    finally:
        lock.unlink()


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination',type=Path,default=Path('pipeline_inputs'))
    parser.add_argument('--cache',type=Path,help='Reuse a local bundle only when each file matches the published SHA256')
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--only',action='append',help='Acquire matching paths only; repeated filters are combined')
    args=parser.parse_args(argv)
    acquire(args.destination,args.cache,args.offline,args.only)


if __name__ == '__main__': main()
