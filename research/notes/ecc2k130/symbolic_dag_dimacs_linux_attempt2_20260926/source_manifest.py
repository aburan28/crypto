#!/usr/bin/env python3
"""Check official CaDiCaL 3.0.1 tracked Git blobs without building or solving."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent


def git(source: Path, *args: str) -> bytes:
    return subprocess.check_output(['git', *args], cwd=source)


def scan(source: Path) -> dict:
    indexed = [entry for entry in git(source, 'ls-files', '-s', '-z').split(b'\0')
               if entry]
    entries = []
    for entry in indexed:
        metadata, name = entry.split(b'\t', 1)
        mode, oid, stage = metadata.split()
        if stage != b'0':
            raise ValueError('unmerged source index')
        entries.append((mode, oid, name))
    batch = subprocess.run(['git', 'cat-file', '--batch'], cwd=source,
                           input=b''.join(oid + b'\n' for _, oid, _ in entries),
                           capture_output=True, check=True).stdout
    stream = io.BytesIO(batch)
    digest = hashlib.sha256(b'cadical-source-manifest-v2\0')
    for mode, oid, name in entries:
        header = stream.readline().split()
        if len(header) != 3 or header[0] != oid or header[1] != b'blob':
            raise ValueError('source Git blob response mismatch')
        raw = stream.read(int(header[2]))
        if stream.read(1) != b'\n':
            raise ValueError('source Git blob framing mismatch')
        digest.update(len(name).to_bytes(4, 'big'))
        digest.update(name)
        digest.update(mode)
        digest.update(hashlib.sha256(raw).digest())
    if stream.read(1):
        raise ValueError('unexpected trailing source Git blob data')
    return {
        'commit': git(source, 'rev-parse', 'HEAD').decode().strip(),
        'tree': git(source, 'rev-parse', 'HEAD^{tree}').decode().strip(),
        'tracked_files': len(entries),
        'source_manifest_sha256': digest.hexdigest(),
        'version': (source / 'VERSION').read_text().strip(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    args = parser.parse_args()
    expected = json.loads((HERE / 'PREPARE.json').read_text())['cadical']
    got = scan(args.source.resolve())
    for key in ('commit', 'tree', 'tracked_files', 'source_manifest_sha256'):
        if got[key] != expected[key]:
            raise SystemExit(f'NOT_ADMITTED: CaDiCaL {key} mismatch')
    if got['version'] != '3.0.1':
        raise SystemExit('NOT_ADMITTED: CaDiCaL version mismatch')
    print(json.dumps({'decision': 'PASS', **got}, sort_keys=True))


if __name__ == '__main__':
    main()
