#!/usr/bin/env python3
"""Rehash both controls and recompute their independent rank proofs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_rank import verify

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / 'evidence'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_archive(folder: str, manifest_name: str) -> dict:
    raw = EVIDENCE / folder
    manifest_path = EVIDENCE / manifest_name
    manifest = json.loads(manifest_path.read_text())
    assert manifest['schema'] == 'compact-orbit-rank-control-archive-v1'
    observed = []
    for path in sorted(raw.rglob('*')):
        if path.is_file():
            observed.append({'path': path.relative_to(raw).as_posix(),
                             'bytes': path.stat().st_size, 'sha256': sha(path)})
    assert observed == manifest['files']
    receipt = json.loads((raw / 'receipt.json').read_text())
    assert receipt['status'] == 'PASS'
    assert receipt['traced_and_control_deterministic_outputs_match'] is True
    listed = {row['path']: row['sha256'] for row in observed if row['path'] != 'receipt.json'}
    assert receipt['files'] == listed
    traced = raw / 'traced'
    command = json.loads((traced / 'command.json').read_text())
    assert command['source_head'].startswith(manifest['source_commit'])
    replay = verify(traced / 'rank.jsonl', traced / 'base.jsonl', traced / 'summary.jsonl')
    assert replay == json.loads((raw / 'independent_rank_replay.json').read_text())
    assert replay['status'] == 'PASS' and replay['rank'] == 2
    return {'raw_files': len(observed), 'manifest_sha256': sha(manifest_path),
            'trace_sha256': replay['trace_sha256'], 'base_sha256': replay['base_sha256'],
            'rank': replay['rank'], 'rank_attempts': replay['attempts']}


def main() -> None:
    local = check_archive('local_control_20260929', 'MANIFEST.json')
    hosted = check_archive('hosted_36550847279', 'HOSTED_MANIFEST.json')
    assert local['trace_sha256'] == hosted['trace_sha256']
    assert local['base_sha256'] == hosted['base_sha256']
    assert local['rank'] == hosted['rank']
    assert local['rank_attempts'] == hosted['rank_attempts']
    print(json.dumps({'status': 'PASS', 'local': local, 'hosted': hosted}, sort_keys=True))


if __name__ == '__main__':
    main()
