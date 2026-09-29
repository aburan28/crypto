#!/usr/bin/env python3
"""Rehash the write-once local control and recompute its independent rank proof."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_rank import verify

HERE = Path(__file__).resolve().parent
RAW = HERE / 'evidence/local_control_20260929'
MANIFEST = HERE / 'evidence/MANIFEST.json'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    assert manifest['schema'] == 'compact-orbit-rank-control-archive-v1'
    observed = []
    for path in sorted(RAW.rglob('*')):
        if path.is_file():
            observed.append({'path': path.relative_to(RAW).as_posix(),
                             'bytes': path.stat().st_size, 'sha256': sha(path)})
    assert observed == manifest['files']
    receipt = json.loads((RAW / 'receipt.json').read_text())
    assert receipt['status'] == 'PASS'
    assert receipt['traced_and_control_deterministic_outputs_match'] is True
    listed = {row['path']: row['sha256'] for row in observed if row['path'] != 'receipt.json'}
    assert receipt['files'] == listed
    traced = RAW / 'traced'
    replay = verify(traced / 'rank.jsonl', traced / 'base.jsonl', traced / 'summary.jsonl')
    assert replay == json.loads((RAW / 'independent_rank_replay.json').read_text())
    assert replay['status'] == 'PASS' and replay['rank'] == 2
    print(json.dumps({'status': 'PASS', 'raw_files': len(observed),
                      'manifest_sha256': sha(MANIFEST), 'rank': replay['rank'],
                      'rank_attempts': replay['attempts']}, sort_keys=True))


if __name__ == '__main__':
    main()
