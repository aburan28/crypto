#!/usr/bin/env python3
"""Hash-only preregistration check, then independent committed evidence replay."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from run import HERE, FILES, freeze, sha, tree_hashes
from verify import Field


def deterministic(value):
    if isinstance(value, dict):
        return {key: deterministic(item) for key, item in value.items()
                if key not in {'wall_seconds', 'cpu_seconds', 'peak_rss_bytes'}}
    if isinstance(value, list):
        return [deterministic(item) for item in value]
    return value


def self_test() -> None:
    f = Field(3, 0b1011)
    for value in range(1, 8):
        inverse = f.inverse(value)
        assert f.mul(value, inverse) == 1
        assert f.square(value) == f.mul(value, value)
    assert f.trace_mask() == 0b001


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path)
    args = parser.parse_args()
    frozen = freeze()
    for name in ('producer.py', 'verify.py', 'run.py', 'ci_replay.py'):
        path = HERE / name
        compile(path.read_bytes(), str(path), 'exec')
    self_test()
    if args.evidence is None:
        print('Frozen unequal-arity protocol/source/input, merged parent hashes, syntax and small-field arithmetic PASS; no unequal outcome read.')
        return
    receipt = json.loads((args.evidence / 'receipt.json').read_text())
    data = json.loads((HERE / 'INPUT.json').read_text())
    assert receipt['status'] == 'success'
    assert receipt['freeze_sha256'] == sha(HERE / 'FROZEN.json')
    assert receipt['input_sha256'] == {label: frozen[label + '_sha256'] for label in FILES}
    assert receipt['parent_balanced_frozen_sha256'] == data['parent_balanced_frozen_sha256']
    assert receipt['parent_balanced_receipt_sha256'] == data['parent_balanced_receipt_sha256']
    expected_commands = [f'm{m}-{stage}' for m in (7, 8, 9, 10) for stage in ('producer', 'verify')]
    assert [row['name'] for row in receipt['commands']] == expected_commands
    assert all(row['exit_code'] == 0 and row['peak_child_rss_bytes_upper'] <= 256 * 1024 * 1024
               for row in receipt['commands'])
    raw = args.evidence / 'raw'
    assert tree_hashes(raw) == receipt['raw_sha256']
    assert {path.name: sha(path) for path in args.evidence.iterdir() if path.is_file()
            and path.name != 'receipt.json'} == receipt['stream_sha256']
    with tempfile.TemporaryDirectory() as directory:
        for m in (7, 8, 9, 10):
            fresh = Path(directory) / f'm{m}-fresh.json'
            subprocess.run([sys.executable, str(HERE / 'verify.py'), '--m', str(m),
                            '--source', str(raw / f'm{m}/producer'), '--out', str(fresh)],
                           check=True, timeout=1230)
            archived = json.loads((raw / f'm{m}/verify.json').read_text())
            assert deterministic(json.loads(fresh.read_text())) == deterministic(archived)
    print('All four unequal n131 low/high F0 mask/column censuses independently replayed.')


if __name__ == '__main__':
    main()
