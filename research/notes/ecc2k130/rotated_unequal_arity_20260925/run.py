#!/usr/bin/env python3
"""Fail-preserving frozen four-arm unequal n131 F0 census runner."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
PARENT_GATE = NOTES / 'rotated_subspace_support_20260925/gate.py'
PARENT_BALANCED = NOTES / 'rotated_arity_screen_20260925'
FILES = {'input': HERE / 'INPUT.json', 'protocol': HERE / 'PROTOCOL.md',
         'parent_gate': PARENT_GATE,
         **{name.removesuffix('.py'): HERE / name for name in
            ('producer.py', 'verify.py', 'run.py', 'ci_replay.py')}}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def peak_child_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return value if sys.platform == 'darwin' else value * 1024


def freeze() -> dict:
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    for label, path in FILES.items():
        assert sha(path) == frozen[label + '_sha256'], label
    data = json.loads((HERE / 'INPUT.json').read_text())
    assert data['domain'] == 'ECC2K130-ROTATED-UNEQUAL-ARITY-20260925-v1'
    assert data['arms'] == [
        {'m': 7, 'd_low': 18, 'd_high': 19, 'r': 5},
        {'m': 8, 'd_low': 16, 'd_high': 17, 'r': 3},
        {'m': 9, 'd_low': 14, 'd_high': 15, 'r': 5},
        {'m': 10, 'd_low': 13, 'd_high': 14, 'r': 1}]
    assert data['field_degree'] == 131 and data['beta'] == 3
    assert data['field_poly'] == (1 << 131) | (1 << 13) | 7
    assert data['group_order'] == 4 * data['q']
    assert data['caps'] == {'producer_wall_seconds': 600,
                            'verifier_wall_seconds': 1200,
                            'rss_bytes': 256 * 1024 * 1024}
    assert data['threshold'] == {'max_nonzero_signed_columns': 16384,
                                 'minimum_ordered_tuple_count_over_q_den': 1,
                                 'minimum_ordered_tuple_count_over_q_num': 1}
    assert sha(PARENT_BALANCED / 'FROZEN.json') == data['parent_balanced_frozen_sha256']
    assert sha(PARENT_BALANCED / 'evidence/receipt.json') == data['parent_balanced_receipt_sha256']
    parent_receipt = json.loads((PARENT_BALANCED / 'evidence/receipt.json').read_text())
    assert parent_receipt['status'] == 'success'
    reference = PARENT_BALANCED / 'evidence/raw/m10/producer/result.json'
    assert sha(reference) == parent_receipt['raw_sha256']['m10/producer/result.json']
    baseline = json.loads(reference.read_text())
    assert baseline['arm'] == {'m': 10, 'd': 13}
    assert baseline['physical_f0_points'] == 7977 and baseline['nonzero_signed_columns'] == 3988
    assert baseline['ordered_physical_tuples'] == 7977 ** 10 and baseline['q'] == data['q']
    return frozen


def tree_hashes(path: Path) -> dict[str, str]:
    return {str(file.relative_to(path)): sha(file) for file in sorted(path.rglob('*')) if file.is_file()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    raw = args.out / 'raw'
    raw.mkdir()
    frozen = freeze()
    data = json.loads((HERE / 'INPUT.json').read_text())
    receipt = {'status': 'started', 'started_utc': utc(),
               'python': sys.version, 'platform': platform.platform(),
               'freeze_sha256': sha(HERE / 'FROZEN.json'),
               'input_sha256': {label: frozen[label + '_sha256'] for label in FILES},
               'parent_balanced_frozen_sha256': data['parent_balanced_frozen_sha256'],
               'parent_balanced_receipt_sha256': data['parent_balanced_receipt_sha256'],
               'commands': []}
    try:
        for arm in data['arms']:
            m = arm['m']
            directory = raw / f'm{m}'
            commands = [(f'm{m}-producer', [sys.executable, str(HERE / 'producer.py'),
                                           '--m', str(m), '--out', str(directory / 'producer')],
                         data['caps']['producer_wall_seconds']),
                        (f'm{m}-verify', [sys.executable, str(HERE / 'verify.py'),
                                         '--m', str(m), '--source', str(directory / 'producer'),
                                         '--out', str(directory / 'verify.json')],
                         data['caps']['verifier_wall_seconds'])]
            directory.mkdir()
            for name, argv, cap in commands:
                item = {'name': name, 'argv': argv, 'started_utc': utc(),
                        'timeout_seconds': cap + 30}
                receipt['commands'].append(item)
                with (args.out / f'{name}.stdout.txt').open('w') as stdout, \
                        (args.out / f'{name}.stderr.txt').open('w') as stderr:
                    try:
                        process = subprocess.run(argv, stdout=stdout, stderr=stderr,
                                                 timeout=cap + 30, check=False)
                        item['exit_code'] = process.returncode
                    except subprocess.TimeoutExpired:
                        item['exit_code'] = 'TIMEOUT'
                item['finished_utc'] = utc()
                item['peak_child_rss_bytes_upper'] = peak_child_rss()
                if item['exit_code'] != 0 or item['peak_child_rss_bytes_upper'] > data['caps']['rss_bytes']:
                    raise RuntimeError(f"{name}: exit/RSS cap {item['exit_code']}, {item['peak_child_rss_bytes_upper']}")
        receipt['status'] = 'success'
    except Exception as error:
        receipt['status'] = 'failed'
        receipt['failure'] = repr(error)
    finally:
        receipt['finished_utc'] = utc()
        receipt['raw_sha256'] = tree_hashes(raw)
        receipt['stream_sha256'] = {path.name: sha(path) for path in args.out.iterdir()
                                    if path.is_file()}
        (args.out / 'receipt.json').write_text(
            json.dumps(receipt, sort_keys=True, separators=(',', ':')) + '\n')
    return int(receipt['status'] != 'success')


if __name__ == '__main__':
    raise SystemExit(main())
