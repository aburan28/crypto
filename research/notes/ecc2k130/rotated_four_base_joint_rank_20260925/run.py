#!/usr/bin/env python3
"""Bounded cold children and complete raw-evidence receipt for joint rank."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import platform
import resource
import subprocess
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n')


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def child(name, argv, cap, evidence, receipt):
    start = utc()
    wall = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        result = subprocess.run(argv, capture_output=True, timeout=cap, check=False)
        code = result.returncode
        stdout, stderr = result.stdout, result.stderr
        status = 'success' if code == 0 else 'producer_failure'
    except subprocess.TimeoutExpired as exc:
        code = None
        stdout, stderr = exc.stdout or b'', exc.stderr or b''
        status = 'timeout'
    (evidence / f'{name}.stdout.txt').write_bytes(stdout)
    (evidence / f'{name}.stderr.txt').write_bytes(stderr)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    record = {'name': name, 'argv': list(map(str, argv)), 'status': status,
              'exit_code': code, 'started_utc': start, 'finished_utc': utc(),
              'wall_seconds': time.monotonic() - wall, 'timeout_seconds': cap,
              'child_cpu_seconds': (after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime),
              'child_peak_rss_bytes_upper': after.ru_maxrss if sys.platform == 'darwin' else after.ru_maxrss * 1024}
    receipt['children'].append(record)
    if status != 'success':
        raise RuntimeError(f'{name}: {status}, exit={code}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path, required=True)
    args = parser.parse_args()
    assert not args.evidence.exists()
    args.evidence.mkdir(parents=True)
    raw = args.evidence / 'raw'
    raw.mkdir()
    receipt = {'schema': 'rotated_joint_run_receipt_v1', 'status': 'running',
               'started_utc': utc(), 'python': sys.version, 'platform': platform.platform(),
               'children': [], 'source_input_sha256': {}}
    for name in ('training.json', 'point_only.json', 'sealed_labels.json', 'construction_receipt.json'):
        receipt['source_input_sha256'][name] = sha(INPUTS / name)
    receipt['source_input_sha256']['FROZEN.json'] = sha(HERE / 'FROZEN.json')
    try:
        child('training', [sys.executable, str(HERE / 'train.py'), '--training',
              str(INPUTS / 'training.json'), '--out', str(raw / 'training')],
              120, args.evidence, receipt)
        child('oracle', [sys.executable, str(HERE / 'oracle.py'), '--points',
              str(INPUTS / 'point_only.json'), '--out', str(raw / 'oracle.json')],
              120, args.evidence, receipt)
        child('point_recovery', [sys.executable, str(HERE / 'recover.py'), '--points',
              str(INPUTS / 'point_only.json'), '--logs', str(raw / 'training/base_logs.json'),
              '--oracle', str(raw / 'oracle.json'), '--out', str(raw / 'recovery.json')],
              120, args.evidence, receipt)
        child('independent_replay', [sys.executable, str(HERE / 'verify.py'), '--raw',
              str(raw), '--inputs', str(INPUTS), '--out', str(raw / 'verification.json')],
              300, args.evidence, receipt)
        verify = json.loads((raw / 'verification.json').read_text())
        assert verify['status'] == 'PASS'
        with tarfile.open(args.evidence / 'raw.tar.gz', 'w:gz') as archive:
            archive.add(raw, arcname='raw')
        receipt['raw_sha256'] = {str(path.relative_to(raw)): sha(path) for path in sorted(raw.rglob('*')) if path.is_file()}
        receipt['raw_archive_sha256'] = sha(args.evidence / 'raw.tar.gz')
        receipt['raw_archive_bytes'] = (args.evidence / 'raw.tar.gz').stat().st_size
        receipt['decision'] = {'full_rank_gate': verify['full_rank_gate'],
                               'semantic_followup_gate': verify['semantic_followup_gate']}
        receipt['status'] = 'success'
    except Exception as error:
        receipt['status'] = 'failure'
        receipt['error'] = repr(error)
        raise
    finally:
        receipt['finished_utc'] = utc()
        save(args.evidence / 'receipt.json', receipt)


if __name__ == '__main__':
    main()
