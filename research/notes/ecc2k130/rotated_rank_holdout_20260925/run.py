#!/usr/bin/env python3
"""Frozen rank/point-only pipeline with durable failure and cost receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import resource
import subprocess
import sys
import tarfile
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'
FROZEN = HERE / 'FROZEN.json'
OLD_ARCHIVE = HERE.parent / 'rotated_row_certificate_20260925/evidence/raw.tar.gz'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preflight():
    match = re.search(r'Frozen rank manifest SHA-256: `([0-9a-f]{64})`',
                      (HERE / 'PROTOCOL.md').read_text())
    assert match and digest(FROZEN) == match.group(1)
    frozen = json.loads(FROZEN.read_text())
    assert digest(OLD_ARCHIVE) == frozen['source_archive_sha256']
    for name, sha in frozen['file_sha256'].items():
        assert digest(HERE / name) == sha, name
    return frozen


def hashes(path):
    return {str(p.relative_to(path)): digest(p) for p in sorted(path.rglob('*')) if p.is_file()}


def stage(argv, name, timeout, receipt):
    start, cpu = time.monotonic(), time.process_time()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    record = {'name': name, 'argv': argv, 'timeout_seconds': timeout}
    receipt['stages'].append(record)
    try:
        result = subprocess.run(argv, cwd=HERE, capture_output=True, text=True, timeout=timeout)
        record['returncode'] = result.returncode
        record['stdout'] = result.stdout[-16000:]
        record['stderr'] = result.stderr[-16000:]
        assert result.returncode == 0, (name, result.returncode, record['stderr'])
    except subprocess.TimeoutExpired as exc:
        record['timeout'] = True
        record['stdout'] = str(exc.stdout)[-16000:]
        record['stderr'] = str(exc.stderr)[-16000:]
        raise
    finally:
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        record['wall_seconds'] = time.monotonic() - start
        record['launcher_cpu_seconds'] = time.process_time() - cpu
        record['child_user_cpu_seconds'] = after.ru_utime - before.ru_utime
        record['child_system_cpu_seconds'] = after.ru_stime - before.ru_stime
        record['child_peak_rss_raw'] = after.ru_maxrss


def archive(raw, output):
    with tarfile.open(output, 'w:gz') as tar:
        for item in sorted(raw.rglob('*')):
            tar.add(item, arcname=str(item.relative_to(raw.parent)), recursive=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path, default=HERE / 'evidence')
    args = parser.parse_args()
    evidence = args.evidence.resolve()
    evidence.mkdir(parents=True, exist_ok=True)
    raw = evidence / 'raw'
    receipt = {'kind': 'rotated_rank_point_only_holdout', 'status': 'failed', 'stages': [],
               'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    try:
        assert not raw.exists() and not (evidence / 'raw.tar.gz').exists()
        frozen = preflight()
        receipt['frozen_manifest_sha256'] = digest(FROZEN)
        receipt['source_archive_sha256'] = frozen['source_archive_sha256']
        receipt['input_sha256'] = {name: digest(INPUTS / name)
                                   for name in ('manifest.json', 'point_only.json',
                                                'sealed_labels.json', 'construction_receipt.json')}
        receipt['input_construction_cost'] = json.loads((INPUTS / 'construction_receipt.json').read_text())
        assert receipt['input_construction_cost']['status'] == 'success'
        raw.mkdir()
        producer_start = time.monotonic()
        def remaining():
            budget = 120 - (time.monotonic() - producer_start)
            assert budget > 1, 'producer global 120-second cap exhausted'
            return int(budget)
        stage([sys.executable, str(HERE / 'rank_train.py'),
               '--points', str(INPUTS / 'point_only.json'), '--out', str(raw / 'training')],
              'rank_training', remaining(), receipt)
        stage([sys.executable, str(HERE / 'oracle.py'),
               '--points', str(INPUTS / 'point_only.json'),
               '--out', str(raw / 'oracle_responses.json')],
              'archive_oracle_full_scan', remaining(), receipt)
        # No sealed labels, archive path, or k-valued training summary in this argv.
        stage([sys.executable, str(HERE / 'recover.py'),
               '--points', str(INPUTS / 'point_only.json'),
               '--logs', str(raw / 'training' / 'base_logs.json'),
               '--oracle', str(raw / 'oracle_responses.json'),
               '--out', str(raw / 'recovery.json')],
              'point_only_recovery', remaining(), receipt)
        receipt['producer_total_wall_seconds'] = time.monotonic() - producer_start
        assert receipt['producer_total_wall_seconds'] <= 120
        stage([sys.executable, str(HERE / 'verify.py'), '--raw', str(raw),
               '--inputs', str(INPUTS), '--report', str(raw / 'verify_report.json')],
              'independent_replay_and_sealed_labels', 300, receipt)
        assert json.loads((raw / 'verify_report.json').read_text())['status'] == 'success'
        receipt['raw_file_sha256'] = hashes(raw)
        archive(raw, evidence / 'raw.tar.gz')
        receipt['archive_sha256'] = digest(evidence / 'raw.tar.gz')
        receipt['archive_bytes'] = (evidence / 'raw.tar.gz').stat().st_size
        receipt['status'] = 'success'
    except Exception as exc:
        receipt['error'] = repr(exc)
        receipt['traceback'] = traceback.format_exc()
        if raw.exists():
            receipt['partial_raw_file_sha256'] = hashes(raw)
        raise
    finally:
        receipt['finished_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        output = evidence / 'receipt.json'
        output.write_text(json.dumps(receipt, sort_keys=True, separators=(',', ':')) + '\n')
        (evidence / 'receipt.json.sha256').write_text(digest(output) + '  receipt.json\n')


if __name__ == '__main__':
    main()
