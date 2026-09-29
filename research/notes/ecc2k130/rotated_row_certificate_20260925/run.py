#!/usr/bin/env python3
"""Execute frozen row producer and independent replay; preserve all failures."""
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
EVIDENCE = HERE / 'evidence'
PROTOCOL = HERE / 'PROTOCOL.md'
FROZEN = HERE / 'FROZEN.json'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_freeze() -> dict:
    text = PROTOCOL.read_text()
    match = re.search(r'Frozen row manifest SHA-256: `([0-9a-f]{64})`', text)
    assert match and sha(FROZEN) == match.group(1)
    manifest = json.loads(FROZEN.read_text())
    for name, expected in manifest['file_sha256'].items():
        assert sha(HERE / name) == expected, name
    old = HERE.parent / 'rotated_subspace_support_20260925/evidence/raw.tar.gz'
    assert sha(old) == manifest['old_archive_sha256']
    return manifest


def hashes(folder: Path):
    return {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}


def run_child(argv, label: str, timeout: int, receipt: dict):
    start_wall, start_cpu = time.monotonic(), time.process_time()
    previous = resource.getrusage(resource.RUSAGE_CHILDREN)
    stage = {'name': label, 'argv': argv, 'timeout_seconds': timeout}
    receipt['stages'].append(stage)
    try:
        proc = subprocess.run(argv, cwd=HERE, capture_output=True, text=True, timeout=timeout)
        stage['returncode'] = proc.returncode
        stage['stdout'] = proc.stdout[-16000:]
        stage['stderr'] = proc.stderr[-16000:]
        assert proc.returncode == 0, (label, proc.returncode, stage['stderr'])
    except subprocess.TimeoutExpired as exc:
        stage['timeout'] = True
        stage['stdout'] = str(exc.stdout)[-16000:]
        stage['stderr'] = str(exc.stderr)[-16000:]
        raise
    finally:
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        stage['wall_seconds'] = time.monotonic() - start_wall
        stage['launcher_cpu_seconds'] = time.process_time() - start_cpu
        stage['child_user_cpu_seconds'] = after.ru_utime - previous.ru_utime
        stage['child_system_cpu_seconds'] = after.ru_stime - previous.ru_stime
        # ru_maxrss is a child high-water value across the launcher lifetime, not additive.
        stage['child_peak_rss_raw'] = after.ru_maxrss


def archive_raw(raw: Path, output: Path):
    with tarfile.open(output, 'w:gz') as tar:
        for item in sorted(raw.rglob('*')):
            tar.add(item, arcname=str(item.relative_to(raw.parent)), recursive=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path, default=EVIDENCE)
    args = parser.parse_args()
    evidence = args.evidence.resolve()
    evidence.mkdir(parents=True, exist_ok=True)
    raw = evidence / 'raw'
    receipt = {'kind': 'rotated_compressed_row_certificate', 'status': 'failed',
               'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               'stages': []}
    try:
        assert not raw.exists() and not (evidence / 'raw.tar.gz').exists()
        manifest = verify_freeze()
        receipt['frozen_manifest_sha256'] = sha(FROZEN)
        receipt['old_archive_sha256'] = manifest['old_archive_sha256']
        receipt['planted_input_sha256'] = manifest['file_sha256']['planted_inputs.json']
        run_child([sys.executable, str(HERE / 'row.py'), '--out', str(raw),
                   '--inputs', str(HERE / 'planted_inputs.json')],
                  'producer_all_cells', 600, receipt)
        assert raw.exists()
        run_child([sys.executable, str(HERE / 'replay.py'), '--raw', str(raw),
                   '--inputs', str(HERE / 'planted_inputs.json'),
                   '--report', str(raw / 'verify_report.json')],
                  'independent_full_replay', 600, receipt)
        report = json.loads((raw / 'verify_report.json').read_text())
        assert report['status'] == 'success' and report['total_wall_seconds'] <= 600
        receipt['raw_file_sha256'] = hashes(raw)
        archive_raw(raw, evidence / 'raw.tar.gz')
        receipt['archive_sha256'] = sha(evidence / 'raw.tar.gz')
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
        receipt_path = evidence / 'receipt.json'
        receipt_path.write_text(json.dumps(receipt, sort_keys=True, separators=(',', ':')) + '\n')
        (evidence / 'receipt.json.sha256').write_text(sha(receipt_path) + '  receipt.json\n')


if __name__ == '__main__':
    main()
