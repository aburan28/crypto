#!/usr/bin/env python3
"""Hash-only prereg CI, then full archived rank/holdout independent replay."""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / 'evidence'
INPUTS = HERE / 'inputs'
SOURCE = HERE.parent / 'rotated_row_certificate_20260925/evidence/raw.tar.gz'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze_check():
    match = re.search(r'Frozen rank manifest SHA-256: `([0-9a-f]{64})`',
                      (HERE / 'PROTOCOL.md').read_text())
    assert match and digest(HERE / 'FROZEN.json') == match.group(1)
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    assert digest(SOURCE) == frozen['source_archive_sha256']
    for name, sha in frozen['file_sha256'].items():
        assert digest(HERE / name) == sha, name
    with tempfile.TemporaryDirectory() as tempdir:
        new_inputs = Path(tempdir) / 'inputs'
        subprocess.run([sys.executable, str(HERE / 'select_inputs.py'),
                        '--out', str(new_inputs)], check=True, timeout=35)
        for name in ('manifest.json', 'point_only.json', 'sealed_labels.json'):
            assert digest(new_inputs / name) == frozen['file_sha256']['inputs/' + name]
    print('frozen source/split and input regeneration: PASS')
    return frozen


def archive_check(frozen):
    archive = EVIDENCE / 'raw.tar.gz'
    receipt_path = EVIDENCE / 'receipt.json'
    assert digest(receipt_path) == (EVIDENCE / 'receipt.json.sha256').read_text().split()[0]
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'success'
    assert receipt['frozen_manifest_sha256'] == digest(HERE / 'FROZEN.json')
    assert receipt['source_archive_sha256'] == frozen['source_archive_sha256']
    assert receipt['archive_sha256'] == digest(archive)
    assert receipt['archive_bytes'] == archive.stat().st_size
    with tempfile.TemporaryDirectory() as tempdir:
        root = Path(tempdir)
        with tarfile.open(archive, 'r:gz') as tar:
            names = tar.getnames()
            assert len(names) == len(set(names))
            assert all(name == 'raw' or name.startswith('raw/') for name in names)
            tar.extractall(root, filter='data')
        raw = root / 'raw'
        assert {str(p.relative_to(raw)): digest(p) for p in sorted(raw.rglob('*')) if p.is_file()} == receipt['raw_file_sha256']
        report = root / 'fresh_verify.json'
        subprocess.run([sys.executable, str(HERE / 'verify.py'), '--raw', str(raw),
                        '--inputs', str(INPUTS), '--report', str(report)],
                       check=True, timeout=300)
        original = json.loads((raw / 'verify_report.json').read_text())
        fresh = json.loads(report.read_text())
        for key in ('status', 'source_archive_sha256', 'input_hashes', 'arms'):
            assert original[key] == fresh[key]
        print('archive-only independent rank, point-log and sealed-label replay: PASS')


def main():
    frozen = freeze_check()
    if (HERE / 'RESULT.md').exists():
        assert (EVIDENCE / 'raw.tar.gz').exists() and (EVIDENCE / 'receipt.json').exists()
    if (EVIDENCE / 'raw.tar.gz').exists():
        archive_check(frozen)
    else:
        print('no outcome archive; hash-only preregistration PASS')


if __name__ == '__main__':
    main()
