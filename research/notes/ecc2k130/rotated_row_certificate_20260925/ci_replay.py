#!/usr/bin/env python3
"""Check preregistration hashes and replay archived rows without producer execution."""
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
OLD = HERE.parent / 'rotated_subspace_support_20260925/evidence/raw.tar.gz'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preflight():
    match = re.search(r'Frozen row manifest SHA-256: `([0-9a-f]{64})`',
                      (HERE / 'PROTOCOL.md').read_text())
    assert match and digest(HERE / 'FROZEN.json') == match.group(1)
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    assert frozen['old_archive_sha256'] == digest(OLD)
    for name, sha in frozen['file_sha256'].items():
        assert digest(HERE / name) == sha, name
    with tempfile.TemporaryDirectory() as tempdir:
        generated = Path(tempdir) / 'inputs.json'
        construction_receipt = Path(tempdir) / 'receipt.json'
        subprocess.run([sys.executable, str(HERE / 'make_inputs.py'),
                        '--out', str(generated), '--receipt', str(construction_receipt)],
                       check=True, timeout=65)
        assert digest(generated) == frozen['file_sha256']['planted_inputs.json']
        assert json.loads(construction_receipt.read_text())['status'] == 'success'
    print('frozen row protocol/source/input and constructor replay: PASS')
    return frozen


def replay_archive(frozen):
    receipt_path = EVIDENCE / 'receipt.json'
    archive = EVIDENCE / 'raw.tar.gz'
    assert digest(receipt_path) == (EVIDENCE / 'receipt.json.sha256').read_text().split()[0]
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'success'
    assert receipt['frozen_manifest_sha256'] == digest(HERE / 'FROZEN.json')
    assert receipt['old_archive_sha256'] == frozen['old_archive_sha256']
    assert receipt['planted_input_sha256'] == frozen['file_sha256']['planted_inputs.json']
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
        actual = {str(path.relative_to(raw)): digest(path)
                  for path in sorted(raw.rglob('*')) if path.is_file()}
        assert actual == receipt['raw_file_sha256']
        report = root / 'independent_replay.json'
        subprocess.run([sys.executable, str(HERE / 'replay.py'), '--raw', str(raw),
                        '--inputs', str(HERE / 'planted_inputs.json'),
                        '--report', str(report)], check=True, timeout=600)
        original = json.loads((raw / 'verify_report.json').read_text())
        independently = json.loads(report.read_text())
        for field in ('status', 'old_archive_sha256', 'planted_input_sha256',
                      'toy', 'planted'):
            if field in ('toy', 'planted'):
                def scrub(rows):
                    return [{key: value for key, value in row.items()
                             if key != 'wall_seconds'} for row in rows]
                assert scrub(original[field]) == scrub(independently[field])
            else:
                assert original[field] == independently[field]
        print('archive-only independent full-row replay: PASS')


def main():
    frozen = preflight()
    if (HERE / 'RESULT.md').exists():
        assert (EVIDENCE / 'raw.tar.gz').exists()
        assert (EVIDENCE / 'receipt.json').exists()
    if (EVIDENCE / 'raw.tar.gz').exists():
        replay_archive(frozen)
    else:
        print('no outcome archive yet; hash-only preregistration PASS')


if __name__ == '__main__':
    main()
