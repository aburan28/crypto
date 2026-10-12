#!/usr/bin/env python3
"""Hash-only preregistration gate, then independently replay committed evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preflight():
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    assert frozen['schema'] == 'rotated_joint_frozen_v1'
    for name, expected in frozen['files_sha256'].items():
        path = HERE / name
        assert path.is_file() and sha(path) == expected, name
        if name.endswith('.py'):
            compile(path.read_bytes(), str(path), 'exec')
    for name, expected in frozen['external_sha256'].items():
        path = NOTES / name
        assert path.is_file() and sha(path) == expected, name
    training = json.loads((HERE / 'inputs/training.json').read_text())
    point_only = json.loads((HERE / 'inputs/point_only.json').read_text())
    sealed = json.loads((HERE / 'inputs/sealed_labels.json').read_text())
    assert len(training['targets']) == 256 and len(point_only['targets']) == len(sealed['targets']) == 64
    assert not ({c['k'] for c in training['targets']} & {c['k'] for c in sealed['targets']})
    assert not any('k' in c for c in point_only['targets'])
    return frozen


def evidence_replay(folder):
    receipt = json.loads((folder / 'receipt.json').read_text())
    assert receipt['status'] == 'success'
    archive = folder / 'raw.tar.gz'
    assert sha(archive) == receipt['raw_archive_sha256']
    assert archive.stat().st_size == receipt['raw_archive_bytes']
    with tempfile.TemporaryDirectory(prefix='joint-rank-ci-') as temp:
        root = Path(temp)
        with tarfile.open(archive, 'r:gz') as tar:
            for member in tar.getmembers():
                path = Path(member.name)
                assert not path.is_absolute() and '..' not in path.parts
                assert member.isfile() or member.isdir()
            tar.extractall(root, filter='data')
        raw = root / 'raw'
        assert {str(p.relative_to(raw)): sha(p) for p in sorted(raw.rglob('*')) if p.is_file()} == receipt['raw_sha256']
        rerun = root / 'ci_verification.json'
        subprocess.run([sys.executable, str(HERE / 'verify.py'), '--raw', str(raw),
                        '--inputs', str(HERE / 'inputs'), '--out', str(rerun)],
                       check=True, timeout=300)
        prior = json.loads((raw / 'verification.json').read_text())
        actual = json.loads(rerun.read_text())
        for key in ('status', 'training_targets', 'holdout_targets', 'beta3_rank',
                    'beta3_columns', 'joint_rank', 'joint_columns', 'beta3_recovered',
                    'joint_recovered', 'full_rank_gate', 'semantic_followup_gate'):
            assert actual[key] == prior[key], key
        assert receipt['decision'] == {key: prior[key] for key in
               ('full_rank_gate', 'semantic_followup_gate')}
    print('Independent n19 four-base joint-rank and sealed holdout replay PASS.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path)
    args = parser.parse_args()
    preflight()
    if args.evidence:
        evidence_replay(args.evidence)
    else:
        print('Frozen protocol, source, archives and disjoint point-only split PASS; no rank outcome read.')


if __name__ == '__main__':
    main()
