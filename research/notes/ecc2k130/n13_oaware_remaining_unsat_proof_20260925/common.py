#!/usr/bin/env python3
"""Frozen inputs for the 22 remaining #785 n13 negative proof queries."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
PARENT = NOTES / 'n13_oaware_unsat_proof_20260925'
SAT = NOTES / 'n13_oaware_sat_benchmark_20260925'
CHECKER = PARENT / 'checker/drat-trim.c'
WORKFLOW = HERE.parents[3] / '.github/workflows/ecc2k130-n13-oaware-remaining-unsat-proof.yml'
FILES = {
    'protocol': HERE / 'PROTOCOL.md',
    'input': HERE / 'INPUT.json',
    'common': HERE / 'common.py',
    'run': HERE / 'run.py',
    'ci_replay': HERE / 'ci_replay.py',
    'workflow': WORKFLOW,
}
EXPECTED_IDS = ('Q0T1', 'Q0T2', 'Q1T2', 'Q1T3', 'Q2T0', 'Q2T1', 'Q2T3',
                'Q3T1', 'Q3T2', 'Q3T3', 'Q5T0', 'Q5T1', 'Q5T2', 'Q5T3',
                'Q6T0', 'Q6T1', 'Q6T2', 'Q6T3', 'Q7T0', 'Q7T1', 'Q7T2', 'Q7T3')
PARENT_HEAD = '3f8e052da9265f88dbcf58035c832140b0a0b689'


def _parent_common():
    spec = importlib.util.spec_from_file_location('n13_parent_proof_common', PARENT / 'common.py')
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


parent_common = _parent_common()
query_bytes = parent_common.query_bytes
checker_selftest = parent_common.checker_selftest


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def file_hashes(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): sha(path) for path in sorted(root.rglob('*'))
            if path.is_file() and path.name != 'receipt.json'}


def freeze() -> tuple[dict, dict, bytes, dict]:
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    for label, path in FILES.items():
        assert sha(path) == frozen[label + '_sha256'], label
    data = json.loads((HERE / 'INPUT.json').read_text())
    assert data['domain'] == 'ECC2K130-N13-OAWARE-REMAINING-UNSAT-PROOF-20260925-v1'
    assert data['parent_proof_head'] == PARENT_HEAD
    merge_commit = data['parent_proof_merge_commit']
    assert merge_commit is None or (isinstance(merge_commit, str) and len(merge_commit) == 40
                                    and all(ch in '0123456789abcdef' for ch in merge_commit))
    assert sha(PARENT / 'FROZEN.json') == data['parent_proof_frozen_sha256']
    assert sha(PARENT / 'evidence/receipt.json') == data['parent_proof_receipt_sha256']
    parent_data, _, base, schema = parent_common.freeze()
    for key, value in parent_data.items():
        if key not in ('domain', 'cases'):
            assert data[key] == value, key
    parent_receipt = json.loads((PARENT / 'evidence/receipt.json').read_text())
    assert parent_receipt['status'] == 'success'
    assert [row['id'] for row in parent_receipt['cases']] == [row['id'] for row in parent_data['cases']]
    assert all(row['status'] == 'PASS' for row in parent_receipt['cases'])
    targets = json.loads((SAT / 'INPUT.json').read_text())['targets']
    negative = [row for row in targets if not row['point_oracle_positive']]
    positive = [row for row in targets if row['point_oracle_positive']]
    selected = {row['id'] for row in parent_data['cases']}
    remaining = [row for row in negative if row['id'] not in selected]
    assert len(negative) == 27 and len(positive) == 5
    assert len(selected) == 5 and len(remaining) == 22
    assert tuple(row['id'] for row in remaining) == EXPECTED_IDS
    assert [row['id'] for row in data['cases']] == list(EXPECTED_IDS)
    for target, case in zip(remaining, data['cases'], strict=True):
        assert case['literal'] == target['literal']
        panel = json.loads((SAT / f"evidence/panel/{case['id']}-cadical.json").read_text())
        assert panel['id'] == case['id'] and panel['assumption_literal'] == case['literal']
        assert panel['verdict'] == 'UNSAT' and panel['point_oracle_positive'] is False
        raw = query_bytes(base, case['literal'], schema['variables'], schema['clauses'])
        assert sha_bytes(raw) == case['prior_query_sha256'] == panel['input_sha256']
    return data, frozen, base, schema
