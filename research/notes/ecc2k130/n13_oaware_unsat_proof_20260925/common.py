#!/usr/bin/env python3
"""Frozen input and independent DRAT-checker preflight for n13 proof gate."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
PARENT = NOTES / 'n13_oaware_sat_benchmark_20260925'
CNF = NOTES / 'rotated_s3_o_branch_20260925/evidence/producer/n13-m5'
CHECKER = HERE / 'checker/drat-trim.c'
FILES = {
    'protocol': HERE / 'PROTOCOL.md',
    'input': HERE / 'INPUT.json',
    'common': HERE / 'common.py',
    'run': HERE / 'run.py',
    'ci_replay': HERE / 'ci_replay.py',
    'checker_source': CHECKER,
    'checker_license': HERE / 'checker/LICENSE',
    'workflow': HERE.parents[3] / '.github/workflows/ecc2k130-n13-oaware-unsat-proof.yml',
    'preoutcome_failure': HERE / 'PREOUTCOME_PREFLIGHT_FAILURE.md',
}
PARENT_FILES = {
    'parent_frozen': PARENT / 'FROZEN.json',
    'parent_input': PARENT / 'INPUT.json',
    'parent_panel_result': PARENT / 'evidence/panel/result.json',
    'parent_manifest': PARENT / 'evidence/MANIFEST.json',
    'base_cnf': CNF / 'base.cnf',
    'schema': CNF / 'schema.json',
}
EXPECTED_CASES = [('Q0T0', 676), ('Q4T0', 629), ('Q4T1', 504),
                  ('Q4T2', 985), ('Q4T3', 690)]


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


def query_bytes(base: bytes, literal: int, variables: int, clauses: int) -> bytes:
    first, header, rest = base.split(b'\n', 2)
    assert first.startswith(b'c ')
    assert header == f'p cnf {variables} {clauses}'.encode()
    assert rest.endswith(b'\n') and 1 <= literal <= variables
    return first + b'\n' + f'p cnf {variables} {clauses + 1}\n'.encode() + \
        rest + f'{literal} 0\n'.encode()


def freeze() -> tuple[dict, dict, bytes, dict]:
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    for label, path in FILES.items():
        assert sha(path) == frozen[label + '_sha256'], label
    data = json.loads((HERE / 'INPUT.json').read_text())
    assert data['domain'] == 'ECC2K130-N13-OAWARE-UNSAT-PROOF-20260925-v1'
    assert data['parent_main_merge'] == '91eb1567ab6d02f4af4e090ff091c7b6e0d1adff'
    assert [(row['id'], row['literal']) for row in data['cases']] == EXPECTED_CASES
    assert data['positive_control']['id'] == 'Q0T3'
    assert data['positive_control']['literal'] == 596
    assert data['caps'] == {'solver_wall_seconds': 60, 'checker_wall_seconds': 180,
                            'rss_bytes': 1024 * 1024 * 1024,
                            'proof_bytes': 128 * 1024 * 1024,
                            'archive_bytes': 512 * 1024 * 1024,
                            'temporary_bytes': 320 * 1024 * 1024}
    assert data['solver_options'] == ['--no-binary']
    assert data['checker_upstream_commit'] == '2e3b2dc0ecf938addbd779d42877b6ed69d9a985'
    assert data['checker_build'] == ['cc', '-std=c99', '-O2',
        'research/notes/ecc2k130/n13_oaware_unsat_proof_20260925/checker/drat-trim.c',
        '-o', '/private/tmp/ecc2k130-drat-trim-vendored-20260925']
    for label, path in PARENT_FILES.items():
        assert sha(path) == data[label + '_sha256'], label
    assert sha(CHECKER) == data['checker_source_sha256']
    assert sha(HERE / 'checker/LICENSE') == data['checker_license_sha256']
    parent_frozen = json.loads((PARENT / 'FROZEN.json').read_text())
    assert data['solver_binary_sha256'] == parent_frozen['binaries']['cadical']['sha256']
    assert data['solver_version'] == parent_frozen['binaries']['cadical']['version']
    assert data['base_cnf_sha256'] == parent_frozen['base_sha256']
    schema = json.loads((CNF / 'schema.json').read_text())
    assert (schema['field_degree'], schema['field_poly'], schema['m']) == (13, 0x201b, 5)
    assert (schema['variables'], schema['clauses']) == (1263, 1195344)
    parent_targets = {row['id']: row for row in json.loads((PARENT / 'INPUT.json').read_text())['targets']}
    base = (CNF / 'base.cnf').read_bytes()
    for case in data['cases'] + [data['positive_control']]:
        target = parent_targets[case['id']]
        assert target['literal'] == case['literal']
        assert target['point_oracle_positive'] == (case['id'] == 'Q0T3')
        prior = json.loads((PARENT / f"evidence/panel/{case['id']}-cadical.json").read_text())
        assert prior['id'] == case['id'] and prior['assumption_literal'] == case['literal']
        assert prior['verdict'] == ('SAT' if case['id'] == 'Q0T3' else 'UNSAT')
        raw = query_bytes(base, case['literal'], schema['variables'], schema['clauses'])
        assert sha_bytes(raw) == case['prior_query_sha256'] == prior['input_sha256']
    return data, frozen, base, schema


def checker_selftest(binary: Path, directory: Path) -> dict:
    """A fixed non-unit UNSAT XOR needs valid derived unit steps before the empty clause."""
    directory.mkdir(parents=True, exist_ok=True)
    cnf = directory / 'xor.cnf'
    good = directory / 'xor-valid.drat'
    bad = directory / 'xor-mutated.drat'
    cnf.write_text('p cnf 2 4\n1 2 0\n-1 2 0\n1 -2 0\n-1 -2 0\n')
    good.write_text('1 0\n2 0\n0\n')
    bad.write_text('0\n')
    results = {}
    for label, proof in (('valid', good), ('mutated', bad)):
        process = subprocess.run([str(binary), str(cnf), str(proof)],
                                 capture_output=True, timeout=15)
        stdout = process.stdout.decode('utf-8', errors='replace')
        (directory / f'{label}.stdout.txt').write_bytes(process.stdout)
        (directory / f'{label}.stderr.txt').write_bytes(process.stderr)
        results[label] = {'exit_code': process.returncode,
                          'verdict_verified': 's VERIFIED' in stdout,
                          'verdict_not_verified': 's NOT VERIFIED' in stdout}
    assert results['valid']['exit_code'] == 0 and results['valid']['verdict_verified']
    assert results['mutated']['exit_code'] != 0 and results['mutated']['verdict_not_verified']
    return {'cnf_sha256': sha(cnf), 'valid_proof_sha256': sha(good),
            'mutated_proof_sha256': sha(bad), 'results': results}
