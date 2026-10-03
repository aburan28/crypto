#!/usr/bin/env python3
"""Hash-only remaining-22 preregistration; then independent archived DRAT verification."""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import tempfile
from pathlib import Path

from common import HERE, CHECKER, FILES, checker_selftest, file_hashes, freeze, \
    query_bytes, sha, sha_bytes


def compile_checker(binary: Path) -> None:
    subprocess.run(['cc', '-std=c99', '-O2', str(CHECKER), '-o', str(binary)],
                   check=True, timeout=30, capture_output=True)
    assert binary.is_file()


def run_check(binary: Path, query: Path, proof: Path, cap: int) -> tuple[int, str]:
    process = subprocess.run([str(binary), str(query), str(proof)],
                             capture_output=True, timeout=cap + 15)
    return process.returncode, process.stdout.decode('utf-8', errors='replace')


def archive_replay(evidence: Path, checker: Path, data: dict, frozen: dict,
                   base: bytes, schema: dict) -> None:
    receipt = json.loads((evidence / 'receipt.json').read_text())
    assert receipt['status'] == 'success'
    assert receipt['freeze_sha256'] == sha(HERE / 'FROZEN.json')
    assert receipt['source_sha256'] == frozen
    assert receipt['parent_main_merge'] == data['parent_main_merge']
    assert receipt['parent_proof_head'] == data['parent_proof_head']
    assert receipt['parent_proof_merge_commit'] == data['parent_proof_merge_commit']
    assert receipt['parent_proof_receipt_sha256'] == data['parent_proof_receipt_sha256']
    assert receipt['checker_upstream_commit'] == data['checker_upstream_commit']
    assert receipt['solver_binary_sha256'] == data['solver_binary_sha256']
    assert receipt['checker_binary_sha256'] == data['checker_mac_binary_sha256']
    assert file_hashes(evidence) == receipt['files_sha256']
    assert receipt['archive_bytes'] <= data['caps']['archive_bytes']
    assert sum(path.stat().st_size for path in evidence.rglob('*') if path.is_file()) <= data['caps']['archive_bytes']
    assert [item['id'] for item in receipt['cases']] == [case['id'] for case in data['cases']]
    assert len(receipt['commands']) == 2 * len(data['cases']) + 1
    assert [entry['name'] for entry in receipt['commands']] == [
        name for case in data['cases'] for name in (f"{case['id']}-solver", f"{case['id']}-checker")
    ] + ['Q0T1-proof-on-Q0T3']
    for entry in receipt['commands']:
        cap = data['caps']['solver_wall_seconds'] if entry['name'].endswith('-solver') \
            else data['caps']['checker_wall_seconds']
        assert entry['wall_seconds'] <= cap
        if entry['name'].endswith('-solver'):
            assert entry['exit_code'] == 20
        elif entry['name'] == 'Q0T1-proof-on-Q0T3':
            assert entry['exit_code'] != 0
        else:
            assert entry['exit_code'] == 0
        assert entry['peak_child_rss_bytes_upper'] <= data['caps']['rss_bytes']
        for stream in ('stdout', 'stderr'):
            assert sha(evidence / f"{entry['name']}.{stream}.txt") == entry[f'{stream}_sha256']
    control = receipt['mutated_proof_control']
    assert control['results']['valid']['exit_code'] == 0
    assert control['results']['valid']['verdict_verified']
    assert control['results']['mutated']['exit_code'] != 0
    assert control['results']['mutated']['verdict_not_verified']
    with tempfile.TemporaryDirectory(prefix='n13-unsat-replay-') as tmp:
        scratch = Path(tmp)
        first_proof = None
        for spec, case in zip(data['cases'], receipt['cases'], strict=True):
            ident = spec['id']
            assert case['status'] == 'PASS' and case['literal'] == spec['literal']
            assert case['point_oracle_positive'] is False
            query_raw = query_bytes(base, spec['literal'], schema['variables'], schema['clauses'])
            assert sha_bytes(query_raw) == case['query_sha256'] == spec['prior_query_sha256']
            assert len(query_raw) == case['query_bytes']
            query = scratch / f'{ident}.cnf'
            query.write_bytes(query_raw)
            compressed = evidence / 'proofs' / f'{ident}.drat.gz'
            assert sha(compressed) == case['gzip_proof_sha256']
            assert compressed.stat().st_size == case['gzip_proof_bytes']
            proof_raw = gzip.decompress(compressed.read_bytes())
            assert 0 < len(proof_raw) == case['raw_proof_bytes'] <= data['caps']['proof_bytes']
            assert sha_bytes(proof_raw) == case['raw_proof_sha256']
            proof = scratch / f'{ident}.drat'
            proof.write_bytes(proof_raw)
            assert sum(path.stat().st_size for path in scratch.rglob('*') if path.is_file()) <= data['caps']['temporary_bytes']
            archived_solver = (evidence / f'{ident}-solver.stdout.txt').read_text(errors='replace')
            archived_checker = (evidence / f'{ident}-checker.stdout.txt').read_text(errors='replace')
            assert 's UNSATISFIABLE' in archived_solver
            assert 's VERIFIED' in archived_checker and 's NOT VERIFIED' not in archived_checker
            code, text = run_check(checker, query, proof, data['caps']['checker_wall_seconds'])
            assert code == 0 and 's VERIFIED' in text and 's NOT VERIFIED' not in text, ident
            query.unlink()
            if ident == 'Q0T1':
                first_proof = proof
            else:
                proof.unlink()
        assert first_proof is not None
        positive = data['positive_control']
        positive_raw = query_bytes(base, positive['literal'], schema['variables'], schema['clauses'])
        assert sha_bytes(positive_raw) == positive['prior_query_sha256']
        positive_query = scratch / 'Q0T3.cnf'
        positive_query.write_bytes(positive_raw)
        assert sum(path.stat().st_size for path in scratch.rglob('*') if path.is_file()) <= data['caps']['temporary_bytes']
        wrong = receipt['wrong_input_control']
        assert wrong == {'proof_id': 'Q0T1', 'sat_query_id': 'Q0T3',
                         'sat_query_sha256': positive['prior_query_sha256'],
                         'command': 'Q0T1-proof-on-Q0T3'}
        code, text = run_check(checker, positive_query, first_proof,
                               data['caps']['checker_wall_seconds'])
        assert code != 0 and 's NOT VERIFIED' in text
    print('All 22 remaining fixed n13 Q+T DRAT proofs externally verified; mutated and wrong-input controls rejected.')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path)
    args = parser.parse_args()
    data, frozen, base, schema = freeze()
    for path in FILES.values():
        if path.suffix == '.py':
            compile(path.read_bytes(), str(path), 'exec')
    with tempfile.TemporaryDirectory(prefix='n13-checker-build-') as tmp:
        checker = Path(tmp) / 'drat-trim'
        compile_checker(checker)
        checker_selftest(checker, Path(tmp) / 'control')
        if args.evidence is None:
            print('Frozen #790 dependency, all 22 #785 inputs, source, external checker and mutated-proof preflight PASS; no selected proof outcome read.')
        else:
            archive_replay(args.evidence, checker, data, frozen, base, schema)


if __name__ == '__main__':
    main()
