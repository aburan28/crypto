#!/usr/bin/env python3
"""Hash-only preflight before release; archive replay if evidence is committed."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = HERE.parent / 'symbolic_dag_fullpoint_20260925'
ARCHIVED_CADICAL = ROOT / 'research/notes/ecc2k130/n13_oaware_sat_benchmark_20260925/evidence/panel/Q0T3-cadical.stdout'
DOMAIN = 'k0-symbolic-dag-dimacs-gate-v1'
SOURCES = ('PROTOCOL.md', 'INPUT.json', 'export.py', 'verify.py', 'bounded.py',
           'produce.py', 'run.py', 'ci_replay.py', 'third_party/drat-trim.c',
           'third_party/LICENSE')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze() -> dict:
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    inputs = json.loads((HERE / 'INPUT.json').read_text())
    assert frozen['domain'] == DOMAIN and inputs['domain'] == DOMAIN
    assert inputs['parent_pr'] == 802
    assert inputs['parent_exact_head'] == 'a7d392229800683636adab7a1a88a670fab74508'
    assert inputs['field'] == {'n': 13, 'modulus': 0x201B}
    assert inputs['n131_single_edge_modulus'] == 2722258935367507707706996859454145699847
    assert len(inputs['pairs']) == 9 and len(inputs['invalid']) == 2
    assert sha(PARENT / 'FROZEN.json') == frozen['parent_freeze_sha256']
    assert sha(PARENT / 'evidence/producer/rows.jsonl.gz') == frozen['parent_rows_sha256']
    assert inputs['parent_rows_sha256'] == frozen['parent_rows_sha256']
    assert sha(ARCHIVED_CADICAL) == frozen['archived_cadical_sat_output_sha256']
    for name, digest in frozen['preoutcome_failure_sha256'].items():
        assert sha(HERE / 'preoutcome_failure_0' / name) == digest
    assert set(frozen['source_sha256']) == set(SOURCES)
    for name in SOURCES:
        assert sha(HERE / name) == frozen['source_sha256'][name], name
    assert sha(ROOT / '.github/workflows/ecc2k130-symbolic-dag-dimacs-gate.yml') == frozen['workflow_sha256']
    assert frozen['toy_cnf_byte_cap'] == 1 << 20
    assert frozen['n13_cnf_byte_cap'] == 8 << 20
    assert frozen['n131_cnf_byte_cap'] == 64 << 20
    assert frozen['n131_node_cap'] == 500000
    assert frozen['n131_external_wall_cap_seconds'] == 90
    assert frozen['n131_rss_cap_bytes'] == 512 << 20
    assert frozen['proof_byte_cap'] == 128 << 20
    if frozen['release_main_head'] is not None:
        assert len(frozen['release_main_head']) == 40
    # Import the actual producer in hash-only CI; this catches parent/local
    # module shadowing before the protected outcome run.
    sys.path.insert(0, str(HERE))
    import produce  # noqa: F401
    from verify import parse_solver_output
    assert parse_solver_output('s SATISFIABLE\nv 1 -2\nv 3 0\n', 3, 10) == ('SAT', [1, 0, 1])
    assert parse_solver_output('s UNSATISFIABLE\n', 3, 20) == ('UNSAT', None)
    archived_status, archived_model = parse_solver_output(ARCHIVED_CADICAL.read_text(), 1263, 10)
    assert archived_status == 'SAT' and archived_model is not None and len(archived_model) == 1263
    for bad in ('s SATISFIABLE\nv 1 -2 3\n',
                's SATISFIABLE\nv 1 -2 3 0 0\n',
                's SATISFIABLE\nv 1 -2 3 0 1\n',
                's SATISFIABLE\nv 1 -2 1 0\n',
                's SATISFIABLE\nv 1 -2 0\n'):
        try:
            parse_solver_output(bad, 3, 10)
        except ValueError:
            pass
        else:
            raise AssertionError('malformed SAT-model control accepted')
    return frozen


def _panel_expected() -> dict[str, dict]:
    inputs = json.loads((HERE / 'INPUT.json').read_text())
    cases = {}
    for pair in inputs['pairs']:
        for suffix, expected, result_key in (('_sat', 'SAT', 'sat_r'),
                                             ('_wrong', 'UNSAT', 'unsat_r')):
            cases[pair['id'] + suffix] = {'expected': expected,
                                          'points': {key: tuple(pair[key]) for key in ('p', 'q')}
                                          | {'r': tuple(pair[result_key])}}
    for row in inputs['invalid']:
        cases[row['id']] = {'expected': 'UNSAT',
                            'points': {key: tuple(row[key]) for key in ('p', 'q', 'r')}}
    assert len(cases) == 20
    return cases


def check_evidence(receipt_path: Path, frozen: dict) -> dict:
    receipt_path = receipt_path.resolve()
    archive = receipt_path.parent
    receipt = json.loads(receipt_path.read_text())
    assert receipt['domain'] == DOMAIN
    assert receipt['freeze_sha256'] == sha(HERE / 'FROZEN.json')
    assert receipt['decision'] in ('PASS', 'FAIL_OR_CENSORED')
    manifest_path = archive / 'MANIFEST.json'
    assert sha(manifest_path) == receipt['manifest_sha256']
    manifest = json.loads(manifest_path.read_text())
    actual = {path.relative_to(archive).as_posix(): path for path in archive.rglob('*')
              if path.is_file() and path.relative_to(archive).as_posix() not in ('receipt.json', 'MANIFEST.json')}
    assert [row['path'] for row in manifest] == sorted(actual)
    for row in manifest:
        artifact = actual[row['path']]
        assert artifact.stat().st_size == row['bytes'] and sha(artifact) == row['sha256']
    attempts = receipt['attempts']
    assert [item['phase'] for item in attempts] == ['toy', 'panel', 'n131'][:len(attempts)]
    for item in attempts:
        phase = item['phase']
        assert sha(archive / f'{phase}.stdout.txt') == item['stdout_sha256']
        assert sha(archive / f'{phase}.stderr.txt') == item['stderr_sha256']
        result = archive / phase / 'result.json'
        if result.is_file():
            assert sha(result) == item['result_sha256']
            parsed_result = json.loads(result.read_text())
            assert item['reported_decision'] == parsed_result['decision']
            assert item['reported_peak_rss_bytes'] == parsed_result['peak_rss_bytes']
            assert item['reported_cpu_seconds'] == parsed_result['cpu_seconds']
    if receipt['decision'] != 'PASS':
        return {'decision': 'ARCHIVED_FAILURE_ONLY', 'phases': len(attempts)}
    caps = {'toy': (frozen['toy_external_wall_cap_seconds'], frozen['toy_rss_cap_bytes']),
            'panel': (frozen['panel_external_wall_cap_seconds'], frozen['panel_rss_cap_bytes']),
            'n131': (frozen['n131_external_wall_cap_seconds'], frozen['n131_rss_cap_bytes'])}
    assert len(attempts) == 3 and all(item['exit_code'] == 0 and
                                     item['stop_reason'] is None and
                                     item['reported_decision'] == 'PASS' and
                                     item['reported_peak_rss_bytes'] is not None and
                                     item['reported_peak_rss_bytes'] <= caps[item['phase']][1] and
                                     item['sampled_peak_rss_bytes'] <= caps[item['phase']][1] and
                                     item['wall_seconds'] <= caps[item['phase']][0]
                                     for item in attempts)
    sys.path.insert(0, str(HERE))
    from verify import (ReferenceField, check_relation_cnf, evaluate_cnf,
                        lift_solver_model, parse_cnf, parse_solver_output,
                        verify_parent_rows)
    from export import build_relation, fixed_point_units
    toy = json.loads((archive / 'toy/result.json').read_text())
    toy_cnfs = {n: archive / 'toy' / f'n{n}.cnf' for n in (2, 3)}
    replay_toy = verify_parent_rows(toy_cnfs)
    assert [row['sha256'] for row in toy['fields']] == [sha(toy_cnfs[n]) for n in (2, 3)]
    panel = json.loads((archive / 'panel/result.json').read_text())
    expected = _panel_expected()
    assert len(panel['queries']) == 20
    relation = build_relation(13, 0x201B)
    field = ReferenceField(13, 0x201B)
    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)
        checker = temporary / 'drat-trim'
        subprocess.run(['cc', '-O2', '-o', str(checker),
                        str(HERE / 'third_party/drat-trim.c')],
                       check=True, capture_output=True, timeout=60)
        for row in panel['queries']:
            name = row['id']
            assert name in expected and row['expected'] == expected[name]['expected']
            points = {key: tuple(row['points'][key]) for key in ('p', 'q', 'r')}
            assert points == expected[name]['points']
            compressed = archive / 'panel/queries' / f'{name}.cnf.gz'
            assert sha(compressed) == row['cnf_gzip_sha256']
            query = temporary / 'query.cnf'
            query_bytes = gzip.decompress(compressed.read_bytes())
            assert len(query_bytes) == row['cnf']['bytes'] <= frozen['n13_cnf_byte_cap']
            query.write_bytes(query_bytes)
            assert sha(query) == row['cnf']['sha256']
            variables, clauses_count, clauses = parse_cnf(query)
            assert variables == row['cnf']['variables'] and clauses_count == row['cnf']['clauses']
            units = fixed_point_units(relation, points['p'], points['q'], points['r'])
            assert row['units'] == units
            check_relation_cnf(relation, clauses, units)
            solver_stdout = (archive / 'panel/logs' / f'{name}.solver.stdout.txt').read_text()
            status, values = parse_solver_output(solver_stdout, variables,
                                                  row['solver']['exit_code'])
            assert status == row['expected']
            if status == 'SAT':
                assert row['status'] == 'SAT_LIFTED' and values is not None
                assert evaluate_cnf(clauses, values)
                lifted = lift_solver_model(relation, values, points['p'],
                                           points['q'], points['r'], field)
                assert tuple(row['lifted']['p']) == lifted['p']
                assert tuple(row['lifted']['q']) == lifted['q']
                assert tuple(row['lifted']['r']) == lifted['r']
                assert row['lifted']['lambda'] == lifted['lambda']
            else:
                assert row['status'] == 'PROVED_UNSAT' and row['proof_text_valid'] is True
                assert row['error'] is None
                packed = archive / 'panel/proofs' / f'{name}.drat.gz'
                assert sha(packed) == row['proof_gzip_sha256']
                proof = temporary / 'proof.drat'
                proof_bytes = gzip.decompress(packed.read_bytes())
                assert 0 < len(proof_bytes) <= frozen['proof_byte_cap']
                assert b'\x00' not in proof_bytes and all(byte < 128 for byte in proof_bytes)
                proof.write_bytes(proof_bytes)
                checked = subprocess.run([str(checker), str(query), str(proof)],
                                         capture_output=True, timeout=180)
                assert checked.returncode == 0
    n131 = json.loads((archive / 'n131/result.json').read_text())
    assert n131['single_edge_only'] is True and n131['n'] == 131
    assert n131['modulus'] == 2722258935367507707706996859454145699847
    n131_packed = archive / 'n131/single_edge.cnf.gz'
    assert sha(n131_packed) == n131['cnf_gzip_sha256']
    assert n131_packed.stat().st_size == n131['cnf_gzip_bytes']
    with tempfile.TemporaryDirectory() as n131_temporary:
        n131_cnf = Path(n131_temporary) / 'single_edge.cnf'
        raw_n131 = gzip.decompress(n131_packed.read_bytes())
        assert len(raw_n131) == n131['cnf']['bytes'] <= frozen['n131_cnf_byte_cap']
        n131_cnf.write_bytes(raw_n131)
        assert sha(n131_cnf) == n131['cnf']['sha256']
        assert n131['cnf']['dag']['variables'] == 920
        assert n131['cnf']['dag']['model_limbs'] == 15
        variables, clauses_count, _ = parse_cnf(n131_cnf, keep_clauses=False)
        assert (variables, clauses_count) == (n131['cnf']['variables'], n131['cnf']['clauses'])
    return {'decision': 'PASS', 'toy': replay_toy,
            'n13_queries': len(panel['queries']),
            'n131_cnf_sha256': n131['cnf']['sha256']}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path)
    args = parser.parse_args()
    frozen = check_freeze()
    evidence = check_evidence(args.evidence, frozen) if args.evidence else None
    print(json.dumps({'decision': 'PASS', 'freeze_sha256': sha(HERE / 'FROZEN.json'),
                      'evidence': evidence}, sort_keys=True))


if __name__ == '__main__':
    main()
