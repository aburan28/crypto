#!/usr/bin/env python3
"""Frozen toy equivalence, n13 local-gate panel, and n131 export smoke."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

from bounded import run_child, sha
from export import build_relation, fixed_point_units, write_cnf
from verify import ReferenceField, evaluate_cnf, lift_solver_model, parse_cnf, parse_solver_output

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'symbolic_dag_fullpoint_20260925'
INPUT = json.loads((HERE / 'INPUT.json').read_text())
N13_POLY = 0x201B


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == 'darwin' else value * 1024


def _text_proof_valid(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    with path.open('rb') as stream:
        while chunk := stream.read(1 << 20):
            if b'\x00' in chunk or any(byte > 127 for byte in chunk):
                return False
    return True


def _gzip_copy(source: Path, dest: Path) -> None:
    with source.open('rb') as src, dest.open('wb') as dst:
        with gzip.GzipFile(filename='', mode='wb', fileobj=dst, mtime=0) as zipped:
            while chunk := src.read(1 << 20):
                zipped.write(chunk)


def toy(out: Path, frozen: dict) -> dict:
    rows = []
    for n, modulus in ((2, 0x7), (3, 0xB)):
        relation = build_relation(n, modulus)
        cnf = out / f'n{n}.cnf'
        meta = write_cnf(relation, cnf, byte_cap=frozen['toy_cnf_byte_cap'])
        rows.append({'n': n, 'modulus': modulus, **meta})
    return {'fields': rows, 'decision': 'PASS'}


def _panel_cases(field: ReferenceField) -> list[dict]:
    assert (INPUT['field']['n'], INPUT['field']['modulus']) == (13, N13_POLY)
    cases = []
    for row in INPUT['pairs']:
        p, q, correct, wrong = (tuple(row[key]) for key in ('p', 'q', 'sat_r', 'unsat_r'))
        assert all(field.point_valid(point) for point in (p, q, correct, wrong))
        assert field.add(p, q) == correct and correct != wrong
        cases.append({'id': row['id'] + '_sat', 'expected': 'SAT', 'p': p, 'q': q, 'r': correct})
        cases.append({'id': row['id'] + '_wrong', 'expected': 'UNSAT', 'p': p, 'q': q, 'r': wrong})
    for row in INPUT['invalid']:
        p, q, r = (tuple(row[key]) for key in ('p', 'q', 'r'))
        assert not all(field.point_valid(point) for point in (p, q, r))
        cases.append({'id': row['id'], 'expected': 'UNSAT', 'p': p, 'q': q, 'r': r})
    assert len(cases) == 20 and len(set(row['id'] for row in cases)) == 20
    return cases


def panel(out: Path, frozen: dict) -> dict:
    field = ReferenceField(13, N13_POLY)
    relation = build_relation(13, N13_POLY)
    base_meta = write_cnf(relation, out / 'base.cnf', byte_cap=frozen['n13_cnf_byte_cap'])
    checker = out / 'drat-trim'
    compile_receipt = run_child(['cc', '-O2', '-o', str(checker),
                                 str(HERE / 'third_party/drat-trim.c')],
                                cwd=HERE, stdout=out / 'checker_compile.stdout.txt',
                                stderr=out / 'checker_compile.stderr.txt',
                                wall_cap=frozen['checker_compile_wall_cap_seconds'],
                                rss_cap=frozen['n13_rss_cap_bytes'])
    if compile_receipt['exit_code'] != 0 or compile_receipt['stop_reason']:
        return {'decision': 'CHECKER_BUILD_FAILURE', 'base': base_meta,
                'checker_compile': compile_receipt, 'queries': []}
    compile_receipt['binary_sha256'] = sha(checker)
    queries = []
    (out / 'queries').mkdir()
    (out / 'proofs').mkdir()
    (out / 'logs').mkdir()
    for case in _panel_cases(field):
        name = case['id']
        query = out / 'query.cnf'
        units = fixed_point_units(relation, case['p'], case['q'], case['r'])
        meta = write_cnf(relation, query, units=units,
                         byte_cap=frozen['n13_cnf_byte_cap'])
        variables, _, clauses = parse_cnf(query)
        proof = out / 'proof.drat'
        proof.unlink(missing_ok=True)
        command = [frozen['cadical_path'], '--no-binary', str(query)]
        if case['expected'] == 'UNSAT':
            command.append(str(proof))
        solver = run_child(command, cwd=HERE,
                           stdout=out / 'logs' / f'{name}.solver.stdout.txt',
                           stderr=out / 'logs' / f'{name}.solver.stderr.txt',
                           wall_cap=frozen['solver_wall_cap_seconds'],
                           rss_cap=frozen['n13_rss_cap_bytes'],
                           watched_file=proof if case['expected'] == 'UNSAT' else None,
                           watched_file_cap=frozen['proof_byte_cap'],
                           file_cap=frozen['proof_byte_cap'] if case['expected'] == 'UNSAT' else None)
        parsed = None
        lifted = None
        error = None
        if solver['stop_reason'] is None and solver['exit_code'] in (10, 20):
            try:
                parsed, values = parse_solver_output(
                    (out / 'logs' / f'{name}.solver.stdout.txt').read_text(),
                    variables, solver['exit_code'])
                if parsed == 'SAT':
                    assert values is not None and evaluate_cnf(clauses, values)
                    lifted = lift_solver_model(relation, values,
                                               case['p'], case['q'], case['r'], field)
            except (AssertionError, ValueError) as exc:
                error = f'{type(exc).__name__}: {exc}'
        else:
            error = 'solver cap, error, or unexpected exit'
        proof_valid = _text_proof_valid(proof) and proof.stat().st_size <= frozen['proof_byte_cap']
        checker_receipt = None
        if parsed == 'UNSAT' and proof_valid:
            checker_receipt = run_child([str(checker), str(query), str(proof)],
                                        cwd=HERE,
                                        stdout=out / 'logs' / f'{name}.checker.stdout.txt',
                                        stderr=out / 'logs' / f'{name}.checker.stderr.txt',
                                        wall_cap=frozen['checker_wall_cap_seconds'],
                                        rss_cap=frozen['n13_rss_cap_bytes'])
        archived_query = out / 'queries' / f'{name}.cnf.gz'
        _gzip_copy(query, archived_query)
        archived_proof = None
        if proof.is_file():
            archived_proof = out / 'proofs' / f'{name}.drat.gz'
            _gzip_copy(proof, archived_proof)
        query.unlink()
        proof.unlink(missing_ok=True)
        if parsed == 'UNSAT' and not proof_valid:
            error = 'missing, empty, non-ASCII, or oversized DRAT proof'
        verified_unsat = bool(checker_receipt and checker_receipt['exit_code'] == 0
                              and checker_receipt['stop_reason'] is None)
        status = ('SAT_LIFTED' if parsed == 'SAT' and lifted is not None and error is None else
                  'PROVED_UNSAT' if parsed == 'UNSAT' and verified_unsat and proof_valid and error is None else
                  'UNPROVEN_UNSAT' if parsed == 'UNSAT' else 'FAIL')
        query_row = {'id': name, 'expected': case['expected'], 'observed': parsed,
                     'status': status, 'points': {k: case[k] for k in ('p', 'q', 'r')},
                     'units': units, 'cnf': meta, 'cnf_gzip_sha256': sha(archived_query),
                     'proof_gzip_sha256': sha(archived_proof) if archived_proof else None,
                     'proof_text_valid': proof_valid,
                     'solver': solver, 'checker': checker_receipt,
                     'lifted': lifted, 'error': error}
        queries.append(query_row)
        (out / 'progress.json').write_text(json.dumps({'completed': len(queries),
                                                        'last': query_row}, sort_keys=True, indent=2) + '\n')
        if (case['expected'] == 'SAT' and status != 'SAT_LIFTED') or (case['expected'] == 'UNSAT' and status != 'PROVED_UNSAT'):
            break
    checker.unlink(missing_ok=True)
    result = {'decision': 'PASS' if len(queries) == 20 and all(
                  q['status'] == ('SAT_LIFTED' if q['expected'] == 'SAT' else 'PROVED_UNSAT')
                  for q in queries) else 'FAIL_OR_CENSORED',
              'base': base_meta, 'checker_compile': compile_receipt,
              'queries': queries}
    return result


def n131(out: Path, frozen: dict) -> dict:
    n = 131
    modulus = INPUT['n131_single_edge_modulus']
    relation = build_relation(n, modulus)
    if len(relation.dag.names) != 920 or relation.dag.counts()['model_limbs'] != 15:
        raise AssertionError('n131 input width changed')
    if len(relation.dag.nodes) > frozen['n131_node_cap']:
        return {'decision': 'NODE_CAP', 'nodes': len(relation.dag.nodes),
                'dag': relation.dag.counts(), 'cnf': None}
    try:
        meta = write_cnf(relation, out / 'single_edge.cnf',
                         byte_cap=frozen['n131_cnf_byte_cap'])
    except Exception as exc:
        return {'decision': 'EXPORT_FAILURE', 'error': f'{type(exc).__name__}: {exc}',
                'nodes': len(relation.dag.nodes), 'dag': relation.dag.counts(),
                'partial_bytes': (out / 'single_edge.cnf').stat().st_size
                if (out / 'single_edge.cnf').is_file() else None}
    archived = out / 'single_edge.cnf.gz'
    _gzip_copy(out / 'single_edge.cnf', archived)
    (out / 'single_edge.cnf').unlink()
    return {'decision': 'PASS', 'n': n, 'modulus': modulus,
            'single_edge_only': True, 'cnf': meta,
            'cnf_gzip_sha256': sha(archived),
            'cnf_gzip_bytes': archived.stat().st_size}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=('toy', 'panel', 'n131'), required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    subprocess.run([sys.executable, str(HERE / 'ci_replay.py')], check=True)
    from run import release_gate
    release_gate(frozen)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    start, cpu = time.monotonic(), time.process_time()
    result = {'toy': toy, 'panel': panel, 'n131': n131}[args.phase](out, frozen)
    result.update({'phase': args.phase, 'domain': frozen['domain'],
                   'wall_seconds': time.monotonic() - start,
                   'cpu_seconds': time.process_time() - cpu,
                   'peak_rss_bytes': _rss_bytes()})
    (out / 'result.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'phase': args.phase, 'decision': result['decision'],
                      'wall_seconds': result['wall_seconds']}, sort_keys=True))
    return 0 if result['decision'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
