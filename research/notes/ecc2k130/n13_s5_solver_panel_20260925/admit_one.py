#!/usr/bin/env python3
"""Bounded child: admit exact SAT model or externally checked text DRAT."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'n13_s5_fullpoint_control_20260925'))
from admission5 import admit_sat, check_cnf  # noqa: E402
from query import circuit  # noqa: E402
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from bounded import run_child  # noqa: E402
from verify import parse_cnf, parse_solver_output  # noqa: E402


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def admit_unsat_receipted(c, row: dict, query: Path, proof: Path,
                          solver_stdout: Path, solver_exit: int,
                          checker: Path, checker_sha256: str,
                          out: Path, frozen: dict) -> dict:
    """Retain the exact external checker receipt, including failed checks."""
    target = tuple(row['target'])
    check_cnf(c, query, c.target_units(target))
    variables, _, _ = parse_cnf(query, keep_clauses=False)
    status, values = parse_solver_output(solver_stdout.read_text(), variables,
                                         solver_exit)
    if status != 'UNSAT' or values is not None:
        raise AssertionError('solver did not emit exact UNSAT status and exit 20')
    if not proof.is_file() or not 0 < proof.stat().st_size <= frozen['proof_byte_max']:
        raise AssertionError('missing, empty, or over-cap proof')
    with proof.open('rb') as stream:
        while chunk := stream.read(1 << 20):
            if b'\x00' in chunk or any(byte > 127 for byte in chunk):
                raise AssertionError('proof is not text DRAT')
    if not checker.is_file() or sha(checker) != checker_sha256:
        raise AssertionError('external DRAT checker binary SHA-256 mismatch')
    logs = out / 'checker_logs'
    logs.mkdir(parents=True, exist_ok=False)
    receipt = run_child([str(checker.resolve()), str(query.resolve()),
                         str(proof.resolve())], cwd=HERE,
                        stdout=logs / 'checker.stdout.txt',
                        stderr=logs / 'checker.stderr.txt',
                        wall_cap=frozen['checker_wall_seconds'],
                        rss_cap=frozen['checker_rss_bytes'],
                        file_cap=frozen['checker_log_byte_max'])
    (out / 'checker_receipt.json').write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    if receipt['exit_code'] != 0 or receipt['stop_reason'] is not None:
        raise AssertionError('external DRAT checker did not succeed')
    return {'status': 'PROVED_UNSAT', 'label': row['label'],
            'query_sha256': sha(query), 'proof_sha256': sha(proof),
            'checker_binary_sha256': checker_sha256,
            'checker': receipt}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--label', required=True)
    p.add_argument('--kind', choices=('SAT', 'UNSAT'), required=True)
    p.add_argument('--query', type=Path, required=True)
    p.add_argument('--solver-stdout', type=Path, required=True)
    p.add_argument('--solver-exit', type=int, required=True)
    p.add_argument('--proof', type=Path)
    p.add_argument('--checker', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--internal-supervised', action='store_true')
    args = p.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    if not args.internal_supervised or frozen['release_main_head'] is None:
        raise SystemExit('HELD: merged-parent release and bounded supervisor required')
    rows = json.loads((HERE / 'QUERY_MANIFEST.json').read_text())['labels']
    matches = [row for row in rows if row['label'] == args.label]
    if len(matches) != 1 or sha(args.query) != matches[0]['sha256']:
        raise SystemExit('exact query does not match frozen label/hash')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    try:
        c = circuit()
        if args.kind == 'SAT':
            verified = admit_sat(c, args.label, args.query,
                                 args.solver_stdout, args.solver_exit)
        else:
            if args.proof is None or args.checker is None:
                raise ValueError('UNSAT admission needs proof and checker')
            verified = admit_unsat_receipted(c, matches[0], args.query,
                                              args.proof, args.solver_stdout,
                                              args.solver_exit, args.checker,
                                              frozen['checker_binary_sha256'],
                                              out, frozen)
        result = {'status': verified['status'], 'label': args.label,
                  'query_sha256': sha(args.query),
                  'solver_stdout_sha256': sha(args.solver_stdout),
                  'verified': verified}
    except Exception as exc:
        result = {'status': 'ADMISSION_FAILURE', 'label': args.label,
                  'query_sha256': sha(args.query),
                  'solver_stdout_sha256': sha(args.solver_stdout),
                  'error': f'{type(exc).__name__}: {exc}',
                  'proof_sha256': sha(args.proof) if args.proof and args.proof.is_file() else None}
        (out / 'admission.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        raise
    (out / 'admission.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'label': args.label}, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
