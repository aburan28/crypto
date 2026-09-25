#!/usr/bin/env python3
"""Fail-closed, all-labelled-tuple n19 semantic admission; no solver calls.

This is intentionally expensive.  A partial/capped run has no admission
status and must be retained as a failure receipt, not counted as coverage.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import resource
import signal
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

from s6 import ARCHIVES, BETAS, O, build, read_factors

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from export import write_cnf  # noqa: E402
from verify import evaluate_cnf, extend_dag, parse_cnf  # noqa: E402
from admission import check_cnf, check_full_model  # noqa: E402

N, MODULUS, Q = 19, 0x80027, 130873
TORSION = (O, (0, 0, 1), (0, 1, 0), (0, 1, 1))
ARCHIVE_SUPPORT = {3: 62389, 338435: 66179, 303097: 66203, 464276: 59323, 42605: 64657}


def point(value):
    return O if value is None else (0, *value)


def fermat_reference():
    """Load the independent #767 replay law, not the Euclid witness law."""
    path = HERE.parent / 'rotated_subspace_support_20260925/verify.py'
    spec = importlib.util.spec_from_file_location('s6_fermat_reference', path)
    if spec is None or spec.loader is None:
        raise AssertionError('missing independent Fermat reference')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.E(module.GF(N, MODULUS))


def bare(p):
    return None if p == O else (p[1], p[2])


def tagged(p):
    return O if p is None else (0, *p)


def negative(p):
    return None if p is None else (p[0], p[0] ^ p[1])


def _archived_histogram(beta: int) -> dict[tuple[int, int, int], int]:
    archive, expected_sha, factors_member = ARCHIVES[beta]
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected_sha:
        raise AssertionError('archive hash drift')
    member = factors_member.replace('factors.json', 'projected_histogram.jsonl')
    with tarfile.open(archive, 'r:gz') as tar:
        entry = tar.getmember(member)
        if not entry.isfile() or entry.size > 40_000_000:
            raise AssertionError('invalid histogram member')
        stream = tar.extractfile(entry)
        if stream is None:
            raise AssertionError('missing histogram')
        hist = {}
        for line in stream:
            row = json.loads(line)
            key = point(row['point'])
            if key in hist or row['count'] < 1:
                raise AssertionError('invalid or duplicate archive row')
            hist[key] = row['count']
    if sum(hist.values()) != 7 ** 6:
        raise AssertionError('archive tuple count drift')
    if beta in ARCHIVE_SUPPORT and len(hist) != ARCHIVE_SUPPORT[beta]:
        raise AssertionError('archive support drift')
    return hist


def run_one(beta: int, out: Path, *, max_seconds: float) -> dict:
    """Every tuple gets exact group, Q/T, DAG, and parsed-CNF validation."""
    factors = read_factors(beta)
    circuit = build(N, MODULUS, factors)
    if len(circuit.dag.names) != 566 or len(circuit.dag.nodes) > 2_000_000:
        raise AssertionError('S6 primary width or DAG node cap exceeded')
    ref = fermat_reference()
    for torsion in TORSION:
        if not ref.on(bare(torsion)) or ref.scalar(bare(torsion), 4) is not None:
            raise AssertionError('frozen torsion point invalid')
    cnf = out / f'beta-{beta}.cnf'
    cnf_meta = write_cnf(circuit, cnf, byte_cap=256 << 20)
    check_cnf(circuit, cnf, [])
    variables, clauses_count, clauses = parse_cnf(cnf)
    if variables != len(circuit.dag.nodes) or clauses_count != cnf_meta['clauses']:
        raise AssertionError('CNF dimensions disagree')
    archived = _archived_histogram(beta)
    observed: Counter = Counter()
    branches: Counter = Counter()
    q_inv4 = pow(4, -1, Q)
    started = time.monotonic()
    rolling = hashlib.sha256()
    completed = 0
    for indices in itertools.product(range(7), repeat=6):
        if time.monotonic() - started > max_seconds:
            raise TimeoutError(f'all-tuple semantic cap after {completed} tuples')
        # An independently calculated sum and the corresponding #802-derived
        # complete Boolean witness must agree.  Neither archive is an oracle.
        total = None
        for i, j in enumerate(indices):
            total = ref.add(total, bare(factors[i][j]))
        projected = ref.scalar(total, 4)
        subgroup_q = ref.scalar(projected, q_inv4)
        torsion_t = ref.add(total, negative(subgroup_q))
        if tagged(torsion_t) not in TORSION or ref.scalar(torsion_t, 4) is not None:
            raise AssertionError('S-Q is not in frozen four-torsion')
        if ref.scalar(subgroup_q, Q) is not None:
            raise AssertionError('Q is outside prime-order subgroup')
        if ref.add(subgroup_q, torsion_t) != total:
            raise AssertionError('full-point Q+T does not equal tuple sum')
        packed, target, path_branches = circuit.witness(indices)
        if target != tagged(total) or not circuit.dag.evaluate(packed, circuit.output):
            raise AssertionError('lifted S6 relation differs from independent sum')
        full = extend_dag(circuit.dag, packed)
        if not evaluate_cnf(clauses, full):
            raise AssertionError('lifted model fails exact parsed CNF')
        lifted = check_full_model(circuit, clauses, full, target)
        if any((full[abs(lit) - 1] == 1) != (lit > 0)
               for lit in circuit.target_units(target)):
            raise AssertionError('exact full-target unit pin fails')
        if lifted['indices'] != list(indices):
            raise AssertionError('selector index drift')
        observed[tagged(projected)] += 1
        branches.update(path_branches)
        # Hash every labelled tuple with its exact (Q,T,S) triple, including O.
        row = [list(indices), list(tagged(subgroup_q)), list(tagged(torsion_t)), list(tagged(total))]
        rolling.update((json.dumps(row, separators=(',', ':')) + '\n').encode())
        completed += 1
        if completed % 1024 == 0:
            (out / 'progress.json').write_text(json.dumps({
                'beta': beta, 'completed': completed, 'rolling_sha256': rolling.hexdigest(),
                'seconds': time.monotonic() - started}, sort_keys=True) + '\n')
    if completed != 7 ** 6 or dict(observed) != archived:
        raise AssertionError('all-tuple count/histogram mismatch')
    if observed[O] < 1 or branches['inverse'] < 1 or branches['copy_q'] < 1:
        raise AssertionError('O and inverse-prefix coverage missing')
    return {'status': 'SEMANTIC_PASS', 'beta': beta, 'tuples': completed,
            'projected_support': len(observed), 'projected_O_multiplicity': observed[O],
            'branch_counts': dict(sorted(branches.items())),
            'tuple_q_t_s_sha256': rolling.hexdigest(), 'cnf': cnf_meta,
            'seconds': time.monotonic() - started}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--beta', required=True, type=int, choices=BETAS)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--internal-supervised', action='store_true')
    args = parser.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    if not args.internal_supervised or not 0 < args.max_seconds <= frozen['semantic_wall_seconds']:
        raise SystemExit('HELD: use bounded run.py with the frozen eight-hour cap')
    # This is a post-merge/re-freeze gate.  No numeric n19 result is produced
    # while the parent #802/#804/#807 PRs remain unmerged.
    if frozen['release_main_head'] is None:
        raise SystemExit('HELD: release_main_head is null; parent merge/re-freeze pending')
    resource.setrlimit(resource.RLIMIT_AS, (frozen['semantic_rss_bytes'], frozen['semantic_rss_bytes']))
    signal.alarm(math.ceil(args.max_seconds))
    if args.out.exists():
        raise SystemExit('refusing to overwrite an earlier evidence directory')
    args.out.mkdir(parents=True)
    try:
        result = run_one(args.beta, args.out, max_seconds=args.max_seconds)
    except Exception as exc:
        result = {'status': 'CENSORED_OR_FAILED', 'beta': args.beta,
                  'error': f'{type(exc).__name__}: {exc}'}
        (args.out / 'result.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        raise
    (args.out / 'result.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'beta', 'tuples', 'projected_support')}, sort_keys=True))


if __name__ == '__main__':
    main()
