#!/usr/bin/env python3
"""Held exhaustive five-slot n13 semantic control; no SAT solver."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path

from admission5 import decode, fermat_reference
from query import circuit
from s5 import read_targets

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from export import write_cnf  # noqa: E402
from verify import evaluate_cnf, extend_dag, parse_cnf  # noqa: E402
sys.path.insert(0, str(HERE.parent / 'n19_s6_fullpoint_exporter_20260925'))
from admission import check_cnf  # noqa: E402


def tagged(p):
    return (1, 0, 0) if p is None else (0, *p)


def run(out: Path, frozen: dict) -> dict:
    c = circuit()
    if len(c.dag.names) != 320 or len(c.dag.nodes) > frozen['dag_node_max']:
        raise AssertionError('n13 primary width or DAG node cap drift')
    cnf = out / 'base.cnf'
    meta = write_cnf(c, cnf, byte_cap=frozen['cnf_byte_max'])
    check_cnf(c, cnf, [])
    _, _, clauses = parse_cnf(cnf)
    ref = fermat_reference(c.n, c.modulus)
    targets = read_targets()
    if len(targets) != 32 or len({tuple(row['full_target']) for row in targets}) != 32:
        raise AssertionError('32 distinct fixed full-target labels required')
    by_point = {tuple(row['full_target']): row for row in targets}
    counts: Counter = Counter()
    branches: Counter = Counter()
    rolling = hashlib.sha256()
    exceptional = None
    started = time.monotonic()
    completed = 0
    for indices in itertools.product(range(5), repeat=5):
        if time.monotonic() - started > frozen['semantic_wall_seconds']:
            raise TimeoutError(f'n13 semantic gate capped after {completed} tuples')
        total = None
        for i, j in enumerate(indices):
            p = c.factors[i][j]
            total = ref.add(total, (p[1], p[2]))
        packed, model_sum, path = c.witness(indices)
        if model_sum != tagged(total) or not c.dag.evaluate(packed, c.output):
            raise AssertionError('S5 model differs from independent Fermat sum')
        full = extend_dag(c.dag, packed)
        if not evaluate_cnf(clauses, full):
            raise AssertionError('complete S5 witness fails parsed DIMACS')
        lifted = decode(c, full, model_sum)
        if lifted['indices'] != list(indices):
            raise AssertionError('selector model drift')
        branches.update(path)
        label = by_point.get(model_sum)
        if label is not None:
            if any((full[abs(lit) - 1] == 1) != (lit > 0)
                   for lit in c.target_units(model_sum)):
                raise AssertionError('full-target unit pin failed')
            counts[label['id']] += 1
            if indices == (0, 0, 0, 4, 1):
                if label['id'] != 'Q3T0' or lifted['prefixes'][0] != [1, 0, 0]:
                    raise AssertionError('exceptional O-prefix witness drift')
                exceptional = {'label': label['id'], 'indices': list(indices),
                               'full_target': list(model_sum),
                               'first_prefix': lifted['prefixes'][0],
                               'branches': lifted['branches']}
        rolling.update((json.dumps([list(indices), list(model_sum)],
                                   separators=(',', ':')) + '\n').encode())
        completed += 1
        if completed % 64 == 0:
            (out / 'progress.json').write_text(json.dumps({
                'completed': completed, 'rolling_sha256': rolling.hexdigest(),
                'wall_seconds': time.monotonic() - started}, sort_keys=True) + '\n')
    if completed != 5 ** 5 or exceptional is None:
        raise AssertionError('incomplete tuple census or exceptional witness')
    for target in targets:
        if counts[target['id']] != target['archived_exact_tuple_count']:
            raise AssertionError(f"archived exact count mismatch for {target['id']}")
    return {'status': 'SEMANTIC_PASS', 'tuples': completed,
            'target_counts': {row['id']: counts[row['id']] for row in targets},
            'branch_counts': dict(sorted(branches.items())),
            'exceptional': exceptional, 'tuple_sum_sha256': rolling.hexdigest(),
            'cnf': meta, 'wall_seconds': time.monotonic() - started}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--internal-supervised', action='store_true')
    args = parser.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    if not args.internal_supervised or frozen['release_main_head'] is None:
        raise SystemExit('HELD: parent merge, exact-head re-freeze and bounded supervisor required')
    resource.setrlimit(resource.RLIMIT_AS,
                       (frozen['semantic_rss_bytes'], frozen['semantic_rss_bytes']))
    signal.alarm(math.ceil(frozen['semantic_wall_seconds']))
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    try:
        result = run(out, frozen)
    except Exception as exc:
        result = {'status': 'CENSORED_OR_FAILED',
                  'error': f'{type(exc).__name__}: {exc}'}
        (out / 'result.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        raise
    (out / 'result.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'tuples': result['tuples']}, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
