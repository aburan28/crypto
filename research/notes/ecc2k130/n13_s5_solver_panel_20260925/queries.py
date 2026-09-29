#!/usr/bin/env python3
"""Pre-outcome exact DIMACS query hashes for all 32 frozen n13 labels."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'n13_s5_fullpoint_control_20260925'))
from query import circuit, write_query  # noqa: E402
from s5 import read_targets  # noqa: E402

MANIFEST = HERE / 'QUERY_MANIFEST.json'
MAX_QUERY_BYTES = 64 * (1 << 20)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_manifest() -> dict:
    c = circuit()
    if len(c.dag.names) != 320 or len(c.dag.nodes) > 2_000_000:
        raise AssertionError('unexpected five-slot circuit dimensions')
    targets = read_targets()
    if len(targets) != 32 or len({row['id'] for row in targets}) != 32:
        raise AssertionError('not the exact 32-label matrix')
    rows = []
    with tempfile.TemporaryDirectory(prefix='s5-query-freeze-') as scratch:
        for index, target in enumerate(targets):
            if target['index'] != index or target['id'] != f'Q{index//4}T{index%4}':
                raise AssertionError('target order drift')
            query = Path(scratch) / f"{target['id']}.cnf"
            meta = write_query(c, target['id'], query, byte_cap=MAX_QUERY_BYTES)
            raw = meta['cnf']
            if raw['bytes'] != query.stat().st_size or raw['sha256'] != sha(query):
                raise AssertionError('query writer/file hash drift')
            if len(meta['units']) != 2 * c.n + 1:
                raise AssertionError('full-target unit count drift')
            rows.append({'index': index, 'label': target['id'],
                         'target': target['full_target'],
                         'archived_exact_tuple_count': target['archived_exact_tuple_count'],
                         'expected_status': 'SAT' if target['archived_exact_tuple_count'] else 'UNSAT',
                         'unit_literals': meta['units'],
                         'variables': raw['variables'], 'clauses': raw['clauses'],
                         'bytes': raw['bytes'], 'sha256': raw['sha256']})
    if sum(row['expected_status'] == 'SAT' for row in rows) != 5:
        raise AssertionError('five-positive/27-negative comparator drift')
    return {'schema': 'ecc2k130_n13_s5_exact_queries_v1',
            'source_input_sha256': sha(HERE.parent / 'n13_s5_fullpoint_control_20260925/INPUT.json'),
            'n': c.n, 'modulus': c.modulus, 'primary_bits': len(c.dag.names),
            'dag_nodes': len(c.dag.nodes), 'labels': rows}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--write', action='store_true')
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    if args.write == args.verify:
        raise SystemExit('choose exactly one of --write or --verify')
    actual = make_manifest()
    if args.write:
        if MANIFEST.exists():
            raise FileExistsError('query manifest already exists')
        MANIFEST.write_text(json.dumps(actual, sort_keys=True, indent=2) + '\n')
    else:
        if actual != json.loads(MANIFEST.read_text()):
            raise AssertionError('frozen exact-query bytes changed')
    print(json.dumps({'decision': 'QUERY_BYTES_MATCH' if args.verify else 'QUERY_BYTES_FROZEN',
                      'labels': len(actual['labels']),
                      'positives': sum(row['expected_status'] == 'SAT'
                                       for row in actual['labels']),
                      'dag_nodes': actual['dag_nodes'],
                      'first_query_sha256': actual['labels'][0]['sha256']}, sort_keys=True))


if __name__ == '__main__':
    main()
