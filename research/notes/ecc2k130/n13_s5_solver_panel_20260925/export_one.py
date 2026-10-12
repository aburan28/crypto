#!/usr/bin/env python3
"""Bounded child: write one exact S5 DIMACS query, without a solver call."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'n13_s5_fullpoint_control_20260925'))
from query import circuit, write_query  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--label', required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--internal-supervised', action='store_true')
    args = p.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    if not args.internal_supervised or frozen['release_main_head'] is None:
        raise SystemExit('HELD: merged-parent release and bounded supervisor required')
    rows = json.loads((HERE / 'QUERY_MANIFEST.json').read_text())['labels']
    matches = [row for row in rows if row['label'] == args.label]
    if len(matches) != 1:
        raise SystemExit('unknown or duplicate frozen label')
    row = matches[0]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    try:
        meta = write_query(circuit(), args.label, out / 'query.cnf',
                           byte_cap=frozen['query_cnf_byte_max'])
        for key in ('sha256', 'bytes', 'variables', 'clauses'):
            if meta['cnf'][key] != row[key]:
                raise AssertionError(f'exact query {key} drift')
        if meta['units'] != row['unit_literals'] or meta['target'] != row['target']:
            raise AssertionError('full-target literal/point drift')
        result = {'status': 'EXACT_QUERY', 'label': args.label,
                  'elapsed_seconds': time.monotonic() - started,
                  'query': meta}
    except Exception as exc:
        result = {'status': 'EXPORT_FAILURE', 'label': args.label,
                  'elapsed_seconds': time.monotonic() - started,
                  'error': f'{type(exc).__name__}: {exc}',
                  'partial_query_bytes': (out / 'query.cnf').stat().st_size
                  if (out / 'query.cnf').is_file() else None}
        (out / 'query_meta.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        raise
    (out / 'query_meta.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'label': args.label}, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
