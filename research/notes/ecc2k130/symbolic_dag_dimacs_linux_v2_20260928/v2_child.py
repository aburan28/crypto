#!/usr/bin/env python3
"""Linux host correction for the unchanged #804 phase functions.

The #804 ``produce.main`` entrypoint has a masked release-gate arity error.
This wrapper imports its frozen semantic functions without calling that
entrypoint.  The only per-phase configuration change is the pinned Linux
CaDiCaL executable path.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from ci_replay import check_static
from run import release_gate

HERE = Path(__file__).resolve().parent
FIRST = HERE.parent / 'symbolic_dag_dimacs_gate_20260925'


def inherited_producer():
    # Deliberately resolve #804's bounded/export/verify modules from its own
    # frozen directory.  Importing produce does not execute produce.main().
    sys.path.insert(0, str(FIRST))
    import produce  # type: ignore[import-not-found]
    return produce


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--gate-only', action='store_true')
    parser.add_argument('--phase', choices=('toy', 'panel', 'n131'))
    parser.add_argument('--out', type=Path)
    parser.add_argument('--expected-head')
    args = parser.parse_args()
    frozen = check_static()
    producer = inherited_producer()
    if args.gate_only:
        if args.phase or args.out:
            parser.error('--gate-only cannot select a measured phase')
        assert all(callable(getattr(producer, name)) for name in ('toy', 'panel', 'n131'))
        print(json.dumps({'decision': 'HASH_ONLY_NO_MEASURED_CHILD',
                          'release_status': frozen['status'],
                          'inherited_functions': ['toy', 'panel', 'n131']},
                         sort_keys=True))
        return 0
    if not args.phase or not args.out or not args.expected_head:
        parser.error('measured child requires --phase, --out and --expected-head')
    release_gate(frozen, args.expected_head)
    original = json.loads((FIRST / 'FROZEN.json').read_text())
    original['cadical_path'] = str((HERE.parent / frozen['linux_binary_path']).resolve())
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    start_wall, start_cpu = time.monotonic(), time.process_time()
    result = {'toy': producer.toy, 'panel': producer.panel,
              'n131': producer.n131}[args.phase](out, original)
    result.update({'phase': args.phase, 'domain': original['domain'],
                   'v2_wrapper': frozen['domain'],
                   'wall_seconds': time.monotonic() - start_wall,
                   'cpu_seconds': time.process_time() - start_cpu,
                   'peak_rss_bytes': producer._rss_bytes()})
    (out / 'result.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'phase': args.phase, 'decision': result['decision'],
                      'wall_seconds': result['wall_seconds']}, sort_keys=True))
    return 0 if result['decision'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
