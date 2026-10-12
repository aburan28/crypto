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
import os
import re
import sys
import time
from pathlib import Path

from ci_replay import check_child_bytes
from run import git, require, sha

HERE = Path(__file__).resolve().parent
FIRST = HERE.parent / 'symbolic_dag_dimacs_gate_20260925'


def inherited_producer():
    # Deliberately resolve #804's bounded/export/verify modules from its own
    # frozen directory.  Importing produce does not execute produce.main().
    sys.path.insert(0, str(FIRST))
    import produce  # type: ignore[import-not-found]
    return produce


def dispatch_gate(frozen: dict, expected_head: str, out: Path,
                  dispatch_path: Path, dispatch_sha256: str) -> dict:
    """Consume only a byte-checked dispatch already issued by the outer gate."""
    require(frozen['status'] == 'RELEASED', 'v2 remains held')
    require(re.fullmatch(r'[0-9a-f]{40}', expected_head or '') is not None,
            'missing reviewed exact head')
    require(re.fullmatch(r'[0-9a-f]{64}', dispatch_sha256 or '') is not None,
            'missing dispatch digest')
    dispatch_path = dispatch_path.resolve()
    out = out.resolve()
    require(dispatch_path == out.parent / 'DISPATCH.json' and
            sha(dispatch_path) == dispatch_sha256,
            'dispatch file path or digest differs from supervisor')
    gate = json.loads(dispatch_path.read_text())
    run_id = int(os.environ.get('GITHUB_RUN_ID', '0'))
    require(os.environ.get('GITHUB_ACTIONS') == 'true' and run_id > 0 and
            sys.platform.startswith('linux') and os.uname().machine == 'x86_64' and
            sys.version_info[:2] == (3, 12),
            'wrong measured child host or event')
    require(git('rev-parse', 'HEAD') == expected_head and
            gate['reviewed_head'] == gate['checkout_head'] ==
            gate['pr_head'] == expected_head and
            gate['freeze_sha256'] == sha(HERE / 'FROZEN.json') and
            gate['release_main_head'] == frozen['release_main_head'] and
            gate['pr_number'] == frozen['release_pr_number'] and
            gate['preparation_freeze_sha256'] ==
            frozen['preparation']['freeze_sha256'] and
            gate['linux_binary_sha256'] == frozen['linux_binary_sha256'] and
            gate['event']['run_id'] == run_id and
            gate['event']['event_head'] == expected_head and
            gate['one_shot']['run_id'] == run_id and
            gate['one_shot']['matching_labeled_runs'] == [run_id] and
            sha(dispatch_path.parent / 'target_rlimit_probe.json') ==
            gate['target_cap_probe']['receipt_sha256'],
            'supervisor dispatch provenance drift')
    return gate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--gate-only', action='store_true')
    parser.add_argument('--phase', choices=('toy', 'panel', 'n131'))
    parser.add_argument('--out', type=Path)
    parser.add_argument('--expected-head')
    parser.add_argument('--dispatch', type=Path)
    parser.add_argument('--dispatch-sha256')
    args = parser.parse_args()
    frozen = check_child_bytes()
    producer = inherited_producer()
    if args.gate_only:
        if args.phase or args.out or args.dispatch or args.dispatch_sha256:
            parser.error('--gate-only cannot select a measured phase')
        assert all(callable(getattr(producer, name)) for name in ('toy', 'panel', 'n131'))
        print(json.dumps({'decision': 'HASH_ONLY_NO_MEASURED_CHILD',
                          'release_status': frozen['status'],
                          'inherited_functions': ['toy', 'panel', 'n131']},
                         sort_keys=True))
        return 0
    if (not args.phase or not args.out or not args.expected_head or
            not args.dispatch or not args.dispatch_sha256):
        parser.error('measured child requires phase, out, head and sealed dispatch')
    dispatch_gate(frozen, args.expected_head, args.out, args.dispatch,
                  args.dispatch_sha256)
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
