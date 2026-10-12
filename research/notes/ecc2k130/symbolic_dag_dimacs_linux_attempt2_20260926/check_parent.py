#!/usr/bin/env python3
"""Verify attempt-1 archive and inherited semantics before Linux preparation."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FIRST = HERE.parent / 'symbolic_dag_dimacs_gate_20260925'
PARENT = HERE.parent / 'symbolic_dag_fullpoint_20260925'
ARCHIVED_SAT = HERE.parent / 'n13_oaware_sat_benchmark_20260925/evidence/panel/Q0T3-cadical.stdout'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    prepare = json.loads((HERE / 'PREPARE.json').read_text())
    if subprocess.run(['git', 'merge-base', '--is-ancestor',
                       prepare['attempt1_record_head'], 'HEAD'], cwd=ROOT,
                      capture_output=True).returncode != 0:
        raise SystemExit('NOT_ADMITTED: first-failure record is not an ancestor')
    pinned = {
        'FROZEN.json': prepare['attempt1_freeze_sha256'],
        'evidence/receipt.json': prepare['attempt1_receipt_sha256'],
        'evidence/MANIFEST.json': prepare['attempt1_manifest_sha256'],
        'INPUT.json': prepare['input_sha256'],
        **prepare['semantic_source_sha256'],
    }
    for path, expected in ((PARENT / 'FROZEN.json', prepare['parent_freeze_sha256']),
                           (PARENT / 'evidence/producer/rows.jsonl.gz',
                            prepare['parent_rows_sha256']),
                           (ARCHIVED_SAT,
                            prepare['archived_cadical_sat_output_sha256'])):
        if sha(path) != expected:
            raise SystemExit(f'NOT_ADMITTED: inherited {path.name} drifted')
    for name, expected in pinned.items():
        if sha(FIRST / name) != expected:
            raise SystemExit(f'NOT_ADMITTED: first-attempt {name} drifted')
    receipt = json.loads((FIRST / 'evidence/receipt.json').read_text())
    if (receipt['freeze_sha256'] != prepare['attempt1_freeze_sha256'] or
            receipt['release_gate']['checkout_head'] !=
            prepare['attempt1_runner_head'] or
            receipt['decision'] != 'FAIL_OR_CENSORED' or
            len(receipt['attempts']) != 1 or
            receipt['attempts'][0]['stop_reason'] != 'LAUNCH_ERROR'):
        raise SystemExit('NOT_ADMITTED: first-attempt classification changed')
    print(json.dumps({'decision': 'PASS',
                      'first_attempt_receipt_sha256':
                      prepare['attempt1_receipt_sha256'],
                      'semantic_sources': len(prepare['semantic_source_sha256'])},
                     sort_keys=True))


if __name__ == '__main__':
    main()
