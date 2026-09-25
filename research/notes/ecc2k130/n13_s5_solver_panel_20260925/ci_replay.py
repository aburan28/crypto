#!/usr/bin/env python3
"""Hash/syntax/exact-query preflight only; no SAT/proof calls."""
from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

from queries import make_manifest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if sys.version_info < (3, 12):
        raise SystemExit('Python 3.12+ required by frozen panel')
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    assert frozen['schema'] == 'ecc2k130_n13_s5_solver_panel_v1'
    assert frozen['required_merged_prs'] == [802, 804, 807, 811, 813]
    assert (frozen['export_wall_seconds'], frozen['export_rss_bytes']) == (60, 2 << 30)
    assert (frozen['solver_wall_seconds'], frozen['solver_rss_bytes']) == (120, 4 << 30)
    assert (frozen['checker_wall_seconds'], frozen['checker_rss_bytes']) == (300, 4 << 30)
    assert (frozen['sat_lift_wall_seconds'], frozen['sat_lift_rss_bytes']) == (60, 1 << 30)
    assert (frozen['total_panel_wall_seconds'], frozen['total_panel_rss_bytes']) == (28800, 16 << 30)
    assert frozen['query_cnf_byte_max'] == 64 << 20
    assert frozen['proof_byte_max'] == 1 << 30
    assert frozen['checker_log_byte_max'] == 64 << 20
    assert frozen['panel_lock_path'] == '/private/tmp/ecc2k130-s5-panel.lock'
    for relative, expected in frozen['source_sha256'].items():
        path = HERE / relative
        assert path.is_file() and sha(path) == expected, relative
        if path.suffix == '.py':
            ast.parse(path.read_text(), filename=str(path))
    for category in ('merged_parent_sha256', 'workflow_sha256'):
        for relative, expected in frozen[category].items():
            path = ROOT / relative
            assert path.is_file() and sha(path) == expected, relative
    checker_source = ROOT / frozen['checker_source_path']
    assert sha(checker_source) == frozen['checker_source_sha256']
    # Regenerate *all* query bytes with the frozen complete full-point exporter.
    # No solver or checker process is called here or in the hosted workflow.
    actual = make_manifest()
    expected = json.loads((HERE / 'QUERY_MANIFEST.json').read_text())
    assert actual == expected
    labels = expected['labels']
    assert len(labels) == 32 and [r['label'] for r in labels] == [f'Q{i//4}T{i%4}' for i in range(32)]
    assert len({r['sha256'] for r in labels}) == 32
    assert sum(r['expected_status'] == 'SAT' for r in labels) == 5
    assert sum(r['expected_status'] == 'UNSAT' for r in labels) == 27
    assert all((r['expected_status'] == 'SAT') == (r['archived_exact_tuple_count'] > 0)
               for r in labels)
    assert labels[12]['label'] == 'Q3T0' and labels[12]['target'] == [0, 7256, 3272]
    assert expected['dag_nodes'] == 19219 and expected['primary_bits'] == 320
    print(json.dumps({'decision': 'HASH_QUERY_PREFLIGHT_PASS',
                      'labels': len(labels), 'sat_expected': 5,
                      'unsat_expected': 27, 'n13_solver_outcome': 'UNRUN'}, sort_keys=True))


if __name__ == '__main__':
    main()
