#!/usr/bin/env python3
"""Hash/syntax-only S6 preflight; never runs n19 or invokes a solver."""
from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if sys.version_info < (3, 10):
        raise SystemExit('Python >=3.10 required by parent zip(strict=True)')
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    assert frozen['schema'] == 'ecc2k130_n19_s6_fullpoint_exporter_v1'
    assert frozen['beta_order'] == [3, 338435, 303097, 464276, 42605]
    assert frozen['required_merged_prs'] == [802, 804, 807]
    assert frozen['semantic_wall_seconds'] == 28800
    assert frozen['semantic_rss_bytes'] == 16 * (1 << 30)
    assert frozen['cnf_byte_max'] == 256 * (1 << 20)
    for relative, expected in frozen['source_sha256'].items():
        path = HERE / relative
        assert path.is_file() and sha(path) == expected, relative
        if path.suffix == '.py':
            ast.parse(path.read_text(), filename=str(path))
    for relative, expected in frozen['merged_parent_sha256'].items():
        path = ROOT / relative
        assert path.is_file() and sha(path) == expected, relative
    for relative, expected in frozen['workflow_sha256'].items():
        path = ROOT / relative
        assert path.is_file() and sha(path) == expected, relative
    for relative, expected in frozen['archive_sha256'].items():
        path = ROOT / relative
        assert path.is_file() and sha(path) == expected, relative
    toy = json.loads((HERE / 'evidence/toy_selftest.json').read_text())
    assert toy['decision'] == 'TOY_SELFTEST_PASS'
    assert toy['cnf']['variables'] == 2662 and toy['cnf']['clauses'] == 9018
    print(json.dumps({'decision': 'HASH_SYNTAX_PASS',
                      'source_files': len(frozen['source_sha256']),
                      'parent_files': len(frozen['merged_parent_sha256']),
                      'n19_outcome': 'UNRUN'}, sort_keys=True))


if __name__ == '__main__':
    main()
