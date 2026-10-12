#!/usr/bin/env python3
"""Check the frozen source and fixed n13 trace-contract input."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def main() -> None:
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    assert frozen['schema'] == 'compact-orbit-rank-evidence-freeze-v1'
    assert frozen['control'] == {
        'n': 13, 'a': 0, 'columns': 2, 'rank_seed': 7,
        'target_scalar': 7, 'scope': 'nonmeasuring_trace_contract',
    }
    assert (HERE / 'targets_n13.txt').read_text() == '7\n'
    for relative, expected in frozen['sha256'].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected, relative
    subprocess.run(['git', 'merge-base', '--is-ancestor', frozen['source_base_main'], 'HEAD'],
                   cwd=ROOT, check=True)
    print(json.dumps({'status': 'PASS', 'source_base_main': frozen['source_base_main'],
                      'freeze_sha256': hashlib.sha256((HERE / 'FROZEN.json').read_bytes()).hexdigest()},
                     sort_keys=True))


if __name__ == '__main__':
    main()
