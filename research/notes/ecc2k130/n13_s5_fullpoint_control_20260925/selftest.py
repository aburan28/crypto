#!/usr/bin/env python3
"""Small-field S5 Boolean/point model regression; no n13 solver outcome."""
from __future__ import annotations

import itertools
import json
import sys
import tempfile
from pathlib import Path

from admission5 import decode
from s5 import ReferenceField, build, input_data

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from export import write_cnf  # noqa: E402
from verify import evaluate_cnf, extend_dag, parse_cnf  # noqa: E402
sys.path.insert(0, str(HERE.parent / 'n19_s6_fullpoint_exporter_20260925'))
from admission import check_cnf  # noqa: E402


def main() -> None:
    ref = ReferenceField(2, 0x7)
    finite = tuple((0, x, y) for x in range(4) for y in range(4)
                   if ref.valid((0, x, y)))[:5]
    assert len(finite) == 5 and (0, 0, 1) in finite
    c = build(2, 0x7, (finite,) * 5)
    assert len(c.dag.names) == 9 * 5 + 4 * 2 + 25 == 78
    representatives = {}
    for indices in itertools.product(range(5), repeat=5):
        _, _, branches = c.witness(indices)
        for branch in branches:
            representatives.setdefault(branch, indices)
        if set(representatives) == {'inverse', 'double', 'generic', 'copy_q'}:
            break
    assert set(representatives) == {'inverse', 'double', 'generic', 'copy_q'}
    checked = set(representatives.values())
    checked.add((0, 0, 0, 0, 0))
    with tempfile.TemporaryDirectory() as temp:
        cnf = Path(temp) / 'toy.cnf'
        meta = write_cnf(c, cnf, byte_cap=16 << 20)
        assert check_cnf(c, cnf, [])['sha256'] == meta['sha256']
        _, _, clauses = parse_cnf(cnf)
        for indices in sorted(checked):
            packed, target, _ = c.witness(indices)
            assert c.dag.evaluate(packed, c.output)
            full = extend_dag(c.dag, packed)
            assert evaluate_cnf(clauses, full)
            lifted = decode(c, full, target)
            assert lifted['indices'] == list(indices)
            assert all((full[abs(lit) - 1] == 1) == (lit > 0)
                       for lit in c.target_units(target))
        wrong = full.copy()
        wrong[c.names_to_nodes['A0_0']] ^= 1
        assert not evaluate_cnf(clauses, wrong)
    frozen = input_data()
    assert len(frozen['targets']) == 32 and len(frozen['factors']) == 5
    assert frozen['exceptional']['factor_indices'] == [0, 0, 0, 4, 1]
    print(json.dumps({'decision': 'TOY_SELFTEST_PASS', 'n': 2,
                      'primary_bits': len(c.dag.names), 'cnf': meta,
                      'lifted_representatives': len(checked),
                      'branches': sorted(representatives),
                      'n13_target_labels_pinned': len(frozen['targets'])},
                     sort_keys=True))


if __name__ == '__main__':
    main()
