#!/usr/bin/env python3
"""Algebra/Boolean toy control only; never launches a PDP solver."""
from __future__ import annotations

import itertools
import json
import tempfile
from pathlib import Path

from admission import check_cnf, check_full_model
from s6 import ReferenceField, build

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE.parent / "symbolic_dag_dimacs_gate_20260925"))
from export import write_cnf  # noqa: E402
from verify import evaluate_cnf, extend_dag, parse_cnf  # noqa: E402


def main() -> None:
    n, modulus = 2, 0x7
    field = ReferenceField(n, modulus)
    finite = tuple((0, x, y) for x in range(1 << n) for y in range(1 << n)
                   if field.valid((0, x, y)))
    assert len(finite) == 7 and (0, 0, 1) in finite
    circuit = build(n, modulus, (finite,) * 6)
    assert len(circuit.dag.names) == 11 * (2 * n + 1) + 5 * n + 42
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "toy.cnf"
        meta = write_cnf(circuit, path, byte_cap=16 << 20)
        assert check_cnf(circuit, path, [])["sha256"] == meta["sha256"]
        _, _, clauses = parse_cnf(path)
        # Exact distinct-slot labelled tuples are census-inspected.  One
        # witness for each branch type and O-prefix is lifted through every
        # Tseitin gate and independently parsed CNF below.
        representatives = {}
        targets = set()
        for indices in itertools.product(range(7), repeat=6):
            _, target, branches = circuit.witness(indices)
            targets.add(target)
            for branch in branches:
                representatives.setdefault(branch, indices)
            # The second summand is always a finite factor, so copy_p is
            # absent from an S6 factor chain. The local gate tests it separately.
            if all(label in representatives for label in
                   ("inverse", "double", "generic", "copy_q")) and len(targets) == 8:
                break
        assert len(targets) == 8 and len(representatives) == 4
        checked = set(representatives.values())
        checked.update(((0, 0, 0, 0, 0, 0), (6, 6, 6, 6, 6, 6),
                        (0, 1, 2, 3, 4, 5)))
        for indices in sorted(checked):
            packed, target, branches = circuit.witness(indices)
            assert circuit.dag.evaluate(packed, circuit.output)
            full = extend_dag(circuit.dag, packed)
            assert evaluate_cnf(clauses, full)
            assert check_full_model(circuit, clauses, full, target)["indices"] == list(indices)
            fixed = circuit.target_units(target)
            assert len(fixed) == 2 * n + 1
            assert all((full[abs(lit) - 1] == 1) == (lit > 0) for lit in fixed)
        wrong = full.copy()
        wrong[circuit.names_to_nodes["A0_0"]] ^= 1
        assert not evaluate_cnf(clauses, wrong)
        print(json.dumps({"decision": "TOY_SELFTEST_PASS", "n": n,
                          "distinct_finite_points": len(finite),
                          "lifted_representatives": len(checked),
                          "branch_types": sorted(representatives),
                          "cnf": meta}, sort_keys=True))


if __name__ == "__main__":
    main()
