#!/usr/bin/env python3
"""Small non-n131 regression for implicit-domain wiring and byte accounting."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "symbolic_dag_fullpoint_20260925"))
sys.path.insert(0, str(HERE.parent / "symbolic_dag_dimacs_gate_20260925"))
sys.path.insert(0, str(HERE))
from dag import PackedModel  # noqa: E402
from export import write_cnf  # noqa: E402
from verify import independent_dimacs_size  # noqa: E402
from capacity import build_chain, exact_dimacs_size  # noqa: E402


def point_bits(label: str, point, n: int):
    o, x, y = point
    return {f"{label}_o": o,
            **{f"{label}_x_{i}": (x >> i) & 1 for i in range(n)},
            **{f"{label}_y_{i}": (y >> i) & 1 for i in range(n)}}


def witness(chain, factors, intermediate, total, slopes):
    values = {}
    for i, (_, x, y) in enumerate(factors):
        assert x in (0, 1)
        values[f"F{i}_a_0"] = x
        values.update({f"F{i}_y_{bit}": (y >> bit) & 1
                       for bit in range(chain.n)})
    values.update(point_bits("S2", intermediate, chain.n))
    values.update(point_bits("SUM", total, chain.n))
    for i, slope in enumerate(slopes):
        values.update({f"L{i}_{bit}": (slope >> bit) & 1
                       for bit in range(chain.n)})
    assert set(values) == set(chain.dag.names)
    return PackedModel.from_bits([values[name] for name in chain.dag.names])


def main() -> dict:
    # These four points on E/GF(2^3) exercise x=0 and x=1 toy slots.
    chain = build_chain(3, 0xb, ((1,), (1,), (1,)), 20000)
    assert chain.counts()["variables"] == 3 + 3 * 3 + 2 * (3 * 3 + 1)
    doubled = witness(chain, ((0, 1, 0), (0, 1, 0), (0, 1, 1)),
                      (0, 0, 1), (0, 1, 0), (1, 0))
    assert chain.dag.evaluate(doubled, chain.output)
    inverse_then_copy = witness(chain, ((0, 0, 1), (0, 0, 1), (0, 0, 1)),
                                (1, 0, 0), (0, 0, 1), (0, 0))
    assert chain.dag.evaluate(inverse_then_copy, chain.output)
    wrong = list(doubled.bit(i) for i in range(doubled.bit_count))
    wrong[chain.dag.names.index("SUM_y_0")] ^= 1
    assert not chain.dag.evaluate(PackedModel.from_bits(wrong), chain.output)
    generic = exact_dimacs_size(chain, fixed_target_O=False)
    target_o = exact_dimacs_size(chain, fixed_target_O=True)
    assert generic == independent_dimacs_size(chain, False)
    assert target_o == independent_dimacs_size(chain, True)
    with tempfile.TemporaryDirectory(prefix="m10_capacity_toy_") as directory:
        p = Path(directory) / "generic.cnf"
        written = write_cnf(chain, p, byte_cap=1 << 20)
        assert (written["variables"], written["clauses"], written["bytes"]) == (
            generic["variables"], generic["clauses"], generic["bytes"])
        o, x, y = chain.roles["SUM"]
        units = [o + 1] + [-(node + 1) for node in (*x, *y)]
        written_o = write_cnf(chain, Path(directory) / "target_O.cnf",
                              units=units, byte_cap=1 << 20)
        assert (written_o["variables"], written_o["clauses"], written_o["bytes"]) == (
            target_o["variables"], target_o["clauses"], target_o["bytes"])
    return {"schema": "ecc2k130-m10-capacity-toy-selftest-v1", "status": "PASS",
            "toy_field_degree": 3, "toy_factor_slots": 3,
            "primary_bits": chain.counts()["variables"],
            "dag_nodes": chain.counts()["total_nodes"],
            "generic_cnf_bytes": generic["bytes"],
            "target_O_cnf_bytes": target_o["bytes"],
            "checked_branches": ["double", "generic", "inverse", "copy"]}


if __name__ == "__main__":
    print(json.dumps(main(), sort_keys=True, separators=(",", ":")))
