#!/usr/bin/env python3
"""Complete three-factor native-leaf chain with implicit finite x domains."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
sys.path.insert(0, str(NOTES / "m10_export_capacity_20260925"))
sys.path.insert(0, str(NOTES / "native_fullpoint_edge_20260929"))
from capacity import (CappedDag, Chain, NodeCapExceeded, _copy_relation,
                      _factor_vars, _point_vars, exact_dimacs_size)  # noqa: E402
from native_relation import build_relation  # noqa: E402
from dag import Field  # noqa: E402
import gate  # noqa: E402 - imported by the frozen source-chain capacity builder


EXACT_MODULUS = 0x800000000000000000000000000002007
LEAF_DIMS = (44, 44, 43)


def leaf_slot_bases() -> tuple[tuple[int, ...], ...]:
    """Three disjoint normal-basis slices spanning all 131 coordinates."""
    field = gate.Field(131, [0, 1, 2, 13])
    conjugates = gate.normal_conjugates(field, 3)
    assert len(conjugates) == 131 and gate.rank(conjugates) == 131
    indices = tuple(tuple(3 * j + i for j in range(d))
                    for i, d in enumerate(LEAF_DIMS))
    assert set(index for slot in indices for index in slot) == set(range(131))
    bases = tuple(tuple(conjugates[index] for index in slot) for slot in indices)
    assert all(gate.rank(list(basis)) == len(basis) for basis in bases)
    assert gate.rank([value for basis in bases for value in basis]) == 131
    return bases


def build_native_chain(n: int, modulus: int, a: int, b: int,
                       bases: tuple[tuple[int, ...], ...], node_cap: int,
                       checkpoint=None) -> Chain:
    if len(bases) != 3 or any(not basis for basis in bases):
        raise ValueError("exactly three nonempty factor domains required")
    if any(value < 0 or value >= 1 << n for basis in bases for value in basis):
        raise ValueError("noncanonical factor-domain vector")
    if any(gate.rank(list(basis)) != len(basis) for basis in bases):
        raise ValueError("dependent factor-domain basis")
    local = build_relation(n, modulus, a, b)
    d = CappedDag(node_cap)
    Field(d, n, modulus)
    roles = {}
    selector_names = []
    for i, basis in enumerate(bases):
        d.stage = f"factor_{i}_domain"
        roles[f"F{i}"], names = _factor_vars(d, f"F{i}", n, basis)
        selector_names.append(names)
    d.stage = "prefix_inputs"
    roles["S2"] = _point_vars(d, "S2", n)
    d.stage = "target_inputs"
    roles["SUM"] = _point_vars(d, "SUM", n)
    slopes = []
    for i in range(2):
        d.stage = f"slope_{i}_inputs"
        slopes.append(tuple(d.var(f"L{i}_{bit}") for bit in range(n)))
    if checkpoint:
        checkpoint({"stage": "primary_wires", "counts": d.counts(),
                    "prefix_sha256": d.prefix_sha256(),
                    "local_relation_counts": local.dag.counts()})
    edges = (("F0", "F1", "S2"), ("S2", "F2", "SUM"))
    edge_outputs = []
    for edge_index, (left, right, out) in enumerate(edges):
        d.stage = f"edge_{edge_index}"
        bindings = {}
        for role, prefix in ((left, "p"), (right, "q"), (out, "r")):
            o, x, y = roles[role]
            bindings[f"{prefix}_o"] = o
            bindings.update({f"{prefix}_x_{bit}": node for bit, node in enumerate(x)})
            bindings.update({f"{prefix}_y_{bit}": node for bit, node in enumerate(y)})
        bindings.update({f"lambda_{bit}": node
                         for bit, node in enumerate(slopes[edge_index])})
        if set(bindings) != set(local.dag.names):
            raise AssertionError("full-point edge input interface drift")
        edge_outputs.append(_copy_relation(d, local.dag, local.output, bindings))
        if checkpoint:
            checkpoint({"stage": d.stage, "counts": d.counts(),
                        "prefix_sha256": d.prefix_sha256()})
    d.stage = "chain_output"
    output = d.all_(edge_outputs)
    if checkpoint:
        checkpoint({"stage": d.stage, "counts": d.counts(),
                    "prefix_sha256": d.prefix_sha256(), "output_node": output})
    return Chain(d, output, n, modulus, tuple(map(len, bases)), bases,
                 tuple(selector_names), roles, edges, tuple(edge_outputs),
                 local.dag.counts())
