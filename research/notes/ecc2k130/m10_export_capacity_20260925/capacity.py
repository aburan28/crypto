#!/usr/bin/env python3
"""Size-only complete-chain full-point DAG with implicit rotated x domains.

Every factor has a free y word. Its x word is *wired* from exactly the
slot's normal-basis selector bits, and the inherited addition relation checks
the full curve equation. No one-hot factor enumeration or solver is used.
"""
from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
sys.path.insert(0, str(NOTES / "symbolic_dag_fullpoint_20260925"))
sys.path.insert(0, str(NOTES / "rotated_subspace_support_20260925"))
from dag import Dag, Field, build_relation  # noqa: E402
import gate  # noqa: E402


class NodeCapExceeded(RuntimeError):
    def __init__(self, stage: str, counts: dict, prefix_sha256: str):
        super().__init__(f"node cap exceeded in {stage}: {counts['total_nodes']} nodes")
        self.stage, self.counts, self.prefix_sha256 = stage, counts, prefix_sha256


class CappedDag(Dag):
    """The parent DAG with a checked append cap and streaming prefix digest."""

    def __init__(self, node_cap: int):
        super().__init__()
        self.node_cap = node_cap
        self.stage = "initialization"
        self._sha = hashlib.sha256(b"zero,0,0\none,0,0\n")
        self._xor_count = 0
        self._and_count = 0

    def _appended(self, old_count: int) -> None:
        if len(self.nodes) == old_count:
            return
        op, a, b = self.nodes[-1]
        self._sha.update(f"{op},{a},{b}\n".encode("ascii"))
        self._xor_count += op == "xor"
        self._and_count += op == "and"
        if len(self.nodes) > self.node_cap:
            raise NodeCapExceeded(self.stage, self.counts(), self._sha.hexdigest())

    def var(self, name: str) -> int:
        old_count = len(self.nodes)
        result = super().var(name)
        self._appended(old_count)
        return result

    def _binary(self, op: str, a: int, b: int) -> int:
        old_count = len(self.nodes)
        result = super()._binary(op, a, b)
        self._appended(old_count)
        return result

    def counts(self) -> dict:
        return {"variables": len(self.names), "xor": self._xor_count,
                "and": self._and_count, "total_nodes": len(self.nodes),
                "model_limbs": (len(self.names) + 63) // 64}

    def prefix_sha256(self) -> str:
        return self._sha.hexdigest()


@dataclass
class Chain:
    dag: CappedDag
    output: int
    n: int
    modulus: int
    dimensions: tuple[int, ...]
    slot_bases: tuple[tuple[int, ...], ...]
    selector_names: tuple[tuple[str, ...], ...]
    roles: dict[str, tuple[int, tuple[int, ...], tuple[int, ...]]]
    edges: tuple[tuple[str, str, str], ...]
    edge_outputs: tuple[int, ...]
    local_relation_counts: dict

    def counts(self) -> dict:
        return self.dag.counts()


def slot_bases_n131(arm: str) -> tuple[tuple[int, ...], ...]:
    if arm not in ("balanced", "unequal"):
        raise ValueError("arm must be balanced or unequal")
    f = gate.Field(131, [0, 1, 2, 13])
    conjugates = gate.normal_conjugates(f, 3)
    dims = [13] * 10 if arm == "balanced" else [14] + [13] * 9
    indices = [[10 * j + i for j in range(d)] for i, d in enumerate(dims)]
    flat = [index for slot in indices for index in slot]
    assert len(flat) == len(set(flat)) == sum(dims)
    assert set(flat) == (set(range(131)) - {130} if arm == "balanced" else set(range(131)))
    bases = tuple(tuple(conjugates[index] for index in slot) for slot in indices)
    assert all(gate.rank(list(basis)) == dimension
               for basis, dimension in zip(bases, dims, strict=True))
    assert gate.rank([value for basis in bases for value in basis]) == sum(dims)
    return bases


def _point_vars(d: CappedDag, label: str, n: int):
    return (d.var(f"{label}_o"),
            tuple(d.var(f"{label}_x_{bit}") for bit in range(n)),
            tuple(d.var(f"{label}_y_{bit}") for bit in range(n)))


def _factor_vars(d: CappedDag, label: str, n: int, basis: tuple[int, ...]):
    selectors = tuple(d.var(f"{label}_a_{j}") for j in range(len(basis)))
    y = tuple(d.var(f"{label}_y_{bit}") for bit in range(n))
    x = []
    for bit in range(n):
        wire = 0
        for j, value in enumerate(basis):
            if (value >> bit) & 1:
                wire = d.xor(wire, selectors[j])
        x.append(wire)
    return (0, tuple(x), y), tuple(f"{label}_a_{j}" for j in range(len(basis)))


def _copy_relation(master: CappedDag, local: Dag, local_output: int,
                   bindings: dict[str, int]) -> int:
    translated = [0, 1]
    for op, a, b in local.nodes[2:]:
        if op == "var":
            translated.append(bindings[local.names[a]])
        elif op == "xor":
            translated.append(master.xor(translated[a], translated[b]))
        elif op == "and":
            translated.append(master.and_(translated[a], translated[b]))
        else:
            raise AssertionError("unexpected parent DAG node")
    return translated[local_output]


def build_chain(n: int, modulus: int, bases: tuple[tuple[int, ...], ...],
                node_cap: int, checkpoint: Callable[[dict], None] | None = None) -> Chain:
    m = len(bases)
    if m < 2 or any(not basis for basis in bases):
        raise ValueError("chain needs at least two nonempty factor slots")
    if any(gate.rank(list(basis)) != len(basis) for basis in bases):
        raise ValueError("dependent slot basis")
    if any(value < 0 or value >= (1 << n) for basis in bases for value in basis):
        raise ValueError("noncanonical slot basis")
    local = build_relation(n, modulus)
    local_counts = local.dag.counts()
    d = CappedDag(node_cap)
    Field(d, n, modulus)
    roles = {}
    selectors = []
    for i, basis in enumerate(bases):
        d.stage = f"factor_{i}_domain"
        roles[f"F{i}"], names = _factor_vars(d, f"F{i}", n, basis)
        selectors.append(names)
    for i in range(2, m):
        d.stage = f"prefix_{i}_inputs"
        roles[f"S{i}"] = _point_vars(d, f"S{i}", n)
    d.stage = "target_inputs"
    roles["SUM"] = _point_vars(d, "SUM", n)
    slopes = []
    for i in range(m - 1):
        d.stage = f"slope_{i}_inputs"
        slopes.append(tuple(d.var(f"L{i}_{bit}") for bit in range(n)))
    if checkpoint:
        checkpoint({"stage": "primary_wires", "counts": d.counts(),
                    "prefix_sha256": d.prefix_sha256(),
                    "local_relation_counts": local_counts})
    edges = []
    edge_outputs = []
    acc = "F0"
    for i in range(1, m):
        out = "SUM" if i == m - 1 else f"S{i+1}"
        right = f"F{i}"
        edges.append((acc, right, out))
        d.stage = f"edge_{i-1}"
        bindings = {}
        for label, prefix in ((acc, "p"), (right, "q"), (out, "r")):
            o, x, y = roles[label]
            bindings[f"{prefix}_o"] = o
            bindings.update({f"{prefix}_x_{bit}": node for bit, node in enumerate(x)})
            bindings.update({f"{prefix}_y_{bit}": node for bit, node in enumerate(y)})
        bindings.update({f"lambda_{bit}": node
                         for bit, node in enumerate(slopes[i-1])})
        if set(bindings) != set(local.dag.names):
            raise AssertionError("parent addition interface drift")
        edge_outputs.append(_copy_relation(d, local.dag, local.output, bindings))
        acc = out
        if checkpoint:
            checkpoint({"stage": d.stage, "counts": d.counts(),
                        "prefix_sha256": d.prefix_sha256()})
    d.stage = "chain_output"
    output = d.all_(edge_outputs)
    if checkpoint:
        checkpoint({"stage": d.stage, "counts": d.counts(),
                    "prefix_sha256": d.prefix_sha256(),
                    "output_node": output})
    return Chain(d, output, n, modulus, tuple(map(len, bases)), bases,
                 tuple(selectors), roles, tuple(edges), tuple(edge_outputs), local_counts)


def _clauses(op: str, a: int, b: int, z: int):
    if op == "xor":
        return ((-a, -b, -z), (a, b, -z), (a, -b, z), (-a, b, z))
    if op == "and":
        return ((-a, -b, z), (a, -z), (b, -z))
    raise ValueError("unexpected gate")


def exact_dimacs_size(chain: Chain, fixed_target_O: bool = True) -> dict:
    """Count exact #804 ASCII DIMACS bytes without writing a large CNF file."""
    d = chain.dag
    counts = d.counts()
    variables = counts["total_nodes"]
    if variables > (1 << 31) - 1:
        raise ValueError("DIMACS signed variable ID overflow")
    units = []
    if fixed_target_O:
        o, x, y = chain.roles["SUM"]
        units = [o + 1] + [-(node + 1) for node in (*x, *y)]
    clauses = 2 + 4 * counts["xor"] + 3 * counts["and"] + 1 + len(units)
    total = len(f"p cnf {variables} {clauses}\n".encode("ascii"))
    total += len("-1 0\n") + len("2 0\n")
    for node_id, (op, a, b) in enumerate(d.nodes[2:], 2):
        if op == "var":
            continue
        for clause in _clauses(op, a + 1, b + 1, node_id + 1):
            total += len((" ".join(map(str, clause)) + " 0\n").encode("ascii"))
    total += len(f"{chain.output + 1} 0\n")
    for unit in units:
        total += len(f"{unit} 0\n")
    return {"variables": variables, "clauses": clauses, "bytes": total,
            "unit_count": len(units), "target": "O" if fixed_target_O else "generic",
            "dag": counts, "dag_prefix_sha256": d.prefix_sha256()}
