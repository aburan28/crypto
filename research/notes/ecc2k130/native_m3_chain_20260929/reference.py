#!/usr/bin/env python3
"""Independent toy group oracle, raw DIMACS checker and exact model lift."""
from __future__ import annotations

import hashlib
import itertools
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "native_fullpoint_edge_20260929"))
from verify import PolynomialCurve  # noqa: E402 - independent long-division point law

O = (1, 0, 0)
TOY_N = 5
TOY_MODULUS = 0x25
TOY_A = 0
TOY_B = 1
TOY_BASES = ((1,), (2,), (4,))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank(vectors: list[int]) -> int:
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def square_mod(value: int, modulus: int, n: int) -> int:
    product = 0
    while value:
        low = value & -value
        product ^= 1 << (2 * (low.bit_length() - 1))
        value ^= low
    while product.bit_length() > n:
        product ^= modulus << (product.bit_length() - n - 1)
    return product


def independent_leaf_bases() -> tuple[tuple[int, ...], ...]:
    n, modulus = 131, 0x800000000000000000000000000002007
    conjugates = [3]
    for _ in range(n - 1):
        conjugates.append(square_mod(conjugates[-1], modulus, n))
    assert square_mod(conjugates[-1], modulus, n) == conjugates[0]
    assert rank(conjugates) == n
    dims = (44, 44, 43)
    bases = tuple(tuple(conjugates[3 * j + i] for j in range(dim))
                  for i, dim in enumerate(dims))
    assert all(rank(list(slot)) == len(slot) for slot in bases)
    assert rank([v for slot in bases for v in slot]) == n
    return bases


def toy_oracle() -> dict:
    curve = PolynomialCurve(TOY_N, TOY_MODULUS, TOY_A, TOY_B)
    points = curve.points()
    factors = []
    for basis in TOY_BASES:
        xs = {0, basis[0]}
        factors.append(tuple(point for point in points
                             if point != O and point[1] in xs))
    assert all(factors)
    sums = {}
    witnesses = {}
    tuples = 0
    for triple in itertools.product(*factors):
        s2 = curve.add(triple[0], triple[1])
        total = curve.add(s2, triple[2])
        sums.setdefault(total, (triple, s2))
        witnesses.setdefault(total, []).append((triple, s2))
        tuples += 1
    positives = sorted(point for point in sums if point != O)
    negatives = sorted(point for point in points if point not in sums)
    assert positives and negatives
    positive, negative = positives[0], negatives[0]
    generic = sorted(point for point, paths in witnesses.items()
                     if point != O and all(branch(path[0][0], path[0][1]) == "generic"
                                           and branch(path[1], path[0][2]) == "generic"
                                           for path in paths))
    assert generic
    return {"curve": curve, "points": points, "factors": factors,
            "support": sums, "witnesses": witnesses, "triples": tuples,
            "positive": positive, "generic_positive": generic[0],
            "negative": negative}


def branch(p, q):
    if p == O:
        return "copy_q"
    if q == O:
        return "copy_p"
    if p[1] == q[1] and p[2] ^ q[2] == p[1]:
        return "inverse"
    if p[1] == q[1]:
        return "double"
    return "generic"


def primary_model(chain, triple, s2, target, slopes, masks):
    from dag import PackedModel  # imported from the frozen edge's DAG
    n = chain.n
    values = {}
    for i, (point, mask) in enumerate(zip(triple, masks, strict=True)):
        _, x, y = point
        calculated = 0
        for j, vector in enumerate(chain.slot_bases[i]):
            values[f"F{i}_a_{j}"] = (mask >> j) & 1
            if mask >> j & 1:
                calculated ^= vector
        assert x == calculated
        values.update({f"F{i}_y_{bit}": (y >> bit) & 1 for bit in range(n)})
    for label, point in (("S2", s2), ("SUM", target)):
        o, x, y = point
        values[f"{label}_o"] = o
        values.update({f"{label}_x_{bit}": (x >> bit) & 1 for bit in range(n)})
        values.update({f"{label}_y_{bit}": (y >> bit) & 1 for bit in range(n)})
    for i, slope in enumerate(slopes):
        values.update({f"L{i}_{bit}": (slope >> bit) & 1 for bit in range(n)})
    assert set(values) == set(chain.dag.names)
    return PackedModel.from_bits([values[name] for name in chain.dag.names])


def decode_model(chain, bits: list[int], curve: PolynomialCurve, target):
    """Lift all primary inputs and independently check every group operation."""
    assert len(bits) == len(chain.dag.nodes) and all(bit in (0, 1) for bit in bits)
    inputs = {chain.dag.names[index]: bits[node_id]
              for node_id, (op, index, _) in enumerate(chain.dag.nodes)
              if op == "var"}
    factors = []
    masks = []
    for i, basis in enumerate(chain.slot_bases):
        mask = sum(inputs[f"F{i}_a_{j}"] << j for j in range(len(basis)))
        x = 0
        for j, vector in enumerate(basis):
            if mask >> j & 1:
                x ^= vector
        y = sum(inputs[f"F{i}_y_{bit}"] << bit for bit in range(chain.n))
        point = (0, x, y)
        assert curve.affine(x, y)
        factors.append(point)
        masks.append(mask)
    def point_role(label):
        return (inputs[f"{label}_o"],
                sum(inputs[f"{label}_x_{bit}"] << bit for bit in range(chain.n)),
                sum(inputs[f"{label}_y_{bit}"] << bit for bit in range(chain.n)))
    s2, total = point_role("S2"), point_role("SUM")
    assert s2 == curve.add(factors[0], factors[1])
    assert total == curve.add(s2, factors[2]) == target
    for i, (left, right) in enumerate(((factors[0], factors[1]), (s2, factors[2]))):
        slope = sum(inputs[f"L{i}_{bit}"] << bit for bit in range(chain.n))
        kind = ("copy" if left == O or right == O else
                "inverse" if left[1] == right[1] and left[2] ^ right[2] == left[1]
                else "double" if left[1] == right[1] else "generic")
        if kind in ("double", "generic"):
            assert slope == curve.slope(left, right)
    return {"factors": [list(point) for point in factors], "masks": masks,
            "s2": list(s2), "target": list(total)}


def parse_cnf_and_model(path: Path, solver_output: str, exit_code: int):
    """Read every raw clause and require a complete satisfying SAT assignment."""
    with path.open("rt", encoding="ascii") as stream:
        header = stream.readline().split()
        assert header[:2] == ["p", "cnf"] and len(header) == 4
        variables, declared = map(int, header[2:])
        clauses = []
        for line in stream:
            literals = [int(token) for token in line.split()]
            assert len(literals) >= 2 and literals[-1] == 0
            assert 0 not in literals[:-1]
            assert all(1 <= abs(lit) <= variables for lit in literals[:-1])
            clauses.append(literals[:-1])
    assert len(clauses) == declared
    status = [line[2:].strip() for line in solver_output.splitlines()
              if line.startswith("s ")]
    assert len(status) == 1
    if status == ["UNSATISFIABLE"]:
        assert exit_code == 20
        return "UNSAT", None, variables, declared
    assert status == ["SATISFIABLE"] and exit_code == 10
    literals = []
    terminators = 0
    for line in solver_output.splitlines():
        if line.startswith("v "):
            for value in map(int, line[2:].split()):
                if value == 0:
                    terminators += 1
                else:
                    assert terminators == 0
                    literals.append(value)
        else:
            assert line.startswith("c ") or line.startswith("s ") or not line.strip()
    assert terminators == 1 and len(literals) == variables
    bits = [None] * variables
    for literal in literals:
        index = abs(literal) - 1
        assert index < variables and bits[index] is None
        bits[index] = int(literal > 0)
    assert all(bit in (0, 1) for bit in bits)
    assert all(any((bits[abs(lit) - 1] == 1) == (lit > 0) for lit in clause)
               for clause in clauses)
    return "SAT", bits, variables, declared
