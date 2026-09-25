#!/usr/bin/env python3
"""Five-edge, six-distinct-slot full-point S6 Boolean circuit.

This is a representation only.  The charged solver campaign has a separate
release gate in PROTOCOL.md; constructing this circuit is not a PDP outcome.
"""
from __future__ import annotations

import hashlib
import json
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
sys.path.insert(0, str(NOTES / "symbolic_dag_fullpoint_20260925"))
from dag import Dag, Field, PackedModel, build_relation  # noqa: E402

O = (1, 0, 0)
ROLES = ("F0", "F1", "F2", "F3", "F4", "F5", "S2", "S3", "S4", "S5", "SUM")
EDGES = (("F0", "F1", "S2"), ("S2", "F2", "S3"),
         ("S3", "F3", "S4"), ("S4", "F4", "S5"), ("S5", "F5", "SUM"))
BETAS = (3, 338435, 303097, 464276, 42605)
ARCHIVES = {
    3: (NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz",
        "39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c",
        "raw/n19-m6/factors.json"),
    **{beta: (NOTES / "rotated_beta_sweep_20260925/evidence/raw.tar.gz",
               "fe84aef6a2cf7f6f4c950245c9c8e870354fb750666f5482b81f9a997d107140",
               f"raw/beta-{beta}/factors.json")
       for beta in BETAS[1:]},
}


def read_factors(beta: int) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    """Read exact finite full-point slots after checking the source archive."""
    archive, expected_sha, member = ARCHIVES[beta]
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected_sha:
        raise ValueError("factor archive SHA-256 mismatch")
    with tarfile.open(archive, "r:gz") as tar:
        item = tar.getmember(member)
        if not item.isfile() or item.size > 10_000:
            raise ValueError("invalid factor member")
        stream = tar.extractfile(item)
        if stream is None:
            raise ValueError("missing factor member")
        data = json.loads(stream.read())
    factors = tuple(tuple((0, *map(int, p)) for p in slot) for slot in data)
    if len(factors) != 6 or any(len(slot) != 7 for slot in factors):
        raise ValueError("expected six seven-point slots")
    return factors


class ReferenceField:
    """Bit-serial field law, independent of the Boolean DAG and exporter."""

    def __init__(self, n: int, modulus: int) -> None:
        self.n, self.modulus, self.bound = n, modulus, 1 << n

    def mul(self, a: int, b: int) -> int:
        if not 0 <= a < self.bound or not 0 <= b < self.bound:
            raise ValueError("noncanonical field operand")
        result = 0
        for _ in range(self.n):
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.bound:
                a ^= self.modulus
        return result

    def inv(self, a: int) -> int:
        if not 0 < a < self.bound:
            raise ZeroDivisionError
        # Binary polynomial Euclid; unlike the parent toy verifier this is
        # viable at n=19 and independent of a solver-provided slope.
        u, v, g, h = a, self.modulus, 1, 0
        while u != 1:
            if not u:
                raise ArithmeticError("noninvertible element")
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v, g, h = v, u, h, g
                shift = -shift
            u ^= v << shift
            g ^= h << shift
        # g may be above degree n after Euclid.
        while g.bit_length() > self.n:
            g ^= self.modulus << (g.bit_length() - self.n - 1)
        if self.mul(a, g) != 1:
            raise AssertionError("Euclid inverse failed")
        return g

    def valid(self, p: tuple[int, int, int]) -> bool:
        o, x, y = p
        if o:
            return p == O
        if not 0 <= x < self.bound or not 0 <= y < self.bound:
            return False
        return self.mul(y, y) ^ self.mul(x, y) == self.mul(self.mul(x, x), x) ^ 1

    def neg(self, p: tuple[int, int, int]) -> tuple[int, int, int]:
        return O if p == O else (0, p[1], p[1] ^ p[2])

    def scalar(self, p: tuple[int, int, int], k: int) -> tuple[int, int, int]:
        if k < 0:
            raise ValueError("negative scalar")
        acc = O
        while k:
            if k & 1:
                acc = self.add_slope(acc, p)[0]
            p = self.add_slope(p, p)[0]
            k >>= 1
        return acc

    def add_slope(self, p: tuple[int, int, int], q: tuple[int, int, int]
                  ) -> tuple[tuple[int, int, int], int, str]:
        if not self.valid(p) or not self.valid(q):
            raise ValueError("invalid source point")
        if p == O:
            return q, 0, "copy_q"
        if q == O:
            return p, 0, "copy_p"
        _, x, y = p
        _, u, v = q
        if x == u and y ^ v == x:
            return O, 0, "inverse"
        if x == u:
            if x == 0 or y != v:
                raise AssertionError("impossible same-x branch")
            lam = x ^ self.mul(y, self.inv(x))
            rx = self.mul(lam, lam) ^ lam
            ry = self.mul(x, x) ^ self.mul(lam ^ 1, rx)
            branch = "double"
        else:
            lam = self.mul(y ^ v, self.inv(x ^ u))
            rx = self.mul(lam, lam) ^ lam ^ x ^ u
            ry = self.mul(lam, x ^ rx) ^ rx ^ y
            branch = "generic"
        result = (0, rx, ry)
        if not self.valid(result):
            raise AssertionError("group addition left curve")
        return result, lam, branch


def _point_names(role: str, n: int) -> tuple[str, ...]:
    return (f"{role}_o", *(f"{role}_x_{i}" for i in range(n)),
            *(f"{role}_y_{i}" for i in range(n)))


def _point_bits(p: tuple[int, int, int], n: int) -> tuple[int, ...]:
    o, x, y = p
    if o not in (0, 1) or not 0 <= x < 1 << n or not 0 <= y < 1 << n:
        raise ValueError("noncanonical point")
    if o and p != O:
        raise ValueError("noncanonical infinity")
    return (o, *((x >> i) & 1 for i in range(n)),
            *((y >> i) & 1 for i in range(n)))


@dataclass
class Circuit:
    dag: Dag
    output: int
    names_to_nodes: dict[str, int]
    n: int
    modulus: int
    factors: tuple[tuple[tuple[int, int, int], ...], ...]
    edge_outputs: tuple[int, ...]

    def target_units(self, target: tuple[int, int, int]) -> list[int]:
        if not ReferenceField(self.n, self.modulus).valid(target):
            raise ValueError("target is not a canonical curve point")
        return [(self.names_to_nodes[name] + 1) * (1 if bit else -1)
                for name, bit in zip(_point_names("SUM", self.n),
                                     _point_bits(target, self.n), strict=True)]

    def witness(self, indices: tuple[int, ...]) -> tuple[PackedModel, tuple[int, int, int], tuple[str, ...]]:
        if len(indices) != 6 or any(not 0 <= j < 7 for j in indices):
            raise ValueError("factor-index tuple outside six slots")
        ref = ReferenceField(self.n, self.modulus)
        roles: dict[str, tuple[int, int, int]] = {
            f"F{i}": self.factors[i][j] for i, j in enumerate(indices)
        }
        slopes: dict[str, int] = {}
        branch_names = []
        for i, (left, right, out) in enumerate(EDGES):
            roles[out], slopes[f"L{i}"], branch = ref.add_slope(roles[left], roles[right])
            branch_names.append(branch)
        assignments: dict[str, int] = {}
        for role in ROLES:
            assignments.update(zip(_point_names(role, self.n),
                                   _point_bits(roles[role], self.n), strict=True))
        for label, slope in slopes.items():
            assignments.update((f"{label}_{bit}", (slope >> bit) & 1)
                               for bit in range(self.n))
        for i, chosen in enumerate(indices):
            assignments.update((f"A{i}_{j}", int(j == chosen)) for j in range(7))
        if set(assignments) != set(self.dag.names):
            raise AssertionError("witness does not cover every primary input")
        return PackedModel.from_bits([assignments[name] for name in self.dag.names]), roles["SUM"], tuple(branch_names)


def _copy_relation(master: Dag, local: Dag, output: int,
                   bindings: dict[str, int]) -> int:
    """Compose the exact parent relation with shared global role variables."""
    translated = [0, 1]
    for op, a, b in local.nodes[2:]:
        if op == "var":
            translated.append(bindings[local.names[a]])
        elif op == "xor":
            translated.append(master.xor(translated[a], translated[b]))
        elif op == "and":
            translated.append(master.and_(translated[a], translated[b]))
        else:
            raise AssertionError("unknown parent DAG node")
    return translated[output]


def build(n: int, modulus: int,
          factors: tuple[tuple[tuple[int, int, int], ...], ...]) -> Circuit:
    if len(factors) != 6 or any(len(slot) != 7 for slot in factors):
        raise ValueError("expected six seven-point factor slots")
    ref = ReferenceField(n, modulus)
    if any(len(set(slot)) != 7 or any(p == O or not ref.valid(p) for p in slot)
           for slot in factors):
        raise ValueError("factor slots must contain seven distinct finite curve points")
    master = Dag()
    Field(master, n, modulus)  # includes a Rabin irreducibility check
    point_vars: dict[str, tuple[int, ...]] = {}
    for role in ROLES:
        point_vars[role] = tuple(master.var(name) for name in _point_names(role, n))
    slope_vars = {f"L{i}": tuple(master.var(f"L{i}_{bit}") for bit in range(n))
                  for i in range(5)}
    selectors = [tuple(master.var(f"A{i}_{j}") for j in range(7)) for i in range(6)]
    names_to_nodes = {name: node_id for node_id, (op, a, _) in enumerate(master.nodes)
                      if op == "var" for name in (master.names[a],)}
    assert len(master.names) == 11 * (2 * n + 1) + 5 * n + 42

    relation = build_relation(n, modulus)
    edge_outputs = []
    for i, (left, right, out) in enumerate(EDGES):
        bindings = {name: names_to_nodes[f"{role}_{name[2:]}"]
                    for role, prefix in ((left, "p"), (right, "q"), (out, "r"))
                    for name in relation.dag.names if name.startswith(prefix + "_")}
        bindings.update({f"lambda_{bit}": slope_vars[f"L{i}"][bit] for bit in range(n)})
        if set(bindings) != set(relation.dag.names):
            raise AssertionError("local addition relation variable drift")
        edge_outputs.append(_copy_relation(master, relation.dag, relation.output, bindings))

    conditions = list(edge_outputs)
    for i, slot in enumerate(factors):
        choices = selectors[i]
        conditions.append(master.any_(list(choices)))
        for j in range(7):
            for k in range(j + 1, 7):
                conditions.append(master.not_(master.and_(choices[j], choices[k])))
            bits = _point_bits(slot[j], n)
            for node_id, bit in zip(point_vars[f"F{i}"], bits, strict=True):
                conditions.append(master.implies(choices[j], master.eq(node_id, bit)))
    output = master.all_(conditions)
    return Circuit(master, output, names_to_nodes, n, modulus, factors,
                   tuple(edge_outputs))
