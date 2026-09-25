#!/usr/bin/env python3
"""Five-distinct-slot n13 full-point chain, stacked on the S6 circuit gate."""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'n19_s6_fullpoint_exporter_20260925'))
from s6 import (Dag, Field, O, PackedModel, ReferenceField, _copy_relation,
                _point_bits, _point_names, build_relation)  # noqa: E402

ROLES = ('F0', 'F1', 'F2', 'F3', 'F4', 'S2', 'S3', 'S4', 'SUM')
EDGES = (('F0', 'F1', 'S2'), ('S2', 'F2', 'S3'),
         ('S3', 'F3', 'S4'), ('S4', 'F4', 'SUM'))
N, MODULUS = 13, 0x201B


def input_data() -> dict:
    value = json.loads((HERE / 'INPUT.json').read_text())
    if value['schema'] != 'ecc2k130_n13_s5_control_input_v1':
        raise ValueError('n13 S5 input schema drift')
    return value


def read_factors() -> tuple[tuple[tuple[int, int, int], ...], ...]:
    data = input_data()
    return tuple(tuple(tuple(p) for p in slot) for slot in data['factors'])


def read_targets() -> tuple[dict, ...]:
    return tuple(input_data()['targets'])


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
        """Pin the exact full target Q+T_torsion; SUM is not T_torsion."""
        if not ReferenceField(self.n, self.modulus).valid(target):
            raise ValueError('target is not a canonical curve point')
        return [(self.names_to_nodes[name] + 1) * (1 if bit else -1)
                for name, bit in zip(_point_names('SUM', self.n),
                                     _point_bits(target, self.n), strict=True)]

    def witness(self, indices: tuple[int, ...]) -> tuple[PackedModel, tuple[int, int, int], tuple[str, ...]]:
        if len(indices) != 5 or any(not 0 <= j < len(self.factors[i])
                                    for i, j in enumerate(indices)):
            raise ValueError('factor-index tuple outside five distinct slots')
        ref = ReferenceField(self.n, self.modulus)
        roles = {f'F{i}': self.factors[i][j] for i, j in enumerate(indices)}
        slopes = {}
        branches = []
        for i, (left, right, out) in enumerate(EDGES):
            roles[out], slopes[f'L{i}'], branch = ref.add_slope(roles[left], roles[right])
            branches.append(branch)
        assigned = {}
        for role in ROLES:
            assigned.update(zip(_point_names(role, self.n),
                                _point_bits(roles[role], self.n), strict=True))
        for role, slope in slopes.items():
            assigned.update((f'{role}_{bit}', (slope >> bit) & 1)
                            for bit in range(self.n))
        for i, chosen in enumerate(indices):
            assigned.update((f'A{i}_{j}', int(j == chosen))
                            for j in range(len(self.factors[i])))
        if set(assigned) != set(self.dag.names):
            raise AssertionError('complete primary model coverage failed')
        return (PackedModel.from_bits([assigned[name] for name in self.dag.names]),
                roles['SUM'], tuple(branches))


def build(n: int, modulus: int,
          factors: tuple[tuple[tuple[int, int, int], ...], ...]) -> Circuit:
    if len(factors) != 5 or any(len(slot) != 5 for slot in factors):
        raise ValueError('expected five five-point factor slots')
    ref = ReferenceField(n, modulus)
    if any(len(set(slot)) != 5 or any(p == O or not ref.valid(p) for p in slot)
           for slot in factors):
        raise ValueError('every factor must be a distinct finite full point')
    d = Dag()
    Field(d, n, modulus)  # Rabin irreducibility check
    point_nodes = {role: tuple(d.var(name) for name in _point_names(role, n))
                   for role in ROLES}
    slopes = {f'L{i}': tuple(d.var(f'L{i}_{bit}') for bit in range(n))
              for i in range(4)}
    selectors = [tuple(d.var(f'A{i}_{j}') for j in range(5)) for i in range(5)]
    names_to_nodes = {d.names[a]: node_id for node_id, (op, a, _) in enumerate(d.nodes)
                      if op == 'var'}
    assert len(d.names) == 9 * (2 * n + 1) + 4 * n + 25
    local = build_relation(n, modulus)
    edge_outputs = []
    for i, (left, right, out) in enumerate(EDGES):
        bindings = {name: names_to_nodes[f'{role}_{name[2:]}']
                    for role, prefix in ((left, 'p'), (right, 'q'), (out, 'r'))
                    for name in local.dag.names if name.startswith(prefix + '_')}
        bindings.update({f'lambda_{bit}': slopes[f'L{i}'][bit]
                         for bit in range(n)})
        if set(bindings) != set(local.dag.names):
            raise AssertionError('local relation variable drift')
        edge_outputs.append(_copy_relation(d, local.dag, local.output, bindings))
    constraints = list(edge_outputs)
    for i, slot in enumerate(factors):
        choices = selectors[i]
        constraints.append(d.any_(list(choices)))
        for j in range(5):
            for k in range(j + 1, 5):
                constraints.append(d.not_(d.and_(choices[j], choices[k])))
            for node_id, bit in zip(point_nodes[f'F{i}'], _point_bits(slot[j], n), strict=True):
                constraints.append(d.implies(choices[j], d.eq(node_id, bit)))
    output = d.all_(constraints)
    return Circuit(d, output, names_to_nodes, n, modulus, factors,
                   tuple(edge_outputs))
