#!/usr/bin/env python3
"""Exact-query n13 S5 SAT and DRAT admission, with independent Fermat replay."""
from __future__ import annotations

import importlib.util
from functools import lru_cache
import sys
from pathlib import Path

from s5 import Circuit, EDGES, ROLES, read_targets

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'n19_s6_fullpoint_exporter_20260925'))
from admission import admit_unsat as _admit_unsat, check_cnf  # noqa: E402
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from verify import evaluate_cnf, extend_dag, parse_cnf, parse_solver_output  # noqa: E402
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_fullpoint_20260925'))
from dag import PackedModel  # noqa: E402


@lru_cache(maxsize=4)
def fermat_reference(n: int, modulus: int):
    path = HERE.parent / 'rotated_subspace_support_20260925/verify.py'
    spec = importlib.util.spec_from_file_location('s5_fermat_reference', path)
    if spec is None or spec.loader is None:
        raise AssertionError('missing independent field verifier')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.E(mod.GF(n, modulus))


def _int(bits: list[int]) -> int:
    return sum(bit << i for i, bit in enumerate(bits))


def _bare(p: tuple[int, int, int]):
    return None if p[0] else (p[1], p[2])


def _tagged(p):
    return (1, 0, 0) if p is None else (0, *p)


def _role(inputs: dict[str, int], name: str, n: int) -> tuple[int, int, int]:
    return (inputs[f'{name}_o'],
            _int([inputs[f'{name}_x_{i}'] for i in range(n)]),
            _int([inputs[f'{name}_y_{i}'] for i in range(n)]))


def decode(circuit: Circuit, values: list[int], target: tuple[int, int, int]) -> dict:
    if len(values) != len(circuit.dag.nodes) or any(v not in (0, 1) for v in values):
        raise ValueError('incomplete/non-Boolean exact model')
    inputs = {circuit.dag.names[a]: values[node_id]
              for node_id, (op, a, _) in enumerate(circuit.dag.nodes) if op == 'var'}
    if len(inputs) != len(circuit.dag.names):
        raise AssertionError('primary input map incomplete')
    roles = {name: _role(inputs, name, circuit.n) for name in ROLES}
    ref = fermat_reference(circuit.n, circuit.modulus)
    if any((p[0] and p != (1, 0, 0)) or p[0] not in (0, 1) or
           (not p[0] and not ref.on(_bare(p))) for p in roles.values()):
        raise AssertionError('noncanonical or off-curve point')
    if roles['SUM'] != target:
        raise AssertionError('full-target pin does not match SUM')
    indices = []
    for i in range(5):
        chosen = [j for j in range(5) if inputs[f'A{i}_{j}']]
        if len(chosen) != 1 or roles[f'F{i}'] != circuit.factors[i][chosen[0]]:
            raise AssertionError('one-hot selector/full-point link failed')
        indices.append(chosen[0])
    branches = []
    for i, (left, right, out) in enumerate(EDGES):
        p, q = _bare(roles[left]), _bare(roles[right])
        expected = ref.add(p, q)
        if roles[out] != _tagged(expected):
            raise AssertionError('independent full-point chain failed')
        if p is None:
            branch, slope = 'copy_q', None
        elif q is None:
            branch, slope = 'copy_p', None
        elif p[0] == q[0] and p[1] ^ q[1] == p[0]:
            branch, slope = 'inverse', None
        elif p[0] == q[0]:
            branch = 'double'
            slope = p[0] ^ ref.f.mul(p[1], ref.f.inv(p[0]))
        else:
            branch = 'generic'
            slope = ref.f.mul(p[1] ^ q[1], ref.f.inv(p[0] ^ q[0]))
        model_slope = _int([inputs[f'L{i}_{bit}'] for bit in range(circuit.n)])
        if slope is not None and model_slope != slope:
            raise AssertionError('model slope differs from independent Fermat law')
        branches.append(branch)
    return {'indices': indices, 'full_target': list(target), 'branches': branches,
            'factors': [list(roles[f'F{i}']) for i in range(5)],
            'prefixes': [list(roles[f'S{i}']) for i in range(2, 5)]}


def _target(label: str) -> tuple[int, int, int]:
    matches = [row for row in read_targets() if row['id'] == label]
    if len(matches) != 1:
        raise ValueError('unknown or duplicate n13 branch label')
    return tuple(matches[0]['full_target'])


def admit_sat(circuit: Circuit, label: str, query: Path,
              solver_stdout: Path, exit_code: int) -> dict:
    target = _target(label)
    check_cnf(circuit, query, circuit.target_units(target))
    variables, _, clauses = parse_cnf(query)
    status, values = parse_solver_output(solver_stdout.read_text(), variables, exit_code)
    if status != 'SAT' or values is None or not evaluate_cnf(clauses, values):
        raise AssertionError('no complete satisfying SAT model')
    primary = [None] * len(circuit.dag.names)
    for node_id, (op, a, _) in enumerate(circuit.dag.nodes):
        if op == 'var':
            primary[a] = values[node_id]
    packed = PackedModel.from_bits(primary)
    if extend_dag(circuit.dag, packed) != values or not circuit.dag.evaluate(packed, circuit.output):
        raise AssertionError('Tseitin or shared-DAG model is inconsistent')
    lifted = decode(circuit, values, target)
    return {'status': 'SAT_LIFTED', 'label': label, **lifted}


def admit_unsat(circuit: Circuit, label: str, query: Path, proof: Path,
                solver_stdout: Path, exit_code: int, checker_binary: Path,
                checker_sha256: str, logs: Path, *, proof_cap: int,
                checker_wall_seconds: float, checker_rss_bytes: int) -> dict:
    target = _target(label)
    result = _admit_unsat(circuit, query, proof, solver_stdout, exit_code,
                          target, checker_binary, checker_sha256, logs,
                          proof_cap=proof_cap,
                          checker_wall_seconds=checker_wall_seconds,
                          checker_rss_bytes=checker_rss_bytes)
    return {'label': label, **result}
