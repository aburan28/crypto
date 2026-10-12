#!/usr/bin/env python3
"""Deterministic, bounded Tseitin DIMACS export of the parent XOR/AND DAG."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'symbolic_dag_fullpoint_20260925'
sys.path.append(str(PARENT))
from dag import Dag, PackedModel, Relation, build_relation  # noqa: E402

DIMACS_MAX_VAR = (1 << 31) - 1


class ExportBudgetExceeded(RuntimeError):
    pass


def expected_clauses(dag: Dag, units: list[int]) -> int:
    counts = dag.counts()
    return 2 + 4 * counts['xor'] + 3 * counts['and'] + 1 + len(units)


def _gate_clauses(op: str, a: int, b: int, z: int) -> tuple[tuple[int, ...], ...]:
    """Input/output arguments are already positive one-based DIMACS IDs."""
    if op == 'xor':
        return ((-a, -b, -z), (a, b, -z), (a, -b, z), (-a, b, z))
    if op == 'and':
        return ((-a, -b, z), (a, -z), (b, -z))
    raise ValueError('not a binary gate')


def write_cnf(relation: Relation, path: Path, *, units: list[int] | None = None,
              byte_cap: int) -> dict:
    """Write canonical DIMACS. Keep partial bytes on budget failure for audit."""
    units = list(units or [])
    dag = relation.dag
    var_count = len(dag.nodes)
    if var_count > DIMACS_MAX_VAR:
        raise ExportBudgetExceeded('DIMACS signed-variable-ID limit')
    if any(not isinstance(lit, int) or lit == 0 or abs(lit) > var_count for lit in units):
        raise ValueError('unit literal outside exported DAG')
    clauses = expected_clauses(dag, units)
    written = 0
    digest = hashlib.sha256()

    def emit(stream, text: str) -> None:
        nonlocal written
        data = text.encode('ascii')
        if written + len(data) > byte_cap:
            raise ExportBudgetExceeded(f'DIMACS byte cap {byte_cap} exceeded')
        stream.write(data)
        digest.update(data)
        written += len(data)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('wb') as stream:
        emit(stream, f'p cnf {var_count} {clauses}\n')
        emit(stream, '-1 0\n')
        emit(stream, '2 0\n')
        for node_id, (op, a, b) in enumerate(dag.nodes[2:], 2):
            if op == 'var':
                continue
            for clause in _gate_clauses(op, a + 1, b + 1, node_id + 1):
                emit(stream, ' '.join(map(str, clause)) + ' 0\n')
        emit(stream, f'{relation.output + 1} 0\n')
        for lit in units:
            emit(stream, f'{lit} 0\n')
    return {'variables': var_count, 'clauses': clauses, 'bytes': written,
            'sha256': digest.hexdigest(), 'dag': dag.counts(),
            'relation_output_literal': relation.output + 1}


def fixed_point_units(relation: Relation, p: tuple[int, int, int],
                      q: tuple[int, int, int], r: tuple[int, int, int]) -> list[int]:
    """Fix P,Q,R only; leave every λ bit existential for the solver."""
    n = relation.field.n
    model = relation.model(p, q, r, 0)
    assert len(relation.dag.names) == 7 * n + 3
    fixed_bits = 3 * (2 * n + 1)
    units = []
    for node_id, (op, var_index, _) in enumerate(relation.dag.nodes):
        if op == 'var' and var_index < fixed_bits:
            units.append((node_id + 1) if model.bit(var_index) else -(node_id + 1))
    assert len(units) == fixed_bits
    return units
