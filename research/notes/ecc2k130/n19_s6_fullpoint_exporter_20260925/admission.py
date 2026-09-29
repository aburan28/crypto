#!/usr/bin/env python3
"""Independent full-point and exact-DIMACS admission for the S6 exporter."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

from s6 import Circuit, EDGES, ROLES, ReferenceField

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "symbolic_dag_dimacs_gate_20260925"))
from verify import (evaluate_cnf, extend_dag, parse_cnf,  # type: ignore  # noqa: E402
                    parse_solver_output, verify_gate_encoding)


def _from_bits(bits: list[int]) -> int:
    return sum(bit << i for i, bit in enumerate(bits))


def _point(inputs: dict[str, int], role: str, n: int) -> tuple[int, int, int]:
    return (inputs[f"{role}_o"],
            _from_bits([inputs[f"{role}_x_{i}"] for i in range(n)]),
            _from_bits([inputs[f"{role}_y_{i}"] for i in range(n)]))


def decode(circuit: Circuit, values: list[int],
           target: tuple[int, int, int]) -> dict:
    """Require a complete Boolean model and independently re-add full points."""
    dag = circuit.dag
    if len(values) != len(dag.nodes) or any(bit not in (0, 1) for bit in values):
        raise ValueError("incomplete/non-Boolean DIMACS model")
    inputs = {dag.names[a]: values[node_id]
              for node_id, (op, a, _) in enumerate(dag.nodes) if op == "var"}
    if len(inputs) != len(dag.names):
        raise AssertionError("primary variable name drift")
    ref = ReferenceField(circuit.n, circuit.modulus)
    roles = {role: _point(inputs, role, circuit.n) for role in ROLES}
    if any(not ref.valid(p) for p in roles.values()) or roles["SUM"] != target:
        raise AssertionError("noncanonical, off-curve, or wrong target point")
    indices = []
    for i in range(6):
        chosen = [j for j in range(7) if inputs[f"A{i}_{j}"]]
        if len(chosen) != 1:
            raise AssertionError("selector group is not exactly one-hot")
        j = chosen[0]
        if roles[f"F{i}"] != circuit.factors[i][j]:
            raise AssertionError("selector-to-full-point link failed")
        indices.append(j)
    branches = []
    for i, (left, right, out) in enumerate(EDGES):
        expected, algebraic_slope, branch = ref.add_slope(roles[left], roles[right])
        if roles[out] != expected:
            raise AssertionError("independent group-law chain failed")
        encoded_slope = _from_bits([inputs[f"L{i}_{bit}"] for bit in range(circuit.n)])
        if branch in ("generic", "double") and encoded_slope != algebraic_slope:
            raise AssertionError("slope is inconsistent with group law")
        branches.append(branch)
    # The circuit/clauses are checked separately.  This equality ensures the
    # fixed target is the exact six-point sum, not just a cofactor projection.
    return {"indices": indices, "target": list(target), "branches": branches,
            "factor_points": [list(roles[f"F{i}"]) for i in range(6)],
            "prefix_points": [list(roles[f"S{i}"]) for i in range(2, 6)]}


def check_cnf(circuit: Circuit, path: Path, units: list[int]) -> dict:
    variables, clause_count, clauses = parse_cnf(path)
    if variables != len(circuit.dag.nodes):
        raise AssertionError("DIMACS variable count drift")
    gate = verify_gate_encoding(circuit.dag, clauses, extra_units=units)
    at = gate["cursor_before_output"]
    if clauses[at:] != [(circuit.output + 1,)] + [(lit,) for lit in units]:
        raise AssertionError("output or target pin clause drift")
    return {"variables": variables, "clauses": clause_count,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def check_full_model(circuit: Circuit, clauses: list[tuple[int, ...]],
                     values: list[int], target: tuple[int, int, int]) -> dict:
    if not evaluate_cnf(clauses, values):
        raise AssertionError("SAT assignment does not satisfy exact CNF")
    # A second DAG evaluation checks that every Tseitin auxiliary agrees with
    # its gate, even if a different exporter wrote the query.
    from dag import PackedModel
    primary = [None] * len(circuit.dag.names)
    for node_id, (op, a, _) in enumerate(circuit.dag.nodes):
        if op == "var":
            primary[a] = values[node_id]
    packed = PackedModel.from_bits(primary)
    if extend_dag(circuit.dag, packed) != values:
        raise AssertionError("incorrect auxiliary variable in SAT model")
    if not circuit.dag.evaluate(packed, circuit.output):
        raise AssertionError("Boolean full-point circuit rejects SAT model")
    return decode(circuit, values, target)


def admit_sat(circuit: Circuit, query: Path, solver_stdout: Path,
              exit_code: int, target: tuple[int, int, int]) -> dict:
    units = circuit.target_units(target)  # full RHS, not torsion alone
    check_cnf(circuit, query, units)
    variables, _, clauses = parse_cnf(query)
    status, values = parse_solver_output(solver_stdout.read_text(), variables, exit_code)
    if status != "SAT" or values is None:
        raise AssertionError("not an exact complete SAT model")
    return check_full_model(circuit, clauses, values, target)


def admit_unsat(circuit: Circuit, query: Path, proof: Path, solver_stdout: Path,
                exit_code: int, target: tuple[int, int, int],
                checker_binary: Path, checker_sha256: str, logs: Path,
                *, proof_cap: int, checker_wall_seconds: float,
                checker_rss_bytes: int) -> dict:
    """Re-run a pinned external checker on this exact, fully audited query."""
    units = circuit.target_units(target)
    check_cnf(circuit, query, units)
    variables, _, _ = parse_cnf(query, keep_clauses=False)
    status, values = parse_solver_output(solver_stdout.read_text(), variables, exit_code)
    if status != "UNSAT" or values is not None:
        raise AssertionError("solver did not emit exact UNSAT status and exit 20")
    if not proof.is_file() or not 0 < proof.stat().st_size <= proof_cap:
        raise AssertionError("missing, empty, or over-cap proof")
    with proof.open("rb") as stream:
        while chunk := stream.read(1 << 20):
            if b"\x00" in chunk or any(byte > 127 for byte in chunk):
                raise AssertionError("proof is not text DRAT")
    if not checker_binary.is_file() or hashlib.sha256(checker_binary.read_bytes()).hexdigest() != checker_sha256:
        raise AssertionError("external DRAT checker binary SHA-256 mismatch")
    sys.path.insert(0, str(HERE.parent / "symbolic_dag_dimacs_gate_20260925"))
    from bounded import run_child
    logs.mkdir(parents=True, exist_ok=False)
    receipt = run_child([str(checker_binary.resolve()), str(query.resolve()), str(proof.resolve())],
                        cwd=HERE, stdout=logs / 'checker.stdout.txt',
                        stderr=logs / 'checker.stderr.txt',
                        wall_cap=checker_wall_seconds, rss_cap=checker_rss_bytes)
    if receipt.get('exit_code') != 0 or receipt.get('stop_reason') is not None:
        raise AssertionError("external DRAT checker did not succeed")
    return {"status": "PROVED_UNSAT", "query_sha256": hashlib.sha256(query.read_bytes()).hexdigest(),
            "proof_sha256": hashlib.sha256(proof.read_bytes()).hexdigest(),
            "checker_binary_sha256": checker_sha256, "checker": receipt}
