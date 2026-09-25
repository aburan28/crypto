#!/usr/bin/env python3
"""Independent DIMACS parser, gate truth-table checker, and model lifter."""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'symbolic_dag_fullpoint_20260925'
sys.path.insert(0, str(PARENT))
from dag import PackedModel, build_relation  # noqa: E402

O = (1, 0, 0)


def parse_cnf(path: Path, *, keep_clauses: bool = True) -> tuple[int, int, list[tuple[int, ...]]]:
    variables = clauses_declared = None
    clauses: list[tuple[int, ...]] = []
    count = 0
    with path.open('rt', encoding='ascii') as stream:
        for line_no, line in enumerate(stream, 1):
            parts = line.split()
            if not parts:
                raise ValueError(f'empty DIMACS line {line_no}')
            if variables is None:
                if len(parts) != 4 or parts[:2] != ['p', 'cnf']:
                    raise ValueError('missing exact DIMACS header')
                variables, clauses_declared = map(int, parts[2:])
                if not 1 <= variables <= (1 << 31) - 1 or clauses_declared < 1:
                    raise ValueError('invalid DIMACS dimensions')
                continue
            literals = [int(part) for part in parts]
            if len(literals) < 2 or literals[-1] != 0 or 0 in literals[:-1]:
                raise ValueError(f'malformed clause on line {line_no}')
            clause = tuple(literals[:-1])
            if any(abs(lit) > variables for lit in clause):
                raise ValueError(f'literal outside header on line {line_no}')
            if keep_clauses:
                clauses.append(clause)
            count += 1
    if variables is None or count != clauses_declared:
        raise ValueError('DIMACS clause count mismatch')
    return variables, count, clauses


def evaluate_cnf(clauses: list[tuple[int, ...]], values: list[int]) -> bool:
    if any(v not in (0, 1) for v in values):
        raise ValueError('incomplete or non-Boolean CNF model')
    return all(any((values[abs(lit) - 1] == 1) == (lit > 0) for lit in clause)
               for clause in clauses)


def extend_dag(dag, model: PackedModel) -> list[int]:
    """Independent direct node evaluation; no exporter clause generator used."""
    model.validate(len(dag.names))
    answer = [0] * len(dag.nodes)
    answer[1] = 1
    for i in range(2, len(dag.nodes)):
        op, a, b = dag.nodes[i]
        if op == 'var':
            answer[i] = model.bit(a)
        elif op == 'and':
            answer[i] = int(answer[a] == 1 and answer[b] == 1)
        elif op == 'xor':
            answer[i] = (answer[a] + answer[b]) % 2
        else:
            raise AssertionError('unknown DAG node')
    return answer


def _local_truth(clauses: list[tuple[int, ...]], a: int, b: int, z: int,
                 expected_op: str) -> None:
    for av in (0, 1):
        for bv in (0, 1):
            for zv in (0, 1):
                local = {a: av, b: bv, z: zv}
                holds = all(any((local[abs(lit)] == 1) == (lit > 0)
                                for lit in clause) for clause in clauses)
                expected = (zv == (av ^ bv if expected_op == 'xor' else av & bv))
                if holds != expected:
                    raise AssertionError(f'{expected_op} local Tseitin truth table failed')


def verify_gate_encoding(dag, clauses: list[tuple[int, ...]],
                         *, extra_units: list[int] | None = None) -> dict:
    extra_units = list(extra_units or [])
    cursor = 0
    assert clauses[cursor] == (-1,)
    cursor += 1
    assert clauses[cursor] == (2,)
    cursor += 1
    binary = 0
    for i, (op, a, b) in enumerate(dag.nodes[2:], 2):
        if op == 'var':
            continue
        width = 4 if op == 'xor' else 3 if op == 'and' else None
        if width is None:
            raise AssertionError('unexpected DAG gate')
        gate_clauses = clauses[cursor:cursor + width]
        assert len(gate_clauses) == width
        assert all(set(abs(lit) for lit in clause) <= {a + 1, b + 1, i + 1}
                   for clause in gate_clauses)
        _local_truth(gate_clauses, a + 1, b + 1, i + 1, op)
        cursor += width
        binary += 1
    # The asserted output literal is provided by the caller's relation, not
    # necessarily the final node in a hash-consed DAG.
    return {'cursor_before_output': cursor, 'binary_gates': binary,
            'extra_units': extra_units}


def check_relation_cnf(relation, clauses: list[tuple[int, ...]],
                       units: list[int] | None = None) -> None:
    units = list(units or [])
    gate = verify_gate_encoding(relation.dag, clauses, extra_units=units)
    cursor = gate['cursor_before_output']
    assert clauses[cursor] == (relation.output + 1,)
    cursor += 1
    assert clauses[cursor:] == [(lit,) for lit in units]


def parse_solver_output(output: str, variables: int, exit_code: int) -> tuple[str, list[int] | None]:
    states = []
    model_literals = []
    for line in output.splitlines():
        if line.startswith('s '):
            states.append(line[2:].strip())
        elif line.startswith('v '):
            tokens = line[2:].split()
            if not tokens or tokens[-1] != '0':
                raise ValueError('unterminated SAT model line')
            model_literals.extend(int(token) for token in tokens[:-1])
        elif line.startswith('c ') or line == 'c' or not line.strip():
            continue
        else:
            raise ValueError('unexpected solver output line')
    if states == ['UNSATISFIABLE'] and exit_code == 20 and not model_literals:
        return 'UNSAT', None
    if states != ['SATISFIABLE'] or exit_code != 10:
        raise ValueError('solver status/exit mismatch')
    assigned = [None] * variables
    for lit in model_literals:
        if lit == 0 or abs(lit) > variables:
            raise ValueError('SAT model literal out of range')
        i = abs(lit) - 1
        if assigned[i] is not None:
            raise ValueError('duplicate or conflicting SAT model literal')
        assigned[i] = int(lit > 0)
    if any(bit is None for bit in assigned):
        raise ValueError('partial SAT model')
    return 'SAT', assigned


class ReferenceField:
    """Polynomial product/remainder and Euclid inverse, independent of DAG."""

    def __init__(self, n: int, modulus: int):
        self.n, self.modulus = n, modulus

    def mul(self, a: int, b: int) -> int:
        product = 0
        for i in range(self.n):
            if a >> i & 1:
                product ^= b << i
        while product.bit_length() > self.n:
            product ^= self.modulus << (product.bit_length() - self.n - 1)
        return product

    def inv(self, a: int) -> int:
        if not a:
            raise ZeroDivisionError
        for candidate in range(1, 1 << self.n):
            if self.mul(a, candidate) == 1:
                return candidate
        raise AssertionError('not a field')

    def point_valid(self, p: tuple[int, int, int]) -> bool:
        o, x, y = p
        if o:
            return p == O
        return self.mul(y, y) ^ self.mul(x, y) == self.mul(self.mul(x, x), x) ^ 1

    def add(self, p: tuple[int, int, int], q: tuple[int, int, int]) -> tuple[int, int, int]:
        if p == O:
            return q
        if q == O:
            return p
        _, x, y = p
        _, u, v = q
        if x == u and y ^ v == x:
            return O
        if x == u:
            slope = x ^ self.mul(y, self.inv(x))
            rx = self.mul(slope, slope) ^ slope
            ry = self.mul(x, x) ^ self.mul(slope ^ 1, rx)
        else:
            slope = self.mul(y ^ v, self.inv(x ^ u))
            rx = self.mul(slope, slope) ^ slope ^ x ^ u
            ry = self.mul(slope, x ^ rx) ^ rx ^ y
        return (0, rx, ry)


def lift_solver_model(relation, values: list[int], p, q, r,
                      field: ReferenceField) -> dict:
    if len(values) != len(relation.dag.nodes):
        raise ValueError('solver model does not cover the exact CNF')
    input_values = [None] * len(relation.dag.names)
    for node_id, (op, index, _) in enumerate(relation.dag.nodes):
        if op == 'var':
            input_values[index] = values[node_id]
    packed = PackedModel.from_bits(input_values)
    n = relation.field.n
    def point(offset):
        bits = input_values[offset:offset + 2 * n + 1]
        o = bits[0]
        x = sum(bit << i for i, bit in enumerate(bits[1:n + 1]))
        y = sum(bit << i for i, bit in enumerate(bits[n + 1:]))
        return (o, x, y)
    recovered = (point(0), point(2 * n + 1), point(4 * n + 2))
    slope = sum(bit << i for i, bit in enumerate(input_values[6 * n + 3:]))
    if recovered != (p, q, r):
        raise AssertionError('solver point model does not match fixed query')
    if not relation.dag.evaluate(packed, relation.output):
        raise AssertionError('solver model does not satisfy point relation')
    if not all(field.point_valid(item) for item in recovered):
        raise AssertionError('lifted solver model has invalid full point')
    if field.add(p, q) != r:
        raise AssertionError('lifted solver model has wrong group sum')
    return {'p': p, 'q': q, 'r': r, 'lambda': slope}


def verify_parent_rows(cnfs: dict[int, Path]) -> dict:
    archive = PARENT / 'evidence/producer/rows.jsonl.gz'
    expected_sha = 'a2b760b72e477cb9ab4d8409e71102cd6a2a09db0b8a5e8f5ad170d1a1102da9'
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == expected_sha
    counts = []
    with gzip.open(archive, 'rt') as stream:
        for n, modulus in ((2, 0x7), (3, 0xB)):
            relation = build_relation(n, modulus)
            variables, clause_count, clauses = parse_cnf(cnfs[n])
            assert variables == len(relation.dag.nodes)
            check_relation_cnf(relation, clauses)
            rows = evaluations = 0
            target_rows = 512 if n == 2 else 64
            for _ in range(target_rows):
                row = json.loads(stream.readline())
                assert row['n'] == n
                p, q, r = (tuple(row[key]) for key in ('p', 'q', 'r'))
                accepted = 0
                for lam in range(1 << n):
                    packed = relation.model(p, q, r, lam)
                    nodes = extend_dag(relation.dag, packed)
                    truth = evaluate_cnf(clauses, nodes)
                    assert truth == relation.dag.evaluate(packed, relation.output)
                    accepted += truth
                    evaluations += 1
                assert accepted == row['witnesses']
                rows += 1
            ref = ReferenceField(n, modulus)
            universe = [(o, x, y) for o in (0, 1) for x in range(1 << n)
                        for y in range(1 << n)]
            invalid = [point for point in universe if not ref.point_valid(point)]
            invalid_evaluations = 0
            for bad in invalid:
                for p, q, r in ((bad, O, O), (O, bad, O), (O, O, bad)):
                    for lam in range(1 << n):
                        packed = relation.model(p, q, r, lam)
                        truth = evaluate_cnf(clauses, extend_dag(relation.dag, packed))
                        assert not truth
                        invalid_evaluations += 1
            archived = json.loads((PARENT / 'evidence/producer/result.json').read_text())
            archived_field = next(item for item in archived['fields'] if item['n'] == n)
            assert invalid_evaluations == archived_field['invalid_model_evaluations']
            counts.append({'n': n, 'rows': rows, 'assignments': evaluations,
                           'invalid_assignments': invalid_evaluations,
                           'variables': variables, 'clauses': clause_count})
        assert stream.readline() == ''
    return {'parent_rows_sha256': expected_sha, 'fields': counts}
