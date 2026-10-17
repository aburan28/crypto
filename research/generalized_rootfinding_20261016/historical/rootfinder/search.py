"""Exact constrained search without Gröbner bases.

The search space is D**n, where D is the explicitly supplied common domain.
The meet-in-the-middle route applies only when the polynomial interaction
graph has at least two connected components. It is not a generic solver for
coupled systems, and the work limit counts assignments rather than field
operations or bytes of memory.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import product
from typing import Any, Sequence


class _PrimeArithmetic:
    """Integer-modulus adapter; callers are responsible for primality checks."""

    def __init__(self, order: int) -> None:
        self.order = order

    def add(self, left: int, right: int) -> int:
        return (left + right) % self.order

    def mul(self, left: int, right: int) -> int:
        return left * right % self.order

    def neg(self, value: int) -> int:
        return -value % self.order

    def pow(self, value: int, exponent: int) -> int:
        return pow(value, exponent, self.order)


def _arithmetic(p: Any) -> tuple[Any, bool]:
    if type(p) is int:
        if p < 2:
            raise ValueError("The prime-field modulus must be at least 2.")
        return _PrimeArithmetic(p), True
    if type(getattr(p, "order", None)) is not int or p.order < 2:
        raise ValueError("A field must expose an integer order of at least 2.")
    if not all(callable(getattr(p, name, None)) for name in ("add", "mul", "neg", "pow")):
        raise ValueError("A field must expose add, mul, neg, and pow operations.")
    return p, False


def _prepare_equations(
    equations: Sequence[Sequence[tuple[int, tuple[int, ...]]]],
    n_variables: int,
    field: Any,
    integer_modulus: bool,
) -> list[list[tuple[int, tuple[int, ...]]]]:
    prepared = []
    for equation in equations:
        coefficients: dict[tuple[int, ...], int] = {}
        for term in equation:
            if len(term) != 2:
                raise ValueError("Each term must contain a coefficient and exponent tuple.")
            coefficient, powers = term
            if type(coefficient) is not int:
                raise ValueError("Polynomial coefficients must be integers.")
            powers = tuple(powers)
            if len(powers) != n_variables or any(type(e) is not int or e < 0 for e in powers):
                raise ValueError("Exponent tuples must contain one nonnegative integer per variable.")
            if integer_modulus:
                coefficient %= field.order
            elif coefficient < 0 or coefficient >= field.order:
                raise ValueError("Polynomial coefficients must be canonical field elements.")
            coefficients[powers] = field.add(coefficients.get(powers, 0), coefficient)
        prepared.append([(c, powers) for powers, c in sorted(coefficients.items()) if c])
    return prepared


def _components(equations: Sequence[Sequence[tuple[int, tuple[int, ...]]]], n: int) -> list[tuple[int, ...]]:
    parents = list(range(n))

    def find(variable: int) -> int:
        while parents[variable] != variable:
            parents[variable] = parents[parents[variable]]
            variable = parents[variable]
        return variable

    for equation in equations:
        for _, powers in equation:
            support = [variable for variable, exponent in enumerate(powers) if exponent]
            if support:
                root = find(support[0])
                for variable in support[1:]:
                    parents[find(variable)] = root
    groups: dict[int, list[int]] = defaultdict(list)
    for variable in range(n):
        groups[find(variable)].append(variable)
    return sorted((tuple(group) for group in groups.values()), key=lambda group: group[0])


def _partition(components: Sequence[tuple[int, ...]]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left: list[int] = []
    right: list[int] = []
    for component in sorted(components, key=lambda group: (-len(group), group[0])):
        (left if len(left) <= len(right) else right).extend(component)
    if len(left) > len(right):
        left, right = right, left
    return tuple(sorted(left)), tuple(sorted(right))


def _compile(equations: Sequence[Sequence[tuple[int, tuple[int, ...]]]]) -> list[list[tuple[int, tuple[tuple[int, int], ...]]]]:
    return [
        [(coefficient, tuple((i, e) for i, e in enumerate(powers) if e)) for coefficient, powers in equation]
        for equation in equations
    ]


def _values(compiled: Sequence[Sequence[tuple[int, tuple[tuple[int, int], ...]]]], assignment: tuple[int, ...], field: Any) -> tuple[int, ...]:
    values = []
    for equation in compiled:
        value = 0
        for coefficient, support in equation:
            term = coefficient
            for variable, exponent in support:
                factor = assignment[variable] if exponent == 1 else field.pow(assignment[variable], exponent)
                term = field.mul(term, factor)
                if term == 0:
                    break
            value = field.add(value, term)
        values.append(value)
    return tuple(values)


def solve_on_domain(
    equations: Sequence[Sequence[tuple[int, tuple[int, ...]]]],
    n_variables: int,
    domain: Sequence[int],
    p: Any,
    *,
    method: str = "auto",
    max_work: int = 1_000_000,
    max_solutions: int = 10_000,
    max_table: int = 100_000,
) -> dict[str, Any]:
    """Find and exactly count all common roots in the explicit domain.

    ``p`` is a prime integer modulus (validated by the caller), or a field
    object exposing ``order``, ``add``, ``mul``, ``neg``, and ``pow``. Extension
    fields must use the latter interface: reduction modulo their order is
    not field arithmetic. ``max_work`` bounds enumerated assignments and
    ``max_table`` bounds stored assignments in the meet-in-the-middle table.
    A truncated result still reports the exact total number of solutions.
    """
    if type(n_variables) is not int or n_variables < 1:
        raise ValueError("n_variables must be a positive integer.")
    if method not in ("auto", "exhaustive", "mitm"):
        raise ValueError("method must be 'auto', 'exhaustive', or 'mitm'.")
    if type(max_work) is not int or max_work < 0:
        raise ValueError("max_work must be a nonnegative integer.")
    if type(max_solutions) is not int or max_solutions < 0:
        raise ValueError("max_solutions must be a nonnegative integer.")
    if type(max_table) is not int or max_table < 0:
        raise ValueError("max_table must be a nonnegative integer.")
    field, integer_modulus = _arithmetic(p)
    domain = tuple(domain)
    if any(type(value) is not int or value < 0 or value >= field.order for value in domain):
        raise ValueError("Domain values must be canonical field elements.")
    if len(set(domain)) != len(domain):
        raise ValueError("Domain values must be distinct.")
    domain = tuple(sorted(domain))
    equations = _prepare_equations(equations, n_variables, field, integer_modulus)
    components = _components(equations, n_variables)
    split = _partition(components) if len(components) >= 2 else None
    exhaustive_work = len(domain) ** n_variables
    mitm_work = sum(len(domain) ** len(part) for part in split) if split else None
    mitm_table = len(domain) ** len(split[0]) if split else 0
    if method == "mitm" and split is None and domain:
        raise ValueError("Meet-in-the-middle requires at least two independent variable blocks.")
    if method == "auto":
        method = "mitm" if mitm_work is not None and mitm_work <= exhaustive_work else "exhaustive"
        if method == "mitm" and mitm_table > max_table and exhaustive_work <= max_work:
            method = "exhaustive"
    work_estimate = mitm_work if method == "mitm" and split is not None else exhaustive_work
    result: dict[str, Any] = {
        "method": method,
        "status": "complete",
        "solutions": [],
        "solution_count_exact": 0,
        "work_estimate": work_estimate,
        "table_assignment_estimate": mitm_table if method == "mitm" else 0,
        "assignments_visited": 0,
        "reason": "Every assignment in the constrained domain has been accounted for.",
    }
    if not domain:
        result["reason"] = "The common domain is empty, so there are no candidate assignments."
        return result
    if work_estimate > max_work:
        result.update(
            status="budget_exceeded",
            solution_count_exact=None,
            reason=f"The estimate of {work_estimate} assignments exceeds max_work={max_work}; search was not started.",
        )
        return result
    if method == "mitm" and mitm_table > max_table:
        result.update(
            status="budget_exceeded",
            solution_count_exact=None,
            reason=f"The table estimate of {mitm_table} stored assignments exceeds max_table={max_table}; search was not started.",
        )
        return result

    original = _compile(equations)
    solutions: list[list[int]] = result["solutions"]
    count = 0
    if method == "exhaustive":
        for assignment in product(domain, repeat=n_variables):
            result["assignments_visited"] += 1
            if not any(_values(original, assignment, field)):
                count += 1
                if len(solutions) < max_solutions:
                    solutions.append(list(assignment))
    else:
        assert split is not None
        left_indices, right_indices = split
        left_set = set(left_indices)
        left_equations: list[list[tuple[int, tuple[int, ...]]]] = []
        right_equations: list[list[tuple[int, tuple[int, ...]]]] = []
        for equation in equations:
            left_terms = []
            right_terms = []
            for coefficient, powers in equation:
                support = [i for i, e in enumerate(powers) if e]
                if not support or support[0] in left_set:
                    left_terms.append((coefficient, tuple(powers[i] for i in left_indices)))
                else:
                    right_terms.append((coefficient, tuple(powers[i] for i in right_indices)))
            left_equations.append(left_terms)
            right_equations.append(right_terms)
        left_compiled, right_compiled = _compile(left_equations), _compile(right_equations)
        buckets: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
        for left in product(domain, repeat=len(left_indices)):
            result["assignments_visited"] += 1
            buckets[_values(left_compiled, left, field)].append(left)
        for right in product(domain, repeat=len(right_indices)):
            result["assignments_visited"] += 1
            target = tuple(field.neg(value) for value in _values(right_compiled, right, field))
            matches = buckets.get(target, ())
            count += len(matches)
            for left in matches[:max_solutions - len(solutions)]:
                assignment = [0] * n_variables
                for variable, value in zip(left_indices, left):
                    assignment[variable] = value
                for variable, value in zip(right_indices, right):
                    assignment[variable] = value
                if any(_values(original, tuple(assignment), field)):
                    raise RuntimeError("A meet-in-the-middle candidate failed exact verification.")
                solutions.append(assignment)
    solutions.sort()
    result["solution_count_exact"] = count
    if count > max_solutions:
        result.update(
            status="output_truncated",
            reason=f"All candidates were accounted for, but only {max_solutions} of {count} solutions were retained.",
        )
    return result
