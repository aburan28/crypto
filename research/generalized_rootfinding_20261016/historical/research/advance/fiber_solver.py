"""Exact finite-domain residual-span decision DAG, without a Gröbner basis.

Research reference. Specialization, coefficient-space Gaussian elimination,
finite-domain polynomial remainders and terminal univariate gcd are the only
algebraic operations; no equation is multiplied by a monomial. Worst-case
complexity remains exponential. See GENERAL_METHOD.md for the precise scope.
"""

from dataclasses import dataclass
from itertools import product

from rootfinder import polynomial as poly


class _Limit(Exception):
    pass


class _Budget:
    def __init__(self, max_states, max_transitions, max_arithmetic, max_cached_terms):
        self.limits = {"states": max_states, "transitions": max_transitions,
                       "arithmetic_primitives": max_arithmetic,
                       "cached_coefficient_entries": max_cached_terms}
        self.counts = {name: 0 for name in self.limits}

    def charge(self, name, amount=1):
        if self.counts[name] + amount > self.limits[name]:
            raise _Limit(f"{name} exceeds configured budget {self.limits[name]}")
        self.counts[name] += amount


class _Field:
    def __init__(self, field, budget):
        self.field, self.budget = field, budget

    def __getattr__(self, name):
        if name.startswith("_raw_"):
            raise AttributeError(name)
        return getattr(self.field, name)

    def _op(self, name, *values):
        self.budget.charge("arithmetic_primitives")
        return getattr(self.field, name)(*values)

    def add(self, a, b): return self._op("add", a, b)
    def sub(self, a, b): return self._op("sub", a, b)
    def neg(self, a): return self._op("neg", a)
    def mul(self, a, b): return self._op("mul", a, b)
    def pow(self, a, b): return self._op("pow", a, b)
    def inv(self, a): return self._op("inv", a)


@dataclass
class _Node:
    kind: str
    variables: tuple
    count: int
    # branch: (selected index, ((value, child), ...)); free: (indices, child)
    # fiber: squarefree common univariate gcd; accept/reject: None.
    data: object = None


def _nonnegative(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _canonical(equations, width, field, max_terms):
    """Canonical coefficient-row span; a constant row proves inconsistency."""
    rows = []
    for equation in equations:
        row = {}
        for coefficient, powers in equation:
            if coefficient:
                row[powers] = field.add(row.get(powers, 0), coefficient)
                if not row[powers]:
                    del row[powers]
        if row:
            rows.append(row)
    columns = sorted({powers for row in rows for powers in row})
    if sum(map(len, rows)) > max_terms:
        raise _Limit(f"residual monomial storage exceeds max_terms={max_terms}")
    rank = 0
    for column in columns:
        pivot = next((i for i in range(rank, len(rows)) if rows[i].get(column)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = field.inv(rows[rank][column])
        rows[rank] = {powers: field.mul(value, inverse)
                      for powers, value in rows[rank].items()}
        for i in range(len(rows)):
            if i == rank or not rows[i].get(column):
                continue
            factor = rows[i][column]
            for powers, value in rows[rank].items():
                updated = field.sub(rows[i].get(powers, 0), field.mul(factor, value))
                if updated:
                    rows[i][powers] = updated
                else:
                    rows[i].pop(powers, None)
        if sum(map(len, rows)) > max_terms:
            raise _Limit(f"row reduction exceeds max_terms={max_terms}")
        rank += 1
        if rank == len(rows):
            break
    zero = (0,) * width
    if any(len(row) == 1 and zero in row for row in rows[:rank]):
        return None
    return tuple(tuple((coefficient, powers) for powers, coefficient in sorted(row.items()))
                 for row in rows[:rank])


def _domain_polynomial(domain, field):
    result = (1,)
    for value in domain:
        result = poly.mul(result, (field.neg(value), 1), field)
    return result


def _normalize(equations, n, domain_poly, field, max_terms):
    """Take the normal form modulo each separate domain polynomial."""
    degree = len(domain_poly) - 1
    remainders = {}
    normalized = []
    for equation in equations:
        terms = {}
        for coefficient, raw_powers in equation:
            field.validate(coefficient)
            powers = tuple(raw_powers)
            if len(powers) != n or any(type(e) is not int or e < 0 for e in powers):
                raise ValueError("one nonnegative integer exponent per variable is required")
            if any(e.bit_length() > 256 for e in powers):
                raise ValueError("exponent bit length exceeds the research implementation limit 256")
            if not coefficient:
                continue
            expansion = {(): coefficient}
            for exponent in powers:
                if exponent < degree:
                    remainder = ((exponent, 1),)
                else:
                    if exponent not in remainders:
                        coefficients = poly.powmod_poly((0, 1), exponent, domain_poly, field)
                        remainders[exponent] = tuple((i, c) for i, c in enumerate(coefficients) if c)
                    remainder = remainders[exponent]
                next_expansion = {}
                for prefix, value in expansion.items():
                    for power, factor in remainder:
                        next_expansion[prefix + (power,)] = field.mul(value, factor)
                expansion = next_expansion
                if len(expansion) > max_terms:
                    raise _Limit(f"domain normal form exceeds max_terms={max_terms}")
            for exponents, value in expansion.items():
                updated = field.add(terms.get(exponents, 0), value)
                if updated:
                    terms[exponents] = updated
                else:
                    terms.pop(exponents, None)
            if len(terms) > max_terms:
                raise _Limit(f"domain normal form exceeds max_terms={max_terms}")
        normalized.append(tuple((c, e) for e, c in sorted(terms.items())))
    return normalized


def _specialize(equations, variable, value, field):
    result = []
    powers = {}
    for equation in equations:
        merged = {}
        for coefficient, exponents in equation:
            exponent = exponents[variable]
            if exponent not in powers:
                powers[exponent] = field.pow(value, exponent)
            coefficient = field.mul(coefficient, powers[exponent])
            reduced = exponents[:variable] + exponents[variable + 1:]
            merged[reduced] = field.add(merged.get(reduced, 0), coefficient)
        result.append(tuple((c, e) for e, c in sorted(merged.items()) if c))
    return result


def _univariate(equation, variable):
    coefficients = {}
    for coefficient, powers in equation:
        if any(power for i, power in enumerate(powers) if i != variable):
            return None
        coefficients[powers[variable]] = coefficient
    if not coefficients:
        return ()
    return tuple(coefficients.get(i, 0) for i in range(max(coefficients) + 1))


def _root_values(polynomial, domain, field, metrics):
    if len(polynomial) <= 1:
        return domain if not polynomial else ()
    if len(polynomial) == 2:
        root = field.neg(field.mul(polynomial[0], field.inv(polynomial[1])))
        return (root,) if root in domain else ()
    roots = []
    for value in domain:
        metrics["domain_root_tests"] += 1
        if poly.evaluate_poly(polynomial, value, field) == 0:
            roots.append(value)
    if len(roots) != len(polynomial) - 1:
        raise ArithmeticError("domain gcd degree and extracted distinct roots disagree")
    return tuple(roots)


def solve_fiber_dag(equations, n_variables, domain, field, *, variable_order=None,
                    max_states=100_000, max_transitions=1_000_000,
                    max_arithmetic=10_000_000, max_terms=50_000,
                    max_cached_terms=500_000, max_solutions=10_000):
    """Return exact all-root count and bounded witnesses, or explicit refusal.

    ``equations`` consists of coefficient/exponent-tuple sparse polynomials.
    ``domain`` is a common explicit set of canonical field elements. Supports
    arbitrary total degree and any exact rootfinder finite-field backend.
    Budgets count states, transitions, scalar arithmetic primitive calls and
    per-state and cached stored terms; they are not wall-clock speed guarantees.
    """
    if type(n_variables) is not int or n_variables < 1 or n_variables > 512:
        raise ValueError("n_variables must be between 1 and 512")
    for name, value in (("max_states", max_states), ("max_transitions", max_transitions),
                        ("max_arithmetic", max_arithmetic), ("max_terms", max_terms),
                        ("max_cached_terms", max_cached_terms), ("max_solutions", max_solutions)):
        _nonnegative(value, name)
    domain = tuple(domain)
    if len(set(domain)) != len(domain):
        raise ValueError("domain must contain distinct elements")
    for value in domain:
        field.validate(value)
    domain = tuple(sorted(domain))
    if variable_order is None:
        variable_order = tuple(range(n_variables))
    else:
        variable_order = tuple(variable_order)
    if (len(variable_order) != n_variables or set(variable_order) != set(range(n_variables))
            or any(type(i) is not int for i in variable_order)):
        raise ValueError("variable_order must be a permutation of variable indices")
    equations = tuple(tuple((c, tuple(e)) for c, e in equation) for equation in equations)
    result = {"method": "finite_domain_residual_span_fiber_dag", "status": "complete",
              "solutions": [], "solution_count_exact": 0, "output_complete": True,
              "count_complete": True, "exhaustive_assignments": len(domain) ** n_variables,
              "variable_order": variable_order, "metrics": {},
              "reason": "All constrained tuples accounted for by exact residual-span recursion."}
    budget = _Budget(max_states, max_transitions, max_arithmetic, max_cached_terms)
    tracked = _Field(field, budget)
    metrics = {"cache_hits": 0, "univariate_fibers": 0, "unary_restrictions": 0,
               "free_variable_collapses": 0, "constant_prunes": 0, "domain_root_tests": 0,
               "domain_root_cache_hits": 0}
    result["metrics"] = metrics
    if not domain:
        metrics.update(budget.counts)
        return result
    memo = {}
    root_values_cache = {}

    def roots_for(polynomial):
        if polynomial in root_values_cache:
            metrics["domain_root_cache_hits"] += 1
            return root_values_cache[polynomial]
        roots = _root_values(polynomial, domain, tracked, metrics)
        root_values_cache[polynomial] = roots
        return roots

    priority = {variable: i for i, variable in enumerate(variable_order)}

    def build(equations, variables):
        canonical = _canonical(equations, len(variables), tracked, max_terms)
        key = (variables, canonical)
        if key in memo:
            metrics["cache_hits"] += 1
            return memo[key]
        budget.charge("states")
        budget.charge("cached_coefficient_entries", sum(map(len, canonical)) if canonical else 0)
        if canonical is None:
            metrics["constant_prunes"] += 1
            node = _Node("reject", variables, 0)
        elif not canonical:
            node = _Node("accept", variables, len(domain) ** len(variables))
        else:
            support = {i for equation in canonical for _, powers in equation
                       for i, exponent in enumerate(powers) if exponent}
            unused = tuple(i for i in range(len(variables)) if i not in support)
            if unused:
                used = tuple(i for i in range(len(variables)) if i in support)
                reduced = tuple(tuple((c, tuple(e[i] for i in used)) for c, e in equation)
                                for equation in canonical)
                child = build(reduced, tuple(variables[i] for i in used))
                metrics["free_variable_collapses"] += len(unused)
                node = _Node("free", variables, child.count * len(domain) ** len(unused),
                             (unused, child))
            elif len(variables) == 1:
                common = domain_poly
                for equation in canonical:
                    common = poly.gcd_poly(common, _univariate(equation, 0), tracked)
                    if len(common) == 1:
                        break
                metrics["univariate_fibers"] += 1
                node = _Node("fiber", variables, len(common) - 1, common)
            else:
                unary = {}
                for equation in canonical:
                    active = {i for _, powers in equation for i, exponent in enumerate(powers) if exponent}
                    if len(active) == 1:
                        index = next(iter(active))
                        common = unary.get(index, domain_poly)
                        unary[index] = poly.gcd_poly(common, _univariate(equation, index), tracked)
                if unary:
                    variable = min(unary, key=lambda i: (len(unary[i]), priority[variables[i]]))
                    candidates = roots_for(unary[variable])
                    if len(candidates) < len(domain):
                        metrics["unary_restrictions"] += 1
                else:
                    variable = min(range(len(variables)), key=lambda i: priority[variables[i]])
                    candidates = domain
                edges = []
                for value in candidates:
                    budget.charge("transitions")
                    child = build(_specialize(canonical, variable, value, tracked),
                                  variables[:variable] + variables[variable + 1:])
                    if child.count:
                        edges.append((value, child))
                node = _Node("branch", variables, sum(child.count for _, child in edges),
                             (variable, tuple(edges)))
        memo[key] = node
        return node

    def witnesses(node):
        if not node.count:
            return
        if node.kind == "accept":
            yield from product(domain, repeat=len(node.variables))
        elif node.kind == "fiber":
            yield from ((root,) for root in roots_for(node.data))
        elif node.kind == "branch":
            index, edges = node.data
            for value, child in edges:
                for row in witnesses(child):
                    yield row[:index] + (value,) + row[index:]
        elif node.kind == "free":
            unused, child = node.data
            used = tuple(i for i in range(len(node.variables)) if i not in unused)
            for row in witnesses(child):
                for free in product(domain, repeat=len(unused)):
                    assignment = [0] * len(node.variables)
                    for index, value in zip(used, row): assignment[index] = value
                    for index, value in zip(unused, free): assignment[index] = value
                    yield tuple(assignment)

    try:
        domain_poly = _domain_polynomial(domain, tracked)
        normalized = _normalize(equations, n_variables, domain_poly, tracked, max_terms)
        root = build(normalized, tuple(range(n_variables)))
        result["solution_count_exact"] = root.count
    except _Limit as error:
        result.update(status="budget_exceeded", solution_count_exact=None,
                      output_complete=False, count_complete=False, reason=str(error))
        metrics.update(budget.counts)
        return result
    try:
        if max_solutions:
            for assignment in witnesses(root):
                # Check retained roots against original, unnormalized input.
                for equation in equations:
                    total = 0
                    for coefficient, powers in equation:
                        term = coefficient
                        for value, exponent in zip(assignment, powers):
                            if exponent:
                                term = tracked.mul(term, tracked.pow(value, exponent))
                        total = tracked.add(total, term)
                    if total:
                        raise ArithmeticError("a DAG witness failed original-equation verification")
                result["solutions"].append(list(assignment))
                if len(result["solutions"]) == max_solutions:
                    break
    except _Limit as error:
        result["reason"] = f"Exact count obtained; witness extraction stopped: {error}"
    result["output_complete"] = len(result["solutions"]) == result["solution_count_exact"]
    if not result["output_complete"]:
        result["status"] = "output_truncated"
        if result["reason"].startswith("All constrained"):
            result["reason"] = "Exact count obtained; retained witnesses are capped by max_solutions."
    result["solutions"].sort()
    metrics.update(budget.counts)
    return result
