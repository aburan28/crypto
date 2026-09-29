#!/usr/bin/env python3
"""Replay a frozen, at-most-eight-variable Boolean algebra corpus.

No curve arithmetic, target generation, decomposition API, or external input
service is provided. This is a bounded evidence checker, not a solver backend.
"""
import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path


class InvalidEvidence(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidEvidence(message)


def integer(value, low, high, name):
    require(type(value) is int and low <= value <= high, f"invalid {name}")
    return value


class Ring:
    def __init__(self, variables):
        self.n = integer(variables, 1, 8, "variables")
        self.size = 1 << self.n
        self.masks = sorted(range(self.size), key=lambda x: (x.bit_count(), x))
        self.index = {mask: i for i, mask in enumerate(self.masks)}

    def decode(self, value):
        require(isinstance(value, str) and value.startswith("0x")
                and 3 <= len(value) <= 2 + (self.size + 3) // 4,
                "invalid polynomial encoding")
        try:
            result = int(value, 16)
        except ValueError as exc:
            raise InvalidEvidence("invalid hexadecimal polynomial") from exc
        require(0 <= result < (1 << self.size), "polynomial exceeds ring")
        return result

    def degree(self, poly):
        return self.masks[poly.bit_length() - 1].bit_count() if poly else -1

    def terms(self, poly):
        while poly:
            bit = poly & -poly
            yield self.masks[bit.bit_length() - 1]
            poly ^= bit

    def evaluate(self, poly, assignment):
        return sum((mask & assignment) == mask for mask in self.terms(poly)) % 2

    def multiply_variable(self, poly, variable):
        result = 0
        for mask in self.terms(poly):
            result ^= 1 << self.index[mask | (1 << variable)]
        return result


def remainder(poly, pivots):
    while poly:
        pivot = poly.bit_length() - 1
        if pivot not in pivots:
            break
        poly ^= pivots[pivot]
    return poly


def closure(ring, equations, bound):
    """Close under variable multiples of rows of degree strictly below D."""
    pivots, pending = {}, deque()

    def insert(poly):
        require(ring.degree(poly) <= bound, "row exceeds degree bound")
        poly = remainder(poly, pivots)
        if poly:
            pivots[poly.bit_length() - 1] = poly
            if ring.degree(poly) < bound:
                pending.append(poly)

    for poly in equations:
        insert(poly)
    while pending:
        poly = pending.popleft()
        for variable in range(ring.n):
            insert(ring.multiply_variable(poly, variable))
    return pivots


def check_reduced_basis(ring, basis, solutions):
    """Check reduced form, zeros, and quotient dimension by enumeration.

    In the finite Boolean ring, ideals correspond to their zero sets. Vanishing
    and a standard-monomial count equal to the zero-set size prove completeness.
    """
    require(all(basis) and len(set(basis)) == len(basis), "zero or duplicate basis row")
    leads = [ring.masks[p.bit_length() - 1] for p in basis]
    for i, poly in enumerate(basis):
        for mask in ring.terms(poly):
            for j, lead in enumerate(leads):
                if i == j and mask == leads[i]:
                    continue
                require(mask & lead != lead, "basis is not reduced")
        require(all(ring.evaluate(poly, a) == 0 for a in solutions),
                "basis fails at an input solution")
    standard = [m for m in ring.masks if not any(m & lead == lead for lead in leads)]
    require(len(standard) == len(solutions), "basis quotient dimension mismatch")


def check_sympy(ring, equations, basis):
    import sympy as sp
    xs = sp.symbols(f"x0:{ring.n}")
    order = tuple(reversed(xs))

    def expression(poly):
        return sum(sp.prod(xs[i] for i in range(ring.n) if mask >> i & 1)
                   for mask in ring.terms(poly))

    actual = sp.groebner([expression(p) for p in equations]
                         + [x*x + x for x in xs], *order, order="grlex", modulus=2)
    leads = [ring.masks[p.bit_length() - 1] for p in basis]
    expected = [expression(p) for p in basis]
    expected += [xs[i]**2 + xs[i] for i in range(ring.n)
                 if not any((1 << i) & lead == lead for lead in leads)]

    def normalized(poly):
        return tuple(sorted(sp.Poly(poly, *order, modulus=2).terms()))

    require({normalized(p.as_expr()) for p in actual.polys}
            == {normalized(p) for p in expected}, "SymPy basis mismatch")


def check_record(record, use_sympy=False):
    require(isinstance(record, dict), "record is not an object")
    require(isinstance(record.get("id"), str) and record["id"], "missing id")
    ring = Ring(record.get("variables"))
    equations_hex = record.get("equations_hex")
    basis_hex = record.get("boolean_gb_hex")
    require(isinstance(equations_hex, list) and 1 <= len(equations_hex) <= 32,
            "invalid equation count")
    require(isinstance(basis_hex, list) and len(basis_hex) <= ring.size,
            "invalid basis count")
    equations = [ring.decode(p) for p in equations_hex]
    basis = [ring.decode(p) for p in basis_hex]
    integer(record.get("input_degree"), 0, ring.n, "input degree")
    input_degree = max(0, max(map(ring.degree, equations)))
    require(record.get("input_degree") == input_degree, "input degree mismatch")
    solutions = [a for a in range(ring.size)
                 if all(ring.evaluate(p, a) == 0 for p in equations)]
    claimed_solutions = record.get("solution_assignments")
    require(isinstance(claimed_solutions, list)
            and all(type(a) is int and 0 <= a < ring.size for a in claimed_solutions),
            "invalid solution assignments")
    recorded_solutions = record.get("solution_assignments")
    require(isinstance(recorded_solutions, list)
            and all(type(a) is int for a in recorded_solutions)
            and recorded_solutions == solutions, "solution set mismatch")
    check_reduced_basis(ring, basis, solutions)
    bound = integer(record.get("truncated_closure_bound"), max(1, input_degree),
                    ring.n + 1, "truncated closure bound")
    profiles = record.get("profiles_through_bound")
    require(isinstance(profiles, list) and len(profiles) == bound - max(1, input_degree) + 1,
            "missing degree profiles")
    for degree, expected in zip(range(max(1, input_degree), bound + 1), profiles):
        require(isinstance(expected, dict)
                and type(expected.get("degree")) is int
                and type(expected.get("rank")) is int
                and type(expected.get("contains_basis")) is bool, "invalid degree profile")
        pivots = closure(ring, equations, degree)
        complete = all(remainder(p, pivots) == 0 for p in basis)
        require(expected == {"degree": degree, "rank": len(pivots),
                             "contains_basis": complete}, "degree profile mismatch")
        require(complete == (degree == bound), "claimed bound is not first completion")
    if use_sympy:
        check_sympy(ring, equations, basis)
    return {"solutions": len(solutions), "bound": bound}


def load_corpus(root):
    raw = (root / "corpus.jsonl").read_bytes()
    require(len(raw) <= 4_000_000, "corpus exceeds checker budget")
    manifest = json.loads((root / "manifest.json").read_text())
    require(manifest.get("schema_version") == 1, "unsupported schema")
    require(hashlib.sha256(raw).hexdigest() == manifest.get("corpus_sha256"),
            "corpus checksum mismatch")
    rows = [json.loads(line) for line in raw.splitlines()]
    require(len(rows) == manifest.get("records") == 1152, "corpus count mismatch")
    require(all(isinstance(r, dict) and isinstance(r.get("id"), str) for r in rows),
            "invalid corpus record")
    require(len({r["id"] for r in rows}) == len(rows), "duplicate record id")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sympy", action="store_true", help="also compare with SymPy")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    rows = load_corpus(root)
    counts, bounds = Counter(), Counter()
    for row in rows:
        try:
            result = check_record(row, args.sympy)
        except (InvalidEvidence, KeyError, TypeError) as exc:
            raise InvalidEvidence(f"{row['id']}: {exc}") from exc
        counts["sat" if result["solutions"] else "unsat"] += 1
        bounds[result["bound"]] += 1
    print(json.dumps({"records_checked": len(rows), "verdicts": dict(counts),
                      "truncated_closure_bounds": dict(sorted(bounds.items())),
                      "sympy_checked": args.sympy, "regularity_claim": False}, sort_keys=True))


if __name__ == "__main__":
    main()
