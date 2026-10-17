"""Independent exact checker for odd-field pole-four coefficient reduction.

Run from any directory:
    python lab/research/odd_algebra_check.py --output \
        lab/research/odd_algebra_check_results.json

This checks algebra and exceptional cases, not performance or novelty.
All polynomial coefficients use ascending order and canonical field encodings.
"""

import argparse
import hashlib
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from rootfinder.fields import ExtensionField, PrimeField
from rootfinder import polynomial as P


class Arithmetic:
    def __init__(self, field):
        self.f = field

    def scalar(self, integer):
        return integer % self.f.characteristic

    def add(self, *values):
        result = 0
        for value in values:
            result = self.f.add(result, value)
        return result

    def mul(self, *values):
        result = 1
        for value in values:
            result = self.f.mul(result, value)
        return result

    def neg(self, value):
        return self.f.neg(value)

    def sub(self, left, right):
        return self.f.sub(left, right)

    def square(self, value):
        return self.mul(value, value)

    def cube(self, value):
        return self.mul(value, value, value)

    def evaluate(self, polynomial, value):
        result = 0
        for coefficient in reversed(polynomial):
            result = self.add(self.mul(result, value), coefficient)
        return result


def formulas(field, a2, a4, z, y, s, A, B):
    """Return E(b), C(b), R(b), and verify the polynomial identity."""
    a = Arithmetic(field)
    add, mul, neg, sub, sq, cb, sc = (
        a.add, a.mul, a.neg, a.sub, a.square, a.cube, a.scalar
    )
    quarter = field.inv(sc(4))
    half = field.inv(sc(2))
    e2 = (
        mul(quarter, add(sq(s), mul(sc(2), s, z), neg(mul(sc(3), sq(z))))),
        neg(mul(sc(2), y)),
        neg(add(a2, mul(half, s), mul(sc(3), half, z))),
        0,
        quarter,
    )
    e3 = (
        mul(quarter, add(mul(sq(s), z), neg(mul(sc(2), s, sq(z))), cb(z))),
        mul(sub(z, s), y),
        add(mul(a2, z), a4, neg(mul(half, s, z)), mul(sc(3), half, sq(z))),
        y,
        mul(z, quarter),
    )
    E = P.sub(e2, (A,), field)
    C_literal = P.sub(P.sub(e3, tuple(mul(z, v) for v in e2), field),
                      (sub(B, mul(z, A)),), field)
    k = add(mul(sc(3), sq(z)), mul(sc(2), a2, z), a4)
    m = add(cb(z), neg(mul(s, sq(z))), neg(B), mul(z, A))
    C = P.trim((m, mul(y, sub(mul(sc(3), z), s)), k, y), field)
    assert C_literal == C
    R = P.trim((
        add(mul(sq(y), add(sq(s), mul(sc(2), s, z), neg(mul(sc(3), sq(z))),
                          neg(mul(sc(4), A)))), mul(k, m)),
        mul(y, add(neg(mul(sc(8), sq(y))), neg(m),
                   mul(k, sub(mul(sc(3), z), s)))),
        sub(sq(k), mul(add(mul(sc(4), a2), s, mul(sc(9), z)), sq(y))),
    ), field)
    literal_R = P.sub(tuple(mul(sc(4), sq(y), v) for v in E),
                      P.mul((neg(k), y), C, field), field)
    assert literal_R == R
    return e2, e3, E, C, R


def literal_norm(field, a2, a4, a6, z, y, b, s):
    """Compute (X²+cX+d)²-b²C(X), then divide exactly by X-z."""
    a = Arithmetic(field)
    c = a.mul(a.sub(a.sub(a.square(b), z), s), field.inv(a.scalar(2)))
    d = a.neg(a.add(a.square(z), a.mul(c, z), a.mul(b, y)))
    F = (d, c, 1)
    curve = (a6, a4, a2, 1)
    norm = P.sub(P.mul(F, F, field), tuple(a.mul(a.square(b), v) for v in curve), field)
    H, remainder = P.divmod_poly(norm, (a.neg(z), 1), field)
    assert not remainder
    return H, c, d


def run_random_checks():
    rng = random.Random(20261012)
    records = []
    for field in (PrimeField(3), PrimeField(5), ExtensionField(3, 2)):
        a = Arithmetic(field)
        counters = dict(random_draws=300, point_parameter_cases=0,
                        norm_divisions_checked=0, reduced_equivalences_checked=0,
                        nonzero_y_zero_remainders=0, zero_y_cases=0)
        for _ in range(counters["random_draws"]):
            a2, a4, a6, z, s, A, B = (rng.randrange(field.order) for _ in range(7))
            cy = a.add(a.cube(z), a.mul(a2, a.square(z)), a.mul(a4, z), a6)
            ys = [value for value in range(field.order) if a.square(value) == cy]
            if not ys:
                continue
            y = next((value for value in ys if value), 0)
            e2, e3, E, C, R = formulas(field, a2, a4, z, y, s, A, B)
            counters["point_parameter_cases"] += 1
            counters["zero_y_cases"] += int(y == 0)
            counters["nonzero_y_zero_remainders"] += int(bool(y) and not R)
            for b in range(field.order):
                H, _, _ = literal_norm(field, a2, a4, a6, z, y, b, s)
                expected_H = (a.neg(a.evaluate(e3, b)), a.evaluate(e2, b), a.neg(s), 1)
                assert H == P.trim(expected_H, field)
                counters["norm_divisions_checked"] += 1
                if y:
                    original = a.evaluate(E, b) == 0 and a.evaluate(C, b) == 0
                    reduced = a.evaluate(R, b) == 0 and a.evaluate(C, b) == 0
                    assert original == reduced
                    counters["reduced_equivalences_checked"] += 1
                else:
                    assert len(C) <= 3
        records.append({"order": field.order, "characteristic": field.characteristic,
                        "degree": field.degree, "modulus": getattr(field, "modulus", None),
                        **counters})
    return records


def exceptional_cases():
    field = PrimeField(3)
    # Nonsingular y²=x³+x+1; derivative is identically 1 in characteristic 3.
    e2, e3, E, C, R = formulas(field, 0, 1, 0, 1, 1, 1, 0)
    assert not R and len(C) == 4 and C == (0, 2, 1, 1)
    assert P.gcd_poly((1, 1, 0, 1), (1,), field) == (1,)
    a = Arithmetic(field)
    roots = [b for b in range(3) if a.evaluate(E, b) == a.evaluate(C, b) == 0]
    cubic_roots = [b for b in range(3) if a.evaluate(C, b) == 0]
    assert roots == cubic_roots == [0]
    fallback = {"field_order": 3, "curve": {"a2": 0, "a4": 1, "a6": 1},
                "target": [0, 1], "s": 1, "A": 1, "B": 0,
                "E": E, "C": C, "R": R, "candidate_roots": cubic_roots,
                "conclusion": "R vanishes identically with nonzero y; cubic fallback is necessary."}

    field = PrimeField(5)
    # A nonsingular curve with three distinct b candidates in the fallback.
    e2, e3, E, C, R = formulas(field, 0, 4, 2, 1, 3, 1, 3)
    assert not R and C == (0, 3, 1, 1)
    assert P.gcd_poly((0, 4, 0, 1), (4, 0, 3), field) == (1,)
    a = Arithmetic(field)
    roots = [b for b in range(5) if a.evaluate(E, b) == a.evaluate(C, b) == 0]
    assert roots == [0, 1, 3]
    for b in roots:
        H, _, _ = literal_norm(field, 0, 4, 0, 2, 1, b, 3)
        assert H == (2, 1, 2, 1)
    fallback["three_candidate_example"] = {
        "field_order": 5, "curve": {"a2": 0, "a4": 4, "a6": 0},
        "target": [2, 1], "s": 3, "A": 1, "B": 3,
        "E": E, "C": C, "R": R, "candidate_roots": roots,
        "conclusion": "Three distinct b candidates disprove a uniform two-candidate bound."
    }
    # H=(X-2)(X²+2), where X²+2 has no GF(5) root.
    H, c, d = literal_norm(field, 0, 0, 1, 0, 1, 1, 2)
    assert H == (1, 2, 3, 1) and c == 2 and d == 4
    support_polynomial = (3, 1)
    irreducible_part = (2, 0, 1)
    assert P.mul(support_polynomial, irreducible_part, field) == H
    frobenius = (0, 4, 0, 0, 0, 1)
    assert P.gcd_poly(H, frobenius, field) == support_polynomial
    assert P.gcd_poly(irreducible_part, frobenius, field) == (1,)
    assert P.gcd_poly((1, 0, 0, 1), (0, 0, 3), field) == (1,)
    counterexample = {
        "field_order": 5, "curve": {"a2": 0, "a4": 0, "a6": 1},
        "target": [0, 1], "b": 1, "s": 2, "c": c, "d": d,
        "H": H, "support": [2], "support_polynomial": support_polynomial,
        "gcd_H_with_Xq_minus_X": support_polynomial,
        "irreducible_remaining_factor": irreducible_part,
        "conclusion": "One support root does not certify three support roots or field splitting."
    }
    return fallback, counterexample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = run_random_checks()
    fallback, counterexample = exceptional_cases()
    result = {"status": "PASS", "seed": 20261012, "checker_sha256":
              hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Independent literal norm division, coefficient formulas, polynomial reduction, and degeneracies; no performance or novelty claim.",
              "fields": records, "zero_remainder_fallback": fallback,
              "support_membership_counterexample": counterexample}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
