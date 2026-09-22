#!/usr/bin/env python3
"""Exhaustive relation geometry on four fixed fields of at most 12 bits.

This educational experiment has no DLP solver, imported targets, actual Oakley
parameters, network use or scalable field configuration.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import time

HERE = Path(__file__).resolve().parent
MODULI = {6: 0x43, 9: 0x211, 10: 0x409, 12: 0x1053}


def poly_mod(a: int, b: int) -> int:
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_mod(a, b)
    return a


def basis(values) -> list[int]:
    pivots = {}
    for value in values:
        while value:
            p = value.bit_length() - 1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
    return [pivots[p] for p in sorted(pivots)]


def span(vectors: list[int]) -> list[int]:
    result = [0]
    for v in vectors:
        result += [x ^ v for x in result]
    return sorted(result)


class ToyField:
    def __init__(self, n: int):
        if n not in MODULI:
            raise ValueError("Only the four fixed toy fields (6, 9, 10, 12 bits) are supported")
        self.n, self.modulus, self.size = n, MODULI[n], 1 << n
        x = 2
        for i in range(1, n + 1):
            x = self.mul(x, x)
            if i <= n // 2 and poly_gcd(x ^ 2, self.modulus) != 1:
                raise ValueError("Reducible toy modulus")
        if x != 2:
            raise ValueError("Invalid toy field")

    def mul(self, a: int, b: int) -> int:
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & self.size:
                a ^= self.modulus
        return out

    def frobenius(self, a: int, k: int) -> int:
        for _ in range(k):
            a = self.mul(a, a)
        return a

    def inv(self, a: int) -> int:
        if not a:
            raise ZeroDivisionError("zero has no inverse")
        u, v, g, h = a, self.modulus, 1, 0
        while u != 1:
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v, g, h = v, u, h, g
                shift = -shift
            u ^= v << shift
            g ^= h << shift
        return poly_mod(g, self.modulus)

    def degree(self, a: int) -> int:
        x = self.frobenius(a, 1)
        d = 1
        while x != a:
            x = self.frobenius(x, 1)
            d += 1
        return d


class ToyCurve:
    def __init__(self, field: ToyField, b: int):
        if not 0 < b < field.size:
            raise ValueError("nonzero toy coefficient required")
        self.f, self.b = field, b
        # z^2 + z has either zero or two roots, differing by 1.
        roots = {}
        for z in range(field.size):
            roots.setdefault(field.mul(z, z) ^ z, z)
        self.by_x = {0: [(0, field.frobenius(b, field.n - 1))]}
        for x in range(1, field.size):
            ix = field.inv(x)
            c = x ^ field.mul(b, field.mul(ix, ix))
            if c in roots:
                y = field.mul(x, roots[c])
                self.by_x[x] = [(x, y), (x, y ^ x)]
        self.affine = {p for points in self.by_x.values() for p in points}
        self.order = len(self.affine) + 1
        if not all(self.on_curve(p) for p in self.affine):
            raise AssertionError("invalid point enumeration")
        if (self.order - field.size - 1) ** 2 > 4 * field.size:
            raise AssertionError("point count violates Hasse")

    def on_curve(self, p) -> bool:
        if p is None:
            return True
        f, (x, y) = self.f, p
        return f.mul(y, y) ^ f.mul(x, y) == f.mul(f.mul(x, x), x) ^ self.b

    def neg(self, p):
        return None if p is None else (p[0], p[0] ^ p[1])

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        f, (x, y), (u, v) = self.f, p, q
        if x == u:
            if y != v or x == 0:
                return None
            t = x ^ f.mul(y, f.inv(x))
            xx = f.mul(t, t) ^ t
            return xx, f.mul(x, x) ^ f.mul(t ^ 1, xx)
        t = f.mul(y ^ v, f.inv(x ^ u))
        xx = f.mul(t, t) ^ t ^ x ^ u
        return xx, f.mul(t, x ^ xx) ^ xx ^ y


def sha_json(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def random_basis(f: ToyField, dimension: int, rng: random.Random) -> list[int]:
    out = []
    while len(out) < dimension:
        out = basis(out + [rng.randrange(1, f.size)])
    return out


def inspect_base(curve: ToyCurve, vectors: list[int], k: int) -> dict:
    started = time.perf_counter()
    f = curve.f
    domain = span(vectors)
    fb = sorted(p for x in domain for p in curve.by_x.get(x, []))
    products = basis(f.mul(x, y) for i, x in enumerate(vectors) for y in vectors[i:])
    setup_seconds = time.perf_counter() - started
    started = time.perf_counter()
    sums = Counter()
    infinity_pairs = 0
    checked = 0
    for i, p in enumerate(fb):
        for q in fb[i + 1:]:
            r = curve.add(p, q)
            checked += 1
            if r is None:
                infinity_pairs += 1
                continue
            # Independent equation check, membership and inverse identity.
            if r not in curve.affine or not curve.on_curve(r):
                raise AssertionError("pair sum fails curve equation")
            if curve.add(r, curve.neg(p)) != q:
                raise AssertionError("pair sum fails inverse identity")
            x, y, z = p[0], q[0], r[0]
            s = f.mul(x, y) ^ f.mul(x, z) ^ f.mul(y, z)
            if f.mul(s, s) ^ f.mul(f.mul(x, y), z) ^ curve.b:
                raise AssertionError("pair sum fails independent S3 identity")
            sums[r] += 1
    distinct = len(sums)
    upper = min(checked, curve.order - 1)
    assert distinct <= upper
    assert len(products) <= min(f.n, len(vectors) * (len(vectors) + 1) // 2)
    all_in_subfield = all(f.frobenius(x, k) == x and f.frobenius(y, k) == y for x, y in fb)
    sums_in_subfield = all(f.frobenius(x, k) == x and f.frobenius(y, k) == y for x, y in sums)
    if all_in_subfield and not sums_in_subfield:
        raise AssertionError("subfield subgroup closure failed")
    return {
        "basis": vectors,
        "x_dimension": len(vectors),
        "x_count": len(domain),
        "factor_base_points": len(fb),
        "product_span_dimension": len(products),
        "unordered_distinct_pairs": checked,
        "infinity_pairs_excluded": infinity_pairs,
        "distinct_nonidentity_targets": distinct,
        "coverage_all_curve_points": distinct / curve.order,
        "distinct_target_upper_bound": upper,
        "target_count_over_upper_bound": distinct / upper if upper else None,
        "all_base_points_in_subfield": all_in_subfield,
        "all_pair_sums_in_subfield": sums_in_subfield,
        "nonidentity_pairs_verified": checked - infinity_pairs,
        "target_multiplicity_sha256": sha_json(sorted((x, y, c) for (x, y), c in sums.items())),
        "geometry_setup_seconds": setup_seconds,
        "enumeration_and_checks_seconds": time.perf_counter() - started,
        "S": None,
        "ratio_to_rho": None,
        "ratio_to_dlp_floor": None,
        "full_dlp_total_operations": None,
        "classification": "structural_diagnostic",
    }


def experiment(contract: dict) -> dict:
    rows = []
    for n, k in contract["field_pairs_n_k"]:
        f = ToyField(n)
        subfield = [x for x in range(f.size) if f.frobenius(x, k) == x]
        assert len(subfield) == 1 << k
        u = basis(subfield)
        gamma = next(x for x in range(2, f.size) if x not in subfield)
        tower = basis(u + [f.mul(gamma, x) for x in u])
        assert len(tower) == 2 * k
        rng = random.Random(contract["seed"] + 100 * n + k)
        holdout_b = rng.randrange(3, f.size)
        while f.degree(holdout_b) != n:
            holdout_b = rng.randrange(3, f.size)
        for label, b in [("subfield_control", 1), ("full_degree_main", 2), ("full_degree_holdout", holdout_b)]:
            curve_start = time.perf_counter()
            curve = ToyCurve(f, b)
            curve_seconds = time.perf_counter() - curve_start
            for multiplier in contract["dimension_multipliers"]:
                dimension = multiplier * k
                structured = u if multiplier == 1 else tower
                candidates = {
                    "generic_linear": random_basis(f, dimension, rng),
                    "subfield_linear": structured,
                    "scaled_subfield_linear": basis(f.mul(gamma, x) for x in structured),
                }
                for family, vectors in candidates.items():
                    row = inspect_base(curve, vectors, k)
                    row.update(n=n, k=k, extension_degree=n // k, curve_case=label,
                               curve_a=0, curve_b=b, coefficient_degree=f.degree(b),
                               modulus=f.modulus, curve_order=curve.order, family=family,
                               shared_curve_enumeration_seconds=curve_seconds)
                    rows.append(row)
    screens = []
    for n, k in contract["field_pairs_n_k"]:
        cell = [r for r in rows if r["n"] == n and r["curve_case"] == "full_degree_holdout" and r["x_dimension"] == k]
        reference = next(r for r in cell if r["family"] == "generic_linear")
        admitted = [r["family"] for r in cell if r["family"] != "generic_linear"
                    and r["product_span_dimension"] < reference["product_span_dimension"]
                    and r["distinct_nonidentity_targets"] >= reference["distinct_nonidentity_targets"]]
        screens.append({"n": n, "k": k, "families_passing": admitted})
    return {
        "scope": contract["scope"], "rows": rows, "holdout_screen": screens,
        "screen_pass_cells": sum(bool(x["families_passing"]) for x in screens),
        "screen_passed": sum(bool(x["families_passing"]) for x in screens) >= 3,
        "limitations": [
            "No DLP or polynomial solver; no full-width Oakley computation.",
            "Product-span rank is a structural proxy, not a solving-degree or runtime prediction.",
            "Coverage uses the whole toy curve, not a selected prime-order subgroup.",
            "Matched x dimensions do not imply matched factor-base point counts.",
            "Duplicate target sums are collapsed; relation-matrix independence is unmeasured.",
            "Wall times include checks and are practicality diagnostics only.",
            "Different towers are separate cells, not a scaling-exponent fit.",
            "All target coverage is exact; one generic subspace per cell limits generality.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite existing evidence")
    contract_bytes = (HERE / "contract.json").read_bytes()
    result = experiment(json.loads(contract_bytes))
    result["contract_sha256"] = hashlib.sha256(contract_bytes).hexdigest()
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, sort_keys=True)
        out.write("\n")
    print(json.dumps({"rows": len(result["rows"]), "screen_pass_cells": result["screen_pass_cells"],
                      "screen_passed": result["screen_passed"]}))


if __name__ == "__main__":
    main()
