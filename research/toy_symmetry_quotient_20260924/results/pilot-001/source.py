#!/usr/bin/env python3
"""Fixed GF(128) symmetry correctness pilot; no configurable cryptographic target.

Enumerates complete point-decomposition catalogues, not Groebner bases or DLPs.
Only the standard library is needed. See README.md for accounting definitions.
"""

import argparse
from collections import Counter, defaultdict
from itertools import combinations_with_replacement, permutations, product
import hashlib
import json
from math import comb
from pathlib import Path
import platform
import random
import statistics
import subprocess
import time

O = (-1, -1)
MODULUS = 0x83  # X^7 + X + 1, in the polynomial basis.
SEEDS = (101, 202, 503, 607)
VARIANTS = ("ordered", "permutation", "invariant_lift")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


class ToyCurve:
    """K_0: y^2+xy=x^3+1 over the fixed seven-bit field."""

    def __init__(self):
        self.counts = Counter()

    def mul(self, a, b):
        self.counts["field_multiply"] += 1
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & 128:
                a ^= MODULUS
        return result

    def square(self, a):
        self.counts["field_square"] += 1
        return self.mul(a, a)

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("zero has no field inverse")
        self.counts["field_invert"] += 1
        result, exponent = 1, 126
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.square(a)
            exponent >>= 1
        return result

    def contains(self, p):
        if p == O:
            return True
        x, y = p
        return (0 <= x < 128 and 0 <= y < 128
                and self.square(y) ^ self.mul(x, y) == self.mul(self.square(x), x) ^ 1)

    def points(self):
        return (O,) + tuple((x, y) for x in range(128) for y in range(128)
                            if self.contains((x, y)))

    @staticmethod
    def neg(p):
        return O if p == O else (p[0], p[0] ^ p[1])

    def add(self, p, q):
        self.counts["ec_add"] += 1  # Includes identity calls and doublings.
        if p == O:
            return q
        if q == O:
            return p
        x, y = p
        u, v = q
        if x == u:
            if y != v or x == 0:
                return O
            self.counts["ec_double"] += 1
            slope = x ^ self.mul(y, self.inv(x))
            out_x = self.square(slope) ^ slope
            out_y = self.square(x) ^ self.mul(slope ^ 1, out_x)
        else:
            slope = self.mul(y ^ v, self.inv(x ^ u))
            out_x = self.square(slope) ^ slope ^ x ^ u
            out_y = self.mul(slope, x ^ out_x) ^ out_x ^ y
        return out_x, out_y

    def total(self, points):
        result = O
        for p in points:
            result = self.add(result, p)
        return result

    def frobenius(self, p, power=1):
        for _ in range(power % 7):
            self.counts["ec_frobenius"] += 1
            if p != O:
                p = self.square(p[0]), self.square(p[1])
        return p


def point_orbit(curve, p):
    return tuple(sorted({curve.frobenius(p, k) for k in range(7)}))


def require_closed(curve, support):
    missing = {curve.frobenius(p) for p in support} - set(support)
    if missing:
        raise ValueError("factor base is not Frobenius invariant")


def fixtures(curve, points):
    """Frozen seed cells, independent holdouts, and a linear-subspace control."""
    signed_orbits = sorted({tuple(sorted(set(point_orbit(curve, p)) |
                                        set(point_orbit(curve, curve.neg(p)))))
                            for p in points if len(point_orbit(curve, p)) == 7})
    cases = []
    for seed in SEEDS:
        rng = random.Random(seed)
        support = tuple(sorted(((0, 1),) + rng.choice(signed_orbits)))
        cases.append({"id": f"orbit-s{seed}", "seed": seed,
                      "split": "frozen" if seed in SEEDS[:2] else "holdout",
                      "kind": "signed_orbit", "support": support})
        ordinary = tuple(sorted(rng.sample(list(points[1:]), len(support))))
        cases.append({"id": f"ordinary-s{seed}", "seed": seed,
                      "split": "frozen" if seed in SEEDS[:2] else "holdout",
                      "kind": "matched_random_points", "support": ordinary})
    # Both degree-three factors give invariant subspaces, but different point yields.
    for middle in (2, 4):
        subspace = []
        for x in range(128):
            x2 = curve.square(x)
            x4 = curve.square(x2)
            if curve.square(x4) ^ (x2 if middle == 2 else x4) ^ x == 0:
                subspace.append(x)
        if len(subspace) != 8:
            raise AssertionError("wrong invariant subspace dimension")
        support = tuple(p for p in points[1:] if p[0] in subspace)
        cases.append({"id": f"linear-d3-x{middle}", "seed": None, "split": "control",
                      "kind": f"linearized_x8_x{middle}_x", "support": support,
                      "subspace": subspace})
    return cases


def root_polynomial(curve, roots):
    """Low-to-high coefficients of product(T+x_i), i.e. symmetric invariants."""
    coefficients = [1]
    for root in roots:
        new = [0] * (len(coefficients) + 1)
        for i, c in enumerate(coefficients):
            new[i] ^= curve.mul(c, root)
            new[i + 1] ^= c
        coefficients = new
    return tuple(coefficients)


def lift_roots(curve, coefficients, allowed_x):
    """Recover a multiset, including multiplicities; reject non-admissible roots."""
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("expected a monic cubic invariant polynomial")
    work, roots = list(coefficients), []
    for x in sorted(set(allowed_x)):
        while len(work) > 1:
            curve.counts["root_tests"] += 1
            quotient = [0] * (len(work) - 1)
            quotient[-1] = work[-1]
            for i in range(len(quotient) - 2, -1, -1):
                quotient[i] = work[i + 1] ^ curve.mul(x, quotient[i + 1])
            remainder = work[0] ^ curve.mul(x, quotient[0])
            if remainder:
                break
            roots.append(x)
            work = quotient
    if work != [1] or len(roots) != 3:
        raise ValueError("invariant polynomial does not split in the allowed support")
    return tuple(roots)


def candidates(curve, support, variant):
    if variant == "ordered":
        yield from product(support, repeat=3)
    elif variant == "permutation":
        yield from combinations_with_replacement(support, 3)
    elif variant == "invariant_lift":
        by_x = defaultdict(list)
        for p in support:
            by_x[p[0]].append(p)
        for xs in combinations_with_replacement(sorted(by_x), 3):
            curve.counts["invariant_polynomials"] += 1
            coefficients = root_polynomial(curve, xs)
            recovered = lift_roots(curve, coefficients, by_x)
            # Opposite points share x; retain all admissible y lifts.
            yield from product(*(by_x[x] for x in recovered))
    else:
        raise ValueError("unknown variant")


def catalogue(curve, support, variant):
    """Complete enumeration, grouped by target; never a first-solution timing."""
    if not support or len(support) > 24 or len(set(support)) != len(support):
        raise ValueError("support must contain 1..24 distinct points")
    if O in support or any(not curve.contains(p) for p in support):
        raise ValueError("support must contain affine points of the fixed toy curve")
    result = defaultdict(set)
    for triple in candidates(curve, support, variant):
        curve.counts["tuple_sum_checks"] += 1
        result[curve.total(triple)].add(tuple(sorted(triple)))
    return dict(result)


def encode_catalogue(catalog, ids):
    return [[ids[target], sorted([list(map(ids.__getitem__, triple)) for triple in triples])]
            for target, triples in sorted(catalog.items())]


def verify_catalogue(curve, catalog, support):
    allowed = set(support)
    for target, triples in catalog.items():
        for triple in triples:
            if (tuple(sorted(triple)) != triple or len(triple) != 3
                    or any(p not in allowed for p in triple) or curve.total(triple) != target):
                raise AssertionError("invalid point decomposition")


def frobenius_audit(curve, catalog, support, points, ids):
    """Whole-target transport, inverse lifting, and fixed-target stabilizers."""
    require_closed(curve, support)
    orbit_reps, solution_orbits, independent_counterexample = set(), set(), None
    fixed_by_power = [0] * 7
    fixed_target_checks = 0
    for target in points:
        orbit_reps.add(min(point_orbit(curve, target)))
        original = catalog.get(target, set())
        for k in range(7):
            mapped_target = curve.frobenius(target, k)
            transformed = {tuple(sorted(curve.frobenius(p, k) for p in t)) for t in original}
            if transformed != catalog.get(mapped_target, set()):
                raise AssertionError("Frobenius changed exact solution coverage")
            restored = {tuple(sorted(curve.frobenius(p, -k) for p in t)) for t in transformed}
            if restored != original:
                raise AssertionError("inverse transport failed")
            if mapped_target == target:
                fixed_target_checks += 1
                if transformed != original:
                    raise AssertionError("target stabilizer mismatch")
        for triple in original:
            orbit = []
            for k in range(7):
                mapped = tuple(sorted(curve.frobenius(p, k) for p in triple))
                mapped_target = curve.frobenius(target, k)
                orbit.append((mapped_target, mapped))
                fixed_by_power[k] += int((mapped_target, mapped) == (target, triple))
            solution_orbits.add(min(orbit))
            bad = (curve.frobenius(triple[0]), triple[1], triple[2])
            bad_target = curve.total(bad)
            if bad_target != target and independent_counterexample is None:
                independent_counterexample = {"target": ids[target], "triple": [ids[p] for p in triple],
                    "rotated_triple": [ids[p] for p in bad], "changed_target": ids[bad_target]}
    if sum(fixed_by_power) != 7 * len(solution_orbits):
        raise AssertionError("Burnside orbit count mismatch")
    if independent_counterexample is None and any(curve.frobenius(p) != p for p in support):
        raise AssertionError("negative control did not expose independent rotation")
    return {"status": "VERIFIED", "target_orbits": len(orbit_reps),
            "target_count": len(points), "solution_orbits": len(solution_orbits),
            "fixed_solution_pairs_by_power": fixed_by_power,
            "fixed_target_checks": fixed_target_checks,
            "negative_control_status": "COUNTEREXAMPLE_FOUND" if independent_counterexample else "TRIVIAL_ACTION",
            "independent_rotation_counterexample": independent_counterexample,
            "solver_calls_saved": None, "independent_relation_rank": None}


def permutation_audit(support, catalog):
    """Burnside and orbit-stabilizer counts, including repeated summands."""
    slots = list(permutations(range(3)))
    fixed = [0] * len(slots)
    for t in product(support, repeat=3):
        for i, perm in enumerate(slots):
            fixed[i] += int(tuple(t[j] for j in perm) == t)
    unordered = sum(map(len, catalog.values()))
    reconstructed_ordered = sum(len(set(permutations(t)))
                                for triples in catalog.values() for t in triples)
    if sum(fixed) != 6 * unordered or reconstructed_ordered != len(support) ** 3:
        raise AssertionError("permutation orbit coverage failed")
    return {"fixed_tuples_by_permutation": fixed, "unordered": unordered,
            "reconstructed_ordered": reconstructed_ordered}


def run(repetitions=3):
    if not 1 <= repetitions <= 5:
        raise ValueError("repetitions must be between one and five")
    setup_curve = ToyCurve()
    start = time.perf_counter()
    points = setup_curve.points()
    cases = fixtures(setup_curve, points)
    setup_seconds = time.perf_counter() - start
    if len(points) != 116:
        raise AssertionError("trace-recurrence group order mismatch")
    ids = {p: i for i, p in enumerate(points)}
    instance_data = {"field_degree": 7, "field_modulus": MODULUS,
        "field_basis": "polynomial", "curve": "y^2+xy=x^3+1", "m": 3,
        "point_encoding": [list(p) for p in points], "identity_encoding": list(O),
        "target_sampling": "all 116 rational points, including identity and unsatisfiable targets",
        "cases": [{**c, "support": [ids[p] for p in c["support"]]} for c in cases]}
    certificates, measurements, summaries = {}, [], []
    for case in cases:
        support = case["support"]
        # Independent ordered ground truth is never passed into a candidate.
        truth = catalogue(ToyCurve(), support, "ordered")
        encoded = encode_catalogue(truth, ids)
        perm_audit = permutation_audit(support, truth)
        checker = ToyCurve()
        try:
            require_closed(checker, support)
        except ValueError:
            frob = {"status": "REJECTED_NON_INVARIANT_SUPPORT"}
        else:
            frob = frobenius_audit(checker, truth, support, points, ids)
        certificates[case["id"]] = {"catalogue": encoded, "sha256": digest(encoded),
                                    "permutation": perm_audit, "frobenius": frob}
        for rep in range(repetitions):
            order = list(VARIANTS)
            random.Random(f"{case['id']}-{rep}").shuffle(order)
            for variant in order:
                curve = ToyCurve()
                start = time.perf_counter()
                found = catalogue(curve, support, variant)
                enumeration_seconds = time.perf_counter() - start
                enumeration_counts = dict(curve.counts)
                curve.counts.clear()
                start = time.perf_counter()
                verify_catalogue(curve, found, support)
                verification_seconds = time.perf_counter() - start
                actual = encode_catalogue(found, ids)
                if actual != encoded:
                    raise AssertionError(f"coverage mismatch: {case['id']}/{variant}")
                measurements.append({"case": case["id"], "variant": variant, "repetition": rep,
                    "status": "VERIFIED_COMPLETE_ENUMERATION", "enumeration_counts": enumeration_counts,
                    "verification_counts": dict(curve.counts),
                    "enumeration_seconds": enumeration_seconds,
                    "verification_seconds": verification_seconds,
                    "certificate_sha256": digest(actual)})
        for variant in VARIANTS:
            rows = [r for r in measurements if r["case"] == case["id"] and r["variant"] == variant]
            checks = rows[0]["enumeration_counts"]["tuple_sum_checks"]
            if any(r["enumeration_counts"] != rows[0]["enumeration_counts"] for r in rows):
                raise AssertionError("non-reproducible operation counts")
            bound = comb(len(support) + 2, 3)
            summaries.append({"case": case["id"], "variant": variant, "class": "accounting",
                "support_size": len(support), "tuple_sum_checks": checks,
                "enumeration_count_floor": bound, "checks_over_enumeration_floor": checks / bound,
                "checks_over_ordered_reference": checks / len(support) ** 3,
                "verified_multisets": sum(map(len, truth.values())),
                "solvable_targets": len(truth), "unsatisfiable_targets": len(points) - len(truth),
                "median_enumeration_and_verification_seconds": statistics.median(
                    r["enumeration_seconds"] + r["verification_seconds"] for r in rows),
                "correct": True, "max_processed_degree": None, "f4_matrix_rows": None,
                "f4_matrix_columns": None, "full_dlp_S": None, "ratio_to_rho": None,
                "full_dlp_speedup": None})
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True,
                                         cwd=Path(__file__).parent).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    report = {"schema_version": 1, "scope": "fixed_toy_complete_enumeration_stage_diagnostic",
        "evidence_type": "measured_local", "classification": "accounting",
        "source_base_commit": commit, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "instance_sha256": digest(instance_data), "certificate_sha256": digest(certificates),
        "python": platform.python_version(), "platform": platform.platform(),
        "repetitions": repetitions, "shared_fixture_setup_seconds": setup_seconds,
        "shared_fixture_setup_counts": dict(setup_curve.counts),
        "shared_setup_accounting": "once per whole suite; not included in per-row medians",
        "ground_truth_and_audit_cost": "test infrastructure; excluded from candidate enumeration timings",
        "groebner_status": "NOT_RUN", "neighbor_status": "NOT_RUN", "full_dlp_status": "NOT_RUN",
        "primary_unit": "tuple sum checks during complete enumeration, not common operations",
        "floor_scope": "each explicit enumeration visits every point multiset; not a universal solver lower bound",
        "counter_note": "field_multiply includes calls inside square/invert; ec_add includes doubles/identity; do not sum nested counters",
        "summary": summaries, "measurements": measurements}
    return instance_data, certificates, report


def write_new_results(directory, repetitions):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("use a new result directory; evidence is never overwritten")
    instances, certificates, report = run(repetitions)
    directory.mkdir(parents=True)
    for name, value in (("instances.json", instances), ("certificates.json", certificates), ("results.json", report)):
        with (directory / name).open("x") as handle:
            json.dump(value, handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
            handle.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3, choices=range(1, 6))
    args = parser.parse_args()
    report = write_new_results(args.output, args.repetitions)
    print(json.dumps({"output": str(args.output), "variants": len(report["summary"]),
                      "runs": len(report["measurements"]), "all_correct": True,
                      "groebner_status": report["groebner_status"]}))


if __name__ == "__main__":
    main()
