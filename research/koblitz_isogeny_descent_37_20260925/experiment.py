#!/usr/bin/env sage -python
"""Explicit degree-73 Koblitz descent and matched Boolean stage diagnostic.

Run with:
    sage -python research/koblitz_isogeny_descent_37_20260925/experiment.py \
      --output /tmp/koblitz-37-run
"""
import argparse
import gzip
import hashlib
import json
import random
import sys
import time
from collections import defaultdict
from itertools import product
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, kronecker_symbol, set_random_seed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "research/toy_f5_neighbors_20260924"))
from matrix import Ring, measure


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def digest(value):
    return digest_bytes(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def integer_representation(value):
    return int(value.integer_representation())


def s3(x, y, z, b):
    e2 = x*y + x*z + y*z
    return e2**2 + x*y*z + b


def s4(x, y, z, w, b):
    # Resultant in u of S3(x,y,u) and S3(z,w,u).
    a, bb, c = (x+y)**2, x*y, (x*y)**2 + b
    d, e, f = (z+w)**2, z*w, (z*w)**2 + b
    return (a*f+c*d)**2 + (a*e+bb*d)*(bb*f+c*e)


def normalize_point(point, r, t):
    if point.is_zero():
        return None
    return point[0] + r, point[1] + t


def normalized_coefficients(curve):
    a1, a2, a3, a4, a6 = curve.a_invariants()
    if a1 == 0:
        raise AssertionError("the transported model has a1=0; binary Semaev normalization unavailable")
    r = a3 / a1
    t = (a4 + r**2) / a1
    a2n = a2 + r
    a6n = r**3 + a2*r**2 + a4*r + a6 + t**2
    if a1 != 1 or a3 + r*a1 != 0 or a4 + a1*t + r**2 != 0:
        raise AssertionError("curve did not normalize to y^2+xy=x^3+A*x^2+b")
    return r, t, a2n, a6n


def verify_normalized_point(curve, point, r, t, a2, b):
    q = normalize_point(point, r, t)
    if q is None:
        return True
    x, y = q
    return y**2 + x*y == x**3 + a2*x**2 + b


def build_isogeny(source):
    trace = -534059
    frobenius_discriminant = trace**2 - 4 * 2**37
    if frobenius_discriminant != -7 * 194399**2 or 194399 != 73 * 2663:
        raise AssertionError("registered Frobenius-order conductor identity failed")
    if not Integer(230603167).is_prime(proof=True):
        raise AssertionError("registered factor-base subgroup order is not prime")
    splitting = int(kronecker_symbol(-7, 73))
    if splitting != -1:
        raise AssertionError("73 is not inert in the source CM field")
    isogenies = source.isogenies_prime_degree(73)
    if len(isogenies) != 74:
        raise AssertionError("inert prime should give 74 rational degree-73 kernels")
    # Sage's exact enumeration order is deterministic for a fixed field presentation.
    phi = isogenies[0]
    target = phi.codomain()
    kernel_poly = phi.kernel_polynomial()
    if phi.degree() != 73 or not phi.is_separable() or kernel_poly.degree() != 36:
        raise AssertionError("degree/separability/kernel polynomial certificate failed")
    if kernel_poly.gcd(kernel_poly.derivative()).degree() != 0:
        raise AssertionError("kernel polynomial is not squarefree")
    division_polynomial = source.division_polynomial(73)
    if division_polynomial % kernel_poly != 0:
        raise AssertionError("kernel polynomial does not divide the 73-division polynomial")
    quotient_degree = (division_polynomial // kernel_poly).degree()
    if quotient_degree != 2628:
        raise AssertionError("unexpected 73-division polynomial quotient degree")
    if source.cardinality() != 137439487532 or target.cardinality() != source.cardinality():
        raise AssertionError("source/target point-count certificate failed")
    set_random_seed(730037)
    homomorphism_pairs = 0
    for _ in range(128):
        p, q = source.random_point(), source.random_point()
        if phi(p + q) != phi(p) + phi(q):
            raise AssertionError("sampled homomorphism check failed")
        homomorphism_pairs += 1
    xmap, ymap = phi.rational_maps()
    nr, nt, na2, nb = normalized_coefficients(target)
    certificate = {
        "inventory_count": len(isogenies),
        "degree": int(phi.degree()),
        "separable": bool(phi.is_separable()),
        "source_coefficients": [str(a) for a in source.a_invariants()],
        "target_coefficients": [str(a) for a in target.a_invariants()],
        "target_normalized_coefficients": {
            "x_shift": str(nr), "y_shift": str(nt), "a2": str(na2), "a6": str(nb)
        },
        "kernel_polynomial": str(kernel_poly),
        "kernel_polynomial_sha256": digest_bytes(str(kernel_poly).encode()),
        "kernel_polynomial_degree": int(kernel_poly.degree()),
        "kernel_division_polynomial_quotient_degree": quotient_degree,
        "homomorphism_checked_pairs": homomorphism_pairs,
        "source_group_order": int(source.cardinality()),
        "target_group_order": int(target.cardinality()),
        "field_modulus": str(source.base_ring().modulus()),
        "x_rational_map": str(xmap),
        "y_rational_map": str(ymap),
        "source_j": str(source.j_invariant()),
        "target_j": str(target.j_invariant()),
        "trace_over_gf_2_37": trace,
        "frobenius_order_discriminant": frobenius_discriminant,
        "frobenius_order_conductor": 194399,
        "source_endomorphism_order_discriminant": -7,
        "target_endomorphism_order_discriminant": -7 * 73**2,
        "legendre_symbol_minus7_mod_73": splitting,
        "direction": "descending; 73 is inert in Q(sqrt(-7)) and divides the Frobenius-order conductor",
        "kernel_rationality": "kernel subgroup is defined over the base field; individual kernel points need not be",
    }
    return phi, certificate


def subgroup_generator(curve, order, seed):
    set_random_seed(seed)
    cardinality = Integer(curve.cardinality())
    if cardinality % order:
        raise AssertionError("registered subgroup order does not divide the curve order")
    cofactor = cardinality // order
    while True:
        point = cofactor * curve.random_point()
        if not point.is_zero():
            if order * point != curve(0):
                raise AssertionError("random subgroup generator has wrong order")
            return point


def factor_base(curve, gen, seed):
    rng = random.Random(seed)
    reps, scalars, seen = [], [], set()
    while len(reps) < 4:
        k = rng.randrange(1, 230603167)
        point = Integer(k) * gen
        x = point[0]
        if x not in seen:
            seen.add(x)
            reps.append(point)
            scalars.append(k)
    return reps, scalars


def signed_sum(curve, points):
    total = curve(0)
    for point, sign in points:
        total += point if sign == 1 else -point
    return total


def target_specs(contract, summands):
    return contract["target_patterns"][str(summands)]


def field_equations(curve, support, target, summands, ring, r, t, b, canonical):
    xs = [normalize_point(p, r, t)[0] for p in support]
    xt = None if target.is_zero() else normalize_point(target, r, t)[0]
    values = []
    slots_for_assignment = []
    for assignment in range(1 << ring.n):
        slots = tuple((assignment >> (2*j)) & 3 for j in range(summands))
        slots_for_assignment.append(slots)
        selected = [xs[k] for k in slots]
        if summands == 2:
            value = selected[0] + selected[1] if xt is None else s3(selected[0], selected[1], xt, b)
        elif xt is None:
            value = s3(selected[0], selected[1], selected[2], b)
        else:
            value = s4(selected[0], selected[1], selected[2], xt, b)
        values.append(integer_representation(value))
    generators = [ring.anf([(value >> bit) & 1 for value in values]) for bit in range(37)]
    if canonical:
        generators.append(ring.anf([int(tuple(sorted(slots)) != slots) for slots in slots_for_assignment]))
    return generators, slots_for_assignment


def truth_set(curve, support, target, slots_for_assignment):
    answer = {}
    for assignment, slots in enumerate(slots_for_assignment):
        valid_signs = []
        for signs in product((1, -1), repeat=len(slots)):
            points = [(support[index], sign) for index, sign in zip(slots, signs)]
            if signed_sum(curve, points) == target:
                valid_signs.append(signs)
        if valid_signs:
            answer[assignment] = valid_signs
    return answer


def run():
    contract_bytes = (HERE / "contract.json").read_bytes()
    contract = json.loads(contract_bytes)
    q = 2**37
    field = GF(q, name="a")
    source = EllipticCurve(field, [1, 0, 0, 0, 1])
    if source.cardinality() != contract["source_group_order"]:
        raise AssertionError("registered source curve order mismatch")
    phi, isogeny_certificate = build_isogeny(source)
    target_curve = phi.codomain()
    curves = [("source", source, None)] + [("degree_73", target_curve, phi)]
    data = []
    started = time.monotonic()
    transported_truth_hashes = {}
    for split, seed in contract["seeds"].items():
        generator = subgroup_generator(source, contract["factor_base_subgroup_order"], seed)
        support_source, support_scalars = factor_base(source, generator, seed ^ 0x5A17)
        targets_by_m = {}
        for summands in contract["summands"]:
            specs = target_specs(contract, summands)
            targets_by_m[summands] = [signed_sum(source, [(support_source[i], sign) for i, sign in spec])
                                      for spec in specs]
        source_mapped = [(name, curve, mapping) for name, curve, mapping in curves]
        transported_support = {
            name: [p if mapping is None else mapping(p) for p in support_source]
            for name, _, mapping in source_mapped
        }
        for name, curve, mapping in source_mapped:
            r, t, a2, b = normalized_coefficients(curve)
            if any(not verify_normalized_point(curve, p, r, t, a2, b) for p in transported_support[name]):
                raise AssertionError("normalization failed on factor-base points")
            if not verify_normalized_point(curve, curve(0), r, t, a2, b):
                raise AssertionError("normalization failed at infinity")
            for summands in contract["summands"]:
                for target_index, source_target in enumerate(targets_by_m[summands]):
                    target = source_target if mapping is None else mapping(source_target)
                    for canonical in (False, True):
                        ring = Ring(2 * summands)
                        generators, slots = field_equations(curve, transported_support[name], target,
                            summands, ring, r, t, b, canonical)
                        truth = truth_set(curve, transported_support[name], target, slots)
                        if canonical:
                            truth = {a: value for a, value in truth.items()
                                     if tuple(sorted(slots[a])) == slots[a]}
                        key = (split, summands, target_index, canonical)
                        truth_hash = digest(sorted((assignment, signs) for assignment, signs in truth.items()))
                        if name == "source":
                            transported_truth_hashes[key] = truth_hash
                        elif truth_hash != transported_truth_hashes[key]:
                            raise AssertionError("isogeny changed the transported signed workload")
                        repetitions = [measure(ring, generators) for _ in range(contract["repetitions"])]
                        if any(digest(rep) != digest(repetitions[0]) for rep in repetitions[1:]):
                            raise AssertionError("matrix solver replay was not deterministic")
                        measured = repetitions[0]
                        roots = set(measured["roots"])
                        if not set(truth).issubset(roots):
                            raise AssertionError("summation equations lost a genuine decomposition")
                        data.append({
                            "split": split, "seed": seed, "model": name, "summands": summands,
                            "target_index": target_index,
                            "target_is_infinity": bool(target.is_zero()),
                            "encoding": "canonical" if canonical else "ordered",
                            "support_scalars": support_scalars,
                            "target_pattern": target_specs(contract, summands)[target_index],
                            "target_x": None if target.is_zero() else str(normalize_point(target, r, t)[0]),
                            "input_sha256": digest(generators),
                            "generators": generators,
                            "truth_sha256": truth_hash,
                            "true_selector_assignments": len(truth),
                            "verified_signed_lifts": sum(len(signs) for signs in truth.values()),
                            "rejected_algebraic_roots": len(roots - set(truth)),
                            "completion_degree": measured["completion_degree"],
                            "f4_xors": sum(row["total_matrix_xors"] for row in measured["traces"]["f4"]),
                            "f5_xors": sum(row["total_matrix_xors"] for row in measured["traces"]["f5"]),
                            "matrix_rows": sum(row["matrix_rows"] for row in measured["traces"]["f5"]),
                            "criterion_rows": sum(row["criterion_rows"] for row in measured["traces"]["f5"]),
                            "repetition_sha256": [digest(rep) for rep in repetitions],
                            "status": measured["status"],
                        })
    return {
        "contract": contract,
        "contract_sha256": digest_bytes(contract_bytes),
        "source_sha256": {p.name: digest_bytes(p.read_bytes()) for p in
            [HERE / "experiment.py", HERE / "contract.json", REPO / "research/toy_f5_neighbors_20260924/matrix.py"]},
        "isogeny_certificate": isogeny_certificate,
        "cases": data,
        "elapsed_seconds": time.monotonic() - started,
        "accounting": {
            "matrix_stage_only": True,
            "full_dlp_total_operations": None,
            "S": None,
            "rho_ratio": None,
            "floor_ratio": None,
            "speedup": None,
            "matrix_solver_repetitions_are_independent_samples": False,
        },
    }


def summarize(result):
    groups = defaultdict(list)
    for case in result["cases"]:
        groups[(case["split"], case["summands"], case["encoding"], case["model"])].append(case)
    rows = []
    for (split, summands, encoding, model), cases in sorted(groups.items()):
        rows.append({
            "split": split, "summands": summands, "encoding": encoding, "model": model,
            "systems": len(cases),
            "completion_degree_max": max(c["completion_degree"] for c in cases),
            "completion_degree_sum": sum(c["completion_degree"] for c in cases),
            "f4_xors": sum(c["f4_xors"] for c in cases),
            "f5_xors": sum(c["f5_xors"] for c in cases),
            "matrix_rows": sum(c["matrix_rows"] for c in cases),
            "criterion_rows": sum(c["criterion_rows"] for c in cases),
            "all_verified": all(c["status"] == "VERIFIED" for c in cases),
        })
    indexed = {(c["split"], c["summands"], c["target_index"], c["encoding"], c["model"]): c
               for c in result["cases"]}
    paired_higher = 0
    for case in result["cases"]:
        if case["model"] == "degree_73":
            base = indexed[(case["split"], case["summands"], case["target_index"], case["encoding"], "source")]
            paired_higher += int(case["completion_degree"] > base["completion_degree"])
    ratios = []
    for row in rows:
        if row["model"] == "degree_73":
            base = next(x for x in rows if (x["split"], x["summands"], x["encoding"], x["model"])
                        == (row["split"], row["summands"], row["encoding"], "source"))
            ratios.append({"split": row["split"], "summands": row["summands"], "encoding": row["encoding"],
                           "f5_xor_ratio_to_source": row["f5_xors"] / base["f5_xors"]})
    passed = bool(ratios) and all(r["f5_xor_ratio_to_source"] <= 0.90 for r in ratios) and paired_higher == 0
    return {"rows": rows, "paired_completion_degree_increases": paired_higher,
            "f5_xor_ratios": ratios, "registered_success_rule_passed": passed,
            "classification": "engineering diagnostic" if passed else "registered criterion failed",
            "full_dlp_speedup": None, "systems": len(result["cases"]),
            "solver_replays": len(result["cases"]) * result["contract"]["repetitions"] * 2}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    result = run()
    (args.output / "raw.json.gz").write_bytes(gzip.compress(
        (json.dumps(result, separators=(",", ":")) + "\n").encode(), mtime=0))
    summary = summarize(result)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=2))
