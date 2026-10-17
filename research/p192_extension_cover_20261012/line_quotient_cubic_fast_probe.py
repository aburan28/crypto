"""Exact prime-field cubic cuts for frozen P-192 pivot quotient points."""

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, matrix, prod, vector

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_prelude(name):
    path = HERE / name
    source = path.read_text()
    marker = "start = time.monotonic_ns()"
    assert source.count(marker) == 1
    namespace = {"__file__": str(path), "__name__": name[:-3] + "_prelude"}
    exec(compile(source.split(marker)[0], str(path), "exec"), namespace)
    return namespace


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=16)
    args = parser.parse_args()
    frozen_path = HERE / "type_i_two_pivot_reduced64.json"
    frozen = json.loads(frozen_path.read_text())
    assert 1 <= args.count <= 64 and frozen["reduced_classes"] == 64
    old = load_prelude("line_quotient_cubic_probe.py")
    dt, E, L, p, r = (old[name] for name in ("dt", "E", "L", "p", "r"))
    K = GF(p)
    R = PolynomialRing(K, "x")
    x = R.gen()
    S = PolynomialRing(R, "y")
    y = S.gen()
    plane_path = HERE / "type_i_plane_quartic.json"
    plane = json.loads(plane_path.read_text())
    curve = sum(S(K(int(c)) * x ** a) * y ** b
                for label, c in plane["quartic_coefficients"].items()
                for a, b, _ in [tuple(map(int, label.split(",")))])
    anchors = [tuple(L(int(v)) for v in item["coordinates_Fp_decimal"])
               for item in plane["known_points_verified"]
               if item["source"] == "base"][:3]
    monomials = old["monomials"]
    assert len(anchors) == 3 and len(monomials) == 10 and L.degree() == 6

    def encode(poly):
        return [[str(int(c))] for c in poly.list()]

    def descend(poly):
        assert all(c ** p == c for c in poly.list())
        return R([K(int(c.polynomial()[0])) for c in poly.list()])

    def fast_residual(point):
        fiber = dt.fiber(point)
        coordinates = []
        for cover in fiber:
            ys = dt.check_model(cover)
            assert all(value != 0 for value in ys)
            coordinates.append(tuple(dt.Vandermonde.solve_right(
                vector(L, [1 / value for value in ys]))))
        coordinates.extend(anchors)
        assert len(coordinates) == 5 and all(c[2] != 0 for c in coordinates)
        rows = []
        for coordinate in coordinates:
            evaluations = [prod(coordinate[k] ** e[k] for k in range(3))
                           for e in monomials]
            for basis_index in range(6):
                row = [K(value.polynomial()[basis_index]) for value in evaluations]
                if any(row):
                    rows.append(row)
        evaluation = matrix(K, rows)
        assert evaluation.rank() == 9
        kernel = list(evaluation.right_kernel().basis())
        assert len(kernel) == 1
        coefficients = list(kernel[0])
        leading = next(value for value in coefficients if value != 0)
        coefficients = [value / leading for value in coefficients]
        cubic = sum(S(coefficients[i] * x ** e[0]) * y ** e[1]
                    for i, e in enumerate(monomials))
        LR = PolynomialRing(L, "x")
        lx = LR.gen()
        known_l = prod((lx - (c[0] / c[2]) ** (p ** power)
                        for c in coordinates[:2] for power in range(3)), LR.one())
        known_l *= prod((lx - c[0] / c[2] for c in anchors), LR.one())
        known = descend(known_l)
        resultant = curve.resultant(cubic)
        quotient, remainder = resultant.quo_rem(known)
        assert resultant.degree() == 12 and remainder == 0
        assert known.degree() == 9 and known.gcd(known.derivative()).degree() == 0
        h = quotient.monic()
        assert h.degree() == 3 and h.gcd(h.derivative()).degree() == 0
        assert h.gcd(known).degree() == 0
        linear = next(s for s in curve.subresultants(cubic) if s.degree() == 1)
        assert linear[1].gcd(h).degree() == 0
        y_value = (-linear[0] * linear[1].inverse_mod(h)) % h
        assert y_value.degree() < 3
        assert curve(y_value) % h == 0 and cubic(y_value) % h == 0
        return {
            "cubic_coefficients": [[str(int(c))] for c in coefficients],
            "residual_x_polynomial": encode(h),
            "residual_y_polynomial": encode(y_value),
        }

    def hybrid_residual(point):
        fiber = dt.fiber(point)
        support = [tuple(value ** (p ** power) for value in cover)
                   for power in range(3) for cover in fiber]
        coordinates = []
        for cover in support:
            ys = dt.check_model(cover)
            assert all(value != 0 for value in ys)
            coordinates.append(tuple(dt.Vandermonde.solve_right(
                vector(L, [1 / value for value in ys]))))
        coordinates.extend(anchors)
        evaluations = [[prod(c[k] ** e[k] for k in range(3))
                        for e in monomials] for c in coordinates]
        M = matrix(L, evaluations)
        assert M.rank() == 9
        coefficients_l = list(M.right_kernel().basis()[0])
        leading = next(value for value in coefficients_l if value != 0)
        coefficients_l = [value / leading for value in coefficients_l]
        assert all(value ** p == value for value in coefficients_l)
        coefficients = [K(int(value.polynomial()[0])) for value in coefficients_l]
        cubic = sum(S(coefficients[i] * x ** e[0]) * y ** e[1]
                    for i, e in enumerate(monomials))
        LR = PolynomialRing(L, "x")
        lx = LR.gen()
        known = descend(prod((lx - c[0] / c[2] for c in coordinates), LR.one()))
        resultant = curve.resultant(cubic)
        quotient, remainder = resultant.quo_rem(known)
        assert resultant.degree() == 12 and remainder == 0
        h = quotient.monic()
        assert h.degree() == 3 and h.gcd(h.derivative()).degree() == 0
        assert h.gcd(known).degree() == 0
        linear = next(s for s in curve.subresultants(cubic) if s.degree() == 1)
        assert linear[1].gcd(h).degree() == 0
        y_value = (-linear[0] * linear[1].inverse_mod(h)) % h
        assert curve(y_value) % h == 0 and cubic(y_value) % h == 0
        return {
            "cubic_coefficients": [[str(int(c))] for c in coefficients],
            "residual_x_polynomial": encode(h),
            "residual_y_polynomial": encode(y_value),
        }

    pairs = []
    for index in range(args.count):
        row = frozen["rows"][index]
        xx, yy = map(int, row["quotient"])
        Q = E(L(xx), L(yy))
        H = old["inv2"] * Q
        assert 2 * H == Q and r * Q == E(0)
        measurements = {}
        order = ("fast", "hybrid", "original") if index % 3 == 0 else (
            ("original", "fast", "hybrid") if index % 3 == 1 else
            ("hybrid", "original", "fast"))
        for label, point in (("positive_input_residual", Q),
                             ("negative_input_residual", H)):
            for method in order:
                start = time.monotonic_ns()
                result = {"fast": fast_residual, "hybrid": hybrid_residual,
                          "original": old["residual"]}[method](point)
                measurements[label + "_" + method + "_ns"] = time.monotonic_ns() - start
                assert result == row[label]
        measurements["order"] = list(order)
        pairs.append({"index": index, **measurements})
        print("FAST_CUBIC_PAIR_PASS", index, flush=True)

    fast_ns = sum(item[key] for item in pairs for key in
                  ("positive_input_residual_fast_ns",
                   "negative_input_residual_fast_ns"))
    hybrid_ns = sum(item[key] for item in pairs for key in
                    ("positive_input_residual_hybrid_ns",
                     "negative_input_residual_hybrid_ns"))
    original_ns = sum(item[key] for item in pairs for key in
                      ("positive_input_residual_original_ns",
                       "negative_input_residual_original_ns"))
    output = {
        "status": "matched_prime_field_cubic_reduction_probe",
        "field_prime": str(p), "subgroup_order": str(r),
        "classes": len(pairs), "residual_cuts": 2 * len(pairs),
        "pairs": pairs, "fast_total_ns_unisolated": fast_ns,
        "hybrid_total_ns_unisolated": hybrid_ns,
        "original_total_ns_unisolated": original_ns,
        "source_sha256": sha(__file__), "frozen_sha256": sha(frozen_path),
        "baseline_source_sha256": sha(HERE / "line_quotient_cubic_probe.py"),
        "plane_sha256": sha(plane_path),
        "runtime_info_sha256": sha(HERE / "type_i_fast_cubic_runtime_info.json"),
        "boundary": "Exact same-input cubic coefficients and degree-three residuals. Sequential shared-host stage times include fiber lifts and are not isolated wall-time speedups.",
    }
    path = HERE / "type_i_fast_cubic_probe.json"
    path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print("FAST_CUBIC_PASS", len(pairs), fast_ns, hybrid_ns, original_ns, flush=True)


if __name__ == "__main__":
    main()
