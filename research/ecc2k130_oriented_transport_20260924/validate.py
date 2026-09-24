#!/usr/bin/env python3
"""Replayable correctness checks; saves a fresh PASS or FAIL receipt."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import traceback

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent / "ecc2k130_direction_review_20260924"))
from oriented_velu import BinaryVeluMap
from fastfield import FastGF2m, IRR131
from relations import (Koblitz, CHALLENGE_ELL, CHALLENGE_PX, CHALLENGE_PY,
                       CHALLENGE_QX, CHALLENGE_QY)
from redteam_velu_replay import Field, Replay


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encode(point):
    return None if point is None else [hex(z) for z in point]


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def points(curve):
    return [None] + [point for x in range(1 << curve.F.deg) for point in curve.points_over(x)]


def largest_odd_generator(curve, all_points):
    cofactor, odd = 1, len(all_points)
    while odd % 2 == 0:
        cofactor, odd = 2*cofactor, odd//2
    best = (1, None)
    for point in all_points[1:]:
        generator = curve.mul(point, cofactor)
        if generator is None:
            continue
        current, degree = generator, 1
        while current is not None:
            current = curve.add(current, generator)
            degree += 1
            check(degree <= odd, "toy subgroup did not close")
        if degree > best[0]:
            best = degree, generator
        if degree == odd:
            break
    return best


def expect_rejection(operation):
    try:
        operation()
    except (ValueError, TypeError):
        return 1
    raise AssertionError("invalid input was accepted")


def toy_checks(receipt, deadline):
    rows = []
    receipt["toy_rows"] = rows
    for n, irr in ((3, 0b1011), (5, 0b100101), (7, 0b10000011)):
        field = FastGF2m(n, irr)
        curves = [Koblitz(field, a=a) for a in (0, 1)]
        all_points = [points(curve) for curve in curves]
        for kernel_a in (0, 1):
            kernel_curve = curves[kernel_a]
            degree, generator = largest_odd_generator(kernel_curve, all_points[kernel_a])
            if generator is None:
                continue
            for source_a in (0, 1):
                check(time.monotonic() < deadline, "validation time limit")
                source = curves[source_a]
                phi = BinaryVeluMap.from_generator(source, kernel_curve, generator, degree)
                inputs = all_points[source_a]
                images = {point: phi(point) for point in inputs}
                check(images[None] is None, "infinity image")
                check(all(phi.codomain.on_curve(image) for image in images.values()), "toy codomain")
                check(all(images[source.neg(p)] == phi.codomain.neg(images[p]) for p in inputs), "toy negation")
                for p in inputs:
                    for q in inputs:
                        check(images[source.add(p, q)] == phi.codomain.add(images[p], images[q]), "toy additivity")
                rational_kernel = [p for p in inputs if p is not None and p[0] in phi.kernel_abscissae]
                check(all(images[p] is None for p in rational_kernel), "rational kernel image")
                check(images[(0, 1)] == (0, field.frobenius(phi.codomain.b, n-1)), "2-torsion image")
                rejected = 0
                for malformed in ([0, 1], (0, 0), (1 << n, 1), (True, 1), (0,), "infinity"):
                    rejected += expect_rejection(lambda p=malformed: phi(p))
                u = next(iter(phi.kernel_abscissae))
                invalid_y = next(y for y in range(1 << n) if not source.on_curve((u, y)))
                rejected += expect_rejection(lambda: phi((u, invalid_y)))
                rejected += expect_rejection(lambda: BinaryVeluMap.from_generator(source, kernel_curve, None, degree))
                rejected += expect_rejection(lambda: BinaryVeluMap.from_generator(source, kernel_curve, generator, 2))
                rejected += expect_rejection(lambda: BinaryVeluMap.from_generator(source, kernel_curve, generator, degree*3))
                rejected += expect_rejection(lambda: BinaryVeluMap.from_generator(source, Koblitz(field, b=0), generator, degree))

                # Where dual kernels also have rational abscissae, construct
                # reverse quotients and check the entire toy group. A normalized
                # dual may differ by negation, which is recorded explicitly.
                compositions, seen = [], set()
                if n <= 5:
                    for reverse_a in (0, 1):
                        reverse_kernel = Koblitz(field, a=reverse_a, b=phi.codomain.b)
                        for candidate in points(reverse_kernel)[1:]:
                            if reverse_kernel.mul(candidate, degree) is not None:
                                continue
                            try:
                                reverse = BinaryVeluMap.from_generator(phi.codomain, reverse_kernel, candidate, degree)
                            except ValueError:
                                continue
                            if reverse.kernel_abscissae in seen:
                                continue
                            seen.add(reverse.kernel_abscissae)
                            if reverse.codomain.b != source.b:
                                continue
                            composed = phi.then(reverse)
                            got = [composed(p) for p in inputs]
                            expected = [source.mul(p, degree) for p in inputs]
                            sign = 1 if got == expected else -1 if got == [source.neg(p) for p in expected] else 0
                            if sign:
                                check(composed.degree == degree*degree, "composition degree")
                                compositions.append({"reverse_kernel_a": reverse_a, "scalar_sign": sign,
                                                     "all_point_equalities": len(inputs)})
                    check(compositions, "expected a normalized dual in tiny control")
                rows.append({"n": n, "modulus": hex(irr), "source_a": source_a,
                             "kernel_a": kernel_a, "degree": degree, "source_order": len(inputs),
                             "codomain_b": hex(phi.codomain.b), "codomain_checks": len(inputs),
                             "full_additivity_checks": len(inputs)**2,
                             "rational_kernel_inputs": len(rational_kernel), "infinity_checks": 1,
                             "two_torsion_checks": 1, "rejected_invalid_inputs": rejected,
                             "normalized_dual_compositions": compositions})


def generic_interface_checks(receipt):
    field = FastGF2m(5, 0b100101)
    for b in range(2, 32):
        source = Koblitz(field, a=0, b=b)
        inputs = points(source)
        degree, generator = largest_odd_generator(source, inputs)
        if degree > 3 and any(degree % d == 0 for d in range(2, degree)):
            break
    else:
        raise AssertionError("no composite cyclic control found")
    phi = BinaryVeluMap.from_generator(source, source, generator, degree)
    kernel, current, total = [], generator, 0
    for _ in range((degree-1)//2):
        kernel.append(current)
        total ^= current[0]
        current = source.add(current, generator)
    images = {p: phi(p) for p in inputs}
    direct_checks = 0
    for p in inputs:
        if p is None or p[0] in phi.kernel_abscissae:
            continue
        x, y = p
        for point in kernel:
            negative = source.neg(point)
            plus, minus = source.add(p, point), source.add(p, negative)
            x ^= plus[0] ^ minus[0]
            y ^= plus[1] ^ minus[1] ^ point[1] ^ negative[1]
        check(images[p] == (x, y ^ total), "generic direct rational Velu")
        direct_checks += 1
    for p in inputs:
        for q in inputs:
            check(images[source.add(p, q)] == phi.codomain.add(images[p], images[q]), "composite additivity")
    invalid = 0
    invalid += expect_rejection(lambda: BinaryVeluMap.from_generator(source, Koblitz(field, b=b ^ 1), generator, degree))
    invalid += expect_rejection(lambda: BinaryVeluMap.from_generator(source, Koblitz(FastGF2m(3, 0b1011)), generator, degree))
    invalid += expect_rejection(lambda: BinaryVeluMap.from_generator(source, source, generator, degree+2))
    other = BinaryVeluMap.from_generator(Koblitz(field, a=1, b=b), source, generator, degree)
    invalid += expect_rejection(lambda: phi.then(other))
    receipt["generic_interface_control"] = {"n": 5, "source_a": 0, "source_b": b,
        "source_order": len(inputs), "degree": degree,
        "direct_rational_coordinate_checks": direct_checks,
        "full_additivity_checks": len(inputs)**2, "parameter_rejections": invalid}


def exact_checks(receipt, frozen, deadline):
    p = CHALLENGE_PX, CHALLENGE_PY
    q = CHALLENGE_QX, CHALLENGE_QY
    rows = []
    receipt["exact_rows"] = rows
    for run in frozen["runs"]:
        field = FastGF2m(131, IRR131)
        source, twist = Koblitz(field, a=0), Koblitz(field, a=1)
        independent = Replay(Field(131, IRR131))
        for saved in run["representative_x_maps"]:
            check(time.monotonic() < deadline, "validation time limit")
            generator = tuple(int(v, 16) for v in saved["kernel_generator_on_twist"])
            phi = BinaryVeluMap.from_generator(source, twist, generator, 263)
            check(phi.kernel_abscissae == {int(x, 16) for x in saved["kernel_abscissae"]}, "saved kernel")
            check(phi.codomain.b == int(saved["codomain_b"], 16), "saved codomain")
            # Replay builds kernel points with a separate field and group implementation.
            kernel, current = [], generator
            for _ in range(131):
                check(independent.on_curve(current, 1), "independent kernel curve")
                kernel.append(current)
                current = independent.add(current, generator, 1)
            check(independent.mul(generator, 263, 1) is None, "independent kernel order")
            targets = {"P": p, "Q": q, "P+Q": source.add(p, q), "2P": source.add(p, p),
                       "2-torsion": (0, 1)}
            images = {name: phi(point) for name, point in targets.items()}
            for name, point in targets.items():
                check(images[name] == independent.direct_velu(point, kernel), "independent full coordinates "+name)
            check(images["P"][0] == int(saved["image_P_x"], 16), "saved P abscissa")
            check(images["Q"][0] == int(saved["image_Q_x"], 16), "saved Q abscissa")
            check(phi(None) is None, "exact infinity")
            for scalar in (0, 1, -1, 2, 3, 17, 263):
                check(phi(source.mul(p, scalar)) == phi.codomain.mul(images["P"], scalar), "oriented scalar transport")
            check(phi.codomain.mul(images["P"], CHALLENGE_ELL) is None, "P image order")
            check(phi.codomain.mul(images["Q"], CHALLENGE_ELL) is None, "Q image order")
            check(images["P+Q"] == phi.codomain.add(images["P"], images["Q"]), "exact additivity")
            check(images["P+Q"] != phi.codomain.add(images["P"], phi.codomain.neg(images["Q"])), "negative sign control")
            # Same kernel on its rational twist exercises all 262 exceptional
            # finite inputs at the actual degree, not only tiny toy kernels.
            twist_phi = BinaryVeluMap.from_generator(twist, twist, generator, 263)
            for point in kernel:
                check(twist_phi(point) is None and twist_phi(twist.neg(point)) is None, "exact rational kernel")
            check(twist_phi(None) is None, "twist infinity")
            bad_y = 0
            u = next(iter(phi.kernel_abscissae))
            while source.on_curve((u, bad_y)):
                bad_y += 1
            expect_rejection(lambda: phi((u, bad_y)))
            rows.append({"seed": run["seed"], "line": saved["line"], "orbit_length": saved["orbit_length"],
                         "codomain_b": hex(phi.codomain.b),
                         "oriented_images": {name: encode(image) for name, image in images.items()},
                         "independent_full_coordinate_checks": len(targets),
                         "planted_scalar_checks": [0, 1, -1, 2, 3, 17, 263],
                         "subgroup_order_checks": 2, "full_additivity_checks": 1,
                         "negative_sign_controls": 1, "infinity_checks": 2,
                         "rational_kernel_inputs_on_twist": 262,
                         "invalid_kernel_abscissa_inputs_rejected": 1})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=ROOT.parent/"ecc2k130_direction_review_20260924"/"twist_torsion_results.json")
    args = parser.parse_args()
    if args.out.exists():
        parser.error("preserve earlier evidence: choose a fresh --out")
    start = time.monotonic()
    files = [ROOT/"PROTOCOL.md", ROOT/"oriented_velu.py", Path(__file__), args.input,
             ROOT.parent/"ecc2k130_relations"/"fastfield.py",
             ROOT.parent/"ecc2k130_relations"/"relations.py",
             ROOT.parent/"ecc2k130_direction_review_20260924"/"redteam_velu_replay.py"]
    receipt = {"schema": "oriented-binary-velu-validation-v1", "status": "RUNNING",
               "timestamp_utc": datetime.now(timezone.utc).isoformat(),
               "command": [sys.executable, *sys.argv], "host": platform.platform(), "python": sys.version,
               "sha256": {str(path.relative_to(ROOT.parent)): sha(path) for path in files},
               "full_ECDLP_cost": None, "PDP_cost": None, "end_to_end_speedup": None}
    try:
        toy_checks(receipt, start+180)
        generic_interface_checks(receipt)
        exact_checks(receipt, json.loads(args.input.read_text()), start+180)
        check(time.monotonic() < start+180, "validation time limit")
        receipt["status"] = "PASS"
    except Exception:
        receipt["status"] = "FAIL"
        receipt["failure"] = traceback.format_exc()
    receipt["elapsed_seconds"] = time.monotonic()-start
    with args.out.open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print(json.dumps({key: receipt[key] for key in ("status", "elapsed_seconds")}, indent=2))
    if receipt["status"] != "PASS":
        print(receipt["failure"], file=sys.stderr)
        return 1
    print(json.dumps({"toy_maps": len(receipt["toy_rows"]), "exact_maps": len(receipt["exact_rows"]),
                      "toy_additivity_equalities": sum(row["full_additivity_checks"] for row in receipt["toy_rows"]),
                      "exact_kernel_inputs": sum(row["rational_kernel_inputs_on_twist"] for row in receipt["exact_rows"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
