#!/usr/bin/env python3
"""Exact ECC2K-130 degree-263 dual certificate with independent full-point replay."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
PRE = RESEARCH / "ecc2k130_direction_review_20260924"
REL = RESEARCH / "ecc2k130_relations"
sys.dont_write_bytecode = True
for directory in (HERE, PRE, REL):
    sys.path.insert(0, str(directory))

from dual_transport import DualTransport
from redteam_velu_replay import Field, Replay
from fastfield import FastGF2m, IRR131
from relations import (Koblitz, CHALLENGE_ELL, CHALLENGE_PX, CHALLENGE_PY,
                       CHALLENGE_QX, CHALLENGE_QY)

DEGREE = 263
FROZEN = {
    "twist_torsion_results.json": "6f1e7b22f3471214d38ec0d3196c88fe9edaf45791e9c05ea67764831fd75368",
    "redteam_velu_replay.py": "c08fbb8de9d54906d783be967a36a6073b2a1d1a3583292badd3adcdd7126790",
    "oriented_velu.py": "a1a656cc32efd612250b1ecf50dff18af79a8f62f970ee79841bfd390603423b",
    "fastfield.py": "b5990fda51700bbfba363251c41f02febfd0fc92edfea4f5ddd647bab472789e",
    "relations.py": "0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f",
}
FILES = {
    "twist_torsion_results.json": PRE / "twist_torsion_results.json",
    "redteam_velu_replay.py": PRE / "redteam_velu_replay.py",
    "oriented_velu.py": RESEARCH / "ecc2k130_oriented_transport_20260924" / "oriented_velu.py",
    "fastfield.py": REL / "fastfield.py",
    "relations.py": REL / "relations.py",
}
P = (CHALLENGE_PX, CHALLENGE_PY)
Q = (CHALLENGE_QX, CHALLENGE_QY)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def encode(point):
    return None if point is None else [hex(point[0]), hex(point[1])]


def neg(point):
    return None if point is None else (point[0], point[1] ^ point[0])


def half_kernel(reference, generator, a, b=1):
    point, result = generator, []
    for _ in range((DEGREE - 1) // 2):
        check(point is not None and reference.on_curve(point, a, b), "half-kernel point")
        result.append(point)
        point = reference.add(point, generator, a)
    check(point == neg(result[-1]), "half-kernel did not close")
    check(len({x for x, _ in result}) == (DEGREE - 1) // 2, "half-kernel duplicate abscissa")
    return result


def direct_rational_velu(reference, point, kernel, a):
    """Independent paired full Vélu sum for a rational same-curve kernel."""
    x, y = point
    total = 0
    for kernel_point in kernel:
        minus = neg(kernel_point)
        plus_image = reference.add(point, kernel_point, a)
        minus_image = reference.add(point, minus, a)
        check(plus_image is not None and minus_image is not None, "reference denominator exception")
        x ^= plus_image[0] ^ minus_image[0]
        y ^= plus_image[1] ^ minus_image[1] ^ kernel_point[1] ^ minus[1]
        total ^= kernel_point[0]
    return x, y ^ total


def delta(after, before):
    return {key: after[key] - before[key] for key in after}


def expected_rejection(operation):
    try:
        operation()
    except (ValueError, TypeError):
        return True
    raise AssertionError("invalid control was accepted")


class CountedField(FastGF2m):
    def __init__(self):
        self.counts = {}
        self._inside_inv = 0
        super().__init__(131, IRR131)
        self.reset()

    def reset(self):
        self.counts = {"mul": 0, "sqr": 0, "inv": 0,
                       "mul_inside_inv": 0, "sqr_inside_inv": 0}

    def mul(self, x, y):
        self.counts["mul"] += 1
        if self._inside_inv:
            self.counts["mul_inside_inv"] += 1
        return super().mul(x, y)

    def sqr(self, x):
        self.counts["sqr"] += 1
        if self._inside_inv:
            self.counts["sqr_inside_inv"] += 1
        return super().sqr(x)

    def inv(self, x):
        self.counts["inv"] += 1
        self._inside_inv += 1
        try:
            return super().inv(x)
        finally:
            self._inside_inv -= 1


def run_row(saved_run, saved, field, deadline):
    check(time.monotonic() < deadline, "300-second stop")
    source, twist = Koblitz(field, a=0), Koblitz(field, a=1)
    generator = tuple(int(z, 16) for z in saved["kernel_generator_on_twist"])
    complement = tuple(int(z, 16) for z in saved_run["basis_on_twist"][1])
    check(saved["line"][0] == 1, "frozen complement rule no longer applies")
    check(source.on_curve(P) and source.on_curve(Q), "public source points")
    field.reset()
    started = time.perf_counter()
    transport = DualTransport(source, twist, generator, complement, DEGREE, P)
    setup_seconds = time.perf_counter() - started
    setup_counts = dict(field.counts)
    check(transport.raw_reverse_scalar_sign == -1, "frozen normalized sign changed")
    check(transport.raw_reverse.codomain.b == source.b == 1, "reverse curve model")
    check(transport.forward.codomain.b == int(saved["codomain_b"], 16), "forward curve model")
    check(transport.forward_twist.codomain.b == transport.forward.codomain.b, "twist quotient b")
    check(transport.source.mul(P, CHALLENGE_ELL) is None, "P subgroup")
    check(transport.source.mul(Q, CHALLENGE_ELL) is None, "Q subgroup")

    field.reset()
    warm = transport.compose(P)
    composition_counts = dict(field.counts)
    check(warm == source.mul(P, DEGREE), "warm composition")
    composition_times = []
    for _ in range(11):
        check(time.monotonic() < deadline, "300-second stop")
        field.reset()
        started = time.perf_counter()
        check(transport.compose(P) == warm, "repeated composition changed")
        composition_times.append(time.perf_counter() - started)
        check(field.counts == composition_counts, "composition count changed")

    # This independent group law imports no production arithmetic. In
    # particular, the complementary-kernel image is reconstructed by direct
    # rational full Vélu sums, not by trusting the production twist map.
    started = time.perf_counter()
    reference_field = Field(131, IRR131)
    reference = Replay(reference_field)
    check(reference.on_curve(P) and reference.on_curve(Q), "reference source membership")
    check(reference.on_curve(generator, 1) and reference.on_curve(complement, 1), "reference twist membership")
    check(reference.mul(generator, DEGREE, 1) is None, "reference forward kernel order")
    check(reference.mul(complement, DEGREE, 1) is None, "reference complement order")
    kernel = half_kernel(reference, generator, 1)
    independent_complement_image = direct_rational_velu(reference, complement, kernel, 1)
    check(independent_complement_image == transport.reverse_kernel_generator,
          "independent reverse kernel generator")
    dest_b = transport.forward.codomain.b
    check(reference.on_curve(independent_complement_image, 1, dest_b), "reference reverse twist membership")
    check(reference.mul(independent_complement_image, DEGREE, 1) is None,
          "reference reverse kernel order")
    reverse_kernel = half_kernel(reference, independent_complement_image, 1, dest_b)
    check({x for x, _ in kernel} == transport.forward.kernel_abscissae, "forward kernel abscissae")
    check({x for x, _ in reverse_kernel} == transport.raw_reverse.kernel_abscissae,
          "reverse kernel abscissae")
    controls = {"P": P, "Q": Q, "P+Q": reference.add(P, Q),
                "2P": reference.add(P, P)}
    full_rows = {}
    for label, point in controls.items():
        check(reference.on_curve(point), "reference source control")
        ref_forward = reference.direct_velu(point, kernel)
        prod_forward = transport.forward(point)
        check(ref_forward == prod_forward, "independent forward full point " + label)
        check(reference.on_curve(ref_forward, 0, dest_b), "forward codomain " + label)
        ref_raw_reverse = reference.direct_velu(ref_forward, reverse_kernel)
        prod_raw_reverse = transport.raw_reverse(prod_forward)
        check(ref_raw_reverse == prod_raw_reverse, "independent reverse full point " + label)
        expected = reference.mul(point, DEGREE)
        check(ref_raw_reverse == neg(expected), "raw normalized reverse sign " + label)
        check(transport.dual(prod_forward) == expected, "positive dual identity " + label)
        check(transport.compose(point) == expected, "composition " + label)
        check(transport.forward(transport.dual(prod_forward)) == reference.mul(ref_forward, DEGREE),
              "other dual identity " + label)
        full_rows[label] = {"forward": encode(ref_forward),
                            "raw_reverse": encode(ref_raw_reverse),
                            "positive_dual": encode(expected)}
    scalar_checks = []
    for scalar in (0, 1, -1, 2, 3, 17, 263):
        point = reference.mul(P, scalar) if scalar >= 0 else neg(reference.mul(P, -scalar))
        expected = reference.mul(point, DEGREE)
        check(transport.compose(point) == expected, "signed scalar composition")
        scalar_checks.append(scalar)
    check(transport.forward(None) is None and transport.forward_twist(None) is None,
          "forward infinity")
    check(transport.raw_reverse(None) is None and transport.raw_reverse_twist(None) is None,
          "reverse infinity")
    check(transport.dual(None) is None and transport.compose(None) is None,
          "dual infinity")
    # The untwisted rational groups have no nonzero 263-kernel points. The
    # twist-side maps provide all actual finite exceptional-input controls.
    exception_counts = {"forward_twist_kernel_points": 0,
                        "reverse_twist_kernel_points": 0,
                        "offcurve_kernel_abscissa_rejections": 0,
                        "invalid_complement_rejections": 0,
                        "invalid_witness_rejections": 0}
    for half, mapped, key in ((kernel, transport.forward_twist, "forward_twist_kernel_points"),
                              (reverse_kernel, transport.raw_reverse_twist, "reverse_twist_kernel_points")):
        for point in half:
            for actual in (point, neg(point)):
                check(mapped(actual) is None, "twist kernel did not map to infinity")
                exception_counts[key] += 1
    for mapped in (transport.forward, transport.raw_reverse):
        abscissa = next(iter(mapped.kernel_abscissae))
        bad_y = 0
        while mapped.source.on_curve((abscissa, bad_y)):
            bad_y += 1
        check(expected_rejection(lambda m=mapped, x=abscissa, y=bad_y: m((x, y))),
              "offcurve kernel-abscissa input")
        exception_counts["offcurve_kernel_abscissa_rejections"] += 1
    check(expected_rejection(lambda: DualTransport(source, twist, generator, generator,
                                                   DEGREE, P)), "non-complement rejection")
    exception_counts["invalid_complement_rejections"] += 1
    check(expected_rejection(lambda: DualTransport(source, twist, generator, complement,
                                                   DEGREE, (0, 1))), "ambiguous witness rejection")
    exception_counts["invalid_witness_rejections"] += 1
    validation_replay_seconds = time.perf_counter() - started
    reference_counts = dict(reference_field.counts)
    return {
        "seed": saved_run["seed"], "line": saved["line"],
        "direction": "horizontal" if saved["orbit_length"] == 1 else "descending",
        "forward_b": hex(dest_b), "reverse_b": hex(transport.raw_reverse.codomain.b),
        "reverse_kernel_generator": encode(independent_complement_image),
        "raw_reverse_scalar_sign": transport.raw_reverse_scalar_sign,
        "formal_dual_correction": "negate raw reverse output",
        "independent_full_coordinate_rows": full_rows,
        "scalar_controls": scalar_checks,
        "exception_counts": exception_counts,
        "setup_from_saved_generators_seconds": setup_seconds,
        "setup_from_saved_generators_field_counts": setup_counts,
        "production_composition_field_counts": composition_counts,
        "production_composition_seconds_11": composition_times,
        "production_composition_seconds_median": statistics.median(composition_times),
        "validation_replay_seconds": validation_replay_seconds,
        "independent_reference_field_counts": reference_counts,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("preserve prior evidence: choose a fresh --out")
    started = time.monotonic()
    receipt = {"schema": "ecc2k130_263_dual_composition_v1", "status": "RUNNING",
               "source_sha256": {"certificate.py": sha(__file__),
                                 "dual_transport.py": sha(HERE / "dual_transport.py")},
               "input_sha256": {name: sha(path) for name, path in FILES.items()},
               "host": platform.platform(), "python": sys.version,
               "command": [sys.executable, *sys.argv],
               "PDP_cost": None, "full_ECDLP_cost": None, "end_to_end_speedup": None,
               "cost_scope": "saved torsion generators supplied; no cold basis discovery or relation/PDP work"}
    try:
        check(receipt["input_sha256"] == FROZEN, "frozen input hash changed")
        frozen = json.loads(FILES["twist_torsion_results.json"].read_text())
        begin = time.perf_counter()
        field = CountedField()
        receipt["production_field_initialization_seconds"] = time.perf_counter() - begin
        rows = []
        receipt["rows"] = rows
        for saved_run in frozen["runs"]:
            for saved in saved_run["representative_x_maps"]:
                rows.append(run_row(saved_run, saved, field, started + 300))
        check(len(rows) == 8, "expected eight saved representative maps")
        stable = [{key: row[key] for key in (
            "seed", "line", "direction", "forward_b", "reverse_b",
            "reverse_kernel_generator", "raw_reverse_scalar_sign",
            "formal_dual_correction", "independent_full_coordinate_rows",
            "scalar_controls", "exception_counts", "setup_from_saved_generators_field_counts",
            "production_composition_field_counts", "independent_reference_field_counts")}
            for row in rows]
        receipt["deterministic_sha256"] = hashlib.sha256(
            json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        check(time.monotonic() < started + 300, "300-second stop")
        receipt["status"] = "PASS"
    except Exception:
        receipt["status"] = "FAIL"
        receipt["failure"] = traceback.format_exc()
    receipt["elapsed_seconds"] = time.monotonic() - started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "seconds": receipt["elapsed_seconds"],
                      "out": str(args.out)}, indent=2))
    if receipt["status"] != "PASS":
        print(receipt["failure"], file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
