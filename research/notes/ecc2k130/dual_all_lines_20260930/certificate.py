#!/usr/bin/env python3
"""Frozen full-line degree-263 dual and exceptional-input certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import resource
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RESEARCH = REPO / "research"
for directory in (HERE, RESEARCH / "ecc2k130_direction_review_20260924",
                  RESEARCH / "ecc2k130_relations"):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
sys.dont_write_bytecode = True

from line_dual import DEGREE, build_line_dual
from fastfield import FastGF2m, IRR131
from redteam_velu_replay import Field, Replay
from relations import (CHALLENGE_PX, CHALLENGE_PY, CHALLENGE_QX,
                       CHALLENGE_QY, Koblitz)

P, Q = (CHALLENGE_PX, CHALLENGE_PY), (CHALLENGE_QX, CHALLENGE_QY)
LINES = [(1, t) for t in range(DEGREE)] + [(0, 1)]
REFERENCE_LINES = [(1, t) for t in (0, 1, 2, 131, 262)] + [(0, 1)]
WALL_LIMIT_SECONDS = 2100
MEMORY_LIMIT_BYTES = 2 * 1024**3


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def encoded(point):
    return None if point is None else [hex(point[0]), hex(point[1])]


def neg(point):
    return None if point is None else (point[0], point[0] ^ point[1])


def peak_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if platform.system() == "Darwin" else value * 1024


def within_budget(deadline):
    if time.monotonic() >= deadline:
        raise TimeoutError("2,100-second campaign wall limit")
    if peak_bytes() > MEMORY_LIMIT_BYTES:
        raise MemoryError("2 GiB campaign peak-memory limit")


def reject(call, label):
    try:
        call()
    except (ValueError, TypeError):
        return label
    raise AssertionError("invalid control accepted: " + label)


def basis_from_run(run):
    return tuple(tuple(int(z, 16) for z in point) for point in run["basis_on_twist"])


def kernel_exceptions(curve, generator, mapped):
    """Exercise every nonzero point, including both signs, on a twist kernel."""
    point, seen, count = generator, set(), 0
    for _ in range((DEGREE - 1) // 2):
        check(point is not None and curve.on_curve(point), "twist kernel membership")
        check(point[0] not in seen, "twist kernel duplicate abscissa")
        seen.add(point[0])
        for actual in (point, curve.neg(point)):
            check(mapped(actual) is None, "twist kernel did not map to infinity")
            count += 1
        previous = point
        point = curve.add(point, generator)
    check(point == curve.neg(previous), "twist kernel did not close")
    check(seen == mapped.kernel_abscissae, "twist kernel abscissae disagree")
    return count


def offcurve_kernel_rejection(mapped):
    x = min(mapped.kernel_abscissae)
    y = 0
    while mapped.source.on_curve((x, y)):
        y += 1
    reject(lambda: mapped((x, y)), "off-curve kernel-abscissa point")


def run_line(source, twist, basis, line, deadline):
    within_budget(deadline)
    begun = time.perf_counter()
    transport = build_line_dual(source, twist, basis, line, P)
    setup_seconds = time.perf_counter() - begun
    check(transport.forward.codomain.b == transport.forward_twist.codomain.b,
          "forward twist quotient model")
    check(transport.raw_reverse.codomain.a == source.a
          and transport.raw_reverse.codomain.b == source.b,
          "reverse quotient model")
    check(transport.forward(None) is None and transport.dual(None) is None
          and transport.compose(None) is None
          and transport.forward_twist(None) is None
          and transport.raw_reverse_twist(None) is None,
          "infinity transport")
    point_rows = {}
    for label, point in (("P", P), ("Q", Q)):
        forward = transport.forward(point)
        expected = source.mul(point, DEGREE)
        check(forward is not None, "source witness killed by twist kernel")
        check(transport.dual(forward) == expected,
              "positive dual identity " + label)
        check(transport.compose(point) == expected,
              "composition identity " + label)
        check(transport.forward(transport.dual(forward))
              == transport.codomain.mul(forward, DEGREE),
              "other dual identity " + label)
        point_rows[label] = {"forward": encoded(forward),
                             "positive_dual": encoded(expected)}
    forward_exceptions = kernel_exceptions(twist, transport.line_generator,
                                           transport.forward_twist)
    reverse_exceptions = kernel_exceptions(transport.forward_twist.codomain,
                                           transport.reverse_kernel_generator,
                                           transport.raw_reverse_twist)
    offcurve_kernel_rejection(transport.forward)
    offcurve_kernel_rejection(transport.raw_reverse)
    check((forward_exceptions, reverse_exceptions) == (262, 262),
          "incomplete twist-kernel exception enumeration")
    row = {
        "line": list(line),
        "generator": encoded(transport.line_generator),
        "complement": encoded(transport.line_complement),
        "forward_b": hex(transport.forward.codomain.b),
        "reverse_b": hex(transport.raw_reverse.codomain.b),
        "raw_reverse_scalar_sign": transport.raw_reverse_scalar_sign,
        "forward_kernel_hash": digest(sorted(transport.forward.kernel_abscissae)),
        "reverse_kernel_hash": digest(sorted(transport.raw_reverse.kernel_abscissae)),
        "point_controls": point_rows,
        "forward_twist_kernel_points": forward_exceptions,
        "reverse_twist_kernel_points": reverse_exceptions,
        "offcurve_kernel_abscissa_rejections": 2,
        "infinity_controls": 5,
    }
    within_budget(deadline)
    return (row, setup_seconds, time.perf_counter() - begun,
            tuple(sorted(transport.forward.kernel_abscissae)))


def reference_half_kernel(reference, generator, a, b):
    point, points, seen = generator, [], set()
    for _ in range((DEGREE - 1) // 2):
        check(point is not None and reference.on_curve(point, a, b),
              "independent twist kernel membership")
        check(point[0] not in seen, "independent duplicate kernel abscissa")
        points.append(point)
        seen.add(point[0])
        previous = point
        point = reference.add(point, generator, a)
    check(point == neg(previous), "independent twist kernel closure")
    return points


def reference_rational_velu(reference, point, kernel, a):
    x, y = point
    shift = 0
    for k in kernel:
        plus = reference.add(point, k, a)
        minus = reference.add(point, neg(k), a)
        check(plus is not None and minus is not None,
              "independent complementary image exception")
        x ^= plus[0] ^ minus[0]
        y ^= plus[1] ^ minus[1] ^ k[1] ^ neg(k)[1]
        shift ^= k[0]
    return x, y ^ shift


def reference_line(source, twist, basis, line, expected_row, deadline):
    within_budget(deadline)
    reference = Replay(Field(131, IRR131))
    u, v = basis
    check(reference.on_curve(u, 1) and reference.on_curve(v, 1),
          "independent twist basis membership")
    if line[0] == 1:
        generator = reference.add(u, reference.mul(v, line[1], 1), 1)
        complement = v
    else:
        generator, complement = v, u
    check(encoded(generator) == expected_row["generator"]
          and encoded(complement) == expected_row["complement"],
          "independent basis-line derivation")
    kernel = reference_half_kernel(reference, generator, 1, 1)
    reverse_generator = reference_rational_velu(reference, complement, kernel, 1)
    check(reference.on_curve(reverse_generator, 1,
                             int(expected_row["forward_b"], 16)),
          "independent reverse generator membership")
    reverse_kernel = reference_half_kernel(reference, reverse_generator, 1,
                                           int(expected_row["forward_b"], 16))
    check(digest(sorted(x for x, _ in kernel))
          == expected_row["forward_kernel_hash"],
          "independent forward kernel fingerprint")
    check(digest(sorted(x for x, _ in reverse_kernel))
          == expected_row["reverse_kernel_hash"],
          "independent reverse kernel fingerprint")
    results = {}
    for label, point in (("P", P), ("Q", Q)):
        expected_forward = reference.direct_velu(point, kernel)
        check(encoded(expected_forward)
              == expected_row["point_controls"][label]["forward"],
              "independent forward full coordinates " + label)
        raw_reverse = reference.direct_velu(expected_forward, reverse_kernel)
        expected_positive = reference.mul(point, DEGREE)
        actual_positive = (raw_reverse if expected_row["raw_reverse_scalar_sign"] == 1
                           else neg(raw_reverse))
        check(actual_positive == expected_positive,
              "independent dual full coordinates " + label)
        check(encoded(expected_positive)
              == expected_row["point_controls"][label]["positive_dual"],
              "independent positive dual receipt " + label)
        results[label] = {"forward": encoded(expected_forward),
                          "raw_reverse": encoded(raw_reverse)}
    within_budget(deadline)
    return {"line": list(line), "points": results,
            "reverse_generator": encoded(reverse_generator)}


def negative_controls(source, twist, basis):
    u, v = basis
    controls = []
    for line in ((0, 0), (2, 0), (1, 263), (-1, 0), (1, True), [1, 0]):
        controls.append(reject(lambda x=line: build_line_dual(source, twist, basis, x, P),
                               "noncanonical line " + repr(line)))
    controls.append(reject(lambda: build_line_dual(source, twist, (u, (0, 0)),
                                                   (1, 0), P), "off-curve basis"))
    controls.append(reject(lambda: build_line_dual(source, twist, (u, (0, 1)),
                                                   (1, 0), P), "wrong-order basis"))
    controls.append(reject(lambda: build_line_dual(source, twist, (u, u),
                                                   (1, 0), P), "dependent basis"))
    controls.append(reject(lambda: build_line_dual(source, twist, basis,
                                                   (1, 0), (0, 1)), "ambiguous witness"))
    return controls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("preserve earlier evidence: choose a fresh --out")
    begun = time.monotonic()
    deadline = begun + WALL_LIMIT_SECONDS
    result = {"schema": "ecc2k130_degree263_dual_all_lines_v1",
              "status": "RUNNING", "base_commit": "dc7fa83bc2f6061f6b62927a4adb6a564a83eeae",
              "host": platform.platform(), "python": sys.version,
              "command": [sys.executable, *sys.argv],
              "cost_scope": "saved torsion bases supplied; no cold discovery, PDP, rank or DLP",
              "PDP_cost": None, "full_ECDLP_cost": None, "end_to_end_speedup": None,
              "runs": []}
    stage = "lock"
    try:
        lock = json.loads((HERE / "LOCK.json").read_text())
        actual = {name: sha(REPO / name) for name in lock["sha256"]}
        check(actual == lock["sha256"], "source/input SHA-256 lock mismatch")
        result["source_input_sha256"] = actual
        frozen = json.loads((RESEARCH / "ecc2k130_direction_review_20260924"
                             / "twist_torsion_results.json").read_text())
        check([run["seed"] for run in frozen["runs"]] == [20260924, 20260925],
              "frozen basis seeds changed")
        all_kernel_sets = []
        for index, saved in enumerate(frozen["runs"]):
            stage = "all lines seed " + str(saved["seed"])
            field = FastGF2m(131, IRR131)
            source, twist = Koblitz(field, a=0), Koblitz(field, a=1)
            basis = basis_from_run(saved)
            check(source.on_curve(P) and source.on_curve(Q), "public point membership")
            seed_result = {"seed": saved["seed"], "basis": [encoded(x) for x in basis],
                           "rows": [], "setup_seconds": 0.0,
                           "all_checks_seconds": 0.0, "reference_rows": []}
            result["runs"].append(seed_result)
            seed_result["negative_controls"] = negative_controls(source, twist, basis)
            kernel_sets = set()
            tested_lines = LINES if index == 0 else REFERENCE_LINES
            for line in tested_lines:
                stage = "line " + str(saved["seed"]) + ":" + str(line)
                row, setup_s, total_s, roots = run_line(
                    source, twist, basis, line, deadline)
                seed_result["rows"].append(row)
                kernel_sets.add(roots)
                seed_result["setup_seconds"] += setup_s
                seed_result["all_checks_seconds"] += total_s
            check(len(seed_result["rows"]) == len(tested_lines),
                  "missing preregistered direction")
            check(len(kernel_sets) == len(tested_lines), "duplicate exact kernel line")
            all_kernel_sets.append(kernel_sets)
            for line in REFERENCE_LINES:
                stage = "independent reference " + str(saved["seed"]) + ":" + str(line)
                row = next(row for row in seed_result["rows"]
                           if row["line"] == list(line))
                seed_result["reference_rows"].append(
                    reference_line(source, twist, basis, line, row, deadline))
        check(len(all_kernel_sets[0]) == 264
              and all_kernel_sets[1].issubset(all_kernel_sets[0]),
              "second basis direction missing from full 263-kernel line set")
        result["distinct_line_fingerprints"] = len(all_kernel_sets[0])
        stable = [{key: run[key] for key in ("seed", "basis", "rows",
                   "reference_rows", "negative_controls")} for run in result["runs"]]
        result["deterministic_sha256"] = digest(stable)
        within_budget(deadline)
        result["status"] = "PASS_ALL_LINES"
    except (TimeoutError, MemoryError) as exc:
        result["status"] = "CENSORED"
        result["failure"] = {"stage": stage, "type": type(exc).__name__,
                             "message": str(exc)}
    except Exception as exc:
        result["status"] = "FAIL"
        result["failure"] = {"stage": stage, "type": type(exc).__name__,
                             "message": str(exc), "traceback": traceback.format_exc()}
    result["elapsed_seconds"] = time.monotonic() - begun
    result["peak_rss_bytes"] = peak_bytes()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "rows": [len(x["rows"])
                      for x in result["runs"]], "elapsed_seconds": result["elapsed_seconds"],
                      "out": str(args.out)}))
    return 0 if result["status"] == "PASS_ALL_LINES" else 1


if __name__ == "__main__":
    raise SystemExit(main())
