#!/usr/bin/env python3
"""Fresh-process replay of the frozen all-line dual certificate.

All 270 line identities and exceptional inputs are recomputed. Twelve fixed
lines also use the separate bit-polynomial reference group law and direct
full-coordinate Vélu sums instead of trusting the production point map.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
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

from line_dual import build_line_dual
from fastfield import FastGF2m, IRR131
from redteam_velu_replay import Field, Replay
from relations import (CHALLENGE_PX, CHALLENGE_PY, CHALLENGE_QX,
                       CHALLENGE_QY, Koblitz)

P, Q = (CHALLENGE_PX, CHALLENGE_PY), (CHALLENGE_QX, CHALLENGE_QY)
ALL = [(1, t) for t in range(263)] + [(0, 1)]
REF = [(1, t) for t in (0, 1, 2, 131, 262)] + [(0, 1)]
WALL_LIMIT_SECONDS = 2100
MEMORY_LIMIT_BYTES = 2 * 1024**3


def require(condition, why):
    if not condition:
        raise AssertionError(why)


def hash_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hash_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def point_json(point):
    return None if point is None else [hex(point[0]), hex(point[1])]


def check_budget(deadline):
    if time.monotonic() >= deadline:
        raise TimeoutError("2,100-second independent replay wall limit")
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if (rss if platform.system() == "Darwin" else rss * 1024) > MEMORY_LIMIT_BYTES:
        raise MemoryError("2 GiB independent replay peak-memory limit")


def negative(point):
    return None if point is None else (point[0], point[0] ^ point[1])


def replay_kernel(curve, generator, velu):
    """Check the finite exceptional path without using certificate helpers."""
    point = generator
    seen = set()
    for _ in range(131):
        require(point is not None and curve.on_curve(point), "replay kernel membership")
        require(point[0] not in seen, "replay duplicate kernel root")
        seen.add(point[0])
        require(velu(point) is None, "positive kernel exception")
        require(velu(curve.neg(point)) is None, "negative kernel exception")
        previous = point
        point = curve.add(point, generator)
    require(point == curve.neg(previous), "replay kernel closure")
    require(seen == velu.kernel_abscissae, "replay kernel root set")
    return 262


def replay_bad_abscissa(velu):
    root = min(velu.kernel_abscissae)
    y = 0
    while velu.source.on_curve((root, y)):
        y += 1
    try:
        velu((root, y))
    except ValueError:
        return 1
    raise AssertionError("off-curve kernel-abscissa input accepted")


def reference_half(reference, generator, a, b):
    half = []
    point = generator
    for _ in range(131):
        require(point is not None and reference.on_curve(point, a, b),
                "reference kernel membership")
        half.append(point)
        previous = point
        point = reference.add(point, generator, a)
    require(point == negative(previous), "reference kernel closure")
    require(len({p[0] for p in half}) == 131, "reference distinct roots")
    return half


def reference_complement_image(reference, complement, kernel):
    x, y = complement
    total = 0
    for point in kernel:
        plus = reference.add(complement, point, 1)
        minus = reference.add(complement, negative(point), 1)
        require(plus is not None and minus is not None,
                "reference complement in kernel")
        x ^= plus[0] ^ minus[0]
        y ^= plus[1] ^ minus[1] ^ point[0]
        total ^= point[0]
    return x, y ^ total


def replay_reference(basis, row, saved):
    reference = Replay(Field(131, IRR131))
    line = tuple(row["line"])
    u, v = basis
    if line[0]:
        generator = reference.add(u, reference.mul(v, line[1], 1), 1)
        complement = v
    else:
        generator, complement = v, u
    require(point_json(generator) == row["generator"],
            "reference line generator")
    forward_half = reference_half(reference, generator, 1, 1)
    reverse_generator = reference_complement_image(reference, complement,
                                                   forward_half)
    require(point_json(reverse_generator) == saved["reverse_generator"],
            "reference reverse generator")
    reverse_half = reference_half(reference, reverse_generator, 1,
                                  int(row["forward_b"], 16))
    require(hash_json(sorted(p[0] for p in forward_half))
            == row["forward_kernel_hash"], "reference forward roots")
    require(hash_json(sorted(p[0] for p in reverse_half))
            == row["reverse_kernel_hash"], "reference reverse roots")
    for label, point in (("P", P), ("Q", Q)):
        forward = reference.direct_velu(point, forward_half)
        raw = reference.direct_velu(forward, reverse_half)
        corrected = raw if row["raw_reverse_scalar_sign"] == 1 else negative(raw)
        require(corrected == reference.mul(point, 263),
                "reference dual scalar identity " + label)
        require(saved["points"][label]
                == {"forward": point_json(forward), "raw_reverse": point_json(raw)},
                "reference point receipt " + label)
        require(row["point_controls"][label]["forward"] == point_json(forward),
                "reference forward receipt " + label)


def replay_negative_controls(source, twist, basis, expected):
    u, _ = basis
    cases = [
        ((0, 0), basis, P), ((2, 0), basis, P), ((1, 263), basis, P),
        ((-1, 0), basis, P), ((1, True), basis, P), ([1, 0], basis, P),
        ((1, 0), (u, (0, 0)), P), ((1, 0), (u, (0, 1)), P),
        ((1, 0), (u, u), P), ((1, 0), basis, (0, 1)),
    ]
    require(len(expected) == len(cases), "negative-control receipt count")
    for line, checked_basis, witness in cases:
        try:
            build_line_dual(source, twist, checked_basis, line, witness)
        except (ValueError, TypeError):
            continue
        raise AssertionError("negative control accepted: " + repr(line))
    return len(cases)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("preserve earlier evidence: choose a fresh --out")
    begun = time.monotonic()
    deadline = begun + WALL_LIMIT_SECONDS
    receipt = {"schema": "ecc2k130_degree263_dual_all_lines_replay_v1",
               "status": "RUNNING", "result_sha256": hash_file(args.result),
               "checked_lines": 0, "checked_kernel_inputs": 0,
               "checked_reference_lines": 0, "checked_negative_controls": 0}
    stage = "lock"
    try:
        lock = json.loads((HERE / "LOCK.json").read_text())
        require({name: hash_file(REPO / name) for name in lock["sha256"]}
                == lock["sha256"], "source/input lock mismatch")
        result = json.loads(args.result.read_text())
        require(result["status"] == "PASS_ALL_LINES", "producer did not pass")
        stable = [{key: run[key] for key in ("seed", "basis", "rows",
                   "reference_rows", "negative_controls")} for run in result["runs"]]
        require(hash_json(stable) == result["deterministic_sha256"],
                "deterministic result digest mismatch")
        source_data = json.loads((RESEARCH / "ecc2k130_direction_review_20260924"
                                  / "twist_torsion_results.json").read_text())
        require([r["seed"] for r in result["runs"]] == [20260924, 20260925],
                "result seed order")
        first_roots = set()
        for index, (saved, run) in enumerate(zip(source_data["runs"], result["runs"])):
            stage = "seed " + str(run["seed"])
            expected_lines = ALL if index == 0 else REF
            require([row["line"] for row in run["rows"]]
                    == [list(line) for line in expected_lines],
                    "line order or coverage")
            basis = tuple(tuple(int(z, 16) for z in point)
                          for point in saved["basis_on_twist"])
            require(run["basis"] == [point_json(x) for x in basis],
                    "basis receipt")
            field = FastGF2m(131, IRR131)
            source, twist = Koblitz(field, a=0), Koblitz(field, a=1)
            receipt["checked_negative_controls"] += replay_negative_controls(
                source, twist, basis, run["negative_controls"])
            by_line = {tuple(row["line"]): row for row in run["rows"]}
            roots = set()
            for line in expected_lines:
                check_budget(deadline)
                stage = "line " + str(run["seed"]) + ":" + str(line)
                row = by_line[line]
                transport = build_line_dual(source, twist, basis, line, P)
                require(point_json(transport.line_generator) == row["generator"]
                        and point_json(transport.line_complement) == row["complement"],
                        "line point receipt")
                require(hex(transport.codomain.b) == row["forward_b"]
                        and hex(transport.raw_reverse.codomain.b) == row["reverse_b"]
                        and transport.raw_reverse_scalar_sign == row["raw_reverse_scalar_sign"],
                        "quotient model or orientation")
                fingerprint = hash_json(sorted(transport.forward.kernel_abscissae))
                require(fingerprint == row["forward_kernel_hash"],
                        "forward kernel fingerprint")
                require(hash_json(sorted(transport.raw_reverse.kernel_abscissae))
                        == row["reverse_kernel_hash"], "reverse kernel fingerprint")
                roots.add(tuple(sorted(transport.forward.kernel_abscissae)))
                for label, point in (("P", P), ("Q", Q)):
                    forward = transport.forward(point)
                    expected = source.mul(point, 263)
                    require(transport.dual(forward) == expected
                            and transport.forward(transport.dual(forward))
                            == transport.codomain.mul(forward, 263),
                            "dual identities " + label)
                    require(row["point_controls"][label]
                            == {"forward": point_json(forward),
                                "positive_dual": point_json(expected)},
                            "point receipt " + label)
                require(transport.forward(None) is None
                        and transport.forward_twist(None) is None
                        and transport.raw_reverse_twist(None) is None
                        and transport.dual(None) is None
                        and transport.compose(None) is None,
                        "infinity exceptions")
                count = replay_kernel(twist, transport.line_generator,
                                      transport.forward_twist)
                count += replay_kernel(transport.forward_twist.codomain,
                                       transport.reverse_kernel_generator,
                                       transport.raw_reverse_twist)
                require(count == row["forward_twist_kernel_points"]
                        + row["reverse_twist_kernel_points"] == 524,
                        "kernel exception counts")
                require(replay_bad_abscissa(transport.forward)
                        + replay_bad_abscissa(transport.raw_reverse)
                        == row["offcurve_kernel_abscissa_rejections"] == 2,
                        "invalid kernel-abscissa receipt")
                require(row["infinity_controls"] == 5, "infinity receipt")
                receipt["checked_lines"] += 1
                receipt["checked_kernel_inputs"] += count
            require(len(roots) == len(expected_lines), "duplicate exact kernel line")
            if index == 0:
                first_roots = roots
            else:
                require(roots.issubset(first_roots), "second-basis line not in first set")
            refs = {tuple(x["line"]): x for x in run["reference_rows"]}
            require(set(refs) == set(REF), "missing independent reference line")
            for line in REF:
                check_budget(deadline)
                stage = "independent reference " + str(run["seed"]) + ":" + str(line)
                replay_reference(basis, by_line[line], refs[line])
                receipt["checked_reference_lines"] += 1
        require(result["distinct_line_fingerprints"] == len(first_roots) == 264,
                "distinct first-basis line count")
        require((receipt["checked_lines"], receipt["checked_kernel_inputs"],
                 receipt["checked_reference_lines"], receipt["checked_negative_controls"])
                == (270, 270 * 524, 12, 20), "final replay coverage")
        receipt["deterministic_sha256"] = result["deterministic_sha256"]
        receipt["status"] = "PASS"
    except Exception as exc:
        receipt["status"] = ("CENSORED" if isinstance(exc, (TimeoutError, MemoryError))
                             else "FAIL")
        receipt["failure"] = {"stage": stage, "type": type(exc).__name__,
                              "message": str(exc), "traceback": traceback.format_exc()}
    receipt["elapsed_seconds"] = time.monotonic() - begun
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "checked_lines": receipt["checked_lines"],
                      "out": str(args.out)}))
    return 0 if receipt["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
