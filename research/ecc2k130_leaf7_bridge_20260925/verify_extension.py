#!/usr/bin/env python3
"""Separate Sage replay of archived full exceptional points over Fq^6.

This script does not import run.py, verify.py, or their curve-map helpers. It
reconstructs the frozen field tower, the 7-isogeny, and its dual from the
archive, then checks all twelve saved extension points. It is never run in CI.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def start_clock():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    peak = int(usage.ru_maxrss)
    if sys.platform.startswith("linux"):
        peak *= 1024
    return time.monotonic(), time.process_time(), peak


def timeout_handler(_signum, _frame):
    raise TimeoutError("frozen 300-second extension replay wall cap reached")


def base_field(spec):
    from sage.all import GF, PolynomialRing

    b = GF(2)
    R = PolynomialRing(b, "z0")
    bits = int(spec["field"]["modulus_hex"], 16)
    modulus = R([b((bits >> i) & 1) for i in range(bits.bit_length())])
    assert bits.bit_length() - 1 == spec["field"]["degree"]
    F = GF(2 ** spec["field"]["degree"], name="z", modulus=modulus)
    assert F.modulus() == modulus
    return F


def to_fq(F, bits: int):
    assert 0 <= bits < (1 << F.degree())
    result = F.zero()
    z = F.gen()
    for i in range(F.degree() - 1, -1, -1):
        result = result * z + F((bits >> i) & 1)
    return result


def from_fq(value):
    return sum(int(coef) << i for i, coef in enumerate(value.polynomial().list()))


def tower_from_polynomial(F, h):
    from sage.all import PolynomialRing

    F3 = F.extension(h, "u")
    R3 = PolynomialRing(F3, "w")
    w = R3.gen()
    F6 = F3.extension(w ** 2 + w + 1, "v")
    assert F3.degree() == 3 and F6.degree() == 2
    return F3, F6


def encode3(value):
    coeffs = value.lift().list()
    assert len(coeffs) <= 3
    return [hex(from_fq(coeffs[i])) if i < len(coeffs) else "0x0"
            for i in range(3)]


def encode6(value):
    coeffs = value.lift().list()
    assert len(coeffs) <= 2
    zero = value.parent().base_ring().zero()
    return [encode3(coeffs[i] if i < len(coeffs) else zero)
            for i in range(2)]


def decode3(F, F3, saved):
    assert isinstance(saved, list) and len(saved) == 3
    result = F3.zero()
    U = F3.gen()
    for i, label in enumerate(saved):
        result += F3(to_fq(F, int(label, 16))) * U ** i
    assert encode3(result) == saved
    return result


def decode6(F, F3, F6, saved):
    assert isinstance(saved, list) and len(saved) == 2
    result = F6.zero()
    V = F6.gen()
    for i, row in enumerate(saved):
        result += F6(decode3(F, F3, row)) * V ** i
    assert encode6(result) == saved
    return result


def canonical(record):
    return json.dumps(record, sort_keys=True)


def replay_kernel(F, F3, F6, curve, isogeny, polynomial, archived, q):
    assert int(polynomial.degree()) == 3
    expected_coeffs = [encode6(polynomial[i]) for i in range(4)]
    assert archived["kernel_polynomial_low_to_high"] == expected_coeffs
    roots = polynomial.roots()
    assert len(roots) == 3 and all(int(multiplicity) == 1 for _, multiplicity in roots)
    actual_roots = sorted([encode6(x) for x, _ in roots], key=canonical)
    assert actual_roots == archived["roots"]
    assert archived["root_count"] == 3 and archived["point_count"] == 6
    assert len(archived["points"]) == 6
    recorded = set()
    counts = {}
    for saved in archived["points"]:
        key = canonical(saved)
        assert key not in recorded
        recorded.add(key)
        x = decode6(F, F3, F6, saved["x"])
        y = decode6(F, F3, F6, saved["y"])
        assert polynomial(x) == 0
        point = curve([x, y])
        assert not point.is_zero() and (7 * point).is_zero()
        five = 5 * point
        assert (point[0] ** q, point[1] ** q) == (five[0], five[1])
        assert isogeny(point).is_zero()
        xkey = canonical(saved["x"])
        counts[xkey] = counts.get(xkey, 0) + 1
    assert len(counts) == 3 and all(count == 2 for count in counts.values())
    # Recompute both lifts at every archived root, independently of their
    # order in the producer's list.
    actual_points = set()
    for x, _ in roots:
        lifts = curve.lift_x(x, all=True)
        assert len(lifts) == 2
        for point in lifts:
            actual_points.add(canonical({"x": encode6(point[0]),
                                         "y": encode6(point[1])}))
    assert recorded == actual_points
    return {"roots": 3, "points": 6, "order_seven": True,
            "q_frobenius_eigenvalue": 5, "mapped_to_infinity": True,
            "matches_all_full_lifts": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--fq-replay", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite extension replay receipt")
    started = start_clock()
    phase = "initialization"
    result = {"schema": "ecc2k130-leaf7-bridge-extension-replay-v1",
              "status": "STOP", "structural_gate_status": "UNRESOLVED",
              "phase": phase}
    try:
        spec = json.loads((HERE / "FROZEN.json").read_text())
        source = json.loads(args.input.read_text())
        fq = json.loads(args.fq_replay.read_text())
        result["freeze_sha256"] = sha(HERE / "FROZEN.json")
        result["source_receipt_sha256"] = sha(args.input)
        result["fq_replay_sha256"] = sha(args.fq_replay)
        if spec["release_main_head"] is None:
            raise RuntimeError("parent note unmerged: post-merge re-freeze and review required")
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.setitimer(signal.ITIMER_REAL, spec["caps"]["child_wall_seconds"])
        try:
            resource.setrlimit(resource.RLIMIT_AS,
                               (spec["caps"]["child_peak_rss_bytes"],
                                spec["caps"]["child_peak_rss_bytes"]))
        except (AttributeError, OSError, ValueError) as exc:
            raise RuntimeError("frozen 2-GiB memory cap could not be enforced") from exc
        phase = "frozen_inputs"
        assert source["status"] == "PRODUCER_PASS"
        assert fq["status"] == "FQ_REPLAY_PASS"
        assert source["structural_gate_status"] == "UNRESOLVED_PENDING_REPLAY"
        assert fq["structural_gate_status"] == "UNRESOLVED_PENDING_EXTENSION_REPLAY"
        assert source["freeze_sha256"] == result["freeze_sha256"]
        assert fq["freeze_sha256"] == result["freeze_sha256"]
        assert fq["source_receipt_sha256"] == result["source_receipt_sha256"]
        assert source["release_main_head"] == spec["release_main_head"]
        assert source["input_sha256"] == spec["input_sha256"]
        assert source["implementation_sha256"] == spec["implementation_sha256"]
        for relative, expected in spec["input_sha256"].items():
            assert sha(REPO / relative) == expected, relative
        for name, expected in spec["implementation_sha256"].items():
            assert sha(HERE / name) == expected, name

        phase = "reconstruct_cubic_and_division_factors"
        from sage.all import EllipticCurve, PolynomialRing
        from sage.version import version as sage_version
        assert str(sage_version).startswith(spec["sage_version"])
        F = base_field(spec)
        R = PolynomialRing(F, "x")
        saved_h = source["kernel_polynomial_low_to_high"]
        assert len(saved_h) == 4 and saved_h[3] == "0x1"
        h = R([to_fq(F, int(label, 16)) for label in saved_h])
        assert h.is_monic() and h.is_irreducible() and int(h.degree()) == 3
        E0 = EllipticCurve(F, [F.one(), F.zero(), F.zero(), F.zero(),
                               to_fq(F, int(source["leaf0_b"], 16))])
        factors = [(poly, int(power)) for poly, power in E0.division_polynomial(7).factor()]
        degrees = sorted([[int(poly.degree()), power] for poly, power in factors])
        assert degrees == spec["gate"]["predicted_7_division_factor_degrees"]
        assert any(poly.monic() == h and power == 1 for poly, power in factors)
        psi = E0.isogeny(h)
        assert int(psi.degree()) == 7 and psi.is_separable()
        assert [hex(from_fq(a)) for a in psi.codomain().a_invariants()] == source["sage_codomain_ainvs"]

        phase = "reconstruct_extension_maps"
        F3, F6 = tower_from_polynomial(F, h)
        E0e = E0.base_extend(F6)
        psie = E0e.isogeny(h.change_ring(F6))
        duale = psie.dual()
        q = 2 ** spec["field"]["degree"]
        result["forward"] = replay_kernel(F, F3, F6, E0e, psie,
                                          h.change_ring(F6), source["forward_kernel"], q)
        result["dual"] = replay_kernel(F, F3, F6, psie.codomain(), duale,
                                       duale.kernel_polynomial(), source["reverse_kernel"], q)
        assert psie(E0e(0)).is_zero() and duale(psie.codomain()(0)).is_zero()
        assert source["cost"]["incremental_mul_equivalent_ratio"] is None
        assert source["cost"]["cold_mul_equivalent_ratio"] is None
        assert source["cost"]["promotion"] is False
        result["cost_promotion"] = False
        result["status"] = "PASS"
        result["structural_gate_status"] = "PASS"
        phase = "complete"
    except BaseException as exc:
        result["status"] = "STOP"
        result["structural_gate_status"] = "UNRESOLVED"
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
        result["traceback"] = traceback.format_exc(limit=12)
    finally:
        end = start_clock()
        result["phase"] = phase
        result["resources"] = {"wall_seconds": end[0] - started[0],
                               "cpu_seconds": end[1] - started[1],
                               "peak_rss_bytes": end[2]}
        if "spec" in locals():
            if result["resources"]["wall_seconds"] > spec["caps"]["child_wall_seconds"]:
                result["status"] = "STOP"
                result["structural_gate_status"] = "UNRESOLVED"
                result.setdefault("error", "wall time exceeded frozen cap")
            if end[2] > spec["caps"]["child_peak_rss_bytes"]:
                result["status"] = "STOP"
                result["structural_gate_status"] = "UNRESOLVED"
                result.setdefault("error", "peak RSS exceeded frozen cap")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        signal.setitimer(signal.ITIMER_REAL, 0)
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
