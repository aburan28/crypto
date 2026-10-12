#!/usr/bin/env python3
"""Frozen truth-table producer for generic and exact-leaf full-point edges."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

from native_relation import build_relation
from dag import PackedModel  # noqa: E402 - native_relation pins the frozen parent DAG

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
from fastfield import FastGF2m  # noqa: E402
from relations import Koblitz  # noqa: E402

class BitField:
    """Shift/reduce reference for the exhaustive small fields."""

    def __init__(self, n: int, modulus: int):
        self.deg, self.modulus = n, modulus

    def mul(self, a: int, b: int) -> int:
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & (1 << self.deg):
                a ^= self.modulus
        return result

    def sqr(self, a: int) -> int:
        return self.mul(a, a)

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError
        for candidate in range(1, 1 << self.deg):
            if self.mul(a, candidate) == 1:
                return candidate
        raise AssertionError("no field inverse")


FIELDS = ((2, 0x7), (3, 0xB))
EXACT_MODULUS = 0x800000000000000000000000000002007
O = (1, 0, 0)
DOMAIN = "ecc2k130-native-fullpoint-edge-20260929-v1"
EXACT = ROOT / "research/ecc2k130_factor_base_replication_20260925/exact_smoke.json"
EXACT_REPLAY = ROOT / "research/ecc2k130_factor_base_replication_20260925/exact_replay.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def raw(point):
    return O if point is None else (0, point[0], point[1])


def ref_points(curve, n: int):
    return [O] + [(0, x, y) for x in range(1 << n) for y in range(1 << n)
                  if curve.on_curve((x, y))]


def branch(p, q):
    if p[0]:
        return "copy_q"
    if q[0]:
        return "copy_p"
    if p[1] == q[1] and p[2] ^ q[2] == p[1]:
        return "inverse"
    if p[1] == q[1]:
        return "double"
    return "generic"


def exact_slope(curve, p, q):
    kind = branch(p, q)
    if kind in ("copy_q", "copy_p", "inverse"):
        return 0
    f = curve.F
    if kind == "double":
        return p[1] ^ f.mul(p[2], f.inv(p[1]))
    return f.mul(p[2] ^ q[2], f.inv(p[1] ^ q[1]))


def check_models(relation):
    model = relation.model(O, O, O, 0)
    checks = (
        PackedModel(model.bit_count - 1, model.limbs),
        PackedModel(model.bit_count, model.limbs[:-1]),
        PackedModel(model.bit_count, model.limbs[:-1] + (model.limbs[-1] | (1 << 63),)),
    )
    for bad in checks:
        try:
            relation.dag.evaluate(bad, relation.output)
        except ValueError:
            continue
        raise AssertionError("malformed model accepted")
    return len(checks)


def check_coefficients(n, modulus):
    for a, b in ((0, 0), (1 << n, 1), (0, 1 << n)):
        try:
            build_relation(n, modulus, a, b)
        except ValueError:
            continue
        raise AssertionError("invalid curve coefficient accepted")
    return 3


def toy_rows(write):
    summaries = []
    first_mismatch = None
    for n, modulus in FIELDS:
        coefficient_controls = check_coefficients(n, modulus)
        field = BitField(n, modulus)
        for a in (0, 1):
            for b in range(1, 1 << n):
                curve = Koblitz(field, a=a, b=b)
                relation = build_relation(n, modulus, a, b)
                points = ref_points(curve, n)
                universe = [(o, x, y) for o in (0, 1) for x in range(1 << n)
                            for y in range(1 << n)]
                invalid = [point for point in universe if point not in points]
                cases = dict.fromkeys(relation.branches, 0)
                rows = valid_evals = invalid_evals = mismatch_count = 0
                for p in points:
                    for q in points:
                        result = raw(curve.add(None if p[0] else p[1:],
                                               None if q[0] else q[1:]))
                        assert result in points
                        kind = branch(p, q)
                        cases[kind] += 1
                        expected_lambdas = (1 << n) if kind in ("copy_q", "copy_p", "inverse") else 1
                        for r in points:
                            accepted = sum(relation.accepts(p, q, r, lam)
                                           for lam in range(1 << n))
                            expected = expected_lambdas if r == result else 0
                            row = {"kind": "toy", "n": n, "modulus": modulus,
                                   "a": a, "b": b, "p": p, "q": q, "r": r,
                                   "branch": kind, "accepted": accepted, "expected": expected}
                            write(row)
                            rows += 1
                            valid_evals += 1 << n
                            if accepted != expected:
                                mismatch_count += 1
                                first_mismatch = first_mismatch or row
                for bad in invalid:
                    for role, (p, q, r) in enumerate(((bad, O, O), (O, bad, O), (O, O, bad))):
                        accepted = sum(relation.accepts(p, q, r, lam)
                                       for lam in range(1 << n))
                        row = {"kind": "invalid", "n": n, "modulus": modulus,
                               "a": a, "b": b, "role": role, "point": bad,
                               "accepted": accepted}
                        write(row)
                        invalid_evals += 1 << n
                        if accepted:
                            mismatch_count += 1
                            first_mismatch = first_mismatch or row
                summaries.append({"n": n, "modulus": modulus, "a": a, "b": b,
                                  "curve_points": len(points), "rows": rows,
                                  "valid_evaluations": valid_evals,
                                  "invalid_points": len(invalid),
                                  "invalid_evaluations": invalid_evals,
                                  "case_counts": cases, "dag": relation.dag.counts(),
                                  "coefficient_controls": coefficient_controls,
                                  "model_controls": check_models(relation),
                                  "mismatches": mismatch_count})
    if any(sum(item["case_counts"][case] for item in summaries) == 0
           for case in ("copy_q", "copy_p", "inverse", "double", "generic")):
        first_mismatch = first_mismatch or {"missing_branch": True}
    return summaries, first_mismatch


def leaf_rows(write):
    archived = json.loads(EXACT.read_text())
    replay = json.loads(EXACT_REPLAY.read_text())
    assert archived["status"] == replay["status"] == "PASS"
    assert replay["input_sha256"] == sha(EXACT)
    field = FastGF2m(131, EXACT_MODULUS)
    summaries = []
    first_mismatch = None
    for record in archived["lines"]:
        line = tuple(record["line"])
        assert line in ((1, 0), (1, 4))
        b = int(record["codomain_b"], 16)
        curve = Koblitz(field, a=0, b=b)
        relation = build_relation(131, EXACT_MODULUS, 0, b)
        base = [tuple(int(value, 16) for value in point)
                for point in record["bases"]["descendant_native"]]
        assert len(base) == 16 and len(set(base)) == 16
        assert all(curve.on_curve(point) for point in base)
        p, q = base[0], base[2]
        assert p[0] != q[0]
        zero_point = curve.points_over(0)[0]
        assert curve.on_curve(zero_point) and curve.add(zero_point, zero_point) is None
        panel = (
            ("copy_q", None, p), ("copy_p", p, None),
            ("inverse", p, curve.neg(p)), ("double", p, p),
            ("generic", p, q), ("zero_self_inverse", zero_point, zero_point),
        )
        counts = dict.fromkeys(("copy_q", "copy_p", "inverse", "double", "generic"), 0)
        rows = mismatches = 0
        for label, left, right in panel:
            left_raw, right_raw = raw(left), raw(right)
            kind = branch(left_raw, right_raw)
            assert kind == ("inverse" if label == "zero_self_inverse" else label)
            counts[kind] += 1
            result = curve.add(left, right)
            result_raw = raw(result)
            slope = exact_slope(curve, left_raw, right_raw)
            wrong_result = next(point for point in base if point != result)
            checks = [("correct", result_raw, slope, True),
                      ("wrong_result", raw(wrong_result), slope, False)]
            if kind in ("copy_q", "copy_p", "inverse"):
                checks.append(("free_slope", result_raw, 1, True))
            else:
                checks.append(("wrong_slope", result_raw, slope ^ 1, False))
            for control, r, lam, expected in checks:
                accepted = relation.accepts(left_raw, right_raw, r, lam)
                row = {"kind": "leaf", "line": line, "b": hex(b),
                       "case": label, "control": control, "p": left_raw,
                       "q": right_raw, "r": r, "lambda": hex(lam),
                       "accepted": accepted, "expected": expected}
                write(row)
                rows += 1
                if accepted != expected:
                    mismatches += 1
                    first_mismatch = first_mismatch or row
        offcurve = None
        for bit in range(131):
            candidate = (p[0], p[1] ^ (1 << bit))
            if not curve.on_curve(candidate):
                offcurve = candidate
                break
        assert offcurve is not None
        for control, bad in (("off_curve", raw(offcurve)),
                             ("noncanonical_infinity", (1, 1, 0))):
            accepted = relation.accepts(bad, O, O, 0)
            row = {"kind": "leaf_invalid", "line": line, "b": hex(b),
                   "control": control, "point": bad,
                   "accepted": accepted, "expected": False}
            write(row)
            rows += 1
            if accepted:
                mismatches += 1
                first_mismatch = first_mismatch or row
        summaries.append({"line": line, "b": hex(b), "rows": rows,
                          "case_counts": counts, "dag": relation.dag.counts(),
                          "model_controls": check_models(relation),
                          "mismatches": mismatches})
    assert len(summaries) == 2
    return summaries, first_mismatch


def run(out: Path) -> dict:
    started, cpu_started = time.monotonic(), time.process_time()
    out.mkdir(parents=True, exist_ok=False)
    rows_path = out / "rows.jsonl.gz"
    row_count = 0
    with rows_path.open("wb") as raw_stream:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_stream, mtime=0) as gz:
            def write(row):
                nonlocal row_count
                gz.write((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode())
                row_count += 1
            toy, toy_mismatch = toy_rows(write)
            leaves, leaf_mismatch = leaf_rows(write)
    first_mismatch = toy_mismatch or leaf_mismatch
    result = {"domain": DOMAIN, "decision": "PASS" if first_mismatch is None else "FAIL",
              "toy": toy, "leaves": leaves, "raw_rows": row_count,
              "first_mismatch": first_mismatch,
              "rows_sha256": sha(rows_path), "exact_smoke_sha256": sha(EXACT),
              "exact_replay_sha256": sha(EXACT_REPLAY),
              "wall_seconds": time.monotonic() - started,
              "cpu_seconds": time.process_time() - cpu_started,
              "peak_rss_bytes": rss_bytes()}
    (out / "result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.out), sort_keys=True))
