#!/usr/bin/env python3
"""Independent long-division group-law replay of every full-point edge row."""
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
from dag import PackedModel  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
EXACT = ROOT / "research/ecc2k130_factor_base_replication_20260925/exact_smoke.json"
EXACT_REPLAY = ROOT / "research/ecc2k130_factor_base_replication_20260925/exact_replay.json"
FIELDS = ((2, 0x7), (3, 0xB))
EXACT_MODULUS = 0x800000000000000000000000000002007
DOMAIN = "ecc2k130-native-fullpoint-edge-20260929-v1"
O = (1, 0, 0)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


class PolynomialCurve:
    """Product-then-long-division field arithmetic and separate point law."""

    def __init__(self, n: int, modulus: int, a: int, b: int):
        self.n, self.modulus, self.a, self.b = n, modulus, a, b

    def mul(self, x: int, y: int) -> int:
        product = 0
        for i in range(self.n):
            for j in range(self.n):
                if (x >> i & 1) and (y >> j & 1):
                    product ^= 1 << (i + j)
        while product.bit_length() > self.n:
            product ^= self.modulus << (product.bit_length() - self.n - 1)
        return product

    def inv(self, value: int) -> int:
        if value == 0:
            raise ZeroDivisionError
        left, right = self.modulus, value
        u, v = 0, 1
        while right:
            shift = left.bit_length() - right.bit_length()
            if shift < 0:
                left, right = right, left
                u, v = v, u
                shift = -shift
            left ^= right << shift
            u ^= v << shift
        if left != 1:
            raise AssertionError("noninvertible field element")
        while u.bit_length() > self.n:
            u ^= self.modulus << (u.bit_length() - self.n - 1)
        assert self.mul(value, u) == 1
        return u

    def affine(self, x: int, y: int) -> bool:
        xx = self.mul(x, x)
        return (self.mul(y, y) ^ self.mul(x, y)
                == self.mul(xx, x) ^ self.mul(self.a, xx) ^ self.b)

    def points(self):
        return [O] + [(0, x, y) for x in range(1 << self.n)
                      for y in range(1 << self.n) if self.affine(x, y)]

    def neg(self, point):
        return O if point == O else (0, point[1], point[2] ^ point[1])

    def add(self, p, q):
        if p == O:
            return q
        if q == O:
            return p
        _, x, y = p
        _, u, v = q
        if x == u:
            if y ^ v == x:
                return O
            assert p == q and x != 0
            slope = x ^ self.mul(y, self.inv(x))
            result_x = self.mul(slope, slope) ^ slope ^ self.a
            result_y = self.mul(x, x) ^ self.mul(slope ^ 1, result_x)
        else:
            slope = self.mul(y ^ v, self.inv(x ^ u))
            result_x = self.mul(slope, slope) ^ slope ^ x ^ u ^ self.a
            result_y = self.mul(slope, x ^ result_x) ^ result_x ^ y
        result = (0, result_x, result_y)
        assert self.affine(result_x, result_y)
        return result

    def sqrt(self, value: int) -> int:
        result = value
        for _ in range(self.n - 1):
            result = self.mul(result, result)
        assert self.mul(result, result) == value
        return result

    def slope(self, p, q):
        kind = branch(p, q)
        if kind in ("copy_q", "copy_p", "inverse"):
            return 0
        if kind == "double":
            return p[1] ^ self.mul(p[2], self.inv(p[1]))
        return self.mul(p[2] ^ q[2], self.inv(p[1] ^ q[1]))


def branch(p, q):
    if p == O:
        return "copy_q"
    if q == O:
        return "copy_p"
    if p[1] == q[1] and p[2] ^ q[2] == p[1]:
        return "inverse"
    if p[1] == q[1]:
        return "double"
    return "generic"


def model_controls(relation):
    model = relation.model(O, O, O, 0)
    broken = (
        PackedModel(model.bit_count - 1, model.limbs),
        PackedModel(model.bit_count, model.limbs[:-1]),
        PackedModel(model.bit_count, model.limbs[:-1] + (model.limbs[-1] | (1 << 63),)),
    )
    for case in broken:
        try:
            relation.dag.evaluate(case, relation.output)
        except ValueError:
            pass
        else:
            raise AssertionError("bad packed model accepted")
    return len(broken)


def coefficient_controls(n, modulus):
    for a, b in ((0, 0), (1 << n, 1), (0, 1 << n)):
        try:
            build_relation(n, modulus, a, b)
        except ValueError:
            continue
        raise AssertionError("invalid curve coefficient accepted")
    return 3


def next_row(rows):
    line = rows.readline()
    assert line, "truncated producer row archive"
    return json.loads(line)


def replay(producer: Path, out: Path) -> dict:
    started, cpu_started = time.monotonic(), time.process_time()
    summary = json.loads((producer / "result.json").read_text())
    assert summary["domain"] == DOMAIN and summary["decision"] == "PASS"
    assert summary["first_mismatch"] is None
    rows_file = producer / "rows.jsonl.gz"
    assert sha(rows_file) == summary["rows_sha256"]
    assert sha(EXACT) == summary["exact_smoke_sha256"]
    assert sha(EXACT_REPLAY) == summary["exact_replay_sha256"]
    toy_counts, leaf_counts, total_rows = [], [], 0
    with gzip.open(rows_file, "rt") as rows:
        for n, modulus in FIELDS:
            invalid_coefficients = coefficient_controls(n, modulus)
            for a in (0, 1):
                for b in range(1, 1 << n):
                    curve = PolynomialCurve(n, modulus, a, b)
                    assert all(curve.mul(z, curve.inv(z)) == 1
                               for z in range(1, 1 << n))
                    relation = build_relation(n, modulus, a, b)
                    points = curve.points()
                    universe = [(o, x, y) for o in (0, 1) for x in range(1 << n)
                                for y in range(1 << n)]
                    invalid = [point for point in universe if point not in points]
                    cases = dict.fromkeys(relation.branches, 0)
                    count = valid_evals = invalid_evals = 0
                    for p in points:
                        for q in points:
                            result = curve.add(p, q)
                            assert result in points
                            kind = branch(p, q)
                            cases[kind] += 1
                            multiplier = (1 << n) if kind in ("copy_q", "copy_p", "inverse") else 1
                            for r in points:
                                accepted = sum(relation.accepts(p, q, r, lam)
                                               for lam in range(1 << n))
                                expected = multiplier if r == result else 0
                                row = next_row(rows)
                                assert row == {"kind": "toy", "n": n, "modulus": modulus,
                                               "a": a, "b": b, "p": list(p),
                                               "q": list(q), "r": list(r),
                                               "branch": kind, "accepted": accepted,
                                               "expected": expected}
                                assert accepted == expected
                                count += 1
                                valid_evals += 1 << n
                    for bad in invalid:
                        for role, (p, q, r) in enumerate(((bad, O, O), (O, bad, O), (O, O, bad))):
                            accepted = sum(relation.accepts(p, q, r, lam)
                                           for lam in range(1 << n))
                            row = next_row(rows)
                            assert row == {"kind": "invalid", "n": n, "modulus": modulus,
                                           "a": a, "b": b, "role": role,
                                           "point": list(bad), "accepted": accepted}
                            assert accepted == 0
                            invalid_evals += 1 << n
                    actual = {"n": n, "modulus": modulus, "a": a, "b": b,
                              "curve_points": len(points), "rows": count,
                              "valid_evaluations": valid_evals,
                              "invalid_points": len(invalid),
                              "invalid_evaluations": invalid_evals,
                              "case_counts": cases, "dag": relation.dag.counts(),
                              "coefficient_controls": invalid_coefficients,
                              "model_controls": model_controls(relation),
                              "mismatches": 0}
                    assert summary["toy"][len(toy_counts)] == actual
                    toy_counts.append(actual)
                    total_rows += count + len(invalid) * 3
        assert all(sum(item["case_counts"][case] for item in toy_counts) > 0
                   for case in ("copy_q", "copy_p", "inverse", "double", "generic"))

        archived = json.loads(EXACT.read_text())
        independent = json.loads(EXACT_REPLAY.read_text())
        assert archived["status"] == independent["status"] == "PASS"
        assert independent["input_sha256"] == sha(EXACT)
        for item in archived["lines"]:
            line = tuple(item["line"])
            assert line in ((1, 0), (1, 4))
            b = int(item["codomain_b"], 16)
            curve = PolynomialCurve(131, EXACT_MODULUS, 0, b)
            relation = build_relation(131, EXACT_MODULUS, 0, b)
            base = [(0, *(int(z, 16) for z in pair))
                    for pair in item["bases"]["descendant_native"]]
            assert len(base) == len(set(base)) == 16
            assert all(curve.affine(point[1], point[2]) for point in base)
            p, q = base[0], base[2]
            assert p[1] != q[1]
            t = (0, 0, curve.sqrt(b))
            assert curve.affine(t[1], t[2]) and curve.add(t, t) == O
            cases = (("copy_q", O, p), ("copy_p", p, O),
                     ("inverse", p, curve.neg(p)), ("double", p, p),
                     ("generic", p, q), ("zero_self_inverse", t, t))
            counts = dict.fromkeys(("copy_q", "copy_p", "inverse", "double", "generic"), 0)
            count = 0
            for label, left, right in cases:
                kind = branch(left, right)
                assert kind == ("inverse" if label == "zero_self_inverse" else label)
                counts[kind] += 1
                true_result = curve.add(left, right)
                lam = curve.slope(left, right)
                wrong_result = next(point for point in base if point != true_result)
                checks = [("correct", true_result, lam, True),
                          ("wrong_result", wrong_result, lam, False)]
                if kind in ("copy_q", "copy_p", "inverse"):
                    checks.append(("free_slope", true_result, 1, True))
                else:
                    checks.append(("wrong_slope", true_result, lam ^ 1, False))
                for control, result, slope, expected in checks:
                    accepted = relation.accepts(left, right, result, slope)
                    row = next_row(rows)
                    assert row == {"kind": "leaf", "line": list(line), "b": hex(b),
                                   "case": label, "control": control,
                                   "p": list(left), "q": list(right), "r": list(result),
                                   "lambda": hex(slope), "accepted": accepted,
                                   "expected": expected}
                    assert accepted == expected
                    count += 1
            offcurve = None
            for bit in range(131):
                candidate = (0, p[1], p[2] ^ (1 << bit))
                if not curve.affine(candidate[1], candidate[2]):
                    offcurve = candidate
                    break
            assert offcurve is not None
            for control, bad in (("off_curve", offcurve),
                                 ("noncanonical_infinity", (1, 1, 0))):
                accepted = relation.accepts(bad, O, O, 0)
                row = next_row(rows)
                assert row == {"kind": "leaf_invalid", "line": list(line),
                               "b": hex(b), "control": control, "point": list(bad),
                               "accepted": accepted, "expected": False}
                assert not accepted
                count += 1
            actual = {"line": list(line), "b": hex(b), "rows": count,
                      "case_counts": counts, "dag": relation.dag.counts(),
                      "model_controls": model_controls(relation),
                      "mismatches": 0}
            assert summary["leaves"][len(leaf_counts)] == actual
            leaf_counts.append(actual)
            total_rows += count
        assert rows.readline() == "", "extra archived rows"
    assert len(leaf_counts) == 2 and total_rows == summary["raw_rows"]
    receipt = {"domain": DOMAIN, "decision": "PASS", "toy": toy_counts,
               "leaves": leaf_counts, "raw_rows": total_rows,
               "rows_sha256": summary["rows_sha256"],
               "producer_sha256": sha(producer / "result.json"),
               "wall_seconds": time.monotonic() - started,
               "cpu_seconds": time.process_time() - cpu_started,
               "peak_rss_bytes": rss_bytes()}
    out.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(replay(args.producer, args.out), sort_keys=True))
