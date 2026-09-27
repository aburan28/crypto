#!/usr/bin/env python3
"""Frozen exhaustive producer for the small-field full-point relation."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

from dag import Dag, PackedModel, build_relation

FIELDS = ((2, 0x7), (3, 0xB))
O = (1, 0, 0)


class RefField:
    """Bit-serial, shift/reduce oracle; deliberately separate from the DAG."""

    def __init__(self, n: int, modulus: int) -> None:
        self.n, self.modulus = n, modulus

    def mul(self, a: int, b: int) -> int:
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & (1 << self.n):
                a ^= self.modulus
        return result

    def inv(self, a: int) -> int:
        if not a:
            raise ZeroDivisionError
        for candidate in range(1, 1 << self.n):
            if self.mul(a, candidate) == 1:
                return candidate
        raise AssertionError("input was not a field")

    def on_curve(self, x: int, y: int) -> bool:
        xx = self.mul(x, x)
        return self.mul(y, y) ^ self.mul(x, y) == self.mul(xx, x) ^ 1

    def points(self) -> list[tuple[int, int, int]]:
        return [O] + [(0, x, y) for x in range(1 << self.n)
                      for y in range(1 << self.n) if self.on_curve(x, y)]

    def add(self, p: tuple[int, int, int], q: tuple[int, int, int]) -> tuple[int, int, int]:
        if p[0]:
            return q
        if q[0]:
            return p
        _, x1, y1 = p
        _, x2, y2 = q
        if x1 == x2 and y1 ^ y2 == x1:
            return O
        if x1 == x2:
            slope = x1 ^ self.mul(y1, self.inv(x1))
            x3 = self.mul(slope, slope) ^ slope
            y3 = self.mul(x1, x1) ^ self.mul(slope ^ 1, x3)
        else:
            slope = self.mul(y1 ^ y2, self.inv(x1 ^ x2))
            x3 = self.mul(slope, slope) ^ slope ^ x1 ^ x2
            y3 = self.mul(slope, x1 ^ x3) ^ x3 ^ y1
        return (0, x3, y3)


def _must_reject(call) -> None:
    try:
        call()
    except ValueError:
        return
    raise AssertionError("malformed model or coordinate was accepted")


def model_controls(relation, n: int) -> int:
    d = relation.dag
    count = 0
    zero = relation.model(O, O, O, 0)
    _must_reject(lambda: d.evaluate(PackedModel(zero.bit_count - 1, zero.limbs), relation.output))
    count += 1
    _must_reject(lambda: d.evaluate(PackedModel(zero.bit_count, ()), relation.output))
    count += 1
    _must_reject(lambda: d.evaluate(PackedModel(zero.bit_count, (1 << 64,)), relation.output))
    count += 1
    _must_reject(lambda: d.evaluate(PackedModel(zero.bit_count, (1 << zero.bit_count,)), relation.output))
    count += 1
    _must_reject(lambda: relation.model(O, O, (0, 1 << n, 0), 0))
    count += 1
    _must_reject(lambda: relation.model(O, O, O, 1 << n))
    count += 1
    _must_reject(lambda: relation.model((2, 0, 0), O, O, 0))
    count += 1
    wide = Dag()
    variables = [wide.var(f"v{i}") for i in range(920)]
    result = wide.xor(variables[0], variables[-1])
    bits = [0] * 920
    bits[-1] = 1
    packed = PackedModel.from_bits(bits)
    assert len(packed.limbs) == 15 and wide.evaluate(packed, result)
    count += 1
    _must_reject(lambda: wide.evaluate(PackedModel(920, packed.limbs[:-1]), result))
    count += 1
    last = list(packed.limbs)
    last[-1] |= 1 << 63
    _must_reject(lambda: wide.evaluate(PackedModel(920, tuple(last)), result))
    count += 1
    return count


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def run(out: Path) -> dict:
    started = time.monotonic()
    cpu_started = time.process_time()
    out.mkdir(parents=True, exist_ok=False)
    rows_path = out / "rows.jsonl.gz"
    summaries = []
    first_mismatch = None
    with rows_path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        for n, modulus in FIELDS:
            relation = build_relation(n, modulus)
            ref = RefField(n, modulus)
            points = ref.points()
            universe = [(o, x, y) for o in (0, 1) for x in range(1 << n)
                        for y in range(1 << n)]
            invalid = [p for p in universe if p not in points]
            case_counts = {case: 0 for case in relation.branches}
            rows = evaluations = invalid_evaluations = mismatches = 0
            for p in points:
                for q in points:
                    result = ref.add(p, q)
                    if result not in points:
                        raise AssertionError("reference group law left the curve")
                    if p[0]:
                        case = "copy_q"
                    elif q[0]:
                        case = "copy_p"
                    elif p[1] == q[1] and p[2] ^ q[2] == p[1]:
                        case = "inverse"
                    elif p[1] == q[1]:
                        case = "double"
                    else:
                        case = "generic"
                    case_counts[case] += 1
                    for r in points:
                        accepted = sum(relation.accepts(p, q, r, lam)
                                       for lam in range(1 << n))
                        evaluations += 1 << n
                        expect = int(r == result)
                        row = {"n": n, "p": p, "q": q, "r": r,
                               "witnesses": accepted, "expected": expect, "case": case}
                        gz.write((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode())
                        rows += 1
                        expected_witnesses = (1 << n) if case in ("copy_q", "copy_p", "inverse") else 1
                        if accepted != expect * expected_witnesses:
                            mismatches += 1
                            if first_mismatch is None:
                                first_mismatch = row
            # Every malformed or off-curve point is tested in each of the
            # three roles, with every slope; no role is silently trusted.
            for bad in invalid:
                for p, q, r in ((bad, O, O), (O, bad, O), (O, O, bad)):
                    for lam in range(1 << n):
                        if relation.accepts(p, q, r, lam):
                            mismatches += 1
                            if first_mismatch is None:
                                first_mismatch = {"n": n, "invalid": bad,
                                                  "p": p, "q": q, "r": r, "lambda": lam}
                        invalid_evaluations += 1
            controls = model_controls(relation, n)
            summaries.append({"n": n, "modulus": modulus, "curve_points": len(points),
                              "rows": rows, "valid_model_evaluations": evaluations,
                              "invalid_points": len(invalid),
                              "invalid_model_evaluations": invalid_evaluations,
                              "model_controls": controls, "case_counts": case_counts,
                              "dag": relation.dag.counts(), "mismatches": mismatches})
    for case in summaries[0]["case_counts"]:
        if sum(field["case_counts"][case] for field in summaries) == 0:
            first_mismatch = first_mismatch or {"missing_case": case}
    result = {"domain": "k0-symbolic-dag-fullpoint-n2n3-v1",
              "fields": summaries, "first_mismatch": first_mismatch,
              "rows_sha256": hashlib.sha256(rows_path.read_bytes()).hexdigest(),
              "wall_seconds": time.monotonic() - started,
              "cpu_seconds": time.process_time() - cpu_started,
              "peak_rss_bytes": _rss_bytes(),
              "decision": "PASS" if first_mismatch is None else "FAIL"}
    (out / "result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.out), sort_keys=True))
