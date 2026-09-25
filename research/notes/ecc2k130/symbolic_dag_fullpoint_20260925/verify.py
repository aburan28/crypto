#!/usr/bin/env python3
"""Independent field/group oracle and full archive replay for the DAG relation."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

from dag import PackedModel, build_relation

FIELDS = ((2, 0x7), (3, 0xB))
O = (1, 0, 0)


class PolynomialField:
    """Product-then-long-division arithmetic, Euclid inverse."""

    def __init__(self, n: int, modulus: int) -> None:
        self.n, self.modulus = n, modulus

    def mul(self, a: int, b: int) -> int:
        product = 0
        for i in range(self.n):
            for j in range(self.n):
                if (a >> i & 1) and (b >> j & 1):
                    product ^= 1 << (i + j)
        while product.bit_length() > self.n:
            product ^= self.modulus << (product.bit_length() - self.n - 1)
        return product

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError
        left, right = self.modulus, a
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
            raise AssertionError("reducible modulus")
        while u.bit_length() > self.n:
            u ^= self.modulus << (u.bit_length() - self.n - 1)
        if self.mul(a, u) != 1:
            raise AssertionError("Euclid inverse failed")
        return u

    def valid_affine(self, x: int, y: int) -> bool:
        return self.mul(y, y) ^ self.mul(x, y) == self.mul(self.mul(x, x), x) ^ 1

    def points(self) -> list[tuple[int, int, int]]:
        affine = []
        for y in range(1 << self.n):
            for x in range(1 << self.n):
                if self.valid_affine(x, y):
                    affine.append((0, x, y))
        return [O] + sorted(affine)

    def sum(self, p: tuple[int, int, int], q: tuple[int, int, int]) -> tuple[int, int, int]:
        if p == O:
            return q
        if q == O:
            return p
        _, x, y = p
        _, u, v = q
        if x == u:
            if y ^ v == x:
                return O
            if p != q or x == 0:
                raise AssertionError("same-x finite pair violates curve fibre")
            slope = x ^ self.mul(y, self.inv(x))
            new_x = self.mul(slope, slope) ^ slope
            new_y = self.mul(x, x) ^ self.mul(slope ^ 1, new_x)
        else:
            slope = self.mul(y ^ v, self.inv(x ^ u))
            new_x = self.mul(slope, slope) ^ slope ^ x ^ u
            new_y = self.mul(slope, x ^ new_x) ^ new_x ^ y
        return (0, new_x, new_y)


def _reject(action) -> None:
    try:
        action()
    except ValueError:
        return
    raise AssertionError("malformed model was accepted")


def check_models(relation, n: int) -> int:
    d = relation.dag
    packed = relation.model(O, O, O, 0)
    assert packed.bit_count == 7 * n + 3
    assert len(packed.limbs) == (packed.bit_count + 63) // 64
    probes = (
        lambda: d.evaluate(PackedModel(packed.bit_count - 1, packed.limbs), relation.output),
        lambda: d.evaluate(PackedModel(packed.bit_count, ()), relation.output),
        lambda: d.evaluate(PackedModel(packed.bit_count, (1 << 64,)), relation.output),
        lambda: d.evaluate(PackedModel(packed.bit_count, (1 << packed.bit_count,)), relation.output),
        lambda: relation.model(O, O, (0, 1 << n, 0), 0),
        lambda: relation.model(O, O, O, 1 << n),
        lambda: relation.model((2, 0, 0), O, O, 0),
    )
    for probe in probes:
        _reject(probe)
    # Synthetic n131-width assignment exercises the highest input, limb 14,
    # without constructing or making any claim about an n131 curve circuit.
    from dag import Dag
    wide = Dag()
    first = wide.var("first")
    for i in range(1, 919):
        wide.var(f"middle_{i}")
    last = wide.var("last")
    output = wide.xor(first, last)
    bits = [0] * 920
    bits[919] = 1
    model = PackedModel.from_bits(bits)
    assert len(model.limbs) == 15 and wide.evaluate(model, output)
    _reject(lambda: wide.evaluate(PackedModel(920, model.limbs[:-1]), output))
    corrupt = list(model.limbs)
    corrupt[14] |= 1 << 63
    _reject(lambda: wide.evaluate(PackedModel(920, tuple(corrupt)), output))
    return len(probes) + 3


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def replay(producer: Path, out: Path) -> dict:
    started = time.monotonic()
    cpu_started = time.process_time()
    summary = json.loads((producer / "result.json").read_text())
    assert summary["domain"] == "k0-symbolic-dag-fullpoint-n2n3-v1"
    assert summary["decision"] == "PASS" and summary["first_mismatch"] is None
    rows_file = producer / "rows.jsonl.gz"
    assert hashlib.sha256(rows_file.read_bytes()).hexdigest() == summary["rows_sha256"]
    counts = []
    with gzip.open(rows_file, "rt") as rows:
        for field_index, (n, modulus) in enumerate(FIELDS):
            f = PolynomialField(n, modulus)
            assert all(f.mul(a, f.inv(a)) == 1 for a in range(1, 1 << n))
            points = f.points()
            relation = build_relation(n, modulus)
            recorded = summary["fields"][field_index]
            assert recorded["n"] == n and recorded["modulus"] == modulus
            assert recorded["curve_points"] == len(points)
            assert recorded["dag"] == relation.dag.counts()
            assert recorded["model_controls"] == check_models(relation, n)
            branch_counts = dict.fromkeys(relation.branches, 0)
            actual_rows = valid_evals = invalid_evals = 0
            for p in points:
                for q in points:
                    expected_point = f.sum(p, q)
                    assert expected_point in points
                    if p == O:
                        case = "copy_q"
                    elif q == O:
                        case = "copy_p"
                    elif p[1] == q[1] and p[2] ^ q[2] == p[1]:
                        case = "inverse"
                    elif p[1] == q[1]:
                        case = "double"
                    else:
                        case = "generic"
                    branch_counts[case] += 1
                    for r in points:
                        line = rows.readline()
                        if not line:
                            raise AssertionError("truncated raw row archive")
                        row = json.loads(line)
                        assert row["n"] == n and tuple(row["p"]) == p
                        assert tuple(row["q"]) == q and tuple(row["r"]) == r
                        assert row["case"] == case
                        assert row["expected"] == int(r == expected_point)
                        truth = [relation.accepts(p, q, r, lam) for lam in range(1 << n)]
                        assert row["witnesses"] == sum(truth)
                        expected_count = (1 << n) if case in ("copy_p", "copy_q", "inverse") else 1
                        assert sum(truth) == (expected_count if r == expected_point else 0)
                        valid_evals += 1 << n
                        actual_rows += 1
            universe = [(o, x, y) for o in (0, 1) for x in range(1 << n)
                        for y in range(1 << n)]
            invalid = [p for p in universe if p not in points]
            for bad in invalid:
                for p, q, r in ((bad, O, O), (O, bad, O), (O, O, bad)):
                    for lam in range(1 << n):
                        assert not relation.accepts(p, q, r, lam)
                        invalid_evals += 1
            assert recorded["rows"] == actual_rows
            assert recorded["valid_model_evaluations"] == valid_evals
            assert recorded["invalid_points"] == len(invalid)
            assert recorded["invalid_model_evaluations"] == invalid_evals
            assert recorded["case_counts"] == branch_counts
            assert recorded["mismatches"] == 0
            counts.append({"n": n, "rows": actual_rows,
                           "valid_model_evaluations": valid_evals,
                           "invalid_model_evaluations": invalid_evals,
                           "curve_points": len(points), "case_counts": branch_counts,
                           "dag": relation.dag.counts()})
        assert rows.readline() == "", "extra raw row archive records"
    assert all(sum(field["case_counts"][case] for field in counts) > 0
               for case in ("copy_q", "copy_p", "inverse", "double", "generic"))
    result = {"domain": summary["domain"], "decision": "PASS", "fields": counts,
              "producer_sha256": hashlib.sha256((producer / "result.json").read_bytes()).hexdigest(),
              "rows_sha256": summary["rows_sha256"],
              "wall_seconds": time.monotonic() - started,
              "cpu_seconds": time.process_time() - cpu_started,
              "peak_rss_bytes": _rss_bytes()}
    out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(replay(args.producer, args.out), sort_keys=True))
