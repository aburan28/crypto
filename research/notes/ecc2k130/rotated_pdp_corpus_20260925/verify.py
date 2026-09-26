#!/usr/bin/env python3
"""Independent bit-serial/Fermat full-tuple replay of the fixed point corpus."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925"
_spec = importlib.util.spec_from_file_location("rotated_parent_verify", PARENT / "verify.py")
assert _spec is not None and _spec.loader is not None
parent_verify = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(parent_verify)

DOMAIN = "ECC2K130-ROTATED-PDP-CORPUS-20260925-v1"
ARMS = {
    "n13-m5": (13, 0x201b, 2003, 89, 5, 2),
    "n19-m6": (19, 0x80027, 130873, 41811, 6, 2),
}
CAP_SECONDS = 600
CAP_RSS = 512 * 1024 * 1024
Point = tuple[int, int] | None


def pjson(point: Point):
    return None if point is None else list(point)


def point(value):
    return None if value is None else tuple(value)


def key(value: Point):
    return (-1, -1) if value is None else value


def load(path: Path):
    return json.loads(path.read_text())


def jsonl(path: Path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def peak_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def is_prime(n: int) -> bool:
    return n > 1 and all(n % p for p in range(2, math.isqrt(n) + 1))


def poly_rem(a: int, modulus: int) -> int:
    while a.bit_length() >= modulus.bit_length():
        a ^= modulus << (a.bit_length() - modulus.bit_length())
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_rem(a, b)
    return a


def rabin_prime_degree(f: parent_verify.GF) -> None:
    assert is_prime(f.n)
    assert f.poly & 1 and poly_gcd(f.poly, 0b110) == 1
    x = 2
    for _ in range(f.n):
        x = f.square(x)
    assert x == 2


def source_order(n: int) -> int:
    previous, current = 2, -1
    for _ in range(2, n + 1):
        previous, current = current, -(current + 2 * previous)
    return (1 << n) + 1 - current


def rank(vectors: list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in vectors:
        while value:
            leading = value.bit_length() - 1
            if leading in pivots:
                value ^= pivots[leading]
            else:
                pivots[leading] = value
                break
    return len(pivots)


def half_trace(f: parent_verify.GF, rhs: int) -> int:
    assert f.n & 1 and f.trace(rhs) == 0
    acc, term = 0, rhs
    for _ in range((f.n + 1) // 2):
        acc ^= term
        term = f.square(f.square(term))
    assert f.square(acc) ^ acc == rhs
    return acc


def lifts(curve: parent_verify.E, x: int) -> list[tuple[int, int]]:
    f = curve.f
    if x == 0:
        return [(0, 1)]
    rhs = x ^ f.square(f.inv(x))
    if f.trace(rhs):
        return []
    z = half_trace(f, rhs)
    points = sorted([(x, f.mul(x, z)), (x, f.mul(x, z ^ 1))])
    assert all(curve.on(p) for p in points) and len(set(points)) == 2
    return points


def generator(curve: parent_verify.E, q: int) -> Point:
    for x in range(1 << curve.f.n):
        for p in lifts(curve, x):
            h = curve.scalar(p, 4)
            if h is not None:
                assert curve.scalar(h, q) is None
                return h
    raise AssertionError("no projected generator")


def neg(p: Point) -> Point:
    return None if p is None else (p[0], p[0] ^ p[1])


def ordinal_indices(ordinal: int, sizes: list[int]) -> list[int]:
    indices = []
    for size in reversed(sizes):
        ordinal, digit = divmod(ordinal, size)
        indices.append(digit)
    assert ordinal == 0
    return list(reversed(indices))


def candidate(arm: str, kind: str, counter: int, modulus: int) -> int:
    payload = f"{DOMAIN}/{arm}/{kind}/{counter}".encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest(), "big") % modulus


def tuple_sum(curve: parent_verify.E, factors, indices: list[int]) -> Point:
    acc = None
    for factor, index in zip(factors, indices):
        acc = curve.add(acc, factor[index])
    return acc


def factor_points(curve: parent_verify.E, beta: int, m: int, d: int):
    f = curve.f
    conjugates = []
    value = beta
    for _ in range(f.n):
        conjugates.append(value)
        value = f.square(value)
    assert value == beta and rank(conjugates) == f.n and f.trace(beta) == 1
    bases = [[conjugates[m * j + i] for j in range(d)] for i in range(m)]
    assert all(rank(b) == d for b in bases)
    assert rank([v for b in bases for v in b]) == m * d
    factors = []
    for basis in bases:
        xs = [0, basis[0], basis[1], basis[0] ^ basis[1]]
        assert len(set(xs)) == 4
        factors.append(sorted(p for x in xs for p in lifts(curve, x)))
    assert all(factors[i + 1] == sorted(curve.tau(p) for p in factors[i])
               for i in range(m - 1))
    return factors


def verify_selection(curve: parent_verify.E, arm: str, q: int, h: Point,
                     factors, full, projected, rows, attempts) -> tuple[int, int]:
    assert len(rows) == 8 and [row["class"] for row in rows] == ["planted"] * 4 + ["negative"] * 4
    sizes = list(map(len, factors))
    total = math.prod(sizes)
    inverse_four = pow(4, -1, q)
    torsion = [None, (0, 1), (1, 0), (1, 1)]
    expected_attempts, expected_rows = [], []
    seen_planted, seen_negative, tried_negative_k = set(), set(), set()
    for counter in range(100000):
        ordinal = candidate(arm, "planted", counter, total)
        indices = ordinal_indices(ordinal, sizes)
        s = tuple_sum(curve, factors, indices)
        r = curve.scalar(s, 4)
        accepted = r is not None and r not in seen_planted
        expected_attempts.append({"class": "planted", "counter": counter,
                                  "ordinal": ordinal, "tuple_indices": indices,
                                  "R": pjson(r), "accepted": accepted})
        if accepted:
            seen_planted.add(r)
            qpoint = curve.scalar(r, inverse_four)
            assert curve.scalar(qpoint, 4) == r
            shift = curve.add(s, neg(qpoint))
            assert shift in torsion and curve.add(qpoint, shift) == s
            expected_rows.append({"class": "planted", "counter": counter,
                                  "tuple_ordinal": ordinal, "tuple_indices": indices,
                                  "S": pjson(s), "Q": pjson(qpoint), "R": pjson(r),
                                  "T": pjson(shift)})
            if len(seen_planted) == 4:
                break
    assert len(seen_planted) == 4
    for counter in range(100000):
        k = 1 + candidate(arm, "negative", counter, q - 1)
        qpoint = curve.scalar(h, k)
        r = curve.scalar(qpoint, 4)
        duplicate = k in tried_negative_k
        tried_negative_k.add(k)
        accepted = not duplicate and r not in projected
        expected_attempts.append({"class": "negative", "counter": counter,
                                  "k": k, "Q": pjson(qpoint), "R": pjson(r),
                                  "duplicate": duplicate, "supported": r in projected,
                                  "accepted": accepted})
        if accepted:
            seen_negative.add(k)
            expected_rows.append({"class": "negative", "counter": counter,
                                  "k": k, "Q": pjson(qpoint), "R": pjson(r)})
            if len(seen_negative) == 4:
                break
    assert len(seen_negative) == 4 and attempts == expected_attempts
    for actual, expected in zip(rows, expected_rows):
        for k, v in expected.items():
            assert actual[k] == v, (arm, k, actual[k], v)
        qpoint, r = point(actual["Q"]), point(actual["R"])
        assert curve.scalar(qpoint, 4) == r
        cosets = [curve.add(qpoint, t) for t in torsion]
        assert len(set(cosets)) == 4
        counts = [full.get(p, 0) for p in cosets]
        assert actual["coset_multiplicities"] == counts
        assert actual["projected_multiplicity"] == projected.get(r, 0) == sum(counts)
        for i, coset in enumerate(cosets):
            witness = actual["coset_witness_indices"][i]
            if counts[i]:
                assert witness is not None and tuple_sum(curve, factors, witness) == coset
            else:
                assert witness is None
        assert (r in projected) == (actual["class"] == "planted")
    return len(seen_planted), len(seen_negative)


def verify_arm(arm: str, archive: Path) -> dict:
    started, cpu_started = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"{arm}: {CAP_SECONDS}s verifier cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        n, poly, q, lam, m, d = ARMS[arm]
        assert source_order(n) == 4 * q and is_prime(q)
        f = parent_verify.GF(n, poly)
        rabin_prime_degree(f)
        curve = parent_verify.E(f)
        torsion = [None, (0, 1), (1, 0), (1, 1)]
        assert all(curve.on(t) and curve.scalar(t, 4) is None for t in torsion)
        h = generator(curve, q)
        if n == 13:
            assert h == (4793, 2429)
        assert (lam * lam + lam + 2) % q == 0 and pow(lam, n, q) == 1
        assert curve.tau(h) == curve.scalar(h, lam)
        factors = factor_points(curve, 3, m, d)
        assert load(archive / "factors.json") == [[pjson(p) for p in b] for b in factors]
        sizes = list(map(len, factors))
        full = Counter()
        for choice in itertools.product(*factors):
            acc = None
            for p in choice:
                acc = curve.add(acc, p)
            full[acc] += 1
        assert sum(full.values()) == math.prod(sizes) <= 7 ** m
        projected = Counter()
        for s, count in full.items():
            projected[curve.scalar(s, 4)] += count
        frozen_full = jsonl(archive / "full_histogram.jsonl")
        assert [point(r["point"]) for r in frozen_full] == sorted(full, key=key)
        assert {point(r["point"]): r["count"] for r in frozen_full} == dict(full)
        for row in frozen_full:
            assert tuple_sum(curve, factors, row["witness_indices"]) == point(row["point"])
        frozen_projected = jsonl(archive / "projected_histogram.jsonl")
        assert [point(r["point"]) for r in frozen_projected] == sorted(projected, key=key)
        assert {point(r["point"]): r["count"] for r in frozen_projected} == dict(projected)
        for row in frozen_projected:
            s = tuple_sum(curve, factors, row["witness_indices"])
            assert s == point(row["full_sum"])
            assert curve.scalar(s, 4) == point(row["point"])
        targets = load(archive / "targets.json")
        attempts = jsonl(archive / "selection_attempts.jsonl")
        planted, negative = verify_selection(curve, arm, q, h, factors,
                                              full, projected, targets, attempts)
        summary = load(archive / "summary.json")
        assert summary["arm"] == arm and summary["n"] == n and summary["poly"] == poly
        assert summary["q"] == q and summary["lambda"] == lam
        assert summary["m"] == m and summary["d"] == d and summary["beta"] == 3
        assert point(summary["generator_H"]) == h
        assert summary["factor_sizes"] == sizes
        assert summary["tuple_count"] == sum(full.values())
        assert summary["distinct_full_sums"] == len(full)
        assert summary["distinct_projected_sums"] == len(projected)
        assert summary["projected_nonzero_support"] == len(projected) - int(None in projected)
        assert summary["projected_negative_guarantee"] == q - min(q, 7 ** m)
        assert summary["selection_attempts"] == len(attempts)
        assert summary["planted"] == planted and summary["negative"] == negative
        assert summary["total_wall_seconds"] <= 300 and summary["peak_rss_bytes"] <= CAP_RSS
        result = {"arm": arm, "full_histogram_replayed": len(full),
                  "projected_histogram_replayed": len(projected),
                  "labelled_tuples_replayed": sum(full.values()),
                  "planted_witnesses": planted, "complete_negative_targets": negative,
                  "selection_attempts_replayed": len(attempts),
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu_started,
                  "peak_rss_bytes": peak_rss(),
                  "operations": {"field": dict(f.ops), "curve": dict(curve.ops)}}
        assert result["wall_seconds"] <= CAP_SECONDS and result["peak_rss_bytes"] <= CAP_RSS
        return result
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=sorted(ARMS), required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = verify_arm(args.arm, args.archive)
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
