#!/usr/bin/env python3
"""Frozen full group-law oracle for the rotated m5/m6 point corpus."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import resource
import signal
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925"
sys.path.insert(0, str(PARENT))
import gate as parent  # noqa: E402 - pinned to PR #762 by run.py

DOMAIN = "ECC2K130-ROTATED-PDP-CORPUS-20260925-v1"
ARMS = {
    "n13-m5": {"n": 13, "low": [0, 1, 3, 4], "q": 2003,
               "lambda": 89, "beta": 3, "m": 5, "d": 2},
    "n19-m6": {"n": 19, "low": [0, 1, 2, 5], "q": 130873,
               "lambda": 41811, "beta": 3, "m": 6, "d": 2},
}
CAP_SECONDS = 300
CAP_RSS = 512 * 1024 * 1024
Point = tuple[int, int] | None


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def pjson(point: Point):
    return None if point is None else list(point)


def key(point: Point):
    return (-1, -1) if point is None else point


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def is_prime(q: int) -> bool:
    return q > 1 and all(q % j for j in range(2, math.isqrt(q) + 1))


def half_trace(f: parent.Field, a: int) -> int:
    assert f.n & 1 and f.trace(a) == 0
    z, term = 0, a
    for _ in range((f.n + 1) // 2):
        z ^= term
        term = f.square(f.square(term))
    assert f.square(z) ^ z == a
    return z


def lifts(curve: parent.Curve, x: int) -> list[tuple[int, int]]:
    f = curve.f
    if x == 0:
        return [(0, 1)]
    rhs = x ^ f.square(f.inverse(x))
    if f.trace(rhs):
        return []
    z = half_trace(f, rhs)
    points = sorted([(x, f.mul(x, z)), (x, f.mul(x, z ^ 1))])
    assert len(set(points)) == 2 and all(curve.on_curve(p) for p in points)
    return points


def generator(curve: parent.Curve, q: int) -> tuple[int, int]:
    for x in range(1 << curve.f.n):
        for p in lifts(curve, x):
            h = curve.scalar(p, 4)
            if h is not None:
                assert curve.scalar(h, q) is None
                return h
    raise AssertionError("no projected generator")


def digest_index(arm: str, kind: str, counter: int, modulus: int) -> int:
    seed = f"{DOMAIN}/{arm}/{kind}/{counter}".encode("ascii")
    return int.from_bytes(hashlib.sha256(seed).digest(), "big") % modulus


def decode_ordinal(ordinal: int, sizes: list[int]) -> list[int]:
    indices = [0] * len(sizes)
    for i in range(len(sizes) - 1, -1, -1):
        ordinal, indices[i] = divmod(ordinal, sizes[i])
    assert ordinal == 0
    return indices


def tuple_sum(curve: parent.Curve, factors, indices: list[int]) -> Point:
    result = None
    for factor, index in zip(factors, indices):
        result = curve.add(result, factor[index])
    return result


def select_targets(curve: parent.Curve, arm: str, q: int, h: Point,
                   factors, full: dict[Point, int], full_witness,
                   projected: dict[Point, int]):
    sizes = list(map(len, factors))
    total = math.prod(sizes)
    inverse_four = pow(4, -1, q)
    torsion = [None, (0, 1), (1, 0), (1, 1)]
    rows, attempts = [], []
    seen_planted, seen_negative, tried_negative_k = set(), set(), set()
    for counter in range(100000):
        ordinal = digest_index(arm, "planted", counter, total)
        indices = decode_ordinal(ordinal, sizes)
        s = tuple_sum(curve, factors, indices)
        r = curve.scalar(s, 4)
        accepted = r is not None and r not in seen_planted
        attempts.append({"class": "planted", "counter": counter,
                         "ordinal": ordinal, "tuple_indices": indices,
                         "R": pjson(r), "accepted": accepted})
        if not accepted:
            continue
        seen_planted.add(r)
        qpoint = curve.scalar(r, inverse_four)
        assert curve.scalar(qpoint, 4) == r
        shift = curve.add(s, curve.neg(qpoint))
        assert shift in torsion and curve.add(qpoint, shift) == s
        rows.append({"class": "planted", "counter": counter,
                     "tuple_ordinal": ordinal, "tuple_indices": indices,
                     "S": pjson(s), "Q": pjson(qpoint), "R": pjson(r),
                     "T": pjson(shift)})
        if len(seen_planted) == 4:
            break
    assert len(seen_planted) == 4, "planted stream exhausted"
    for counter in range(100000):
        k = 1 + digest_index(arm, "negative", counter, q - 1)
        qpoint = curve.scalar(h, k)
        r = curve.scalar(qpoint, 4)
        duplicate = k in tried_negative_k
        tried_negative_k.add(k)
        accepted = not duplicate and r not in projected
        attempts.append({"class": "negative", "counter": counter,
                         "k": k, "Q": pjson(qpoint), "R": pjson(r),
                         "duplicate": duplicate, "supported": r in projected,
                         "accepted": accepted})
        if not accepted:
            continue
        seen_negative.add(k)
        rows.append({"class": "negative", "counter": counter, "k": k,
                     "Q": pjson(qpoint), "R": pjson(r)})
        if len(seen_negative) == 4:
            break
    assert len(seen_negative) == 4, "negative stream exhausted"
    for row in rows:
        qpoint = None if row["Q"] is None else tuple(row["Q"])
        r = None if row["R"] is None else tuple(row["R"])
        assert curve.scalar(qpoint, 4) == r
        cosets = [curve.add(qpoint, t) for t in torsion]
        assert len(set(cosets)) == 4
        row["coset_multiplicities"] = [full.get(p, 0) for p in cosets]
        row["coset_witness_indices"] = [list(full_witness[p]) if p in full_witness else None
                                         for p in cosets]
        row["projected_multiplicity"] = projected.get(r, 0)
        assert sum(row["coset_multiplicities"]) == row["projected_multiplicity"]
        assert (row["projected_multiplicity"] > 0) == (row["class"] == "planted")
    return rows, attempts


def run_arm(arm: str, out: Path) -> None:
    cfg = ARMS[arm]
    started, cpu_started = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"{arm}: {CAP_SECONDS}s wall acceptance cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        n, q, m, d = cfg["n"], cfg["q"], cfg["m"], cfg["d"]
        f = parent.Field(n, cfg["low"])
        f.rabin_prime_degree()
        assert f.poly == (0x201b if n == 13 else 0x80027)
        assert is_prime(q) and parent.source_group_order(n) == 4 * q
        c = parent.Curve(f)
        torsion = [None, (0, 1), (1, 0), (1, 1)]
        assert all(c.on_curve(t) and c.scalar(t, 4) is None for t in torsion)
        h = generator(c, q)
        if n == 13:
            assert h == (4793, 2429)
        lam = cfg["lambda"]
        assert (lam * lam + lam + 2) % q == 0 and pow(lam, n, q) == 1
        assert c.tau(h) == c.scalar(h, lam)
        conjugates = parent.normal_conjugates(f, cfg["beta"])
        bases = parent.subspace_basis(conjugates, m, d)
        assert parent.rank([v for basis in bases for v in basis]) == m * d
        factors = []
        for basis in bases:
            xs = parent.all_x(basis)
            assert len(set(xs)) == 4
            factors.append(sorted(p for x in xs for p in lifts(c, x)))
        assert all(factors[i + 1] == sorted(c.tau(p) for p in factors[i])
                   for i in range(m - 1))
        sizes = list(map(len, factors))
        assert all(1 <= size <= 7 for size in sizes)
        setup_wall = time.perf_counter() - started
        setup_ops = {"field": dict(f.operations), "curve": dict(c.operations)}

        full, witnesses = parent.full_histogram(c, factors)
        tuple_count = math.prod(sizes)
        assert sum(full.values()) == tuple_count <= 7 ** m
        histogram_wall = time.perf_counter() - started - setup_wall
        projected, projected_witness = {}, {}
        for s, count in full.items():
            r = c.scalar(s, 4)
            projected[r] = projected.get(r, 0) + count
            projected_witness.setdefault(r, (s, witnesses[s]))
        assert sum(projected.values()) == tuple_count
        projection_wall = time.perf_counter() - started - setup_wall - histogram_wall
        rows, attempts = select_targets(c, arm, q, h, factors, full,
                                        witnesses, projected)
        assert [row["class"] for row in rows] == ["planted"] * 4 + ["negative"] * 4
        selection_wall = time.perf_counter() - started - setup_wall - histogram_wall - projection_wall
        save(out / "factors.json", [[pjson(p) for p in factor] for factor in factors])
        save(out / "targets.json", rows)
        with (out / "selection_attempts.jsonl").open("w") as stream:
            for row in attempts:
                stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        with (out / "full_histogram.jsonl").open("w") as stream:
            for p in sorted(full, key=key):
                stream.write(json.dumps({"point": pjson(p), "count": full[p],
                                         "witness_indices": list(witnesses[p])},
                                        sort_keys=True, separators=(",", ":")) + "\n")
        with (out / "projected_histogram.jsonl").open("w") as stream:
            for p in sorted(projected, key=key):
                s, witness = projected_witness[p]
                stream.write(json.dumps({"point": pjson(p), "count": projected[p],
                                         "full_sum": pjson(s), "witness_indices": list(witness)},
                                        sort_keys=True, separators=(",", ":")) + "\n")
        summary = {"arm": arm, "n": n, "poly": f.poly, "q": q, "beta": cfg["beta"],
                   "m": m, "d": d, "lambda": lam, "generator_H": pjson(h),
                   "factor_sizes": sizes, "tuple_count": tuple_count,
                   "distinct_full_sums": len(full), "distinct_projected_sums": len(projected),
                   "projected_nonzero_support": len(projected) - int(None in projected),
                   "projected_negative_guarantee": q - min(q, 7 ** m),
                   "planted": 4, "negative": 4, "selection_attempts": len(attempts),
                   "setup_wall_seconds": setup_wall, "histogram_wall_seconds": histogram_wall,
                   "projection_wall_seconds": projection_wall,
                   "selection_wall_seconds": selection_wall,
                   "total_wall_seconds": time.perf_counter() - started,
                   "total_cpu_seconds": time.process_time() - cpu_started,
                   "peak_rss_bytes": rss(),
                   "setup_operations": setup_ops,
                   "total_operations": {"field": dict(f.operations), "curve": dict(c.operations)}}
        save(out / "summary.json", summary)
        assert summary["total_wall_seconds"] <= CAP_SECONDS and rss() <= CAP_RSS
    except Exception as error:
        save(out / "producer_failure.json", {"arm": arm, "error": repr(error),
                                               "wall_seconds": time.perf_counter() - started,
                                               "cpu_seconds": time.process_time() - cpu_started,
                                               "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=sorted(ARMS), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    run_arm(args.arm, args.out)


if __name__ == "__main__":
    main()
