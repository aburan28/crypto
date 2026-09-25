#!/usr/bin/env python3
"""Independent Euclid-field/group replay of finite-fibre S3 controls."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import itertools
import json
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925/gate.py"
FIELDS = ((1, 0x3), (2, 0x7), (3, 0xb), (4, 0x13), (5, 0x25), (6, 0x43), (7, 0x83))
PANELS = {3: 4, 4: 4, 5: 3}
WALL_CAP = 180
RSS_CAP = 512 * 1024 * 1024


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def parent_module():
    spec = importlib.util.spec_from_file_location("theorem_independent_curve", PARENT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def poly_rem(a: int, b: int) -> int:
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_rem(a, b)
    return a


def verify_field(f):
    n, poly = f.n, f.poly
    assert poly.bit_length() == n + 1 and poly & 1
    x = poly_rem(2, poly)
    powers = [x]
    for _ in range(n):
        powers.append(f.square(powers[-1]))
    assert powers[n] == x
    prime_divisors = [p for p in range(2, n + 1)
                      if n % p == 0 and all(p % d for d in range(2, int(p ** .5) + 1))]
    assert all(poly_gcd(poly, powers[n // p] ^ x) == 1 for p in prime_divisors)


def expected_order(n: int) -> int:
    t = [2, -1]
    for _ in range(2, n + 1):
        t.append(-t[-1] - 2 * t[-2])
    return (1 << n) + 1 - t[n]


def rational_fibres(curve):
    f = curve.f
    rows = {}
    for a in range(1 << f.n):
        points = tuple((a, y) for y in range(1 << f.n) if curve.on_curve((a, y)))
        if points:
            assert len(points) == (1 if a == 0 else 2)
            assert {curve.neg(p) for p in points} == set(points)
            rows[a] = points
    assert rows[0] == ((0, 1),)
    assert 1 + sum(len(row) for row in rows.values()) == expected_order(f.n)
    return rows


def polynomial(f, a: int, b: int, c: int) -> int:
    symmetric = f.mul(a, b) ^ f.mul(a, c) ^ f.mul(b, c)
    return f.square(symmetric) ^ f.mul(f.mul(a, b), c) ^ 1


def label(a: int, b: int) -> str:
    if a == 0 and b == 0:
        return "zero_zero_no_root"
    if a == b:
        return "equal_nonzero_linear"
    if 0 in (a, b):
        return "distinct_one_zero_double"
    return "distinct_nonzero_two"


def check_pairs(curve, fibres):
    f = curve.f
    xs = sorted(fibres)
    rows, roots_by_pair = [], {}
    counts = Counter()
    first = {}
    total_roots = 0
    for a in xs:
        for b in xs:
            algebraic = [c for c in range(1 << f.n) if polynomial(f, a, b, c) == 0]
            signed_sums = [curve.add(p, q) for p in fibres[a] for q in fibres[b]]
            geometric = sorted({point[0] for point in signed_sums if point is not None})
            assert algebraic == geometric
            case = label(a, b)
            counts[case] += 1
            total_roots += len(algebraic)
            row = {"n": f.n, "a": a, "b": b, "case": case,
                   "roots": algebraic, "signed_sum_x": geometric,
                   "has_infinity_sum": None in signed_sums}
            first.setdefault(case, {"a": a, "b": b, "roots": algebraic,
                                    "has_infinity_sum": None in signed_sums})
            rows.append(row)
            roots_by_pair[a, b] = tuple(algebraic)
    summary = {"n": f.n, "field_poly": f.poly,
               "group_order": expected_order(f.n),
               "affine_points": sum(map(len, fibres.values())),
               "rational_x_fibres": len(fibres),
               "ordered_rational_pairs": len(rows),
               "polynomial_evaluations": len(rows) * (1 << f.n),
               "distinct_root_occurrences": total_roots,
               "case_counts": dict(sorted(counts.items())),
               "first_cases": first}
    return rows, roots_by_pair, summary


def find_controls(curve, fibres, roots_by_pair):
    f = curve.f
    formal = None
    for a in range(1 << f.n):
        if a in fibres:
            continue
        roots = [c for c in range(1 << f.n) if polynomial(f, a, a, c) == 0]
        assert a and len(roots) == 1
        formal = {"n": f.n, "a": a, "b": a,
                  "formal_roots": roots, "rational_lift_count": 0}
        break
    restricted = None
    nonzero = [a for a in sorted(fibres) if a and len(fibres[a]) == 2]
    for a, b in itertools.combinations(nonzero, 2):
        permitted = curve.add(fibres[a][0], fibres[b][0])
        assert permitted is not None
        roots = roots_by_pair[a, b]
        lost = sorted(set(roots) - {permitted[0]})
        if lost:
            restricted = {"n": f.n, "a": a, "b": b,
                          "permitted_points": [list(fibres[a][0]), list(fibres[b][0])],
                          "permitted_sum_x": permitted[0],
                          "formal_roots": list(roots),
                          "first_excluded_sign_root": lost[0]}
            break
    return formal, restricted


def check_chain(curve, fibres, roots_by_pair, m: int):
    xs = tuple(sorted(fibres))
    algebraic = set()
    for factors in itertools.product(xs, repeat=m):
        paths = [tuple([root]) for root in roots_by_pair[factors[0], factors[1]]]
        for x in factors[2:-1]:
            fresh = []
            for path in paths:
                fresh.extend((*path, next_x) for next_x in roots_by_pair[path[-1], x])
            paths = fresh
        for path in paths:
            for endpoint in roots_by_pair[path[-1], factors[-1]]:
                algebraic.add((factors, path, endpoint))
    from_points = {}
    all_points = tuple(p for fibre in fibres.values() for p in fibre)
    total = o_prefix = terminal_o = all_affine = 0
    for chosen in itertools.product(all_points, repeat=m):
        total += 1
        running = chosen[0]
        prefixes = []
        for next_point in chosen[1:-1]:
            running = curve.add(running, next_point)
            prefixes.append(None if running is None else running[0])
        missing_prefix = any(value is None for value in prefixes)
        o_prefix += missing_prefix
        running = curve.add(running, chosen[-1])
        terminal_o += running is None
        if missing_prefix or running is None:
            continue
        all_affine += 1
        key = (tuple(p[0] for p in chosen), tuple(prefixes), running[0])
        from_points.setdefault(key, set()).add(running)
    assert set(from_points) == algebraic
    assert all(from_points[key] == set(fibres[key[2]]) for key in algebraic)
    rows = [{"n": curve.f.n, "m": m, "factor_x": list(key[0]),
             "prefix_x": list(key[1]), "terminal_x": key[2],
             "realised_full_points": [list(p) for p in sorted(from_points[key])]}
            for key in sorted(algebraic)]
    summary = {"n": curve.f.n, "m": m,
               "rational_x_tuples": len(xs) ** m,
               "signed_point_tuples": total,
               "all_affine_signed_tuples": all_affine,
               "signed_tuples_with_O_prefix": o_prefix,
               "signed_tuples_with_O_terminal": terminal_o,
               "distinct_affine_candidate_paths": len(algebraic),
               "candidate_only_paths": 0,
               "candidate_paths_missing_target_sign": 0}
    return rows, summary


def read_rows(path: Path):
    with gzip.open(path, "rt") as stream:
        return [json.loads(line) for line in stream]


def replay(producer_dir: Path):
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == frozen["verify_sha256"]
    assert sha(PARENT) == frozen["parent_sha256"]
    producer = json.loads((producer_dir / "result.json").read_text())
    assert producer["domain"] == frozen["domain"]
    assert sha(producer_dir / "pair_rows.jsonl.gz") == producer["pair_rows_sha256"]
    assert sha(producer_dir / "chain_rows.jsonl.gz") == producer["chain_rows_sha256"]
    parent = parent_module()
    all_pair_rows, all_chain_rows = [], []
    pair_summaries, chain_summaries = [], []
    first_nonlift = first_restricted = None
    field_ops, curve_ops = Counter(), Counter()
    for n, poly in FIELDS:
        low = [bit for bit in range(n) if poly >> bit & 1]
        field = parent.Field(n, low)
        assert field.poly == poly
        verify_field(field)
        curve = parent.Curve(field)
        by_x = rational_fibres(curve)
        pair_rows, roots, pair_summary = check_pairs(curve, by_x)
        all_pair_rows.extend(pair_rows)
        pair_summaries.append(pair_summary)
        nonlift, restricted = find_controls(curve, by_x, roots)
        if first_nonlift is None and nonlift is not None:
            first_nonlift = nonlift
        if first_restricted is None and restricted is not None:
            first_restricted = restricted
        if n in PANELS:
            chain_rows, chain_summary = check_chain(curve, by_x, roots, PANELS[n])
            all_chain_rows.extend(chain_rows)
            chain_summaries.append(chain_summary)
        field_ops.update(field.operations)
        curve_ops.update(curve.operations)
    assert all_pair_rows == read_rows(producer_dir / "pair_rows.jsonl.gz")
    assert all_chain_rows == read_rows(producer_dir / "chain_rows.jsonl.gz")
    assert producer["fields"] == pair_summaries
    assert producer["chain_panels"] == chain_summaries
    assert producer["first_nonlift_formal_root"] == first_nonlift
    assert producer["first_sign_restriction_loss"] == first_restricted
    return {"decision": "PASS", "domain": frozen["domain"],
            "pair_rows_checked": len(all_pair_rows),
            "chain_rows_checked": len(all_chain_rows),
            "producer_sha256": sha(producer_dir / "result.json"),
            "field_operations": dict(sorted(field_ops.items())),
            "curve_operations": dict(sorted(curve_ops.items()))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"{WALL_CAP}s independent theorem replay cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, WALL_CAP)
    try:
        result = replay(args.producer)
        result.update({"wall_seconds": time.perf_counter() - started,
                       "cpu_seconds": time.process_time() - cpu,
                       "peak_rss_bytes": rss()})
        assert result["wall_seconds"] <= WALL_CAP and result["peak_rss_bytes"] <= RSS_CAP
        save(args.out, result)
    except Exception as error:
        save(args.out, {"decision": "FAIL", "error": repr(error),
                        "wall_seconds": time.perf_counter() - started,
                        "cpu_seconds": time.process_time() - cpu,
                        "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    main()
