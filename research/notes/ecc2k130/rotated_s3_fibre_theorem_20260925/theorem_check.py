#!/usr/bin/env python3
"""Exhaustive finite-fibre and recursive-S3 controls; see PROTOCOL.md."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIELDS = {1: 0x3, 2: 0x7, 3: 0xb, 4: 0x13, 5: 0x25, 6: 0x43, 7: 0x83}
PANELS = ((3, 4), (4, 4), (5, 3))
WALL_CAP = 180
RSS_CAP = 512 * 1024 * 1024


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def remainder(a: int, modulus: int) -> int:
    while a.bit_length() >= modulus.bit_length():
        a ^= modulus << (a.bit_length() - modulus.bit_length())
    return a


def irreducible(poly: int, n: int) -> bool:
    if poly.bit_length() != n + 1 or poly & 1 == 0:
        return False
    for degree in range(1, n // 2 + 1):
        for low in range(1, 1 << degree, 2):
            if remainder(poly, (1 << degree) | low) == 0:
                return False
    return True


class GF:
    def __init__(self, n: int, poly: int):
        self.n, self.poly, self.size = n, poly, 1 << n
        self.ops = Counter()
        self.inverse_cache = {}

    def mul(self, a: int, b: int) -> int:
        self.ops["field_mul"] += 1
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & self.size:
                a ^= self.poly
        return out

    def square(self, a: int) -> int:
        self.ops["field_square"] += 1
        return self.mul(a, a)

    def inv(self, a: int) -> int:
        assert a
        if a not in self.inverse_cache:
            self.ops["field_inverse"] += 1
            exponent, base, result = self.size - 2, a, 1
            while exponent:
                if exponent & 1:
                    result = self.mul(result, base)
                base = self.square(base)
                exponent >>= 1
            assert self.mul(a, result) == 1
            self.inverse_cache[a] = result
        return self.inverse_cache[a]


class E:
    def __init__(self, f: GF):
        self.f = f
        self.ops = Counter()

    def on(self, p) -> bool:
        if p is None:
            return True
        x, y = p
        f = self.f
        return f.square(y) ^ f.mul(x, y) == f.mul(f.square(x), x) ^ 1

    def neg(self, p):
        return None if p is None else (p[0], p[0] ^ p[1])

    def add(self, p, q):
        self.ops["point_add"] += 1
        if p is None:
            return q
        if q is None:
            return p
        x, y = p
        u, v = q
        f = self.f
        if x == u:
            if y ^ v == x:
                return None
            assert y == v and x != 0
            slope = x ^ f.mul(y, f.inv(x))
            nx = f.square(slope) ^ slope
            ny = f.square(x) ^ f.mul(slope ^ 1, nx)
        else:
            slope = f.mul(y ^ v, f.inv(x ^ u))
            nx = f.square(slope) ^ slope ^ x ^ u
            ny = f.mul(slope, x ^ nx) ^ nx ^ y
        return nx, ny


def order(n: int) -> int:
    previous, current = 2, -1
    for _ in range(2, n + 1):
        previous, current = current, -(current + 2 * previous)
    return (1 << n) + 1 - current


def fibres(curve: E):
    f = curve.f
    result = {}
    for x in range(f.size):
        row = tuple((x, y) for y in range(f.size) if curve.on((x, y)))
        if row:
            result[x] = row
            assert (row == ((0, 1),)) if x == 0 else (len(row) == 2)
            assert {curve.neg(p) for p in row} == set(row)
    assert sum(map(len, result.values())) + 1 == order(f.n)
    assert result[0] == ((0, 1),)
    return result


def s3(f: GF, a: int, b: int, c: int) -> int:
    product = f.mul(a, b)
    return f.mul(f.square(a ^ b), f.square(c)) ^ f.mul(product, c) ^ f.square(product) ^ 1


def pair_case(a: int, b: int) -> str:
    if a == b == 0:
        return "zero_zero_no_root"
    if a == b:
        return "equal_nonzero_linear"
    if a * b == 0:
        return "distinct_one_zero_double"
    return "distinct_nonzero_two"


def json_line(row) -> bytes:
    return (json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode()


def write_rows(path: Path, rows):
    with path.open("wb") as file:
        with gzip.GzipFile(filename="", mode="wb", fileobj=file, mtime=0) as zipped:
            for row in rows:
                zipped.write(json_line(row))


def pair_control(curve: E, fibres_by_x: dict):
    f = curve.f
    rows = []
    cases = Counter()
    root_count = 0
    first_cases = {}
    pair_root_cache = {}
    for a, b in itertools.product(sorted(fibres_by_x), repeat=2):
        roots = tuple(c for c in range(f.size) if s3(f, a, b, c) == 0)
        finite = {curve.add(p, q) for p in fibres_by_x[a] for q in fibres_by_x[b]}
        geometry = tuple(sorted(p[0] for p in finite if p is not None))
        # Multiple point sums may share x; compare distinct x values.
        geometry = tuple(sorted(set(geometry)))
        assert roots == geometry, (f.n, a, b, roots, geometry)
        label = pair_case(a, b)
        cases[label] += 1
        root_count += len(roots)
        if label not in first_cases:
            first_cases[label] = {"a": a, "b": b, "roots": list(roots),
                                  "has_infinity_sum": None in finite}
        pair_root_cache[a, b] = roots
        rows.append({"n": f.n, "a": a, "b": b, "case": label,
                     "roots": list(roots), "signed_sum_x": list(geometry),
                     "has_infinity_sum": None in finite})
    return rows, {"n": f.n, "field_poly": f.poly, "group_order": order(f.n),
                  "affine_points": sum(map(len, fibres_by_x.values())),
                  "rational_x_fibres": len(fibres_by_x),
                  "ordered_rational_pairs": len(rows),
                  "polynomial_evaluations": len(rows) * f.size,
                  "distinct_root_occurrences": root_count,
                  "case_counts": dict(sorted(cases.items())),
                  "first_cases": first_cases}, pair_root_cache


def negative_controls(curve: E, fibres_by_x: dict, pair_cache: dict):
    f = curve.f
    nonlift = None
    for x in range(f.size):
        if x in fibres_by_x:
            continue
        roots = [c for c in range(f.size) if s3(f, x, x, c) == 0]
        assert x != 0 and len(roots) == 1
        nonlift = {"n": f.n, "a": x, "b": x, "formal_roots": roots,
                   "rational_lift_count": 0}
        break
    restricted = None
    xs = sorted(x for x, row in fibres_by_x.items() if x != 0 and len(row) == 2)
    for a in xs:
        for b in xs:
            if a >= b:
                continue
            only = curve.add(fibres_by_x[a][0], fibres_by_x[b][0])
            assert only is not None
            roots = pair_cache[a, b]
            lost = sorted(set(roots) - {only[0]})
            if lost:
                restricted = {"n": f.n, "a": a, "b": b,
                              "permitted_points": [list(fibres_by_x[a][0]),
                                                   list(fibres_by_x[b][0])],
                              "permitted_sum_x": only[0],
                              "formal_roots": list(roots),
                              "first_excluded_sign_root": lost[0]}
                break
        if restricted is not None:
            break
    return nonlift, restricted


def chain_control(curve: E, fibres_by_x: dict, m: int, root_cache: dict):
    f = curve.f
    xs = tuple(sorted(fibres_by_x))
    algebraic = set()
    for x_tuple in itertools.product(xs, repeat=m):
        paths = [(c,) for c in root_cache[x_tuple[0], x_tuple[1]]]
        for factor_x in x_tuple[2:-1]:
            paths = [path + (c,) for path in paths
                     for c in root_cache[path[-1], factor_x]]
        for path in paths:
            for terminal_x in root_cache[path[-1], x_tuple[-1]]:
                algebraic.add((x_tuple, path, terminal_x))
    geometry = {}
    total, o_prefix, terminal_o, affine = 0, 0, 0, 0
    points = tuple(p for row in fibres_by_x.values() for p in row)
    for tuple_points in itertools.product(points, repeat=m):
        total += 1
        acc = tuple_points[0]
        path = []
        for factor in tuple_points[1:-1]:
            acc = curve.add(acc, factor)
            path.append(None if acc is None else acc[0])
        has_o = None in path
        if has_o:
            o_prefix += 1
        acc = curve.add(acc, tuple_points[-1])
        if acc is None:
            terminal_o += 1
        if has_o or acc is None:
            continue
        affine += 1
        key = (tuple(p[0] for p in tuple_points), tuple(path), acc[0])
        geometry.setdefault(key, set()).add(acc)
    assert set(geometry) == algebraic, (f.n, m, len(geometry), len(algebraic))
    assert all(points == set(fibres_by_x[key[2]]) for key, points in geometry.items())
    rows = [{"n": f.n, "m": m, "factor_x": list(key[0]),
             "prefix_x": list(key[1]), "terminal_x": key[2],
             "realised_full_points": [list(p) for p in sorted(geometry[key])]}
            for key in sorted(algebraic)]
    return rows, {"n": f.n, "m": m, "rational_x_tuples": len(xs) ** m,
                  "signed_point_tuples": total,
                  "all_affine_signed_tuples": affine,
                  "signed_tuples_with_O_prefix": o_prefix,
                  "signed_tuples_with_O_terminal": terminal_o,
                  "distinct_affine_candidate_paths": len(algebraic),
                  "candidate_only_paths": 0,
                  "candidate_paths_missing_target_sign": 0}


def run(out: Path):
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_sig, _frame):
        raise TimeoutError(f"{WALL_CAP}s theorem producer cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, WALL_CAP)
    out.mkdir(parents=True, exist_ok=False)
    try:
        freeze = json.loads((HERE / "FROZEN.json").read_text())
        assert sha(Path(__file__)) == freeze["producer_sha256"]
        assert sha(HERE / "PROTOCOL.md") == freeze["protocol_sha256"]
        assert sha(HERE / "PROOF.md") == freeze["proof_sha256"]
        pair_rows, chain_rows = [], []
        pair_summaries, chain_summaries = [], []
        first_nonlift = first_restricted = None
        field_ops, curve_ops = Counter(), Counter()
        for n, poly in FIELDS.items():
            assert irreducible(poly, n)
            f = GF(n, poly)
            curve = E(f)
            by_x = fibres(curve)
            rows, summary, cache = pair_control(curve, by_x)
            pair_rows.extend(rows)
            pair_summaries.append(summary)
            nonlift, restricted = negative_controls(curve, by_x, cache)
            if first_nonlift is None and nonlift is not None:
                first_nonlift = nonlift
            if first_restricted is None and restricted is not None:
                first_restricted = restricted
            for panel_n, m in PANELS:
                if panel_n != n:
                    continue
                rows, chain_summary = chain_control(curve, by_x, m, cache)
                chain_rows.extend(rows)
                chain_summaries.append(chain_summary)
            field_ops.update(f.ops)
            curve_ops.update(curve.ops)
        assert first_nonlift is not None and first_restricted is not None
        write_rows(out / "pair_rows.jsonl.gz", pair_rows)
        write_rows(out / "chain_rows.jsonl.gz", chain_rows)
        result = {"domain": freeze["domain"], "fields": pair_summaries,
                  "chain_panels": chain_summaries,
                  "first_nonlift_formal_root": first_nonlift,
                  "first_sign_restriction_loss": first_restricted,
                  "pair_rows_sha256": sha(out / "pair_rows.jsonl.gz"),
                  "chain_rows_sha256": sha(out / "chain_rows.jsonl.gz"),
                  "field_operations": dict(sorted(field_ops.items())),
                  "curve_operations": dict(sorted(curve_ops.items())),
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu,
                  "peak_rss_bytes": rss()}
        assert result["wall_seconds"] <= WALL_CAP and result["peak_rss_bytes"] <= RSS_CAP
        save(out / "result.json", result)
    except Exception as error:
        save(out / "failure.json", {"error": repr(error),
                                    "wall_seconds": time.perf_counter() - started,
                                    "cpu_seconds": time.process_time() - cpu,
                                    "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)


if __name__ == "__main__":
    main()
