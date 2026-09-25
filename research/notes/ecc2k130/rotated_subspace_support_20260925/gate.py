#!/usr/bin/env python3
"""Frozen rotated normal-basis support gate; see PROTOCOL.md before running."""
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
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOMAIN = "ECC2K130-ROTATED-SUBSPACE-20260925-v1"
Q131 = 680564733841876926932320129493409985129
LAMBDA131 = 196511074115861092422032515080945363956
MODELS = {
    13: {"low": [0, 1, 3, 4], "order": 2003},
    131: {"low": [0, 1, 2, 13], "order": Q131},
}
GRID = [(5, 24), (5, 25), (5, 26), (6, 20), (6, 21)]
Point = tuple[int, int] | None
ACTIVE_STATE: dict = {}


class Field:
    def __init__(self, n: int, low: list[int]):
        self.n = n
        self.low = tuple(low)
        self.poly = (1 << n) | sum(1 << i for i in low)
        self.mask = (1 << n) - 1
        self.operations = Counter()

    def _reduce(self, v: int) -> int:
        while v.bit_length() > self.n:
            v ^= self.poly << (v.bit_length() - self.n - 1)
        return v

    def _mul(self, a: int, b: int) -> int:
        out = 0
        while b:
            if b & 1:
                out ^= a
            a <<= 1
            b >>= 1
        return self._reduce(out)

    def mul(self, a: int, b: int) -> int:
        self.operations["field_mul"] += 1
        return self._mul(a, b)

    def square(self, a: int) -> int:
        self.operations["field_square"] += 1
        return self._mul(a, a)

    def inverse(self, a: int) -> int:
        assert 0 < a <= self.mask
        self.operations["field_inverse"] += 1
        u, v, g1, g2 = a, self.poly, 1, 0
        while u != 1:
            if u == 0:
                raise ArithmeticError("noninvertible element")
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v = v, u
                g1, g2 = g2, g1
                shift = -shift
            u ^= v << shift
            g1 ^= g2 << shift
        return self._reduce(g1)

    def trace(self, a: int) -> int:
        acc = 0
        v = a
        for _ in range(self.n):
            acc ^= v
            v = self.square(v)
        assert acc in (0, 1) and v == a
        return acc

    def rabin_prime_degree(self) -> None:
        assert self.poly & 1
        assert bin(self.poly).count("1") & 1  # no X+1 factor
        a = 2
        for _ in range(self.n):
            a = self.square(a)
        assert a == 2
        assert is_prime_by_trial(self.n)


class Curve:
    def __init__(self, field: Field):
        self.f = field
        self.operations = Counter()

    def neg(self, p: Point) -> Point:
        return None if p is None else (p[0], p[0] ^ p[1])

    def tau(self, p: Point) -> Point:
        return None if p is None else (self.f.square(p[0]), self.f.square(p[1]))

    def on_curve(self, p: Point) -> bool:
        if p is None:
            return True
        x, y = p
        f = self.f
        return f.square(y) ^ f.mul(x, y) == f.mul(f.square(x), x) ^ 1

    def add(self, p: Point, q: Point) -> Point:
        self.operations["point_add"] += 1
        if p is None:
            return q
        if q is None:
            return p
        x1, y1 = p
        x2, y2 = q
        f = self.f
        if x1 == x2:
            if y1 ^ y2 == x1:
                return None
            assert y1 == y2 and x1 != 0
            lam = x1 ^ f.mul(y1, f.inverse(x1))
            x3 = f.square(lam) ^ lam
            y3 = f.square(x1) ^ f.mul(lam ^ 1, x3)
            return (x3, y3)
        lam = f.mul(y1 ^ y2, f.inverse(x1 ^ x2))
        x3 = f.square(lam) ^ lam ^ x1 ^ x2
        y3 = f.mul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)

    def scalar(self, p: Point, k: int) -> Point:
        self.operations["scalar_call"] += 1
        assert k >= 0
        acc = None
        while k:
            if k & 1:
                acc = self.add(acc, p)
            p = self.add(p, p)
            k >>= 1
        return acc


def source_group_order(n: int) -> int:
    # E(F_2)={O,(0,1),(1,0),(1,1)}, hence t_1=2+1-4=-1.
    previous, current = 2, -1
    if n == 1:
        return 4
    for _ in range(2, n + 1):
        previous, current = current, -current - 2 * previous
    return (1 << n) + 1 - current


def is_prime_by_trial(n: int) -> bool:
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def rank(vectors: list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in vectors:
        v = value
        while v:
            bit = v.bit_length() - 1
            if bit in pivots:
                v ^= pivots[bit]
            else:
                pivots[bit] = v
                break
    return len(pivots)


def normal_conjugates(f: Field, beta: int) -> list[int]:
    out = []
    cur = beta
    for _ in range(f.n):
        out.append(cur)
        cur = f.square(cur)
    assert cur == beta and rank(out) == f.n
    assert f.trace(beta) == 1
    return out


def first_alternate_normal_beta(f: Field, beta_a: int) -> int:
    source_orbit = set(normal_conjugates(f, beta_a))
    for candidate in range(beta_a + 1, 1 << f.n):
        if candidate in source_orbit:
            continue
        conjugates = []
        value = candidate
        for _ in range(f.n):
            conjugates.append(value)
            value = f.square(value)
        if rank(conjugates) == f.n:
            return candidate
    raise AssertionError("no second normal generator")


def subspace_basis(conjugates: list[int], m: int, d: int) -> list[list[int]]:
    n = len(conjugates)
    assert m * d <= n
    bases = [[conjugates[m * j + i] for j in range(d)] for i in range(m)]
    assert all(rank(b) == d for b in bases)
    assert rank([v for b in bases for v in b]) == m * d
    return bases


def x_from_mask(basis: list[int], mask: int) -> int:
    value = 0
    for j, v in enumerate(basis):
        if (mask >> j) & 1:
            value ^= v
    return value


def all_x(basis: list[int]) -> list[int]:
    return [x_from_mask(basis, mask) for mask in range(1 << len(basis))]


def artin_schreier_roots(f: Field) -> dict[int, list[int]]:
    assert f.n == 13
    roots: dict[int, list[int]] = {}
    for z in range(1 << f.n):
        rhs = f.square(z) ^ z
        roots.setdefault(rhs, []).append(z)
    assert len(roots) == 1 << (f.n - 1)
    assert all(len(v) == 2 for v in roots.values())
    return roots


def lift_x(f: Field, roots: dict[int, list[int]], x: int) -> list[Point]:
    if x == 0:
        return [(0, 1)]
    rhs = x ^ f.square(f.inverse(x))
    return [(x, f.mul(x, z)) for z in roots.get(rhs, [])]


def whole_curve(curve: Curve) -> tuple[list[tuple[int, int]], dict[int, list[tuple[int, int]]]]:
    f = curve.f
    roots = artin_schreier_roots(f)
    by_x: dict[int, list[tuple[int, int]]] = {}
    all_points = []
    for x in range(1 << f.n):
        points = sorted(lift_x(f, roots, x))
        assert all(curve.on_curve(p) for p in points)
        by_x[x] = points
        all_points.extend(points)
    assert len(all_points) + 1 == source_group_order(13) == 8012
    return all_points, by_x


def subgroup_generator(curve: Curve, points: list[tuple[int, int]]) -> tuple[int, int]:
    for p in points:
        h = curve.scalar(p, 4)
        if h is not None:
            assert curve.scalar(h, 2003) is None
            return h
    raise AssertionError("no subgroup generator")


def four_torsion(curve: Curve) -> list[Point]:
    torsion = [None, (0, 1), (1, 0), (1, 1)]
    assert all(curve.on_curve(t) for t in torsion)
    assert curve.scalar(torsion[1], 2) is None
    assert curve.scalar(torsion[2], 2) == torsion[1]
    assert curve.scalar(torsion[3], 2) == torsion[1]
    return torsion


def full_histogram(curve: Curve, bases: list[list[tuple[int, int]]]):
    counts: dict[Point, int] = {None: 1}
    witnesses: dict[Point, tuple[int, ...]] = {None: ()}
    for base in bases:
        new_counts: dict[Point, int] = {}
        new_witnesses: dict[Point, tuple[int, ...]] = {}
        for partial, multiplicity in counts.items():
            for index, point in enumerate(base):
                target = curve.add(partial, point)
                new_counts[target] = new_counts.get(target, 0) + multiplicity
                if target not in new_witnesses:
                    new_witnesses[target] = witnesses[partial] + (index,)
        counts, witnesses = new_counts, new_witnesses
    assert sum(counts.values()) == math.prod(len(b) for b in bases)
    return counts, witnesses


def pjson(p: Point):
    return None if p is None else list(p)


def tau_orbit(curve: Curve, p: Point) -> int:
    v = curve.tau(p)
    for length in range(1, curve.f.n + 1):
        if v == p:
            assert curve.f.n % length == 0
            return length
        v = curve.tau(v)
    raise AssertionError("Frobenius orbit did not close")


def hard_deadline(seconds: int, label: str, field: Field, curve: Curve | None = None) -> None:
    ACTIVE_STATE.update({"label": label, "field": field, "curve": curve,
                         "started_perf": time.perf_counter()})
    def expire(_signum, _frame):
        raise TimeoutError(f"{label}: {seconds}s hard wall deadline")
    signal.signal(signal.SIGALRM, expire)
    signal.setitimer(signal.ITIMER_REAL, seconds)


def clear_deadline() -> None:
    signal.setitimer(signal.ITIMER_REAL, 0)
    ACTIVE_STATE.clear()


def peak_rss_bytes() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def phase_delta(f: Field, c: Curve, before_f: Counter, before_c: Counter) -> dict[str, int]:
    values = {key: f.operations[key] - before_f[key] for key in sorted(f.operations)}
    values.update({key: c.operations[key] - before_c[key] for key in sorted(c.operations)})
    return values


def save_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n")


def toy_arm(curve: Curve, by_x, beta: int, m: int, repeated: bool,
            subgroup: list[Point], torsion: list[Point], features: list[dict]):
    f = curve.f
    before_field, before_curve = f.operations.copy(), curve.operations.copy()
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    snapshots = [("start", wall_start, cpu_start, before_field, before_curve)]
    conjugates = normal_conjugates(f, beta)
    spaces = subspace_basis(conjugates, m, 2)
    factors = []
    sign_pairs = []
    for basis in spaces:
        xs = all_x(basis)
        assert len(set(xs)) == 4
        factor = sorted(p for x in xs for p in by_x[x])
        assert all(curve.on_curve(p) for p in factor)
        factors.append(factor)
        sign_pairs.append(sum(len(by_x[x]) == 2 for x in xs))
    for i in range(1, m):
        assert factors[i] == sorted(curve.tau(p) for p in factors[i - 1])
    actual = [factors[0]] * m if repeated else factors
    snapshots.append(("factor_setup", time.perf_counter(), time.process_time(),
                      f.operations.copy(), curve.operations.copy()))
    histogram, witnesses = full_histogram(curve, actual)
    snapshots.append(("histogram", time.perf_counter(), time.process_time(),
                      f.operations.copy(), curve.operations.copy()))
    projected_histogram: dict[Point, int] = {}
    for full_sum, multiplicity in histogram.items():
        projected = curve.scalar(full_sum, 4)
        projected_histogram[projected] = projected_histogram.get(projected, 0) + multiplicity
        witness = witnesses[full_sum]
        assert len(witness) == m
        replay = None
        for factor, index in zip(actual, witness):
            replay = curve.add(replay, factor[index])
        assert replay == full_sum
    assert sum(projected_histogram.values()) == sum(histogram.values())
    all_projected = [set(curve.scalar(p, 4) for p in factor) for factor in actual]
    for factor in factors:
        for p in factor:
            assert curve.scalar(curve.tau(p), 4) == curve.tau(curve.scalar(p, 4))
    snapshots.append(("projection_and_witness_self_check", time.perf_counter(),
                      time.process_time(), f.operations.copy(), curve.operations.copy()))
    rows = []
    coset_rhs_charges = 0
    for k, q in enumerate(subgroup):
        shifted = [curve.add(q, t) for t in torsion]
        coset_rhs_charges += 4
        assert len(set(shifted)) == 4
        multiplicities = [histogram.get(p, 0) for p in shifted]
        indexed_witnesses = [list(witnesses[p]) if p in witnesses else None for p in shifted]
        projection = curve.scalar(q, 4)
        projected_count = projected_histogram.get(projection, 0)
        assert projected_count == sum(multiplicities), (k, projected_count, multiplicities)
        rows.append({
            "k": k, "point": pjson(q), "projected_point": pjson(projection),
            **features[k],
            "coset_multiplicities": multiplicities,
            "coset_witness_indices": indexed_witnesses,
            "projected_multiplicity": projected_count,
            "raw_hit": multiplicities[0] > 0,
            "projected_hit": projected_count > 0,
        })
    assert coset_rhs_charges == 4 * len(subgroup)
    snapshots.append(("all_four_rhs", time.perf_counter(), time.process_time(),
                      f.operations.copy(), curve.operations.copy()))
    stages = []
    for prev, cur in zip(snapshots, snapshots[1:]):
        name, wall, cpu, field_ops, curve_ops = cur
        stages.append({"name": name, "wall_seconds": wall - prev[1],
                       "cpu_seconds": cpu - prev[2],
                       "operations": {key: field_ops[key] - prev[3][key]
                                      for key in sorted(set(field_ops) | set(prev[3]))} | {
                           key: curve_ops[key] - prev[4][key]
                           for key in sorted(set(curve_ops) | set(prev[4]))}})
    summary = {
        "n": 13, "m": m, "d": 2, "beta": beta,
        "policy": "repeated" if repeated else "rotated",
        "factor_sizes": list(map(len, actual)),
        "sign_pair_x_counts": sign_pairs if not repeated else [sign_pairs[0]] * m,
        "physical_point_choices": sum(map(len, actual)),
        "projected_unique_each_factor": [len(s) for s in all_projected],
        "projected_duplicates_each_factor": [len(b) - len(s) for b, s in zip(actual, all_projected)],
        "compressed_log_columns_upper": len(all_projected[0]),
        "tuple_count": sum(histogram.values()),
        "distinct_full_sums": len(histogram),
        "distinct_projected_sums": len(projected_histogram),
        "tuple_collisions": sum(histogram.values()) - len(histogram),
        "projected_tuple_collisions": sum(histogram.values()) - len(projected_histogram),
        "raw_H_hits": sum(r["raw_hit"] for r in rows),
        "projected_H_hits": sum(r["projected_hit"] for r in rows),
        "raw_nonzero_H_hits": sum(r["raw_hit"] for r in rows[1:]),
        "projected_nonzero_H_hits": sum(r["projected_hit"] for r in rows[1:]),
        "coset_hits": [sum(r["coset_multiplicities"][j] > 0 for r in rows) for j in range(4)],
        "coset_rhs_charges": coset_rhs_charges,
        "stages": stages,
        "total_wall_seconds": snapshots[-1][1] - wall_start,
        "total_cpu_seconds": snapshots[-1][2] - cpu_start,
        "peak_rss_bytes": peak_rss_bytes(),
        "total_operations": phase_delta(f, curve, before_field, before_curve),
    }
    assert summary["total_wall_seconds"] <= 180, (beta, m, repeated, summary["total_wall_seconds"])
    assert summary["peak_rss_bytes"] <= 512 * 1024 * 1024
    return summary, rows, actual


def run_toy(out: Path) -> None:
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    f = Field(13, MODELS[13]["low"])
    f.rabin_prime_degree()
    c = Curve(f)
    points, by_x = whole_curve(c)
    torsion = four_torsion(c)
    h = subgroup_generator(c, points)
    assert c.tau(h) == c.scalar(h, 89)
    assert c.scalar(h, 2003) is None
    assert 2003 % 2 and is_prime_by_trial(2003)
    assert first_alternate_normal_beta(f, 3) == 7
    subgroup = []
    q = None
    for k in range(2003):
        subgroup.append(q)
        q = c.add(q, h)
    assert q is None and len(set(subgroup)) == 2003
    features = []
    for q in subgroup:
        torsion_projection = c.scalar(q, 2003)
        assert torsion_projection in torsion
        cofactor_class = torsion.index(torsion_projection)
        features.append({
            "trace_x": None if q is None else f.trace(q[0]),
            "tau_orbit_size": tau_orbit(c, q),
            "cofactor_class": cofactor_class,
            "torsion_projection": pjson(torsion_projection),
        })
    assert all(row["cofactor_class"] == 0 for row in features)
    assert all(row["tau_orbit_size"] == 13 and row["trace_x"] == 0 for row in features[1:])
    global_setup = {"wall_seconds": time.perf_counter() - wall_start,
                    "cpu_seconds": time.process_time() - cpu_start,
                    "operations": dict(f.operations) | dict(c.operations),
                    "peak_rss_bytes": peak_rss_bytes()}
    assert global_setup["peak_rss_bytes"] <= 512 * 1024 * 1024
    summary_rows = []
    for beta in (3, 7):
        for m in (5, 6):
            for repeated in (False, True):
                hard_deadline(180, f"n13 beta={beta} m={m} repeated={repeated}", f, c)
                summary, rows, factors = toy_arm(c, by_x, beta, m, repeated,
                                                subgroup, torsion, features)
                stem = f"n13-b{beta}-m{m}-{'repeated' if repeated else 'rotated'}"
                save_json(out / (stem + "-summary.json"), summary)
                with (out / (stem + "-targets.jsonl")).open("w") as stream:
                    for row in rows:
                        stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
                save_json(out / (stem + "-factors.json"), [[pjson(p) for p in factor] for factor in factors])
                summary_rows.append(summary)
                assert peak_rss_bytes() <= 512 * 1024 * 1024
                clear_deadline()
    save_json(out / "toy_summary.json", {"generator": pjson(h), "torsion": list(map(pjson, torsion)),
                                          "curve_order": 8012, "global_setup": global_setup,
                                          "rows": summary_rows})


def half_trace(f: Field, a: int) -> int:
    assert f.n & 1 and f.trace(a) == 0
    z, term = 0, a
    for _ in range((f.n + 1) // 2):
        z ^= term
        term = f.square(f.square(term))
    assert f.square(z) ^ z == a
    return z


def n131_projected_lambda_point(curve: Curve, tmask: int) -> dict:
    f = curve.f
    assert f.n == 131
    lam = LAMBDA131
    assert (lam * lam + lam + 2) % Q131 == 0
    assert pow(lam, 131, Q131) == 1
    for x in range(2, 1 << 12):
        rhs = x ^ f.square(f.inverse(x))
        if trace_fast(tmask, rhs) != 0:
            continue
        z = half_trace(f, rhs)
        p = (x, f.mul(x, z))
        assert curve.on_curve(p)
        h = curve.scalar(p, 4)
        if h is None:
            continue
        assert curve.scalar(h, Q131) is None
        assert curve.tau(h) == curve.scalar(h, lam)
        assert curve.scalar(curve.tau(p), 4) == curve.tau(h)
        return {"x": x, "point": pjson(p), "projected_H": pjson(h),
                "lambda": lam, "group_order": 4 * Q131,
                "checks": ["point_on_curve", "nonzero_[4]projection",
                           "[q]H=O", "tau(H)=[lambda]H",
                           "[4]tau(P)=tau([4]P)"]}
    raise AssertionError("no deterministic projected n131 point")


def trace_mask(f: Field) -> int:
    mask = 0
    for bit in range(f.n):
        mask |= f.trace(1 << bit) << bit
    return mask


def trace_fast(mask: int, x: int) -> int:
    return (mask & x).bit_count() & 1


def deterministic_masks(kind: str, m: int, d: int, count: int,
                        excluded: set[int] | None = None) -> list[int]:
    assert kind in ("cov", "density")
    excluded = set() if excluded is None else set(excluded)
    result = []
    used = excluded.copy()
    if kind == "cov":
        for mask in (0, 1, 1 << (d - 1), (1 << d) - 1):
            assert mask not in used
            used.add(mask)
            result.append(mask)
    counter = 0
    while len(result) < count:
        source = f"{DOMAIN}/{kind}/{m}/{d}/{counter}".encode("ascii")
        candidate = int.from_bytes(hashlib.sha256(source).digest(), "big") % (1 << d)
        counter += 1
        if candidate not in used and not (kind == "density" and candidate == 0):
            used.add(candidate)
            result.append(candidate)
    return result


def freeze_inputs(out: Path) -> None:
    f = Field(131, MODELS[131]["low"])
    f.rabin_prime_degree()
    conjugates = normal_conjugates(f, 3)
    payload = {"domain": DOMAIN, "model": 131, "polynomial": f.poly,
               "beta": 3, "generator": "SHA256 big-endian reduced modulo 2^d; skip duplicate/zero/excluded; counter starts zero",
               "cells": []}
    for m, d in GRID:
        bases = subspace_basis(conjugates, m, d)
        for i in range(m - 1):
            assert [f.square(x) for x in bases[i]] == bases[i + 1]
        cov = deterministic_masks("cov", m, d, 256)
        density = deterministic_masks("density", m, d, 1 << 14, set(cov))
        assert len(set(cov) | set(density)) == len(cov) + len(density)
        prefix = f"n131-m{m}-d{d}"
        entries = {}
        for kind, masks in (("cov", cov), ("density", density)):
            path = out / f"{prefix}-{kind}-masks.json"
            save_json(path, masks)
            entries[kind] = {"path": path.name,
                             "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                             "count": len(masks)}
        payload["cells"].append({"m": m, "d": d, **entries})
    save_json(out / "input_manifest.json", payload)


def wilson_bounds(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = success / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    radius = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return center - radius, center + radius


def load_frozen_mask(path: Path, expected_sha: str, count: int) -> list[int]:
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == expected_sha
    masks = json.loads(raw)
    assert len(masks) == count and len(set(masks)) == count
    return masks


def run_density(out: Path, inputs: Path) -> None:
    global_wall_start, global_cpu_start = time.perf_counter(), time.process_time()
    manifest = json.loads((inputs / "input_manifest.json").read_text())
    assert manifest["domain"] == DOMAIN and manifest["model"] == 131
    f = Field(131, MODELS[131]["low"])
    f.rabin_prime_degree()
    assert manifest["polynomial"] == f.poly and manifest["beta"] == 3
    assert source_group_order(131) == 4 * Q131
    conjugates = normal_conjugates(f, 3)
    tmask = trace_mask(f)
    assert trace_fast(tmask, 3) == 1
    c = Curve(f)
    lambda_point = n131_projected_lambda_point(c, tmask)
    global_setup = {"wall_seconds": time.perf_counter() - global_wall_start,
                    "cpu_seconds": time.process_time() - global_cpu_start,
                    "field_operations": dict(f.operations),
                    "curve_operations": dict(c.operations),
                    "peak_rss_bytes": peak_rss_bytes()}
    assert global_setup["peak_rss_bytes"] <= 512 * 1024 * 1024
    summaries = []
    for cell in manifest["cells"]:
        m, d = cell["m"], cell["d"]
        assert (m, d) in GRID
        hard_deadline(120, f"n131 m={m} d={d}", f)
        wall_start, cpu_start = time.perf_counter(), time.process_time()
        bf = f.operations.copy()
        bases = subspace_basis(conjugates, m, d)
        for i in range(m - 1):
            assert [f.square(x) for x in bases[i]] == bases[i + 1]
        cov = load_frozen_mask(inputs / cell["cov"]["path"], cell["cov"]["sha256"], 256)
        density = load_frozen_mask(inputs / cell["density"]["path"], cell["density"]["sha256"], 1 << 14)
        assert cov == deterministic_masks("cov", m, d, 256)
        assert density == deterministic_masks("density", m, d, 1 << 14, set(cov))
        cov_rows = []
        for mask in cov:
            xs = [x_from_mask(basis, mask) for basis in bases]
            for i in range(1, m):
                assert f.square(xs[i - 1]) == xs[i]
            parity = mask.bit_count() & 1
            assert all(trace_fast(tmask, x) == parity for x in xs)
            lifts = []
            for x in xs:
                solvable = True if x == 0 else trace_fast(tmask, x ^ f.square(f.inverse(x))) == 0
                lifts.append(1 if x == 0 else 2 * int(solvable))
            assert len(set(lifts)) == 1
            cov_rows.append({"mask": mask, "x": xs, "lifts": lifts})
        density_rows = []
        for index, mask in enumerate(density):
            x = x_from_mask(bases[0], mask)
            assert x != 0 and trace_fast(tmask, x) == (mask.bit_count() & 1)
            rhs = x ^ f.square(f.inverse(x))
            solvable = trace_fast(tmask, rhs) == 0
            if index < 32 or index % 1024 == 0:
                assert f.trace(x) == (mask.bit_count() & 1)
                assert f.trace(rhs) == int(not solvable)
            density_rows.append({"mask": mask, "x": x, "solvable": solvable})
        success = sum(row["solvable"] for row in density_rows)
        lower, upper = wilson_bounds(success, len(density_rows))
        # These are descriptive extrapolations, never exact full-space bounds.
        estimate_points = 1 + 2 * ((1 << d) - 1) * success / len(density_rows)
        upper_points = 1 + 2 * ((1 << d) - 1) * upper
        ideal_points = 1 << d
        raw_ratio = min(1, ideal_points**m / (4 * Q131))
        projected_ratio = min(1, ideal_points**m / Q131)
        summary = {"m": m, "d": d, "md": m * d, "beta": 3,
                   "covariance_masks": len(cov), "density_masks": len(density),
                   "density_solvable": success, "density_rate": success / len(density),
                   "wilson_95_model_interval": [lower, upper],
                   "conditional_point_size_estimate": estimate_points,
                   "conditional_point_size_upper": upper_points,
                   "ideal_raw_tuple_ratio_ceiling": raw_ratio,
                   "ideal_projected_tuple_ratio_ceiling": projected_ratio,
                   "conditional_sampled_raw_ratio": min(1, estimate_points**m / (4 * Q131)),
                   "conditional_sampled_projected_ratio": min(1, estimate_points**m / Q131),
                   "conditional_wilson_upper_projected_ratio": min(1, upper_points**m / Q131),
                   "physical_ideal_point_choices": m * ideal_points,
                   "compressed_ideal_point_choices_proxy": ideal_points,
                   "compressed_rigorous_point_choices_upper": 2 * ideal_points - 1,
                   "wall_seconds": time.perf_counter() - wall_start,
                   "cpu_seconds": time.process_time() - cpu_start,
                   "field_operations": {k: f.operations[k] - bf[k] for k in sorted(f.operations)},
                   "peak_rss_bytes": peak_rss_bytes()}
        prefix = f"n131-m{m}-d{d}"
        save_json(out / f"{prefix}-summary.json", summary)
        with (out / f"{prefix}-covariance.jsonl").open("w") as stream:
            for row in cov_rows:
                stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        with (out / f"{prefix}-density.jsonl").open("w") as stream:
            for row in density_rows:
                stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        summaries.append(summary)
        assert summary["wall_seconds"] <= 120, (m, d, summary["wall_seconds"])
        assert peak_rss_bytes() <= 512 * 1024 * 1024
        clear_deadline()
    save_json(out / "density_summary.json", {"n": 131, "rows": summaries,
                                            "trace_mask": tmask,
                                            "projected_lambda_point": lambda_point,
                                            "global_setup": global_setup})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["freeze-inputs", "toy", "density"])
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--inputs", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    try:
        if args.mode == "freeze-inputs":
            freeze_inputs(args.out)
        elif args.mode == "toy":
            run_toy(args.out)
        else:
            assert args.inputs is not None
            run_density(args.out, args.inputs)
    except Exception as error:
        state = ACTIVE_STATE.copy()
        field = state.get("field")
        curve = state.get("curve")
        save_json(args.out / "producer_failure.json", {
            "error": repr(error), "active_label": state.get("label"),
            "active_wall_seconds": time.perf_counter() - state["started_perf"] if state else None,
            "field_operations_so_far": dict(field.operations) if field else None,
            "curve_operations_so_far": dict(curve.operations) if curve else None,
            "peak_rss_bytes": peak_rss_bytes(),
        })
        raise


if __name__ == "__main__":
    main()
