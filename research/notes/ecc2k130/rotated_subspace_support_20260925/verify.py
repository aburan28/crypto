#!/usr/bin/env python3
"""Independent arithmetic and exhaustive replay of the rotated-space gate."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import resource
import sys
import time
from collections import Counter
from functools import lru_cache
from pathlib import Path

DOMAIN = "ECC2K130-ROTATED-SUBSPACE-20260925-v1"
P13 = (1 << 13) | (1 << 4) | (1 << 3) | (1 << 1) | 1
P131 = (1 << 131) | (1 << 13) | (1 << 2) | (1 << 1) | 1
Q131 = 680564733841876926932320129493409985129
LAMBDA131 = 196511074115861092422032515080945363956
GRID = [(5, 24), (5, 25), (5, 26), (6, 20), (6, 21)]


class GF:
    def __init__(self, n: int, poly: int):
        self.n, self.poly = n, poly
        self.mask = (1 << n) - 1
        self.ops = Counter()

    def mul(self, a: int, b: int) -> int:
        self.ops["mul"] += 1
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a >> self.n:
                a ^= self.poly
        return result

    def square(self, a: int) -> int:
        self.ops["square"] += 1
        return self.mul(a, a)

    @lru_cache(maxsize=1 << 16)
    def inv(self, a: int) -> int:
        assert a != 0
        self.ops["inverse"] += 1
        # Independent Fermat inverse; producer uses extended Euclid.
        exponent = (1 << self.n) - 2
        result = 1
        base = a
        while exponent:
            if exponent & 1:
                result = self.mul(result, base)
            base = self.square(base)
            exponent >>= 1
        assert self.mul(a, result) == 1
        return result

    def trace(self, a: int) -> int:
        value, current = 0, a
        for _ in range(self.n):
            value ^= current
            current = self.square(current)
        assert current == a and value in (0, 1)
        return value

    def trace_mask(self) -> int:
        return sum(self.trace(1 << bit) << bit for bit in range(self.n))


class E:
    def __init__(self, f: GF):
        self.f = f
        self.ops = Counter()

    def on(self, p):
        if p is None:
            return True
        x, y = p
        f = self.f
        return f.square(y) ^ f.mul(x, y) == f.mul(f.square(x), x) ^ 1

    def tau(self, p):
        return None if p is None else (self.f.square(p[0]), self.f.square(p[1]))

    def add(self, p, q):
        self.ops["add"] += 1
        if p is None:
            return q
        if q is None:
            return p
        x, y = p
        u, v = q
        f = self.f
        if x == u:
            if (y ^ v) == x:
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

    def scalar(self, p, k):
        result = None
        while k:
            if k & 1:
                result = self.add(result, p)
            p = self.add(p, p)
            k >>= 1
        return result


def peak_rss_bytes() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def source_order_by_recurrence(n: int) -> int:
    # Four F_2-rational points give t_1=-1 and the Weil trace recurrence.
    a, b = 2, -1
    for _ in range(2, n + 1):
        a, b = b, -(b + 2 * a)
    return (1 << n) + 1 - (b if n > 1 else -1)


def ptr(value):
    return None if value is None else tuple(value)


def read_json(path: Path):
    return json.loads(path.read_text())


def roots_and_points(curve: E):
    f = curve.f
    roots = {}
    for z in range(1 << 13):
        key = f.square(z) ^ z
        roots.setdefault(key, []).append(z)
    assert len(roots) == 1 << 12
    by_x = {}
    points = []
    for x in range(1 << 13):
        if x == 0:
            lifts = [(0, 1)]
        else:
            key = x ^ f.square(f.inv(x))
            lifts = [(x, f.mul(x, z)) for z in roots.get(key, [])]
        assert all(curve.on(p) for p in lifts)
        by_x[x] = sorted(lifts)
        points.extend(lifts)
    assert len(points) == 8011
    return by_x, points


def factor_points(curve: E, by_x, beta: int, m: int, repeated: bool):
    f = curve.f
    conjugates = []
    v = beta
    for _ in range(13):
        conjugates.append(v)
        v = f.square(v)
    assert v == beta
    assert len(set(conjugates)) == 13
    factors = []
    for i in range(m):
        basis = [conjugates[m * j + i] for j in range(2)]
        xs = [0, basis[0], basis[1], basis[0] ^ basis[1]]
        assert len(set(xs)) == 4
        factors.append(sorted(p for x in xs for p in by_x[x]))
    return [factors[0]] * m if repeated else factors


def toy(archive: Path):
    start = time.perf_counter()
    start_cpu = time.process_time()
    f = GF(13, P13)
    curve = E(f)
    by_x, points = roots_and_points(curve)
    summary = read_json(archive / "toy_summary.json")
    assert summary["curve_order"] == len(points) + 1 == source_order_by_recurrence(13) == 8012
    torsion = [None, (0, 1), (1, 0), (1, 1)]
    assert [ptr(t) for t in summary["torsion"]] == torsion
    assert all(curve.on(t) for t in torsion)
    assert curve.scalar(torsion[1], 2) is None
    assert curve.scalar(torsion[2], 2) == curve.scalar(torsion[3], 2) == torsion[1]
    h = next(curve.scalar(p, 4) for p in points if curve.scalar(p, 4) is not None)
    assert ptr(summary["generator"]) == h
    assert curve.scalar(h, 2003) is None and curve.tau(h) == curve.scalar(h, 89)
    subgroup = []
    q = None
    for _ in range(2003):
        subgroup.append(q)
        q = curve.add(q, h)
    assert q is None and len(set(subgroup)) == 2003
    expected_arms = [(beta, m, policy) for beta in (3, 7)
                     for m in (5, 6) for policy in ("rotated", "repeated")]
    assert [(r["beta"], r["m"], r["policy"]) for r in summary["rows"]] == expected_arms
    arm_reports = []
    for meta in summary["rows"]:
        arm_start = time.perf_counter()
        arm_cpu = time.process_time()
        m, beta, repeated = meta["m"], meta["beta"], meta["policy"] == "repeated"
        stem = f"n13-b{beta}-m{m}-{'repeated' if repeated else 'rotated'}"
        factors = factor_points(curve, by_x, beta, m, repeated)
        frozen_factors = read_json(archive / (stem + "-factors.json"))
        assert [[list(p) for p in factor] for factor in factors] == frozen_factors
        histogram = Counter()
        # This enumerates each labelled tuple explicitly; producer uses a DP.
        for choice in itertools.product(*factors):
            acc = None
            for p in choice:
                acc = curve.add(acc, p)
            histogram[acc] += 1
        assert sum(histogram.values()) == math.prod(map(len, factors)) == meta["tuple_count"]
        assert len(histogram) == meta["distinct_full_sums"]
        assert meta["tuple_collisions"] == sum(histogram.values()) - len(histogram)
        assert meta["factor_sizes"] == list(map(len, factors))
        assert meta["physical_point_choices"] == sum(map(len, factors))
        sign_pair_counts = [sum(n == 2 for n in Counter(p[0] for p in factor).values())
                            for factor in factors]
        assert meta["sign_pair_x_counts"] == sign_pair_counts
        projected_factors = [set(curve.scalar(p, 4) for p in factor) for factor in factors]
        assert meta["projected_unique_each_factor"] == list(map(len, projected_factors))
        assert meta["projected_duplicates_each_factor"] == [len(b)-len(s)
                                                                for b,s in zip(factors,projected_factors)]
        assert meta["compressed_log_columns_upper"] == len(projected_factors[0])
        projected = Counter()
        for p, count in histogram.items():
            projected[curve.scalar(p, 4)] += count
        assert len(projected) == meta["distinct_projected_sums"]
        assert meta["projected_tuple_collisions"] == sum(histogram.values()) - len(projected)
        rows = [json.loads(line) for line in (archive / (stem + "-targets.jsonl")).read_text().splitlines()]
        assert len(rows) == 2003
        all_rhss = set()
        witness_replays = 0
        miss_replays = 0
        for k, (q, row) in enumerate(zip(subgroup, rows)):
            assert row["k"] == k and ptr(row["point"]) == q
            proj = curve.scalar(q, 4)
            assert ptr(row["projected_point"]) == proj
            tau = q
            orbit = 0
            while True:
                orbit += 1
                tau = curve.tau(tau)
                if tau == q:
                    break
                assert orbit < 14
            assert row["tau_orbit_size"] == orbit
            assert row["trace_x"] == (None if q is None else f.trace(q[0]))
            torsion_image = curve.scalar(q, 2003)
            assert row["cofactor_class"] == torsion.index(torsion_image)
            assert ptr(row["torsion_projection"]) == torsion_image
            shifted = [curve.add(q, t) for t in torsion]
            assert len(set(shifted)) == 4
            counts = [histogram[p] for p in shifted]
            assert counts == row["coset_multiplicities"]
            assert sum(counts) == projected[proj] == row["projected_multiplicity"]
            assert row["raw_hit"] == (counts[0] > 0)
            assert row["projected_hit"] == (projected[proj] > 0)
            for j, p in enumerate(shifted):
                all_rhss.add(p)
                witness = row["coset_witness_indices"][j]
                if histogram[p]:
                    assert witness is not None and len(witness) == m
                    actual = None
                    for factor, index in zip(factors, witness):
                        actual = curve.add(actual, factor[index])
                    assert actual == p
                    witness_replays += 1
                else:
                    assert witness is None
                    miss_replays += 1
        assert len(all_rhss) == 8012  # Full E(F_2^13), including O and all T cosets.
        assert meta["raw_H_hits"] == sum(row["raw_hit"] for row in rows)
        assert meta["projected_H_hits"] == sum(row["projected_hit"] for row in rows)
        assert meta["raw_nonzero_H_hits"] == sum(row["raw_hit"] for row in rows[1:])
        assert meta["projected_nonzero_H_hits"] == sum(row["projected_hit"] for row in rows[1:])
        assert meta["coset_hits"] == [sum(row["coset_multiplicities"][j] > 0 for row in rows)
                                       for j in range(4)]
        assert meta["coset_rhs_charges"] == 4 * len(rows)
        arm_reports.append({"stem": stem, "tuples": sum(histogram.values()),
                            "distinct_full_sums": len(histogram),
                            "witness_replays": witness_replays,
                            "complete_negative_point_support_replays": miss_replays,
                            "replay_wall_seconds": time.perf_counter() - arm_start,
                            "replay_cpu_seconds": time.process_time() - arm_cpu,
                            "peak_rss_bytes": peak_rss_bytes()})
    return {"toy_arithmetic": "bit-serial reduction and Fermat inverses",
            "toy_replay_wall_seconds": time.perf_counter() - start,
            "toy_replay_cpu_seconds": time.process_time() - start_cpu,
            "toy_peak_rss_bytes": peak_rss_bytes(),
            "toy_operations": dict(f.ops) | dict(curve.ops), "arms": arm_reports}


def batch_inverses(f: GF, values: list[int]) -> list[int]:
    assert all(v != 0 for v in values)
    prefixes = [1]
    for value in values:
        prefixes.append(f.mul(prefixes[-1], value))
    inverse = f.inv(prefixes[-1])
    outputs = [0] * len(values)
    for i in range(len(values) - 1, -1, -1):
        outputs[i] = f.mul(inverse, prefixes[i])
        inverse = f.mul(inverse, values[i])
    assert inverse == 1
    return outputs


def frozen_masks(inputs: Path, entry: dict, count: int) -> list[int]:
    path = inputs / entry["path"]
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == entry["sha256"]
    masks = json.loads(raw)
    assert len(masks) == count and len(set(masks)) == count
    return masks


def independent_masks(kind: str, m: int, d: int, count: int,
                      excluded: set[int] | None = None) -> list[int]:
    values, used = [], set() if excluded is None else set(excluded)
    if kind == "cov":
        for v in (0, 1, 1 << (d - 1), (1 << d) - 1):
            values.append(v)
            used.add(v)
    counter = 0
    while len(values) < count:
        seed = f"{DOMAIN}/{kind}/{m}/{d}/{counter}".encode("ascii")
        v = int.from_bytes(hashlib.sha256(seed).digest(), "big") & ((1 << d) - 1)
        counter += 1
        if v not in used and (kind != "density" or v != 0):
            values.append(v)
            used.add(v)
    return values


def independent_half_trace(f: GF, a: int) -> int:
    assert f.trace(a) == 0
    result = 0
    term = a
    for _ in range(66):
        result ^= term
        term = f.square(f.square(term))
    assert f.square(result) ^ result == a
    return result


def independent_lambda_point(curve: E, tmask: int) -> dict:
    f = curve.f
    lam = LAMBDA131
    assert (lam * lam + lam + 2) % Q131 == 0 and pow(lam, 131, Q131) == 1
    for x in range(2, 1 << 12):
        rhs = x ^ f.square(f.inv(x))
        if ((rhs & tmask).bit_count() & 1) != 0:
            continue
        z = independent_half_trace(f, rhs)
        p = (x, f.mul(x, z))
        assert curve.on(p)
        h = curve.scalar(p, 4)
        if h is None:
            continue
        assert curve.scalar(h, Q131) is None
        assert curve.tau(h) == curve.scalar(h, lam)
        assert curve.scalar(curve.tau(p), 4) == curve.tau(h)
        return {"x": x, "point": list(p), "projected_H": list(h),
                "lambda": lam, "group_order": 4 * Q131,
                "checks": ["point_on_curve", "nonzero_[4]projection",
                           "[q]H=O", "tau(H)=[lambda]H",
                           "[4]tau(P)=tau([4]P)"]}
    raise AssertionError("independent n131 lambda point not found")


def density(archive: Path, inputs: Path):
    start = time.perf_counter()
    start_cpu = time.process_time()
    f = GF(131, P131)
    manifest = read_json(inputs / "input_manifest.json")
    assert manifest["domain"] == DOMAIN and manifest["polynomial"] == P131
    assert source_order_by_recurrence(131) == 4 * Q131
    beta = manifest["beta"]
    assert beta == 3
    conjugates = []
    cur = beta
    for _ in range(131):
        conjugates.append(cur)
        cur = f.square(cur)
    assert cur == beta
    tmask = f.trace_mask()
    summary = read_json(archive / "density_summary.json")
    assert summary["trace_mask"] == tmask
    lambda_point = independent_lambda_point(E(f), tmask)
    assert summary["projected_lambda_point"] == lambda_point
    assert [(row["m"], row["d"]) for row in summary["rows"]] == GRID
    reports = []
    for cell, metadata in zip(manifest["cells"], summary["rows"]):
        m, d = cell["m"], cell["d"]
        assert (m, d) == (metadata["m"], metadata["d"])
        basis = [[conjugates[m * j + i] for j in range(d)] for i in range(m)]
        cov = frozen_masks(inputs, cell["cov"], 256)
        masks = frozen_masks(inputs, cell["density"], 1 << 14)
        assert cov == independent_masks("cov", m, d, 256)
        assert masks == independent_masks("density", m, d, 1 << 14, set(cov))
        rows_cov = [json.loads(line) for line in (archive / f"n131-m{m}-d{d}-covariance.jsonl").read_text().splitlines()]
        rows_density = [json.loads(line) for line in (archive / f"n131-m{m}-d{d}-density.jsonl").read_text().splitlines()]
        assert len(rows_cov) == len(cov) and len(rows_density) == len(masks)
        xs = []
        for mask, row in zip(masks, rows_density):
            x = 0
            for j in range(d):
                if (mask >> j) & 1:
                    x ^= basis[0][j]
            assert x == row["x"] and x != 0 and mask == row["mask"]
            assert ((x & tmask).bit_count() & 1) == (mask.bit_count() & 1)
            xs.append(x)
        inverses = batch_inverses(f, xs)
        count = 0
        for x, inv, row in zip(xs, inverses, rows_density):
            rhs = x ^ f.square(inv)
            solvable = ((rhs & tmask).bit_count() & 1) == 0
            assert solvable == row["solvable"]
            count += int(solvable)
        assert count == metadata["density_solvable"]
        assert metadata["density_masks"] == len(masks)
        assert metadata["covariance_masks"] == len(cov)
        assert math.isclose(metadata["density_rate"], count / len(masks), rel_tol=1e-12)
        z = 1.959963984540054
        p_hat = count / len(masks)
        den = 1 + z*z/len(masks)
        center = (p_hat + z*z/(2*len(masks)))/den
        radius = z*math.sqrt(p_hat*(1-p_hat)/len(masks) + z*z/(4*len(masks)**2))/den
        lower, upper = center-radius, center+radius
        assert all(math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)
                   for a,b in zip(metadata["wilson_95_model_interval"],[lower,upper]))
        ideal = 1 << d
        point_est = 1 + 2 * (ideal-1) * p_hat
        point_up = 1 + 2 * (ideal-1) * upper
        expected_ratios = {
            "conditional_point_size_estimate":point_est,
            "conditional_point_size_upper":point_up,
            "ideal_raw_tuple_ratio_ceiling":min(1,ideal**m/(4*Q131)),
            "ideal_projected_tuple_ratio_ceiling":min(1,ideal**m/Q131),
            "conditional_sampled_raw_ratio":min(1,point_est**m/(4*Q131)),
            "conditional_sampled_projected_ratio":min(1,point_est**m/Q131),
            "conditional_wilson_upper_projected_ratio":min(1,point_up**m/Q131),
        }
        for key,value in expected_ratios.items():
            assert math.isclose(metadata[key],value,rel_tol=1e-12,abs_tol=1e-12), key
        assert metadata["physical_ideal_point_choices"] == m*ideal
        assert metadata["compressed_ideal_point_choices_proxy"] == ideal
        assert metadata["compressed_rigorous_point_choices_upper"] == 2*ideal-1
        cov_xs = []
        for mask, row in zip(cov, rows_cov):
            values = []
            for b in basis:
                x = 0
                for j in range(d):
                    if (mask >> j) & 1:
                        x ^= b[j]
                values.append(x)
            assert mask == row["mask"] and values == row["x"]
            assert all(f.square(values[i]) == values[i + 1] for i in range(m - 1))
            assert all(((x & tmask).bit_count() & 1) == (mask.bit_count() & 1) for x in values)
            cov_xs.append(values[0])
        nonzero_cov = [x for x in cov_xs if x != 0]
        cov_inv = iter(batch_inverses(f, nonzero_cov))
        for x, row in zip(cov_xs, rows_cov):
            if x == 0:
                expected = 1
            else:
                rhs = x ^ f.square(next(cov_inv))
                expected = 2 if ((rhs & tmask).bit_count() & 1) == 0 else 0
            assert row["lifts"] == [expected] * m
        reports.append({"m": m, "d": d, "density_masks_replayed": len(masks),
                        "covariance_masks_replayed": len(cov), "solvable": count})
    return {"density_arithmetic": "independent bit-serial field; batch inversion with Fermat final inverse",
            "density_replay_wall_seconds": time.perf_counter() - start,
            "density_replay_cpu_seconds": time.process_time() - start_cpu,
            "density_peak_rss_bytes": peak_rss_bytes(),
            "density_operations": dict(f.ops),
            "n131_lambda_point": lambda_point, "cells": reports}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--inputs", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = toy(args.archive / "toy")
    if args.inputs is not None:
        report["density"] = density(args.archive / "density", args.inputs)
    args.output.write_text(json.dumps(report, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
