#!/usr/bin/env python3
"""Independent group, S3, target-stream and rank replay of compact receipts."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_oriented_transport_20260924"))
from fastfield import FastGF2m  # noqa: E402
from relations import Koblitz, solve_s3_last  # noqa: E402
from oriented_velu import BinaryVeluMap  # noqa: E402

R, H, N, SEED = 421, 4988, 21, 2026092503


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hashed(s):
    return int.from_bytes(hashlib.sha256(s.encode()).digest(), "big")


def point(value):
    return None if value is None else tuple(value)


def encode(P):
    return None if P is None else [P[0], P[1]]


def orbit_key(F, P):
    if P is None:
        return None
    x, seen = P[0], set()
    while x not in seen:
        seen.add(x)
        x = F.sqr(x)
    return min(seen)


def signed_orbit(E, P):
    positive, Q = [], P
    while Q not in positive:
        positive.append(Q)
        Q = E.frobenius(Q)
    assert Q == P and len(positive) == N
    return positive + [E.neg(Q) for Q in positive]


def rank2(values):
    pivots = {}
    for value in values:
        while value:
            p = value.bit_length() - 1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
    return len(pivots)


def source_hash_check(data):
    here = Path(__file__).resolve().parent
    expected = data["source_sha256"]
    source_map = {
        "prior_degree7_runner": ROOT / "research/ecc2k130_factor_base_pilot_20260924/run.py",
        "fastfield": ROOT / "research/ecc2k130_relations/fastfield.py",
        "relations": ROOT / "research/ecc2k130_relations/relations.py",
        "oriented_velu": ROOT / "research/ecc2k130_oriented_transport_20260924/oriented_velu.py",
        "projected_runner": here / "run.py",
        "raw_runner": here / "run_raw.py",
        "historical_512_runner": here / "run_mask_512_v1.py",
        "protocol_mask": here / "PROTOCOL_MASK.md",
    }
    schema = data["schema"]
    if "power" in schema:
        source_map["protocol"] = here / "PROTOCOL_POWER.md"
        source_map["runner"] = here / "run_mask.py"
    elif "masked" in schema:
        source_map["protocol"] = here / "PROTOCOL_MASK.md"
        source_map["runner"] = here / "run_mask_512_v1.py"
    elif "raw" in schema:
        source_map["protocol"] = here / "PROTOCOL_RAW.md"
        source_map["runner"] = here / "run_raw.py"
    else:
        source_map["protocol"] = here / "PROTOCOL.md"
        source_map["runner"] = here / "run.py"
    for key, old_hash in expected.items():
        path = source_map.get(key)
        assert path is not None and path.is_file(), (key, path)
        assert sha(path) == old_hash, (key, sha(path), old_hash)
    assert sha(here / "compact.py") == data["compactor_sha256"]


def verify_selection(E, data):
    schema = data["schema"]
    if "power" in schema or "masked" in schema:
        trials, rank_cap, seed, domain = (
            200_000, 7, 2026092504, "pdp-chain-mask-orbit-v1")
        reported = data["geometry"]["eligible_orbits"]
    elif "raw" in schema:
        trials, rank_cap, seed, domain = (
            20_000, 9, SEED, "pdp-chain-raw-orbit-v1")
        reported = data["geometry"]["eligible_orbits"]
    else:
        return
    rng, seen, eligible = random.Random(seed), set(), []
    for trial in range(trials):
        x = rng.randrange(1, 1 << N)
        orbit_x, q = [], x
        for _ in range(N):
            orbit_x.append(q)
            q = E.F.sqr(q)
        dim = rank2(orbit_x)
        if dim > rank_cap:
            continue
        lifts = E.points_over(x)
        if not lifts:
            continue
        P = min(lifts)
        positive, q = [], P
        while q not in positive:
            positive.append(q)
            q = E.frobenius(q)
        if q != P or len(positive) != N:
            continue
        points = positive + [E.neg(Q) for Q in positive]
        if len(set(points)) != 42:
            continue
        projected = E.mul(P, H)
        if projected is None:
            continue
        key = frozenset(points)
        if key in seen:
            continue
        seen.add(key)
        identifier = hashlib.sha256(
            f"{domain}|{P[0]}|{P[1]}".encode()).hexdigest()
        eligible.append({
            "trial": trial, "representative": encode(P),
            "projected_representative": encode(projected),
            "x_span_rank": dim, "hash": identifier})
    assert eligible == reported, "geometry-selected candidates changed"
    selected = sorted(eligible, key=lambda c: (
        c["x_span_rank"], c["hash"]))[:2]
    selected.sort(key=lambda c: c["hash"])
    for kind, item in zip(("train", "held"), selected):
        base = data["geometry"]["bases"][kind]
        for key in ("trial", "representative", "projected_representative",
                    "hash"):
            assert base[key] == item[key], (kind, key)


def verify_mask(E, i, item):
    for retry in range(item["retries"] + 1):
        x = hashed(f"pdp-chain-mask-v1|{SEED}|{i}|{retry}|x") & ((1 << N) - 1)
        sign = hashed(f"pdp-chain-mask-v1|{SEED}|{i}|{retry}|sign") & 1
        points = E.points_over(x)
        accepted = bool(points) and not (x == 0 and sign)
        assert accepted == (retry == item["retries"])
        if accepted:
            W = sorted(points)[sign]
            M = E.mul(W, R)
            assert encode(W) == item["source_point"]
            assert encode(M) == item["mask"]
            assert E.mul(M, H) is None
            return M
    raise AssertionError("unreachable")


def pair_table(E, points):
    table = defaultdict(list)
    for a in range(len(points)):
        for b in range(a, len(points)):
            table[E.add(points[a], points[b])].append((a, b))
    return table


def finite_s3_root_set(E, points):
    xs = sorted({P[0] for P in points})
    pairs = [(xs[i], xs[j]) for i in range(len(xs))
             for j in range(i, len(xs))]
    roots = solve_s3_last(E, pairs)
    result = set()
    for (x, y), values in zip(pairs, roots):
        for t in values:
            assert E.s3(x, y, t) == 0
            result.add(t)
    return result


class Rank:
    def __init__(self):
        self.pivots = {}

    def add(self, a, b, rhs):
        row = [a % R, b % R]
        rhs %= R
        for p in sorted(self.pivots):
            lead, const = self.pivots[p]
            factor = row[p]
            if factor:
                row = [(x - factor * y) % R for x, y in zip(row, lead)]
                rhs = (rhs - factor * const) % R
        p = next((i for i, x in enumerate(row) if x), None)
        if p is None:
            assert rhs == 0
            return False
        inv = pow(row[p], -1, R)
        self.pivots[p] = ([x * inv % R for x in row], rhs * inv % R)
        return True

    def solve(self):
        assert len(self.pivots) == 2
        last = self.pivots[1][1]
        first_row, first_rhs = self.pivots[0]
        first = (first_rhs - first_row[1] * last) % R
        return [first, last]


def verify_arm(E, G, Q, points, weights, targets, keys, held_keys, arm,
               *, projected, masked):
    table = pair_table(E, points)
    assert sum(map(len, table.values())) == arm["pair_entries"] == 903
    assert len(table) == arm["distinct_pair_sums"]
    root_x = finite_s3_root_set(E, points)
    rank = Rank()
    cell_ranks = defaultdict(Rank)
    first_by_cell = {}
    first_rank = None
    cache = set()
    counts = Counter()
    for i, (target, key, case) in enumerate(zip(targets, keys, arm["cases"])):
        u, v = target[0], target[1]
        T = target[2]
        witness, tested, skip_mask, lookups, misses, gained, after = case
        assert skip_mask >> tested == 0, (i, skip_mask, tested)
        if arm["policy"] == "baseline":
            assert skip_mask == 0 and misses == 0
        found, observed_lookups, observed_misses = None, 0, 0
        for k, P3 in enumerate(points):
            U = E.add(T, E.neg(P3))
            skip = bool(skip_mask & (1 << k))
            if arm["policy"] == "affine_screen" and U is not None:
                if U[0] not in cache:
                    cache.add(U[0])
                    observed_misses += 1
            if skip:
                assert U is not None and U not in table
                assert U[0] not in root_x, (
                    "certified skip had a finite-base algebraic S3 root", i, k)
                continue
            observed_lookups += 1
            if U in table:
                a, b = table[U][0]
                found = [a, b, k]
                assert E.add(E.add(points[a], points[b]), P3) == T
                if masked:
                    actual = E.mul(T, H)
                    expected = E.add(E.mul(G, H * u), E.mul(Q, H * v))
                    assert actual == expected
                break
        assert found == witness and observed_lookups == lookups
        assert observed_misses == misses, (i, observed_misses, misses)
        assert (found[2] + 1 if found else len(points)) == tested
        counts["hits"] += int(found is not None)
        counts["skips"] += skip_mask.bit_count()
        counts["lookups"] += observed_lookups
        cell = ("infinity" if key is None else
                "orbit_holdout" if arm["base_kind"] == "train" and key in held_keys else
                "train" if arm["base_kind"] == "train" else
                "confirmation" if key in held_keys else "base_holdout")
        counts[f"{cell}_attempts"] += 1
        counts[f"{cell}_hits"] += int(found is not None)
        counts[f"{cell}_skips"] += skip_mask.bit_count()
        increase = False
        if found is not None:
            scale = 1 if projected else H
            a = sum(weights[j] for j in found) % R
            cell_ranks[cell].add(a, -(scale * v), scale * u)
            if len(cell_ranks[cell].pivots) == 2 and cell not in first_by_cell:
                first_by_cell[cell] = i + 1
            if first_rank is None:
                increase = rank.add(a, -(scale * v), scale * u)
                if len(rank.pivots) == 2:
                    first_rank = i + 1
        assert bool(gained) == increase and after == len(rank.pivots)
    assert arm["screen_cache_entries"] == len(cache)
    assert arm["rank"] == len(rank.pivots)
    assert arm["first_full_rank_attempt"] == first_rank
    solution = rank.solve() if first_rank is not None else None
    assert arm["solution"] == solution
    assert bool(arm["verified"]) == (solution is not None)
    rep = points[0] if projected else E.mul(points[0], H)
    if solution is not None:
        assert E.mul(G, solution[0]) == rep
        assert E.mul(G, solution[1]) == Q
    for cell, cell_rank in cell_ranks.items():
        counts[f"{cell}_rank"] = len(cell_rank.pivots)
        if len(cell_rank.pivots) == 2:
            local = cell_rank.solve()
            assert E.mul(G, local[0]) == rep
            assert E.mul(G, local[1]) == Q
            counts[f"{cell}_first_rank_attempt"] = first_by_cell[cell]
            counts[f"{cell}_recovered_k"] = local[1]
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--skip-selection", action="store_true")
    args = parser.parse_args()
    data = json.loads(args.receipt.read_text())
    source_hash_check(data)
    param = data["parameters"]
    assert param["degree"] == N and param["subgroup_order"] == R
    assert param["irreducible"] == 0x200005
    attempts = param["attempts"]
    assert attempts in (256, 512, 2048)
    F = FastGF2m(N, 0x200005)
    E0, Twist = Koblitz(F, 0, 1), Koblitz(F, 1, 1)
    kernel = point(data["geometry"]["kernel_generator"])
    assert Twist.on_curve(kernel) and Twist.mul(kernel, 7) is None
    phi = BinaryVeluMap.from_generator(E0, Twist, kernel, 7)
    E1 = Koblitz(F, phi.codomain.a, phi.codomain.b)
    assert E1.b == data["geometry"]["codomain_b"] == 0x11584f
    G, G1 = (point(data["geometry"][k]) for k in
             ("generator", "transported_generator"))
    Q, Q1 = (point(data["geometry"][k]) for k in
             ("challenge", "transported_challenge"))
    assert E0.on_curve(G) and E0.mul(G, R) is None
    assert phi(G) == G1 and phi(Q) == Q1
    secret = 1 + hashed(f"pdp-chain-secret-v1|{SEED}") % (R - 1)
    assert data["geometry"]["secret_audit_only"] == secret
    assert E0.mul(G, secret) == Q and E1.mul(G1, secret) == Q1
    if not args.skip_selection:
        verify_selection(E0, data)
    projected = "pdp-chained-orbit" in data["schema"]
    masked = "masked" in data["schema"] or "power" in data["schema"]
    bases = {}
    for kind in ("train", "held"):
        row = data["geometry"]["bases"][kind]
        original = [point(x) for x in row["original"]]
        mapped = [point(x) for x in row["transported"]]
        assert len(set(original)) == len(set(mapped)) == 42
        assert original == signed_orbit(E0, point(row["representative"]))
        assert [phi(P) for P in original] == mapped
        assert all(E0.on_curve(P) for P in original)
        assert all(E1.on_curve(P) for P in mapped)
        assert len(row["weights"]) == 42
        for E, ps, gen in ((E0, original, G), (E1, mapped, G1)):
            rep = ps[0] if projected else E.mul(ps[0], H)
            for P, w in zip(ps, row["weights"]):
                label = P if projected else E.mul(P, H)
                assert E.mul(rep, w) == label
        bases[kind] = {"original": original, "transported": mapped,
                       "weights": row["weights"]}
    assert not set(bases["train"]["original"]) & set(bases["held"]["original"])
    source_targets, mapped_targets, keys = [], [], []
    coeffs = data["target_coefficients"]
    assert len(coeffs) == attempts
    masks = data.get("masks")
    if masked:
        assert len(masks) == attempts
    for i, (u, v) in enumerate(coeffs):
        assert u == hashed(f"pdp-chain-target-v1|{SEED}|{i}|u") % R
        assert v == 1 + hashed(f"pdp-chain-target-v1|{SEED}|{i}|v") % (R - 1)
        M = verify_mask(E0, i, masks[i]) if masked else None
        T0 = E0.add(M, E0.add(E0.mul(G, u), E0.mul(Q, v)))
        T1 = E1.add(phi(M), E1.add(E1.mul(G1, u), E1.mul(Q1, v)))
        assert encode(T0) == data["source_targets"][i]
        assert encode(T1) == data["transported_targets"][i]
        assert phi(T0) == T1
        key = orbit_key(F, T0)
        assert key == data["target_keys"][i]
        source_targets.append((u, v, T0))
        mapped_targets.append((u, v, T1))
        keys.append(key)
    held_keys = set(data["geometry"]["held_target_orbits"])
    if masked:
        assert held_keys == {key for key in keys if key is not None and
            int.from_bytes(hashlib.sha256(
                f"pdp-chain-mask-target-orbit-v1|{key}".encode()).digest()[:8],
                "big") % 3 == 0}
    else:
        assert all(key in data["geometry"]["all_target_orbits"]
                   for key in keys if key is not None)
    arm_counts = {}
    for kind in ("train", "held"):
        for geometry, E, gen, challenge, targets in (
            ("original", E0, G, Q, source_targets),
            ("transported", E1, G1, Q1, mapped_targets),
        ):
            for mode in ("baseline", "screen"):
                name = f"{kind}_{geometry}_{mode}"
                arm = data["arms"][name]
                assert len(arm["cases"]) == attempts
                arm_counts[name] = verify_arm(
                    E, gen, challenge, bases[kind][geometry],
                    bases[kind]["weights"], targets, keys, held_keys,
                    arm, projected=projected, masked=masked)
    for kind in ("train", "held"):
        for geometry in ("original", "transported"):
            a, b = (data["arms"][f"{kind}_{geometry}_{mode}"]
                    for mode in ("baseline", "screen"))
            assert [(r[0], r[5], r[6]) for r in a["cases"]] == [
                (r[0], r[5], r[6]) for r in b["cases"]]
        for mode in ("baseline", "screen"):
            a, b = (data["arms"][f"{kind}_{geometry}_{mode}"]
                    for geometry in ("original", "transported"))
            assert [(r[0], r[5], r[6]) for r in a["cases"]] == [
                (r[0], r[5], r[6]) for r in b["cases"]]
    print(json.dumps({
        "status": "PASS", "schema": data["schema"], "attempts": attempts,
        "hits": {k: c["hits"] for k, c in arm_counts.items()},
        "certified_skips": {k: c["skips"] for k, c in arm_counts.items()},
        "first_rank": {k: a["first_full_rank_attempt"]
                       for k, a in data["arms"].items()},
        "confirmation_cell_rank": {
            k: {"rank": c["confirmation_rank"],
                "first_global_attempt": c.get("confirmation_first_rank_attempt"),
                "recovered_k": c.get("confirmation_recovered_k")}
            for k, c in arm_counts.items() if k.startswith("held_")},
    }, indent=2))


if __name__ == "__main__":
    main()
