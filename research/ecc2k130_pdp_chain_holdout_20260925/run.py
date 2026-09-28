#!/usr/bin/env python3
"""Frozen, charged m=3 PDP affine-screen pilot. See PROTOCOL.md."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_factor_base_pilot_20260924"))
import run as prior  # noqa: E402

SEED = 2026092503
ATTEMPTS = 256
R = prior.R
N = prior.N


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hash_int(label: str) -> int:
    return int.from_bytes(hashlib.sha256(label.encode()).digest(), "big")


def enc(P):
    return None if P is None else [P[0], P[1]]


def echelon(values, work: Counter):
    pivots = {}
    for value in values:
        while value:
            p = value.bit_length() - 1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
            work["row_xor"] += 1
    return pivots


def orbit(E, P):
    positive, Q = [], P
    while Q not in positive:
        positive.append(Q)
        Q = E.frobenius(Q)
    assert Q == P and len(positive) == N
    return list(dict.fromkeys(positive + [E.neg(Q) for Q in positive]))


def x_orbit_key(F, P):
    if P is None:
        return None
    x, seen = P[0], set()
    while x not in seen:
        seen.add(x)
        x = F.sqr(x)
    assert len(seen) <= N
    return min(seen)


def phase(meter, name, fn, costs):
    before = meter.snapshot()
    result = fn()
    costs[name] = meter.delta(before, meter.snapshot())
    return result


def sum_costs(*costs):
    total = Counter()
    for cost in costs:
        total.update(cost)
    return dict(sorted(total.items()))


def choose_bases(E, B0):
    seen, candidates = set(), []
    for P in B0:
        if P in seen:
            continue
        points = orbit(E, P)
        assert len(points) == 42
        seen.update(points)
        work = Counter()
        rank = len(echelon([Q[0] for Q in points], work))
        key = hashlib.sha256(
            f"pdp-chain-orbit-v1|{P[0]}|{P[1]}".encode()).hexdigest()
        candidates.append({"representative": P, "points": points,
                           "x_span_rank": rank, "hash": key})
    assert len(candidates) >= 2
    chosen = sorted(candidates, key=lambda x: (x["x_span_rank"], x["hash"]))[:2]
    chosen.sort(key=lambda x: x["hash"])
    return candidates, {"train": chosen[0], "held": chosen[1]}


def target_coefficients():
    secret = 1 + hash_int(f"pdp-chain-secret-v1|{SEED}") % (R - 1)
    out = []
    for i in range(ATTEMPTS):
        u = hash_int(f"pdp-chain-target-v1|{SEED}|{i}|u") % R
        v = 1 + hash_int(f"pdp-chain-target-v1|{SEED}|{i}|v") % (R - 1)
        out.append((u, v))
    return secret, out


def target_partition(E, G):
    all_keys = set()
    for d in range(1, R):
        all_keys.add(x_orbit_key(E.F, E.mul(G, d)))
    ranked = sorted(all_keys, key=lambda x: hashlib.sha256(
        f"pdp-chain-target-orbit-v1|{x}".encode()).digest())
    assert len(ranked) >= 3
    return ranked, set(ranked[:math.ceil(len(ranked) / 3)])


def frobenius_scalar(E, G):
    target = E.frobenius(G)
    roots = [x for x in range(R) if (x * x + x + 2) % R == 0]
    matches = [x for x in roots if E.mul(G, x) == target]
    assert len(matches) == 1
    return matches[0]


def point_weights(E, representative, points, tau):
    power, weights = 1, {}
    for _ in range(N):
        P = E.mul(representative, power)
        weights[P] = power
        weights[E.neg(P)] = -power % R
        power = power * tau % R
    assert power == 1 and len(weights) == len(points) == 42
    assert all(P in weights and E.mul(representative, weights[P]) == P
               for P in points)
    return [weights[P] for P in points]


def pair_table(E, points):
    table = defaultdict(list)
    for i, P in enumerate(points):
        for j in range(i, len(points)):
            table[E.add(P, points[j])].append((i, j))
    return dict(table)


class AffineScreen:
    """The #706 S3 row-span test, with curve b and cached base products."""

    def __init__(self, E, points):
        self.F, self.b = E.F, E.b
        self.work = Counter()
        self.basis = list(echelon([P[0] for P in points], self.work).values())
        self.products = [self.F.mul(a, b) for a in self.basis for b in self.basis]
        self.cache = {}

    def contradicts(self, t):
        if t in self.cache:
            return self.cache[t], False
        F, basis, products = self.F, self.basis, self.products
        linear = [F.mul(F.sqr(t), F.sqr(u)) for u in basis]
        cross = [F.sqr(z) ^ F.mul(t, z) for z in products]
        offset = 1 + 2 * len(basis)
        rows = []
        for bit in range(F.deg):
            row = (self.b >> bit) & 1
            for j, value in enumerate(linear):
                row |= ((value >> bit) & 1) << (j + 1)
                row |= ((value >> bit) & 1) << (j + 1 + len(basis))
            for j, value in enumerate(cross):
                row |= ((value >> bit) & 1) << (offset + j)
            rows.append(row)
        inconsistent = 1 in echelon(rows, self.work).values()
        self.cache[t] = inconsistent
        return inconsistent, True


def split_name(base_kind, target_key, held_keys):
    if target_key is None:
        return "infinity"
    held = target_key in held_keys
    if base_kind == "train":
        return "orbit_holdout" if held else "train"
    return "confirmation" if held else "base_holdout"


def algorithm(E, points, weights, G, Q, targets, target_keys, held_keys,
              base_kind, *, screen, meter):
    phase_costs = {}
    table = phase(meter, "pair_table", lambda: pair_table(E, points), phase_costs)
    screen_obj = (phase(meter, "screen_setup",
                        lambda: AffineScreen(E, points), phase_costs)
                  if screen else None)
    tracker = prior.RankTracker(width=2, modulus=R)
    records = []
    first_rank = None
    skipped = []
    for i, ((u, v, T), key) in enumerate(zip(targets, target_keys)):
        before = meter.snapshot()
        before_xor = screen_obj.work["row_xor"] if screen_obj else 0
        witness, tested, skips, lookups, cache_misses = None, 0, 0, 0, 0
        for k, P3 in enumerate(points):
            tested += 1
            U = E.add(T, E.neg(P3))
            if screen_obj is not None and U is not None:
                contradiction, miss = screen_obj.contradicts(U[0])
                cache_misses += int(miss)
                if contradiction:
                    skips += 1
                    skipped.append([i, k, U[0]])
                    continue
            lookups += 1
            pairs = table.get(U)
            if pairs:
                a, b = pairs[0]
                witness = [a, b, k]
                assert E.add(E.add(points[a], points[b]), P3) == T
                break
        gained = False
        if witness is not None and first_rank is None:
            row = [sum(weights[j] for j in witness) % R, -v % R]
            gained = tracker.add(row, u)
            if len(tracker.pivots) == 2:
                first_rank = i + 1
        rec = {
            "i": i, "u": u, "v": v, "target": enc(T),
            "target_orbit": key, "split": split_name(base_kind, key, held_keys),
            "first_witness": witness, "hit": witness is not None,
            "third_candidates_tested": tested, "certified_skips": skips,
            "pair_lookups": lookups, "profile_cache_misses": cache_misses,
            "independent_before_rank_stop": gained,
            "rank_after": len(tracker.pivots),
            "scan_cost": meter.delta(before, meter.snapshot()),
            "scan_row_xor": (
                screen_obj.work["row_xor"] - before_xor if screen_obj else 0),
        }
        records.append(rec)
    solution = None
    verified = False
    if first_rank is not None:
        before = meter.snapshot()
        solution = tracker.solve()
        verified = (solution[1] == target_coefficients()[0]
                    and E.mul(G, solution[1]) == Q
                    and E.mul(G, solution[0]) == points[0])
        # points[0] is the selected orbit representative by construction.
        phase_costs["recovery_verify"] = meter.delta(before, meter.snapshot())
    else:
        phase_costs["recovery_verify"] = {}
    assert verified, f"{base_kind} did not recover a verified scalar"
    phase_costs["scan_all"] = sum_costs(*(r["scan_cost"] for r in records))
    phase_costs["scan_to_rank"] = sum_costs(
        *(r["scan_cost"] for r in records[:first_rank]))
    pair_entries = sum(map(len, table.values()))
    assert pair_entries == 42 * 43 // 2
    return {
        "base_kind": base_kind, "policy": "affine_screen" if screen else "baseline",
        "base_size": len(points), "x_span_rank": (
            len(screen_obj.basis) if screen_obj else
            len(echelon([P[0] for P in points], Counter()))),
        "pair_entries": pair_entries, "distinct_pair_sums": len(table),
        "screen_cache_entries": len(screen_obj.cache) if screen_obj else 0,
        "screen_row_xor": dict(screen_obj.work) if screen_obj else {},
        "skipped_residuals": skipped,
        "first_full_rank_attempt": first_rank, "rank": len(tracker.pivots),
        "solution": solution, "verified": verified,
        "rank_ops": dict(tracker.ops), "phase_costs": phase_costs,
        "records": records,
    }


def summarize(variants):
    out = {}
    for key, arm in variants.items():
        cells = {}
        for name in ("train", "base_holdout", "orbit_holdout",
                     "confirmation", "infinity"):
            rows = [r for r in arm["records"] if r["split"] == name]
            cells[name] = {
                "attempts": len(rows), "hits": sum(r["hit"] for r in rows),
                "certified_skips": sum(r["certified_skips"] for r in rows),
                "pair_lookups": sum(r["pair_lookups"] for r in rows),
                "independent_gain_before_rank_stop": sum(
                    r["independent_before_rank_stop"] for r in rows),
                "scan_cost": sum_costs(*(r["scan_cost"] for r in rows)),
                "scan_row_xor": sum(r["scan_row_xor"] for r in rows),
            }
        out[key] = cells
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "preserve prior evidence"
    costs = {}
    meter = prior.Meter()
    def setup_field():
        field = prior.CountedField(meter)
        field.frobenius(1, N - 1)
        return field

    F = phase(meter, "field_setup", setup_field, costs)
    E0 = prior.CountedCurve(meter, F, 0, 1)
    twist = prior.CountedCurve(meter, F, 1, 1)
    lines, selected, trials = phase(
        meter, "kernel_line_search",
        lambda: prior.order_seven_lines(F, twist), costs)
    kernel = lines[selected]
    phi = phase(
        meter, "map_construction",
        lambda: prior.BinaryVeluMap.from_generator(
            E0, twist, kernel["generator"], prior.ELL), costs)
    E1 = prior.CountedCurve(meter, F, phi.codomain.a, phi.codomain.b)
    assert E1.b == kernel["codomain_b"] == 0x11584f
    G = phase(meter, "generator_search",
              lambda: prior.choose_generator(E0)[0], costs)
    B0 = phase(meter, "candidate_base_scan",
               lambda: prior.candidate_base(E0)[0], costs)
    candidates, selected_bases = phase(
        meter, "signed_orbit_selection",
        lambda: choose_bases(E0, B0), costs)
    tau = phase(meter, "frobenius_scalar",
                lambda: frobenius_scalar(E0, G), costs)
    base = {}
    for kind, record in selected_bases.items():
        points = record["points"]
        weights = phase(
            meter, f"orbit_weights_{kind}",
            lambda points=points: point_weights(
                E0, record["representative"], points, tau), costs)
        images = phase(meter, f"base_transport_{kind}",
                       lambda points=points: [phi(P) for P in points], costs)
        assert len(set(images)) == 42
        assert all(P is not None and E1.on_curve(P)
                   and E1.mul(P, R) is None for P in images)
        base[kind] = {"original": points, "transported": images,
                      "weights": weights, "representative": record["representative"],
                      "hash": record["hash"]}
    ranked_keys, held_keys = phase(
        meter, "target_orbit_partition_audit",
        lambda: target_partition(E0, G), costs)
    secret, coefficients = target_coefficients()
    Q = phase(meter, "challenge_generation_audit",
              lambda: E0.mul(G, secret), costs)
    assert Q is not None and E0.mul(Q, R) is None
    G1, Q1 = phase(meter, "generator_challenge_transport",
                   lambda: (phi(G), phi(Q)), costs)
    assert G1 is not None and Q1 == E1.mul(G1, secret)
    targets0, targets1, source_target_costs, codomain_target_costs = [], [], [], []
    for u, v in coefficients:
        before = meter.snapshot()
        T0 = E0.add(E0.mul(G, u), E0.mul(Q, v))
        targets0.append((u, v, T0))
        source_target_costs.append(meter.delta(before, meter.snapshot()))
        before = meter.snapshot()
        T1 = E1.add(E1.mul(G1, u), E1.mul(Q1, v))
        targets1.append((u, v, T1))
        codomain_target_costs.append(meter.delta(before, meter.snapshot()))
    target_keys = phase(
        meter, "target_key_audit",
        lambda: [x_orbit_key(F, T) for _, _, T in targets0], costs)
    assert all(k is None or k in ranked_keys for k in target_keys)
    parity = phase(
        meter, "target_map_audit",
        lambda: all(phi(T0) == T1 for (_, _, T0), (_, _, T1)
                    in zip(targets0, targets1)), costs)
    assert parity
    arms = {}
    for kind in ("train", "held"):
        for geometry, E, gen, challenge, targets in [
            ("original", E0, G, Q, targets0),
            ("transported", E1, G1, Q1, targets1),
        ]:
            for screen in (False, True):
                key = f"{kind}_{geometry}_{'screen' if screen else 'baseline'}"
                arms[key] = algorithm(
                    E, base[kind][geometry], base[kind]["weights"],
                    gen, challenge, targets, target_keys, held_keys,
                    kind, screen=screen, meter=meter)
    for kind in ("train", "held"):
        for geometry in ("original", "transported"):
            a = arms[f"{kind}_{geometry}_baseline"]
            b = arms[f"{kind}_{geometry}_screen"]
            assert [(r["first_witness"], r["rank_after"]) for r in a["records"]] == [
                (r["first_witness"], r["rank_after"]) for r in b["records"]]
            assert all(not a["records"][i]["hit"]
                       for i, _, _ in b["skipped_residuals"]
                       if not b["records"][i]["hit"])
            assert a["solution"] == b["solution"]
        for mode in ("baseline", "screen"):
            a = arms[f"{kind}_original_{mode}"]
            b = arms[f"{kind}_transported_{mode}"]
            assert [(r["first_witness"], r["rank_after"]) for r in a["records"]] == [
                (r["first_witness"], r["rank_after"]) for r in b["records"]]
            assert a["solution"] == b["solution"]
    summary = summarize(arms)
    cold = {}
    for name, arm in arms.items():
        kind, geometry, mode = name.split("_")
        prefix = arm["first_full_rank_attempt"]
        assert prefix is not None
        shared = [costs["field_setup"], costs["generator_search"],
                  costs["candidate_base_scan"], costs["signed_orbit_selection"],
                  costs["frobenius_scalar"], costs[f"orbit_weights_{kind}"],
                  costs["challenge_generation_audit"]]
        if geometry == "transported":
            shared += [costs["kernel_line_search"], costs["map_construction"],
                       costs[f"base_transport_{kind}"],
                       costs["generator_challenge_transport"]]
            shared += codomain_target_costs[:prefix]
        else:
            shared += source_target_costs[:prefix]
        phase_c = arm["phase_costs"]
        shared += [phase_c["pair_table"], phase_c["scan_to_rank"],
                   phase_c["recovery_verify"]]
        if mode == "screen":
            shared.append(phase_c["screen_setup"])
        cold[name] = sum_costs(*shared)
        cold[name]["row_xor"] = sum(
            r["scan_row_xor"] for r in arm["records"][:prefix])
        cold[name]["pair_lookups"] = sum(
            r["pair_lookups"] for r in arm["records"][:prefix])
    artifact = {
        "schema": "ecc2k130-pdp-chained-orbit-holdout-v1",
        "status": "TOY_DIAGNOSTIC",
        "parameters": {"seed": SEED, "degree": N,
                       "irreducible": prior.IRR, "subgroup_order": R,
                       "attempts": ATTEMPTS, "isogeny_degree": prior.ELL,
                       "source_order": prior.ORDER,
                       "source_cofactor": prior.COFACTOR},
        "source_sha256": {
            "protocol": digest(Path(__file__).with_name("PROTOCOL.md")),
            "runner": digest(Path(__file__)),
            "prior_degree7_runner": digest(ROOT / "research/ecc2k130_factor_base_pilot_20260924/run.py"),
            "fastfield": digest(ROOT / "research/ecc2k130_relations/fastfield.py"),
            "relations": digest(ROOT / "research/ecc2k130_relations/relations.py"),
            "oriented_velu": digest(ROOT / "research/ecc2k130_oriented_transport_20260924/oriented_velu.py"),
        },
        "geometry": {
            "kernel_line_search_trials": trials,
            "selected_kernel_x": list(selected),
            "kernel_generator": enc(kernel["generator"]),
            "codomain_b": E1.b, "generator": enc(G),
            "transported_generator": enc(G1), "challenge": enc(Q),
            "transported_challenge": enc(Q1), "secret_audit_only": secret,
            "frobenius_scalar_mod_r": tau,
            "candidate_orbits": [
                {"representative": enc(c["representative"]),
                 "x_span_rank": c["x_span_rank"], "hash": c["hash"]}
                for c in candidates],
            "bases": {
                k: {"representative": enc(v["representative"]),
                    "hash": v["hash"], "weights": v["weights"],
                    "original": [enc(P) for P in v["original"]],
                    "transported": [enc(P) for P in v["transported"]]}
                for k, v in base.items()},
            "all_target_orbits": ranked_keys,
            "held_target_orbits": sorted(held_keys),
        },
        "target_coefficients": [list(x) for x in coefficients],
        "source_targets": [enc(T) for _, _, T in targets0],
        "transported_targets": [enc(T) for _, _, T in targets1],
        "target_keys": target_keys,
        "setup_costs": costs,
        "source_target_costs": source_target_costs,
        "codomain_target_costs": codomain_target_costs,
        "arms": arms, "summary": summary,
        "cold_cost_to_verified_rank": cold,
        "limitations": [
            "toy order-421 subgroup, 42-point one-orbit bases, m=3 exact pair table",
            "transferred source orbit is not a native codomain Frobenius orbit",
            "affine screen applies to each fixed-third residual, not the whole S4 chain",
            "lookup, row XOR, field and curve counters are separate native units",
            "no degree-131 or rho crossover inference",
        ],
        "platform": platform.platform(),
        "python": sys.version,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.out), "secret_verified": all(
            arm["verified"] for arm in arms.values()),
        "base_x_span_ranks": {
            k: {g: arms[f"{k}_{g}_baseline"]["x_span_rank"]
                for g in ("original", "transported")}
            for k in ("train", "held")},
        "first_rank": {k: v["first_full_rank_attempt"]
                       for k, v in arms.items()},
        "certified_skips": {k: sum(r["certified_skips"] for r in v["records"])
                            for k, v in arms.items()},
    }, indent=2))


if __name__ == "__main__":
    main()
