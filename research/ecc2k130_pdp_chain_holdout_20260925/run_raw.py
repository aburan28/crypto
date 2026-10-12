#!/usr/bin/env python3
"""Predeclared sparse raw-orbit m=3 PDP control; see PROTOCOL_RAW.md."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
pilot_path = Path(__file__).with_name("run.py")
spec = importlib.util.spec_from_file_location("pdp_projected_runner", pilot_path)
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)

prior = pilot.prior
N, R, H, SEED, ATTEMPTS = prior.N, prior.R, prior.COFACTOR, pilot.SEED, 512


def raw_bases(E):
    rng = random.Random(SEED)
    stats = Counter()
    seen, eligible = set(), []
    for trial in range(20_000):
        x = rng.randrange(1, 1 << N)
        stats["x_trials"] += 1
        x_orbit, q = [], x
        for _ in range(N):
            x_orbit.append(q)
            q = E.F.sqr(q)
        rank = len(pilot.echelon(x_orbit, stats))
        if rank > 9:
            stats["rank_rejected"] += 1
            continue
        stats["lift_checked"] += 1
        lifts = E.points_over(x)
        if not lifts:
            stats["lift_rejected"] += 1
            continue
        P = min(lifts)
        try:
            points = pilot.orbit(E, P)
        except AssertionError:
            stats["short_orbit_rejected"] += 1
            continue
        if len(points) != 42:
            stats["short_orbit_rejected"] += 1
            continue
        projected = E.mul(P, H)
        if projected is None:
            stats["zero_projection_rejected"] += 1
            continue
        assert E.mul(projected, R) is None
        key = frozenset(points)
        if key in seen:
            stats["duplicate_orbit_rejected"] += 1
            continue
        seen.add(key)
        identifier = hashlib.sha256(
            f"pdp-chain-raw-orbit-v1|{P[0]}|{P[1]}".encode()).hexdigest()
        eligible.append({"trial": trial, "x": x, "representative": P,
                         "projected_representative": projected,
                         "points": points, "x_span_rank": rank,
                         "hash": identifier})
    assert len(eligible) >= 2
    selected = sorted(eligible, key=lambda c: (
        c["x_span_rank"], c["hash"]))[:2]
    selected.sort(key=lambda c: c["hash"])
    assert not set(selected[0]["points"]) & set(selected[1]["points"])
    return stats, eligible, {"train": selected[0], "held": selected[1]}


def target_coefficients():
    secret = 1 + pilot.hash_int(
        f"pdp-chain-secret-v1|{SEED}") % (R - 1)
    return secret, [
        (pilot.hash_int(f"pdp-chain-target-v1|{SEED}|{i}|u") % R,
         1 + pilot.hash_int(f"pdp-chain-target-v1|{SEED}|{i}|v") % (R - 1))
        for i in range(ATTEMPTS)]


def raw_algorithm(E, points, weights, G, Q, targets, target_keys,
                  held_keys, base_kind, *, screen, meter):
    phase_costs = {}
    table = pilot.phase(
        meter, "pair_table", lambda: pilot.pair_table(E, points), phase_costs)
    aff = (pilot.phase(
        meter, "screen_setup", lambda: pilot.AffineScreen(E, points), phase_costs)
        if screen else None)
    tracker = prior.RankTracker(width=2, modulus=R)
    first_rank, records, skipped = None, [], []
    for i, ((u, v, T), key) in enumerate(zip(targets, target_keys)):
        before = meter.snapshot()
        xor_before = aff.work["row_xor"] if aff else 0
        witness = None
        tested = skips = lookups = misses = 0
        for k, P3 in enumerate(points):
            tested += 1
            U = E.add(T, E.neg(P3))
            if aff is not None and U is not None:
                contradiction, miss = aff.contradicts(U[0])
                misses += int(miss)
                if contradiction:
                    skips += 1
                    skipped.append([i, k, pilot.enc(U)])
                    assert U not in table, "screen discarded a group pair"
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
            row = [sum(weights[j] for j in witness) % R,
                   -(H % R) * v % R]
            gained = tracker.add(row, (H % R) * u % R)
            if len(tracker.pivots) == 2:
                first_rank = i + 1
        records.append({
            "i": i, "u": u, "v": v, "target": pilot.enc(T),
            "target_orbit": key,
            "split": pilot.split_name(base_kind, key, held_keys),
            "first_witness": witness, "hit": witness is not None,
            "third_candidates_tested": tested,
            "certified_skips": skips, "pair_lookups": lookups,
            "profile_cache_misses": misses,
            "independent_before_rank_stop": gained,
            "rank_after": len(tracker.pivots),
            "scan_cost": meter.delta(before, meter.snapshot()),
            "scan_row_xor": aff.work["row_xor"] - xor_before if aff else 0,
        })
    solution, verified = None, False
    if first_rank is not None:
        before = meter.snapshot()
        solution = tracker.solve()
        projected_rep = E.mul(points[0], H)
        verified = (E.mul(G, solution[0]) == projected_rep
                    and E.mul(G, solution[1]) == Q)
        phase_costs["recovery_verify"] = meter.delta(before, meter.snapshot())
        assert verified, "full-rank system recovered wrong group labels"
    else:
        phase_costs["recovery_verify"] = {}
    phase_costs["scan_all"] = pilot.sum_costs(
        *(r["scan_cost"] for r in records))
    phase_costs["scan_to_rank"] = (
        pilot.sum_costs(*(r["scan_cost"] for r in records[:first_rank]))
        if first_rank is not None else None)
    return {
        "base_kind": base_kind,
        "policy": "affine_screen" if screen else "baseline",
        "base_size": len(points),
        "x_span_rank": len(aff.basis) if aff else len(
            pilot.echelon([P[0] for P in points], Counter())),
        "pair_entries": sum(map(len, table.values())),
        "distinct_pair_sums": len(table),
        "screen_cache_entries": len(aff.cache) if aff else 0,
        "screen_row_xor": dict(aff.work) if aff else {},
        "skipped_residuals": skipped,
        "first_full_rank_attempt": first_rank, "rank": len(tracker.pivots),
        "solution": solution, "verified": verified,
        "rank_ops": dict(tracker.ops),
        "phase_costs": phase_costs, "records": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "preserve previous evidence"
    costs, meter = {}, prior.Meter()

    def setup_field():
        field = prior.CountedField(meter)
        field.frobenius(1, N - 1)
        return field

    F = pilot.phase(meter, "field_setup", setup_field, costs)
    E0 = prior.CountedCurve(meter, F, 0, 1)
    twist = prior.CountedCurve(meter, F, 1, 1)
    lines, selected, trials = pilot.phase(
        meter, "kernel_line_search",
        lambda: prior.order_seven_lines(F, twist), costs)
    kernel = lines[selected]
    phi = pilot.phase(
        meter, "map_construction",
        lambda: prior.BinaryVeluMap.from_generator(
            E0, twist, kernel["generator"], prior.ELL), costs)
    E1 = prior.CountedCurve(meter, F, phi.codomain.a, phi.codomain.b)
    assert E1.b == kernel["codomain_b"] == 0x11584f
    G = pilot.phase(meter, "generator_search",
                    lambda: prior.choose_generator(E0)[0], costs)
    search_stats, eligible, selected_bases = pilot.phase(
        meter, "raw_orbit_search",
        lambda: raw_bases(E0), costs)
    tau = pilot.phase(meter, "frobenius_scalar",
                      lambda: pilot.frobenius_scalar(E0, G), costs)
    bases = {}
    for kind, record in selected_bases.items():
        points = record["points"]

        def make_weights():
            projected = [E0.mul(P, H) for P in points]
            assert len(set(projected)) == 42
            return pilot.point_weights(
                E0, record["projected_representative"], projected, tau)

        weights = pilot.phase(
            meter, f"projected_weights_{kind}", make_weights, costs)
        images = pilot.phase(
            meter, f"base_transport_{kind}",
            lambda points=points: [phi(P) for P in points], costs)
        assert len(set(images)) == 42
        assert all(P is not None and E1.on_curve(P) for P in images)
        assert all(phi(E0.mul(P, H)) == E1.mul(phi(P), H)
                   for P in points)
        bases[kind] = {
            "record": record, "weights": weights,
            "original": points, "transported": images}
    ranked_keys, held_keys = pilot.phase(
        meter, "target_orbit_partition_audit",
        lambda: pilot.target_partition(E0, G), costs)
    secret, coefficients = target_coefficients()
    Q = pilot.phase(
        meter, "challenge_generation_audit",
        lambda: E0.mul(G, secret), costs)
    assert Q is not None and E0.mul(Q, R) is None
    G1, Q1 = pilot.phase(
        meter, "generator_challenge_transport",
        lambda: (phi(G), phi(Q)), costs)
    assert G1 is not None and E1.mul(G1, secret) == Q1
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
    target_keys = pilot.phase(
        meter, "target_key_audit",
        lambda: [pilot.x_orbit_key(F, T) for _, _, T in targets0], costs)
    assert all(k is None or k in ranked_keys for k in target_keys)
    parity = pilot.phase(
        meter, "target_map_audit",
        lambda: all(phi(A[2]) == B[2]
                    for A, B in zip(targets0, targets1)), costs)
    assert parity
    arms = {}
    for kind in ("train", "held"):
        for geometry, E, gen, challenge, targets in [
            ("original", E0, G, Q, targets0),
            ("transported", E1, G1, Q1, targets1),
        ]:
            for screen in (False, True):
                name = f"{kind}_{geometry}_{'screen' if screen else 'baseline'}"
                arms[name] = raw_algorithm(
                    E, bases[kind][geometry], bases[kind]["weights"],
                    gen, challenge, targets, target_keys, held_keys, kind,
                    screen=screen, meter=meter)
    for kind in ("train", "held"):
        for geometry in ("original", "transported"):
            baseline = arms[f"{kind}_{geometry}_baseline"]
            screened = arms[f"{kind}_{geometry}_screen"]
            assert [(r["first_witness"], r["rank_after"]) for r in baseline["records"]] == [
                (r["first_witness"], r["rank_after"]) for r in screened["records"]]
            assert baseline["solution"] == screened["solution"]
        for policy in ("baseline", "screen"):
            source = arms[f"{kind}_original_{policy}"]
            mapped = arms[f"{kind}_transported_{policy}"]
            assert [(r["first_witness"], r["rank_after"]) for r in source["records"]] == [
                (r["first_witness"], r["rank_after"]) for r in mapped["records"]]
            assert source["solution"] == mapped["solution"]
    summary = pilot.summarize(arms)
    cold = {}
    for name, arm in arms.items():
        kind, geometry, mode = name.split("_")
        prefix = arm["first_full_rank_attempt"]
        if prefix is None:
            cold[name] = None
            continue
        components = [
            costs["field_setup"], costs["generator_search"],
            costs["raw_orbit_search"], costs["frobenius_scalar"],
            costs[f"projected_weights_{kind}"],
            costs["challenge_generation_audit"],
            arm["phase_costs"]["pair_table"],
            arm["phase_costs"]["scan_to_rank"],
            arm["phase_costs"]["recovery_verify"],
        ]
        if geometry == "transported":
            components += [
                costs["kernel_line_search"], costs["map_construction"],
                costs[f"base_transport_{kind}"],
                costs["generator_challenge_transport"]]
            components += codomain_target_costs[:prefix]
        else:
            components += source_target_costs[:prefix]
        if mode == "screen":
            components.append(arm["phase_costs"]["screen_setup"])
        total = pilot.sum_costs(*components)
        total["row_xor"] = sum(
            r["scan_row_xor"] for r in arm["records"][:prefix])
        total["pair_lookups"] = sum(
            r["pair_lookups"] for r in arm["records"][:prefix])
        cold[name] = total
    artifact = {
        "schema": "ecc2k130-pdp-raw-orbit-holdout-v1",
        "status": "TOY_DIAGNOSTIC",
        "parameters": {
            "seed": SEED, "degree": N, "irreducible": prior.IRR,
            "subgroup_order": R, "cofactor": H,
            "source_order": prior.ORDER, "attempts": ATTEMPTS,
            "isogeny_degree": prior.ELL, "raw_x_trials": 20_000},
        "source_sha256": {
            "protocol": pilot.digest(Path(__file__).with_name("PROTOCOL_RAW.md")),
            "runner": pilot.digest(Path(__file__)),
            "projected_runner": pilot.digest(Path(__file__).with_name("run.py")),
            "prior_degree7_runner": pilot.digest(
                ROOT / "research/ecc2k130_factor_base_pilot_20260924/run.py"),
            "fastfield": pilot.digest(
                ROOT / "research/ecc2k130_relations/fastfield.py"),
            "relations": pilot.digest(
                ROOT / "research/ecc2k130_relations/relations.py"),
            "oriented_velu": pilot.digest(
                ROOT / "research/ecc2k130_oriented_transport_20260924/oriented_velu.py"),
        },
        "geometry": {
            "kernel_line_search_trials": trials,
            "selected_kernel_x": list(selected),
            "kernel_generator": pilot.enc(kernel["generator"]),
            "codomain_b": E1.b, "generator": pilot.enc(G),
            "transported_generator": pilot.enc(G1),
            "challenge": pilot.enc(Q), "transported_challenge": pilot.enc(Q1),
            "secret_audit_only": secret,
            "frobenius_scalar_mod_r": tau,
            "raw_search_counts": dict(search_stats),
            "eligible_orbits": [
                {"trial": c["trial"], "representative": pilot.enc(c["representative"]),
                 "projected_representative": pilot.enc(c["projected_representative"]),
                 "x_span_rank": c["x_span_rank"], "hash": c["hash"]}
                for c in eligible],
            "bases": {
                kind: {
                    "trial": b["record"]["trial"],
                    "representative": pilot.enc(b["record"]["representative"]),
                    "projected_representative": pilot.enc(
                        b["record"]["projected_representative"]),
                    "hash": b["record"]["hash"], "weights": b["weights"],
                    "original": [pilot.enc(P) for P in b["original"]],
                    "transported": [pilot.enc(P) for P in b["transported"]]}
                for kind, b in bases.items()},
            "all_target_orbits": ranked_keys,
            "held_target_orbits": sorted(held_keys),
        },
        "target_coefficients": [list(x) for x in coefficients],
        "source_targets": [pilot.enc(T) for _, _, T in targets0],
        "transported_targets": [pilot.enc(T) for _, _, T in targets1],
        "target_keys": target_keys,
        "setup_costs": costs,
        "source_target_costs": source_target_costs,
        "codomain_target_costs": codomain_target_costs,
        "arms": arms, "summary": summary,
        "cold_cost_to_verified_rank": cold,
        "limitations": [
            "toy order-421 challenge subgroup and degree-7 isogeny",
            "raw full-group base points require cofactor-projected log labels",
            "transported source orbit is not a native codomain Frobenius orbit",
            "exact pair-table lookup is the complete finite-base comparator",
            "no degree-131 or rho crossover inference"],
        "platform": platform.platform(), "python": sys.version,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.out),
        "base_trials": {k: b["record"]["trial"] for k, b in bases.items()},
        "base_x_span_ranks": {
            k: {g: arms[f"{k}_{g}_baseline"]["x_span_rank"]
                for g in ("original", "transported")}
            for k in ("train", "held")},
        "hits": {k: sum(r["hit"] for r in a["records"])
                 for k, a in arms.items()},
        "first_rank": {k: a["first_full_rank_attempt"] for k, a in arms.items()},
        "certified_skips": {
            k: sum(r["certified_skips"] for r in a["records"])
            for k, a in arms.items()},
    }, indent=2))


if __name__ == "__main__":
    main()
