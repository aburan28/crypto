#!/usr/bin/env python3
"""Preregistered disjoint-orbit degree-7 factor-base replication."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_factor_base_pilot_20260924"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_dual_transport_20260925"))
import run as pilot  # noqa: E402
from dual_transport import DualTransport  # noqa: E402

SEEDS = (2026092511, 2026092512)
HOLDOUTS = ("A", "B")
VARIANTS = ("original", "transported", "descendant_native", "pullback")
MAX_BASE_X = 4096
MAX_TARGET_DRAWS = 8192


def hash_int(s: str) -> int:
    return int.from_bytes(hashlib.sha256(s.encode()).digest(), "big")


def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def phase(meter, ledger, name, fn):
    before = meter.snapshot()
    result = fn()
    ledger[name] = meter.delta(before, meter.snapshot())
    return result


def canonical(E, point):
    if point is None:
        return None
    choices = []
    q = point
    for _ in range(pilot.N):
        choices.extend((q, E.neg(q)))
        q = (E.F.sqr(q[0]), E.F.sqr(q[1]))
    assert q == point
    return min(choices)


def orbit_universe(E, G):
    by_point, members = {}, {}
    q = None
    for _ in range(pilot.R):
        if q is not None:
            representative = canonical(E, q)
            by_point[q] = representative
            members.setdefault(representative, []).append(q)
        q = E.add(q, G)
    assert q is None and len(by_point) == pilot.R - 1
    representatives = sorted(members)
    assert len(representatives) == 10
    assert all(len(members[k]) == 42 for k in representatives)
    return representatives, by_point, members


def construct_base(curve, seed, role, by_point, quota, source, dual):
    points, back_points = [], []
    seen_x, seen_point = set(), set()
    accepted = Counter()
    meta = Counter()
    for trial in range(MAX_BASE_X):
        x = hash_int(f"degree7-base-v1|{seed}|{role}|{trial}") & ((1 << pilot.N) - 1)
        if x in seen_x:
            meta["repeated_x"] += 1
            continue
        seen_x.add(x)
        meta["x_scanned"] += 1
        for p in curve.points_over(x):
            meta["raw_projected"] += 1
            q = curve.mul(p, pilot.COFACTOR)
            if q is None:
                meta["infinity"] += 1
                continue
            if q in seen_point:
                meta["duplicate"] += 1
                continue
            seen_point.add(q)
            if role == "source":
                back = q
            else:
                back = source.mul(dual.dual(q), pow(pilot.ELL, -1, pilot.R))
                assert back is not None
                meta["dual_pullbacks"] += 1
            rep = by_point.get(back)
            assert rep is not None, "candidate failed subgroup/orbit audit"
            meta["classified"] += 1
            if rep in quota and accepted[rep] < quota[rep]:
                points.append(q)
                back_points.append(back)
                accepted[rep] += 1
                if len(points) == pilot.BASE_SIZE:
                    meta["candidate_trials"] = trial + 1
                    assert all(accepted[k] == quota[k] for k in quota)
                    return points, back_points, dict(meta), dict(accepted)
            else:
                meta["out_of_quota"] += 1
    meta["candidate_trials"] = MAX_BASE_X
    return points, back_points, dict(meta), dict(accepted)


def make_targets(E, G, Q, by_point, allowed, label, meter):
    accepted, costs = [], []
    draws, rejections = 0, Counter()
    previous = meter.snapshot()
    for trial in range(MAX_TARGET_DRAWS):
        draws += 1
        u = hash_int(f"degree7-target-v1|{label}|{trial}|u") % pilot.R
        v = 1 + hash_int(f"degree7-target-v1|{label}|{trial}|v") % (pilot.R - 1)
        T = E.add(E.mul(G, u), E.mul(Q, v))
        rep = by_point.get(T)
        if rep is None:
            rejections["infinity"] += 1
        elif rep not in allowed:
            rejections["other_orbit"] += 1
        else:
            accepted.append((u, v, T))
            now = meter.snapshot()
            costs.append(meter.delta(previous, now))
            previous = now
            if len(accepted) == pilot.ATTEMPTS:
                break
    return accepted, costs, {"draws": draws, "rejections": dict(rejections),
                             "accepted": len(accepted)}


def codomain_targets(E, G, Q, targets, meter):
    output, costs = [], []
    for u, v, _ in targets:
        before = meter.snapshot()
        output.append((u, v, E.add(E.mul(G, u), E.mul(Q, v))))
        costs.append(meter.delta(before, meter.snapshot()))
    return output, costs


def slim_variant(v):
    return {key: value for key, value in v.items() if key != "cases"} | {
        "cases": [[c["witness_count"], c["first_witness"],
                   c["independent_before_rank_stop"]] for c in v["cases"]]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    assert not args.out.exists(), "receipt must be immutable"
    started = time.monotonic()
    meter, ledger = pilot.Meter(), {}
    before = meter.snapshot()
    F = pilot.CountedField(meter)
    F.frobenius(1, pilot.N - 1)
    ledger["field_setup"] = meter.delta(before, meter.snapshot())
    E0 = pilot.CountedCurve(meter, F, 0, 1)
    Twist = pilot.CountedCurve(meter, F, 1, 1)
    lines, selected, line_trials = phase(
        meter, ledger, "kernel_search", lambda: pilot.order_seven_lines(F, Twist))
    old = json.loads((ROOT / "research/ecc2k130_factor_base_pilot_20260924/results_final.json").read_text())
    assert list(selected) == old["geometry"]["selected_kernel_x"]
    record = lines[selected]
    phi = phase(meter, ledger, "forward_setup", lambda: pilot.BinaryVeluMap.from_generator(
        E0, Twist, record["generator"], pilot.ELL))
    complement = next(lines[k]["generator"] for k in sorted(lines) if k != selected)
    G, gen_trials = phase(meter, ledger, "generator", lambda: pilot.choose_generator(E0))
    secret = hash_int("degree7-replication-secret-v1") % pilot.R
    assert 0 < secret < pilot.R
    Q = phase(meter, ledger, "challenge", lambda: E0.mul(G, secret))
    assert Q is not None and E0.mul(Q, pilot.R) is None
    D = phase(meter, ledger, "dual_setup", lambda: DualTransport(
        E0, Twist, record["generator"], complement, pilot.ELL, G))
    assert (D.codomain.a, D.codomain.b) == (phi.codomain.a, phi.codomain.b)
    assert D.compose(G) == E0.mul(G, pilot.ELL)
    E1 = pilot.CountedCurve(meter, F, phi.codomain.a, phi.codomain.b)
    assert E1.b == record["codomain_b"]
    representatives, by_point, members = phase(
        meter, ledger, "source_orbit_universe", lambda: orbit_universe(E0, G))
    quotas = {k: 4 for k in representatives[:4]}
    bases, base_meta, occupancy = {}, {}, {}
    for seed in SEEDS:
        key = str(seed)
        B0, _, meta0, occ0 = phase(meter, ledger, f"source_base_{key}",
            lambda seed=seed: construct_base(E0, seed, "source", by_point,
                                              quotas, E0, D))
        B1, back, meta1, occ1 = phase(meter, ledger, f"native_base_{key}",
            lambda seed=seed: construct_base(E1, seed, "leaf", by_point,
                                              quotas, E0, D))
        assert len(B0) == len(B1) == pilot.BASE_SIZE
        assert occ0 == occ1 == quotas
        Btransport = phase(meter, ledger, f"transport_base_{key}",
            lambda B0=B0: [phi(p) for p in B0])
        assert len(set(Btransport)) == pilot.BASE_SIZE
        assert len(set(back)) == pilot.BASE_SIZE
        bases[key] = {"original": B0, "transported": Btransport,
                      "descendant_native": B1, "pullback": back}
        base_meta[key] = {"source": meta0, "native": meta1}
        occupancy[key] = {"source": {str(k): occ0[k] for k in quotas},
                          "native": {str(k): occ1[k] for k in quotas}}
    assert bases[str(SEEDS[0])]["original"] != bases[str(SEEDS[1])]["original"]
    assert bases[str(SEEDS[0])]["descendant_native"] != bases[str(SEEDS[1])]["descendant_native"]
    G1, Q1 = phase(meter, ledger, "codomain_generator_challenge",
        lambda: (phi(G), phi(Q)))
    assert E1.mul(G1, secret) == Q1
    target_data, target_cost, target_meta, leaf_target_cost = {}, {}, {}, {}
    for label, allowed in (("A", set(representatives[:5])),
                           ("B", set(representatives[5:]))):
        targets, costs, meta = make_targets(E0, G, Q, by_point, allowed, label, meter)
        assert len(targets) == pilot.ATTEMPTS
        target_data[label], target_cost[label], target_meta[label] = targets, costs, meta
        leaf, leaf_cost = codomain_targets(E1, G1, Q1, targets, meter)
        target_data[f"{label}_leaf"] = leaf
        leaf_target_cost[label] = leaf_cost
        # Full covariance is an audit and is excluded from candidate costs.
        phase(meter, ledger, f"covariance_audit_{label}",
              lambda targets=targets, leaf=leaf: [
                  (_ for _ in ()).throw(AssertionError("map target mismatch"))
                  if phi(src[2]) != dst[2] else None
                  for src, dst in zip(targets, leaf)])
    results, cold, controls = {}, {}, {}
    for seed in SEEDS:
        key = str(seed)
        results[key], cold[key], controls[key] = {}, {}, {}
        for label in HOLDOUTS:
            results[key][label], cold[key][label] = {}, {}
            for variant in VARIANTS:
                leaf = variant in ("transported", "descendant_native")
                curve, gen, challenge = (E1, G1, Q1) if leaf else (E0, G, Q)
                targets = target_data[f"{label}_leaf"] if leaf else target_data[label]
                v = pilot.relation_run(curve, gen, challenge, targets,
                    bases[key][variant], expected_secret=secret, meter=meter)
                results[key][label][variant] = slim_variant(v)
                limit = v["first_full_rank_attempt"] or pilot.ATTEMPTS
                shared = [ledger["field_setup"], ledger["generator"], ledger["challenge"],
                          ledger["source_orbit_universe"],
                          *target_cost[label][:limit]]
                specific = [v["costs"]["pair_table"], v["costs"]["scan_to_rank"],
                            v["costs"]["recovery_verify"]]
                if variant == "original":
                    specific.append(ledger[f"source_base_{key}"])
                elif variant == "transported":
                    specific.extend((ledger["kernel_search"], ledger["forward_setup"],
                                     ledger[f"source_base_{key}"],
                                     ledger[f"transport_base_{key}"],
                                     ledger["codomain_generator_challenge"],
                                     *leaf_target_cost[label][:limit]))
                elif variant == "descendant_native":
                    specific.extend((ledger["kernel_search"], ledger["dual_setup"],
                                     ledger[f"native_base_{key}"],
                                     ledger["codomain_generator_challenge"],
                                     *leaf_target_cost[label][:limit]))
                else:
                    specific.extend((ledger["kernel_search"], ledger["dual_setup"],
                                     ledger[f"native_base_{key}"]))
                cold[key][label][variant] = pilot.add_costs(*shared, *specific)
            a, b, c, d = (results[key][label][name] for name in VARIANTS)
            controls[key][label] = {
                "forward_pair_invariant": [x[0] for x in a["cases"]] == [x[0] for x in b["cases"]],
                "dual_pair_invariant": [x[0] for x in c["cases"]] == [x[0] for x in d["cases"]],
                "forward_rank_equal": a["first_full_rank_attempt"] == b["first_full_rank_attempt"],
                "dual_rank_equal": c["first_full_rank_attempt"] == d["first_full_rank_attempt"],
                "all_verified": all(results[key][label][name]["verified"] for name in VARIANTS),
            }
    sources = {name: digest(path) for name, path in {
        "protocol": Path(__file__).with_name("PROTOCOL.md"), "runner": Path(__file__),
        "pilot": ROOT / "research/ecc2k130_factor_base_pilot_20260924/run.py",
        "pilot_result": ROOT / "research/ecc2k130_factor_base_pilot_20260924/results_final.json",
        "dual": ROOT / "research/ecc2k130_dual_transport_20260925/dual_transport.py",
        "field": ROOT / "research/ecc2k130_relations/fastfield.py",
        "curve": ROOT / "research/ecc2k130_relations/relations.py",
        "velu": ROOT / "research/ecc2k130_oriented_transport_20260924/oriented_velu.py",
    }.items()}
    raw = {"schema": "degree7-disjoint-orbit-factor-base-replication-v1",
           "timestamp_utc": datetime.now(timezone.utc).isoformat(),
           "platform": platform.platform(), "python": sys.version,
           "source_sha256": sources,
           "parameters": {"base_seeds": SEEDS, "holdouts": HOLDOUTS,
                          "r": pilot.R, "cofactor": pilot.COFACTOR,
                          "attempts_per_holdout": pilot.ATTEMPTS},
           "geometry": {"selected_kernel_x": list(selected),
                        "selected_codomain_b": E1.b, "line_trials": line_trials,
                        "generator_trials": gen_trials,
                        "dual_raw_sign": D.raw_reverse_scalar_sign,
                        "orbit_representatives": [pilot.encode_point(k) for k in representatives],
                        "orbit_sizes": [len(members[k]) for k in representatives]},
           "challenge": {"G": pilot.encode_point(G), "Q": pilot.encode_point(Q),
                         "secret_audit_only": secret},
           "bases": {k: {n: [pilot.encode_point(p) for p in v]
                         for n, v in group.items()} for k, group in bases.items()},
           "base_meta": base_meta, "occupancy": occupancy,
           "target_streams": {label: {"coefficients": [[u, v] for u, v, _ in target_data[label]],
                                       "meta": target_meta[label],
                                       "point_orbits": [representatives.index(by_point[T])
                                                        for _, _, T in target_data[label]]}
                              for label in HOLDOUTS},
           "phase_costs": ledger, "target_prefix_costs": target_cost,
           "codomain_target_prefix_costs": leaf_target_cost,
           "variants": results, "cold_cost_to_rank_or_512": cold,
           "controls": controls, "wall_seconds": time.monotonic() - started,
           "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
           "limitations": ["Toy degree-7 ramified map; no n131 PDP or DLP inference",
                           "Signed source tau orbits used only as matched controls",
                           "Audit covariance excluded from candidate costs"]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(raw, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"receipt": str(args.out), "wall_seconds": raw["wall_seconds"],
                      "full_rank": {k: {h: {v: results[k][h][v]["first_full_rank_attempt"]
                                          for v in VARIANTS} for h in HOLDOUTS}
                                    for k in results},
                      "controls": controls}, indent=2))


if __name__ == "__main__":
    main()
