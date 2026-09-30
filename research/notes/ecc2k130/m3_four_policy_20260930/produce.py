#!/usr/bin/env python3
"""Frozen, charged degree-7 four-policy m3 held-out relation producer."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/ecc2k130_factor_base_pilot_20260924"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_dual_transport_20260925"))
import run as pilot  # noqa: E402
from dual_transport import DualTransport  # noqa: E402

REP_PATH = ROOT / "research/ecc2k130_factor_base_replication_20260925/run.py"
spec = importlib.util.spec_from_file_location("m3_replication_helpers", REP_PATH)
assert spec is not None and spec.loader is not None
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)

CONFIG = HERE / "CONFIG.json"
VARIANTS = ("original", "transported", "descendant_native", "pullback")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def checked_config(require_lock: bool = True) -> dict:
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-degree7-m3-four-policy-v1"
    assert config["field_degree"] == pilot.N == 21
    assert int(config["field_modulus_hex"], 16) == pilot.IRR
    assert config["subgroup_order"] == pilot.R == 421
    assert config["cofactor"] == pilot.COFACTOR == 4988
    assert config["isogeny_degree"] == pilot.ELL == 7
    assert config["base_size"] == 8 and config["arity"] == 3
    for relative, expected in config["inputs_sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    if require_lock:
        lock = json.loads((HERE / "FROZEN.json").read_text())
        assert lock["schema"] == "ecc2k130-degree7-m3-four-policy-source-lock-v1"
        for relative, expected in lock["sha256"].items():
            assert sha(ROOT / relative) == expected, relative
    return config


def construct_base(E, seed: int, role: str, source, dual,
                   by_point: dict, quotas: dict, config: dict):
    points, pullbacks = [], []
    used_x, used_points = set(), set()
    accepted, meta = Counter(), Counter()
    inv_ell = pow(pilot.ELL, -1, pilot.R)
    for trial in range(config["max_base_x_trials_per_curve_seed"]):
        label = f"degree7-m3-four-policy-base-v1|{seed}|{role}|{trial}"
        x = rep.hash_int(label) & ((1 << pilot.N) - 1)
        if x in used_x:
            meta["repeated_x"] += 1
            continue
        used_x.add(x)
        meta["x_scanned"] += 1
        for point in E.points_over(x):
            meta["raw_projected"] += 1
            projected = E.mul(point, pilot.COFACTOR)
            if projected is None:
                meta["infinity"] += 1
                continue
            if projected in used_points:
                meta["duplicate"] += 1
                continue
            used_points.add(projected)
            back = (projected if role == "source" else
                    source.mul(dual.dual(projected), inv_ell))
            assert back is not None
            if role == "leaf":
                meta["dual_pullbacks"] += 1
            orbit = by_point.get(back)
            assert orbit is not None, "candidate failed source subgroup audit"
            meta["classified"] += 1
            if orbit in quotas and accepted[orbit] < quotas[orbit]:
                points.append(projected)
                pullbacks.append(back)
                accepted[orbit] += 1
                if len(points) == config["base_size"]:
                    meta["candidate_trials"] = trial + 1
                    assert all(accepted[k] == quotas[k] for k in quotas)
                    return points, pullbacks, dict(meta), dict(accepted)
            else:
                meta["out_of_quota"] += 1
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
    raise AssertionError(f"incomplete {role} base: {len(points)} points after cap")


def make_targets(E, G, Q, by_point, allowed, label: str, meter,
                 config: dict):
    accepted, costs, rejections = [], [], Counter()
    previous = meter.snapshot()
    for trial in range(config["max_target_draws_per_holdout"]):
        u = rep.hash_int(f"degree7-m3-four-policy-target-v1|{label}|{trial}|u") % pilot.R
        v = 1 + rep.hash_int(
            f"degree7-m3-four-policy-target-v1|{label}|{trial}|v") % (pilot.R - 1)
        T = E.add(E.mul(G, u), E.mul(Q, v))
        orbit = by_point.get(T)
        if orbit is None:
            rejections["infinity"] += 1
        elif orbit not in allowed:
            rejections["other_orbit"] += 1
        else:
            accepted.append((u, v, T))
            now = meter.snapshot()
            costs.append(meter.delta(previous, now))
            previous = now
            if len(accepted) == config["accepted_targets_per_holdout"]:
                return accepted, costs, {
                    "draws": trial + 1, "accepted": len(accepted),
                    "rejections": dict(rejections)}
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
    raise AssertionError(f"incomplete target holdout {label}")


def relation_run(E, G, Q, targets, base, meter, secret: int) -> dict:
    size = len(base)
    assert size == 8
    before_table = meter.snapshot()
    table = pilot.pair_table(E, base)
    after_table = meter.snapshot()
    tracker = pilot.RankTracker(width=size + 1, modulus=pilot.R)
    cases = []
    first_rank = None
    solution = None
    verified = False
    at_rank = after_verify = None
    hits = dependent = 0
    for index, (u, v, T) in enumerate(targets):
        before_target = meter.snapshot()
        witness, multiplicity, tested = None, 0, 0
        for k, P3 in enumerate(base):
            tested += 1
            residual = E.add(T, E.neg(P3))
            pairs = table.get(residual, ())
            if pairs:
                a, b = pairs[0]
                witness = [a, b, k]
                multiplicity = len(pairs)
                break
        if witness is not None:
            hits += 1
            assert E.add(E.add(base[witness[0]], base[witness[1]]),
                         base[witness[2]]) == T
        independent = False
        if witness is not None and first_rank is None:
            row = [0] * (size + 1)
            for j in witness:
                row[j] += 1
            row[-1] = -v % pilot.R
            independent = tracker.add(row, u)
            if not independent:
                dependent += 1
            if len(tracker.pivots) == size + 1:
                first_rank = index + 1
        after_target_scan = meter.snapshot()
        if first_rank == index + 1:
            at_rank = after_target_scan
            solution = tracker.solve()
            verified = (solution[-1] == secret and
                        E.mul(G, solution[-1]) == Q and
                        all(E.mul(G, d) == P
                            for d, P in zip(solution[:-1], base)))
            after_verify = meter.snapshot()
        cases.append({"i": index, "u": u, "v": v,
                      "target": pilot.encode_point(T),
                      "first_witness": witness,
                      "selected_pair_multiplicity": multiplicity,
                      "third_candidates_tested": tested,
                      "scan_cost": meter.delta(before_target, after_target_scan),
                      "independent_before_rank_stop": independent,
                      "rank_after": len(tracker.pivots)})
    after_scan = meter.snapshot()
    if first_rank is None:
        at_rank = after_scan
        after_verify = at_rank
    assert at_rank is not None and after_verify is not None
    costs = {"pair_table": meter.delta(before_table, after_table),
             "scan_to_rank": meter.delta(after_table, at_rank),
             "recovery_verify": meter.delta(at_rank, after_verify),
             "audit_tail": meter.delta(after_verify, after_scan)}
    return {"attempts": len(targets), "hits": hits,
            "dependent_before_stop": dependent,
            "first_full_rank_attempt": first_rank,
            "rank": len(tracker.pivots), "solution": solution,
            "verified": bool(verified), "mod_r_ops": dict(tracker.ops),
            "pair_entries": sum(map(len, table.values())),
            "distinct_pair_sums": len(table),
            "costs": costs, "cases": cases}


def triple_support_audit(E, base) -> int:
    supported = set()
    for i in range(len(base)):
        for j in range(i, len(base)):
            pair = E.add(base[i], base[j])
            for k in range(j, len(base)):
                supported.add(E.add(pair, base[k]))
    assert len(supported) <= 120
    return len(supported)


def run(output: Path) -> dict:
    config = checked_config()
    started, cpu_started = time.perf_counter(), time.process_time()
    meter, ledger = pilot.Meter(), {}
    def make_field():
        field = pilot.CountedField(meter)
        field.frobenius(1, pilot.N - 1)
        return field
    F = rep.phase(meter, ledger, "field_setup", make_field)
    E0 = pilot.CountedCurve(meter, F, 0, 1)
    twist = pilot.CountedCurve(meter, F, 1, 1)
    lines, selected, line_trials = rep.phase(
        meter, ledger, "kernel_search", lambda: pilot.order_seven_lines(F, twist))
    prior = json.loads((ROOT / "research/ecc2k130_factor_base_pilot_20260924/results_final.json").read_text())
    assert list(selected) == prior["geometry"]["selected_kernel_x"]
    record = lines[selected]
    phi = rep.phase(meter, ledger, "forward_setup", lambda:
        pilot.BinaryVeluMap.from_generator(E0, twist, record["generator"], pilot.ELL))
    complement = next(lines[k]["generator"] for k in sorted(lines) if k != selected)
    G, generator_trials = rep.phase(meter, ledger, "generator", lambda:
        pilot.choose_generator(E0))
    secret = 1 + rep.hash_int(config["labels"]["secret"]) % (pilot.R - 1)
    Q = rep.phase(meter, ledger, "challenge", lambda: E0.mul(G, secret))
    rep.phase(meter, ledger, "challenge_order_audit", lambda:
              (Q is not None and E0.mul(Q, pilot.R) is None) or
              (_ for _ in ()).throw(AssertionError("challenge order")))
    D = rep.phase(meter, ledger, "dual_setup", lambda:
        DualTransport(E0, twist, record["generator"], complement, pilot.ELL, G))
    rep.phase(meter, ledger, "dual_identity_audit", lambda:
              (D.compose(G) == E0.mul(G, pilot.ELL)) or
              (_ for _ in ()).throw(AssertionError("dual identity")))
    E1 = pilot.CountedCurve(meter, F, phi.codomain.a, phi.codomain.b)
    assert (E1.a, E1.b) == (D.codomain.a, D.codomain.b)
    assert E1.b == record["codomain_b"]
    reps, by_point, members = rep.phase(meter, ledger, "source_orbit_universe",
                                       lambda: rep.orbit_universe(E0, G))
    quotas = {orbit: config["base_points_per_source_signed_orbit"]
              for orbit in reps[:4]}
    bases, base_meta, occupancy = {}, {}, {}
    for seed in config["base_seeds"]:
        key = str(seed)
        B0, _, meta0, occ0 = rep.phase(meter, ledger, f"source_base_{key}",
            lambda seed=seed: construct_base(E0, seed, "source", E0, D,
                                              by_point, quotas, config))
        B1, back, meta1, occ1 = rep.phase(meter, ledger, f"native_base_{key}",
            lambda seed=seed: construct_base(E1, seed, "leaf", E0, D,
                                              by_point, quotas, config))
        assert len(B0) == len(B1) == config["base_size"]
        assert occ0 == occ1 == quotas
        transport = rep.phase(meter, ledger, f"transport_base_{key}",
                              lambda B0=B0: [phi(P) for P in B0])
        bases[key] = {"original": B0, "transported": transport,
                      "descendant_native": B1, "pullback": back}
        base_meta[key] = {"source": meta0, "native": meta1}
        occupancy[key] = {"source": {str(k): occ0[k] for k in quotas},
                          "native": {str(k): occ1[k] for k in quotas}}
        assert all(len(set(points)) == config["base_size"]
                   for points in bases[key].values())
        rep.phase(meter, ledger, f"base_covariance_audit_{key}",
                  lambda B1=B1, back=back: [
                      (_ for _ in ()).throw(AssertionError("base map mismatch"))
                      if phi(src) != dst else None
                      for src, dst in zip(back, B1)])
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
    assert bases[str(config["base_seeds"][0])]["original"] != bases[
        str(config["base_seeds"][1])]["original"]
    assert bases[str(config["base_seeds"][0])]["descendant_native"] != bases[
        str(config["base_seeds"][1])]["descendant_native"]
    G1, Q1 = rep.phase(meter, ledger, "codomain_generator_challenge",
                       lambda: (phi(G), phi(Q)))
    rep.phase(meter, ledger, "codomain_challenge_audit", lambda:
              (E1.mul(G1, secret) == Q1) or
              (_ for _ in ()).throw(AssertionError("codomain challenge")))
    target_data, target_cost, target_meta, leaf_target_cost = {}, {}, {}, {}
    for label, allowed in (("A", set(reps[:5])), ("B", set(reps[5:]))):
        targets, costs, meta = make_targets(E0, G, Q, by_point, allowed,
                                            label, meter, config)
        target_data[label], target_cost[label], target_meta[label] = targets, costs, meta
        leaf, leaf_cost = rep.codomain_targets(E1, G1, Q1, targets, meter)
        target_data[f"{label}_leaf"] = leaf
        leaf_target_cost[label] = leaf_cost
        rep.phase(meter, ledger, f"target_covariance_audit_{label}",
                  lambda targets=targets, leaf=leaf: [
                      (_ for _ in ()).throw(AssertionError("target map mismatch"))
                      if phi(src[2]) != dst[2] else None
                      for src, dst in zip(targets, leaf)])
    variants, cold, support, controls, cell_hashes = {}, {}, {}, {}, {}
    for seed in config["base_seeds"]:
        key = str(seed)
        variants[key], cold[key], support[key], controls[key], cell_hashes[key] = {}, {}, {}, {}, {}
        for label in config["holdouts"]:
            variants[key][label], cold[key][label], support[key][label], cell_hashes[key][label] = {}, {}, {}, {}
            trajectories = {}
            for name in VARIANTS:
                leaf = name in ("transported", "descendant_native")
                E, gen, challenge = (E1, G1, Q1) if leaf else (E0, G, Q)
                targets = target_data[f"{label}_leaf"] if leaf else target_data[label]
                result = relation_run(E, gen, challenge, targets,
                                      bases[key][name], meter, secret)
                variants[key][label][name] = {
                    k: v for k, v in result.items() if k != "cases"}
                trajectories[name] = [
                    (case["first_witness"], case["rank_after"])
                    for case in result["cases"]]
                path = output / "cells" / key / label / f"{name}.json"
                save(path, result)
                cell_hashes[key][label][name] = sha(path)
                support[key][label][name] = rep.phase(meter, ledger,
                    f"support_audit_{key}_{label}_{name}",
                    lambda E=E, base=bases[key][name]: triple_support_audit(E, base))
                limit = result["first_full_rank_attempt"] or len(targets)
                shared = [ledger["field_setup"], ledger["generator"],
                          ledger["challenge"], ledger["source_orbit_universe"],
                          *target_cost[label][:limit]]
                specific = [result["costs"]["pair_table"],
                            result["costs"]["scan_to_rank"],
                            result["costs"]["recovery_verify"]]
                if name == "original":
                    specific.append(ledger[f"source_base_{key}"])
                elif name == "transported":
                    specific.extend((ledger["kernel_search"], ledger["forward_setup"],
                                     ledger[f"source_base_{key}"],
                                     ledger[f"transport_base_{key}"],
                                     ledger["codomain_generator_challenge"],
                                     *leaf_target_cost[label][:limit]))
                elif name == "descendant_native":
                    specific.extend((ledger["kernel_search"], ledger["dual_setup"],
                                     ledger[f"native_base_{key}"],
                                     ledger["codomain_generator_challenge"],
                                     *leaf_target_cost[label][:limit]))
                else:
                    specific.extend((ledger["kernel_search"], ledger["dual_setup"],
                                     ledger[f"native_base_{key}"]))
                cold[key][label][name] = pilot.add_costs(*shared, *specific)
                assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
            a, b, c, d = (variants[key][label][name] for name in VARIANTS)
            controls[key][label] = {
                "forward_hit_parity": a["hits"] == b["hits"],
                "dual_hit_parity": c["hits"] == d["hits"],
                "forward_rank_parity": a["first_full_rank_attempt"] == b["first_full_rank_attempt"],
                "dual_rank_parity": c["first_full_rank_attempt"] == d["first_full_rank_attempt"],
                "forward_trajectory_parity": trajectories["original"] == trajectories["transported"],
                "dual_trajectory_parity": trajectories["descendant_native"] == trajectories["pullback"],
                "all_verified": all(variants[key][label][name]["verified"] for name in VARIANTS)}
            assert all(value for name, value in controls[key][label].items()
                       if name != "all_verified"), controls[key][label]
    status = "PASS_PANEL" if all(variants[key][label][name]["verified"]
        for key in variants for label in variants[key] for name in VARIANTS) else "NEGATIVE_PANEL"
    result = {"schema": "ecc2k130-degree7-m3-four-policy-result-v1",
              "status": status, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "platform": platform.platform(), "python": sys.version,
              "source_head": subprocess.check_output(
                  ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "config_sha256": sha(CONFIG), "frozen_sha256": sha(HERE / "FROZEN.json"),
              "field_degree": pilot.N, "subgroup_order": pilot.R,
              "geometry": {"selected_kernel_x": list(selected),
                           "codomain_b": E1.b, "line_trials": line_trials,
                           "generator_trials": generator_trials,
                           "orbit_representatives": [pilot.encode_point(k) for k in reps],
                           "orbit_sizes": [len(members[k]) for k in reps]},
              "challenge": {"G": pilot.encode_point(G), "Q": pilot.encode_point(Q),
                            "secret_audit_only": secret},
              "bases": {key: {name: [pilot.encode_point(P) for P in points]
                              for name, points in group.items()}
                        for key, group in bases.items()},
              "base_meta": base_meta, "occupancy": occupancy,
              "target_streams": {label: {"coefficients": [[u, v] for u, v, _ in target_data[label]],
                                         "meta": target_meta[label],
                                         "point_orbits": [reps.index(by_point[T]) for _, _, T in target_data[label]]}
                                 for label in config["holdouts"]},
              "phase_costs": ledger, "target_prefix_costs": target_cost,
              "codomain_target_prefix_costs": leaf_target_cost,
              "variants": variants, "cell_sha256": cell_hashes,
              "triple_support_cardinality": support,
              "cold_cost_to_rank_or_512": cold, "controls": controls,
              "wall_seconds": time.perf_counter() - started,
              "cpu_seconds": time.process_time() - cpu_started,
              "peak_rss_bytes": peak_rss_bytes(),
              "PDP_yield_n131": None, "full_ECDLP_cost": None,
              "method_crossover": None}
    save(output / "result.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a frozen panel"
    args.out.mkdir(parents=True)
    config = json.loads(CONFIG.read_text())
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(
        TimeoutError("frozen producer wall cap")))
    signal.alarm(config["producer_wall_cap_seconds"])
    try:
        result = run(args.out)
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
        print(json.dumps({"status": result["status"],
                          "cells": sum(len(v) for seeds in result["variants"].values()
                                       for v in seeds.values())}, sort_keys=True))
    except BaseException as error:
        save(args.out / "failure.json", {
            "status": "CENSORED" if isinstance(error, (TimeoutError, MemoryError)) else "FAIL",
            "error_type": type(error).__name__, "error": str(error),
            "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
