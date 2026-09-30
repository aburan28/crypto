#!/usr/bin/env python3
"""Frozen two-candidate m3 selector; never run before FROZEN.json is merged."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PRIOR = ROOT / "research/notes/ecc2k130/m3_four_policy_20260930/produce.py"
prior_spec = importlib.util.spec_from_file_location("m3_prior_producer", PRIOR)
assert prior_spec is not None and prior_spec.loader is not None
previous = importlib.util.module_from_spec(prior_spec)
prior_spec.loader.exec_module(previous)
from score import choose, score_candidate  # noqa: E402

pilot, rep, DualTransport = previous.pilot, previous.rep, previous.DualTransport
save, sha, peak_rss_bytes = previous.save, previous.sha, previous.peak_rss_bytes
CONFIG, FROZEN = HERE / "CONFIG.json", HERE / "FROZEN.json"
POLICIES = ("original", "transported", "descendant_native", "pullback")
ARMS = ("control", "selected")


def checked_config() -> dict:
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-degree7-m3-base-selector-protocol-v1"
    assert config["source_lock_required_before_fixture_generation"] is True
    assert (config["field_degree"], int(config["field_modulus_hex"], 16),
            config["subgroup_order"], config["cofactor"],
            config["isogeny_degree"], config["base_size"]) == (
                pilot.N, pilot.IRR, pilot.R, pilot.COFACTOR, pilot.ELL, 8)
    lock = json.loads(FROZEN.read_text())
    assert lock["schema"] == "ecc2k130-degree7-m3-base-selector-source-lock-v1"
    assert lock["protocol_sha256"] == sha(HERE / "PROTOCOL.md")
    for relative, digest in lock["sha256"].items():
        assert sha(ROOT / relative) == digest, relative
    return config


class IncompleteCandidate(AssertionError):
    def __init__(self, receipt: dict):
        self.receipt = receipt
        super().__init__(f"incomplete {receipt['role']} candidate "
                         f"{receipt['index']}: {receipt['points']} points")


def construct_candidate(E, seed: int, role: str, index: int, source, dual,
                        by_point: dict, quotas: dict, config: dict):
    points, pullbacks = [], []
    used_x, used_points = set(), set()
    accepted, meta = Counter(), Counter()
    inv_ell = pow(pilot.ELL, -1, pilot.R)
    for trial in range(config["max_base_x_trials_per_candidate"]):
        label = (f"degree7-m3-base-selector-base-v1|{seed}|{role}|"
                 f"{index}|{trial}")
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
                    assert accepted == quotas
                    return points, pullbacks, dict(meta), dict(accepted)
            else:
                meta["out_of_quota"] += 1
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
    raise IncompleteCandidate({"seed": seed, "role": role, "index": index,
        "points": len(points), "base": [list(P) for P in points],
        "pullback": [list(P) for P in pullbacks], "meta": dict(meta),
        "occupancy": {str(k): accepted[k] for k in quotas}})


def make_targets(E, G, Q, by_point, allowed, holdout: str, meter, config):
    targets, costs, rejections = [], [], Counter()
    previous_mark = meter.snapshot()
    for trial in range(config["max_target_draws_per_holdout"]):
        stem = f"degree7-m3-base-selector-target-v1|{holdout}|{trial}"
        u = rep.hash_int(stem + "|u") % pilot.R
        v = 1 + rep.hash_int(stem + "|v") % (pilot.R - 1)
        T = E.add(E.mul(G, u), E.mul(Q, v))
        orbit = by_point.get(T)
        if orbit is None:
            rejections["infinity"] += 1
        elif orbit not in allowed:
            rejections["other_orbit"] += 1
        else:
            targets.append((u, v, T))
            now = meter.snapshot()
            costs.append(meter.delta(previous_mark, now))
            previous_mark = now
            if len(targets) == config["accepted_relation_targets_per_holdout"]:
                return targets, costs, {"draws": trial + 1,
                    "accepted": len(targets), "rejections": dict(rejections)}
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
    raise AssertionError(f"incomplete target holdout {holdout}")


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
    lines, selected, line_trials = rep.phase(meter, ledger, "kernel_search",
        lambda: pilot.order_seven_lines(F, twist))
    prior = json.loads((ROOT / "research/ecc2k130_factor_base_pilot_20260924/"
                        "results_final.json").read_text())
    assert list(selected) == prior["geometry"]["selected_kernel_x"]
    kernel = lines[selected]
    phi = rep.phase(meter, ledger, "forward_setup", lambda:
        pilot.BinaryVeluMap.from_generator(E0, twist, kernel["generator"], pilot.ELL))
    complement = next(lines[k]["generator"] for k in sorted(lines) if k != selected)
    G, generator_trials = rep.phase(meter, ledger, "generator", lambda:
        pilot.choose_generator(E0))
    D = rep.phase(meter, ledger, "dual_setup", lambda:
        DualTransport(E0, twist, kernel["generator"], complement, pilot.ELL, G))
    rep.phase(meter, ledger, "dual_identity_audit", lambda:
        (D.compose(G) == E0.mul(G, pilot.ELL)) or
        (_ for _ in ()).throw(AssertionError("dual identity")))
    E1 = pilot.CountedCurve(meter, F, phi.codomain.a, phi.codomain.b)
    assert (E1.a, E1.b) == (D.codomain.a, D.codomain.b)
    assert E1.b == kernel["codomain_b"]
    reps, by_point, members = rep.phase(meter, ledger, "source_orbit_universe",
        lambda: rep.orbit_universe(E0, G))
    assert len(by_point) == 420 and len(reps) == 10
    quotas = {orbit: config["source_signed_orbit_quota"] for orbit in reps[:4]}

    # Both scores and the choice are completed before constructing the new Q.
    candidates, choices, scores = {}, {}, {}
    for seed in config["base_seeds"]:
        key = str(seed)
        candidates[key], choices[key], scores[key] = {}, {}, {}
        for role, E in (("source", E0), ("leaf", E1)):
            candidates[key][role], scores[key][role] = [], []
            for index in range(config["candidates_per_role_seed"]):
                name = f"{role}_{key}_{index}"
                try:
                    base, back, meta, occupancy = rep.phase(meter, ledger,
                        f"candidate_{name}", lambda E=E, seed=seed,
                        role=role, index=index: construct_candidate(
                            E, seed, role, index, E0, D, by_point, quotas, config))
                except IncompleteCandidate as error:
                    save(output / "candidate_failures" / f"{name}.json", error.receipt)
                    raise
                score = rep.phase(meter, ledger, f"score_{name}", lambda back=back:
                    score_candidate(E0, back, set(by_point), meter,
                                    pilot.RankTracker, pilot.R))
                assert score["score_cost"] == {k: v for k, v in
                    ledger[f"score_{name}"].items() if k != "cpu_ns"}
                record = {"base": [pilot.encode_point(P) for P in base],
                          "pullback": [pilot.encode_point(P) for P in back],
                          "meta": meta,
                          "occupancy": {str(k): occupancy[k] for k in quotas},
                          "score": score}
                candidates[key][role].append(record)
                scores[key][role].append(score)
            choice = choose(scores[key][role])
            choices[key][role] = choice
            if choice is None:
                save(output / "selector_failures" / f"{key}_{role}.json",
                     {"seed": seed, "role": role, "scores": scores[key][role]})
                raise AssertionError(f"no rank-eight {role} candidate for {key}")

    secret = 1 + rep.hash_int(config["labels"]["secret"]) % (pilot.R - 1)
    Q = rep.phase(meter, ledger, "challenge", lambda: E0.mul(G, secret))
    old_result = json.loads((ROOT / "research/notes/ecc2k130/"
        "m3_four_policy_20260930/evidence_run_36722040881/result.json").read_text())
    if list(Q) == old_result["challenge"]["Q"]:
        save(output / "collision.json", {"new_Q": list(Q),
             "prior_Q": old_result["challenge"]["Q"]})
        raise AssertionError("new Q collides with prior panel")
    rep.phase(meter, ledger, "challenge_order_audit", lambda:
        (Q is not None and E0.mul(Q, pilot.R) is None) or
        (_ for _ in ()).throw(AssertionError("challenge order")))
    G1, Q1 = rep.phase(meter, ledger, "codomain_generator_challenge",
                       lambda: (phi(G), phi(Q)))
    assert E1.mul(G1, secret) == Q1

    target_data, target_cost, target_meta, leaf_target_cost = {}, {}, {}, {}
    for label, allowed in (("A", set(reps[:5])), ("B", set(reps[5:]))):
        targets, costs, meta = make_targets(E0, G, Q, by_point, allowed,
                                            label, meter, config)
        target_data[label], target_cost[label], target_meta[label] = targets, costs, meta
        leaf, leaf_cost = rep.codomain_targets(E1, G1, Q1, targets, meter)
        target_data[f"{label}_leaf"], leaf_target_cost[label] = leaf, leaf_cost
        assert all(phi(source[2]) == mapped[2]
                   for source, mapped in zip(targets, leaf))

    bases, base_map_cost = {}, {}
    for key in candidates:
        bases[key], base_map_cost[key] = {}, {}
        for arm in ARMS:
            source_index = 0 if arm == "control" else choices[key]["source"]
            leaf_index = 0 if arm == "control" else choices[key]["leaf"]
            source_base = [tuple(P) for P in candidates[key]["source"][source_index]["base"]]
            native_base = [tuple(P) for P in candidates[key]["leaf"][leaf_index]["base"]]
            pullback = [tuple(P) for P in candidates[key]["leaf"][leaf_index]["pullback"]]
            transport = rep.phase(meter, ledger, f"transport_{key}_{arm}",
                                  lambda source_base=source_base:
                                  [phi(P) for P in source_base])
            base_map_cost[key][arm] = ledger[f"transport_{key}_{arm}"]
            assert all(phi(P) == R for P, R in zip(pullback, native_base))
            bases[key][arm] = {"original": source_base,
                "transported": transport, "descendant_native": native_base,
                "pullback": pullback}

    cells, cold, cold_mod_r, cell_timings, cell_hashes, controls = (
        {}, {}, {}, {}, {}, {})
    for key in candidates:
        cells[key], cold[key], cold_mod_r[key], cell_timings[key] = {}, {}, {}, {}
        cell_hashes[key], controls[key] = {}, {}
        for label in config["holdouts"]:
            cells[key][label], cold[key][label] = {}, {}
            cold_mod_r[key][label], cell_timings[key][label] = {}, {}
            cell_hashes[key][label], controls[key][label] = {}, {}
            for arm in ARMS:
                cells[key][label][arm], cold[key][label][arm] = {}, {}
                cold_mod_r[key][label][arm], cell_timings[key][label][arm] = {}, {}
                cell_hashes[key][label][arm] = {}
                trajectories = {}
                for policy in POLICIES:
                    leaf = policy in ("transported", "descendant_native")
                    E, gen, challenge = (E1, G1, Q1) if leaf else (E0, G, Q)
                    targets = target_data[f"{label}_leaf"] if leaf else target_data[label]
                    cell_started, cell_cpu_started = time.perf_counter(), time.process_time()
                    cell = previous.relation_run(E, gen, challenge, targets,
                        bases[key][arm][policy], meter, secret)
                    cell_timings[key][label][arm][policy] = {
                        "wall_seconds_including_audit_tail": time.perf_counter() - cell_started,
                        "cpu_seconds_including_audit_tail": time.process_time() - cell_cpu_started}
                    cells[key][label][arm][policy] = {
                        field: value for field, value in cell.items() if field != "cases"}
                    trajectories[policy] = [(case["first_witness"], case["rank_after"])
                                            for case in cell["cases"]]
                    path = output / "cells" / key / label / arm / f"{policy}.json"
                    save(path, cell)
                    cell_hashes[key][label][arm][policy] = sha(path)
                    limit = cell["first_full_rank_attempt"] or len(targets)
                    shared = [ledger[name] for name in ("field_setup", "generator",
                               "challenge", "source_orbit_universe")]
                    shared += target_cost[label][:limit]
                    specific = [cell["costs"][name] for name in (
                        "pair_table", "scan_to_rank", "recovery_verify")]
                    role = "source" if policy in ("original", "transported") else "leaf"
                    chosen_index = 0 if arm == "control" else choices[key][role]
                    candidate_indices = (0,) if arm == "control" else (0, 1)
                    specific += [ledger[f"candidate_{role}_{key}_{i}"]
                                 for i in candidate_indices]
                    if arm == "selected":
                        specific += [ledger[f"score_{role}_{key}_{i}"] for i in (0, 1)]
                    if policy == "transported":
                        specific += [ledger[name] for name in ("kernel_search",
                            "forward_setup", "codomain_generator_challenge")]
                        specific += [base_map_cost[key][arm],
                                     *leaf_target_cost[label][:limit]]
                    elif policy == "descendant_native":
                        specific += [ledger[name] for name in ("kernel_search",
                            "forward_setup", "dual_setup", "codomain_generator_challenge")]
                        specific += leaf_target_cost[label][:limit]
                    elif policy == "pullback":
                        specific += [ledger[name] for name in (
                            "kernel_search", "forward_setup", "dual_setup")]
                    cold[key][label][arm][policy] = pilot.add_costs(*shared, *specific)
                    row_ops = Counter(cell["mod_r_ops"])
                    if arm == "selected":
                        for i in (0, 1):
                            row_ops.update(scores[key][role][i]["score_mod_r_ops"])
                    cold_mod_r[key][label][arm][policy] = dict(sorted(row_ops.items()))
                    assert chosen_index in (0, 1)
                    assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
                a, b, c, d = (cells[key][label][arm][name] for name in POLICIES)
                controls[key][label][arm] = {
                    "forward_hit_parity": a["hits"] == b["hits"],
                    "dual_hit_parity": c["hits"] == d["hits"],
                    "forward_rank_parity": a["first_full_rank_attempt"] == b["first_full_rank_attempt"],
                    "dual_rank_parity": c["first_full_rank_attempt"] == d["first_full_rank_attempt"],
                    "forward_trajectory_parity": trajectories["original"] == trajectories["transported"],
                    "dual_trajectory_parity": trajectories["descendant_native"] == trajectories["pullback"],
                    "all_verified": all(cells[key][label][arm][name]["verified"]
                                        for name in POLICIES)}
                assert all(value for name, value in controls[key][label][arm].items()
                           if name != "all_verified")

    status = "PASS_PANEL" if all(cells[key][label][arm][policy]["verified"]
        for key in cells for label in cells[key] for arm in ARMS
        for policy in POLICIES) else "NEGATIVE_PANEL"
    result = {"schema": "ecc2k130-degree7-m3-base-selector-result-v1",
        "status": status, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "python": sys.version,
        "host": platform.node(), "command": sys.argv,
        "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                               cwd=ROOT, text=True).strip(),
        "config_sha256": sha(CONFIG), "frozen_sha256": sha(FROZEN),
        "geometry": {"selected_kernel_x": list(selected), "codomain_b": E1.b,
            "line_trials": line_trials, "generator_trials": generator_trials,
            "orbit_representatives": [list(P) for P in reps],
            "orbit_sizes": [len(members[P]) for P in reps]},
        "challenge": {"G": list(G), "Q": list(Q), "secret_audit_only": secret},
        "candidates": candidates, "choices": choices,
        "bases": {key: {arm: {policy: [list(P) for P in base]
            for policy, base in bases[key][arm].items()} for arm in ARMS}
            for key in bases},
        "target_streams": {label: {"coefficients": [[u, v] for u, v, _ in target_data[label]],
            "meta": target_meta[label],
            "point_orbits": [reps.index(by_point[T]) for _, _, T in target_data[label]]}
            for label in config["holdouts"]},
        "phase_costs": ledger, "target_prefix_costs": target_cost,
        "codomain_target_prefix_costs": leaf_target_cost,
        "cells": cells, "cell_sha256": cell_hashes,
        "cold_cost_to_rank_or_512": cold,
        "cold_mod_r_ops_to_rank_or_512": cold_mod_r,
        "cell_timings": cell_timings, "controls": controls,
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
    config = checked_config()
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(
        TimeoutError("frozen producer wall cap")))
    signal.alarm(config["producer_wall_cap_seconds"])
    try:
        result = run(args.out)
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
        print(json.dumps({"status": result["status"], "cells": 64}, sort_keys=True))
    except BaseException as error:
        save(args.out / "failure.json", {"status": "CENSORED" if isinstance(
            error, (TimeoutError, MemoryError, IncompleteCandidate)) else "FAIL",
            "error_type": type(error).__name__, "error": str(error),
            "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
