#!/usr/bin/env python3
"""Independent replay of the frozen two-candidate m3 selector evidence."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PRIOR = ROOT / "research/notes/ecc2k130/m3_four_policy_20260930/verify.py"
prior_spec = importlib.util.spec_from_file_location("m3_prior_verifier", PRIOR)
assert prior_spec is not None and prior_spec.loader is not None
reference = importlib.util.module_from_spec(prior_spec)
prior_spec.loader.exec_module(reference)

pilot, DualTransport = reference.pilot, reference.DualTransport
save, sha = reference.save, reference.sha
R = 421
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


def replay_candidate(E, seed, role, index, source, dual, by_point, quotas, config):
    """Regenerate a candidate without the producer's generator or score code."""
    base, back, scanned_x, unique = [], [], set(), set()
    counts, meta = Counter(), Counter()
    for trial in range(config["max_base_x_trials_per_candidate"]):
        label = (f"degree7-m3-base-selector-base-v1|{seed}|{role}|"
                 f"{index}|{trial}")
        x = reference.hash_int(label) & ((1 << 21) - 1)
        if x in scanned_x:
            meta["repeated_x"] += 1
            continue
        scanned_x.add(x)
        meta["x_scanned"] += 1
        for P in E.points_over(x):
            meta["raw_projected"] += 1
            projected = E.mul(P, pilot.COFACTOR)
            if projected is None:
                meta["infinity"] += 1
                continue
            if projected in unique:
                meta["duplicate"] += 1
                continue
            unique.add(projected)
            pulled = (projected if role == "source" else source.mul(
                dual.dual(projected), pow(pilot.ELL, -1, R)))
            assert pulled is not None
            if role == "leaf":
                meta["dual_pullbacks"] += 1
            orbit = by_point.get(pulled)
            assert orbit is not None
            meta["classified"] += 1
            if orbit in quotas and counts[orbit] < quotas[orbit]:
                base.append(projected)
                back.append(pulled)
                counts[orbit] += 1
                if len(base) == 8:
                    meta["candidate_trials"] = trial + 1
                    assert counts == quotas
                    return base, back, dict(meta), {str(k): counts[k] for k in quotas}
            else:
                meta["out_of_quota"] += 1
    raise AssertionError(f"independent {role} candidate {seed}/{index} incomplete")


def independent_score(E, base, universe):
    """Brute all triples, then eliminate rows using the separate verifier."""
    triples, _ = reference.brute_triples(E, base)
    triples.pop(None, None)
    assert set(triples) <= universe
    rows, digest_rows = [], []
    for target in sorted(triples):
        k, i, j = triples[target][0]
        coeffs = [0] * 8
        for index in (i, j, k):
            coeffs[index] += 1
        rows.append((coeffs, 0))
        digest_rows.append([target[0], target[1], k, i, j])
    rank, _ = reference.independent_linear_system(rows, 8)
    return {"distinct_support": len(triples),
        "first_witness_base_row_rank": rank,
        "distinct_first_witness_rows": len({tuple(row) for row, _ in rows}),
        "first_witness_sha256": hashlib.sha256(json.dumps(
            digest_rows, separators=(",", ":")).encode()).hexdigest()}


def independent_choice(scores):
    eligible = [i for i in range(2)
                if scores[i]["first_witness_base_row_rank"] == 8]
    assert eligible, "no rank-eight candidate"
    return sorted(eligible, key=lambda i: (
        -scores[i]["distinct_support"],
        -scores[i]["distinct_first_witness_rows"], i))[0]


def target_stream(E, G, Q, by_point, allowed, label, config):
    targets, rejects = [], Counter()
    for trial in range(config["max_target_draws_per_holdout"]):
        stem = f"degree7-m3-base-selector-target-v1|{label}|{trial}"
        u = reference.hash_int(stem + "|u") % R
        v = 1 + reference.hash_int(stem + "|v") % (R - 1)
        T = E.add(E.mul(G, u), E.mul(Q, v))
        orbit = by_point.get(T)
        if orbit is None:
            rejects["infinity"] += 1
        elif orbit not in allowed:
            rejects["other_orbit"] += 1
        else:
            targets.append((u, v, T))
            if len(targets) == config["accepted_relation_targets_per_holdout"]:
                return targets, {"draws": trial + 1, "accepted": len(targets),
                                 "rejections": dict(rejects)}
    raise AssertionError(f"independent target holdout {label} incomplete")


def sum_costs(parts):
    total = Counter()
    for part in parts:
        assert all(type(value) is int and value >= 0 for value in part.values())
        total.update(part)
    return dict(sorted(total.items()))


def replay_cell(path, expected_sha, summary, E, ref, G, Q, targets, base,
                secret):
    assert sha(path) == expected_sha
    cell = json.loads(path.read_text())
    assert {k: v for k, v in cell.items() if k != "cases"} == summary
    assert len(cell["cases"]) == len(targets) == 512
    triples, pair_sums = reference.brute_triples(E, base)
    assert cell["pair_entries"] == 36
    assert cell["distinct_pair_sums"] == pair_sums
    rows = []
    rank = hits = dependent = 0
    first_rank = solution = None
    for index, ((u, v, T), case) in enumerate(zip(targets, cell["cases"])):
        witnesses = triples.get(T, ())
        first = witnesses[0] if witnesses else None
        witness = [first[1], first[2], first[0]] if first else None
        multiplicity = sum(1 for item in witnesses if item[0] == first[0]) if first else 0
        tested = first[0] + 1 if first else 8
        assert (case["i"], case["u"], case["v"], case["target"]) == (
            index, u, v, list(T))
        assert case["first_witness"] == witness
        assert case["selected_pair_multiplicity"] == multiplicity
        assert case["third_candidates_tested"] == tested
        assert all(type(value) is int and value >= 0
                   for value in case["scan_cost"].values())
        hits += bool(first)
        independent = False
        if first and first_rank is None:
            coeffs = [0] * 9
            for j in witness:
                coeffs[j] += 1
            coeffs[-1] = -v % R
            rows.append((coeffs, u))
            next_rank, candidate = reference.independent_linear_system(rows, 9)
            independent = next_rank > rank
            dependent += not independent
            rank = next_rank
            if rank == 9:
                first_rank, solution = index + 1, candidate
        assert case["independent_before_rank_stop"] is independent
        assert case["rank_after"] == rank
    assert cell["hits"] == hits and cell["dependent_before_stop"] == dependent
    assert cell["rank"] == rank and cell["first_full_rank_attempt"] == first_rank
    assert cell["solution"] == solution
    assert cell["verified"] == (solution is not None)
    for metric in ("group_add", "mul", "sqr", "inv"):
        prefix = sum(case["scan_cost"].get(metric, 0)
                     for case in cell["cases"][:first_rank or 512])
        tail = sum(case["scan_cost"].get(metric, 0)
                   for case in cell["cases"][first_rank or 512:])
        assert prefix == cell["costs"]["scan_to_rank"].get(metric, 0)
        assert tail == cell["costs"]["audit_tail"].get(metric, 0)
    if solution is not None:
        assert solution[-1] == secret
        assert E.mul(G, solution[-1]) == Q
        assert ref.mul(G, solution[-1]) == Q
        assert all(E.mul(G, d) == P and ref.mul(G, d) == P
                   for d, P in zip(solution[:-1], base))
    return cell


def replay(folder: Path, manifest: Path) -> dict:
    config = checked_config()
    manifest_data = json.loads(manifest.read_text())
    assert manifest_data["schema"] == "ecc2k130-degree7-m3-base-selector-evidence-v1"
    assert sha(folder / "result.json") == manifest_data["result_sha256"]
    result = json.loads((folder / "result.json").read_text())
    assert result["schema"] == "ecc2k130-degree7-m3-base-selector-result-v1"
    assert result["status"] in ("PASS_PANEL", "NEGATIVE_PANEL")
    assert result["config_sha256"] == sha(CONFIG)
    assert result["frozen_sha256"] == sha(FROZEN)
    reference.restore_bare_curve()
    F = pilot.FastGF2m(21, pilot.IRR)
    E0, twist = pilot.Koblitz(F, 0, 1), pilot.Koblitz(F, 1, 1)
    lines, selected, line_trials = pilot.order_seven_lines(F, twist)
    saved = json.loads((ROOT / "research/ecc2k130_factor_base_pilot_20260924/"
                        "results_final.json").read_text())
    assert list(selected) == saved["geometry"]["selected_kernel_x"]
    kernel = lines[selected]
    phi = pilot.BinaryVeluMap.from_generator(E0, twist, kernel["generator"], 7)
    complement = next(lines[key]["generator"] for key in sorted(lines) if key != selected)
    G, generator_trials = pilot.choose_generator(E0)
    D = DualTransport(E0, twist, kernel["generator"], complement, 7, G)
    E1 = pilot.Koblitz(F, phi.codomain.a, phi.codomain.b)
    assert E1.b == kernel["codomain_b"]
    assert D.compose(G) == E0.mul(G, 7)
    reps, by_point, members = reference.independent_orbits(E0, G)
    assert result["geometry"] == {"selected_kernel_x": list(selected),
        "codomain_b": E1.b, "line_trials": line_trials,
        "generator_trials": generator_trials,
        "orbit_representatives": [list(P) for P in reps],
        "orbit_sizes": [len(members[P]) for P in reps]}
    secret = 1 + reference.hash_int(config["labels"]["secret"]) % (R - 1)
    Q = E0.mul(G, secret)
    old = json.loads((ROOT / "research/notes/ecc2k130/m3_four_policy_20260930/"
                      "evidence_run_36722040881/result.json").read_text())
    assert list(Q) != old["challenge"]["Q"]
    assert result["challenge"] == {"G": list(G), "Q": list(Q),
                                   "secret_audit_only": secret}
    G1, Q1 = phi(G), phi(Q)
    assert E1.mul(G1, secret) == Q1
    ref0, ref1 = reference.ReferenceCurve(1), reference.ReferenceCurve(E1.b)
    assert ref0.on_curve(G) and ref0.on_curve(Q)
    assert ref0.mul(G, secret) == Q and ref1.mul(G1, secret) == Q1

    quotas = {rep: config["source_signed_orbit_quota"] for rep in reps[:4]}
    candidates, choices, bases = {}, {}, {}
    reference_checks = 4
    for seed in config["base_seeds"]:
        key = str(seed)
        candidates[key], choices[key], bases[key] = {}, {}, {}
        for role, E in (("source", E0), ("leaf", E1)):
            candidates[key][role], scores = [], []
            for index in range(2):
                base, back, meta, occupancy = replay_candidate(E, seed, role,
                    index, E0, D, by_point, quotas, config)
                score = independent_score(E0, back, set(by_point))
                saved_candidate = result["candidates"][key][role][index]
                assert saved_candidate["base"] == [list(P) for P in base]
                assert saved_candidate["pullback"] == [list(P) for P in back]
                assert saved_candidate["meta"] == meta
                assert saved_candidate["occupancy"] == occupancy
                assert {k: saved_candidate["score"][k] for k in score} == score
                assert saved_candidate["score"]["score_cost"]["group_add"] == 324
                assert saved_candidate["score"]["score_cost"] == result[
                    "phase_costs"][f"score_{role}_{key}_{index}"]
                for P in base:
                    ref = ref0 if role == "source" else ref1
                    assert ref.on_curve(P) and ref.mul(P, R) is None
                    reference_checks += 1
                candidates[key][role].append((base, back))
                scores.append(score)
            choices[key][role] = independent_choice(scores)
        assert result["choices"][key] == choices[key]
        for arm in ARMS:
            source_index = 0 if arm == "control" else choices[key]["source"]
            leaf_index = 0 if arm == "control" else choices[key]["leaf"]
            source_base = candidates[key]["source"][source_index][0]
            native, pullback = candidates[key]["leaf"][leaf_index]
            expected = {"original": source_base,
                "transported": [phi(P) for P in source_base],
                "descendant_native": native, "pullback": pullback}
            assert all(phi(P) == Q for P, Q in zip(pullback, native))
            assert result["bases"][key][arm] == {name: [list(P) for P in points]
                for name, points in expected.items()}
            bases[key][arm] = expected

    targets = {}
    for label, allowed in (("A", set(reps[:5])), ("B", set(reps[5:]))):
        source, meta = target_stream(E0, G, Q, by_point, allowed, label, config)
        leaf = [(u, v, E1.add(E1.mul(G1, u), E1.mul(Q1, v)))
                for u, v, _ in source]
        assert all(phi(T[2]) == U[2] for T, U in zip(source, leaf))
        assert result["target_streams"][label] == {
            "coefficients": [[u, v] for u, v, _ in source],
            "meta": meta,
            "point_orbits": [reps.index(by_point[T]) for _, _, T in source]}
        targets[label], targets[label + "_leaf"] = source, leaf
        for cases, ref, gen, challenge in ((source, ref0, G, Q),
                                           (leaf, ref1, G1, Q1)):
            for u, v, T in cases[:2] + cases[-2:]:
                assert ref.on_curve(T)
                assert ref.add(ref.mul(gen, u), ref.mul(challenge, v)) == T
                reference_checks += 1

    replayed, trajectories = {}, {}
    cell_count = total_cases = 0
    ledger = result["phase_costs"]
    for key in bases:
        replayed[key], trajectories[key] = {}, {}
        for label in config["holdouts"]:
            replayed[key][label], trajectories[key][label] = {}, {}
            for arm in ARMS:
                replayed[key][label][arm], trajectories[key][label][arm] = {}, {}
                for policy in POLICIES:
                    leaf = policy in ("transported", "descendant_native")
                    E, ref, gen, challenge = ((E1, ref1, G1, Q1) if leaf else
                                              (E0, ref0, G, Q))
                    target_cases = targets[label + "_leaf"] if leaf else targets[label]
                    path = folder / "cells" / key / label / arm / f"{policy}.json"
                    cell = replay_cell(path, result["cell_sha256"][key][label][arm][policy],
                        result["cells"][key][label][arm][policy], E, ref, gen,
                        challenge, target_cases, bases[key][arm][policy], secret)
                    cell_count += 1
                    total_cases += len(cell["cases"])
                    replayed[key][label][arm][policy] = {
                        "hits": cell["hits"], "rank": cell["rank"],
                        "first_full_rank_attempt": cell["first_full_rank_attempt"]}
                    trajectories[key][label][arm][policy] = [
                        (case["first_witness"], case["rank_after"])
                        for case in cell["cases"]]
                    limit = cell["first_full_rank_attempt"] or 512
                    shared = [ledger[name] for name in ("field_setup", "generator",
                               "challenge", "source_orbit_universe")]
                    shared += result["target_prefix_costs"][label][:limit]
                    specific = [cell["costs"][name] for name in (
                        "pair_table", "scan_to_rank", "recovery_verify")]
                    role = "source" if policy in ("original", "transported") else "leaf"
                    specific += [ledger[f"candidate_{role}_{key}_{i}"]
                                 for i in ((0,) if arm == "control" else (0, 1))]
                    if arm == "selected":
                        specific += [ledger[f"score_{role}_{key}_{i}"] for i in (0, 1)]
                    if policy == "transported":
                        specific += [ledger[name] for name in ("kernel_search",
                            "forward_setup", "codomain_generator_challenge",
                            f"transport_{key}_{arm}")]
                        specific += result["codomain_target_prefix_costs"][label][:limit]
                    elif policy == "descendant_native":
                        specific += [ledger[name] for name in ("kernel_search",
                            "forward_setup", "dual_setup", "codomain_generator_challenge")]
                        specific += result["codomain_target_prefix_costs"][label][:limit]
                    elif policy == "pullback":
                        specific += [ledger[name] for name in ("kernel_search",
                            "forward_setup", "dual_setup")]
                    assert sum_costs(shared + specific) == result[
                        "cold_cost_to_rank_or_512"][key][label][arm][policy]
                    row_ops = [cell["mod_r_ops"]]
                    if arm == "selected":
                        row_ops += [result["candidates"][key][role][i]["score"][
                            "score_mod_r_ops"] for i in (0, 1)]
                    assert sum_costs(row_ops) == result[
                        "cold_mod_r_ops_to_rank_or_512"][key][label][arm][policy]
                    timing = result["cell_timings"][key][label][arm][policy]
                    assert set(timing) == {"wall_seconds_including_audit_tail",
                                           "cpu_seconds_including_audit_tail"}
                    assert all(type(value) in (int, float) and value >= 0
                               for value in timing.values())
                a, b, c, d = (replayed[key][label][arm][name] for name in POLICIES)
                traj = trajectories[key][label][arm]
                expected_controls = {
                    "forward_hit_parity": a["hits"] == b["hits"],
                    "dual_hit_parity": c["hits"] == d["hits"],
                    "forward_rank_parity": a["first_full_rank_attempt"] == b["first_full_rank_attempt"],
                    "dual_rank_parity": c["first_full_rank_attempt"] == d["first_full_rank_attempt"],
                    "forward_trajectory_parity": traj["original"] == traj["transported"],
                    "dual_trajectory_parity": traj["descendant_native"] == traj["pullback"],
                    "all_verified": all(result["cells"][key][label][arm][name]["verified"]
                                        for name in POLICIES)}
                assert result["controls"][key][label][arm] == expected_controls
    assert cell_count == 64 and total_cases == 64 * 512
    expected_status = "PASS_PANEL" if all(result["cells"][key][label][arm][name]["verified"]
        for key in result["cells"] for label in result["cells"][key]
        for arm in ARMS for name in POLICIES) else "NEGATIVE_PANEL"
    assert result["status"] == expected_status
    return {"schema": "ecc2k130-degree7-m3-base-selector-replay-v1",
        "status": "PASS", "panel_status": expected_status,
        "result_sha256": sha(folder / "result.json"),
        "evidence_manifest_sha256": sha(manifest),
        "config_sha256": sha(CONFIG), "frozen_sha256": sha(FROZEN),
        "cells_replayed": cell_count, "cases_replayed": total_cases,
        "bit_polynomial_point_checks": reference_checks,
        "choices": choices, "cells": replayed}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a replay receipt"
    config = checked_config()
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(
        TimeoutError("frozen verifier wall cap")))
    signal.alarm(config["verifier_wall_cap_seconds"])
    try:
        receipt = replay(args.evidence.resolve(), args.manifest.resolve())
        assert reference.peak_rss_bytes() <= config["process_rss_cap_bytes"]
        save(args.out, receipt)
        print(json.dumps({"status": receipt["status"],
            "panel_status": receipt["panel_status"],
            "cases_replayed": receipt["cases_replayed"]}, sort_keys=True))
    except BaseException as error:
        save(args.out, {"status": "CENSORED" if isinstance(error, (
            TimeoutError, MemoryError)) else "FAIL",
            "error_type": type(error).__name__, "error": str(error),
            "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
            "peak_rss_bytes": reference.peak_rss_bytes()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
