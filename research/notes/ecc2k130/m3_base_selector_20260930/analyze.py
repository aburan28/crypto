#!/usr/bin/env python3
"""Derive the frozen selector decision from pinned raw and replay evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import verify

HERE = Path(__file__).resolve().parent
POLICIES = {"source": "original", "leaf": "descendant_native"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(evidence: Path) -> dict:
    verify.checked_config()
    manifest_path = evidence / "EVIDENCE.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema"] == "ecc2k130-degree7-m3-base-selector-evidence-v1"
    assert sha(evidence / "result.json") == manifest["result_sha256"]
    result = json.loads((evidence / "result.json").read_text())
    replay = json.loads((evidence / "replay.json").read_text())
    assert result["status"] == "PASS_PANEL"
    assert replay["status"] == "PASS" and replay["panel_status"] == "PASS_PANEL"
    assert replay["result_sha256"] == manifest["result_sha256"]
    assert replay["evidence_manifest_sha256"] == sha(manifest_path)
    assert replay["cells_replayed"] == 64 and replay["cases_replayed"] == 32768
    assert result["source_head"] == manifest["source_head"]
    assert manifest["raw_cell_files"] == 64 and manifest["cases_per_cell"] == 512
    assert result["method_crossover"] is None

    # This per-holdout breakdown is exploratory: it was not the frozen score.
    verify.reference.restore_bare_curve()
    field = verify.pilot.FastGF2m(21, verify.pilot.IRR)
    curve = verify.pilot.Koblitz(field, 0, 1)
    generator = tuple(result["challenge"]["G"])
    reps, by_point, _ = verify.reference.independent_orbits(curve, generator)
    allowed = {"A": set(reps[:5]), "B": set(reps[5:])}
    support_by_holdout = {}
    for seed in sorted(result["candidates"]):
        support_by_holdout[seed] = {}
        for role in POLICIES:
            support_by_holdout[seed][role] = []
            for candidate in result["candidates"][seed][role]:
                base = [tuple(P) for P in candidate["pullback"]]
                triples, _ = verify.reference.brute_triples(curve, base)
                triples.pop(None, None)
                assert len(triples) == candidate["score"]["distinct_support"]
                counts = {label: sum(by_point[T] in orbits for T in triples)
                          for label, orbits in allowed.items()}
                assert sum(counts.values()) == len(triples)
                support_by_holdout[seed][role].append(counts)

    pairs = []
    covariance_pairs = 0
    for seed in sorted(result["choices"]):
        for holdout in ("A", "B"):
            for arm in ("control", "selected"):
                cells = result["cells"][seed][holdout][arm]
                assert cells["original"]["hits"] == cells["transported"]["hits"]
                assert cells["descendant_native"]["hits"] == cells["pullback"]["hits"]
                assert cells["original"]["first_full_rank_attempt"] == (
                    cells["transported"]["first_full_rank_attempt"])
                assert cells["descendant_native"]["first_full_rank_attempt"] == (
                    cells["pullback"]["first_full_rank_attempt"])
                assert all(result["controls"][seed][holdout][arm].values())
                covariance_pairs += 2
            for role, policy in POLICIES.items():
                control = result["cells"][seed][holdout]["control"][policy]
                selected = result["cells"][seed][holdout]["selected"][policy]
                cold = result["cold_cost_to_rank_or_512"][seed][holdout]
                mod = result["cold_mod_r_ops_to_rank_or_512"][seed][holdout]
                first = cold["control"][policy]
                second = cold["selected"][policy]
                assert control["verified"] and selected["verified"]
                assert control["rank"] == selected["rank"] == 9
                pair = {"seed": int(seed), "holdout": holdout,
                    "role": role, "policy": policy,
                    "selected_candidate": result["choices"][seed][role],
                    "control_hits_out_of_512": control["hits"],
                    "selected_hits_out_of_512": selected["hits"],
                    "control_first_rank": control["first_full_rank_attempt"],
                    "selected_first_rank": selected["first_full_rank_attempt"],
                    "control_cold_field_mul": first["mul"],
                    "selected_cold_field_mul": second["mul"],
                    "delta_cold_field_mul": second["mul"] - first["mul"],
                    "selected_div_control_mul": second["mul"] / first["mul"],
                    "control_cold_field_sqr": first["sqr"],
                    "selected_cold_field_sqr": second["sqr"],
                    "control_cold_inversion_calls": first["inv"],
                    "selected_cold_inversion_calls": second["inv"],
                    "control_cold_group_add": first["group_add"],
                    "selected_cold_group_add": second["group_add"],
                    "control_mod_r_ops": sum(mod["control"][policy].values()),
                    "selected_mod_r_ops": sum(mod["selected"][policy].values()),
                    "strictly_lower_mul": second["mul"] < first["mul"],
                    "no_hit_regression": selected["hits"] >= control["hits"],
                    "no_rank_delay": selected["first_full_rank_attempt"] <= (
                        control["first_full_rank_attempt"]),
                    "verified_rank9_Q": True}
                pairs.append(pair)
    assert covariance_pairs == 32 and len(pairs) == 16
    counts = {role: {"pairs": sum(p["role"] == role for p in pairs),
        "strictly_lower_mul": sum(p["role"] == role and p["strictly_lower_mul"]
                                  for p in pairs),
        "hit_regressions": sum(p["role"] == role and not p["no_hit_regression"]
                               for p in pairs),
        "rank_delays": sum(p["role"] == role and not p["no_rank_delay"]
                           for p in pairs)} for role in POLICIES}
    assert counts["source"]["pairs"] == counts["leaf"]["pairs"] == 8
    advantage = all(p["strictly_lower_mul"] and p["no_hit_regression"] and
                    p["no_rank_delay"] and p["verified_rank9_Q"] for p in pairs)
    decision = "CHARGED_SELECTOR_ADVANTAGE" if advantage else "NO_CHARGED_SELECTOR_ADVANTAGE"
    preflight = HERE / manifest["preflight_without_fixtures"]
    assert json.loads(preflight.read_text())["status"] == "PREFLIGHT_NO_FIXTURES"
    return {"schema": "ecc2k130-degree7-m3-base-selector-analysis-v1",
        "decision": decision, "panel_status": result["status"],
        "replay_status": replay["status"],
        "source_lock_pr": manifest["source_lock_pr"],
        "producer_source_head": result["source_head"],
        "config_sha256": result["config_sha256"],
        "source_lock_sha256": result["frozen_sha256"],
        "result_sha256": manifest["result_sha256"],
        "manifest_sha256": sha(manifest_path),
        "replay_sha256": sha(evidence / "replay.json"),
        "preflight_failure_sha256": sha(preflight),
        "run_host": {"platform": result["platform"],
            "python": result["python"], "hostname": result["host"],
            "wall_seconds": result["wall_seconds"],
            "cpu_seconds": result["cpu_seconds"],
            "peak_rss_bytes": result["peak_rss_bytes"]},
        "selector_choices": result["choices"],
        "candidate_support_by_holdout_posthoc": support_by_holdout,
        "primary_pairs": pairs, "gate_counts": counts,
        "mapped_covariance_pairs": covariance_pairs,
        "ECC2K_130_PDP_yield": None,
        "full_ECDLP_cost": None, "rho_crossover": None}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite analysis output"
    result = analyze(args.evidence.resolve())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"decision": result["decision"],
        "pairs": len(result["primary_pairs"])}, sort_keys=True))


if __name__ == "__main__":
    main()
