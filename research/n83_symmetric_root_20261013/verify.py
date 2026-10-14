#!/usr/bin/env python3
"""Independently replay the N83 symmetric-S3 target and complete base."""

import argparse
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "research/n83_indexed5_20261013/analyze.py"
ANSWER = ROOT / "research/n83_indexed5_20261013/answer.jsonl"
PRIOR = ROOT / "research/n83_dense_root_20261013/runs/primary"
ARMS = ((2543, True), (2200, True), (1798, True), (1798, False))

spec = importlib.util.spec_from_file_location("n83_independent_group", SOURCE)
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)


def one_json(path):
    return json.loads(path.read_text())


def verify_arm(output, k, symmetric, answer):
    arm = f"k{k}_{'symmetric' if symmetric else 'full'}"
    folder = output / "runs" / "primary" / arm
    if not (folder / "usage.json").exists():
        return {"arm": arm, "status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"arm": arm, "status": usage["status"], "usage": usage}
    summary = one_json(folder / "stdout.jsonl")
    target = one_json(folder / "targets_result.jsonl")
    dump = one_json(folder / "base_dump.jsonl")
    if (summary["n"] != 83 or summary["a"] != 1
            or summary["targets"] != 1 or summary["targets_solved"] != 1
            or summary["targets_failed"] != 0 or summary["orbit_columns"] != k
            or summary["rank"] != k
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_table_indexed5_enabled"] is not True
            or summary["root_table_dense_enabled"] is not False
            or summary["root_table_load75_enabled"] is not False
            or summary["root_state_layout"] != "packed83_11byte_roots"
            or summary["root_state_slot_bytes"] != 28
            or summary["symmetric_s3_states_enabled"] is not symmetric):
        raise ValueError(f"{arm}: source layout, base, rank or target mismatch")
    state_bound = (83 * k * (k - 1) // 2 + 42 * k if symmetric
                   else 83 * k * k)
    if summary["regular_states"] > state_bound:
        raise ValueError(f"{arm}: regular-state count exceeds source-loop bound")
    base = group.replay_base(dump, summary)
    q = group.point_from_json(answer["point"])
    scalar = int(target["recovered_scalar"])
    if (group.point_from_json(target["target"]) != q
            or target["fixture_index"] != 0
            or scalar != int(answer["scalar"])
            or not target["group_verified"]
            or not group.on_curve(q)
            or group.multiply(q, group.R) is not None
            or group.multiply(group.G, scalar) != q):
        raise ValueError(f"{arm}: target scalar or subgroup replay failed")
    indices = target["point_indices"]
    if len(indices) != 4 or any(not isinstance(i, int) or i < 0 or i >= len(base)
                                for i in indices):
        raise ValueError(f"{arm}: invalid four-point indices")
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != q:
        raise ValueError(f"{arm}: four-point relation replay failed")
    phase_ms = sum(float(target[key]) for key in
                   ("target_query_ms", "target_pdp_ms",
                    "target_relation_check_ms", "target_descent_ms",
                    "target_recovery_check_ms"))
    online_ms = float(target["target_ms"])
    if abs(online_ms - phase_ms) > 1.0:
        raise ValueError(f"{arm}: exclusive online phases do not sum")
    parity = None
    if k == 1798 and not symmetric:
        prior = PRIOR / "k1798_dense"
        prior_summary = one_json(prior / "stdout.jsonl")
        prior_target = one_json(prior / "targets_result.jsonl")
        parity = (summary["base_hash"] == prior_summary["base_hash"]
                  and summary["regular_states"] == prior_summary["regular_states"]
                  and indices == prior_target["point_indices"]
                  and target["x_codes"] == prior_target["x_codes"]
                  and target["probes"] == prior_target["probes"]
                  and scalar == prior_target["recovered_scalar"])
        if not parity:
            raise ValueError(f"{arm}: archived full-state first-witness parity failed")
    root_key_count_parity = None
    if k == 1798 and symmetric:
        full_folder = output / "runs" / "primary" / "k1798_full"
        full_usage_path = full_folder / "usage.json"
        if (full_usage_path.exists()
                and one_json(full_usage_path)["status"] == "PASS_PRODUCER"):
            try:
                full_summary = one_json(full_folder / "stdout.jsonl")
            except (OSError, json.JSONDecodeError):
                full_summary = None
            if full_summary is not None:
                root_key_count_parity = (
                    summary["base_hash"] == full_summary["base_hash"]
                    and summary["root_table_entries"] == full_summary["root_table_entries"])
                if not root_key_count_parity:
                    raise ValueError(f"{arm}: exact root key or base count differs from full index")
    return {
        "arm": arm, "status": "PASS_INDEPENDENT_REPLAY",
        "actual_factor_base_points": len(base), "effective_columns": k,
        "base_hash": summary["base_hash"],
        "regular_states": summary["regular_states"],
        "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
        "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
        "target_probes": target["probes"],
        "single_target_online_ms_exploratory": online_ms,
        "setup_ms": summary["timing_ms"]["setup_total"],
        "rank_ms": summary["timing_ms"]["rank_stage"],
        "final_la_ms": summary["timing_ms"]["linear_algebra"],
        "peak_rss_bytes": usage["peak_rss_bytes"],
        "first_witness_parity": parity,
        "root_key_count_parity": root_key_count_parity,
        "symmetric_s3_states_enabled": symmetric,
        "root_table_entries": summary["root_table_entries"],
    }


def verify_rho(output, answer):
    folder = output / "runs" / "rho_primary"
    if not (folder / "usage.json").exists():
        return {"status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"status": usage["status"], "usage": usage}
    receipt = one_json(folder / "receipt.json")
    q = group.point_from_json(answer["point"])
    scalar = int(receipt["recovered_scalar"])
    if (receipt["status"] != "PASS_VERIFIED_PUBLIC_TARGET"
            or group.point_from_json(receipt["public_target_q"]) != q
            or scalar != int(answer["scalar"])
            or not receipt["generic_scalar_replay_verified"]
            or group.multiply(group.G, scalar) != q):
        raise ValueError("strong-rho same-point scalar replay failed")
    return {"status": "PASS_INDEPENDENT_REPLAY",
            "single_target_online_ms_exploratory": receipt["online_ms"],
            "peak_rss_bytes": usage["peak_rss_bytes"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    answer = one_json(ANSWER)
    rows = []
    for k, symmetric in ARMS:
        try:
            rows.append(verify_arm(output, k, symmetric, answer))
        except Exception as error:
            rows.append({"arm": f"k{k}_{'symmetric' if symmetric else 'full'}",
                         "status": "VERIFICATION_FAILED",
                         "error": str(error)})
    try:
        rho = verify_rho(output, answer)
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    result = {"schema": "n83.symmetric-s3-tier1-analysis/v1", "arms": rows,
              "rho": rho, "hardware_tier": "cloud_instance_exploratory",
              "controlled_speedup": None}
    (output / "analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    complete = (all(row["status"] == "PASS_INDEPENDENT_REPLAY" for row in rows)
                and rho["status"] == "PASS_INDEPENDENT_REPLAY")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
