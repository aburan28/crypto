#!/usr/bin/env python3
"""Independent GF(2^73) replay for the extra paired runs (R2, R3).

Same checks as verify_single_target.py, parameterized by run tag so each
additional frozen run dir gets its own replay receipt.  Usage:

    python3 verify_single_target_r2_r3.py R2
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_single_target as base  # noqa: E402


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ("R2", "R3"):
        raise SystemExit("usage: verify_single_target_r2_r3.py R2|R3")
    tag = sys.argv[1]
    root = Path(__file__).resolve().parent
    protocol = base.load_json(root / "protocol.json")
    run_dir = root / "runs" / (protocol["run_id"].removesuffix("R1") + tag)
    if not run_dir.exists():
        raise SystemExit(f"missing run dir {run_dir}")

    # Reuse the frozen-module checks with this run's directory.
    fixture = base.load_json(root / "frozen" / "fixture.json")
    target_rows = [
        line
        for line in (root / "frozen" / "target_points.jsonl")
        .read_text()
        .splitlines()
        if line.strip()
    ]
    if len(target_rows) != 1:
        raise ValueError("frozen target input must contain exactly one point")
    ic = base.load_one_jsonl(run_dir / "ic.jsonl")
    ic_summary = base.load_one_jsonl(run_dir / "ic_summary.jsonl")
    rho = base.load_one_jsonl(run_dir / "rho.jsonl")
    ic_execution = base.load_json(run_dir / "ic_execution.json")
    rho_execution = base.load_json(run_dir / "rho_execution.json")
    rho_peak_rss_bytes = int(rho_execution["peak_rss_bytes"])
    ic_peak_rss_bytes = (
        int(ic_summary["peak_rss_bytes"])
        if isinstance(ic_summary.get("peak_rss_bytes"), int)
        else int(ic_execution["peak_rss_bytes"])
    )

    n = int(fixture["n"])
    low_terms = [int(term) for term in fixture["field_modulus_low_terms"]]
    modulus = (1 << n) | 1
    for term in low_terms:
        modulus |= 1 << term
    generator = tuple(int(value) for value in fixture["generator"])
    target = tuple(int(value) for value in fixture["public_target"])
    scalar = int(fixture["fixture_scalar_validation_only"])
    order = int(fixture["subgroup_order"])
    cofactor = int(fixture["cofactor"])
    group_order = int(fixture["curve_order"])
    curve_trace = int(fixture["curve_trace_over_f2n"])

    checks = base.curve_replay(n, modulus, int(fixture["a"]), 1, generator, target, scalar, order)
    checks["fixture_order_factorization_matches"] = order * cofactor == group_order
    checks["fixture_trace_matches_order"] = (1 << n) + 1 - group_order == curve_trace
    checks["ic_uses_frozen_Q"] = tuple(int(v) for v in ic["published_q"]) == target
    checks["rho_uses_frozen_Q"] = tuple(int(v) for v in rho["published_q"]) == target
    checks["rho_uses_frozen_G"] = tuple(int(v) for v in rho["generator"]) == generator
    checks["ic_recovers_fixture_scalar"] = int(ic["recovered_scalar"]) == scalar
    checks["rho_recovers_fixture_scalar"] = int(rho["recovered_fixture_scalar"]) == scalar
    checks["ic_in_process_verification_passed"] = ic.get("group_verified") is True
    checks["rho_in_process_verification_passed"] = rho.get("verified") is True
    checks["one_ic_online_target"] = ic_summary["targets"] == 1 and ic_summary["targets_solved"] == 1
    checks["one_rho_online_target"] = rho["fixture_index"] == 0
    checks["known_answer_not_sent_to_ic"] = ic.get("published_fixture_scalar") is None
    checks["known_answer_not_sent_to_rho"] = (
        rho.get("published_fixture_scalar") is None
        and not rho_execution["known_answer_environment_variable_present"]
    )
    checks["ic_process_succeeded"] = ic_execution["return_code"] == 0
    checks["rho_process_succeeded"] = rho_execution["return_code"] == 0
    checks["no_explicit_pair_table"] = ic_summary["pair_table_entries"] == 0
    checks["no_edge_selectors"] = ic_summary["edge_selectors"] == 0

    ic_phases = {
        "target_query": float(ic["target_query_ms"]),
        "target_PDP": float(ic["target_pdp_ms"]),
        "target_relation_check": float(ic["target_relation_check_ms"]),
        "target_descent": float(ic["target_descent_ms"]),
        "target_recovery_check": float(ic["target_recovery_check_ms"]),
    }
    ic_online = float(ic["target_ms"])
    ic_phase_sum = sum(ic_phases.values())
    rho_phases = {
        "target_specific_jump_setup": float(rho["setup_ms"]),
        "target_walk": float(rho["walk_ms"]),
        "scalar_recovery_check": float(rho["validation_ms"]),
    }
    rho_online = sum(rho_phases.values())
    checks["ic_online_phase_sum_matches"] = abs(ic_phase_sum - ic_online) <= max(
        0.02, ic_online * 1e-8
    )
    checks["rho_used_supplied_public_point"] = rho["target_input_kind"] == "public_point"
    checks["rho_fixture_generation_excluded"] = float(rho["target_generation_ms_excluded"]) == 0.0

    online_speedup = rho_online / ic_online
    passed = all(checks.values())
    receipt = {
        "kind": "independent_n73_single_target_replay",
        "status": "PASS" if passed else "FAIL",
        "run_id": run_dir.name,
        "run_tag": tag,
        "rho_walk_seed": base.load_json(run_dir / "run.json")["rho_walk_seed"],
        "ic_rank_seed": base.load_json(run_dir / "run.json")["ic_rank_seed"],
        "scope": "Standalone standard-library GF(2^73) curve replay, same-point pairing, one-target accounting, and phase sums; does not independently remeasure producer runtimes.",
        "candidate_id": protocol["candidate_id"],
        "workload_id": protocol["workload_id"],
        "checks": checks,
        "timing_ms": {
            "ic_online": ic_online,
            "rho_online": rho_online,
            "online_speedup": online_speedup,
            "ic_phases": ic_phases,
            "rho_phases": rho_phases,
            "ic_peak_rss_bytes": ic_peak_rss_bytes,
            "rho_peak_rss_bytes": rho_peak_rss_bytes,
        },
        "inputs_sha256": {
            "ic": base.sha256(run_dir / "ic.jsonl"),
            "ic_summary": base.sha256(run_dir / "ic_summary.jsonl"),
            "rho": base.sha256(run_dir / "rho.jsonl"),
            "rho_execution": base.sha256(run_dir / "rho_execution.json"),
            "ic_execution": base.sha256(run_dir / "ic_execution.json"),
        },
    }
    out = run_dir / "independent_replay.json"
    out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"run": tag, "status": receipt["status"], "speedup": online_speedup}, sort_keys=True))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
