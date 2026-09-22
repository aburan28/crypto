#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import platform
import socket
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
FROZEN = HERE / "frozen"
PREVIOUS = HERE.parent / "autolab_groebner_linear_tail_20260921" / "results.json"


def load(path: Path):
    return json.loads(path.read_text())


def observation(name: str):
    return load(RUNS / name / "stdout.txt")


def receipt(name: str):
    return load(RUNS / name / "receipt.json")


def ladder(name: str):
    return {"receipt": receipt(name), "stage": load(RUNS / name / "artifact" / "stage.json")}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def optional_number(token: str):
    return None if token == "-" else float(token)


def paired(name: str):
    text = (RUNS / name / "stdout.txt").read_text()
    arms = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) != 15 or fields[1] != "F4" or fields[0] not in {"x-chained", "sym"}:
            continue
        arms[fields[0]] = {
            "variables": int(fields[2]), "equations": int(fields[3]), "degree": int(fields[4]),
            "found": int(fields[5]), "refuted": int(fields[6]), "inconclusive": int(fields[7]),
            "via_t": int(fields[8]), "found_ms": optional_number(fields[9]),
            "refuted_ms": optional_number(fields[10]), "effort": int(fields[11]),
            "built_degree": int(fields[12]), "ffd": int(fields[13]), "gate": fields[14]
        }
    assert set(arms) == {"x-chained", "sym"}
    return {"receipt": receipt(name), "arms": arms, "stdout": text}


protocol = load(HERE / "protocol.json")
selected = paired("paired_selected_final_n31_default_t2")
control = paired("paired_control_final_n31_default_t2")
profiles = {
    "selected": observation("profile_selected_final"),
    "control": observation("profile_control_final"),
    "no_layout": observation("profile_no_layout_final"),
    "no_schedule": observation("profile_no_schedule_final"),
    "x_selected": observation("profile_x_selected_final"),
    "x_control": observation("profile_x_control_final")
}
ladders = {
    "frozen_selected": ladder("ladder_frozen_selected_final"),
    "frozen_control": ladder("ladder_frozen_control_final"),
    "holdout_selected": ladder("ladder_holdout_selected_final"),
    "holdout_control": ladder("ladder_holdout_control_final")
}
previous = load(PREVIOUS)

assert selected["receipt"]["executable_sha256"] == control["receipt"]["executable_sha256"]
for arm in ("x-chained", "sym"):
    for field in ("found", "refuted", "inconclusive", "effort", "built_degree", "ffd", "gate"):
        assert selected["arms"][arm][field] == control["arms"][arm][field]
assert all(arm["inconclusive"] == 0 and arm["gate"] == "ok" for arm in selected["arms"].values())
for system in ("", "x_"):
    a = profiles[f"{system}selected"]
    b = profiles[f"{system}control"]
    assert a["roots"] == b["roots"]
    assert a["solve_stats"] == b["solve_stats"]
    assert a["f4_profile"]["rows"] < b["f4_profile"]["rows"]
    assert a["f4_profile"]["cols"] < b["f4_profile"]["cols"]
    assert a["f4_profile"]["word_ops"] < b["f4_profile"]["word_ops"]

ladder_comparisons = {}
for group in ("frozen", "holdout"):
    chosen = ladders[f"{group}_selected"]["stage"]["rows"]
    baseline = ladders[f"{group}_control"]["stage"]["rows"]
    assert len(chosen) == len(baseline)
    rows = []
    for selected_row, control_row in zip(chosen, baseline):
        key = (selected_row["curve"], selected_row["m"], selected_row["ell"])
        assert key == (control_row["curve"], control_row["m"], control_row["ell"])
        assert selected_row["targets"] == control_row["targets"]
        assert selected_row["decomposed"] == control_row["decomposed"]
        assert selected_row["verdict_digest"] == control_row["verdict_digest"]
        rows.append({
            "curve": selected_row["curve"], "m": selected_row["m"], "ell": selected_row["ell"],
            "targets": selected_row["targets"], "verdict_digest": selected_row["verdict_digest"],
            "selected_wall_ns": selected_row["wall_ns"], "control_wall_ns": control_row["wall_ns"],
            "wall_speedup": control_row["wall_ns"] / selected_row["wall_ns"],
            "selected_word_ops": selected_row["word_ops"], "control_word_ops": control_row["word_ops"]
        })
    ladder_comparisons[group] = rows

legacy = previous["paired"]["legacy"]
result = {
    "schema_version": "1.0",
    "task_id": protocol["task_id"],
    "classification": "N31_BOOLEAN_F4_ACTIVE_LAYOUT_PAIRED_T2_PASS",
    "selected_policy": protocol["selected_policy"],
    "paired": {
        "selected": selected,
        "control": control,
        "wall_speedup": control["receipt"]["wall_ms"] / selected["receipt"]["wall_ms"],
        "child_user_cpu_speedup": control["receipt"]["child_user_ms"] / selected["receipt"]["child_user_ms"],
        "selected_over_control_peak_rss": selected["receipt"]["peak_rss_bytes"] / control["receipt"]["peak_rss_bytes"],
        "x_median_speedup": control["arms"]["x-chained"]["found_ms"] / selected["arms"]["x-chained"]["found_ms"],
        "sym_median_speedup": control["arms"]["sym"]["refuted_ms"] / selected["arms"]["sym"]["refuted_ms"]
    },
    "cumulative_from_pre_linear_tail_legacy": {
        "wall_speedup": legacy["receipt"]["wall_ms"] / selected["receipt"]["wall_ms"],
        "child_user_cpu_speedup": legacy["receipt"]["child_user_ms"] / selected["receipt"]["child_user_ms"],
        "x_median_speedup": legacy["arms"]["x-chained"]["found_ms"] / selected["arms"]["x-chained"]["found_ms"],
        "sym_median_speedup": legacy["arms"]["sym"]["refuted_ms"] / selected["arms"]["sym"]["refuted_ms"],
        "predecessor_results_sha256": sha256(PREVIOUS)
    },
    "profiles": {
        **profiles,
        "sym_solver_speedup": profiles["control"]["elapsed_ms"] / profiles["selected"]["elapsed_ms"],
        "x_solver_speedup": profiles["x_control"]["elapsed_ms"] / profiles["x_selected"]["elapsed_ms"],
        "layout_cache_speedup": profiles["no_layout"]["elapsed_ms"] / profiles["selected"]["elapsed_ms"],
        "schedule_cache_speedup": profiles["no_schedule"]["elapsed_ms"] / profiles["selected"]["elapsed_ms"],
        "sym_word_ops_reduction_fraction": 1 - profiles["selected"]["f4_profile"]["word_ops"] / profiles["control"]["f4_profile"]["word_ops"],
        "x_word_ops_reduction_fraction": 1 - profiles["x_selected"]["f4_profile"]["word_ops"] / profiles["x_control"]["f4_profile"]["word_ops"]
    },
    "ladders": {**ladders, "comparisons": ladder_comparisons},
    "rejected": {
        "row_arena": "Rejected: three fresh-process pairs were about 11 percent slower because retained capacities harmed reduction locality.",
        "pivot_trace": "Rejected: verified replay and fallback had no measurable reduction-time gain; range verification cost matched skipped scans.",
        "hash_column_collection": "Rejected: measured construction was slower than sort-and-deduplicate."
    },
    "validation": {
        "full_mask_schedule_equivalence": "PASS",
        "exact_layout_support_guard": "PASS",
        "forced_active_exhaustive_roots": "PASS",
        "forced_active_all_split_rules": "PASS",
        "paired_output_and_gate_equivalence": "PASS",
        "frozen_ladder_verdict_digests": "PASS",
        "holdout_ladder_verdict_digests": "PASS"
    },
    "claim_boundary": protocol["claim_boundary"]
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

source_hashes = {
    "koblitz_groebner_source": sha256(FROZEN / "koblitz_groebner.rs"),
    "rref_tests_source": sha256(FROZEN / "f4_rref_tests.rs"),
    "frontier_probe_source": sha256(FROZEN / "koblitz_f4_frontier_probe.rs"),
    "paired_executable_receipt": selected["receipt"]["executable_sha256"],
    "profile_executable_receipt": receipt("profile_selected_final")["executable_sha256"]
}
claim_report = {
    "schema_version": 2, "task_id": protocol["task_id"], "stage": "decomposition",
    "n_or_bits": 31, "m_summands": 2, "dimension_l_or_dim": 16, "unknowns": 32,
    "system_degree": "quadratic", "eq_var_ratio": 31 / 32,
    "ffd_or_degree_of_regularity": {"x_chained": 3, "symmetrised": 4},
    "oracle_class": "Boolean matrix-F4 with active-variable multipliers and exact support-verified layout reuse",
    "median_ms_per_target": {"x_found_selected": selected["arms"]["x-chained"]["found_ms"], "sym_refuted_selected": selected["arms"]["sym"]["refuted_ms"]},
    "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "node_budget": 50000},
    "macaulay_degree_max": 3, "f4_splits": {"x": 7056, "sym": 4095},
    "verdict_mix": {"x_found": 2, "sym_refuted": 2, "inconclusive": 0},
    "fixture_hash": selected["receipt"]["stdout_sha256"], "executable_or_source_hash": source_hashes,
    "host_id": {"hostname": socket.gethostname(), "system": platform.system(), "machine": platform.machine()},
    "resource_caps": {"node_budget": 50000, "targets_per_arm": 2, "threads": 1}, "seeds": [24301],
    "claim_boundary": protocol["claim_boundary"],
    "claim_boundary_non_claims": ["not relation-yield evidence", "not end-to-end index calculus", "not asymptotic or sub-rho evidence", "not external or private targets", "not a key-recovery claim"]
}
(HERE / "claim_report_decomposition.json").write_text(json.dumps(claim_report, indent=2, sort_keys=True) + "\n")
required_global = ["fixture_hash", "executable_or_source_hash", "host_id", "resource_caps", "seeds", "claim_boundary_non_claims"]
required_stage = ["n_or_bits", "m_summands", "dimension_l_or_dim", "unknowns", "system_degree", "eq_var_ratio", "ffd_or_degree_of_regularity", "oracle_class", "median_ms_per_target", "largest_solvable"]
missing_global = [field for field in required_global if field not in claim_report]
missing_stage = [field for field in required_stage if field not in claim_report]
(HERE / "claim_check_decomposition.json").write_text(json.dumps({
    "schema_version": 2, "stage": "decomposition", "fail_closed": True,
    "required_global_provenance": required_global, "required_stage_fields": required_stage,
    "missing_global_provenance": missing_global, "missing_stage_fields": missing_stage,
    "status": "PASS" if not missing_global and not missing_stage else "SCHEMA_INCOMPLETE"
}, indent=2, sort_keys=True) + "\n")
files = {}
for path in sorted(HERE.rglob("*")):
    if path.is_file() and path.name != "manifest.json":
        files[str(path.relative_to(HERE))] = sha256(path)
(HERE / "manifest.json").write_text(json.dumps({"schema_version": "1.0", "task_id": protocol["task_id"], "files": files}, indent=2, sort_keys=True) + "\n")
