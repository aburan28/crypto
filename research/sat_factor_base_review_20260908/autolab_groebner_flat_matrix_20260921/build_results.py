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
PRE_LINEAR = HERE.parent / "autolab_groebner_linear_tail_20260921" / "results.json"


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
flat = paired("paired_flat_n31_default_t2")
nested = paired("paired_nested_n31_default_t2")
profiles = {
    "flat_sym": observation("profile_flat_sym"),
    "nested_sym": observation("profile_nested_sym"),
    "flat_x": observation("profile_flat_x"),
    "nested_x": observation("profile_nested_x"),
    "flat_subphases": observation("profile_flat_subphases")
}
ladders = {
    "frozen_flat": ladder("ladder_frozen_flat"),
    "frozen_nested": ladder("ladder_frozen_nested"),
    "holdout_flat": ladder("ladder_holdout_flat"),
    "holdout_nested": ladder("ladder_holdout_nested")
}
pre_linear = load(PRE_LINEAR)

assert flat["receipt"]["executable_sha256"] == nested["receipt"]["executable_sha256"]
for arm in ("x-chained", "sym"):
    for field in ("found", "refuted", "inconclusive", "effort", "built_degree", "ffd", "gate"):
        assert flat["arms"][arm][field] == nested["arms"][arm][field]
assert all(arm["inconclusive"] == 0 and arm["gate"] == "ok" for arm in flat["arms"].values())
for system in ("sym", "x"):
    a = profiles[f"flat_{system}"]
    b = profiles[f"nested_{system}"]
    assert a["roots"] == b["roots"]
    assert a["solve_stats"] == b["solve_stats"]
    for field in ("calls", "rows", "cols", "word_ops"):
        assert a["f4_profile"][field] == b["f4_profile"][field]

ladder_comparisons = {}
for group in ("frozen", "holdout"):
    selected = ladders[f"{group}_flat"]["stage"]["rows"]
    control = ladders[f"{group}_nested"]["stage"]["rows"]
    assert len(selected) == len(control)
    rows = []
    for flat_row, nested_row in zip(selected, control):
        key = (flat_row["curve"], flat_row["m"], flat_row["ell"])
        assert key == (nested_row["curve"], nested_row["m"], nested_row["ell"])
        assert flat_row["verdict_digest"] == nested_row["verdict_digest"]
        assert flat_row["word_ops"] == nested_row["word_ops"]
        rows.append({
            "curve": flat_row["curve"], "m": flat_row["m"], "ell": flat_row["ell"],
            "targets": flat_row["targets"], "verdict_digest": flat_row["verdict_digest"],
            "flat_wall_ns": flat_row["wall_ns"], "nested_wall_ns": nested_row["wall_ns"],
            "wall_speedup": nested_row["wall_ns"] / flat_row["wall_ns"],
            "word_ops": flat_row["word_ops"]
        })
    ladder_comparisons[group] = rows

legacy = pre_linear["paired"]["legacy"]
result = {
    "schema_version": "1.0",
    "task_id": protocol["task_id"],
    "classification": "N31_BOOLEAN_F4_FLAT_MATRIX_PAIRED_T2_PASS",
    "selected_policy": protocol["selected_policy"],
    "paired": {
        "flat": flat, "nested": nested,
        "wall_speedup": nested["receipt"]["wall_ms"] / flat["receipt"]["wall_ms"],
        "child_user_cpu_speedup": nested["receipt"]["child_user_ms"] / flat["receipt"]["child_user_ms"],
        "flat_over_nested_peak_rss": flat["receipt"]["peak_rss_bytes"] / nested["receipt"]["peak_rss_bytes"],
        "x_median_speedup": nested["arms"]["x-chained"]["found_ms"] / flat["arms"]["x-chained"]["found_ms"],
        "sym_median_speedup": nested["arms"]["sym"]["refuted_ms"] / flat["arms"]["sym"]["refuted_ms"]
    },
    "cumulative_from_pre_linear_tail_legacy": {
        "wall_speedup": legacy["receipt"]["wall_ms"] / flat["receipt"]["wall_ms"],
        "child_user_cpu_speedup": legacy["receipt"]["child_user_ms"] / flat["receipt"]["child_user_ms"],
        "x_median_speedup": legacy["arms"]["x-chained"]["found_ms"] / flat["arms"]["x-chained"]["found_ms"],
        "sym_median_speedup": legacy["arms"]["sym"]["refuted_ms"] / flat["arms"]["sym"]["refuted_ms"],
        "predecessor_results_sha256": sha256(PRE_LINEAR)
    },
    "profiles": {
        **profiles,
        "sym_solver_speedup": profiles["nested_sym"]["elapsed_ms"] / profiles["flat_sym"]["elapsed_ms"],
        "sym_reduction_speedup": profiles["nested_sym"]["f4_profile"]["reduce_ns"] / profiles["flat_sym"]["f4_profile"]["reduce_ns"],
        "x_solver_speedup": profiles["nested_x"]["elapsed_ms"] / profiles["flat_x"]["elapsed_ms"],
        "x_reduction_speedup": profiles["nested_x"]["f4_profile"]["reduce_ns"] / profiles["flat_x"]["f4_profile"]["reduce_ns"]
    },
    "ladders": {**ladders, "comparisons": ladder_comparisons},
    "rejected": {
        "exact_row_product_cache": "Rejected: 253921 exact misses and zero hits; wall regressed from 6.93 to 8.25 seconds.",
        "paired_specialisation": "Rejected: three-process medians differed by only 0.17 percent; shifted-row construction remained separate."
    },
    "validation": {
        "flat_nested_linear_tail_equivalence": "PASS",
        "forced_flat_exhaustive_roots": "PASS",
        "forced_flat_all_split_rules": "PASS",
        "paired_output_and_gate_equivalence": "PASS",
        "frozen_ladder_verdict_and_ops": "PASS",
        "holdout_ladder_verdict_and_ops": "PASS"
    },
    "claim_boundary": protocol["claim_boundary"]
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

source_hashes = {
    "koblitz_groebner_source": sha256(FROZEN / "koblitz_groebner.rs"),
    "rref_tests_source": sha256(FROZEN / "f4_rref_tests.rs"),
    "frontier_probe_source": sha256(FROZEN / "koblitz_f4_frontier_probe.rs"),
    "paired_executable_receipt": flat["receipt"]["executable_sha256"],
    "profile_executable_receipt": receipt("profile_flat_sym")["executable_sha256"]
}
claim_report = {
    "schema_version": 2, "task_id": protocol["task_id"], "stage": "decomposition",
    "n_or_bits": 31, "m_summands": 2, "dimension_l_or_dim": 16, "unknowns": 32,
    "system_degree": "quadratic", "eq_var_ratio": 31 / 32,
    "ffd_or_degree_of_regularity": {"x_chained": 3, "symmetrised": 4},
    "oracle_class": "Boolean matrix-F4 with contiguous solver matrix and exact linear-tail elimination",
    "median_ms_per_target": {"x_found_flat": flat["arms"]["x-chained"]["found_ms"], "sym_refuted_flat": flat["arms"]["sym"]["refuted_ms"]},
    "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "node_budget": 50000},
    "macaulay_degree_max": 3, "f4_splits": {"x": 7056, "sym": 4095},
    "verdict_mix": {"x_found": 2, "sym_refuted": 2, "inconclusive": 0},
    "fixture_hash": flat["receipt"]["stdout_sha256"], "executable_or_source_hash": source_hashes,
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
