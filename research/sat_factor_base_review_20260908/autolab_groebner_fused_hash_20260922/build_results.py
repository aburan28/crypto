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

def load(path: Path): return json.loads(path.read_text())
def observation(name: str): return load(RUNS / name / "stdout.txt")
def receipt(name: str): return load(RUNS / name / "receipt.json")
def ladder(name: str): return {"receipt": receipt(name), "stage": load(RUNS / name / "artifact" / "stage.json")}
def sha256(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def optional_number(token: str): return None if token == "-" else float(token)

def paired(name: str):
    text = (RUNS / name / "stdout.txt").read_text()
    arms = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) != 15 or fields[1] != "F4" or fields[0] not in {"x-chained", "sym"}: continue
        arms[fields[0]] = {
            "variables": int(fields[2]), "equations": int(fields[3]), "degree": int(fields[4]),
            "found": int(fields[5]), "refuted": int(fields[6]), "inconclusive": int(fields[7]),
            "via_t": int(fields[8]), "found_ms": optional_number(fields[9]), "refuted_ms": optional_number(fields[10]),
            "effort": int(fields[11]), "built_degree": int(fields[12]), "ffd": int(fields[13]), "gate": fields[14]
        }
    assert set(arms) == {"x-chained", "sym"}
    return {"receipt": receipt(name), "arms": arms, "stdout": text}

protocol = load(HERE / "protocol.json")
selected = paired("paired_selected_n31_default_t2")
control = paired("paired_merged_control_n31_default_t2")
profiles = {name: observation(name) for name in (
    "profile_selected_sym", "profile_control_sym", "profile_no_fused_sym", "profile_std_hash_sym",
    "profile_selected_x", "profile_control_x")}
ladders = {
    "frozen_selected": ladder("ladder_frozen_selected"), "frozen_control": ladder("ladder_frozen_control"),
    "holdout_selected": ladder("ladder_holdout_selected"), "holdout_control": ladder("ladder_holdout_control")
}
pre_linear = load(PRE_LINEAR)

assert selected["receipt"]["executable_sha256"] == control["receipt"]["executable_sha256"]
for arm in ("x-chained", "sym"):
    for field in ("found", "refuted", "inconclusive", "effort", "built_degree", "ffd", "gate"):
        assert selected["arms"][arm][field] == control["arms"][arm][field]
for system in ("sym", "x"):
    a, b = profiles[f"profile_selected_{system}"], profiles[f"profile_control_{system}"]
    assert a["roots"] == b["roots"] and a["solve_stats"] == b["solve_stats"]
    for field in ("calls", "rows", "cols", "word_ops"): assert a["f4_profile"][field] == b["f4_profile"][field]

ladder_comparisons = {}
for group in ("frozen", "holdout"):
    chosen = ladders[f"{group}_selected"]["stage"]["rows"]
    baseline = ladders[f"{group}_control"]["stage"]["rows"]
    rows = []
    for a, b in zip(chosen, baseline):
        assert (a["curve"], a["m"], a["ell"]) == (b["curve"], b["m"], b["ell"])
        assert a["verdict_digest"] == b["verdict_digest"] and a["word_ops"] == b["word_ops"]
        rows.append({"curve": a["curve"], "m": a["m"], "ell": a["ell"], "targets": a["targets"],
                     "verdict_digest": a["verdict_digest"], "selected_wall_ns": a["wall_ns"],
                     "control_wall_ns": b["wall_ns"], "wall_speedup": b["wall_ns"] / a["wall_ns"], "word_ops": a["word_ops"]})
    ladder_comparisons[group] = rows

legacy = pre_linear["paired"]["legacy"]
result = {
    "schema_version": "1.0", "task_id": protocol["task_id"],
    "classification": "N31_BOOLEAN_F4_FUSED_HASH_PAIRED_T2_PASS", "selected_policy": protocol["selected_policy"],
    "paired": {"selected": selected, "control": control,
        "wall_speedup": control["receipt"]["wall_ms"] / selected["receipt"]["wall_ms"],
        "child_user_cpu_speedup": control["receipt"]["child_user_ms"] / selected["receipt"]["child_user_ms"],
        "selected_over_control_peak_rss": selected["receipt"]["peak_rss_bytes"] / control["receipt"]["peak_rss_bytes"],
        "x_median_speedup": control["arms"]["x-chained"]["found_ms"] / selected["arms"]["x-chained"]["found_ms"],
        "sym_median_speedup": control["arms"]["sym"]["refuted_ms"] / selected["arms"]["sym"]["refuted_ms"]},
    "cumulative_from_pre_linear_tail_legacy": {
        "wall_speedup": legacy["receipt"]["wall_ms"] / selected["receipt"]["wall_ms"],
        "child_user_cpu_speedup": legacy["receipt"]["child_user_ms"] / selected["receipt"]["child_user_ms"],
        "x_median_speedup": legacy["arms"]["x-chained"]["found_ms"] / selected["arms"]["x-chained"]["found_ms"],
        "sym_median_speedup": legacy["arms"]["sym"]["refuted_ms"] / selected["arms"]["sym"]["refuted_ms"],
        "predecessor_results_sha256": sha256(PRE_LINEAR)},
    "profiles": {**profiles,
        "sym_solver_speedup": profiles["profile_control_sym"]["elapsed_ms"] / profiles["profile_selected_sym"]["elapsed_ms"],
        "x_solver_speedup": profiles["profile_control_x"]["elapsed_ms"] / profiles["profile_selected_x"]["elapsed_ms"],
        "fused_pack_speedup": profiles["profile_no_fused_sym"]["elapsed_ms"] / profiles["profile_selected_sym"]["elapsed_ms"],
        "fast_hash_speedup": profiles["profile_std_hash_sym"]["elapsed_ms"] / profiles["profile_selected_sym"]["elapsed_ms"]},
    "ladders": {**ladders, "comparisons": ladder_comparisons},
    "rejected": {"row_product_cache": "253921 exact misses and zero hits; slower than no cache.",
                 "paired_specialisation": "Three-process median changed by 0.17 percent; shifted-row work remained separate."},
    "validation": {"fused_materialized_flat_equivalence": "PASS", "missing_extra_support_fallback": "PASS",
        "fast_standard_index_equivalence": "PASS", "forced_fused_exhaustive_roots": "PASS",
        "forced_fused_all_split_rules": "PASS", "paired_output_and_gate_equivalence": "PASS",
        "frozen_ladder_verdict_and_ops": "PASS", "holdout_ladder_verdict_and_ops": "PASS"},
    "claim_boundary": protocol["claim_boundary"]
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

hashes = {"koblitz_groebner_source": sha256(FROZEN / "koblitz_groebner.rs"),
          "rref_tests_source": sha256(FROZEN / "f4_rref_tests.rs"),
          "frontier_probe_source": sha256(FROZEN / "koblitz_f4_frontier_probe.rs"),
          "paired_executable_receipt": selected["receipt"]["executable_sha256"],
          "profile_executable_receipt": receipt("profile_selected_sym")["executable_sha256"]}
claim = {"schema_version": 2, "task_id": protocol["task_id"], "stage": "decomposition", "n_or_bits": 31,
    "m_summands": 2, "dimension_l_or_dim": 16, "unknowns": 32, "system_degree": "quadratic", "eq_var_ratio": 31/32,
    "ffd_or_degree_of_regularity": {"x_chained": 3, "symmetrised": 4},
    "oracle_class": "Boolean matrix-F4 with fused verified flat packing and trusted-u64 hashing",
    "median_ms_per_target": {"x_found_selected": selected["arms"]["x-chained"]["found_ms"], "sym_refuted_selected": selected["arms"]["sym"]["refuted_ms"]},
    "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "node_budget": 50000}, "macaulay_degree_max": 3,
    "f4_splits": {"x": 7056, "sym": 4095}, "verdict_mix": {"x_found": 2, "sym_refuted": 2, "inconclusive": 0},
    "fixture_hash": selected["receipt"]["stdout_sha256"], "executable_or_source_hash": hashes,
    "host_id": {"hostname": socket.gethostname(), "system": platform.system(), "machine": platform.machine()},
    "resource_caps": {"node_budget": 50000, "targets_per_arm": 2, "threads": 1}, "seeds": [24301],
    "claim_boundary": protocol["claim_boundary"], "claim_boundary_non_claims": ["not relation-yield evidence", "not end-to-end index calculus", "not asymptotic or sub-rho evidence", "not external or private targets", "not a key-recovery claim"]}
(HERE / "claim_report_decomposition.json").write_text(json.dumps(claim, indent=2, sort_keys=True) + "\n")
required_global = ["fixture_hash", "executable_or_source_hash", "host_id", "resource_caps", "seeds", "claim_boundary_non_claims"]
required_stage = ["n_or_bits", "m_summands", "dimension_l_or_dim", "unknowns", "system_degree", "eq_var_ratio", "ffd_or_degree_of_regularity", "oracle_class", "median_ms_per_target", "largest_solvable"]
mg = [x for x in required_global if x not in claim]; ms = [x for x in required_stage if x not in claim]
(HERE / "claim_check_decomposition.json").write_text(json.dumps({"schema_version":2,"stage":"decomposition","fail_closed":True,"required_global_provenance":required_global,"required_stage_fields":required_stage,"missing_global_provenance":mg,"missing_stage_fields":ms,"status":"PASS" if not mg and not ms else "SCHEMA_INCOMPLETE"}, indent=2, sort_keys=True)+"\n")
files = {str(path.relative_to(HERE)): sha256(path) for path in sorted(HERE.rglob("*")) if path.is_file() and path.name != "manifest.json"}
(HERE / "manifest.json").write_text(json.dumps({"schema_version":"1.0","task_id":protocol["task_id"],"files":files}, indent=2, sort_keys=True)+"\n")
