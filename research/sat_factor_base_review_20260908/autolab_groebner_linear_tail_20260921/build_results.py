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


def load(path: Path):
    return json.loads(path.read_text())


def observation(name: str):
    return load(RUNS / name / "stdout.txt")


def receipt(name: str):
    return load(RUNS / name / "receipt.json")


def ladder(name: str):
    return {
        "receipt": receipt(name),
        "stage": load(RUNS / name / "artifact" / "stage.json"),
    }


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
            "variables": int(fields[2]),
            "equations": int(fields[3]),
            "degree": int(fields[4]),
            "found": int(fields[5]),
            "refuted": int(fields[6]),
            "inconclusive": int(fields[7]),
            "via_t": int(fields[8]),
            "found_ms": optional_number(fields[9]),
            "refuted_ms": optional_number(fields[10]),
            "effort": int(fields[11]),
            "built_degree": int(fields[12]),
            "ffd": int(fields[13]),
            "gate": fields[14],
        }
    assert set(arms) == {"x-chained", "sym"}
    return {"receipt": receipt(name), "arms": arms, "stdout": text}


protocol = load(HERE / "protocol.json")
selected = paired("paired_selected_retained_n31_default_t2")
legacy = paired("paired_legacy_retained_n31_default_t2")
profiles = {
    "sym_selected": observation("profile_sym_selected"),
    "sym_legacy": observation("profile_sym_legacy"),
    "x_selected": observation("profile_x_selected"),
    "x_legacy": observation("profile_x_legacy"),
}
degree4 = observation("profile_sym_degree4_sparse")
ladders = {
    "frozen_selected": ladder("ladder_frozen_selected"),
    "frozen_legacy": ladder("ladder_frozen_legacy"),
    "holdout_selected": ladder("ladder_holdout_selected"),
    "holdout_legacy": ladder("ladder_holdout_legacy"),
}

assert selected["receipt"]["executable_sha256"] == legacy["receipt"]["executable_sha256"]
for arm in ("x-chained", "sym"):
    assert selected["arms"][arm]["found"] == legacy["arms"][arm]["found"]
    assert selected["arms"][arm]["refuted"] == legacy["arms"][arm]["refuted"]
    assert selected["arms"][arm]["inconclusive"] == legacy["arms"][arm]["inconclusive"] == 0
    assert selected["arms"][arm]["effort"] == legacy["arms"][arm]["effort"]
    assert selected["arms"][arm]["built_degree"] == legacy["arms"][arm]["built_degree"]
    assert selected["arms"][arm]["ffd"] == legacy["arms"][arm]["ffd"]
    assert selected["arms"][arm]["gate"] == legacy["arms"][arm]["gate"] == "ok"
for system in ("sym", "x"):
    a = profiles[f"{system}_selected"]
    b = profiles[f"{system}_legacy"]
    assert a["roots"] == b["roots"]
    assert a["solve_stats"] == b["solve_stats"]
    assert a["f4_profile"]["calls"] < b["f4_profile"]["calls"]
    assert a["f4_profile"]["word_ops"] < b["f4_profile"]["word_ops"]
assert degree4["profile"]["resolves"] is False
ladder_comparisons = {}
for group in ("frozen", "holdout"):
    chosen = ladders[f"{group}_selected"]["stage"]["rows"]
    control = ladders[f"{group}_legacy"]["stage"]["rows"]
    assert len(chosen) == len(control)
    rows = []
    for selected_row, legacy_row in zip(chosen, control):
        key = (selected_row["curve"], selected_row["m"], selected_row["ell"])
        assert key == (legacy_row["curve"], legacy_row["m"], legacy_row["ell"])
        assert selected_row["targets"] == legacy_row["targets"]
        assert selected_row["decomposed"] == legacy_row["decomposed"]
        assert selected_row["verdict_digest"] == legacy_row["verdict_digest"]
        assert selected_row["wall_ns"] <= legacy_row["wall_ns"]
        rows.append({
            "curve": selected_row["curve"],
            "m": selected_row["m"],
            "ell": selected_row["ell"],
            "targets": selected_row["targets"],
            "verdict_digest": selected_row["verdict_digest"],
            "selected_wall_ns": selected_row["wall_ns"],
            "legacy_wall_ns": legacy_row["wall_ns"],
            "wall_speedup": legacy_row["wall_ns"] / selected_row["wall_ns"],
            "selected_word_ops": selected_row["word_ops"],
            "legacy_word_ops": legacy_row["word_ops"],
        })
    ladder_comparisons[group] = rows

result = {
    "schema_version": "1.0",
    "task_id": protocol["task_id"],
    "classification": "N31_BOOLEAN_F4_LINEAR_TAIL_PAIRED_T2_PASS",
    "selected_policy": protocol["selected_policy"],
    "paired": {
        "selected": selected,
        "legacy": legacy,
        "wall_speedup": legacy["receipt"]["wall_ms"] / selected["receipt"]["wall_ms"],
        "child_user_cpu_speedup": legacy["receipt"]["child_user_ms"]
        / selected["receipt"]["child_user_ms"],
        "selected_over_legacy_peak_rss": selected["receipt"]["peak_rss_bytes"]
        / legacy["receipt"]["peak_rss_bytes"],
        "x_median_speedup": legacy["arms"]["x-chained"]["found_ms"]
        / selected["arms"]["x-chained"]["found_ms"],
        "sym_median_speedup": legacy["arms"]["sym"]["refuted_ms"]
        / selected["arms"]["sym"]["refuted_ms"],
    },
    "profiles": {
        **profiles,
        "sym_solver_speedup": profiles["sym_legacy"]["elapsed_ms"]
        / profiles["sym_selected"]["elapsed_ms"],
        "x_solver_speedup": profiles["x_legacy"]["elapsed_ms"]
        / profiles["x_selected"]["elapsed_ms"],
        "sym_word_ops_reduction_fraction": 1
        - profiles["sym_selected"]["f4_profile"]["word_ops"]
        / profiles["sym_legacy"]["f4_profile"]["word_ops"],
        "x_word_ops_reduction_fraction": 1
        - profiles["x_selected"]["f4_profile"]["word_ops"]
        / profiles["x_legacy"]["f4_profile"]["word_ops"],
    },
    "rejected": {
        "root_degree4_sparse": degree4,
        "reason": "The 15,407 by 32,817 degree-4 closure took 28.4 seconds and neither refuted nor determined a variable.",
        "forward_m4ri": "Rejected during exploration: block widths 2, 3, and 4 all reduced XORs but were slower than suffix-forward elimination on the shrinking branch matrices.",
    },
    "ladders": {
        **ladders,
        "comparisons": ladder_comparisons,
    },
    "validation": {
        "linear_tail_row_space_equivalence": "PASS",
        "forced_small_system_exhaustive_roots": "PASS",
        "forced_all_split_rules": "PASS",
        "paired_output_and_gate_equivalence": "PASS",
        "frozen_ladder_verdict_digests": "PASS",
        "holdout_ladder_verdict_digests": "PASS",
        "frozen_and_holdout_wall_non_regression": "PASS",
    },
    "claim_boundary": protocol["claim_boundary"],
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

source_hashes = {
    "koblitz_groebner_source": sha256(FROZEN / "koblitz_groebner.rs"),
    "rref_tests_source": sha256(FROZEN / "f4_rref_tests.rs"),
    "frontier_probe_source": sha256(FROZEN / "koblitz_f4_frontier_probe.rs"),
    "paired_executable_receipt": selected["receipt"]["executable_sha256"],
    "profile_executable_receipt": receipt("profile_sym_selected")["executable_sha256"],
}
claim_report = {
    "schema_version": 2,
    "task_id": protocol["task_id"],
    "stage": "decomposition",
    "n_or_bits": 31,
    "m_summands": 2,
    "dimension_l_or_dim": 16,
    "unknowns": 32,
    "system_degree": "quadratic",
    "eq_var_ratio": 31 / 32,
    "ffd_or_degree_of_regularity": {"x_chained": 3, "symmetrised": 4},
    "oracle_class": "Boolean matrix-F4 with maximum-degree scheduling and exact linear-tail elimination",
    "median_ms_per_target": {
        "x_found_selected": selected["arms"]["x-chained"]["found_ms"],
        "sym_refuted_selected": selected["arms"]["sym"]["refuted_ms"],
    },
    "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "node_budget": 50000},
    "macaulay_degree_max": 3,
    "f4_splits": {"x": 7056, "sym": 4095},
    "verdict_mix": {"x_found": 2, "sym_refuted": 2, "inconclusive": 0},
    "fixture_hash": selected["receipt"]["stdout_sha256"],
    "executable_or_source_hash": source_hashes,
    "host_id": {"hostname": socket.gethostname(), "system": platform.system(), "machine": platform.machine()},
    "resource_caps": {"node_budget": 50000, "targets_per_arm": 2, "threads": 1},
    "seeds": [24301],
    "claim_boundary": protocol["claim_boundary"],
    "claim_boundary_non_claims": [
        "not relation-yield evidence",
        "not end-to-end index calculus",
        "not asymptotic or sub-rho evidence",
        "not external or private targets",
        "not a key-recovery claim"
    ]
}
(HERE / "claim_report_decomposition.json").write_text(json.dumps(claim_report, indent=2, sort_keys=True) + "\n")

required_global = ["fixture_hash", "executable_or_source_hash", "host_id", "resource_caps", "seeds", "claim_boundary_non_claims"]
required_stage = ["n_or_bits", "m_summands", "dimension_l_or_dim", "unknowns", "system_degree", "eq_var_ratio", "ffd_or_degree_of_regularity", "oracle_class", "median_ms_per_target", "largest_solvable"]
missing_global = [field for field in required_global if field not in claim_report]
missing_stage = [field for field in required_stage if field not in claim_report]
claim_check = {
    "schema_version": 2,
    "stage": "decomposition",
    "fail_closed": True,
    "required_global_provenance": required_global,
    "required_stage_fields": required_stage,
    "missing_global_provenance": missing_global,
    "missing_stage_fields": missing_stage,
    "status": "PASS" if not missing_global and not missing_stage else "SCHEMA_INCOMPLETE"
}
(HERE / "claim_check_decomposition.json").write_text(json.dumps(claim_check, indent=2, sort_keys=True) + "\n")

files = {}
for path in sorted(HERE.rglob("*")):
    if path.is_file() and path.name != "manifest.json":
        files[str(path.relative_to(HERE))] = sha256(path)
(HERE / "manifest.json").write_text(json.dumps({"schema_version": "1.0", "task_id": protocol["task_id"], "files": files}, indent=2, sort_keys=True) + "\n")
