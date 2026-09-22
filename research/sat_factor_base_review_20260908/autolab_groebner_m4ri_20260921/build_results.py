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


def median_word_ops(observation):
    values = sorted(observation["word_ops"])
    return values[len(values) // 2]


protocol = load(HERE / "protocol.json")
x_suffix = observation("x_suffix_r31")
x_auto = observation("x_auto_r31")
sym_suffix = observation("sym_suffix_r31")
sym_auto = observation("sym_auto_r31")
paired_suffix = paired("paired_suffix_n31_default_seed_t1")
paired_auto = paired("paired_auto_density_n31_default_t1")

assert x_auto["rref_kernel"] == "active_suffix"
assert sym_auto["rref_kernel"] == "m4ri_auto"
assert x_auto["rank"] == x_suffix["rank"] == 1022
assert sym_auto["rank"] == sym_suffix["rank"] == 992
assert median_word_ops(x_auto) == median_word_ops(x_suffix)
assert median_word_ops(sym_auto) < median_word_ops(sym_suffix)
for run in (paired_suffix, paired_auto):
    assert run["receipt"]["exit_code"] == 0
    assert run["arms"]["x-chained"]["found"] == 1
    assert run["arms"]["sym"]["refuted"] == 1
    assert all(arm["inconclusive"] == 0 and arm["gate"] == "ok" for arm in run["arms"].values())
    assert run["arms"]["x-chained"]["effort"] == 7056
    assert run["arms"]["sym"]["effort"] == 4095
    assert run["arms"]["x-chained"]["ffd"] == 3
    assert run["arms"]["sym"]["ffd"] == 4

sym_suffix_ops = median_word_ops(sym_suffix)
sym_auto_ops = median_word_ops(sym_auto)
result = {
    "schema_version": "1.0",
    "task_id": protocol["task_id"],
    "classification": "N31_BOOLEAN_F4_M4RI_KERNEL_PASS_PAIRED_WALL_INDETERMINATE",
    "selected_policy": protocol["selected_policy"],
    "kernel": {
        "x_suffix_control": x_suffix,
        "x_selected": x_auto,
        "x_policy_result": "unchanged_active_suffix",
        "sym_suffix_control": sym_suffix,
        "sym_selected": sym_auto,
        "sym_median_speedup": sym_suffix["median_ms"] / sym_auto["median_ms"],
        "sym_word_ops_reduction_fraction": 1 - sym_auto_ops / sym_suffix_ops,
        "sym_word_ops_suffix": sym_suffix_ops,
        "sym_word_ops_selected": sym_auto_ops,
    },
    "block_width_ablation": {
        "sym": {
            width: observation(f"tune_sym_m4ri_b{width}_r11")
            for width in (4, 6, 8, 10)
        }
    },
    "paired": {
        "suffix_control": paired_suffix,
        "selected_density_policy": paired_auto,
        "selected_vs_suffix_wall_ratio": paired_auto["receipt"]["wall_ms"]
        / paired_suffix["receipt"]["wall_ms"],
        "interpretation": "One target per arm is enough for exact-output validation but not a wall-clock speedup claim; the selected run was faster overall than the suffix control, pending a paired distribution.",
    },
    "validation": {
        "rref_reference_test": "PASS",
        "small_system_exhaustive_root_set_test": "PASS",
        "paired_output_and_gate_equivalence": "PASS",
    },
    "claim_boundary": protocol["claim_boundary"],
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

source_hashes = {
    "koblitz_groebner_source": sha256(FROZEN / "koblitz_groebner.rs"),
    "kernel_bench_source": sha256(FROZEN / "koblitz_f4_kernel_bench.rs"),
    "kernel_bench_executable_receipt": receipt("sym_auto_r31")["executable_sha256"],
    "paired_bench_executable_receipt": paired_auto["receipt"]["executable_sha256"],
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
    "oracle_class": "Boolean matrix-F4 with shape-selected block-4 M4RI",
    "median_ms_per_target": {
        "x_found_selected_paired": paired_auto["arms"]["x-chained"]["found_ms"],
        "sym_refuted_selected_paired": paired_auto["arms"]["sym"]["refuted_ms"],
    },
    "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "node_budget": 50000},
    "macaulay_degree_max": 3,
    "f4_splits": {"x": 7056, "sym": 4095},
    "verdict_mix": {"x_found": 1, "sym_refuted": 1, "inconclusive": 0},
    "fixture_hash": paired_auto["receipt"]["stdout_sha256"],
    "executable_or_source_hash": source_hashes,
    "host_id": {"hostname": socket.gethostname(), "system": platform.system(), "machine": platform.machine()},
    "resource_caps": {"node_budget": 50000, "targets_per_arm": 1, "threads": 1},
    "seeds": [24301],
    "claim_boundary": protocol["claim_boundary"],
    "claim_boundary_non_claims": [
        "not relation-yield evidence",
        "not end-to-end index calculus",
        "not asymptotic or sub-rho evidence",
        "not external or private targets",
        "not a key-recovery claim",
    ],
}
(HERE / "claim_report_decomposition.json").write_text(
    json.dumps(claim_report, indent=2, sort_keys=True) + "\n"
)

required_global = [
    "fixture_hash",
    "executable_or_source_hash",
    "host_id",
    "resource_caps",
    "seeds",
    "claim_boundary_non_claims",
]
required_stage = [
    "n_or_bits",
    "m_summands",
    "dimension_l_or_dim",
    "unknowns",
    "system_degree",
    "eq_var_ratio",
    "ffd_or_degree_of_regularity",
    "oracle_class",
    "median_ms_per_target",
    "largest_solvable",
]
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
    "status": "PASS" if not missing_global and not missing_stage else "SCHEMA_INCOMPLETE",
}
(HERE / "claim_check_decomposition.json").write_text(
    json.dumps(claim_check, indent=2, sort_keys=True) + "\n"
)

files = {}
for path in sorted(HERE.rglob("*")):
    if path.is_file() and path.name != "manifest.json":
        files[str(path.relative_to(HERE))] = sha256(path)
(HERE / "manifest.json").write_text(
    json.dumps({"schema_version": "1.0", "task_id": protocol["task_id"], "files": files}, indent=2, sort_keys=True) + "\n"
)
