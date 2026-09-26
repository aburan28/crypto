#!/usr/bin/env python3
"""One-shot four-arm n53 rank-stage panel; requires reviewed release gate."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PARENT = HERE.parent / "autolab_combined_l384_coldbase_n53_20260925"
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
BASE_GZ = ORBIT / "independent_replay_20260924_codex/base_header.jsonl.gz"
EXE = REPO / "target/release/examples/koblitz_s5_sat_instance"
MAX_RSS = 2 * 1024**3
GLOBAL_S = 1800
ARM_S = 240
AUDIT_S = 600
ARMS = (
    ("native_lex", "native", "lex"),
    ("certified_cyclic", "certified", "target_cyclic_v1"),
    ("certified_lex", "certified", "lex"),
    ("native_cyclic", "native", "target_cyclic_v1"),
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_measure():
    path = PARENT / "run_panel.py"
    spec = importlib.util.spec_from_file_location("n53_rank_rotation_resource_measure", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.measure, module.complete


def run(out: Path):
    import check_protocol
    frozen = check_protocol.preflight(require_release=True)
    assert EXE.is_file()
    out.mkdir(parents=True, exist_ok=False)
    measure, complete = load_measure()
    started = time.monotonic()
    panel = {
        "schema": "n53_target_cyclic_rank_factorial_panel_v1",
        "classification": "RUNNING",
        "checkout_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "protocol_sha256": sha(HERE / "PROTOCOL.md"),
        "source_sha256": sha(REPO / "examples/koblitz_s5_sat_instance.rs"),
        "binary_sha256": sha(EXE),
        "cargo_lock_sha256": sha(PARENT / "Cargo.lock"),
        "target_scalars_sha256": sha(HERE / "target_scalars.txt"),
        "target_points_sha256": sha(HERE / "target_points.jsonl"),
        "certified_base_gzip_sha256": sha(BASE_GZ),
        "host": platform.platform(), "machine": platform.machine(),
        "github_run_id": os.getenv("GITHUB_RUN_ID"),
        "github_run_attempt": os.getenv("GITHUB_RUN_ATTEMPT"),
        "stage_receipts": {},
        "attack_speed_crossover": None,
        "common_operation_unit": None,
    }
    panel_file = out / "panel.json"

    def save():
        panel["elapsed_wall_s"] = time.monotonic() - started
        panel_file.write_text(json.dumps(panel, indent=2, sort_keys=True) + "\n")

    def cap(stage: str, limit: int) -> float:
        remaining = GLOBAL_S - (time.monotonic() - started)
        if remaining <= 0:
            panel["classification"] = "CENSORED_GLOBAL_WALL_BEFORE_" + stage.upper()
            save()
            return 0
        return min(limit, remaining)

    save()
    materialize_started = time.monotonic_ns()
    base_bytes = gzip.decompress(BASE_GZ.read_bytes())
    assert base_bytes.count(b"\n") == 1
    certified = out / "certified_base.jsonl"
    certified.write_bytes(base_bytes)
    panel["certified_base_materialization_ms"] = (time.monotonic_ns() - materialize_started) / 1e6
    panel["certified_base_materialized_sha256"] = sha(certified)
    save()
    executable = [str(EXE), "53", "0", "1", "10", "natural", "1", "2000", "1", "internal"]
    common = {key: val for key, val in os.environ.items() if not key.startswith("KIC_")}
    common.update({
        "KIC_ALGEBRA_ENCODING": "orbit_factorized",
        "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
        "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
        "KIC_ORBIT_REP_ENCODING": "one_hot",
        "KIC_ORBIT_BATCH_ONLY": "1",
        "KIC_ORBIT_INCLUDE_BASE_HEADER": "1",
        "KIC_ORBIT_TARGET_SCALARS": str(HERE / "target_scalars.txt"),
        "KIC_TASK_ID": "TASK-IC-N53-TARGET-CYCLIC-RANK-20260925",
        "RAYON_NUM_THREADS": "1",
    })
    assert common["KIC_ALGEBRA_ENCODING"] == "orbit_factorized"
    for name, base, policy in ARMS:
        seconds = cap(name, ARM_S)
        if not seconds:
            return panel
        env = common | {"KIC_ORBIT_REGULAR_SCAN_POLICY": policy}
        inputs = {"targets": HERE / "target_scalars.txt", "points_validator": HERE / "target_points.jsonl"}
        if base == "certified":
            env["KIC_FACTOR_BASE_JSONL"] = str(certified)
            inputs["certified_base"] = certified
        receipt = measure(executable, env, out / name, "producer", seconds, inputs)
        panel["stage_receipts"][name] = receipt
        save()
        if not complete(receipt):
            panel["classification"] = "CENSORED_OR_INVALID_" + name.upper()
            save()
            return panel
    seconds = cap("audit", AUDIT_S)
    if not seconds:
        return panel
    audit_dir = out / "audit_stage"
    audit_command = [sys.executable, str(HERE / "audit.py"), "--out", str(out)]
    audit_inputs = {name: out / name / "producer.stdout.jsonl" for name, _, _ in ARMS}
    audit_inputs |= {"targets": HERE / "target_scalars.txt", "points": HERE / "target_points.jsonl"}
    receipt = measure(audit_command, common, audit_dir, "audit", seconds, audit_inputs)
    panel["stage_receipts"]["audit"] = receipt
    panel["classification"] = "RANK_STAGE_REPLAYED" if complete(receipt) and (out / "audit.json").is_file() else "CENSORED_OR_INVALID_AUDIT"
    if panel["classification"] == "RANK_STAGE_REPLAYED":
        panel["audit_sha256"] = sha(out / "audit.json")
    panel["measured_arm_wall_ms"] = {name: panel["stage_receipts"][name]["wall_ms"]
                                     for name, _, _ in ARMS}
    save()
    return panel


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.out)
    print(json.dumps({"classification": result["classification"], "output": str(args.out)}, sort_keys=True))
    if result["classification"] != "RANK_STAGE_REPLAYED":
        raise SystemExit(2)
