#!/usr/bin/env python3
"""Run fixed cold same-Q blocks with exact wait4 CPU and per-arm raw traces."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from run_panel import run as run_ir, sha, write_json  # noqa: E402

CONFIG = HERE / "CONFIG.json"


def load_cell(cell_id: str) -> tuple[dict, dict]:
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-compact-ir-cold-gap-v1"
    assert sha(ROOT / config["instruction_config"]) == config["instruction_config_sha256"]
    assert sha(ROOT / config["source_freeze"]) == config["source_freeze_sha256"]
    old = json.loads((ROOT / config["instruction_config"]).read_text())
    matches = [cell for cell in config["cells"] if cell["id"] == cell_id]
    assert len(matches) == 1
    cell = matches[0]
    parent, = [item for item in old["cells"] if item["id"] == cell_id]
    for field, old_field in (("n", "n"), ("L", "L"), ("K", "k"),
                             ("points_sha256", "points_sha256"),
                             ("fixture_sha256", "fixture_sha256")):
        assert cell[field] == parent[old_field], (cell_id, field)
    assert config["arm_order"] == ["off_a", "rho", "off_b"]
    assert config["arm_order_offsets_rule"] == "block_index_mod_3"
    assert (config["rank_seed"], config["s3_batch_window"],
            config["s3_prefilter"], config["rho_dp_bits"],
            config["rho_canonicalization_backend"],
            config["rho_quotient_mode"]) == (
                7, 64, "off", 4, "normal_basis", "signed_frobenius")
    assert config["per_arm_timeout_seconds"] == 900
    assert config["per_arm_address_space_limit_bytes"] == 5 * 1024**3
    assert cell["blocks"] == (20 if cell_id == "n37_L1" else 5)
    return config, cell


def schedule(config: dict, cell: dict, mode: str) -> list[tuple[int, str]]:
    if mode == "smoke":
        return [(0, "off_a")]
    assert mode == "measure"
    order = config["arm_order"]
    return [(block, arm)
            for block in range(cell["blocks"])
            for arm in (order[block % 3:] + order[:block % 3])]


def run(cell_id: str, frozen_root: Path, batch: Path, rho: Path,
        materialization: Path, output: Path, cpu: int | None, mode: str) -> dict:
    config, cell = load_cell(cell_id)
    if mode == "measure":
        assert sys.platform == "linux" and cpu is not None
        assert cpu in os.sched_getaffinity(0)
    else:
        assert mode == "smoke"
    assert not output.exists(), "never overwrite a cold evaluation run"
    output.mkdir(parents=True)
    plan = schedule(config, cell, mode)
    report = {
        "schema": "ecc2k130-compact-ir-cold-gap-run-v1",
        "status": "RUNNING", "mode": mode, "cell": cell,
        "config_sha256": sha(CONFIG),
        "instruction_config_sha256": config["instruction_config_sha256"],
        "source_freeze_sha256": config["source_freeze_sha256"],
        "materialization_sha256": sha(materialization),
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "git_head": subprocess.check_output(
                     ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                 "reserved_cpu": cpu},
        "binaries": {"compact_sha256": sha(batch), "rho_sha256": sha(rho)},
        "plan": [{"block": block, "arm": arm} for block, arm in plan],
        "runs": [],
    }
    write_json(output / "cold_run.json", report)
    try:
        for block, arm in plan:
            policy = "rho" if arm == "rho" else "off"
            subdir = output / f"b{block}_{arm}"
            child = run_ir(cell_id, frozen_root, batch, rho, materialization,
                           subdir, "native", policy, cpu,
                           config["per_arm_timeout_seconds"],
                           config["per_arm_address_space_limit_bytes"],
                           (config["per_arm_address_space_limit_bytes"]
                            if sys.platform == "linux" else None))
            item, = child["runs"]
            assert item["policy"] == policy
            record = {"block": block, "arm": arm, "policy": policy,
                      "subdir": subdir.name, "status": child["status"],
                      "exit_code": item["exit_code"],
                      "stopped_for": item["stopped_for"],
                      "cpu_seconds": (item["child_user_cpu_seconds"] +
                                      item["child_system_cpu_seconds"]),
                      "wall_seconds": item["elapsed_under_backend_seconds_not_a_cost"],
                      "max_rss_kib_linux": item["child_max_rss_kib_linux"]}
            report["runs"].append(record)
            write_json(output / "cold_run.json", report)
            if child["status"] != "SMOKE_PASS":
                report["status"] = "FAIL"
                break
        else:
            report["status"] = "SMOKE_PASS" if mode == "smoke" else "PASS"
    except BaseException as error:
        report["status"] = "FAIL"
        report["error_type"] = type(error).__name__
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
        write_json(output / "cold_run.json", report)
        raise
    write_json(output / "cold_run.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", required=True)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--materialization", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--mode", choices=("smoke", "measure"), required=True)
    args = parser.parse_args()
    result = run(args.cell, args.frozen_root.resolve(), args.batch.resolve(),
                 args.rho.resolve(), args.materialization.resolve(),
                 args.out.resolve(), args.cpu, args.mode)
    if result["status"] not in ("PASS", "SMOKE_PASS"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
