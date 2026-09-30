#!/usr/bin/env python3
"""Replay all cold ranks and logs, then adjudicate fixed paired CPU controls."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from run_panel import sha  # noqa: E402
from verify_panel import identity, rows, verify as verify_ir  # noqa: E402
from run_cold import CONFIG, load_cell, schedule  # noqa: E402

T_CRIT_95 = {5: 2.7764451051977987, 20: 2.093024054408263}


def paired_interval(values: list[float]) -> tuple[float, list[float]]:
    assert len(values) in T_CRIT_95 and all(value > 0 for value in values)
    logs = [math.log(value) for value in values]
    mean = statistics.mean(logs)
    half = T_CRIT_95[len(values)] * statistics.stdev(logs) / math.sqrt(len(values))
    return statistics.median(values), [math.exp(mean - half), math.exp(mean + half)]


def compact_identity(subdir: Path) -> tuple:
    base, = rows(subdir / "off.base.jsonl")
    ranks = rows(subdir / "off.rank.jsonl")
    targets = rows(subdir / "off.target.jsonl")
    return base, ranks, tuple(identity(target) for target in targets)


def verify(cell_id: str, run_dir: Path, mode: str) -> dict:
    config, cell = load_cell(cell_id)
    report = json.loads((run_dir / "cold_run.json").read_text())
    assert report["schema"] == "ecc2k130-compact-ir-cold-gap-run-v1"
    assert report["mode"] == mode and report["cell"] == cell
    assert report["config_sha256"] == sha(CONFIG)
    assert report["instruction_config_sha256"] == config["instruction_config_sha256"]
    assert report["source_freeze_sha256"] == config["source_freeze_sha256"]
    assert report["status"] == ("SMOKE_PASS" if mode == "smoke" else "PASS")
    expected = schedule(config, cell, mode)
    assert report["plan"] == [{"block": block, "arm": arm} for block, arm in expected]
    assert [(item["block"], item["arm"]) for item in report["runs"]] == expected
    reserved_cpu = report["host"]["reserved_cpu"]
    if mode == "measure":
        assert report["host"]["machine"] == "x86_64" and reserved_cpu is not None
    measurements: dict[int, dict] = {}
    checked: list[dict] = []
    witness: dict[tuple[int, str], tuple] = {}
    for record in report["runs"]:
        block, arm = record["block"], record["arm"]
        policy = "rho" if arm == "rho" else "off"
        assert record["policy"] == policy and record["status"] == "SMOKE_PASS"
        assert record["subdir"] == f"b{block}_{arm}"
        subdir = run_dir / record["subdir"]
        child = json.loads((subdir / "run.json").read_text())
        assert child["status"] == "SMOKE_PASS" and child["backend"] == "native"
        assert child["sequence"] == [policy]
        assert child["host"]["git_head"] == report["host"]["git_head"]
        assert child["host"]["reserved_cpu"] == reserved_cpu
        assert child["binaries"] == report["binaries"]
        assert child["limits"] == {
            "timeout_seconds_per_arm": config["per_arm_timeout_seconds"],
            "rss_limit_bytes_per_arm": config["per_arm_address_space_limit_bytes"],
            "address_space_limit_bytes_per_arm": (
                config["per_arm_address_space_limit_bytes"]
                if child["host"]["platform"].startswith("Linux") else None)}
        item, = child["runs"]
        assert item["exit_code"] == record["exit_code"] == 0
        assert item["stopped_for"] == record["stopped_for"] is None
        assert item["elapsed_under_backend_seconds_not_a_cost"] == record["wall_seconds"]
        assert (item["child_user_cpu_seconds"] + item["child_system_cpu_seconds"]
                == record["cpu_seconds"] > 0)
        assert item["child_max_rss_kib_linux"] == record["max_rss_kib_linux"]
        if mode == "measure":
            assert item["child_max_rss_kib_linux"] > 0
            assert item["child_max_rss_kib_linux"] * 1024 < config[
                "per_arm_address_space_limit_bytes"]
            assert item["command"][:3] == ["taskset", "-c", str(reserved_cpu)]
        replay = verify_ir(cell_id, subdir, smoke=True)
        assert replay["status"] == "SMOKE_PASS"
        assert replay["details"][policy]["targets_verified"] == cell["L"]
        checked.append({"block": block, "arm": arm,
                        "targets_verified": cell["L"],
                        "run_sha256": sha(subdir / "run.json")})
        if policy == "off":
            witness[(block, arm)] = compact_identity(subdir)
        measurements.setdefault(block, {})[arm] = {
            "cpu_seconds": record["cpu_seconds"],
            "wall_seconds": record["wall_seconds"],
            "rss_kib": record["max_rss_kib_linux"]}
    if mode == "smoke":
        return {"status": "SMOKE_PASS", "cell": cell_id, "checks": checked,
                "timing_eligible": None, "paired": None}
    isolation, = rows(run_dir / "isolation.jsonl")
    assert isolation["schema"] == "isolated-bench/1" and isolation["mode"] == "reserve"
    assert isolation["exit_status"] == 0
    assert reserved_cpu in isolation["reserved_cpus"]
    for block in range(cell["blocks"]):
        assert set(measurements[block]) == {"off_a", "rho", "off_b"}
        assert witness[(block, "off_a")] == witness[(block, "off_b")]
    aa_values, ic_values, wall_values = [], [], []
    for block in range(cell["blocks"]):
        data = measurements[block]
        aa_values.append(data["off_b"]["cpu_seconds"] / data["off_a"]["cpu_seconds"])
        off_geo = math.sqrt(data["off_a"]["cpu_seconds"] *
                            data["off_b"]["cpu_seconds"])
        ic_values.append(off_geo / data["rho"]["cpu_seconds"])
        wall_geo = math.sqrt(data["off_a"]["wall_seconds"] *
                             data["off_b"]["wall_seconds"])
        wall_values.append(wall_geo / data["rho"]["wall_seconds"])
    aa_median, aa_interval = paired_interval(aa_values)
    ic_median, ic_interval = paired_interval(ic_values)
    wall_median, wall_interval = paired_interval(wall_values)
    aa_valid = 0.9 <= aa_median <= 1.1 and aa_interval[0] <= 1 <= aa_interval[1]
    uncontended = isolation["contended_samples"] == 0
    return {"status": "PASS", "cell": cell_id, "n": cell["n"], "L": cell["L"],
            "K": cell["K"], "blocks": cell["blocks"],
            "checked": checked, "measurements": measurements,
            "paired": {"off_b_over_off_a": {"median": aa_median,
                                            "interval_95pct": aa_interval},
                       "off_geo_over_rho_cpu": {"median": ic_median,
                                                 "interval_95pct": ic_interval},
                       "off_geo_over_rho_wall": {"median": wall_median,
                                                  "interval_95pct": wall_interval}},
            "aa_valid": aa_valid, "uncontended": uncontended,
            "timing_eligible": aa_valid and uncontended,
            "method_crossover": None,
            "evidence_class": "same_Q_full_process_cold_CPU_not_disjoint_holdout"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--mode", choices=("smoke", "measure"), required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a cold replay receipt"
    try:
        receipt = verify(args.cell, args.run_dir.resolve(), args.mode)
    except BaseException as error:
        receipt = {"status": "FAIL", "error_type": type(error).__name__,
                   "error": str(error), "traceback": traceback.format_exc()}
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in receipt.items()
                      if key not in ("checked", "measurements")}, sort_keys=True))


if __name__ == "__main__":
    main()
