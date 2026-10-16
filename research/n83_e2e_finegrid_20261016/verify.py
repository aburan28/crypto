#!/usr/bin/env python3
"""Independently replay the complete N83 adjacent-K cold screen."""

import argparse
import importlib.util
import json
import statistics
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_finegrid_runner", HERE / "run.py")
holdout = module("n83_finegrid_holdout_verify",
                 ROOT / "research/n83_e2e_holdout_20261016/verify.py")
dynamic = holdout.dynamic


def read(path):
    return json.loads(path.read_text())


def repeatability(by_arm, k):
    a, b = (by_arm.get(f"k{k}_{suffix}") for suffix in ("a", "b"))
    if a is None or b is None:
        return {"status": "INCOMPLETE"}
    if (a["base_hash"] != b["base_hash"] or
            a["root_table_entries"] != b["root_table_entries"]):
        return {"status": "BASE_OR_INDEX_MISMATCH"}
    times = [a["cold_process_ms_exploratory"],
             b["cold_process_ms_exploratory"]]
    spread = abs(times[0] - times[1]) / statistics.mean(times)
    return {"status": "PASS_5_PERCENT" if spread <= 0.05
            else "NOISY_5_PERCENT", "cold_process_ms": times,
            "relative_spread": spread}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, point, scalar = dynamic.coldver.frozen_point()
    if hashes != runner.frozen_inputs():
        raise ValueError("frozen source or public point differs")
    host = read(output / "host.json")
    source = read(output / "source.json")
    preflight = read(output / "preflight.json")
    if (host["schema"] != "n83.e2e-finegrid-host/v1"
            or host["git_sha"] is None
            or preflight["status"] != "PASS_CAPACITY_PREFLIGHT"
            or source["schema"] != "n83.e2e-finegrid-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(HERE / "run.py")
            or source["target_sha256"] != runner.TARGET_SHA256
            or source["arm_order"] != [name for name, _ in runner.ARMS]
            or source["k_grid"] != list(runner.K_GRID)
            or source["control_ks"] != list(runner.CONTROL_KS)
            or source["rho_seed"] != runner.RHO_SEED):
        raise ValueError("fine-grid source, host or schedule differs")
    ordered = source["arm_order"] + ["rho"]
    topology = host.get("topology", {})
    topology_status = ("PASS" if topology.get("lscpu_json")
                       and len(topology.get("visible_cpu_topology", {})) >= 4
                       and (output / "host_telemetry_after.json").exists()
                       else "MISSING")
    telemetry = {name: holdout.telemetry_receipt(output, name)
                 for name in ordered}
    telemetry_status = ("PASS" if topology_status == "PASS" and all(
        item["status"] == "PASS" for item in telemetry.values())
                        else "INCOMPLETE")
    finished = []
    for name in ordered:
        path = output / "runs" / name / "usage.json"
        if not path.exists():
            break
        finished.append(datetime.fromisoformat(read(path)["finished_utc"]))
    order_status = ("PASS" if len(finished) == len(ordered)
                    and all(a < b for a, b in zip(finished, finished[1:]))
                    else "INCOMPLETE_OR_OUT_OF_ORDER")
    arms = []
    for name, k in runner.ARMS:
        try:
            arms.append(dynamic.verify_ic(output, name, k, "root11_dynamic",
                                          point, scalar))
        except Exception as error:
            arms.append({"arm": name, "k": k,
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    try:
        rho = dynamic.coldver.verify_rho(output, point, scalar)
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = {arm["arm"]: arm for arm in arms
                if arm["status"] == "PASS_INDEPENDENT_REPLAY"
                and arm["resources"]["status"] == "PASS"}
    aa = {str(k): repeatability(accepted, k) for k in runner.CONTROL_KS}
    evidence_valid = (order_status == "PASS" and telemetry_status == "PASS"
                      and all(arm["status"] in
                              ("PASS_INDEPENDENT_REPLAY", "TIMEOUT")
                              for arm in arms)
                      and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                      and all(item["status"] != "BASE_OR_INDEX_MISMATCH"
                              for item in aa.values()))
    admitted = (evidence_valid
                and len(accepted) == len(runner.ARMS)
                and rho["resources"]["status"] == "PASS"
                and all(item["status"] == "PASS_5_PERCENT"
                        for item in aa.values()))
    by_k = {}
    if admitted:
        for k in runner.K_GRID:
            names = (f"k{k}_a", f"k{k}_b") if k in runner.CONTROL_KS else (f"k{k}",)
            by_k[str(k)] = statistics.mean(
                accepted[name]["cold_process_ms_exploratory"] for name in names)
    analysis = {
        "schema": "n83.e2e-finegrid-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "exposed_point_complete_adjacent_k_cold_screen",
        "target_count": 1, "k_grid": list(runner.K_GRID),
        "arm_order": source["arm_order"],
        "status": "SCREEN_COMPLETE" if admitted else "BOUNDED_UNKNOWN",
        "evidence_status": "PASS_INDEPENDENT_REPLAY" if evidence_valid else
                           "VERIFICATION_FAILED",
        "selection_admitted": admitted, "selected_candidate_id": None,
        "descriptive_cold_process_ms_by_k": by_k,
        "descriptive_order": sorted(by_k, key=by_k.get),
        "aa_controls": aa,
        "execution_order_status": order_status,
        "topology_status": topology_status,
        "telemetry_status": telemetry_status,
        "host_telemetry_by_arm": telemetry,
        "arms": arms, "rho": rho,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if admitted else 2


if __name__ == "__main__":
    raise SystemExit(main())
