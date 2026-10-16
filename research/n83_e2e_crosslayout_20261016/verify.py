#!/usr/bin/env python3
"""Replay the N83 cold cross-layout screen and its matched controls."""

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


runner = module("n83_crosslayout_runner", HERE / "run.py")
holdout = module("n83_crosslayout_holdout_verify",
                 ROOT / "research/n83_e2e_holdout_20261016/verify.py")
dynamic = holdout.dynamic


def read(path):
    return json.loads(path.read_text())


def aa_control(by_arm, prefix):
    a, b = (by_arm.get(prefix + suffix) for suffix in ("_a", "_b"))
    if a is None or b is None:
        return {"status": "INCOMPLETE"}
    if (a["base_hash"] != b["base_hash"] or
            a["root_table_entries"] != b["root_table_entries"]):
        return {"status": "BASE_OR_INDEX_MISMATCH"}
    times = [a["cold_process_ms_exploratory"],
             b["cold_process_ms_exploratory"]]
    relative = abs(times[0] - times[1]) / statistics.mean(times)
    return {"status": "PASS_5_PERCENT" if relative <= 0.05
            else "NOISY_5_PERCENT", "cold_process_ms": times,
            "relative_spread": relative}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, point, scalar = dynamic.coldver.frozen_point()
    if hashes != runner.frozen_inputs():
        raise ValueError("frozen source or exposed target differs")
    host = read(output / "host.json")
    source = read(output / "source.json")
    preflight = read(output / "preflight.json")
    if (host["schema"] != "n83.e2e-crosslayout-host/v1"
            or host["git_sha"] is None
            or preflight["status"] != "PASS_CAPACITY_PREFLIGHT"
            or source["schema"] != "n83.e2e-crosslayout-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(HERE / "run.py")
            or source["target_sha256"] != runner.TARGET_SHA256
            or source["arm_order"] != [name for name, _, _ in runner.ARMS]
            or source["arm_layouts"] != {name: layout for name, _, layout
                                        in runner.ARMS}
            or source["rho_seed"] != runner.RHO_SEED):
        raise ValueError("cross-layout source, host or schedule differs")
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
    for name, k, layout in runner.ARMS:
        try:
            arms.append(dynamic.verify_ic(output, name, k, layout,
                                          point, scalar))
        except Exception as error:
            arms.append({"arm": name, "k": k, "layout": layout,
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    try:
        rho = dynamic.coldver.verify_rho(output, point, scalar)
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = {arm["arm"]: arm for arm in arms
                if arm["status"] == "PASS_INDEPENDENT_REPLAY"
                and arm["resources"]["status"] == "PASS"}
    aa = {layout: aa_control(accepted, "k1104" + layout)
          for layout in ("root11", "root22")}
    identity = {}
    for k in (1104, 1120):
        rows = [arm for arm in accepted.values() if arm["k"] == k]
        hashes_by_layout = {layout: {arm["base_hash"] for arm in rows
                                     if arm["layout"] == layout}
                            for layout in ("root11_dynamic", "root22_dynamic")}
        identity[str(k)] = {"base_hashes": {key: sorted(value)
                                            for key, value in hashes_by_layout.items()},
                            "status": "PASS" if len(set().union(
                                *hashes_by_layout.values())) == 1
                            and all(hashes_by_layout.values()) else "MISMATCH"}
    evidence_valid = (order_status == "PASS" and telemetry_status == "PASS"
                      and all(arm["status"] in
                              ("PASS_INDEPENDENT_REPLAY", "TIMEOUT")
                              for arm in arms)
                      and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                      and all(item["status"] != "MISMATCH"
                              for item in identity.values()))
    admitted = (evidence_valid
                and len(accepted) == len(runner.ARMS)
                and rho["resources"]["status"] == "PASS"
                and all(item["status"] == "PASS_5_PERCENT"
                        for item in aa.values())
                and all(item["status"] == "PASS" for item in identity.values()))
    by_pipeline = {}
    if admitted:
        for arm in arms:
            key = f"k{arm['k']}{arm['layout']}"
            by_pipeline.setdefault(key, []).append(
                arm["cold_process_ms_exploratory"])
        by_pipeline = {key: statistics.mean(values)
                       for key, values in by_pipeline.items()}
    analysis = {
        "schema": "n83.e2e-crosslayout-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "exposed_point_current_source_cold_layout_screen",
        "target_count": 1, "arm_order": source["arm_order"],
        "status": "SCREEN_COMPLETE" if admitted else "BOUNDED_UNKNOWN",
        "evidence_status": "PASS_INDEPENDENT_REPLAY" if evidence_valid else
                           "VERIFICATION_FAILED",
        "selection_admitted": admitted,
        "selected_candidate_id": None,
        "descriptive_cold_process_ms_by_pipeline": by_pipeline,
        "descriptive_order": sorted(by_pipeline, key=by_pipeline.get),
        "aa_controls": aa, "cross_layout_base_identity": identity,
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
