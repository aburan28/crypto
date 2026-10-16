#!/usr/bin/env python3
"""Independent N83 K3114/K3072 A/A, relation and scalar replay."""

import argparse
import importlib.util
import json
import math
import statistics
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNNER = HERE / "run.py"
DYNAMIC_VERIFIER = ROOT / "research/n83_root11_dynamic_20261016/verify.py"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_online_range_runner", RUNNER)
dynamic = module("n83_online_range_dynamic_verifier", DYNAMIC_VERIFIER)
group = dynamic.group


def one_json(path):
    return json.loads(path.read_text())


def frozen_point():
    hashes = runner.frozen_inputs()
    point_data = one_json(runner.cold.TARGETS)
    answer = one_json(runner.cold.ANSWERS)
    point = group.point_from_json(point_data)
    scalar = int(answer["scalar"])
    if (answer["fixture_index"] != 0 or answer["point"] != point_data
            or not group.on_curve(point)
            or group.multiply(point, group.R) is not None
            or group.multiply(group.G, scalar) != point):
        raise ValueError("frozen exposed point fails independent subgroup replay")
    return hashes, point, scalar


def verify_rho(output, point, scalar):
    result = dynamic.coldver.verify_rho(output, point, scalar)
    if result["status"] != "PASS_INDEPENDENT_REPLAY":
        return result
    receipt = one_json(output / "runs" / "rho" / "receipt.json")
    if (receipt["walk_seed"] != int(runner.cold.RHO_SEED)
            or receipt["target_count"] != 1
            or receipt["known_answer_input"] is not False
            or receipt["validation_sidecar_opened"] is not False
            or receipt["group_verified"] is not True
            or receipt["public_input_sha256"] != runner.parent.digest(
                runner.cold.TARGETS)
            or receipt["compiled_source_sha256"][
                "examples/n83_rho_optimization_eval.rs"] != runner.frozen_inputs()[
                    "examples/n83_rho_optimization_eval.rs"]):
        raise ValueError("same-point rho input or source receipt differs")
    return result


def online_control(rows, prefix):
    names = [f"{prefix}_{suffix}" for suffix in "ab"]
    if any(name not in rows for name in names):
        return {"status": "INCOMPLETE", "arms": names}
    identity = {(rows[name]["base_hash"], rows[name]["root_table_entries"])
                for name in names}
    if len(identity) != 1:
        return {"status": "BASE_OR_INDEX_MISMATCH", "arms": names}
    times = [rows[name]["one_target_online_ms"] for name in names]
    difference = abs(times[0] - times[1])
    spread = difference / statistics.mean(times)
    return {"status": "PASS_10_PERCENT" if spread <= 0.10 else "NOISY",
            "arms": names, "online_ms": times,
            "mean_online_ms": statistics.mean(times),
            "pair_relative_differences": [spread],
            "pair_absolute_differences_ms": [difference],
            "base_hash": next(iter(identity))[0],
            "root_table_entries": next(iter(identity))[1]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, point, scalar = frozen_point()
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    if (host["schema"] != "n83.online-k3114-host/v1"
            or host["git_sha"] is None
            or source["schema"] != "n83.online-k3114-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(RUNNER)
            or source["arm_order"] != [name for name, _, _, _ in runner.ARMS]
            or source["rho_seed"] != runner.cold.RHO_SEED
            or source["binary_sha256"] != runner.parent.digest(
                output.parent / "bin" / "koblitz_orbit_dlp_fast")
            or source["rho_binary_sha256"] != runner.parent.digest(
                output.parent / "bin" / "n83_rho_optimization_eval")):
        raise ValueError("host, source, binary or arm order differs")
    ordered = [name for name, _, _, _ in runner.ARMS] + ["rho"]
    finished = []
    for name in ordered:
        path = output / "runs" / name / "usage.json"
        if not path.exists():
            break
        finished.append(datetime.fromisoformat(one_json(path)["finished_utc"]))
    order_status = ("PASS" if len(finished) == len(ordered)
                    and all(a < b for a, b in zip(finished, finished[1:]))
                    else "INCOMPLETE_OR_OUT_OF_ORDER")
    arms = []
    for name, k, layout, _ in runner.ARMS:
        try:
            item = dynamic.verify_ic(output, name, k, layout, point, scalar)
        except Exception as error:
            item = {"arm": name, "k": k, "layout": layout,
                    "status": "VERIFICATION_FAILED", "error": str(error)}
        arms.append(item)
    try:
        rho = verify_rho(output, point, scalar)
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = {item["arm"]: item for item in arms
                if item["status"] == "PASS_INDEPENDENT_REPLAY"
                and item["resources"]["status"] == "PASS"}
    controls = {prefix: online_control(accepted, prefix)
                for prefix in ("k3072", "k3114")}
    noise_ms = None
    separation_ms = None
    winner = None
    if all(item["status"] == "PASS_10_PERCENT" for item in controls.values()):
        noise_ms = max(value for item in controls.values()
                       for value in item["pair_absolute_differences_ms"])
        lower_3114 = min(controls["k3114"]["online_ms"])
        upper_3072 = max(controls["k3072"]["online_ms"])
        lower_3072 = min(controls["k3072"]["online_ms"])
        upper_3114 = max(controls["k3114"]["online_ms"])
        if lower_3114 - upper_3072 > noise_ms:
            winner = "k3072root11"
            separation_ms = lower_3114 - upper_3072
        elif lower_3072 - upper_3114 > noise_ms:
            winner = "k3114root11"
            separation_ms = lower_3072 - upper_3114
    statuses = {"PASS_INDEPENDENT_REPLAY", "TIMEOUT", "STOPPED_RESOURCE_GUARD",
                "FAILED"}
    evidence_valid = (order_status == "PASS"
                      and all(item["status"] in statuses for item in arms)
                      and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                      and all(item["status"] != "BASE_OR_INDEX_MISMATCH"
                              for item in controls.values()))
    bounded = [item["arm"] for item in arms
               if item["status"] != "PASS_INDEPENDENT_REPLAY"
               or item.get("resources", {}).get("status") != "PASS"]
    analysis = {
        "schema": "n83.online-k3114-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "one_exposed_public_point_k3114_capacity_frontier",
        "status": "SCREEN_SELECTED" if evidence_valid and winner else
                  "BOUNDED_UNKNOWN",
        "evidence_status": ("PASS_INDEPENDENT_REPLAY" if evidence_valid
                            else "VERIFICATION_FAILED"),
        "selected_pipeline": winner if evidence_valid else None,
        "selected_candidate_id": None,
        "online_aa_controls": controls,
        "winner_separation_ms": separation_ms,
        "max_adjacent_aa_difference_ms": noise_ms,
        "bounded_arms": bounded,
        "execution_order_status": order_status,
        "arms": arms,
        "rho": rho,
        "controlled_online_speedup": None,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if evidence_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
