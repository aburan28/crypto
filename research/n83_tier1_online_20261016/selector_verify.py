#!/usr/bin/env python3
"""Independent N83 mirrored finalist A/A, relation and scalar replay."""

import argparse
import importlib.util
import json
import math
import statistics
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNNER = HERE / "selector_run.py"
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


def telemetry_status(output, arm, cpus):
    telemetry = one_json(output / "runs" / arm / "telemetry.json")
    steal = telemetry["steal_fraction_of_wall"]
    status = ("PASS" if telemetry["selected_cpus"] == cpus
              and isinstance(steal, (int, float)) and math.isfinite(steal)
              and 0 <= steal <= 0.01 else "REJECTED_STEAL_OR_AFFINITY")
    return status, telemetry


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
    pair_differences = [abs(times[0] - times[1])]
    pair_spreads = [pair_differences[0] / statistics.mean(times)]
    return {"status": ("PASS_10_PERCENT" if all(value <= 0.10 for value
                       in pair_spreads) else "NOISY"),
            "arms": names, "online_ms": times,
            "mean_online_ms": statistics.mean(times),
            "pair_relative_differences": pair_spreads,
            "pair_absolute_differences_ms": pair_differences,
            "base_hash": next(iter(identity))[0],
            "root_table_entries": next(iter(identity))[1]}


def unsupported_cells():
    cells = []
    working_reserve_bytes = 1 << 30
    for k in (3450, 4096):
        states, request, slots, _, _ = dynamic.expected_layout(
            k, "root11_dynamic")
        core_bytes = 11 * states + 5 * slots
        if core_bytes + working_reserve_bytes <= runner.cold.RSS_CAP_BYTES:
            raise ValueError(f"K{k} dynamic capacity exclusion no longer applies")
        if (k == 3450 and not states < (1 << 29)) or (
                k == 4096 and not states >= (1 << 29)):
            raise ValueError(f"K{k} packed state-ID boundary differs")
        cells.append({"k": k, "layout": "root11_dynamic",
                      "state_count": states, "state_id_limit": 1 << 29,
                      "request_entries": request, "allocated_slots": slots,
                      "root_state_bytes": 11 * states,
                      "root_table_bytes": 5 * slots,
                      "core_bytes_before_base": core_bytes,
                      "working_reserve_bytes": working_reserve_bytes,
                      "rss_cap_bytes": runner.cold.RSS_CAP_BYTES,
                      "status": ("UNSUPPORTED_STATE_ID_AND_CAPACITY" if k == 4096
                                 else "UNSUPPORTED_CAPACITY")})
    return cells


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, point, scalar = frozen_point()
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    envelope = source.get("cpu_envelope")
    cpu_count = {"github_two_core": 2, "aws_four_core": 4}.get(envelope)
    if (host["schema"] != "n83.tier1-selector-host/v1"
            or host["git_sha"] is None
            or cpu_count is None
            or len(host["cpu_affinity"]) != cpu_count
            or host["visible_affinity"] != host["cpu_affinity"]
            or host.get("cpu_envelope") != envelope
            or source["schema"] != "n83.tier1-selector-source/v1"
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
            if item["status"] == "PASS_INDEPENDENT_REPLAY":
                status, telemetry = telemetry_status(
                    output, name, host["cpu_affinity"])
                item["telemetry"] = telemetry
                if status != "PASS":
                    item["resources"]["status"] = status
        except Exception as error:
            item = {"arm": name, "k": k, "layout": layout,
                    "status": "VERIFICATION_FAILED", "error": str(error)}
        arms.append(item)
    try:
        rho = verify_rho(output, point, scalar)
        if rho["status"] == "PASS_INDEPENDENT_REPLAY":
            status, telemetry = telemetry_status(
                output, "rho", host["cpu_affinity"])
            rho["telemetry"] = telemetry
            if status != "PASS":
                rho["status"] = status
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = {item["arm"]: item for item in arms
                if item["status"] == "PASS_INDEPENDENT_REPLAY"
                and item["resources"]["status"] == "PASS"}
    controls = {prefix: online_control(accepted, prefix)
                for prefix in ("k3072", "k2543")}
    noise_ms = None
    separation_ms = None
    winner = None
    if all(item["status"] == "PASS_10_PERCENT" for item in controls.values()):
        noise_ms = max(value for item in controls.values()
                       for value in item["pair_absolute_differences_ms"])
        lower_2543 = min(controls["k2543"]["online_ms"])
        upper_3072 = max(controls["k3072"]["online_ms"])
        lower_3072 = min(controls["k3072"]["online_ms"])
        upper_2543 = max(controls["k2543"]["online_ms"])
        if lower_2543 - upper_3072 > noise_ms:
            winner = "k3072root11"
            separation_ms = lower_2543 - upper_3072
        elif lower_3072 - upper_2543 > noise_ms:
            winner = "k2543root11"
            separation_ms = lower_3072 - upper_2543
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
        "schema": "n83.tier1-selector-analysis/v1",
        "hardware_tier": "tier1_cloud_instance_selector",
        "workload_role": "one_exposed_public_point_same_instance_selector",
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
