#!/usr/bin/env python3
"""Independently replay one N83 finalist Tier 1 host stratum."""

import argparse
import hashlib
import importlib.util
import json
import math
import statistics
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DOMAIN = b"n83-e2e-final-holdout-sixteen-20261016-v1"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_finalist_runner", HERE / "run.py")
finegrid = module("n83_finalist_finegrid_identity",
                  ROOT / "research/n83_e2e_finegrid_20261016/materialize.py")
prior = module("n83_finalist_holdout_verify",
               ROOT / "research/n83_e2e_holdout_20261016/verify.py")
dynamic = prior.dynamic
group = dynamic.group


def read(path):
    return json.loads(path.read_text())


def suffix(path, count):
    return Path(path).parts[-count:]


def fixture_points(fixture_dir):
    runner.fixture(fixture_dir)
    encoded = [json.loads(line) for line in
               (fixture_dir / "targets.jsonl").read_text().splitlines()]
    answers = [json.loads(line) for line in
               (fixture_dir / "answers.jsonl").read_text().splitlines()]
    if (len(encoded) != 16 or len(answers) != 16
            or len({tuple(item) for item in encoded}) != 16):
        raise ValueError("sealed finalist point set differs")
    points, scalars = [], []
    for index, (item, answer) in enumerate(zip(encoded, answers)):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.point_from_json(item)
        if (answer["fixture_index"] != index or answer["point"] != item
                or answer["seed_sha256"] != seed.hex()
                or int(answer["scalar"]) != scalar
                or not group.on_curve(point)
                or group.multiply(point, group.R) is not None
                or group.multiply(group.G, scalar) != point):
            raise ValueError(f"sealed finalist point {index} fails replay")
        points.append(point)
        scalars.append(scalar)
    return encoded, points, scalars


def arm_folder(output, event):
    index, arm = event["point_index"], event["name"]
    point_output = output / "points" / f"p{index:02d}"
    return point_output, arm, point_output / "runs" / arm


def replay_ic(output, event, point, scalar, by_name, field, curve):
    point_output, arm, folder = arm_folder(output, event)
    name = event["candidate"]
    _, k, layout, expected_id = by_name[name]
    try:
        result = dynamic.verify_ic(point_output, arm, k, layout, point, scalar)
        result["candidate_id"] = expected_id
        if result["status"] == "PASS_INDEPENDENT_REPLAY":
            summary = read(folder / "stdout.jsonl")
            manifest = finegrid.dynamic.candidate(
                k, layout, summary, folder / "base_dump.jsonl", field, curve)
            if manifest["candidate_id"] != expected_id:
                raise ValueError("exact IC1 candidate identity differs")
        return result
    except Exception as error:
        return {"arm": arm, "k": k, "layout": layout,
                "candidate_id": expected_id,
                "status": "VERIFICATION_FAILED", "error": str(error)}


def replay_rho(output, event, point, scalar, seed, target, source_hashes):
    _, _, folder = arm_folder(output, event)
    try:
        if not (folder / "usage.json").exists():
            return {"status": "NOT_RUN"}
        usage = read(folder / "usage.json")
        if usage["status"] != "PASS_PRODUCER":
            return {"status": usage["status"], "usage": usage}
        receipt = read(folder / "receipt.json")
        recovered = int(receipt["recovered_scalar"])
        online_ms = float(receipt["online_ms"])
        if (receipt["status"] != "PASS_VERIFIED_PUBLIC_TARGET"
                or group.point_from_json(receipt["public_target_q"]) != point
                or recovered != scalar or receipt["walk_seed"] != seed
                or receipt["target_count"] != 1
                or receipt["known_answer_input"] is not False
                or receipt["validation_sidecar_opened"] is not False
                or receipt["group_verified"] is not True
                or receipt["generic_scalar_replay_verified"] is not True
                or receipt["public_input_sha256"] != runner.parent.digest(target)
                or receipt["compiled_source_sha256"]["examples/n83_rho_optimization_eval.rs"] !=
                source_hashes["examples/n83_rho_optimization_eval.rs"]
                or group.multiply(group.G, recovered) != point
                or not math.isfinite(online_ms) or online_ms <= 0):
            raise ValueError("same-point rho or calibration scalar replay failed")
        return {"status": "PASS_INDEPENDENT_REPLAY", "online_ms": online_ms,
                "resources": dynamic.coldver.resource(usage),
                "full_child_wall_seconds": usage["wall_seconds"]}
    except Exception as error:
        return {"status": "VERIFICATION_FAILED", "error": str(error)}


def control(primary, repeat):
    if primary is None or repeat is None:
        return {"status": "INCOMPLETE"}
    if (primary.get("status") != "PASS_INDEPENDENT_REPLAY"
            or repeat.get("status") != "PASS_INDEPENDENT_REPLAY"):
        return {"status": "REPLAY_FAILED"}
    if (primary["candidate_id"] != repeat["candidate_id"]
            or primary["base_hash"] != repeat["base_hash"]
            or primary["root_table_entries"] != repeat["root_table_entries"]):
        return {"status": "BASE_OR_INDEX_MISMATCH"}
    times = (primary["cold_process_ms_exploratory"],
             repeat["cold_process_ms_exploratory"])
    spread = abs(times[0] - times[1]) / statistics.mean(times)
    return {"status": "PASS_5_PERCENT" if spread <= 0.05
            else "NOISY_5_PERCENT", "cold_process_ms": times,
            "relative_spread": spread}


def calibration(before, after):
    if before is None or after is None:
        return {"status": "INCOMPLETE"}
    if any(item.get("status") != "PASS_INDEPENDENT_REPLAY" or
           item.get("resources", {}).get("status") != "PASS"
           for item in (before, after)):
        return {"status": "REPLAY_OR_RESOURCE_FAILED"}
    times = (before["online_ms"], after["online_ms"])
    spread = abs(times[0] - times[1]) / statistics.mean(times)
    return {"status": "PASS_10_PERCENT" if spread <= 0.10
            else "NOISY_10_PERCENT", "online_ms": times,
            "relative_spread": spread}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--half", type=int, choices=(0, 1), required=True)
    parser.add_argument("--fixture-dir", type=Path, default=HERE)
    args = parser.parse_args()
    output, fixture_dir = args.output.resolve(), args.fixture_dir.resolve()
    hashes, lines = runner.fixture(fixture_dir)
    _, points, scalars = fixture_points(fixture_dir)
    host = read(output / "host.json")
    source = read(output / "source.json")
    preflight = read(output / "preflight.json")
    events = runner.schedule(args.half)
    if (host["schema"] != "n83.e2e-finalist-host/v1"
            or host["git_sha"] is None
            or host["stratum_half"] != args.half
            or host["instance_id"] !=
            f"{host['parent_instance_id']}-half{args.half}"
            or preflight["status"] != "PASS_CAPACITY_PREFLIGHT"
            or source["schema"] != "n83.e2e-finalist-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(HERE / "run.py")
            or source["targets_sha256"] != runner.TARGETS_SHA256
            or source["answers_sha256"] != runner.ANSWERS_SHA256
            or source["fixture_manifest_sha256"] != runner.parent.digest(
                fixture_dir / "fixture_manifest.json")
            or source["half"] != args.half
            or source["arms"] != [list(arm) for arm in runner.ARMS]
            or source["events"] != events
            or source["rho_seed_rule"] != "2026101700+8*point_index+fixed_candidate_index"
            or source["calibration_seed"] != runner.CALIBRATION_RHO_SEED
            or source["calibration_target_sha256"] != runner.parent.digest(runner.cold.TARGETS)
            or source["calibration_answer_sha256"] != runner.parent.digest(runner.cold.ANSWERS)
            or source["global_producer_cap_seconds"] !=
            runner.GLOBAL_PRODUCER_CAP_SECONDS):
        raise ValueError("finalist source, host, fixture or schedule differs")
    binaries_match = (
        runner.parent.digest(output / "executed_bins/koblitz_orbit_dlp_fast") ==
        source["binary_sha256"] and
        runner.parent.digest(output / "executed_bins/n83_rho_optimization_eval") ==
        source["rho_binary_sha256"])
    topology = host.get("topology", {})
    topology_status = ("PASS" if topology.get("lscpu_json")
                       and len(topology.get("visible_cpu_topology", {})) >= 4
                       and (output / "host_telemetry_after.json").exists()
                       else "MISSING")
    finished = []
    telemetry = []
    source_targets = []
    launch_checks = []
    required_events = [event for event in events if event["kind"] != "control"]
    for event in events:
        point_output, arm, folder = arm_folder(output, event)
        usage_path = folder / "usage.json"
        if usage_path.exists():
            finished.append(datetime.fromisoformat(read(usage_path)["finished_utc"]))
        telemetry.append(prior.telemetry_receipt(point_output, arm))
        target = point_output / "target.jsonl"
        source_targets.append(target.exists() and
                              target.read_text() == lines[event["point_index"]])
        try:
            launch = read(folder / "launch.json")
            command = launch["command"]
            if event["kind"] in ("rho", "calibration"):
                expected_target = (runner.cold.TARGETS if event["kind"] == "calibration"
                                   else target)
                seed = (runner.CALIBRATION_RHO_SEED if event["kind"] == "calibration"
                        else runner.PAIR_RHO_SEED_BASE + 8 * event["point_index"] +
                        next(i for i, arm in enumerate(runner.ARMS)
                             if arm[0] == event["candidate"]))
                valid = (suffix(command[1], 3) == suffix(expected_target, 3)
                         and command[2] == str(seed)
                         and suffix(command[3], 3) ==
                         suffix(runner.cold.RHO_CONFIG, 3)
                         and suffix(command[4], 5) ==
                         suffix(folder / "receipt.json", 5)
                         and launch["environment_flags"] == {}
                         and launch["cap_seconds"] == runner.RHO_CAP_SECONDS
                         and Path(command[0]).parts[-2:] ==
                         ("executed_bins", "n83_rho_optimization_eval"))
            else:
                _, k, layout, _ = {item[0]: item for item in runner.ARMS}[
                    event["candidate"]]
                expected_flags = runner.dynamic.environment(folder, layout)
                expected_flags = {key: value for key, value in
                                  expected_flags.items() if key.startswith("KIC_")}
                observed_flags = dict(launch["environment_flags"])
                observed_dump = observed_flags.pop("KIC_DUMP_BASE", "")
                expected_dump = expected_flags.pop("KIC_DUMP_BASE")
                valid = (command[1] == f"construct:83:1:{k}"
                         and suffix(command[2], 3) == suffix(target, 3)
                         and command[3] == "29"
                         and suffix(command[4], 5) ==
                         suffix(folder / "targets_result.jsonl", 5)
                         and observed_flags == expected_flags
                         and suffix(observed_dump, 5) ==
                         suffix(expected_dump, 5)
                         and launch["cap_seconds"] == runner.IC_CAP_SECONDS
                         and Path(command[0]).parts[-2:] ==
                         ("executed_bins", "koblitz_orbit_dlp_fast"))
            launch_checks.append(valid)
        except (KeyError, IndexError, OSError, ValueError, TypeError):
            launch_checks.append(False)
    primary_finished = [datetime.fromisoformat(read(
        arm_folder(output, event)[2] / "usage.json")["finished_utc"])
        for event in required_events if
        (arm_folder(output, event)[2] / "usage.json").exists()]
    order_status = ("PASS" if len(primary_finished) == len(required_events)
                    and all(a < b for a, b in
                            zip(primary_finished, primary_finished[1:]))
                    else "INCOMPLETE_OR_OUT_OF_ORDER")
    primary_indices = [i for i, event in enumerate(events)
                       if event["kind"] != "control"]
    telemetry_status = ("PASS" if topology_status == "PASS"
                        and all(telemetry[i]["status"] == "PASS"
                                and telemetry[i]["steal_jiffies_delta"] == 0
                                for i in primary_indices)
                        else "INCOMPLETE_OR_STEAL")
    field, curve = finegrid.panel.curve_record()
    by_name = {arm[0]: arm for arm in runner.ARMS}
    point_rows = {index: {"point_index": index, "arms": [],
                          "rho_by_candidate": {}, "calibration_by_candidate": {}}
                  for index in range(args.half * 8, args.half * 8 + 8)}
    controls = {}
    cal_encoded = json.loads(runner.cold.TARGETS.read_text().splitlines()[0])
    cal_answer = json.loads(runner.cold.ANSWERS.read_text().splitlines()[0])
    cal_point = group.point_from_json(cal_encoded)
    cal_scalar = int(cal_answer["scalar"])
    if (cal_answer["fixture_index"] != 0 or cal_answer["point"] != cal_encoded
            or group.multiply(group.G, cal_scalar) != cal_point
            or group.multiply(cal_point, group.R) is not None):
        raise ValueError("calibration fixture scalar replay failed")
    calibrations = {}
    for event in events:
        index, name, candidate, kind = (event[key] for key in
                                         ("point_index", "name", "candidate", "kind"))
        if kind == "rho":
            seed = runner.PAIR_RHO_SEED_BASE + 8 * index + next(
                i for i, arm in enumerate(runner.ARMS) if arm[0] == candidate)
            point_rows[index]["rho_by_candidate"][by_name[candidate][3]] = replay_rho(
                output, event, points[index], scalars[index], seed,
                output / "points" / f"p{index:02d}" / "target.jsonl", hashes)
        elif kind == "calibration":
            calibrations[(index, candidate, name.startswith("cal_before_"))] = replay_rho(
                output, event, cal_point, cal_scalar, runner.CALIBRATION_RHO_SEED,
                runner.cold.TARGETS, hashes)
        else:
            result = replay_ic(output, event, points[index], scalars[index],
                               by_name, field, curve)
            if kind == "control":
                controls[candidate] = result
            else:
                point_rows[index]["arms"].append(result)
    for index, point in point_rows.items():
        for name, _, _, candidate_id in runner.ARMS:
            point["calibration_by_candidate"][candidate_id] = calibration(
                calibrations.get((index, name, True)),
                calibrations.get((index, name, False)))
    primary = {item["arm"]: item for item in point_rows[args.half * 8]["arms"]}
    aa = {name: control(primary.get(name), controls.get(name))
          for name in runner.CONTROL_KEYS}
    base_status = "PASS"
    for k in (1104, 1107, 1112, 1120):
        hashes_for_k = {item.get("base_hash") for point in point_rows.values()
                        for item in point["arms"] if item["k"] == k
                        and item["status"] == "PASS_INDEPENDENT_REPLAY"}
        if len(hashes_for_k) != 1:
            base_status = "MISMATCH"
    arms = [item for point in point_rows.values() for item in point["arms"]]
    rho = [item for point in point_rows.values() for item in
           point["rho_by_candidate"].values()]
    calibration_rows = list(calibrations.values())
    calibration_status = ("PASS_10_PERCENT_PER_PAIR" if
                          len(calibration_rows) == 128 and all(
                              len(point["calibration_by_candidate"]) == 8 and
                              all(item["status"] == "PASS_10_PERCENT" for item in
                                  point["calibration_by_candidate"].values())
                              for point in point_rows.values()) else
                          "NOISY_OR_INCOMPLETE")
    replay_status = ("PASS_INDEPENDENT_REPLAY" if binaries_match
                     and all(source_targets[i] for i in primary_indices)
                     and all(launch_checks[i] for i in primary_indices)
                     and len(arms) == 64 and len(rho) == 64
                     and len(calibration_rows) == 128
                     and all(item["status"] == "PASS_INDEPENDENT_REPLAY"
                             for item in arms + rho + calibration_rows)
                     and base_status == "PASS"
                     else "VERIFICATION_FAILED")
    resource_status = ("PASS" if replay_status == "PASS_INDEPENDENT_REPLAY"
                       and all(item["resources"]["status"] == "PASS"
                               for item in arms + rho + calibration_rows)
                       else "FAILED_OR_INCOMPLETE")
    admitted = (replay_status == "PASS_INDEPENDENT_REPLAY"
                and resource_status == "PASS"
                and order_status == "PASS" and telemetry_status == "PASS"
                and not (output / "producer_stop.json").exists()
                and calibration_status == "PASS_10_PERCENT_PER_PAIR")
    validated = []
    for index in sorted(point_rows):
        point = point_rows[index]
        converted = []
        for item in point["arms"]:
            converted.append({
                "candidate_id": item.get("candidate_id"),
                "status": item["status"],
                "resource_status": item.get("resources", {}).get("status"),
                "one_target_online_ms": item.get("one_target_online_ms"),
                "cold_process_ms": item.get("cold_process_ms_exploratory"),
                "raw_replay": item,
            })
        validated.append({"point_index": index,
                          "rho_by_candidate": {
                              candidate_id: {"status": item["status"],
                                             "resource_status": item.get("resources", {}).get("status"),
                                             "online_ms": item.get("online_ms"),
                                             "raw_replay": item}
                              for candidate_id, item in point["rho_by_candidate"].items()},
                          "calibration_by_candidate": point["calibration_by_candidate"],
                          "arms": converted})
    analysis = {
        "schema": "n83.e2e-finalist-half-analysis/v1", "half": args.half,
        "point_indices": list(sorted(point_rows)),
        "source_sha256": runner.SOURCE_SHA256,
        "targets_sha256": runner.TARGETS_SHA256,
        "answers_sha256": runner.ANSWERS_SHA256,
        "hardware_tier": "cloud_instance_tier1" if admitted else
                         "cloud_instance_exploratory",
        "status": "HALF_ADMITTED" if admitted else "BOUNDED_UNKNOWN",
        "evidence_status": replay_status,
        "tier_status": "PASS_TIER1" if admitted else "NOT_ADMITTED",
        "aa_status": calibration_status,
        "aa_controls": aa, "control_arms": controls,
        "raw_calibrations": {f"p{index:02d}/{name}/{side}": item for
                             (index, name, side), item in calibrations.items()},
        "base_status": base_status,
        "resource_status": resource_status,
        "execution_order_status": order_status,
        "topology_status": topology_status,
        "telemetry_status": telemetry_status,
        "host_telemetry_by_event": telemetry,
        "producer_stop": (read(output / "producer_stop.json") if
                          (output / "producer_stop.json").exists() else None),
        "all_per_point_targets_match": all(source_targets),
        "all_launch_commands_match": all(launch_checks),
        "archived_executed_binaries_match": binaries_match,
        "points": validated,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if admitted else 2


if __name__ == "__main__":
    raise SystemExit(main())
