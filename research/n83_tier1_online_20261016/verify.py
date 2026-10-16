#!/usr/bin/env python3
"""Independent five-point N83 IC/rho replay and Tier 1 admission."""

import argparse
import importlib.util
import json
import math
import statistics
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNNER = HERE / "run.py"
DYNAMIC_VERIFIER = ROOT / "research/n83_root11_dynamic_20261016/verify.py"
T_CRITICAL_DF4 = 2.7764451051977987


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_tier1_runner", RUNNER)
dynamic = module("n83_tier1_dynamic_verifier", DYNAMIC_VERIFIER)
group = dynamic.group


def one_json(path):
    return json.loads(path.read_text())


def frozen_answers():
    hashes, encoded_points = runner.fixture_preflight()
    answers = [json.loads(line) for line in
               (HERE / "answers.jsonl").read_text().splitlines()]
    if len(answers) != 5:
        raise ValueError("five answer sidecar rows required")
    scalars = []
    for index, (encoded, answer) in enumerate(zip(encoded_points, answers)):
        point = group.point_from_json(encoded)
        scalar = int(answer["scalar"])
        if (answer["fixture_index"] != index or answer["point"] != encoded
                or not group.on_curve(point)
                or group.multiply(point, group.R) is not None
                or group.multiply(group.G, scalar) != point):
            raise ValueError(f"point {index + 1} fails independent subgroup replay")
        scalars.append(scalar)
    return hashes, encoded_points, scalars


def telemetry(output, name, cpus):
    path = output / "runs" / name / "telemetry.json"
    if not path.exists():
        return {"status": "MISSING"}
    receipt = one_json(path)
    steal = receipt.get("steal_fraction_of_wall")
    if (receipt.get("selected_cpus") != cpus
            or not isinstance(steal, (int, float))
            or not math.isfinite(steal) or steal < 0
            or receipt.get("interrupt_delta", -1) < 0):
        return {"status": "INVALID", "receipt": receipt}
    return {"status": "PASS" if steal <= 0.01 else "STEAL_EXCEEDED_1_PERCENT",
            "receipt": receipt}


def verify_rho(output, name, point, scalar, seed, input_path, cpus):
    folder = output / "runs" / name
    if not (folder / "usage.json").exists():
        return {"status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"status": usage["status"], "usage": usage,
                "telemetry": telemetry(output, name, cpus)}
    receipt = one_json(folder / "receipt.json")
    recovered = int(receipt["recovered_scalar"])
    online_ms = float(receipt["online_ms"])
    if (receipt["status"] != "PASS_VERIFIED_PUBLIC_TARGET"
            or group.point_from_json(receipt["public_target_q"]) != point
            or recovered != scalar or group.multiply(group.G, recovered) != point
            or receipt["walk_seed"] != seed or receipt["target_count"] != 1
            or receipt["known_answer_input"] is not False
            or receipt["validation_sidecar_opened"] is not False
            or receipt["group_verified"] is not True
            or receipt["generic_scalar_replay_verified"] is not True
            or receipt["public_input_sha256"] != runner.parent.digest(input_path)
            or receipt["compiled_source_sha256"][
                "examples/n83_rho_optimization_eval.rs"] != runner.fixture_preflight()[
                    0]["examples/n83_rho_optimization_eval.rs"]
            or not math.isfinite(online_ms) or online_ms <= 0):
        raise ValueError(f"{name}: rho source, point or scalar replay failed")
    return {"status": "PASS_INDEPENDENT_REPLAY", "online_ms": online_ms,
            "resources": dynamic.coldver.resource(usage),
            "telemetry": telemetry(output, name, cpus),
            "full_child_wall_seconds": usage["wall_seconds"]}


def checked_ic(output, name, k, layout, point, scalar, cpus):
    try:
        result = dynamic.verify_ic(output, name, k, layout, point, scalar)
        result["telemetry"] = telemetry(output, name, cpus)
        return result
    except Exception as error:
        return {"arm": name, "k": k, "layout": layout,
                "status": "VERIFICATION_FAILED", "error": str(error)}


def checked_rho(output, name, point, scalar, seed, input_path, cpus):
    try:
        return verify_rho(output, name, point, scalar, seed, input_path, cpus)
    except Exception as error:
        return {"status": "VERIFICATION_FAILED", "error": str(error)}


def admitted(row):
    return (row["status"] == "PASS_INDEPENDENT_REPLAY"
            and row["resources"]["status"] == "PASS"
            and row["telemetry"]["status"] == "PASS")


def paired_interval(ratios):
    logs = [math.log(value) for value in ratios]
    mean = statistics.mean(logs)
    halfwidth = T_CRITICAL_DF4 * statistics.stdev(logs) / math.sqrt(len(logs))
    return {"mean_log_ratio": mean, "sample_sd_log_ratio": statistics.stdev(logs),
            "lower_95_log_ratio": mean - halfwidth,
            "upper_95_log_ratio": mean + halfwidth,
            "lower_95_ratio": math.exp(mean - halfwidth),
            "upper_95_ratio": math.exp(mean + halfwidth)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    start = time.perf_counter()
    hashes, encoded_points, scalars = frozen_answers()
    selection, candidate_path = runner.selection_preflight()
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    cpus = host["cpu_affinity"]
    identity = host["instance_identity"]
    provider = host["provider"]
    aws_host = (provider == "aws_ec2_spot"
                and identity.get("availabilityZone") == runner.EXPECTED_AZ
                and identity.get("instanceType") == runner.EXPECTED_TYPE)
    github_host = (provider == "github_actions"
                   and identity.get("repository") == "aburan28/crypto"
                   and identity.get("runnerOS") == "Linux"
                   and identity.get("sha") == host["git_sha"]
                   and isinstance(identity.get("runId"), int)
                   and identity["runId"] > 0
                   and isinstance(identity.get("runAttempt"), int)
                   and identity["runAttempt"] > 0
                   and identity.get("instanceId") == host["instance_id"]
                   and host.get("run_url") ==
                   f"https://github.com/aburan28/crypto/actions/runs/{identity['runId']}")
    if github_host:
        seal = runner.github_seal_preflight(selection, identity)
        github_host = (seal["selector_cpus"] == cpus
                       and host["visible_affinity"] == cpus)
    schedule = [{"point_index": i, "pair_order": ["ic", "rho"] if i % 2 else
                 ["rho", "ic"]} for i in range(1, 6)]
    if (host["schema"] != "n83.tier1-online-host/v1"
            or not (aws_host or github_host)
            or len(cpus) != 4 or len(set(cpus)) != 4
            or source["schema"] != "n83.tier1-online-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(RUNNER)
            or source["selection_sha256"] != runner.parent.digest(
                HERE / "selection.json")
            or source["candidate_manifest_sha256"] != runner.parent.digest(
                candidate_path)
            or source["candidate_id"] != selection["candidate_id"]
            or source["binary_sha256"] != runner.parent.digest(
                output.parent / "bin" / "koblitz_orbit_dlp_fast")
            or source["rho_binary_sha256"] != runner.parent.digest(
                output.parent / "bin" / "n83_rho_optimization_eval")
            or source["schedule"] != schedule
            or source["rho_seed_base"] != runner.RHO_SEED_BASE
            or source["calibration_seed"] != runner.CALIBRATION_SEED):
        raise ValueError("cloud allocation, selected candidate, source or binary differs")
    point_inputs = []
    for index, encoded in enumerate(encoded_points, start=1):
        path = output / "inputs" / f"point_{index}.jsonl"
        lines = path.read_text().splitlines()
        if (len(lines) != 1 or json.loads(lines[0]) != encoded
                or source["point_input_sha256"][str(index)] !=
                runner.parent.digest(path)):
            raise ValueError(f"point {index}: public producer input differs")
        point_inputs.append(path)
    order = []
    for row in schedule:
        index = row["point_index"]
        order.extend([f"p{index}_cal_before"] +
                     [f"p{index}_{arm}" for arm in row["pair_order"]] +
                     [f"p{index}_cal_after"])
    finished = []
    for name in order:
        path = output / "runs" / name / "usage.json"
        if not path.exists():
            break
        finished.append(datetime.fromisoformat(one_json(path)["finished_utc"]))
    order_status = ("PASS" if len(finished) == len(order)
                    and all(a < b for a, b in zip(finished, finished[1:]))
                    else "INCOMPLETE_OR_OUT_OF_ORDER")
    calibration_point = group.point_from_json(one_json(runner.cold.TARGETS))
    calibration_scalar = int(one_json(runner.cold.ANSWERS)["scalar"])
    point_rows = []
    base_identity = set()
    for index, (encoded, scalar) in enumerate(zip(encoded_points, scalars), start=1):
        point = group.point_from_json(encoded)
        checks = {}
        for suffix in ("before", "after"):
            name = f"p{index}_cal_{suffix}"
            checks[suffix] = checked_rho(
                output, name, calibration_point, calibration_scalar,
                int(runner.CALIBRATION_SEED), runner.cold.TARGETS, cpus)
        ic = checked_ic(output, f"p{index}_ic", selection["k"],
                        selection["layout"], point, scalar, cpus)
        rho = checked_rho(output, f"p{index}_rho", point, scalar,
                          runner.RHO_SEED_BASE + index,
                          point_inputs[index - 1], cpus)
        if ic["status"] == "PASS_INDEPENDENT_REPLAY":
            base_identity.add((ic["base_hash"], ic["root_table_entries"]))
        aa = None
        if admitted(checks["before"]) and admitted(checks["after"]):
            times = (checks["before"]["online_ms"],
                     checks["after"]["online_ms"])
            spread = abs(times[0] - times[1]) / statistics.mean(times)
            aa = {"online_ms": list(times), "relative_difference": spread,
                  "status": "PASS_10_PERCENT" if spread <= 0.10 else "NOISY"}
        pair_admitted = (admitted(ic) and admitted(rho) and aa is not None
                         and aa["status"] == "PASS_10_PERCENT")
        ratio = rho["online_ms"] / ic["one_target_online_ms"] if pair_admitted else None
        point_rows.append({"point_index": index, "public_point": encoded,
                           "calibration": checks, "calibration_aa": aa,
                           "ic": ic, "rho": rho,
                           "pair_admitted": pair_admitted, "online_ratio": ratio})
    evidence_valid = (order_status == "PASS"
                      and len(base_identity) <= 1
                      and all(row["ic"]["status"] != "VERIFICATION_FAILED"
                              and row["rho"]["status"] != "VERIFICATION_FAILED"
                              and all(control["status"] != "VERIFICATION_FAILED"
                                      for control in row["calibration"].values())
                              for row in point_rows))
    full_panel = evidence_valid and all(row["pair_admitted"] for row in point_rows)
    ratios = [row["online_ratio"] for row in point_rows] if full_panel else []
    interval = paired_interval(ratios) if full_panel else None
    win = (full_panel and all(ratio > 1 for ratio in ratios)
           and interval["lower_95_log_ratio"] > 0)
    status = ("TIER1_ONLINE_SPEEDUP_PASS" if win else
              "VALID_PAIRS_NO_SPEEDUP" if full_panel else "BOUNDED_UNKNOWN")
    analysis = {
        "schema": "n83.tier1-online-analysis/v1",
        "hardware_tier": ("tier1_github_hosted_runner" if github_host
                          else "tier1_aws_ec2_spot"),
        "instance_id": host["instance_id"],
        "candidate_id": selection["candidate_id"],
        "workload_role": "five_fresh_single_target_online_pairs",
        "status": status,
        "evidence_status": ("PASS_INDEPENDENT_REPLAY" if evidence_valid
                            else "VERIFICATION_FAILED"),
        "execution_order_status": order_status,
        "cross_point_base_identity": [list(pair) for pair in sorted(base_identity)],
        "pair_count": len(point_rows), "admitted_pair_count": sum(
            row["pair_admitted"] for row in point_rows),
        "online_ratios": ratios,
        "median_online_ratio": statistics.median(ratios) if ratios else None,
        "mean_ic_online_ms": statistics.mean(
            row["ic"]["one_target_online_ms"] for row in point_rows)
            if full_panel else None,
        "mean_rho_online_ms": statistics.mean(
            row["rho"]["online_ms"] for row in point_rows)
            if full_panel else None,
        "paired_log_ratio_interval": interval,
        "point_rows": point_rows,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    (output / "replay_timing.json").write_text(json.dumps({
        "independent_replay_wall_ms": (time.perf_counter() - start) * 1000,
        "outside_producer_online_intervals": True}, indent=2,
        sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if evidence_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
