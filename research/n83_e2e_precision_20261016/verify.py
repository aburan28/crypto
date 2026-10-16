#!/usr/bin/env python3
"""Independently replay the sealed N83 K1104/K1120 precision panel."""

import argparse
import importlib.util
import json
import math
import statistics
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T_CRITICAL_DF15 = 2.131449545559323


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_precision_run", HERE / "run.py")
previous = module("n83_precision_prior_verify",
                  ROOT / "research/n83_e2e_holdout_20261016/verify.py")
dynamic = previous.dynamic
group = dynamic.group
one_json = previous.one_json
telemetry_receipt = previous.telemetry_receipt


def points_and_answers():
    hashes, manifest = runner.fixture()
    points = [json.loads(line) for line in runner.TARGETS.read_text().splitlines()]
    answers = [json.loads(line) for line in runner.ANSWERS.read_text().splitlines()]
    if (len(points) != 16 or len(answers) != 16
            or len({tuple(point) for point in points}) != 16):
        raise ValueError("sixteen distinct frozen precision points required")
    scalars = []
    for i, (encoded, answer) in enumerate(zip(points, answers)):
        point = group.point_from_json(encoded)
        scalar = int(answer["scalar"])
        if (answer["fixture_index"] != i or answer["point"] != encoded
                or not group.on_curve(point)
                or group.multiply(point, group.R) is not None
                or group.multiply(group.G, scalar) != point):
            raise ValueError(f"precision point {i} fails scalar replay")
        scalars.append(scalar)
    return hashes, manifest, points, scalars


def verify_rho(output, i, point, scalar, hashes):
    folder = output / "runs" / f"p{i}_rho"
    if not (folder / "usage.json").exists():
        return {"status": "NOT_RUN"}
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"status": usage["status"], "usage": usage}
    receipt = one_json(folder / "receipt.json")
    recovered = int(receipt["recovered_scalar"])
    online_ms = float(receipt["online_ms"])
    if (receipt["status"] != "PASS_VERIFIED_PUBLIC_TARGET"
            or group.point_from_json(receipt["public_target_q"]) != point
            or recovered != scalar
            or receipt["walk_seed"] != runner.RHO_SEED_BASE + i
            or receipt["target_count"] != 1
            or receipt["known_answer_input"] is not False
            or receipt["validation_sidecar_opened"] is not False
            or receipt["group_verified"] is not True
            or receipt["generic_scalar_replay_verified"] is not True
            or receipt["public_input_sha256"] != runner.parent.digest(
                output / "inputs" / f"point_{i}.jsonl")
            or receipt["compiled_source_sha256"][
                "examples/n83_rho_optimization_eval.rs"] !=
            hashes["examples/n83_rho_optimization_eval.rs"]
            or group.multiply(group.G, recovered) != point
            or not math.isfinite(online_ms) or online_ms <= 0):
        raise ValueError(f"point {i}: same-point rho replay failed")
    return {"status": "PASS_INDEPENDENT_REPLAY", "online_ms": online_ms,
            "resources": dynamic.coldver.resource(usage),
            "full_child_wall_seconds": usage["wall_seconds"]}


def interval(differences):
    mean = statistics.mean(differences)
    sd = statistics.stdev(differences)
    halfwidth = T_CRITICAL_DF15 * sd / math.sqrt(16)
    return {"differences_ms": differences, "mean_ms": mean,
            "sample_sd_ms": sd, "lower_95_ms": mean - halfwidth,
            "upper_95_ms": mean + halfwidth}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, manifest, encoded_points, scalars = points_and_answers()
    schedule = runner.schedule()
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    preflight = one_json(output / "preflight.json")
    if (host["schema"] != "n83.e2e-precision-host/v1"
            or host["git_sha"] is None
            or preflight["status"] != "PASS_CAPACITY_PREFLIGHT"
            or source["schema"] != "n83.e2e-precision-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.screen.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(HERE / "run.py")
            or source["generator_sha256"] != runner.parent.digest(
                HERE / "generate.py")
            or source["fixture_manifest_sha256"] != runner.parent.digest(
                runner.MANIFEST)
            or source["targets_sha256"] != manifest["targets_sha256"]
            or source["answers_sha256"] != manifest["answers_sha256"]
            or source["k_values"] != list(runner.K_VALUES)
            or source["schedule"] != schedule
            or source["rho_seed_base"] != runner.RHO_SEED_BASE):
        raise ValueError("precision host, source, fixture or schedule differs")
    for i, encoded in enumerate(encoded_points, start=1):
        path = output / "inputs" / f"point_{i}.jsonl"
        if (path.read_text().splitlines() !=
                [json.dumps(encoded, separators=(",", ":"))]
                or source["point_input_sha256"][str(i)] !=
                runner.parent.digest(path)):
            raise ValueError(f"point {i}: executed one-point input differs")
    ordered_events = []
    for row in schedule:
        i = row["point_index"]
        rho = f"p{i}_rho"
        ic = [f"p{i}_{arm}" for arm in row["ic_order"]]
        ordered_events.extend(([rho] if row["rho_first"] else []) + ic
                              + ([] if row["rho_first"] else [rho]))
    topology = host.get("topology", {})
    topology_status = ("PASS" if topology.get("lscpu_json")
                       and len(topology.get("visible_cpu_topology", {})) >= 4
                       and (output / "host_telemetry_after.json").exists()
                       else "MISSING")
    telemetry = {name: telemetry_receipt(output, name)
                 for name in ordered_events}
    telemetry_status = ("PASS" if topology_status == "PASS"
                        and all(item["status"] == "PASS"
                                for item in telemetry.values())
                        else "INCOMPLETE")
    finished = []
    for name in ordered_events:
        path = output / "runs" / name / "usage.json"
        if not path.exists():
            break
        finished.append(datetime.fromisoformat(one_json(path)["finished_utc"]))
    order_status = ("PASS" if len(finished) == len(ordered_events)
                    and all(a < b for a, b in zip(finished, finished[1:]))
                    else "INCOMPLETE_OR_OUT_OF_ORDER")
    point_rows = []
    for row in schedule:
        i = row["point_index"]
        point = group.point_from_json(encoded_points[i - 1])
        scalar = scalars[i - 1]
        arms = []
        for arm in row["ic_order"]:
            k = int(arm.split("_", 1)[0][1:])
            name = f"p{i}_{arm}"
            try:
                item = dynamic.verify_ic(output, name, k, "root11_dynamic",
                                         point, scalar)
                item["pipeline"] = f"k{k}"
                arms.append(item)
            except Exception as error:
                arms.append({"arm": name, "pipeline": f"k{k}", "k": k,
                             "status": "VERIFICATION_FAILED", "error": str(error)})
        try:
            rho = verify_rho(output, i, point, scalar, hashes)
        except Exception as error:
            rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
        accepted = {item["arm"]: item for item in arms
                    if item["status"] == "PASS_INDEPENDENT_REPLAY"
                    and item["resources"]["status"] == "PASS"}
        aa = {}
        for k in runner.K_VALUES:
            a, b = (accepted.get(f"p{i}_k{k}_{suffix}")
                    for suffix in ("a", "b"))
            if a is None or b is None:
                aa[f"k{k}"] = {"status": "INCOMPLETE"}
            elif (a["base_hash"] != b["base_hash"]
                  or a["root_table_entries"] != b["root_table_entries"]):
                aa[f"k{k}"] = {"status": "BASE_OR_INDEX_MISMATCH"}
            else:
                times = [a["cold_process_ms_exploratory"],
                         b["cold_process_ms_exploratory"]]
                spread = abs(times[0] - times[1]) / statistics.mean(times)
                aa[f"k{k}"] = {"status": "PASS_2_PERCENT" if spread <= 0.02
                                else "NOISY_2_PERCENT",
                                "cold_process_ms": times,
                                "relative_spread": spread}
        point_rows.append({"point_index": i, "public_point": encoded_points[i - 1],
                           "arms": arms, "rho": rho, "aa_control": aa})
    identities = {}
    for k in runner.K_VALUES:
        pipeline = f"k{k}"
        pairs = sorted({(item["base_hash"], item["root_table_entries"])
                        for row in point_rows for item in row["arms"]
                        if item["pipeline"] == pipeline
                        and item["status"] == "PASS_INDEPENDENT_REPLAY"})
        identities[pipeline] = {"identities": [list(pair) for pair in pairs],
                                "status": "PASS" if len(pairs) == 1 else
                                "MISSING" if not pairs else "MISMATCH"}
    evidence_valid = (order_status == "PASS" and telemetry_status == "PASS"
                      and all(all(item["status"] in
                                  ("PASS_INDEPENDENT_REPLAY", "TIMEOUT")
                                  for item in row["arms"])
                              and row["rho"]["status"] == "PASS_INDEPENDENT_REPLAY"
                              for row in point_rows)
                      and all(item["status"] != "MISMATCH"
                              for item in identities.values()))
    admitted = (evidence_valid and all(
        all(item["status"] == "PASS_INDEPENDENT_REPLAY"
            and item["resources"]["status"] == "PASS" for item in row["arms"])
        and row["rho"]["resources"]["status"] == "PASS"
        and all(control["status"] == "PASS_2_PERCENT"
                for control in row["aa_control"].values())
        for row in point_rows))
    times_by_point, paired = [], None
    selected = None
    if admitted:
        for row in point_rows:
            i = row["point_index"]
            by_arm = {item["arm"]: item for item in row["arms"]}
            times = {f"k{k}": statistics.mean(
                by_arm[f"p{i}_k{k}_{suffix}"]["cold_process_ms_exploratory"]
                for suffix in ("a", "b")) for k in runner.K_VALUES}
            times_by_point.append({"point_index": i, "cold_process_ms": times})
        paired = interval([row["cold_process_ms"]["k1104"] -
                           row["cold_process_ms"]["k1120"]
                           for row in times_by_point])
        if paired["lower_95_ms"] > 0:
            selected = "k1120"
        elif paired["upper_95_ms"] < 0:
            selected = "k1104"
    analysis = {
        "schema": "n83.e2e-precision-analysis/v1",
        "hardware_tier": "cloud_instance_tier1",
        "workload_role": "sixteen_fresh_one_point_cold_precision_pairs",
        "point_count": 16, "ic_arm_count": 64, "rho_arm_count": 16,
        "status": "COLD_WINNER_CONFIRMED_WITHIN_PAIR" if selected else
                  "BOUNDED_UNKNOWN",
        "evidence_status": "PASS_INDEPENDENT_REPLAY" if evidence_valid else
                           "VERIFICATION_FAILED",
        "selection_admitted": admitted,
        "selected_pipeline": selected,
        "selected_candidate_id": None,
        "execution_order_status": order_status,
        "topology_status": topology_status,
        "telemetry_status": telemetry_status,
        "host_telemetry_by_arm": telemetry,
        "ordered_events": ordered_events,
        "cross_point_base_identity": identities,
        "cold_process_ms_by_point": times_by_point,
        "paired_interval": paired,
        "point_rows": point_rows,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if admitted else 2


if __name__ == "__main__":
    raise SystemExit(main())
