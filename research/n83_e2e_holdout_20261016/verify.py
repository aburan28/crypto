#!/usr/bin/env python3
"""Independent five-point cold N83 factor-base holdout replay and decision."""

import argparse
import importlib.util
import json
import math
import statistics
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T_CRITICAL_DF4 = 2.7764451051977987


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_e2e_holdout_run", HERE / "run.py")
dynamic = module("n83_e2e_holdout_dynamic_verify",
                 ROOT / "research/n83_root11_dynamic_20261016/verify.py")
dynamic.runner.SOURCE_SHA256 = runner.screen.SOURCE_SHA256
dynamic.runner.parent.EXPECTED[runner.screen.SOURCE] = runner.screen.SOURCE_SHA256
dynamic.coldver.runner.parent.EXPECTED[runner.screen.SOURCE] = runner.screen.SOURCE_SHA256
group = dynamic.group


def one_json(path):
    return json.loads(path.read_text())


def points_and_answers():
    hashes = runner.fixture()
    points = [json.loads(line) for line in runner.TARGETS.read_text().splitlines()]
    answers = [json.loads(line) for line in runner.ANSWERS.read_text().splitlines()]
    if len(points) != 5 or len(answers) != 5 or len({tuple(p) for p in points}) != 5:
        raise ValueError("five distinct frozen points required")
    scalars = []
    for i, (encoded, answer) in enumerate(zip(points, answers)):
        point = group.point_from_json(encoded)
        scalar = int(answer["scalar"])
        if (answer["fixture_index"] != i or answer["point"] != encoded
                or not group.on_curve(point)
                or group.multiply(point, group.R) is not None
                or group.multiply(group.G, scalar) != point):
            raise ValueError(f"holdout point {i} fails scalar replay")
        scalars.append(scalar)
    return hashes, points, scalars


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


def paired_interval(differences):
    mean = statistics.mean(differences)
    sd = statistics.stdev(differences)
    halfwidth = T_CRITICAL_DF4 * sd / math.sqrt(len(differences))
    return {"differences_ms": differences, "mean_ms": mean,
            "sample_sd_ms": sd, "lower_95_ms": mean - halfwidth,
            "upper_95_ms": mean + halfwidth}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, encoded_points, scalars = points_and_answers()
    challengers, selection = runner.finalists()
    schedule = runner.schedule(challengers)
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    preflight = one_json(output / "preflight.json")
    if (host["schema"] != "n83.e2e-holdout-host/v1"
            or host["git_sha"] is None
            or preflight["status"] != "PASS_CAPACITY_PREFLIGHT"
            or source["schema"] != "n83.e2e-holdout-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.screen.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(HERE / "run.py")
            or source["selection_sha256"] != runner.parent.digest(runner.SELECTION)
            or source["screen_analysis_sha256"] != runner.parent.digest(
                runner.SCREEN_ANALYSIS)
            or source["screen_archive_sha256"] !=
            selection["screen_archive_sha256"]
            or source["challenger_ks"] != challengers
            or source["schedule"] != schedule
            or source["rho_seed_base"] != runner.RHO_SEED_BASE):
        raise ValueError("host, source, selection or schedule differs")
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
            k = runner.screen.CONTROL_K if arm.startswith("k1152_") else int(arm[1:])
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
        aa = dynamic.repeatability(accepted, f"p{i}_k1152_a", f"p{i}_k1152_b")
        point_rows.append({"point_index": i, "public_point": encoded_points[i - 1],
                           "arms": arms, "rho": rho, "aa_control": aa})
    pipelines = ("k1152", *(f"k{k}" for k in challengers))
    identities = {}
    for pipeline in pipelines:
        pairs = sorted({(item["base_hash"], item["root_table_entries"])
                        for row in point_rows for item in row["arms"]
                        if item["pipeline"] == pipeline
                        and item["status"] == "PASS_INDEPENDENT_REPLAY"})
        identities[pipeline] = {"identities": [list(pair) for pair in pairs],
                                "status": "PASS" if len(pairs) == 1 else
                                "MISSING" if not pairs else "MISMATCH"}
    evidence_valid = (order_status == "PASS" and all(
        all(item["status"] in ("PASS_INDEPENDENT_REPLAY", "TIMEOUT")
            for item in row["arms"])
        and row["rho"]["status"] == "PASS_INDEPENDENT_REPLAY"
        and (row["aa_control"] is None
             or row["aa_control"]["status"] != "BASE_OR_INDEX_MISMATCH")
        for row in point_rows)
        and all(item["status"] != "MISMATCH" for item in identities.values()))
    admitted = (evidence_valid and all(
        all(item["status"] == "PASS_INDEPENDENT_REPLAY"
            and item["resources"]["status"] == "PASS" for item in row["arms"])
        and row["rho"]["resources"]["status"] == "PASS"
        and row["aa_control"] is not None
        and row["aa_control"]["status"] == "PASS_10_PERCENT"
        for row in point_rows))
    times_by_point, mean_by_pipeline, paired = [], {}, {}
    winner = None
    if admitted:
        for row in point_rows:
            i = row["point_index"]
            by_arm = {item["arm"].split(f"p{i}_", 1)[1]: item
                      for item in row["arms"]}
            times = {f"k{k}": by_arm[f"k{k}"]["cold_process_ms_exploratory"]
                     for k in challengers}
            times["k1152"] = statistics.mean(
                by_arm[f"k1152_{suffix}"]["cold_process_ms_exploratory"]
                for suffix in ("a", "b"))
            times_by_point.append({"point_index": i, "cold_process_ms": times})
        mean_by_pipeline = {pipeline: statistics.mean(
            row["cold_process_ms"][pipeline] for row in times_by_point)
            for pipeline in pipelines}
        for candidate in pipelines:
            for other in pipelines:
                if candidate != other:
                    paired[f"{candidate}|{other}"] = paired_interval([
                        row["cold_process_ms"][other] -
                        row["cold_process_ms"][candidate]
                        for row in times_by_point])
        fastest = min(pipelines, key=mean_by_pipeline.get)
        if all(
            paired[f"{fastest}|{other}"]["lower_95_ms"] > 0
            and all(difference > abs(row["aa_control"]["cold_process_ms"][0]
                                     - row["aa_control"]["cold_process_ms"][1]) / 2
                    for difference, row in zip(
                        paired[f"{fastest}|{other}"]["differences_ms"], point_rows))
            for other in pipelines if other != fastest):
            winner = fastest
    analysis = {
        "schema": "n83.e2e-holdout-analysis/v1",
        "hardware_tier": "cloud_instance_tier1",
        "workload_role": "five_fresh_one_point_cold_confirmations",
        "point_count": 5, "arm_count": 25,
        "status": ("COLD_WINNER_CONFIRMED_WITHIN_GRID" if winner else
                   "BOUNDED_UNKNOWN"),
        "evidence_status": ("PASS_INDEPENDENT_REPLAY" if evidence_valid
                            else "VERIFICATION_FAILED"),
        "selection_admitted": admitted,
        "selected_pipeline": winner, "selected_candidate_id": None,
        "execution_order_status": order_status,
        "ordered_events": ordered_events,
        "cross_point_base_identity": identities,
        "cold_process_ms_by_point": times_by_point,
        "mean_cold_process_ms_by_pipeline": mean_by_pipeline,
        "paired_intervals": paired,
        "point_rows": point_rows,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if evidence_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
