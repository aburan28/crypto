#!/usr/bin/env python3
"""Replay every N83 cold K-grid answer and retain bounded selection data."""

import argparse
import importlib.util
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_e2e_run", HERE / "run.py")
dynamic = module("n83_e2e_dynamic_verifier",
                 ROOT / "research/n83_root11_dynamic_20261016/verify.py")
dynamic.runner.SOURCE_SHA256 = runner.SOURCE_SHA256
dynamic.runner.parent.EXPECTED[runner.SOURCE] = runner.SOURCE_SHA256
dynamic.coldver.runner.parent.EXPECTED[runner.SOURCE] = runner.SOURCE_SHA256
group = dynamic.group


def one_json(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes = runner.frozen_inputs()
    encoded = one_json(runner.cold.TARGETS)
    answer = one_json(runner.cold.ANSWERS)
    point = group.point_from_json(encoded)
    scalar = int(answer["scalar"])
    if (answer["point"] != encoded or answer["fixture_index"] != 0
            or not group.on_curve(point)
            or group.multiply(point, group.R) is not None
            or group.multiply(group.G, scalar) != point):
        raise ValueError("frozen public point or answer fails independent replay")
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    preflight = one_json(output / "preflight.json")
    if (host["schema"] != "n83.e2e-factor-base-host/v1"
            or host["git_sha"] is None
            or preflight["status"] != "PASS_CAPACITY_PREFLIGHT"
            or source["schema"] != "n83.e2e-factor-base-source/v1"
            or source["hashes"] != hashes
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.parent.digest(HERE / "run.py")
            or source["target_sha256"] != runner.parent.digest(runner.cold.TARGETS)
            or source["arm_order"] != [arm for arm, _ in runner.ARMS]
            or source["k_grid"] != list(runner.K_GRID)
            or source["rho_seed"] != runner.RHO_SEED):
        raise ValueError("source, host, target or grid differs")
    arms = []
    for arm, k in runner.ARMS:
        try:
            arms.append(dynamic.verify_ic(output, arm, k, "root11_dynamic",
                                          point, scalar))
        except Exception as error:
            arms.append({"arm": arm, "k": k, "status": "VERIFICATION_FAILED",
                         "error": str(error)})
    try:
        rho = dynamic.coldver.verify_rho(output, point, scalar)
        if rho["status"] == "PASS_INDEPENDENT_REPLAY":
            receipt = one_json(output / "runs" / "rho" / "receipt.json")
            if (receipt["walk_seed"] != int(runner.RHO_SEED)
                    or receipt["target_count"] != 1
                    or receipt["known_answer_input"] is not False
                    or receipt["validation_sidecar_opened"] is not False
                    or receipt["group_verified"] is not True
                    or receipt["public_input_sha256"] !=
                    runner.parent.digest(runner.cold.TARGETS)
                    or receipt["compiled_source_sha256"]["examples/n83_rho_optimization_eval.rs"]
                    != hashes["examples/n83_rho_optimization_eval.rs"]):
                raise ValueError("same-point rho source, seed or input differs")
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = {row["arm"]: row for row in arms
                if row["status"] == "PASS_INDEPENDENT_REPLAY"
                and row["resources"]["status"] == "PASS"}
    aa = dynamic.repeatability(accepted, "k1152_a", "k1152_b")
    complete = (len(accepted) == len(runner.ARMS)
                and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                and rho["resources"]["status"] == "PASS"
                and aa is not None and aa["status"] == "PASS_10_PERCENT")
    evidence_valid = (all(row["status"] in
                          ("PASS_INDEPENDENT_REPLAY", "TIMEOUT") for row in arms)
                      and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                      and aa is not None
                      and aa["status"] != "BASE_OR_INDEX_MISMATCH")
    cold_ms = {}
    if complete:
        for k in runner.K_GRID:
            cold_ms[str(k)] = (statistics.mean(
                accepted[f"k1152_{suffix}"]["cold_process_ms_exploratory"]
                for suffix in ("a", "b")) if k == runner.CONTROL_K else
                accepted[f"k{k}"]["cold_process_ms_exploratory"])
    descriptive_order = sorted(cold_ms, key=cold_ms.get)
    analysis = {
        "schema": "n83.e2e-factor-base-analysis/v1",
        "hardware_tier": "cloud_instance_tier1_screen",
        "workload_role": "exposed_one_public_point_cold_k_grid",
        "status": "SCREEN_COMPLETE" if complete else "BOUNDED_UNKNOWN",
        "evidence_status": ("PASS_INDEPENDENT_REPLAY" if evidence_valid
                            else "VERIFICATION_FAILED"),
        "selected_candidate_id": None,
        "descriptive_order": descriptive_order,
        "cold_process_ms_by_k": cold_ms,
        "aa_control": aa, "arms": arms, "rho": rho,
        "controlled_online_speedup": None,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if evidence_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
