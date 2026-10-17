#!/usr/bin/env python3
"""Replay every completed cold-neighborhood base, relation and answer."""

import argparse
import importlib.util
import json
import statistics
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


runner = module("n83_neighborhood_runner", RUNNER)
dynamic = module("n83_dynamic_verifier", DYNAMIC_VERIFIER)


def one_json(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, point, scalar = dynamic.coldver.frozen_point()
    if hashes != runner.frozen_inputs():
        raise ValueError("frozen source or target changed")
    host = one_json(output / "host.json")
    source = one_json(output / "source.json")
    if (host["schema"] != "n83.cold-neighborhood-host/v1"
            or source["schema"] != "n83.cold-neighborhood-source/v1"
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["hashes"] != hashes
            or source["runner_sha256"] != runner.parent.digest(RUNNER)
            or source["arm_order"] != [name for name, _, _, _ in runner.ARMS]
            or source["rho_seed"] != runner.cold.RHO_SEED
            or host["git_sha"] is None):
        raise ValueError("executed source, host or arm order differs")
    arms = []
    for name, k, layout, _ in runner.ARMS:
        try:
            arms.append(dynamic.verify_ic(output, name, k, layout, point, scalar))
        except Exception as error:
            arms.append({"arm": name, "k": k, "layout": layout,
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    try:
        rho = dynamic.coldver.verify_rho(output, point, scalar)
    except Exception as error:
        rho = {"status": "VERIFICATION_FAILED", "error": str(error)}
    accepted = [row for row in arms if row["status"] == "PASS_INDEPENDENT_REPLAY"]
    by_arm = {row["arm"]: row for row in accepted
              if row["resources"]["status"] == "PASS"}
    aa = dynamic.repeatability(by_arm, "k1024dynamic_a", "k1024dynamic_b")
    evidence_valid = (len(accepted) > 0
                      and all(row["status"] in ("PASS_INDEPENDENT_REPLAY", "TIMEOUT")
                              for row in arms)
                      and rho["status"] == "PASS_INDEPENDENT_REPLAY"
                      and (aa is None or aa["status"] != "BASE_OR_INDEX_MISMATCH"))
    admitted = (evidence_valid
                and all(row["resources"]["status"] == "PASS" for row in accepted)
                and rho["resources"]["status"] == "PASS"
                and aa is not None and aa["status"] == "PASS_10_PERCENT")
    bounded = [row["arm"] for row in arms if row["status"] == "TIMEOUT"]
    cold_by_pipeline = {}
    if admitted:
        for row in accepted:
            key = f"k{row['k']}{row['layout']}"
            cold_by_pipeline.setdefault(key, []).append(
                row["cold_process_ms_exploratory"])
        cold_by_pipeline = {key: statistics.mean(times)
                            for key, times in cold_by_pipeline.items()}
    ranked = sorted(cold_by_pipeline, key=cold_by_pipeline.get)
    analysis = {
        "schema": "n83.cold-neighborhood-analysis/v1",
        "hardware_tier": "cloud_instance_exploratory",
        "workload_role": "exposed_public_point_cold_neighborhood_screen",
        "target_count": 1, "arm_order": source["arm_order"],
        "status": (("SCREEN_ADMITTED_WITH_TIMEOUT" if bounded else "SCREEN_COMPLETE")
                   if admitted else "BOUNDED_UNKNOWN"),
        "evidence_status": ("PASS_INDEPENDENT_REPLAY" if evidence_valid
                            else "VERIFICATION_FAILED"),
        "selected_factor_base": None, "bounded_arms": bounded,
        "descriptive_cold_process_ms_by_pipeline": cold_by_pipeline,
        "confirmation_candidates": ranked[:2] if admitted else [],
        "aa_control": aa, "arms": arms, "rho": rho,
        "controlled_online_speedup": None,
    }
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if evidence_valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
