#!/usr/bin/env python3
"""Independent exact-base, relation and scalar replay for K4096 ID30."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module("n83_id30_k4096_runner", HERE / "run_k4096.py")
control = module("n83_id30_k4096_verifier", HERE / "verify_control.py")
group = control.group


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    source = control.read(output / "source.json")
    host = control.read(output / "host.json")
    preflight = control.read(output / "preflight.json")
    usage = control.read(output / "runs/opt_in_id30/usage.json")
    dynamic = control.runner.dynamic_runner()
    if (source["schema"] != "n83.id30-k4096-source/v1"
            or source["source_sha256"] != control.runner.SOURCE_SHA256
            or source["runner_sha256"] != control.runner.digest(
                HERE / "run_k4096.py")
            or source["control_module_sha256"] != control.runner.digest(
                HERE / "run_control.py")
            or source["binary_sha256"] != control.runner.digest(
                output / "koblitz_orbit_dlp_fast")
            or source["public_input_sha256"] != control.runner.digest(
                dynamic.cold.TARGETS)
            or source["hashes"] != dynamic.frozen_inputs()
            or source["k"] != runner.K
            or source["cap_seconds"] != runner.CAP_SECONDS
            or source["rss_cap_bytes"] != runner.RSS_CAP_BYTES
            or host["physical_memory_bytes"] < 40 * 1024**3
            or host["free_percent_before"] < runner.MIN_START_FREE_PERCENT
            or preflight["status"] != "PASS_CAPACITY_PREFLIGHT"):
        raise ValueError("source, binary, public input or host preflight differs")
    public = control.read(dynamic.cold.TARGETS)
    answer = control.read(dynamic.cold.ANSWERS)
    point = group.point_from_json(public)
    scalar = int(answer["scalar"])
    if (answer["point"] != public or not group.on_curve(point)
            or group.multiply(point, group.R) is not None
            or group.multiply(group.G, scalar) != point):
        raise ValueError("frozen point or answer fails subgroup replay")
    evidence = {"schema": "n83.id30-k4096-analysis/v1",
                "hardware_tier": "local_correctness_diagnostic",
                "k": runner.K,
                "source_sha256": source["source_sha256"],
                "producer_status": usage["status"],
                "peak_sampled_rss_bytes": usage["peak_sampled_rss_bytes"],
                "min_sampled_free_percent": usage["min_sampled_free_percent"],
                "cap_seconds": runner.CAP_SECONDS,
                "performance_claim": None,
                "selected_factor_base": None}
    if usage["status"] == "PASS_PRODUCER":
        control.runner.K = runner.K
        arm = control.verify_arm(output, "opt_in_id30", point, scalar)
        evidence["arm"] = arm
        evidence["status"] = (
            "PASS_EXACT_K4096" if arm["status"] == "PASS_INDEPENDENT_REPLAY"
            and usage["peak_sampled_rss_bytes"] <= runner.RSS_CAP_BYTES
            and usage["wall_seconds"] <= runner.CAP_SECONDS + 2
            and usage["min_sampled_free_percent"] >= runner.MIN_RUNNING_FREE_PERCENT
            else "BOUNDED_UNKNOWN")
    else:
        evidence["status"] = "BOUNDED_UNKNOWN"
        evidence["arm"] = {"status": usage["status"], "usage": usage}
    (output / "analysis.json").write_text(
        json.dumps(evidence, sort_keys=True, indent=2) + "\n")
    print(json.dumps(evidence, sort_keys=True))
    return 0 if evidence["status"] == "PASS_EXACT_K4096" else 1


if __name__ == "__main__":
    raise SystemExit(main())
