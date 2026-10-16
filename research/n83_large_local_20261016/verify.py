#!/usr/bin/env python3
"""Independent K4096 base, target relation and scalar replay, with local limits."""

import argparse
import importlib.util
import json
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


runner = module("n83_large_local_runner", RUNNER)
dynamic = module("n83_large_local_dynamic_verifier", DYNAMIC_VERIFIER)
dynamic.runner.SOURCE_SHA256 = runner.SOURCE_SHA256
dynamic.coldver.runner.parent.EXPECTED[
    dynamic.coldver.runner.parent.SOURCE] = runner.SOURCE_SHA256
dynamic.coldver.runner.RSS_CAP_BYTES = runner.RSS_CAP_BYTES


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    source = read(output / "source.json")
    host = read(output / "host.json")
    usage = read(output / "runs" / "k4096root11dynamic" / "usage.json")
    if (source["schema"] != "n83.large-local-source/v1"
            or host["schema"] != "n83.large-local-host/v1"
            or host["tier"] != "local_diagnostic"
            or source["source_sha256"] != runner.SOURCE_SHA256
            or source["runner_sha256"] != runner.digest(RUNNER)
            or source["k"] != runner.K
            or source["layout"] != "root11_dynamic"
            or source["target_sha256"] != runner.digest(runner.dynamic.cold.TARGETS)
            or source["hashes"] != runner.dynamic.frozen_inputs()
            or host["rss_cap_bytes"] != runner.RSS_CAP_BYTES
            or host["cap_seconds"] != runner.CAP_SECONDS):
        raise ValueError("frozen local source, input or resource limit differs")
    if (usage["peak_rss_bytes"] > runner.RSS_CAP_BYTES
            and usage["status"] != "STOPPED_RSS_GUARD"):
        raise ValueError("RSS guard status differs")
    point_data = read(runner.dynamic.cold.TARGETS)
    answer = read(runner.dynamic.cold.ANSWERS)
    point = dynamic.group.point_from_json(point_data)
    scalar = int(answer["scalar"])
    if (answer["point"] != point_data
            or not dynamic.group.on_curve(point)
            or dynamic.group.multiply(point, dynamic.group.R) is not None
            or dynamic.group.multiply(dynamic.group.G, scalar) != point):
        raise ValueError("frozen public point fails independent subgroup replay")
    state_count = 83 * runner.K * (runner.K - 1) // 2 + 42 * runner.K
    stderr = (output / "runs" / "k4096root11dynamic" / "stderr.log").read_text()
    structural_stop = (usage["status"] == "FAILED"
                       and usage["exit_code"] == 101
                       and state_count >= 1 << 29
                       and "one-root complete key bound exceeds fixed indexed slots"
                           in stderr)
    if structural_stop:
        arm = {"status": "UNSUPPORTED_STATE_ID", "k": runner.K,
               "state_count": state_count, "state_id_limit": 1 << 29,
               "source_guard": "one-root complete key bound exceeds fixed indexed slots"}
    else:
        try:
            arm = dynamic.verify_ic(output, "k4096root11dynamic", runner.K,
                                    "root11_dynamic", point, scalar)
        except Exception as error:
            arm = {"status": "VERIFICATION_FAILED", "error": str(error)}
    admitted = (arm["status"] == "PASS_INDEPENDENT_REPLAY"
                and usage["status"] == "PASS_PRODUCER"
                and usage["peak_rss_bytes"] <= runner.RSS_CAP_BYTES)
    analysis = {
        "schema": "n83.large-local-analysis/v1", "tier": "local_diagnostic",
        "status": ("PASS_LOCAL_DIAGNOSTIC" if admitted else
                   "UNSUPPORTED_STATE_ID" if structural_stop else "BOUNDED_UNKNOWN"),
        "selection_eligible": False,
        "source_sha256": runner.SOURCE_SHA256,
        "binary_sha256": source["binary_sha256"],
        "resource_receipt": usage,
        "arm": arm,
    }
    (output / "analysis.json").write_text(json.dumps(
        analysis, indent=2, sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if admitted else 1


if __name__ == "__main__":
    raise SystemExit(main())
