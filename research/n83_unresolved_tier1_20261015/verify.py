#!/usr/bin/env python3
"""Independently replay both same-point N83 IC/rho diagnostic pairs."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
RUNNER = HERE / "run.py"
BASE = REPO / "research/n83_final_tier1_20261015/verify.py"

rspec = importlib.util.spec_from_file_location("n83_unresolved_runner", RUNNER)
runner = importlib.util.module_from_spec(rspec)
rspec.loader.exec_module(runner)
bspec = importlib.util.spec_from_file_location("n83_tier1_base_verifier", BASE)
base = importlib.util.module_from_spec(bspec)
bspec.loader.exec_module(base)


def verify_pair(output, point, scalar, host, selected, index):
    pair = output / "pairs" / selected["pipeline"]
    if not (pair / "order.json").exists():
        return {"pipeline": selected["pipeline"], "status": "NOT_RUN"}
    order = base.one_json(pair / "order.json")
    expected_order = ["ic", "rho"] if index == 0 else ["rho", "ic"]
    point_file = output / "inputs" / "point_1.jsonl"
    if (order["pipeline"] != selected["pipeline"]
            or order["candidate_id"] != selected["candidate_id"]
            or order["order"] != expected_order
            or order["point_sha256"] != runner.parent.digest(point_file)
            or base.json_lines(point_file) != base.json_lines(runner.TARGETS)):
        raise ValueError("same-point pairing or frozen arm order differs")
    ic = base.verify_ic(pair / "ic", point, scalar, selected)
    rho = base.verify_rho(pair / "rho", point, scalar)
    if ic["status"] == "PASS_INDEPENDENT_REPLAY":
        target = base.one_json(pair / "ic" / "targets_result.jsonl")
        phase_keys = ("target_query_ms", "target_pdp_ms",
                      "target_relation_check_ms", "target_descent_ms",
                      "target_recovery_check_ms")
        ic["phases_ms"] = {key: target[key] for key in phase_keys}
    noise = base.pair_noise(pair, host, ic, rho)
    passed = (ic["status"] == rho["status"] == "PASS_INDEPENDENT_REPLAY"
              and noise["status"] == "PASS")
    return {"pipeline": selected["pipeline"],
            "candidate_id": selected["candidate_id"],
            "status": "PASS_TIER1_PAIR" if passed else "PAIR_INCOMPLETE_OR_NOISY",
            "order": expected_order, "ic": ic, "rho": rho, "noise": noise,
            "ratio_rho_over_ic": (rho["online_ms"] / ic["online_ms"] if passed else None)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, encoded, answer, candidates = runner.frozen()
    point = base.group.point_from_json(encoded)
    scalar = int(answer["scalar"])
    if (answer["fixture_index"] != 0
            or base.group.point_from_json(answer["point"]) != point
            or not base.group.on_curve(point)
            or base.group.multiply(point, base.group.R) is not None
            or base.group.multiply(base.group.G, scalar) != point):
        raise ValueError("one public point or scalar sidecar fails independent replay")
    host = base.one_json(output / "host.json")
    source = base.one_json(output / "source.json")
    if (source["schema"] != "n83.unresolved-tier1-source/v1"
            or source["hashes"] != hashes
            or source["candidate_order"] != candidates
            or source["confirmation_selection_sha256"] != runner.parent.digest(runner.RESULT)
            or source["runner_sha256"] != runner.parent.digest(RUNNER)
            or host["schema"] != "n83.unresolved-tier1-host/v1"):
        raise ValueError("executed source, host or candidate order differs")
    rows = []
    for index, selected in enumerate(candidates):
        try:
            rows.append(verify_pair(output, point, scalar, host, selected, index))
        except Exception as error:
            rows.append({"pipeline": selected["pipeline"],
                         "candidate_id": selected["candidate_id"],
                         "status": "VERIFICATION_FAILED", "error": str(error)})
    end_path = output / "allocation_end.json"
    allocation_end = base.one_json(end_path) if end_path.exists() else None
    same_allocation = (allocation_end is not None
                       and allocation_end["instance_id"] == host["instance_id"]
                       and allocation_end["visible_affinity"] ==
                           sorted(host["execution_cpus"]))
    complete = same_allocation and all(row["status"] == "PASS_TIER1_PAIR" for row in rows)
    results = {
        "schema": "n83.unresolved-tier1-analysis/v1",
        "hardware_tier": "tier1_cloud_instance",
        "status": "PASS_TIER1_SAME_POINT_PAIR_SET" if complete else "TIER1_GATE_INCOMPLETE",
        "target_count": 1,
        "factor_base_selection_status": "BOUNDED_UNKNOWN",
        "candidate_order": [row["pipeline"] for row in candidates],
        "allocation_end": allocation_end,
        "same_allocation": same_allocation,
        "pairs": rows,
        "online_ratio_rho_over_ic": {
            row["pipeline"]: row["ratio_rho_over_ic"] if row["status"] == "PASS_TIER1_PAIR"
            else None for row in rows},
        "descriptive_ic_online_ms": {
            row["pipeline"]: row["ic"]["online_ms"] if row["status"] == "PASS_TIER1_PAIR"
            else None for row in rows},
    }
    (output / "analysis.json").write_text(json.dumps(results, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps({"status": results["status"],
                      "online_ratio_rho_over_ic": results["online_ratio_rho_over_ic"]},
                     sort_keys=True))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
