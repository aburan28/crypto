#!/usr/bin/env python3
"""Replay the base, four-point relation and scalar in each compact-table arm."""

import argparse
import importlib.util
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location("n83_compact_local_run", HERE / "local_pair.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
group_spec = importlib.util.spec_from_file_location(
    "n83_compact_local_group", ROOT / "research/n83_indexed5_20261013/analyze.py")
group = importlib.util.module_from_spec(group_spec)
group_spec.loader.exec_module(group)


def one_json(path):
    return json.loads(path.read_text())


def verify_arm(output, arm, point, scalar):
    folder = output / "runs" / arm
    compact = arm.startswith("compact")
    usage = one_json(folder / "usage.json")
    if usage["status"] != "PASS_PRODUCER":
        return {"arm": arm, "status": usage["status"], "usage": usage}
    summary = one_json(folder / "stdout.jsonl")
    target = one_json(folder / "targets_result.jsonl")
    dump = one_json(folder / "base_dump.jsonl")
    command = one_json(folder / "command.json")
    expected_bytes = 5 * (1 << (27 if compact else 28))
    if (command["compact27"] is not compact
            or command["env_flags"]["KIC_INDEXED5_COMPACT27"] !=
            ("1" if compact else "0")
            or summary["n"] != 83 or summary["a"] != 1
            or summary["orbit_columns"] != runner.K
            or summary["rank"] != runner.K
            or summary["factor_base_points"] != 166 * runner.K
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_allocated_bytes"] != expected_bytes
            or summary["indexed5_compact27_enabled"] is not compact
            or summary["root_table_entries"] >= (1 << (27 if compact else 28))
            or summary["targets_solved"] != 1
            or target["group_verified"] is not True
            or usage["peak_rss_bytes"] is None
            or usage["peak_rss_bytes"] <= 0):
        raise ValueError(f"{arm}: exact layout, rank or telemetry differs")
    base = group.replay_base(dump, summary)
    if len(base) != 166 * runner.K:
        raise ValueError(f"{arm}: base replay differs")
    relation = None
    indices = target["point_indices"]
    if len(indices) != 4:
        raise ValueError(f"{arm}: four-point relation expected")
    for index in indices:
        if not isinstance(index, int) or not 0 <= index < len(base):
            raise ValueError(f"{arm}: relation index outside base")
        relation = group.add(relation, base[index])
    phases = ("target_query_ms", "target_pdp_ms",
              "target_relation_check_ms", "target_descent_ms",
              "target_recovery_check_ms")
    recovered = int(target["recovered_scalar"])
    timing = summary["timing_ms"]
    total = float(timing["process_total"])
    if (relation != point or recovered != scalar
            or group.multiply(group.G, recovered) != point
            or not math.isclose(sum(target[key] for key in phases),
                                target["target_ms"], abs_tol=1.0)
            or abs(total - timing["setup_total"] - timing["rank_stage"]
                   - timing["targets_total"]) > 1000
            or usage["outer_wall_seconds"] * 1000 < total):
        raise ValueError(f"{arm}: scalar, relation or exclusive costs differ")
    return {"arm": arm, "status": "PASS_INDEPENDENT_REPLAY",
            "compact27": compact, "base_hash": summary["base_hash"],
            "root_table_entries": summary["root_table_entries"],
            "root_table_allocated_bytes": expected_bytes,
            "factor_base_points": len(base), "process_total_ms": total,
            "setup_ms": timing["setup_total"],
            "rank_stage_ms": timing["rank_stage"],
            "target_ms": target["target_ms"],
            "peak_rss_bytes": usage["peak_rss_bytes"],
            "outer_wall_seconds": usage["outer_wall_seconds"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    source = one_json(output / "source.json")
    if (source["source_sha256"] != runner.SOURCE_SHA256
            or source["binary_sha256"] != runner.sha256(
                Path(one_json(output / "runs/off_a/command.json")["argv"][0]))
            or source["target_sha256"] != runner.sha256(runner.TARGETS)
            or source["answer_sidecar_sha256"] != runner.sha256(runner.ANSWERS)
            or source["runner_sha256"] != runner.sha256(HERE / "local_pair.py")
            or source["arm_order"] != list(runner.ARMS)
            or source["k"] != runner.K
            or source["rank_seed"] != 29):
        raise ValueError("source, binary or frozen inputs differ")
    point = group.point_from_json(one_json(runner.TARGETS))
    scalar = int(one_json(runner.ANSWERS)["scalar"])
    if group.multiply(group.G, scalar) != point:
        raise ValueError("frozen answer sidecar differs")
    arms = []
    for arm in runner.ARMS:
        try:
            arms.append(verify_arm(output, arm, point, scalar))
        except Exception as error:
            arms.append({"arm": arm, "status": "VERIFICATION_FAILED",
                         "error": str(error)})
    admitted = all(row["status"] == "PASS_INDEPENDENT_REPLAY" for row in arms)
    aa = None
    comparison = None
    if admitted:
        if (len({row["base_hash"] for row in arms}) != 1
                or len({row["root_table_entries"] for row in arms}) != 1):
            admitted = False
        else:
            times = {row["arm"]: row["process_total_ms"] for row in arms}
            off = (times["off_a"], times["off_b"])
            compact = (times["compact_a"], times["compact_b"])
            aa = {"off_spread": abs(off[0] - off[1]) / statistics.mean(off),
                  "compact_spread": abs(compact[0] - compact[1]) /
                  statistics.mean(compact)}
            off_mean = statistics.mean(off)
            compact_mean = statistics.mean(compact)
            comparison = {"off_mean_ms": off_mean,
                          "compact_mean_ms": compact_mean,
                          "off_over_compact_ratio": off_mean / compact_mean,
                          "advantage_ms": off_mean - compact_mean,
                          "max_aa_difference_ms": max(abs(off[0] - off[1]),
                                                      abs(compact[0] - compact[1]))}
    analysis = {"schema": "n83.indexed5-compact-local-analysis/v1",
                "tier": "local_diagnostic", "k": runner.K,
                "status": "PASS_INDEPENDENT_REPLAY" if admitted
                else "BOUNDED_UNKNOWN", "arms": arms, "aa": aa,
                "comparison": comparison}
    (output / "analysis.json").write_text(json.dumps(analysis, indent=2,
                                                    sort_keys=True) + "\n")
    print(json.dumps(analysis, sort_keys=True))
    return 0 if admitted else 1


if __name__ == "__main__":
    raise SystemExit(main())
