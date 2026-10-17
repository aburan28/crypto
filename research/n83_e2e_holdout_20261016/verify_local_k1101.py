#!/usr/bin/env python3
"""Independently replay the archived K1101 local diagnostic."""

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast.rs"
SOURCE_SHA256 = "8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12"
COLD = ROOT / "research/n83_cold_factor_base_20261016"
spec = importlib.util.spec_from_file_location("n83_k1101_group", COLD / "verify.py")
coldver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coldver)
group = coldver.group


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one_json(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    folder = args.folder.resolve()
    command = one_json(folder / "command.json")
    usage = one_json(folder / "usage.json")
    summary = one_json(folder / "stdout.jsonl")
    target = one_json(folder / "targets_result.jsonl")
    dump = one_json(folder / "base_dump.jsonl")
    binary = folder / "koblitz_orbit_dlp_fast"
    if (sha256(SOURCE) != SOURCE_SHA256
            or command["source_sha256"] != SOURCE_SHA256
            or sha256(binary) != command["binary_sha256"]
            or usage["exit_code"] != 0
            or usage["tier"] != "local_diagnostic"
            or (folder / "stderr.txt").read_text()
            or summary["n"] != 83 or summary["a"] != 1
            or summary["orbit_columns"] != 1101
            or summary["rank"] != 1101
            or summary["factor_base_points"] != 166 * 1101
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_allocated_bytes"] != 671088640
            or target["group_verified"] is not True):
        raise ValueError("local source, binary, base, index or execution differs")
    base = group.replay_base(dump, summary)
    point = group.point_from_json(one_json(COLD / "targets.jsonl"))
    scalar = int(one_json(COLD / "answers.jsonl")["scalar"])
    relation = None
    indices = target["point_indices"]
    if len(indices) != 4 or len(base) != 166 * 1101:
        raise ValueError("base or relation width differs")
    for index in indices:
        if not isinstance(index, int) or not 0 <= index < len(base):
            raise ValueError("relation index outside enumerated base")
        relation = group.add(relation, base[index])
    recovered = int(target["recovered_scalar"])
    phases = ("target_query_ms", "target_pdp_ms",
              "target_relation_check_ms", "target_descent_ms",
              "target_recovery_check_ms")
    timing = summary["timing_ms"]
    if (relation != point or recovered != scalar
            or group.multiply(group.G, recovered) != point
            or not math.isclose(sum(target[key] for key in phases),
                                target["target_ms"], abs_tol=1.0)
            or abs(timing["process_total"] - timing["setup_total"]
                   - timing["rank_stage"] - timing["targets_total"]) > 1000
            or usage["outer_wall_seconds"] * 1000 < timing["process_total"]):
        raise ValueError("point, scalar or exclusive phase replay failed")
    record = {
        "status": "PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY",
        "tier": "local_diagnostic", "k": 1101,
        "actual_factor_base_points": len(base),
        "folded_columns": summary["rank"],
        "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
        "base_hash": summary["base_hash"],
        "recovered_scalar": str(recovered), "point_indices": indices,
        "process_total_ms": timing["process_total"],
        "setup_ms": timing["setup_total"],
        "rank_stage_ms": timing["rank_stage"],
        "target_ms": target["target_ms"],
        "source_sha256": SOURCE_SHA256,
    }
    if record != one_json(folder / "independent_replay.json"):
        raise ValueError("archived independent replay receipt differs")
    print(json.dumps(record, sort_keys=True))


if __name__ == "__main__":
    main()
