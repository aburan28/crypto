"""Bind the old and new ledgers, then bootstrap 33 matched curve indices."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import statistics
import tarfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
FREEZE = HERE / "FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric_from_logs(logs: list[float]) -> float:
    return math.exp(statistics.mean(logs))


def interval(logs: list[float], rng: random.Random) -> list[float]:
    draws = sorted(geometric_from_logs([logs[rng.randrange(len(logs))]
                                        for _ in logs]) for _ in range(20000))
    return [draws[499], draws[19499]]


def analyze(extension: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen_raw = FREEZE.read_bytes()
    frozen = json.loads(frozen_raw)
    assert frozen["kind"] == "s3_crosscurve_extension33_freeze_v1"
    assert frozen["analyzer_sha256"] == sha(Path(__file__).read_bytes())
    old_archive = OLD / "local/raw_pilot24.tar.gz"
    assert sha(old_archive.read_bytes()) == frozen["old_raw_archive_sha256"]
    with tarfile.open(old_archive, "r:gz") as tape:
        old_raw = tape.extractfile("raw.jsonl").read()
    old_analysis = OLD / "local/analysis_pilot24.json"
    assert sha(old_analysis.read_bytes()) == frozen["old_analysis_sha256"]
    assert json.loads(old_analysis.read_bytes())["ledger_sha256"] == sha(old_raw)
    old_rows = [json.loads(line) for line in old_raw.splitlines()]
    assert len(old_rows) == 192 and all(row["verified_native_and_fixture"] for row in old_rows)

    completion_raw = (extension / "completion.json").read_bytes()
    completion = json.loads(completion_raw)
    assert completion["freeze_sha256"] == sha(frozen_raw)
    new_raw = (extension / "raw.jsonl").read_bytes()
    assert completion["ledger_sha256"] == sha(new_raw)
    new_rows = [json.loads(line) for line in new_raw.splitlines()]
    expected = [(row["n"], row["label"], row["arm"], row["run_id"])
                for row in frozen["run_order"]]
    actual = [(row["n"], row["label"], row["arm"], row["run_id"])
              for row in new_rows]
    assert actual == expected[:len(actual)]
    assert all(row["binary_sha256"] == frozen["binary_sha256"][row["arm"]]
               for row in old_rows + new_rows)
    assert len({row["run_id"] for row in old_rows + new_rows}) == len(old_rows + new_rows)
    complete = (completion["complete"] and len(new_rows) == 72
                and all(row["verified_native_and_fixture"] for row in new_rows))
    output.mkdir(parents=True)
    result = {
        "kind": "s3_crosscurve_extension33_analysis_v1",
        "freeze_sha256": sha(frozen_raw),
        "old_ledger_sha256": sha(old_raw),
        "extension_ledger_sha256": sha(new_raw),
        "completion_sha256": sha(completion_raw),
        "old_run_count": len(old_rows),
        "extension_run_count": len(new_rows),
        "complete": complete,
        "matched_index_count": None,
        "operation_count": None,
        "online_wall_exploratory": None,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    if complete:
        pooled = old_raw + new_raw
        assert old_raw.endswith(b"\n") and new_raw.endswith(b"\n")
        (output / "raw.jsonl").write_bytes(pooled)
        spec = importlib.util.spec_from_file_location(
            "crosscurve_old_analyzer", OLD / "analyze_local_pilot.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        old_result = module.analyze(output / "raw.jsonl", output / "prior_analysis.json")
        assert old_result["matched_index_count"] == 33
        assert old_result["run_count"] == 264
        assert old_result["binary_sha256"] == frozen["binary_sha256"]
        blocks = {(block["n"], block["label"]): block for block in old_result["blocks"]}
        assert len(blocks) == 66 and all(block["all_verified"] for block in blocks.values())
        per_curve = {"41": {"root_logs": [], "wall_logs": []},
                     "53": {"root_logs": [], "wall_logs": []}}
        root_difference = []
        wall_difference = []
        for index in range(1, 34):
            label = f"T{index:03}"
            a, b = blocks[41, label], blocks[53, label]
            ar, br = math.log(a["online_root_call_ratio"]), math.log(b["online_root_call_ratio"])
            aw, bw = math.log(a["exploratory_online_ratio"]), math.log(b["exploratory_online_ratio"])
            per_curve["41"]["root_logs"].append(ar)
            per_curve["53"]["root_logs"].append(br)
            per_curve["41"]["wall_logs"].append(aw)
            per_curve["53"]["wall_logs"].append(bw)
            root_difference.append(ar - br)
            wall_difference.append(aw - bw)
        rng = random.Random(20261010)
        operation = {}
        online = {}
        for degree in ("41", "53"):
            roots = per_curve[degree]["root_logs"]
            walls = per_curve[degree]["wall_logs"]
            operation[degree] = {"reference_over_prepared_geomean": geometric_from_logs(roots),
                                 "target_bootstrap95": interval(roots, rng),
                                 "min_point_ratio": math.exp(min(roots))}
            online[degree] = {"reference_over_prepared_geomean": geometric_from_logs(walls),
                              "target_bootstrap95": interval(walls, rng),
                              "min_point_ratio": math.exp(min(walls))}
        operation["n41_over_n53_gain"] = {
            "geomean_ratio": geometric_from_logs(root_difference),
            "matched_index_bootstrap95": interval(root_difference, rng),
        }
        online["n41_over_n53_gain"] = {
            "geomean_ratio": geometric_from_logs(wall_difference),
            "matched_index_bootstrap95": interval(wall_difference, rng),
        }
        result.update({"matched_index_count": 33,
                       "pooled_ledger_sha256": sha(pooled),
                       "prior_analysis_sha256": sha((output / "prior_analysis.json").read_bytes()),
                       "operation_count": operation,
                       "online_wall_exploratory": online,
                       "bootstrap_draws": 20000,
                       "bootstrap_seed": 20261010})
    (output / "analysis.json").write_text(json.dumps(result, indent=2,
                                                      sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in (
        "old_run_count", "extension_run_count", "complete", "matched_index_count",
        "operation_count", "online_wall_exploratory")}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.extension.resolve(), args.output.resolve())
