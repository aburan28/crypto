"""Audit the 136-run N53 confirmation panel and its single-target ratios."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
PANEL = HERE / "PANEL_FREEZE.json"
EXECUTION = HERE / "EXECUTION_FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def archive(run_path: Path, receipt: dict, name: str) -> bytes:
    packed = (run_path / (name + ".gz")).read_bytes()
    expected = receipt["archive"][name]
    assert sha(packed) == expected["gzip_sha256"]
    raw = gzip.decompress(packed)
    assert sha(raw) == expected["raw_sha256"]
    assert len(raw) == expected["raw_bytes"]
    return raw


def geomean(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(math.log(value) for value in values))


def quantile(sorted_values: list[float], fraction: float) -> float:
    position = (len(sorted_values) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    return sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * (
        position - lower)


def analyze(run_dir: Path, output: Path, sage_replay: Path | None) -> dict:
    if output.exists():
        raise FileExistsError(output)
    panel_raw = PANEL.read_bytes()
    panel = json.loads(panel_raw)
    execution_raw = EXECUTION.read_bytes()
    execution = json.loads(execution_raw)
    assert panel["kind"] == "s3_n53_c1800_confirmation34_panel_freeze_v1"
    assert execution["panel_freeze_sha256"] == sha(panel_raw)
    assert execution["analyzer_sha256"] == sha(Path(__file__).read_bytes())
    completion_raw = (run_dir / "completion.json").read_bytes()
    completion = json.loads(completion_raw)
    assert completion["panel_freeze_sha256"] == sha(panel_raw)
    assert completion["execution_freeze_sha256"] == sha(execution_raw)
    assert completion["planned_run_count"] == 136
    ledger_raw = (run_dir / "ledger.jsonl").read_bytes()
    ledger = [json.loads(line) for line in ledger_raw.splitlines()]
    assert len(ledger) == completion["completed_run_count"] <= 136
    schedule = [(item["label"], item["arm"], item["repetition"])
                for item in panel["run_order"]]
    assert len(schedule) == 136
    targets = {target["label"]: target for target in panel["targets"]}
    rows: dict[tuple[str, str], list[dict]] = {}
    run_rows = []
    failures = []
    for index, entry in enumerate(ledger):
        label, arm, repeat = schedule[index]
        assert (entry["label"], entry["arm"], entry["repetition"]) == (
            label, arm, repeat)
        assert entry["index"] == index
        receipt_path = run_dir / entry["relative_receipt"]
        receipt_raw = receipt_path.read_bytes()
        assert sha(receipt_raw) == entry["receipt_sha256"]
        receipt = json.loads(receipt_raw)
        target = targets[label]
        frozen_arm = panel["arms"][arm]
        assert receipt["panel_freeze_sha256"] == sha(panel_raw)
        assert receipt["execution_freeze_sha256"] == sha(execution_raw)
        assert receipt["candidate_id"] == entry["candidate_id"] == frozen_arm["candidate_id"]
        assert receipt["workload_id"] == entry["workload_id"] == target["workload_id"]
        assert receipt["run_id"] == entry["run_id"] == (
            frozen_arm["candidate_id"] + "W" + target["workload_id"] + "R" + str(repeat))
        assert receipt["point_sha256"] == target["point_sha256"]
        assert receipt["source_sha256"] == frozen_arm["source_sha256"]
        assert receipt["binary_sha256"] == execution["binaries"][arm]["sha256"]
        report_raw = archive(receipt_path.parent, receipt, "report.json")
        stdout_raw = archive(receipt_path.parent, receipt, "stdout.bin")
        stderr_raw = archive(receipt_path.parent, receipt, "stderr.txt")
        verified = bool(entry["verified"] and entry["memory_within_envelope"]
                        and receipt["verified_native_and_fixture"]
                        and receipt["memory_within_envelope"]
                        and receipt["phase_accounting_ok"]
                        and receipt["exit_code"] == 0
                        and stdout_raw == report_raw
                        and stderr_raw == b"")
        try:
            report = json.loads(report_raw) if report_raw else None
        except ValueError:
            report = None
        if verified:
            assert report is not None
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["factor_base_digest"] == frozen_arm["factor_base_digest"]
            assert report["factor_base_points"] == frozen_arm["factor_base_points"]
            assert report["prepared_rank"] == frozen_arm["prepare_rank"]
            assert report["target_s3_calls"] == receipt["target_s3_calls"]
            assert receipt["online_ms"] == entry["online_ms"] > 0
            assert abs(sum(receipt["target_phase_ms"][phase] for phase in PHASES) -
                       receipt["online_ms"]) < .02
        else:
            failures.append({"index": index, "label": label, "arm": arm,
                             "repetition": repeat, "receipt_sha256": sha(receipt_raw),
                             "exit_code": receipt["exit_code"], "error": receipt["error"]})
        row = {"index": index, "label": label, "arm": arm, "repetition": repeat,
               "candidate_id": receipt["candidate_id"],
               "workload_id": receipt["workload_id"],
               "run_id": receipt["run_id"],
               "verified": verified,
               "online_ms": receipt["online_ms"],
               "target_phase_ms": receipt["target_phase_ms"],
               "target_s3_calls": receipt["target_s3_calls"],
               "reusable_setup_ms": receipt["reusable_setup_ms"],
               "peak_rss_bytes": receipt["peak_rss_bytes"],
               "receipt_sha256": sha(receipt_raw),
               "report_raw_sha256": sha(report_raw)}
        run_rows.append(row)
        rows.setdefault((label, arm), []).append(row)

    sage_verified = None
    sage_sha = None
    if sage_replay is not None:
        sage_raw = sage_replay.read_bytes()
        sage = json.loads(sage_raw)
        assert sage["panel_freeze_sha256"] == sha(panel_raw)
        assert sage["execution_freeze_sha256"] == sha(execution_raw)
        sage_rows = {(row["label"], row["arm"], row["repetition"]): row
                     for row in sage["runs"]}
        assert len(sage_rows) == len(sage["runs"])
        sage_verified = bool(sage["generator_order_verified"]
                             and sage["verified_run_count"] == 136
                             and len(sage_rows) == 136
                             and sage["replay_source_sha256"] == execution["sage_replay_sha256"]
                             and all((row["label"], row["arm"], row["repetition"]) in sage_rows
                                     and sage_rows[(row["label"], row["arm"], row["repetition"])][
                                         "run_id"] == row["run_id"]
                                     and sage_rows[(row["label"], row["arm"], row["repetition"])][
                                         "report_raw_sha256"] == row["report_raw_sha256"]
                                     and sage_rows[(row["label"], row["arm"], row["repetition"])][
                                         "scalar_times_generator_equals_public_target"]
                                     for row in run_rows))
        sage_sha = sha(sage_raw)
    complete = bool(len(ledger) == 136 and not failures
                    and completion["stopped_reason"] is None
                    and completion["verified_run_count"] == 136
                    and completion["memory_pass_count"] == 136)
    points = []
    ratio = None
    interval = None
    local_gate = False
    if complete:
        assert set(rows) == {(label, arm) for label in panel["target_labels"]
                             for arm in ("baseline", "combined")}
        for label in panel["target_labels"]:
            baseline = rows[(label, "baseline")]
            combined = rows[(label, "combined")]
            assert len(baseline) == len(combined) == 2
            bms = [row["online_ms"] for row in baseline]
            cms = [row["online_ms"] for row in combined]
            points.append({"label": label, "workload_id": targets[label]["workload_id"],
                           "baseline_online_ms": bms, "combined_online_ms": cms,
                           "baseline_over_combined_online_ratio": geomean(bms) / geomean(cms),
                           "baseline_target_s3_calls": [row["target_s3_calls"] for row in baseline],
                           "combined_target_s3_calls": [row["target_s3_calls"] for row in combined]})
        logs = [math.log(point["baseline_over_combined_online_ratio"]) for point in points]
        ratio = math.exp(statistics.mean(logs))
        rng = random.Random(20261010)
        bootstrap = sorted(math.exp(sum(rng.choice(logs) for _ in logs) / len(logs))
                           for _ in range(20000))
        interval = [quantile(bootstrap, .025), quantile(bootstrap, .975)]
        local_gate = bool(sage_verified and ratio >= 2.0 and interval[0] > 2.0)
    result = {
        "kind": "s3_n53_c1800_confirmation34_analysis_v1",
        "panel_freeze_sha256": sha(panel_raw),
        "execution_freeze_sha256": sha(execution_raw),
        "completion_sha256": sha(completion_raw),
        "ledger_sha256": sha(ledger_raw),
        "sage_replay_sha256": sage_sha,
        "independent_sage_replay_verified": sage_verified,
        "planned_run_count": 136,
        "completed_run_count": len(ledger),
        "verified_run_count": sum(row["verified"] for row in run_rows),
        "distinct_target_count": len(points),
        "complete": complete,
        "failed_runs": failures,
        "unstarted_run_count": 136 - len(ledger),
        "online_geomean_baseline_over_combined": ratio,
        "online_bootstrap95_target_resample": interval,
        "local_twofold_gate_pass": local_gate,
        "points": points,
        "runs": run_rows,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("completed_run_count", "verified_run_count", "complete",
                       "online_geomean_baseline_over_combined",
                       "online_bootstrap95_target_resample", "local_twofold_gate_pass")},
                     sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sage-replay", type=Path)
    args = parser.parse_args()
    analyze(args.run_dir.resolve(), args.output.resolve(),
            args.sage_replay.resolve() if args.sage_replay else None)
