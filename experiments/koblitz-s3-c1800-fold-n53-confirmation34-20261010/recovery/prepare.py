"""Freeze an explicit N53 recovery after the T220 baseline report-write failure."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def prepare(partial_dir: Path) -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "Recovery freeze is immutable"
    original_path = PARENT / "FREEZE.json"
    original = json.loads(original_path.read_bytes())
    partial_raw = (partial_dir / "raw.jsonl").read_bytes()
    rows = [json.loads(line) for line in partial_raw.splitlines()]
    partial_analysis_path = PARENT / "local/partial_analysis.json"
    partial_analysis = json.loads(partial_analysis_path.read_bytes())
    assert len(rows) == 28 and partial_analysis["verified_runs"] == 27
    assert partial_analysis["complete"] is False
    assert partial_analysis["ledger_sha256"] == sha(partial_raw)
    assert [(row["label"], row["arm"], row["repetition"]) for row in rows] == [
        (entry["label"], entry["arm"], entry["repetition"])
        for entry in original["run_order"][:28]]
    assert all(row["verified_native_and_fixture"] for row in rows[:27])
    failed = rows[27]
    assert (failed["label"], failed["arm"], failed["repetition"]) == (
        "T220", "baseline", 1)
    assert failed["exit_code"] == 101 and failed["error"] == "invalid_json"
    assert failed["verified_native_and_fixture"] is False
    failed_report = partial_dir / "T220_baseline_R1.json"
    assert sha(failed_report.read_bytes()) == failed["report_sha256"]
    assert len(failed_report.read_bytes()) == 122880
    assert (partial_dir / "T220_baseline_R1.stdout").read_bytes() == b""
    assert (partial_dir / "T220_baseline_R1.stderr").read_bytes() == b""
    partial_archive = PARENT / "local/partial_raw.tar.zst"
    assert partial_archive.exists()
    recovery = copy.deepcopy(original)
    for record in recovery["arms"].values():
        manifest = (PARENT / record["candidate_manifest_relative"]).resolve()
        record["candidate_manifest_relative"] = os.path.relpath(manifest, HERE)
    smoke = (PARENT / recovery["smoke_archive_relative"]).resolve()
    recovery["smoke_archive_relative"] = os.path.relpath(smoke, HERE)
    new_order = copy.deepcopy(original["run_order"])
    new_order[27]["repetition"] = 2
    recovery["run_order"] = new_order
    recovery["kind"] = "s3_n53_c1200_c1800_pmullfold_confirmation34_recovery_v1"
    recovery["recovery_run_order"] = new_order[27:]
    recovery["recovery_process_count"] = len(new_order[27:])
    recovery["original_verified_prefix_count"] = 27
    recovery["original_freeze_sha256"] = sha(original_path.read_bytes())
    recovery["original_runner_sha256"] = sha((PARENT / "run.py").read_bytes())
    recovery["original_partial_ledger_sha256"] = sha(partial_raw)
    recovery["original_partial_analysis_sha256"] = sha(partial_analysis_path.read_bytes())
    recovery["original_partial_archive_sha256"] = sha(partial_archive.read_bytes())
    recovery["original_failure"] = {
        "run_id": failed["run_id"], "report_sha256": failed["report_sha256"],
        "exit_code": failed["exit_code"], "error": failed["error"],
        "report_bytes": len(failed_report.read_bytes()),
    }
    recovery["minimum_output_filesystem_free_bytes"] = 2 * 1024**3
    recovery["runner_sha256"] = sha((HERE / "run.py").read_bytes())
    recovery["analyzer_sha256"] = sha((HERE / "analyze.py").read_bytes())
    recovery["replay_sha256"] = sha((HERE / "replay_sage.py").read_bytes())
    recovery["preparer_sha256"] = sha(Path(__file__).read_bytes())
    assert len(new_order) == recovery["planned_run_count"] == 68
    assert recovery["recovery_process_count"] == 41
    freeze_path.write_text(json.dumps(recovery, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "verified_prefix": 27, "recovery_runs": 41,
                      "failed_original_run_id": failed["run_id"]}))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partial-dir", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.partial_dir.resolve())
