"""Resume the N53 panel with an explicit R2 retry and preserved failure."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def original_runner():
    spec = importlib.util.spec_from_file_location("n53_original_run", PARENT / "run.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def preflight(binaries: Path, probe: Path, partial_dir: Path, output: Path):
    freeze = json.loads((HERE / "FREEZE.json").read_bytes())
    assert freeze["kind"] == "s3_n53_c1200_c1800_pmullfold_confirmation34_recovery_v1"
    assert freeze["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert freeze["analyzer_sha256"] == sha((HERE / "analyze.py").read_bytes())
    assert freeze["replay_sha256"] == sha((HERE / "replay_sage.py").read_bytes())
    assert freeze["original_freeze_sha256"] == sha((PARENT / "FREEZE.json").read_bytes())
    assert freeze["original_runner_sha256"] == sha((PARENT / "run.py").read_bytes())
    assert freeze["original_partial_analysis_sha256"] == sha(
        (PARENT / "local/partial_analysis.json").read_bytes())
    assert freeze["original_partial_archive_sha256"] == sha(
        (PARENT / "local/partial_raw.tar.zst").read_bytes())
    partial_raw = (partial_dir / "raw.jsonl").read_bytes()
    assert freeze["original_partial_ledger_sha256"] == sha(partial_raw)
    rows = [json.loads(line) for line in partial_raw.splitlines()]
    assert len(rows) == freeze["original_verified_prefix_count"] + 1 == 28
    assert all(row["verified_native_and_fixture"] for row in rows[:-1])
    failed = rows[-1]
    assert failed["run_id"] == freeze["original_failure"]["run_id"]
    assert failed["error"] == freeze["original_failure"]["error"]
    assert sha((partial_dir / "T220_baseline_R1.json").read_bytes()) == (
        freeze["original_failure"]["report_sha256"])
    original = original_runner()
    original_freeze = original.preflight(binaries, probe)
    assert len(freeze["run_order"]) == freeze["planned_run_count"] == 68
    assert len(freeze["recovery_run_order"]) == freeze["recovery_process_count"] == 41
    assert freeze["run_order"][:27] == original_freeze["run_order"][:27]
    assert freeze["run_order"][28:] == original_freeze["run_order"][28:]
    assert freeze["run_order"][27] == {
        "label": "T220", "arm": "baseline", "repetition": 2}
    assert not output.exists()
    output.parent.mkdir(parents=True, exist_ok=True)
    stat = shutil.disk_usage(output.parent)
    assert stat.free >= freeze["minimum_output_filesystem_free_bytes"]
    return freeze, original_freeze, original, partial_raw


def run(binaries: Path, probe: Path, partial_dir: Path, output: Path) -> None:
    freeze, original_freeze, original, partial_raw = preflight(
        binaries, probe, partial_dir, output)
    output.mkdir(parents=True, exist_ok=False)
    prefix = freeze["original_verified_prefix_count"]
    for entry in freeze["run_order"][:prefix]:
        stem = f"{entry['label']}_{entry['arm']}_R{entry['repetition']}"
        for suffix in (".json", ".stdout", ".stderr"):
            shutil.copyfile(partial_dir / (stem + suffix), output / (stem + suffix))
    ledger = output / "raw.jsonl"
    ledger.write_bytes(b"".join(partial_raw.splitlines(keepends=True)[:prefix]))
    targets = {target["label"]: target for target in original_freeze["selected_targets"]}
    for entry in freeze["recovery_run_order"]:
        row = original.run_one(original_freeze, binaries, output,
                               targets[entry["label"]], entry["arm"],
                               entry["repetition"])
        assert (row["label"], row["arm"], row["repetition"]) == (
            entry["label"], entry["arm"], entry["repetition"])
        with ledger.open("a") as handle:
            handle.write(json.dumps(row, separators=(",", ":")) + "\n")
        print(json.dumps({"label": row["label"], "arm": row["arm"],
                          "repetition": row["repetition"],
                          "verified": row["verified_native_and_fixture"],
                          "online_ms": row["online_ms"],
                          "peak_rss_bytes": row["peak_rss_bytes"],
                          "error": row["error"]}), flush=True)
        if not row["verified_native_and_fixture"]:
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--partial-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binaries.resolve(), args.probe.resolve(),
        args.partial_dir.resolve(), args.output.resolve())
