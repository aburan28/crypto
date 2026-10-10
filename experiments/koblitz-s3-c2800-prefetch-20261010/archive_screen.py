"""Preserve the complete paired screen with deterministic compressed evidence."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


HERE = Path(__file__).resolve().parent
FREEZE = HERE / "SCREEN_FREEZE.json"
ANALYSIS = HERE / "analysis.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def archive(run_dir: Path) -> None:
    output = HERE / "local/paired_screen"
    if output.exists():
        raise FileExistsError(output)
    plan = json.loads(FREEZE.read_bytes())
    analysis = json.loads(ANALYSIS.read_bytes())
    ledger = (run_dir / "raw.jsonl").read_bytes()
    assert analysis["freeze_sha256"] == sha(FREEZE.read_bytes())
    assert analysis["ledger_sha256"] == sha(ledger)
    assert analysis["complete"]
    assert analysis["run_count"] == analysis["verified_and_memory_ok_count"] == 16
    assert analysis["target_count"] == len(plan["target_labels"]) == 8
    assert analysis["same_first_witness_all"] is True
    assert not analysis["failed_runs"]
    assert all(not differences for differences in
               analysis["unexpected_report_differences"].values())
    reports = sorted(run_dir.glob("T*_R1.json"))
    assert len(reports) == 16
    output.mkdir(parents=True)
    ledger_gzip = output / "raw.jsonl.gz"
    ledger_gzip.write_bytes(gzip.compress(ledger, mtime=0))
    reports_archive = output / "raw_reports.tar.gz"
    with reports_archive.open("wb") as raw_file:
        with gzip.GzipFile(filename="", fileobj=raw_file, mode="wb", mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w|") as tar:
                for path in reports:
                    info = tar.gettarinfo(str(path), arcname=path.name)
                    info.mtime = info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    with path.open("rb") as source:
                        tar.addfile(info, source)
    shutil.copyfile(run_dir / "host.json", output / "host.json")
    receipt = {
        "kind": "s3_n53_c2800_prefetch_paired_screen_archive_v1",
        "freeze_sha256": sha(FREEZE.read_bytes()),
        "analysis_sha256": sha(ANALYSIS.read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze_screen.py").read_bytes()),
        "ledger_decompressed_sha256": sha(ledger),
        "ledger_gzip_sha256": sha(ledger_gzip.read_bytes()),
        "raw_reports_archive_sha256": sha(reports_archive.read_bytes()),
        "host_sha256": sha((output / "host.json").read_bytes()),
        "report_count": len(reports),
        "run_count": analysis["run_count"],
        "controlled_speedup": None,
    }
    (output / "archive_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    archive(args.run_dir)
