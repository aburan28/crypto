#!/usr/bin/env python3
"""Rehash a durable L384 archive and rerun complete independent group replay."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
SHARED = HERE.parent / "autolab_shared_log_n53_20260925"
RHO_STUDY = HERE.parent / "autolab_matched_point_rho_n53_20260925"
HEX = re.compile(r"^[0-9a-f]{64}$")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    return {
        "ic": REPO / "examples/koblitz_s5_sat_instance.rs",
        "rho": REPO / "examples/koblitz_rho_batch_ks.rs",
        "operational": HERE / "operational.py",
        "audit": HERE / "audit.py",
        "runner": HERE / "run_panel.py",
        "archive_sealer": HERE / "archive.py",
        "archive_verifier": HERE / "verify_archive.py",
        "training_schedule_reference": ORBIT / "cold_batch_rank.py",
        "independent_group_reference": ORBIT / "independent_replay_20260924_codex/replay.py",
        "rho_audit_math": RHO_STUDY / "analyze.py",
        "point_generation_math": SHARED / "generate_targets.py",
    }


def unpack_checked(archive: Path, destination: Path) -> int:
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        names = [member.name for member in members]
        assert len(names) == len(set(names))
        assert names[-1] == "SHA256SUMS"
        assert all(member.isfile() for member in members)
        for member in members:
            name = Path(member.name)
            assert not name.is_absolute() and ".." not in name.parts
            assert member.name == "SHA256SUMS" or name.parts[0] == "panel"
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            extracted = tar.extractfile(member)
            assert extracted is not None
            target.write_bytes(extracted.read())
    checks = (destination / "SHA256SUMS").read_text().splitlines()
    assert len(checks) == len(names) - 1
    expected = {}
    for line in checks:
        digest, name = line.split("  ", 1)
        assert HEX.fullmatch(digest) and name.startswith("panel/")
        assert name not in expected
        expected[name] = digest
    assert set(expected) == set(names) - {"SHA256SUMS"}
    for name, digest in expected.items():
        assert sha(destination / name) == digest, name
    return len(expected)


def check_receipts(panel: Path, summary: dict):
    locations = {
        "training_producer": ("training", "producer"),
        "operational_rank": ("training/operational_rank", "rank"),
        "rho": ("rho", "rho"),
        "ic": ("ic", "ic"),
        "operational_recovery": ("ic/operational_recovery", "recovery"),
        "independent_audit": ("audit_driver", "audit"),
    }
    for key, receipt in summary["steps"].items():
        directory, basename = locations[key]
        root = panel / directory
        recorded = json.loads((root / "resource_receipt.json").read_text())
        assert receipt == recorded, key
        assert receipt["manifest_sha256"] == sha(root / "manifest.json"), key
        assert receipt["stdout_sha256"] == sha(root / f"{basename}.stdout.jsonl"), key
        assert receipt["stderr_sha256"] == sha(root / f"{basename}.stderr.txt"), key
        manifest = json.loads((root / "manifest.json").read_text())
        assert HEX.fullmatch(manifest["checkout_head"])
        assert all(HEX.fullmatch(value) for value in manifest["input_sha256"].values())
        assert receipt["wall_ms"] >= 0
        assert receipt["user_cpu_s"] >= 0 and receipt["system_cpu_s"] >= 0
        assert receipt["peak_rss_bytes"] >= 0 and receipt["observed_group_peak_rss_bytes"] >= 0


def check_complete(panel: Path, summary: dict):
    assert summary["classification"] in (
        "COMPLETE_FIXED_STREAM_OPERATIONAL_WALL_WIN",
        "COMPLETE_FIXED_STREAM_NO_CROSSOVER",
        "COMPLETE_FIXED_STREAM_INCONCLUSIVE",
    )
    assert set(summary["steps"]) == {
        "training_producer", "operational_rank", "rho", "ic", "operational_recovery", "independent_audit"
    }
    for receipt in summary["steps"].values():
        assert receipt["returncode"] == 0
        assert not receipt["timed_out"] and not receipt["rss_gate"]
    replay_file = panel / "rechecked_audit.json"
    process = subprocess.run(
        [sys.executable, str(HERE / "audit.py"), "--panel", str(panel), "--out", str(replay_file)],
        cwd=REPO, capture_output=True, text=True, timeout=620,
    )
    assert process.returncode == 0, process.stderr[-2000:]
    assert json.loads(replay_file.read_text()) == json.loads((panel / "audit.json").read_text())
    costs = summary["costs"]
    steps = summary["steps"]
    expected_operational = sum(steps[name]["wall_ms"] for name in (
        "training_producer", "operational_rank", "ic", "operational_recovery"
    ))
    expected_lower = steps["training_producer"]["wall_ms"] + steps["ic"]["wall_ms"]
    rho = steps["rho"]["wall_ms"]
    audit = steps["independent_audit"]["wall_ms"]
    assert math.isclose(costs["ic_operational_wall_ms"], expected_operational)
    assert math.isclose(costs["ic_two_child_lower_wall_ms"], expected_lower)
    assert math.isclose(costs["ic_audit_inclusive_wall_ms"], expected_operational + audit)
    assert math.isclose(costs["rho_operational_wall_ms"], rho)
    assert math.isclose(costs["operational_ratio_to_rho"], expected_operational / rho)
    assert math.isclose(costs["two_child_lower_ratio_to_rho"], expected_lower / rho)
    assert math.isclose(costs["audit_inclusive_ratio_to_rho"], (expected_operational + audit) / rho)
    expected_class = (
        "COMPLETE_FIXED_STREAM_OPERATIONAL_WALL_WIN" if expected_operational < rho
        else "COMPLETE_FIXED_STREAM_NO_CROSSOVER" if expected_lower > rho
        else "COMPLETE_FIXED_STREAM_INCONCLUSIVE"
    )
    assert summary["classification"] == expected_class
    return json.loads(replay_file.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle
    manifest = json.loads((bundle / "archive_manifest.json").read_text())
    archive = bundle / "evidence.tar.gz"
    assert sha(archive) == manifest["archive_sha256"]
    assert archive.stat().st_size == manifest["archive_bytes"]
    freeze = json.loads((HERE / "SOURCE_FREEZE.json").read_text())
    assert freeze["target_sha256"] == sha(HERE / "points_L384.jsonl")
    assert freeze["protocol_sha256"] == sha(HERE / "PROTOCOL.md")
    assert all(sha(path) == freeze["source_sha256"][key] for key, path in sources().items())
    with tempfile.TemporaryDirectory() as temp:
        unpacked = Path(temp)
        count = unpack_checked(archive, unpacked)
        assert count == manifest["files"]
        panel = unpacked / "panel"
        summary = json.loads((panel / "panel.json").read_text())
        assert sha(panel / "panel.json") == manifest["panel_sha256"]
        if (bundle / "panel.json").exists():
            assert (bundle / "panel.json").read_bytes() == (panel / "panel.json").read_bytes()
        assert summary["classification"] == manifest["classification"]
        assert summary["protocol_sha256"] == freeze["protocol_sha256"]
        assert summary["points_sha256"] == freeze["target_sha256"]
        assert summary["source_sha256"] == freeze["source_sha256"]
        assert summary["base_gzip_sha256"] == freeze["base_gzip_sha256"]
        assert summary["validator_manifest_sha256"] == freeze["validator_manifest_sha256"]
        assert HEX.fullmatch(summary["checkout_head"])
        assert all(HEX.fullmatch(value) for value in summary["binary_sha256"].values())
        check_receipts(panel, summary)
        if summary["classification"].startswith("COMPLETE"):
            replay = check_complete(panel, summary)
            print(json.dumps({
                "verdict": "PASS_COMPLETE_INDEPENDENT_REPLAY",
                "training_relations": replay["training"]["training_relations_replayed"],
                "ic_logs": replay["ic"]["ic_relations_replayed"],
                "rho_logs": replay["rho"]["rho_scalars_replayed"],
                "archive_sha256": manifest["archive_sha256"],
            }, sort_keys=True))
        else:
            assert summary["classification"].startswith(("CENSORED", "INVALID"))
            print(json.dumps({
                "verdict": "PASS_CENSORED_OR_INVALID_RAW_INTEGRITY_ONLY",
                "classification": summary["classification"],
                "archive_sha256": manifest["archive_sha256"],
            }, sort_keys=True))


if __name__ == "__main__":
    main()
