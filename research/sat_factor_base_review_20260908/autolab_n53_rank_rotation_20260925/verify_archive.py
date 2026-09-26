#!/usr/bin/env python3
"""Verify raw provenance and replay a completed n53 rotation rank result."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
HEX64 = re.compile(r"^[0-9a-f]{64}$")
HEX40 = re.compile(r"^[0-9a-f]{40}$")
ARMS = (
    ("native_lex", "native", "lex"),
    ("certified_cyclic", "certified", "target_cyclic_v1"),
    ("certified_lex", "certified", "lex"),
    ("native_cyclic", "native", "target_cyclic_v1"),
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha(path: Path) -> str:
    return sha(path.read_bytes())


def verify_current_freeze():
    import check_protocol
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["schema"] == "n53_target_cyclic_rank_factorial_freeze_v1"
    assert frozen["status"] == "released_for_one_outcome"
    assert HEX40.fullmatch(frozen["release_main_head"])
    assert file_sha(HERE / "PROTOCOL.md") == frozen["protocol_sha256"]
    for name, path in check_protocol.sources().items():
        assert file_sha(path) == frozen["source_sha256"][name], name
    for name, path in check_protocol.inputs().items():
        assert file_sha(path) == frozen["input_sha256"][name], name
    return frozen


def unpack_checked(raw: bytes, destination: Path, manifest: dict) -> int:
    with tarfile.open(fileobj=io.BytesIO(gzip.decompress(raw)), mode="r:") as tar:
        infos = tar.getmembers()
        names = [info.name for info in infos]
        assert len(names) == len(set(names)) and names[-1] == "SHA256SUMS"
        assert all(info.isfile() for info in infos)
        for info in infos:
            name = Path(info.name)
            assert not name.is_absolute() and ".." not in name.parts
            assert info.name == "SHA256SUMS" or name.parts[0] == "panel"
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            payload = tar.extractfile(info)
            assert payload is not None
            path.write_bytes(payload.read())
    expected = {}
    for line in (destination / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert HEX64.fullmatch(digest) and name.startswith("panel/")
        assert name not in expected
        expected[name] = digest
    assert set(expected) == set(names) - {"SHA256SUMS"}
    assert len(expected) == manifest["files"]
    for name, digest in expected.items():
        assert file_sha(destination / name) == digest, name
    return len(expected)


def check_stage(panel: Path, summary: dict, frozen: dict, stage: str, name: str | None,
                base: str | None, policy: str | None):
    is_arm = name is not None
    root = panel / (name if is_arm else "audit_stage")
    basename = "producer" if is_arm else "audit"
    receipt = summary["stage_receipts"][stage]
    recorded = json.loads((root / "resource_receipt.json").read_text())
    assert receipt == recorded, stage
    manifest_path = root / "manifest.json"
    stdout = root / f"{basename}.stdout.jsonl"
    stderr = root / f"{basename}.stderr.txt"
    assert receipt["manifest_sha256"] == file_sha(manifest_path), stage
    assert receipt["stdout_sha256"] == file_sha(stdout), stage
    assert receipt["stderr_sha256"] == file_sha(stderr), stage
    assert isinstance(receipt["returncode"], int)
    assert isinstance(receipt["timed_out"], bool)
    assert isinstance(receipt["rss_gate"], bool)
    assert receipt["wall_ms"] >= 0
    assert receipt["user_cpu_s"] >= 0 and receipt["system_cpu_s"] >= 0
    assert receipt["peak_rss_bytes"] >= 0 and receipt["observed_group_peak_rss_bytes"] >= 0
    manifest = json.loads(manifest_path.read_text())
    assert manifest["checkout_head"] == summary["checkout_head"]
    assert manifest["host"] == summary["host"] and manifest["machine"] == summary["machine"]
    assert 0 < manifest["timeout_s"] <= frozen["caps_seconds"]["per_arm" if is_arm else "audit"]
    source = frozen["source_sha256"]
    inputs = frozen["input_sha256"]
    expected = {
        "targets": inputs["target_scalars"],
        "points_validator": inputs["target_points"],
        "producer_source": source["rust_producer"],
        "producer_binary": summary["binary_sha256"],
        "cargo_lock": source["cargo_lock"],
        "protocol": frozen["protocol_sha256"],
    }
    if is_arm:
        assert receipt["rayon_num_threads"] == summary["rayon_num_threads"] == 1
        assert manifest["command"][1:] == ["53", "0", "1", "10", "natural", "1", "2000", "1", "internal"]
        assert Path(manifest["command"][0]).name == "koblitz_s5_sat_instance"
        env = manifest["environment"]
        expected_env_keys = {
            "KIC_ALGEBRA_ENCODING", "KIC_ORBIT_LAZY_RELATIVE_SUPPORT",
            "KIC_ORBIT_BRANCH_ORDER", "KIC_ORBIT_REP_ENCODING",
            "KIC_ORBIT_BATCH_ONLY", "KIC_ORBIT_INCLUDE_BASE_HEADER",
            "KIC_ORBIT_TARGET_SCALARS", "KIC_TASK_ID",
            "KIC_ORBIT_REGULAR_SCAN_POLICY",
        }
        if base == "certified":
            expected_env_keys.add("KIC_FACTOR_BASE_JSONL")
        assert set(env) == expected_env_keys
        assert env["KIC_ORBIT_REGULAR_SCAN_POLICY"] == policy
        assert env["KIC_ALGEBRA_ENCODING"] == "orbit_factorized"
        assert env["KIC_ORBIT_LAZY_RELATIVE_SUPPORT"] == "1"
        assert env["KIC_ORBIT_BRANCH_ORDER"] == "pair_then_pair"
        assert env["KIC_ORBIT_REP_ENCODING"] == "one_hot"
        assert env["KIC_ORBIT_BATCH_ONLY"] == "1"
        assert env["KIC_ORBIT_INCLUDE_BASE_HEADER"] == "1"
        assert env["KIC_TASK_ID"] == "TASK-IC-N53-TARGET-CYCLIC-RANK-20260925"
        assert Path(env["KIC_ORBIT_TARGET_SCALARS"]).name == "target_scalars.txt"
        if base == "certified":
            assert Path(env["KIC_FACTOR_BASE_JSONL"]).name == "certified_base.jsonl"
            expected["certified_base"] = file_sha(panel / "certified_base.jsonl")
        else:
            assert "KIC_FACTOR_BASE_JSONL" not in env
        assert manifest["input_sha256"] == expected, stage
    else:
        assert Path(manifest["command"][0]).name.startswith("python")
        assert Path(manifest["command"][1]).name == "audit.py"
        assert manifest["command"][2] == "--out"
        expected = {arm: file_sha(panel / arm / "producer.stdout.jsonl") for arm, _, _ in ARMS}
        expected |= {
            "targets": inputs["target_scalars"],
            "points": inputs["target_points"],
            "audit_source": source["audit"],
            "independent_group_source": source["independent_group"],
            "independent_rank_source": source["independent_rank"],
            "producer_binary": summary["binary_sha256"],
            "protocol": frozen["protocol_sha256"],
        }
        assert manifest["input_sha256"] == expected


def verify(bundle: Path) -> dict:
    frozen = verify_current_freeze()
    archive = bundle / "evidence.tar.gz"
    outer = json.loads((bundle / "archive_manifest.json").read_text())
    raw = archive.read_bytes()
    assert sha(raw) == outer["archive_sha256"]
    assert len(raw) == outer["archive_bytes"]
    with tempfile.TemporaryDirectory(prefix="n53-rotation-replay-") as scratch:
        unpacked = Path(scratch)
        count = unpack_checked(raw, unpacked, outer)
        panel = unpacked / "panel"
        summary_path = panel / "panel.json"
        summary = json.loads(summary_path.read_text())
        assert file_sha(summary_path) == outer["panel_sha256"]
        assert (bundle / "panel.json").read_bytes() == summary_path.read_bytes()
        assert summary["schema"] == "n53_target_cyclic_rank_factorial_panel_v1"
        assert summary["classification"] == outer["classification"]
        assert summary["protocol_sha256"] == frozen["protocol_sha256"]
        assert summary["source_sha256"] == frozen["source_sha256"]["rust_producer"]
        assert summary["cargo_lock_sha256"] == frozen["source_sha256"]["cargo_lock"]
        assert summary["target_scalars_sha256"] == frozen["input_sha256"]["target_scalars"]
        assert summary["target_points_sha256"] == frozen["input_sha256"]["target_points"]
        assert summary["certified_base_gzip_sha256"] == frozen["input_sha256"]["certified_base_gzip"]
        base_path = panel / "certified_base.jsonl"
        if "certified_base_materialized_sha256" in summary:
            base_bytes = gzip.decompress(check_base_path(frozen).read_bytes())
            assert base_path.read_bytes() == base_bytes
            assert summary["certified_base_materialized_sha256"] == file_sha(base_path)
        else:
            # A terminal wrapper/setup failure can occur before the certified
            # fixture is materialized. No base or rank claim is admitted.
            assert summary["classification"] == "RUNNING" or summary["classification"].startswith(
                "CENSORED_OR_INVALID_WRAPPER_"
            )
            assert summary.get("active_stage") in ("setup", "base_materialization", "unknown")
        assert HEX40.fullmatch(summary["checkout_head"])
        assert HEX40.fullmatch(summary["dispatch_main_head"])
        subprocess.run(["git", "merge-base", "--is-ancestor", frozen["release_main_head"],
                        summary["dispatch_main_head"]], cwd=REPO, check=True)
        assert HEX64.fullmatch(summary["binary_sha256"])
        assert summary["rayon_num_threads"] == 1
        assert summary["process_group_rss_cap_bytes"] == frozen["process_group_rss_cap_bytes"]
        assert summary["arm_wall_cap_s"] == frozen["caps_seconds"]["per_arm"]
        assert summary["audit_wall_cap_s"] == frozen["caps_seconds"]["audit"]
        assert summary["global_wall_cap_s"] == frozen["caps_seconds"]["global"]
        assert summary["attack_speed_crossover"] is None and summary["common_operation_unit"] is None
        subprocess.run(["git", "merge-base", "--is-ancestor", frozen["parent_pr_815_head"],
                        summary["checkout_head"]], cwd=REPO, check=True)
        present = set(summary["stage_receipts"])
        arm_names = [name for name, _, _ in ARMS]
        # panel.json uses sort_keys=True, so dictionary insertion order is
        # not the process order. The completed-stage set must be a prefix.
        assert (present == set(arm_names[:len(present)])
                or present == set(arm_names + ["audit"]))
        for name, base, policy in ARMS:
            if name in summary["stage_receipts"]:
                check_stage(panel, summary, frozen, name, name, base, policy)
        if "audit" in summary["stage_receipts"]:
            check_stage(panel, summary, frozen, "audit", None, None, None)
        classification = summary["classification"]
        wrapper_failure = classification.startswith("CENSORED_OR_INVALID_WRAPPER_")
        if wrapper_failure:
            exception_path = panel / "exception_receipt.json"
            assert summary["exception_receipt_sha256"] == file_sha(exception_path)
            exception = json.loads(exception_path.read_text())
            assert exception["stage"] == summary["active_stage"]
            assert exception["kind"] and exception["traceback"]
            assert exception["completed_stage_keys"] == sorted(summary["stage_receipts"])
        elif classification == "RUNNING":
            # An abrupt OS/workflow termination may leave the last RUNNING
            # snapshot. Preserve its raw bytes as incomplete/invalid, never
            # as a rank or process-cost result.
            assert "exception_receipt_sha256" not in summary
        status = "RAW_PROVENANCE_PASS"
        if classification == "RANK_STAGE_REPLAYED":
            assert present == set(arm_names + ["audit"])
            assert all(receipt["returncode"] == 0 and not receipt["timed_out"]
                       and not receipt["rss_gate"] and
                       receipt["peak_rss_bytes"] < frozen["process_group_rss_cap_bytes"] and
                       receipt["observed_group_peak_rss_bytes"] < frozen["process_group_rss_cap_bytes"]
                       for receipt in summary["stage_receipts"].values())
            expected_audit = (panel / "audit.json").read_bytes()
            assert summary["audit_sha256"] == sha(expected_audit)
            assert summary["measured_arm_wall_ms"] == {
                name: summary["stage_receipts"][name]["wall_ms"] for name in arm_names
            }
            process = subprocess.run([sys.executable, str(HERE / "audit.py"), "--out", str(panel)],
                                     cwd=REPO, capture_output=True, text=True, timeout=620)
            assert process.returncode == 0, process.stderr[-2000:]
            assert (panel / "audit.json").read_bytes() == expected_audit
            status = "RAW_PROVENANCE_AND_INDEPENDENT_RANK_REPLAY_PASS"
        elif wrapper_failure:
            status = "RAW_PROVENANCE_INVALID_WRAPPER_ONLY"
        elif classification == "RUNNING":
            status = "RAW_PROVENANCE_INCOMPLETE_RUNNER_ONLY"
        else:
            assert classification.startswith(("CENSORED", "INVALID"))
        return {"status": status, "classification": classification,
                "files": count, "archive_sha256": outer["archive_sha256"]}


def check_base_path(frozen: dict) -> Path:
    import check_protocol
    path = check_protocol.inputs()["certified_base_gzip"]
    assert file_sha(path) == frozen["input_sha256"]["certified_base_gzip"]
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.bundle), sort_keys=True))
