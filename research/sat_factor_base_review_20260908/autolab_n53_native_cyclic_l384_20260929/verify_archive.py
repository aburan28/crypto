#!/usr/bin/env python3
"""Rehash first-attempt evidence and independently replay full or partial outputs."""
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
HEX = re.compile(r"^[0-9a-f]{64}$")
SHA40 = re.compile(r"^[0-9a-f]{40}$")
CAPS = {"rho": 600, "training": 360, "base": 90, "rank": 180,
        "ic": 300, "recovery": 180, "independent_audit": 1200}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unpack_checked(archive: Path, destination: Path) -> int:
    total = 0
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        names = [member.name for member in members]
        assert len(names) == len(set(names)) and names[-1] == "SHA256SUMS"
        assert len(names) <= 160
        for member in members:
            name = Path(member.name)
            assert member.isfile() and not name.is_absolute() and ".." not in name.parts
            assert member.name == "SHA256SUMS" or name.parts[0] == "panel"
            assert member.size <= 128 * 1024**2
            total += member.size
            assert total <= 512 * 1024**2
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            stream = tar.extractfile(member)
            assert stream is not None
            target.write_bytes(stream.read())
    rows = (destination / "SHA256SUMS").read_text().splitlines()
    assert len(rows) == len(names) - 1
    expected = {}
    for row in rows:
        digest, name = row.split("  ", 1)
        assert HEX.fullmatch(digest) and name.startswith("panel/") and name not in expected
        expected[name] = digest
    assert set(expected) == set(names) - {"SHA256SUMS"}
    for name, digest in expected.items():
        assert sha(destination / name) == digest, name
    return len(expected)


def stage_location(key: str):
    if key == "rho":
        return Path("rho"), "rho", "rho"
    if key == "independent_audit":
        return Path("audit_driver"), "audit", "independent_audit"
    arm, stage = key.split("_", 1)
    assert arm in ("A", "B") and stage in ("training", "base", "rank", "ic", "recovery")
    folders = {"training": ("training", "producer"),
               "base": ("training/base_materialization", "base"),
               "rank": ("training/operational_rank", "rank"),
               "ic": ("ic", "ic"),
               "recovery": ("ic/operational_recovery", "recovery")}
    directory, basename = folders[stage]
    return Path(arm) / directory, basename, stage


def check_build(panel: Path, frozen: dict, summary: dict | None):
    path = panel / "build_receipt.json"
    if not path.is_file():
        assert not (panel / "build.stdout.txt").exists()
        assert not (panel / "build.stderr.txt").exists()
        return None
    build = json.loads(path.read_text())
    assert build["schema"] == "n53_native_cyclic_l384_build_v1"
    assert build["command"] == frozen["build_command"]
    for key, source_key in (("cargo_lock_sha256", "cargo_lock"),
                            ("cargo_toml_sha256", "cargo_toml"),
                            ("ic_source_sha256", "ic"),
                            ("rho_source_sha256", "rho")):
        assert build[key] == frozen["source_sha256"][source_key]
    for key, file in (("stdout_sha256", "build.stdout.txt"),
                      ("stderr_sha256", "build.stderr.txt")):
        if key in build:
            assert build[key] == sha(panel / file)
    if build["status"] == "SUCCESS":
        assert build["returncode"] == 0 and 0 <= build["wall_s"] <= 900
        assert build["toolchain"]["cargo"].startswith("cargo 1.93.1 ")
        assert build["toolchain"]["rustc"].startswith("rustc 1.93.1 ")
        assert set(build["binary_sha256"]) == {"ic", "rho"}
        assert all(HEX.fullmatch(value) for value in build["binary_sha256"].values())
    else:
        assert build["status"] in ("STARTED", "FAILED", "FAILED_OR_CENSORED")
        assert summary is None or summary["classification"] == "INVALID_RUNNER_EXCEPTION"
    if summary is not None and "build_receipt_sha256" in summary:
        assert build["status"] == "SUCCESS"
        assert summary["build_receipt_sha256"] == sha(path)
        assert summary["binary_sha256"] == build["binary_sha256"]
    return build


def expected_inputs(panel: Path, key: str, summary: dict, frozen: dict) -> dict[str, str]:
    sources = frozen["source_sha256"]
    binaries = summary["binary_sha256"]
    points = frozen["input_sha256"]["public_points"]
    if key == "rho":
        return {"rho_source": sources["rho"], "rho_binary": binaries["rho"],
                "points": points}
    if key == "independent_audit":
        result = {"audit_source": sources["audit"], "points": points,
                  "validator_scalars": frozen["input_sha256"]["validator_scalars"]}
        for arm in ("A", "B"):
            for label, path in {
                f"{arm}_training_raw": panel / arm / "training/producer.stdout.jsonl",
                f"{arm}_rank": panel / arm / "training/operational_solution.json",
                f"{arm}_ic_raw": panel / arm / "ic/ic.stdout.jsonl",
                f"{arm}_recovery": panel / arm / "ic/operational_recovery.json",
            }.items():
                if path.is_file():
                    result[label] = sha(path)
        rho_raw = panel / "rho/rho.stdout.jsonl"
        if rho_raw.is_file():
            result["rho_raw"] = sha(rho_raw)
        return result
    arm, stage = key.split("_", 1)
    training = panel / arm / "training"
    if stage == "training":
        return {"ic_source": sources["ic"], "ic_binary": binaries["ic"],
                "training_scalars": sha(training / "target_scalars.txt")}
    if stage == "base":
        return {"materializer_source": sources["materializer"],
                "training_raw": sha(training / "producer.stdout.jsonl")}
    if stage == "rank":
        return {"operational_source": sources["operational"],
                "producer_raw": sha(training / "producer.stdout.jsonl"),
                "base_header": sha(training / "base_header.jsonl"),
                "schedule": sha(training / "target_scalars.txt")}
    if stage == "ic":
        return {"ic_source": sources["ic"], "ic_binary": binaries["ic"],
                "base_header": sha(training / "base_header.jsonl"), "points": points}
    assert stage == "recovery"
    return {"operational_source": sources["operational"],
            "ic_raw": sha(panel / arm / "ic/ic.stdout.jsonl"), "points": points,
            "training_solution": sha(training / "operational_solution.json")}


def check_receipts(panel: Path, summary: dict, frozen: dict):
    ordered = ["rho"] + [f"{arm}_{stage}" for arm in ("B", "A")
                         for stage in ("training", "base", "rank", "ic", "recovery")] + ["independent_audit"]
    sequence = summary["stage_sequence"]
    assert isinstance(sequence, list) and len(sequence) == len(set(sequence))
    assert sequence and sequence[0] == "rho"
    exception = summary["classification"] == "INVALID_RUNNER_EXCEPTION"
    recorded_keys = set(summary["steps"])
    unreceipted = set(sequence) - recorded_keys
    assert recorded_keys <= set(sequence)
    assert not unreceipted or (exception and len(unreceipted) <= 2)
    if unreceipted:
        assert (panel / "failure.json").is_file()
    assert sequence == [key for key in ordered if key in sequence], "stage order changed"
    assert sequence[-1] == "independent_audit" or exception
    assert summary["stage_budget_seconds"] == {"operational": 3000, "independent_audit": 1200}
    prior_operational_ms = 0.0
    for key in sequence:
        directory, basename, stage = stage_location(key)
        expected_timeout = (1200.0 if key == "independent_audit" else
                            max(0.0, min(CAPS[stage], 3000.0 - prior_operational_ms / 1000)))
        if key not in summary["steps"]:
            assert exception
            continue
        recorded = summary["steps"][key]
        if "not_run" in recorded:
            assert key != "independent_audit" and expected_timeout == 0.0
            assert recorded == {"not_run": "GLOBAL_WALL_CAP"}
            continue
        root = panel / directory
        receipt = json.loads((root / "resource_receipt.json").read_text())
        manifest = json.loads((root / "manifest.json").read_text())
        assert receipt == recorded
        assert receipt["manifest_sha256"] == sha(root / "manifest.json")
        assert receipt["stdout_sha256"] == sha(root / f"{basename}.stdout.jsonl")
        assert receipt["stderr_sha256"] == sha(root / f"{basename}.stderr.txt")
        assert manifest["checkout_head"] == summary["checkout_head"]
        assert 0 < manifest["timeout_s"] <= CAPS[stage]
        assert math.isclose(manifest["timeout_s"], expected_timeout, rel_tol=0, abs_tol=1e-6)
        assert receipt["wall_ms"] >= 0
        if key != "independent_audit":
            prior_operational_ms += receipt["wall_ms"]
        assert receipt["user_cpu_s"] >= 0 and receipt["system_cpu_s"] >= 0
        assert receipt["peak_rss_bytes"] >= 0 and receipt["observed_group_peak_rss_bytes"] >= 0
        if summary["classification"].startswith("COMPLETE"):
            assert receipt["returncode"] == 0 and not receipt["timed_out"] and not receipt["rss_gate"]
            assert max(receipt["peak_rss_bytes"], receipt["observed_group_peak_rss_bytes"]) < 2 * 1024**3
        assert manifest["input_sha256"] == expected_inputs(panel, key, summary, frozen), key
        command = manifest["command"]
        env = manifest["environment"]
        assert "validator_scalars_L384.txt" not in json.dumps(command)
        assert "validator_scalars_L384.txt" not in json.dumps(env)
        if stage in ("training", "ic"):
            assert Path(command[0]).name == "koblitz_s5_sat_instance"
            assert command[1:] == ["53", "0", "1", "10", "natural", "1", "2000", "1", "internal"]
            assert env["KIC_ALGEBRA_ENCODING"] == "orbit_factorized"
            assert env["KIC_ORBIT_LAZY_RELATIVE_SUPPORT"] == "1"
            assert env["KIC_ORBIT_BRANCH_ORDER"] == "pair_then_pair"
            assert env["KIC_ORBIT_REP_ENCODING"] == "one_hot"
            assert env["KIC_ORBIT_BATCH_ONLY"] == "1"
            assert env["KIC_ORBIT_REGULAR_SCAN_POLICY"] == "target_cyclic_v1"
            assert env["RAYON_NUM_THREADS"] == "1"
            assert env["KIC_TASK_ID"].endswith("-" + key[0])
            if stage == "training":
                assert set(env) == {"KIC_ALGEBRA_ENCODING", "KIC_ORBIT_LAZY_RELATIVE_SUPPORT",
                                    "KIC_ORBIT_BRANCH_ORDER", "KIC_ORBIT_REP_ENCODING",
                                    "KIC_ORBIT_BATCH_ONLY", "KIC_ORBIT_INCLUDE_BASE_HEADER",
                                    "KIC_FACTOR_BASE_SELECTION", "KIC_ORBIT_REGULAR_SCAN_POLICY",
                                    "KIC_ORBIT_TARGET_SCALARS", "KIC_TASK_ID", "RAYON_NUM_THREADS"}
                assert env["KIC_FACTOR_BASE_SELECTION"] == "ascending_x_v1"
                assert env["KIC_ORBIT_INCLUDE_BASE_HEADER"] == "1"
                assert "KIC_FACTOR_BASE_JSONL" not in env
                assert Path(env["KIC_ORBIT_TARGET_SCALARS"]).name == "target_scalars.txt"
                assert manifest["input_sha256"]["training_scalars"] == frozen["input_sha256"][f"{key[0]}_scalars"]
            else:
                assert set(env) == {"KIC_ALGEBRA_ENCODING", "KIC_ORBIT_LAZY_RELATIVE_SUPPORT",
                                    "KIC_ORBIT_BRANCH_ORDER", "KIC_ORBIT_REP_ENCODING",
                                    "KIC_ORBIT_BATCH_ONLY", "KIC_FACTOR_BASE_JSONL",
                                    "KIC_ORBIT_REGULAR_SCAN_POLICY", "KIC_ORBIT_TARGET_POINTS_JSONL",
                                    "KIC_TASK_ID", "RAYON_NUM_THREADS"}
                assert Path(env["KIC_FACTOR_BASE_JSONL"]).name == "base_header.jsonl"
                assert Path(env["KIC_ORBIT_TARGET_POINTS_JSONL"]).name == "points_L384.jsonl"
        elif stage == "rho":
            assert Path(command[0]).name == "koblitz_rho_batch_ks"
            assert command[1:] == ["53", "0", "signed_frobenius", "384", "531929"]
            assert set(env) == {"KIC_RHO_TARGET_POINTS_JSONL", "KIC_RHO_BATCH_CORPUS",
                                "KIC_RHO_DP_BITS", "KIC_RHO_PRECOMPUTE_WALKS", "RAYON_NUM_THREADS"}
            assert Path(env["KIC_RHO_TARGET_POINTS_JSONL"]).name == "points_L384.jsonl"
            assert env["KIC_RHO_BATCH_CORPUS"] == "n53-native-cyclic-L384-20260929-v1"
            assert env["KIC_RHO_DP_BITS"] == "4" and env["KIC_RHO_PRECOMPUTE_WALKS"] == "0"
            assert env["RAYON_NUM_THREADS"] == "1"
        else:
            assert Path(command[0]).name.startswith("python")
            assert Path(command[1]).name == {"base": "cold_base.py", "rank": "operational.py",
                                                    "recovery": "operational.py", "independent_audit": "audit.py"}[stage]
            assert env == {}
        assert all(HEX.fullmatch(digest) for digest in manifest["input_sha256"].values())


def check_complete(panel: Path, summary: dict):
    assert summary["classification"] in ("COMPLETE_PRIMARY_B_WALL_WIN",
                                         "COMPLETE_PRIMARY_B_NO_CROSSOVER")
    assert set(summary["steps"]) == {"rho", "independent_audit"} | {
        f"{arm}_{stage}" for arm in ("A", "B")
        for stage in ("training", "base", "rank", "ic", "recovery")}
    assert all(summary["arms"][arm]["status"] == "OPERATIONAL_384_LOGS" for arm in ("A", "B"))
    steps = summary["steps"]
    rho = steps["rho"]["wall_ms"]
    costs = summary["costs"]
    assert math.isclose(costs["rho_operational_wall_ms"], rho)
    assert math.isclose(costs["rho_cpu_s"],
                        steps["rho"]["user_cpu_s"] + steps["rho"]["system_cpu_s"])
    for arm in ("A", "B"):
        keys = [f"{arm}_{stage}" for stage in ("training", "base", "rank", "ic", "recovery")]
        operational = sum(steps[key]["wall_ms"] for key in keys)
        lower = sum(steps[f"{arm}_{stage}"]["wall_ms"] for stage in ("training", "base", "ic"))
        c = costs[arm]
        assert math.isclose(c["ic_operational_wall_ms"], operational)
        assert math.isclose(c["ic_two_child_lower_wall_ms"], lower)
        assert math.isclose(c["ic_to_same_rho_ratio"], operational / rho)
        assert math.isclose(c["lower_to_same_rho_ratio"], lower / rho)
        assert math.isclose(c["ic_operational_cpu_s"],
                            sum(steps[key]["user_cpu_s"] + steps[key]["system_cpu_s"] for key in keys))
        assert c["classification"] == ("FIXED_STREAM_WALL_WIN" if operational < rho
                                       else "FIXED_STREAM_NO_CROSSOVER")
    assert summary["classification"] == (
        "COMPLETE_PRIMARY_B_WALL_WIN" if costs["B"]["ic_operational_wall_ms"] < rho
        else "COMPLETE_PRIMARY_B_NO_CROSSOVER")


def replay_audit(panel: Path, summary: dict):
    audit_receipt = summary["steps"].get("independent_audit", {})
    if (audit_receipt.get("returncode") != 0 or audit_receipt.get("timed_out")
            or audit_receipt.get("rss_gate") or not (panel / "audit.json").is_file()
            or "replay" not in summary):
        return "RAW_INTEGRITY_ONLY"
    archived = json.loads((panel / "audit.json").read_text())
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "audit.json"
        process = subprocess.run([sys.executable, str(HERE / "audit.py"),
                                  "--panel", str(panel), "--out", str(out)],
                                 cwd=REPO, capture_output=True, text=True, timeout=1240)
        assert process.returncode == 0, process.stderr[-3000:]
        assert json.loads(out.read_text()) == archived
    assert summary["replay"] == archived
    if summary["classification"].startswith("COMPLETE"):
        assert archived["classification"] == "INDEPENDENT_FULL_REPLAY"
    elif summary["classification"] == "PARTIAL_REPLAY_NO_END_TO_END_VERDICT":
        assert archived["classification"] == "INDEPENDENT_PARTIAL_REPLAY"
    return archived["classification"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle
    manifest = json.loads((bundle / "archive_manifest.json").read_text())
    archive = bundle / "evidence.tar.gz"
    assert sha(archive) == manifest["archive_sha256"]
    assert archive.stat().st_size == manifest["archive_bytes"]
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    import check_protocol
    assert frozen["protocol_sha256"] == sha(HERE / "PROTOCOL.md")
    assert all(sha(path) == frozen["source_sha256"][key]
               for key, path in check_protocol.sources().items())
    with tempfile.TemporaryDirectory() as tmp:
        unpacked = Path(tmp)
        assert unpack_checked(archive, unpacked) == manifest["files"]
        panel = unpacked / "panel"
        gate = json.loads((panel / "predispatch.json").read_text())
        assert gate["schema"] == "n53_native_cyclic_l384_dispatch_gate_v1"
        assert gate["status"] in ("ADMITTED", "REFUSED")
        assert gate["label"] == frozen["one_shot_pr_label"]
        summary_path = panel / "panel.json"
        summary = json.loads(summary_path.read_text()) if summary_path.is_file() else None
        build = check_build(panel, frozen, summary)
        if summary is None:
            assert manifest["classification"] == "NO_PANEL_SUMMARY"
            if gate["status"] == "REFUSED":
                assert build is None
            else:
                assert gate["run_attempt"] == "1"
            print(json.dumps({"verdict": "PASS_PREDISPATCH_OR_BUILD_RAW_INTEGRITY_ONLY",
                              "gate": gate["status"], "archive_sha256": sha(archive)}, sort_keys=True))
            return
        assert gate["status"] == "ADMITTED" and gate["run_attempt"] == "1"
        assert gate["reviewed_head"] == summary.get("checkout_head", gate["reviewed_head"])
        assert sha(summary_path) == manifest["panel_sha256"]
        if (bundle / "panel.json").is_file():
            assert (bundle / "panel.json").read_bytes() == summary_path.read_bytes()
        exception = summary["classification"] == "INVALID_RUNNER_EXCEPTION"
        if exception:
            assert summary["failure_sha256"] == sha(panel / "failure.json")
            assert manifest["classification"] == "INVALID_RUNNER_EXCEPTION"
            if "audit_failure_sha256" in summary:
                assert summary["audit_failure_sha256"] == sha(panel / "audit_failure.json")
            if "steps" not in summary:
                assert "stage_sequence" not in summary
                print(json.dumps({"verdict": "PASS_RUNNER_FAILURE_RAW_INTEGRITY_ONLY",
                                  "archive_sha256": sha(archive)}, sort_keys=True))
                return
        assert build is not None and build["status"] == "SUCCESS"
        assert summary["github_run_id"] == gate["run_id"]
        assert summary["github_run_attempt"] == gate["run_attempt"]
        assert gate["release_main_head"] == frozen["release_main_head"]
        assert summary["classification"] == manifest["classification"]
        assert summary["protocol_sha256"] == frozen["protocol_sha256"]
        assert summary["frozen_sha256"] == sha(HERE / "FROZEN.json")
        assert summary["points_sha256"] == frozen["input_sha256"]["public_points"]
        assert summary["source_sha256"] == frozen["source_sha256"]
        assert summary["input_sha256"] == frozen["input_sha256"]
        assert SHA40.fullmatch(summary["checkout_head"])
        assert summary["release_main_head"] == frozen["release_main_head"]
        assert SHA40.fullmatch(summary["dispatch_main_head"])
        subprocess.run(["git", "merge-base", "--is-ancestor",
                        summary["release_main_head"], summary["dispatch_main_head"]],
                       cwd=REPO, check=True)
        assert all(HEX.fullmatch(value) for value in summary["binary_sha256"].values())
        check_receipts(panel, summary, frozen)
        verdict = replay_audit(panel, summary)
        if summary["classification"].startswith("COMPLETE"):
            check_complete(panel, summary)
            assert verdict == "INDEPENDENT_FULL_REPLAY"
        else:
            assert summary["classification"].startswith(("PARTIAL", "INVALID", "CENSORED"))
        if exception:
            verdict_label = ("PASS_EXCEPTION_PREFIX_REPLAY" if verdict != "RAW_INTEGRITY_ONLY"
                             else "PASS_RUNNER_FAILURE_RAW_INTEGRITY_ONLY")
        else:
            verdict_label = ("PASS_COMPLETE_INDEPENDENT_REPLAY" if verdict == "INDEPENDENT_FULL_REPLAY"
                             else "PASS_PARTIAL_REPLAY" if verdict == "INDEPENDENT_PARTIAL_REPLAY"
                             else "PASS_RAW_INTEGRITY_ONLY")
        print(json.dumps({"verdict": verdict_label, "audit_classification": verdict,
                          "classification": summary["classification"],
                          "archive_sha256": sha(archive)}, sort_keys=True))


if __name__ == "__main__":
    main()
