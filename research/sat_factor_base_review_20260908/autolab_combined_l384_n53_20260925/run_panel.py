#!/usr/bin/env python3
"""Cold n53 L384 one-index/one-table matched panel with split operational work."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
SHARED = HERE.parent / "autolab_shared_log_n53_20260925"
RHO_STUDY = HERE.parent / "autolab_matched_point_rho_n53_20260925"
BASE_GZ = ORBIT / "independent_replay_20260924_codex/base_header.jsonl.gz"
IC_SOURCE = REPO / "examples/koblitz_s5_sat_instance.rs"
RHO_SOURCE = REPO / "examples/koblitz_rho_batch_ks.rs"
IC_EXE = REPO / "target/release/examples/koblitz_s5_sat_instance"
RHO_EXE = REPO / "target/release/examples/koblitz_rho_batch_ks"
POINTS = HERE / "points_L384.jsonl"
BASE_HASH = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"
POINTS_HASH = "d5185187014a12516aeef29306b65d4d20864293bb2c60667a9684fa8e51ac97"
MANIFEST_HASH = "f1843670a169d65645bf83886aeffee2e60e25627faa0002763eb1e7c13d8362"
IC_SOURCE_HASH = "c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09"
RHO_SOURCE_HASH = "fedacb54e441979c8c32860e7b5639e43799234741d49677010164564536d2c8"
BASE_GZ_HASH = "23397af2ef668aed0775bcb409e1ae19555357ded635452818c9a3812f679d08"
SOURCE_POINTS_HASHES = (
    "1ec6ca63fb81d6e4e07ec9d2a8e5389d02f25fb73f0fc56e3d5ccab85a59dd6f",
    "d4cd94396d04d755e4eb6a3f4f550da329769ee6b42a1020fc95615c80c9df55",
    "1b3d3a0d9ed59f6f5f96ec0eff80776ee00f75a47ce07e4abfeb95cbf9d61a22",
)
MAX_RSS = 2 * 1024**3
MAX_WALL = 2700


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean_env() -> dict[str, str]:
    return {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}


def linux_process_group_rss(pgid: int) -> int:
    if not Path("/proc").is_dir():
        return 0
    total = 0
    for entry in Path("/proc").iterdir():
        if not entry.name.isdecimal():
            continue
        try:
            stat = (entry / "stat").read_text()
            fields = stat[stat.rfind(")") + 2:].split()
            if int(fields[2]) != pgid:
                continue
            for line in (entry / "status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    total += int(line.split()[1]) * 1024
                    break
        except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
            continue
    return total


def measure(command: list[str], env: dict[str, str], directory: Path,
            basename: str, timeout: float, input_files: dict[str, Path]) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    stdout = directory / f"{basename}.stdout.jsonl"
    stderr = directory / f"{basename}.stderr.txt"
    assert not stdout.exists() and not stderr.exists(), "measurement output already exists"
    manifest = {
        "command": command,
        "environment": {key: env[key] for key in sorted(env) if key.startswith("KIC_")},
        "timeout_s": timeout,
        "input_sha256": {key: sha(path) for key, path in input_files.items()},
        "checkout_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "host": platform.platform(), "machine": platform.machine(),
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    before_load = os.getloadavg()
    with stdout.open("wb") as out_stream, stderr.open("wb") as err_stream:
        started = time.monotonic_ns()
        child = subprocess.Popen(command, cwd=REPO, env=env, stdout=out_stream, stderr=err_stream, start_new_session=True)
        timed_out = False
        rss_gate = False
        observed_group_peak_rss = 0
        while True:
            group_rss = linux_process_group_rss(child.pid)
            observed_group_peak_rss = max(observed_group_peak_rss, group_rss)
            if group_rss >= MAX_RSS:
                rss_gate = True
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _, status, usage = os.wait4(child.pid, 0)
                break
            pid, status, usage = os.wait4(child.pid, os.WNOHANG)
            if pid:
                break
            if (time.monotonic_ns() - started) / 1e9 >= timeout:
                timed_out = True
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _, status, usage = os.wait4(child.pid, 0)
                break
            time.sleep(0.02)
        child.returncode = os.waitstatus_to_exitcode(status)
        wall_ms = (time.monotonic_ns() - started) / 1e6
    peak_rss = usage.ru_maxrss if platform.system() == "Darwin" else usage.ru_maxrss * 1024
    receipt = {
        "returncode": child.returncode, "timed_out": timed_out, "rss_gate": rss_gate,
        "wall_ms": wall_ms, "user_cpu_s": usage.ru_utime, "system_cpu_s": usage.ru_stime,
        "peak_rss_bytes": peak_rss, "observed_group_peak_rss_bytes": observed_group_peak_rss,
        "load_average_before": before_load, "load_average_after": os.getloadavg(),
        "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr),
        "manifest_sha256": sha(directory / "manifest.json"),
    }
    (directory / "resource_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(basename, child.returncode, round(wall_ms, 1), peak_rss, flush=True)
    return receipt


def complete(receipt: dict) -> bool:
    return (
        receipt["returncode"] == 0 and not receipt["timed_out"] and not receipt["rss_gate"]
        and receipt["peak_rss_bytes"] < MAX_RSS and receipt["observed_group_peak_rss_bytes"] < MAX_RSS
    )


def preflight():
    assert sha(IC_SOURCE) == IC_SOURCE_HASH and sha(RHO_SOURCE) == RHO_SOURCE_HASH
    assert sha(BASE_GZ) == BASE_GZ_HASH and sha(POINTS) == POINTS_HASH
    assert sha(SHARED / "TARGET_MANIFEST.json") == MANIFEST_HASH
    parts = []
    for block, expected in enumerate(SOURCE_POINTS_HASHES):
        file = SHARED / f"points_b{block}_L128.jsonl"
        assert sha(file) == expected
        parts.append(file.read_bytes())
    assert POINTS.read_bytes() == b"".join(parts)
    points = [tuple(json.loads(line)) for line in POINTS.read_text().splitlines()]
    assert len(points) == len(set(points)) == 384
    assert subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip() == ""


def run(args):
    preflight()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    root = {
        "classification": "RUNNING",
        "protocol_sha256": sha(HERE / "PROTOCOL.md"),
        "points_sha256": sha(POINTS),
        "validator_manifest_sha256": MANIFEST_HASH,
        "base_hash": BASE_HASH,
        "base_gzip_sha256": sha(BASE_GZ),
        "checkout_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "source_sha256": {
            "ic": sha(IC_SOURCE), "rho": sha(RHO_SOURCE),
            "operational": sha(HERE / "operational.py"),
            "audit": sha(HERE / "audit.py"),
            "runner": sha(HERE / "run_panel.py"),
            "archive_sealer": sha(HERE / "archive.py"),
            "archive_verifier": sha(HERE / "verify_archive.py"),
            "training_schedule_reference": sha(ORBIT / "cold_batch_rank.py"),
            "independent_group_reference": sha(ORBIT / "independent_replay_20260924_codex/replay.py"),
            "rho_audit_math": sha(RHO_STUDY / "analyze.py"),
            "point_generation_math": sha(SHARED / "generate_targets.py"),
        },
        "binary_sha256": {"ic": sha(IC_EXE), "rho": sha(RHO_EXE)},
        "host": platform.platform(), "machine": platform.machine(),
        "steps": {},
    }
    path = args.out / "panel.json"

    def save():
        root["elapsed_wall_s"] = time.monotonic() - started
        path.write_text(json.dumps(root, indent=2, sort_keys=True) + "\n")

    def budget(stage):
        if time.monotonic() - started >= MAX_WALL:
            root["classification"] = "CENSORED_GLOBAL_WALL_CAP_BEFORE_" + stage.upper()
            save()
            return False
        return True

    def stage(key, command, env, directory, basename, timeout, files):
        if not budget(key):
            return None
        receipt = measure(command, env, directory, basename, timeout, files)
        root["steps"][key] = receipt
        save()
        return receipt

    save()
    training = args.out / "training"
    training.mkdir()
    with gzip.open(BASE_GZ, "rb") as stream:
        base = stream.read()
    assert base.count(b"\n") == 1
    header = json.loads(base)
    assert header["base_hash"] == BASE_HASH
    (training / "base_header.jsonl").write_bytes(base)
    # The training targets are intentionally known-scalar inputs. The L384
    # point validators are a separate file and are never opened here.
    import operational
    schedule = operational.training_schedule()
    (training / "target_scalars.txt").write_text("".join(f"{value}\n" for value in schedule))
    training_env = clean_env() | {
        "KIC_ALGEBRA_ENCODING": "orbit_factorized",
        "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
        "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
        "KIC_ORBIT_REP_ENCODING": "one_hot",
        "KIC_ORBIT_BATCH_ONLY": "1",
        "KIC_FACTOR_BASE_JSONL": str(training / "base_header.jsonl"),
        "KIC_ORBIT_TARGET_SCALARS": str(training / "target_scalars.txt"),
        "KIC_TASK_ID": "TASK-IC-N53-COMBINED-L384-20260925",
    }
    ic_command = [str(IC_EXE), "53", "0", "1", "10", "natural", "1", "2000", "1", "internal"]
    training_child = stage(
        "training_producer", ic_command, training_env, training, "producer", 360,
        {"ic_source": IC_SOURCE, "ic_binary": IC_EXE, "base_header": training / "base_header.jsonl",
         "training_scalars": training / "target_scalars.txt"},
    )
    if training_child is None or not complete(training_child):
        root["classification"] = "CENSORED_TRAINING_PRODUCER"
        save()
        return
    rank_receipt = stage(
        "operational_rank",
        [sys.executable, str(HERE / "operational.py"), "--stage", "training",
         "--training", str(training), "--out", str(training / "operational_solution.json")],
        clean_env(), training / "operational_rank", "rank", 180,
        {"operational_source": HERE / "operational.py", "producer_raw": training / "producer.stdout.jsonl",
         "base_header": training / "base_header.jsonl", "schedule": training / "target_scalars.txt"},
    )
    if rank_receipt is None or not complete(rank_receipt):
        root["classification"] = "INVALID_OR_CENSORED_OPERATIONAL_RANK"
        save()
        return
    rank_report = json.loads((training / "operational_solution.json").read_text())
    if rank_report["rank"] != 220 or rank_report["relations_used"] != 512:
        root["classification"] = "INVALID_TRAINING_RANK_OR_RELATIONS"
        save()
        return
    rho_env = clean_env() | {
        "KIC_RHO_TARGET_POINTS_JSONL": str(POINTS),
        "KIC_RHO_BATCH_CORPUS": "n53-combined-L384-20260925-v1",
        "KIC_RHO_DP_BITS": "4", "KIC_RHO_PRECOMPUTE_WALKS": "0",
    }
    rho = stage(
        "rho", [str(RHO_EXE), "53", "0", "signed_frobenius", "384", "531384"],
        rho_env, args.out / "rho", "rho", 600,
        {"rho_source": RHO_SOURCE, "rho_binary": RHO_EXE, "points": POINTS},
    )
    ic_env = clean_env() | {
        "KIC_ALGEBRA_ENCODING": "orbit_factorized",
        "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
        "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
        "KIC_ORBIT_REP_ENCODING": "one_hot",
        "KIC_ORBIT_BATCH_ONLY": "1",
        "KIC_FACTOR_BASE_JSONL": str(training / "base_header.jsonl"),
        "KIC_ORBIT_TARGET_POINTS_JSONL": str(POINTS),
        "KIC_TASK_ID": "TASK-IC-N53-COMBINED-L384-20260925",
    }
    ic = stage(
        "ic", ic_command, ic_env, args.out / "ic", "ic", 300,
        {"ic_source": IC_SOURCE, "ic_binary": IC_EXE, "base_header": training / "base_header.jsonl",
         "points": POINTS},
    )
    if ic is not None and complete(ic):
        recovery = stage(
            "operational_recovery",
            [sys.executable, str(HERE / "operational.py"), "--stage", "point",
             "--training", str(training), "--raw", str(args.out / "ic"),
             "--points", str(POINTS), "--out", str(args.out / "ic/operational_recovery.json")],
            clean_env(), args.out / "ic/operational_recovery", "recovery", 180,
            {"operational_source": HERE / "operational.py",
             "ic_raw": args.out / "ic/ic.stdout.jsonl", "points": POINTS,
             "training_solution": training / "operational_solution.json"},
        )
    else:
        recovery = None
    if rho is None or ic is None or recovery is None or not all(map(complete, (rho, ic, recovery))):
        root["classification"] = "CENSORED_CHILD_OR_OPERATIONAL_RECOVERY"
        save()
        return
    audit = stage(
        "independent_audit",
        [sys.executable, str(HERE / "audit.py"), "--panel", str(args.out),
         "--out", str(args.out / "audit.json")],
        clean_env(), args.out / "audit_driver", "audit", 600,
        {"audit_source": HERE / "audit.py", "points": POINTS,
         "training_raw": training / "producer.stdout.jsonl",
         "training_solution": training / "operational_solution.json",
         "ic_raw": args.out / "ic/ic.stdout.jsonl",
         "ic_solution": args.out / "ic/operational_recovery.json",
         "rho_raw": args.out / "rho/rho.stdout.jsonl"},
    )
    if audit is None or not complete(audit):
        root["classification"] = "INVALID_INDEPENDENT_AUDIT"
        save()
        return
    replay = json.loads((args.out / "audit.json").read_text())
    assert replay["classification"] == "PASS_ALL_384_SAME_Q_LOGS_AND_512_TRAINING_RELATIONS"
    operational_ms = sum(root["steps"][name]["wall_ms"] for name in (
        "training_producer", "operational_rank", "ic", "operational_recovery"))
    lower_ms = training_child["wall_ms"] + ic["wall_ms"]
    rho_ms = rho["wall_ms"]
    audit_inclusive_ms = operational_ms + audit["wall_ms"]
    root["costs"] = {
        "ic_operational_wall_ms": operational_ms,
        "ic_two_child_lower_wall_ms": lower_ms,
        "ic_audit_inclusive_wall_ms": audit_inclusive_ms,
        "rho_operational_wall_ms": rho_ms,
        "operational_ratio_to_rho": operational_ms / rho_ms,
        "two_child_lower_ratio_to_rho": lower_ms / rho_ms,
        "audit_inclusive_ratio_to_rho": audit_inclusive_ms / rho_ms,
        "ic_operational_cpu_s": sum(
            root["steps"][name]["user_cpu_s"] + root["steps"][name]["system_cpu_s"]
            for name in ("training_producer", "operational_rank", "ic", "operational_recovery")
        ),
        "rho_cpu_s": rho["user_cpu_s"] + rho["system_cpu_s"],
        "ic_operational_peak_rss_bytes": max(
            max(root["steps"][name]["peak_rss_bytes"], root["steps"][name]["observed_group_peak_rss_bytes"])
            for name in ("training_producer", "operational_rank", "ic", "operational_recovery")
        ),
        "rho_peak_rss_bytes": max(rho["peak_rss_bytes"], rho["observed_group_peak_rss_bytes"]),
    }
    root["native_operations"] = {
        "training": replay["training"], "ic": replay["ic"], "rho": replay["rho"],
    }
    root["classification"] = (
        "COMPLETE_FIXED_STREAM_OPERATIONAL_WALL_WIN" if operational_ms < rho_ms
        else "COMPLETE_FIXED_STREAM_NO_CROSSOVER" if lower_ms > rho_ms
        else "COMPLETE_FIXED_STREAM_INCONCLUSIVE"
    )
    save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args)
    classification = json.loads((args.out / "panel.json").read_text())["classification"]
    if classification.startswith("INVALID"):
        raise SystemExit(classification)
