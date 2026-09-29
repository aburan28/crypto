#!/usr/bin/env python3
"""Two-arm, cold native-cyclic n53 L384 IC versus one same-Q batched rho."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
IC_SOURCE = REPO / "examples/koblitz_s5_sat_instance.rs"
RHO_SOURCE = REPO / "examples/koblitz_rho_batch_ks.rs"
IC_EXE = REPO / "target/release/examples/koblitz_s5_sat_instance"
RHO_EXE = REPO / "target/release/examples/koblitz_rho_batch_ks"
POINTS = HERE / "points_L384.jsonl"
MAX_RSS = 2 * 1024**3
MAX_WALL = 4200
MAX_AUDIT_WALL = 1200
MAX_OPERATIONAL_WALL = MAX_WALL - MAX_AUDIT_WALL

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
        "environment": {key: env[key] for key in sorted(env) if key.startswith("KIC_") or key == "RAYON_NUM_THREADS"},
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


def stage_timeout(key: str, requested: float, steps: dict) -> float:
    """Reserve the full audit child slot even after operational cap exhaustion."""
    if key == "independent_audit":
        return min(requested, MAX_AUDIT_WALL)
    used = sum(receipt.get("wall_ms", 0) for name, receipt in steps.items()
               if name != "independent_audit") / 1000
    return max(0.0, min(requested, MAX_OPERATIONAL_WALL - used))


def audit_input_files(panel: Path) -> dict[str, Path]:
    files = {"audit_source": HERE / "audit.py", "points": POINTS,
             "validator_scalars": HERE / "validator_scalars_L384.txt"}
    for arm in ("A", "B"):
        for label, path in {
            f"{arm}_training_raw": panel / arm / "training/producer.stdout.jsonl",
            f"{arm}_rank": panel / arm / "training/operational_solution.json",
            f"{arm}_ic_raw": panel / arm / "ic/ic.stdout.jsonl",
            f"{arm}_recovery": panel / arm / "ic/operational_recovery.json",
        }.items():
            if path.is_file():
                files[label] = path
    rho_raw = panel / "rho/rho.stdout.jsonl"
    if rho_raw.is_file():
        files["rho_raw"] = rho_raw
    return files


def audit_command(panel: Path) -> list[str]:
    return [sys.executable, str(HERE / "audit.py"), "--panel", str(panel),
            "--out", str(panel / "audit.json")]


def has_complete_raw_prefix(panel: Path) -> bool:
    children = [(panel / "rho", "rho")]
    children.extend((panel / arm / stage, name)
                    for arm in ("A", "B")
                    for stage, name in (("training", "producer"), ("ic", "ic")))
    for directory, basename in children:
        receipt_path = directory / "resource_receipt.json"
        raw_path = directory / f"{basename}.stdout.jsonl"
        if receipt_path.is_file() and raw_path.is_file():
            try:
                if complete(json.loads(receipt_path.read_text())):
                    return True
            except (KeyError, ValueError):
                continue
    return False


def audit_after_exception(panel: Path, summary: dict) -> dict:
    """Attempt one independent replay of successful raw prefixes after a runner error."""
    if not has_complete_raw_prefix(panel):
        return summary
    sequence = summary.setdefault("stage_sequence", [])
    if "independent_audit" in sequence:
        receipt_path = panel / "audit_driver/resource_receipt.json"
        if receipt_path.is_file() and "independent_audit" not in summary.get("steps", {}):
            receipt = json.loads(receipt_path.read_text())
            summary.setdefault("steps", {})["independent_audit"] = receipt
            if complete(receipt) and (panel / "audit.json").is_file():
                summary["replay"] = json.loads((panel / "audit.json").read_text())
        return summary  # An audit was already attempted; never retry it.
    sequence.append("independent_audit")
    try:
        receipt = measure(audit_command(panel), clean_env(), panel / "audit_driver",
                          "audit", MAX_AUDIT_WALL, audit_input_files(panel))
        summary.setdefault("steps", {})["independent_audit"] = receipt
        if complete(receipt):
            summary["replay"] = json.loads((panel / "audit.json").read_text())
    except Exception as exc:
        failure = {"error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()}
        failure_path = panel / "audit_failure.json"
        failure_path.write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n")
        summary["audit_failure_sha256"] = sha(failure_path)
    return summary


def preflight(panel: Path):
    import check_protocol
    frozen = check_protocol.preflight(require_release=True)
    assert sha(IC_SOURCE) == frozen["source_sha256"]["ic"]
    assert sha(RHO_SOURCE) == frozen["source_sha256"]["rho"]
    assert sha(POINTS) == frozen["input_sha256"]["public_points"]
    assert not (HERE / "evidence/archive_manifest.json").exists()
    build = json.loads((panel / "build_receipt.json").read_text())
    assert build["status"] == "SUCCESS"
    assert build["binary_sha256"] == {"ic": sha(IC_EXE), "rho": sha(RHO_EXE)}
    return frozen


def run(args):
    frozen = preflight(args.out)
    assert sha(IC_EXE) and sha(RHO_EXE)
    dispatch_main_head = subprocess.check_output(
        ["git", "rev-parse", "origin/main"], cwd=REPO, text=True).strip()
    subprocess.run(["git", "merge-base", "--is-ancestor",
                    frozen["release_main_head"], dispatch_main_head], cwd=REPO, check=True)
    assert args.out.is_dir() and {item.name for item in args.out.iterdir()} == {
        "predispatch.json", "build_receipt.json", "build.stdout.txt", "build.stderr.txt"}
    gate = json.loads((args.out / "predispatch.json").read_text())
    assert gate["status"] == "ADMITTED"
    assert gate["reviewed_head"] == os.environ["KIC_NATIVE_L384_EXPECTED_HEAD"]
    assert gate["run_id"] == os.environ["GITHUB_RUN_ID"]
    started = time.monotonic()
    root = {
        "schema": "n53_native_cyclic_l384_panel_v1",
        "classification": "RUNNING",
        "protocol_sha256": sha(HERE / "PROTOCOL.md"),
        "frozen_sha256": sha(HERE / "FROZEN.json"),
        "points_sha256": sha(POINTS),
        "checkout_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "release_main_head": frozen["release_main_head"],
        "dispatch_main_head": dispatch_main_head,
        "source_sha256": frozen["source_sha256"],
        "input_sha256": frozen["input_sha256"],
        "binary_sha256": {"ic": sha(IC_EXE), "rho": sha(RHO_EXE)},
        "build_receipt_sha256": sha(args.out / "build_receipt.json"),
        "host": platform.platform(), "machine": platform.machine(),
        "github_run_id": os.environ.get("GITHUB_RUN_ID"),
        "github_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "steps": {}, "stage_sequence": [], "arms": {},
        "stage_budget_seconds": {"operational": MAX_OPERATIONAL_WALL,
                                 "independent_audit": MAX_AUDIT_WALL},
    }
    path = args.out / "panel.json"

    def save():
        root["elapsed_wall_s"] = time.monotonic() - started
        path.write_text(json.dumps(root, indent=2, sort_keys=True) + "\n")

    def stage(key, command, env, directory, basename, timeout, files):
        root["stage_sequence"].append(key)
        save()  # Preserve the attempted stage even if the child cannot start.
        allowed = stage_timeout(key, timeout, root["steps"])
        if allowed <= 0:
            assert key != "independent_audit", "independent audit has a reserved slot"
            root["steps"][key] = {"not_run": "GLOBAL_WALL_CAP"}
            save()
            return None
        receipt = measure(command, env, directory, basename, allowed, files)
        root["steps"][key] = receipt
        save()
        return receipt

    ic_command = [str(IC_EXE), "53", "0", "1", "10", "natural", "1", "2000", "1", "internal"]
    rho_env = clean_env() | {
        "KIC_RHO_TARGET_POINTS_JSONL": str(POINTS),
        "KIC_RHO_BATCH_CORPUS": "n53-native-cyclic-L384-20260929-v1",
        "KIC_RHO_DP_BITS": "4", "KIC_RHO_PRECOMPUTE_WALKS": "0",
        "RAYON_NUM_THREADS": "1",
    }
    save()
    rho = stage("rho", [str(RHO_EXE), "53", "0", "signed_frobenius", "384", "531929"],
                rho_env, args.out / "rho", "rho", 600,
                {"rho_source": RHO_SOURCE, "rho_binary": RHO_EXE, "points": POINTS})
    # B is the fresh primary transfer arm, run first. A is the registered
    # #823 stream and remains a positive control even if B fails.
    for arm in ("B", "A"):
        arm_dir = args.out / arm
        training = arm_dir / "training"
        training.mkdir(parents=True)
        schedule = HERE / f"training_{arm}_scalars.txt"
        (training / "target_scalars.txt").write_bytes(schedule.read_bytes())
        state = {"training_source_sha256": sha(schedule), "status": "RUNNING"}
        root["arms"][arm] = state
        save()
        training_env = clean_env() | {
            "KIC_ALGEBRA_ENCODING": "orbit_factorized",
            "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
            "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
            "KIC_ORBIT_REP_ENCODING": "one_hot",
            "KIC_ORBIT_BATCH_ONLY": "1",
            "KIC_ORBIT_INCLUDE_BASE_HEADER": "1",
            "KIC_FACTOR_BASE_SELECTION": "ascending_x_v1",
            "KIC_ORBIT_REGULAR_SCAN_POLICY": "target_cyclic_v1",
            "KIC_ORBIT_TARGET_SCALARS": str(training / "target_scalars.txt"),
            "KIC_TASK_ID": f"TASK-IC-N53-NATIVE-CYCLIC-L384-20260929-{arm}",
            "RAYON_NUM_THREADS": "1",
        }
        training_child = stage(
            f"{arm}_training", ic_command, training_env, training, "producer", 360,
            {"ic_source": IC_SOURCE, "ic_binary": IC_EXE,
             "training_scalars": training / "target_scalars.txt"},
        )
        if training_child is None or not complete(training_child):
            state["status"] = "CENSORED_TRAINING_CHILD"
            save()
            continue
        base_receipt = stage(
            f"{arm}_base", [sys.executable, str(HERE / "cold_base.py"),
                "--training-raw", str(training / "producer.stdout.jsonl"),
                "--header", str(training / "base_header.jsonl"),
                "--receipt", str(training / "base_materialization.json")],
            clean_env(), training / "base_materialization", "base", 90,
            {"materializer_source": HERE / "cold_base.py",
             "training_raw": training / "producer.stdout.jsonl"},
        )
        if base_receipt is None or not complete(base_receipt):
            state["status"] = "INVALID_OR_CENSORED_NATIVE_BASE"
            save()
            continue
        base_report = json.loads((training / "base_materialization.json").read_text())
        assert base_report["classification"] == "FRESH_NATIVE_BASE_MATCHES_PINNED_ORDER"
        state["base_header_sha256"] = sha(training / "base_header.jsonl")
        rank_receipt = stage(
            f"{arm}_rank", [sys.executable, str(HERE / "operational.py"), "--stage", "training",
                "--training", str(training), "--out", str(training / "operational_solution.json")],
            clean_env(), training / "operational_rank", "rank", 180,
            {"operational_source": HERE / "operational.py",
             "producer_raw": training / "producer.stdout.jsonl",
             "base_header": training / "base_header.jsonl",
             "schedule": training / "target_scalars.txt"},
        )
        if rank_receipt is None or not complete(rank_receipt):
            state["status"] = "INVALID_OR_CENSORED_RANK"
            save()
            continue
        rank_report = json.loads((training / "operational_solution.json").read_text())
        state["rank"] = rank_report["rank"]
        state["training_failed_targets"] = len(rank_report["failed_target_scalars"])
        if rank_report["classification"] != "OPERATIONAL_TRAINING_LOGS_RECOVERED":
            state["status"] = "RANK_DEFICIENT"
            save()
            continue
        ic_env = clean_env() | {
            "KIC_ALGEBRA_ENCODING": "orbit_factorized",
            "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
            "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
            "KIC_ORBIT_REP_ENCODING": "one_hot",
            "KIC_ORBIT_BATCH_ONLY": "1",
            "KIC_FACTOR_BASE_JSONL": str(training / "base_header.jsonl"),
            "KIC_ORBIT_REGULAR_SCAN_POLICY": "target_cyclic_v1",
            "KIC_ORBIT_TARGET_POINTS_JSONL": str(POINTS),
            "KIC_TASK_ID": f"TASK-IC-N53-NATIVE-CYCLIC-L384-20260929-{arm}",
            "RAYON_NUM_THREADS": "1",
        }
        ic = stage(
            f"{arm}_ic", ic_command, ic_env, arm_dir / "ic", "ic", 300,
            {"ic_source": IC_SOURCE, "ic_binary": IC_EXE,
             "base_header": training / "base_header.jsonl", "points": POINTS},
        )
        if ic is None or not complete(ic):
            state["status"] = "CENSORED_POINT_CHILD"
            save()
            continue
        recovery = stage(
            f"{arm}_recovery", [sys.executable, str(HERE / "operational.py"), "--stage", "point",
                "--training", str(training), "--raw", str(arm_dir / "ic"),
                "--points", str(POINTS), "--out", str(arm_dir / "ic/operational_recovery.json")],
            clean_env(), arm_dir / "ic/operational_recovery", "recovery", 180,
            {"operational_source": HERE / "operational.py",
             "ic_raw": arm_dir / "ic/ic.stdout.jsonl", "points": POINTS,
             "training_solution": training / "operational_solution.json"},
        )
        if recovery is None or not complete(recovery):
            state["status"] = "INVALID_OR_CENSORED_RECOVERY"
            save()
            continue
        recovery_report = json.loads((arm_dir / "ic/operational_recovery.json").read_text())
        state["point_logs_recovered"] = recovery_report["count"]
        state["status"] = ("OPERATIONAL_384_LOGS" if recovery_report["count"] == 384
                           else "INCOMPLETE_POINT_LOGS")
        save()
    audit = stage("independent_audit", audit_command(args.out), clean_env(),
                  args.out / "audit_driver", "audit", MAX_AUDIT_WALL,
                  audit_input_files(args.out))
    if audit is None or not complete(audit):
        root["classification"] = "INVALID_OR_CENSORED_INDEPENDENT_AUDIT"
        save()
        return
    replay = json.loads((args.out / "audit.json").read_text())
    if replay["classification"] != "INDEPENDENT_FULL_REPLAY":
        root["classification"] = "PARTIAL_REPLAY_NO_END_TO_END_VERDICT"
        root["replay"] = replay
        save()
        return
    assert rho is not None and complete(rho)
    assert all(root["arms"][arm]["status"] == "OPERATIONAL_384_LOGS" for arm in ("A", "B"))
    rho_ms = rho["wall_ms"]
    costs = {"rho_operational_wall_ms": rho_ms,
             "rho_cpu_s": rho["user_cpu_s"] + rho["system_cpu_s"]}
    for arm in ("A", "B"):
        keys = [f"{arm}_{stage}" for stage in ("training", "base", "rank", "ic", "recovery")]
        assert all(complete(root["steps"][key]) for key in keys)
        operational = sum(root["steps"][key]["wall_ms"] for key in keys)
        lower = sum(root["steps"][f"{arm}_{stage}"]["wall_ms"]
                    for stage in ("training", "base", "ic"))
        costs[arm] = {
            "ic_operational_wall_ms": operational,
            "ic_two_child_lower_wall_ms": lower,
            "ic_to_same_rho_ratio": operational / rho_ms,
            "lower_to_same_rho_ratio": lower / rho_ms,
            "ic_operational_cpu_s": sum(root["steps"][key]["user_cpu_s"] +
                                        root["steps"][key]["system_cpu_s"] for key in keys),
            "classification": ("FIXED_STREAM_WALL_WIN" if operational < rho_ms
                               else "FIXED_STREAM_NO_CROSSOVER"),
        }
    root["costs"] = costs
    root["replay"] = replay
    root["classification"] = ("COMPLETE_PRIMARY_B_WALL_WIN" if costs["B"]["ic_operational_wall_ms"] < rho_ms
                              else "COMPLETE_PRIMARY_B_NO_CROSSOVER")
    save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--gate-only", action="store_true")
    args = parser.parse_args()
    if args.gate_only:
        import check_protocol
        try:
            check_protocol.preflight(require_release=True)
        except AssertionError as exc:
            if "HELD" not in str(exc):
                raise
            print(json.dumps({"gate": "HELD_NO_CHILDREN", "reason": str(exc)}, sort_keys=True))
            return
        raise AssertionError("release gate unexpectedly open; gate-only check must not run a child")
    assert args.out is not None
    try:
        run(args)
    except Exception as exc:
        if args.out.is_dir():
            failure = {"classification": "INVALID_RUNNER_EXCEPTION",
                       "error": f"{type(exc).__name__}: {exc}",
                       "traceback": traceback.format_exc()}
            (args.out / "failure.json").write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n")
            panel_path = args.out / "panel.json"
            panel = json.loads(panel_path.read_text()) if panel_path.is_file() else {}
            panel["classification"] = "INVALID_RUNNER_EXCEPTION"
            panel["failure_sha256"] = sha(args.out / "failure.json")
            panel = audit_after_exception(args.out, panel)
            panel_path.write_text(json.dumps(panel, indent=2, sort_keys=True) + "\n")
        raise
    classification = json.loads((args.out / "panel.json").read_text())["classification"]
    if classification.startswith("INVALID"):
        raise SystemExit(classification)


if __name__ == "__main__":
    main()
