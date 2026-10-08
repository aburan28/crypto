#!/usr/bin/env python3
"""One-host, one-training, paired n53 compact IC / batched-rho panel."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORBIT = REPO / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"
BASE_GZ = ORBIT / "independent_replay_20260924_codex/base_header.jsonl.gz"
BASE_HASH = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"
IC_SOURCE = REPO / "examples/koblitz_s5_sat_instance.rs"
RHO_SOURCE = REPO / "examples/koblitz_rho_batch_ks.rs"
IC_EXE = REPO / "target/release/examples/koblitz_s5_sat_instance"
RHO_EXE = REPO / "target/release/examples/koblitz_rho_batch_ks"
MAX_RSS = 2 * 1024**3
MAX_WALL = 3000


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean_env() -> dict[str, str]:
    return {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}


def linux_process_group_rss(pgid: int) -> int:
    """Current summed RSS for the fresh measured process group on Linux."""
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
    directory.mkdir(parents=True, exist_ok=False)
    stdout = directory / f"{basename}.stdout.jsonl"
    stderr = directory / f"{basename}.stderr.txt"
    source_hashes = {key: sha(path) for key, path in input_files.items()}
    manifest = {
        "command": command, "environment": {key: env[key] for key in sorted(env) if key.startswith("KIC_")},
        "timeout_s": timeout, "input_sha256": source_hashes,
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
            elapsed = (time.monotonic_ns() - started) / 1e9
            if elapsed >= timeout:
                timed_out = True
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _, status, usage = os.wait4(child.pid, 0)
                break
            time.sleep(0.02)
        wall_ms = (time.monotonic_ns() - started) / 1e6
        child.returncode = os.waitstatus_to_exitcode(status)
    peak_rss = usage.ru_maxrss if platform.system() == "Darwin" else usage.ru_maxrss * 1024
    receipt = {
        "returncode": child.returncode, "timed_out": timed_out, "rss_gate": rss_gate,
        "wall_ms": wall_ms, "user_cpu_s": usage.ru_utime,
        "system_cpu_s": usage.ru_stime, "peak_rss_bytes": peak_rss,
        "observed_group_peak_rss_bytes": observed_group_peak_rss,
        "load_average_before": before_load, "load_average_after": os.getloadavg(),
        "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr),
        "manifest_sha256": sha(directory / "manifest.json"),
    }
    (directory / "resource_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(basename, directory.name, child.returncode, round(wall_ms, 1), peak_rss, flush=True)
    return receipt


def complete(receipt: dict) -> bool:
    return (receipt["returncode"] == 0 and not receipt["timed_out"]
            and not receipt["rss_gate"] and receipt["peak_rss_bytes"] < MAX_RSS
            and receipt["observed_group_peak_rss_bytes"] < MAX_RSS)


def checked_targets():
    frozen = json.loads((HERE / "TARGET_MANIFEST.json").read_text())
    with tempfile.TemporaryDirectory() as temp:
        generated = Path(temp) / "generated"
        generated_process = subprocess.run(
            [sys.executable, str(HERE / "generate_targets.py"), "--out", str(generated)],
            cwd=REPO, capture_output=True, text=True,
        )
        assert generated_process.returncode == 0, generated_process.stderr
        assert (generated / "TARGET_MANIFEST.json").read_bytes() == (HERE / "TARGET_MANIFEST.json").read_bytes()
        for row in frozen["blocks"]:
            for count in (32, 128):
                name = row[f"L{count}_file"]
                assert (generated / name).read_bytes() == (HERE / name).read_bytes()
                assert sha(HERE / name) == row[f"L{count}_sha256"]
                if count == 32:
                    assert (HERE / name).read_bytes() == b"".join(
                        (HERE / row["L128_file"]).read_bytes().splitlines(keepends=True)[:32]
                    )
    return frozen


def run(args):
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    frozen = checked_targets()
    git_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    git_status = subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip()
    assert not git_status, git_status
    root = {
        "protocol_sha256": sha(HERE / "PROTOCOL.md"),
        "target_manifest_sha256": sha(HERE / "TARGET_MANIFEST.json"),
        "generator_sha256": sha(HERE / "generate_targets.py"),
        "checkout_head": git_head, "source_sha256": {
            "ic": sha(IC_SOURCE), "rho": sha(RHO_SOURCE),
            "training_driver": sha(ORBIT / "cold_batch_rank.py"),
            "independent_replay": sha(ORBIT / "independent_replay_20260924_codex/replay.py"),
            "pair_verifier": sha(HERE / "verify_pair.py"),
            "panel_runner": sha(HERE / "run_panel.py"),
            "archive_sealer": sha(HERE / "archive.py"),
        },
        "binary_sha256": {"ic": sha(IC_EXE), "rho": sha(RHO_EXE)},
        "base_gzip_sha256": sha(BASE_GZ), "base_hash": BASE_HASH,
        "host": platform.platform(), "machine": platform.machine(),
        "steps": [], "classification": "RUNNING",
    }
    root_path = args.out / "panel.json"

    def save():
        root["elapsed_wall_s"] = time.monotonic() - started
        root_path.write_text(json.dumps(root, indent=2, sort_keys=True) + "\n")

    save()
    training = args.out / "training"
    rank_script = ORBIT / "cold_batch_rank.py"
    driver = measure(
        [sys.executable, str(rank_script), "--out", str(training), "--targets", "512", "--timeout", "330"],
        clean_env(), args.out / "training_driver", "training_driver", 360,
        {"training_driver": rank_script, "ic_binary": IC_EXE, "base_gzip": BASE_GZ},
    )
    root["training_driver"] = driver
    if not complete(driver) or not (training / "validation.json").exists():
        root["classification"] = "CENSORED_TRAINING_PRODUCER_OR_REPLAY"
        save()
        return
    validation = json.loads((training / "validation.json").read_text())
    if not (validation["rank"] == validation["columns"] == 220
            and validation["factor_base_log_solution_verified"]
            and validation["targets_extracted"] == 512):
        root["classification"] = "INVALID_TRAINING_RANK_OR_RELATIONS"
        save()
        return
    root["training_child_wall_ms"] = validation["process_wall_ms"]
    training_observation = json.loads((training / "producer.stdout.jsonl").read_text())
    training_batch = training_observation["compact_orbit_batch"]
    training_queries = training_batch["query_observations"]
    training_scalars = [int(line) for line in (training / "target_scalars.txt").read_text().splitlines()]
    assert len(training_queries) == len(training_scalars) == 512
    assert all(row["scalar"] == scalar and row["hit"]
               for row, scalar in zip(training_queries, training_scalars))
    root["training_operations"] = {
        "regular_states": training_batch["regular_states"],
        "index_entries": training_batch["index_entries"],
        "query_observations": training_queries,
        "s3_calls_sum": sum(row["s3_calls"] for row in training_queries),
        "partner_roots_sum": sum(row["partner_roots"] for row in training_queries),
        "indexed_partner_hits_sum": sum(row["indexed_partner_hits"] for row in training_queries),
        "group_lift_attempts_sum": sum(row["group_lift_attempts"] for row in training_queries),
        "query_ms_sum": training_batch["query_ms_sum"],
        "rank": validation["rank"],
        "first_full_rank_at_extracted": validation["full_rank_at_extracted"],
    }
    root["training_child_peak_rss_bytes"] = (
        validation["child_peak_rss_raw"] if platform.system() == "Darwin"
        else validation["child_peak_rss_raw"] * 1024
    )
    if root["training_child_peak_rss_bytes"] >= MAX_RSS:
        root["classification"] = "CENSORED_TRAINING_RSS"
        save()
        return
    save()

    def pair(block: int, count: int) -> bool:
        row = frozen["blocks"][block]
        point_file = HERE / row[f"L{count}_file"]
        assert sha(point_file) == row[f"L{count}_sha256"]
        pair_dir = args.out / f"L{count}_b{block}"
        pair_dir.mkdir()
        common = clean_env()
        ic_env = common | {
            "KIC_ALGEBRA_ENCODING": "orbit_factorized",
            "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
            "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
            "KIC_ORBIT_REP_ENCODING": "one_hot",
            "KIC_ORBIT_BATCH_ONLY": "1",
            "KIC_FACTOR_BASE_JSONL": str(training / "base_header.jsonl"),
            "KIC_ORBIT_TARGET_POINTS_JSONL": str(point_file),
            "KIC_TASK_ID": "TASK-IC-N53-TRUE-SHARED-LOG-20260925",
        }
        rho_env = common | {
            "KIC_RHO_TARGET_POINTS_JSONL": str(point_file),
            "KIC_RHO_BATCH_CORPUS": f"n53-shared-log-20260925-b{block}-L{count}",
            "KIC_RHO_DP_BITS": "4", "KIC_RHO_PRECOMPUTE_WALKS": "0",
        }
        commands = {
            "ic": ([str(IC_EXE), "53", "0", "1", "10", "natural", "1", "2000", "1", "internal"], ic_env, IC_SOURCE),
            "rho": ([str(RHO_EXE), "53", "0", "signed_frobenius", str(count), str(row["rho_seed"])], rho_env, RHO_SOURCE),
        }
        order = ("ic", "rho") if block != 1 else ("rho", "ic")
        record = {"block": block, "count": count, "order": order, "points_sha256": sha(point_file), "arms": {}}
        root["steps"].append(record)
        save()
        for arm in order:
            if time.monotonic() - started >= MAX_WALL:
                record["classification"] = "CENSORED_GLOBAL_WALL_CAP"
                save()
                return False
            command, env, source = commands[arm]
            receipt = measure(command, env, pair_dir / arm, arm,
                              180 if count == 32 else 300,
                              {"points": point_file, "source": source,
                               "binary": IC_EXE if arm == "ic" else RHO_EXE})
            record["arms"][arm] = receipt
            save()
        if not all(complete(record["arms"][arm]) for arm in ("ic", "rho")):
            record["classification"] = "CENSORED_CHILD_FAILURE_TIMEOUT_OR_RSS"
            save()
            return False
        for arm in ("ic", "rho"):
            if time.monotonic() - started >= MAX_WALL:
                record["classification"] = "CENSORED_GLOBAL_WALL_CAP"
                save()
                return False
            receipt = measure(
                [sys.executable, str(HERE / "verify_pair.py"), "--arm", arm,
                 "--training", str(training), "--block", str(block), "--count", str(count),
                 "--raw", str(pair_dir / arm), "--out", str(pair_dir / arm / "validation.json")],
                clean_env(), pair_dir / f"verify_{arm}", f"verify_{arm}", 300,
                {"verifier": HERE / "verify_pair.py", "points": point_file,
                 "raw": pair_dir / arm / f"{arm}.stdout.jsonl"},
            )
            record[f"verify_{arm}"] = receipt
            save()
            if not complete(receipt):
                record["classification"] = "INVALID_INDEPENDENT_REPLAY"
                save()
                return False
        ic = json.loads((pair_dir / "ic/validation.json").read_text())
        rho = json.loads((pair_dir / "rho/validation.json").read_text())
        assert ic["recovered_scalars"] == rho["recovered_scalars"]
        record["classification"] = "PASS"
        record["ic_operations"] = {
            "regular_states": ic["regular_states"],
            "index_entries": ic["index_entries"],
            "partner_trials_sum": ic["partner_trials_sum"],
            "query_ms_sum": ic["query_ms_sum"],
            "batch_loop_wall_ms": ic["batch_loop_wall_ms"],
            "s3_calls_sum": sum(row["s3_calls"] for row in ic["query_observations"]),
            "partner_roots_sum": sum(row["partner_roots"] for row in ic["query_observations"]),
            "indexed_partner_hits_sum": sum(row["indexed_partner_hits"] for row in ic["query_observations"]),
            "group_lift_attempts_sum": sum(row["group_lift_attempts"] for row in ic["query_observations"]),
        }
        record["rho_operations"] = {
            "walk_steps": rho["walk_steps"], "table_entries": rho["table_entries"],
            "cross_target_solves": rho["cross_target_solves"], "charges": rho["charges"],
        }
        record["ic_lower_wall_ms"] = root["training_child_wall_ms"] + record["arms"]["ic"]["wall_ms"]
        record["ic_conservative_upper_wall_ms"] = (driver["wall_ms"] + record["arms"]["ic"]["wall_ms"]
                                                    + record["verify_ic"]["wall_ms"])
        record["rho_wall_ms"] = record["arms"]["rho"]["wall_ms"]
        record["ic_over_rho_lower"] = record["ic_lower_wall_ms"] / record["rho_wall_ms"]
        record["ic_over_rho_upper"] = record["ic_conservative_upper_wall_ms"] / record["rho_wall_ms"]
        record["ic_conservative_upper_cpu_s"] = sum(
            item["user_cpu_s"] + item["system_cpu_s"] for item in
            (driver, record["arms"]["ic"], record["verify_ic"])
        )
        record["rho_cpu_s"] = (record["arms"]["rho"]["user_cpu_s"]
                               + record["arms"]["rho"]["system_cpu_s"])
        record["ic_peak_rss_bytes"] = max(
            root["training_child_peak_rss_bytes"], driver["peak_rss_bytes"],
            driver["observed_group_peak_rss_bytes"],
            record["arms"]["ic"]["peak_rss_bytes"],
            record["arms"]["ic"]["observed_group_peak_rss_bytes"],
            record["verify_ic"]["peak_rss_bytes"],
        )
        record["rho_peak_rss_bytes"] = max(
            record["arms"]["rho"]["peak_rss_bytes"],
            record["arms"]["rho"]["observed_group_peak_rss_bytes"],
        )
        save()
        return True

    for block in range(3):
        if not pair(block, 32):
            root["classification"] = (
                "INVALID_L32_REPLAY" if root["steps"][-1]["classification"].startswith("INVALID")
                else "CENSORED_L32_CHILD_OR_CAP"
            )
            save()
            return
    l32 = [x for x in root["steps"] if x["count"] == 32]
    max_query_ms = max(x["arms"]["ic"]["wall_ms"] for x in l32)
    max_rss = max(x["arms"]["ic"]["peak_rss_bytes"] for x in l32)
    root["L128_gate"] = {
        "four_times_max_L32_ic_wall_ms": 4 * max_query_ms,
        "max_L32_ic_rss_bytes": max_rss,
        "pass": 4 * max_query_ms < 300000 and max_rss < int(1.75 * 1024**3),
    }
    if not root["L128_gate"]["pass"]:
        root["classification"] = "CENSORED_L128_PREDICTED_COST_OR_RSS"
        save()
        return
    for block in range(3):
        if not pair(block, 128):
            root["classification"] = (
                "INVALID_L128_REPLAY" if root["steps"][-1]["classification"].startswith("INVALID")
                else "CENSORED_L128_CHILD_OR_CAP"
            )
            save()
            return
    root["classification"] = "COMPLETE_ALL_SIX_PAIRS"
    for count in (32, 128):
        rows = [step for step in root["steps"] if step["count"] == count]
        root[f"L{count}_portfolio"] = {
            "ic_upper_wall_ms_training_once": driver["wall_ms"] + sum(
                step["arms"]["ic"]["wall_ms"] + step["verify_ic"]["wall_ms"] for step in rows
            ),
            "rho_wall_ms": sum(step["arms"]["rho"]["wall_ms"] for step in rows),
        }
    save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments)
    classification = json.loads((arguments.out / "panel.json").read_text())["classification"]
    if classification.startswith("INVALID") or "INVALID" in classification:
        raise SystemExit(f"independent validation failed: {classification}")
