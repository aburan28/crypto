#!/usr/bin/env python3
"""One-shot, fresh-process nine-arm n41 base-window/rho screen.

Run only after FROZEN.json and INPUT_RECEIPT.json have been reviewed and
merged. The runner consumes public point files; it never opens scalar fixtures.
Use tools/isolated_bench.py reserve around this command on the timing host.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from run_panel import rss_bytes, sha, write_json  # noqa: E402

FROZEN = HERE / "FROZEN.json"
INPUT_RECEIPT = HERE / "INPUT_RECEIPT.json"
SOURCE_FREEZE = ROOT / "research/notes/ecc2k130/compact_s3_prefilter_20260930/FROZEN.json"
SOURCE_FREEZE_SHA256 = "3e9f67cc2cd6de5a8458badb3525983d561d9c3449c05118819fb424c5093b2b"
FROZEN_LOCK_SHA256 = "7d671f48c2da93f133d98802e80f858d1d9ea3b86996f7037f758990e1566627"
FROZEN_FAST_ARITH_SHA256 = "2a5c54bdedb6a0b1ffcf95dd7badd32c41d1482b53055e31fb9f28d514ac8c9e"


def config() -> dict:
    result = json.loads((HERE / "CONFIG.json").read_text())
    assert result["schema"] == "ecc2k130-base-window-screen-protocol-v1"
    return result


def schedule(cfg: dict) -> list[tuple[int, str]]:
    arms = tuple(cfg["arm_order_before_rotation"])
    assert arms == ("w0_a", "w1_a", "w2_a", "w3_a", "rho",
                    "w0_b", "w1_b", "w2_b", "w3_b")
    return [
        (block, arm)
        for block in range(cfg["blocks"])
        for arm in arms[block % len(arms):] + arms[:block % len(arms)]
    ]


def checked_inputs(cfg: dict) -> dict:
    assert FROZEN.is_file() and INPUT_RECEIPT.is_file()
    frozen = json.loads(FROZEN.read_text())
    receipt = json.loads(INPUT_RECEIPT.read_text())
    assert frozen["schema"] == "ecc2k130-base-window-input-freeze-v1"
    assert receipt["status"] == "PASS"
    assert receipt["frozen_sha256"] == sha(FROZEN)
    assert receipt["point_equations_verified"] == cfg["blocks"] * cfg["public_targets_per_block"]
    assert frozen["curve_slug"] == cfg["curve_slug"]
    assert (frozen["n"], frozen["a"], frozen["subgroup_order"],
            frozen["cofactor"]) == (
                cfg["n"], cfg["a"], cfg["subgroup_order"], cfg["cofactor"]
            )
    assert frozen["protocol_config_sha256"] == sha(HERE / "CONFIG.json")
    assert frozen["protocol_sha256"] == sha(HERE / "PROTOCOL.md")
    assert len(frozen["blocks"]) == cfg["blocks"]
    for block, item in enumerate(frozen["blocks"]):
        assert item["block"] == block
        points = HERE / item["points_file"]
        assert points.is_file() and sha(points) == item["points_sha256"]
        assert len(points.read_bytes().splitlines()) == cfg["public_targets_per_block"]
    return frozen


def checked_source(cfg: dict, frozen: dict, source_root: Path,
                   generator: Path, compact: Path, rho: Path,
                   materialization: Path, build_receipt: Path) -> dict:
    assert sha(SOURCE_FREEZE) == SOURCE_FREEZE_SHA256
    source_freeze = json.loads(SOURCE_FREEZE.read_text())
    record = json.loads(materialization.read_text())
    assert record["schema"] == "compact-frozen-source-materialization-v1"
    assert record["pinned_files"] == len(source_freeze["source_sha256"]) == 20
    assert record["freezes"] == [str(SOURCE_FREEZE.relative_to(ROOT))]
    for name, expected in source_freeze["source_sha256"].items():
        assert sha(source_root / name) == expected, name
    assert sha(source_root / "examples/koblitz_orbit_dlp_s3_batch.rs") == cfg["compact_source_sha256"]
    assert sha(source_root / "examples/koblitz_rho_batch_ks_v3.rs") == cfg["rho_source_sha256"]
    assert sha(source_root / "src/cryptanalysis/koblitz_fast_arith.rs") == FROZEN_FAST_ARITH_SHA256
    assert sha(source_root / "Cargo.lock") == FROZEN_LOCK_SHA256
    source_map = frozen["source_lock"]["sha256"]
    for name, expected in source_map.items():
        assert sha(ROOT / name) == expected, name
    generator_source = "examples/koblitz_base_window.rs"
    assert sha(source_root / generator_source) == source_map[generator_source]
    assert cfg["compact_query_backend"] == "s3"  # The frozen executable has only S3.
    for binary in (generator, compact, rho):
        assert binary.is_file()
    built = json.loads(build_receipt.read_text())
    assert built["schema"] == "ecc2k130-base-window-frozen-build-v1"
    assert built["status"] == "PASS" and built["profile"] == "release"
    assert build_receipt.parent / "source" == source_root
    assert built["source_freeze_sha256"] == SOURCE_FREEZE_SHA256
    assert built["materialization_sha256"] == sha(materialization)
    assert built["generator_source_sha256"] == source_map[generator_source]
    assert built["compact_source_sha256"] == cfg["compact_source_sha256"]
    assert built["rho_source_sha256"] == cfg["rho_source_sha256"]
    assert built["fast_arith_sha256"] == FROZEN_FAST_ARITH_SHA256
    assert built["cargo_lock_sha256"] == FROZEN_LOCK_SHA256
    for name, binary in (
        ("koblitz_base_window", generator),
        ("koblitz_orbit_dlp_s3_batch", compact),
        ("koblitz_rho_batch_ks_v3", rho),
    ):
        assert Path(built["binaries"][name]["path"]) == binary
        assert built["binaries"][name]["sha256"] == sha(binary)
    return {
        "source_freeze_sha256": SOURCE_FREEZE_SHA256,
        "materialization_sha256": sha(materialization),
        "materialization": record,
        "build_receipt_sha256": sha(build_receipt),
        "generator_source_sha256": source_map[generator_source],
        "compact_source_sha256": cfg["compact_source_sha256"],
        "rho_source_sha256": cfg["rho_source_sha256"],
        "frozen_fast_arith_sha256": FROZEN_FAST_ARITH_SHA256,
        "cargo_lock_sha256": FROZEN_LOCK_SHA256,
        "generator_binary_sha256": sha(generator),
        "compact_binary_sha256": sha(compact),
        "rho_binary_sha256": sha(rho),
        "rustc_version_verbose": subprocess.check_output(
            ["rustc", "--version", "--verbose"], text=True
        ).strip(),
        "cargo_version": subprocess.check_output(["cargo", "--version"], text=True).strip(),
    }


def child(command: list[str], env: dict[str, str], stdout: Path, stderr: Path,
          timeout: float, rss_limit: int) -> dict:
    """wait4 CPU/RSS accounting with a separate process-launch wall receipt."""
    started = time.monotonic()
    stopped_for = None
    peak_rss = 0

    def limits() -> None:
        resource.setrlimit(resource.RLIMIT_AS, (rss_limit, rss_limit))

    with stdout.open("xb") as out, stderr.open("xb") as err:
        process = subprocess.Popen(
            command, stdout=out, stderr=err, env=env,
            start_new_session=True, preexec_fn=limits,
        )
        launch_wall = time.monotonic() - started
        while True:
            waited_pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if waited_pid:
                break
            peak_rss = max(peak_rss, rss_bytes(process.pid))
            if time.monotonic() - started > timeout:
                stopped_for = "timeout"
            elif peak_rss > rss_limit:
                stopped_for = "rss_limit"
            if stopped_for:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                waited_pid, status, usage = os.wait4(process.pid, 0)
                break
            time.sleep(0.2)
    assert waited_pid == process.pid
    return {
        "exit_code": os.waitstatus_to_exitcode(status),
        "stopped_for": stopped_for,
        "elapsed_wall_seconds_not_primary_cost": time.monotonic() - started,
        "process_launch_wall_seconds_not_primary_cost": launch_wall,
        "observed_peak_rss_bytes": peak_rss,
        "child_user_cpu_seconds": usage.ru_utime,
        "child_system_cpu_seconds": usage.ru_stime,
        "child_max_rss_kib_linux": usage.ru_maxrss,
    }


def file_meta(paths: dict[str, Path]) -> dict:
    return {
        name: {"name": path.name, "sha256": sha(path), "bytes": path.stat().st_size}
        for name, path in paths.items() if path.is_file()
    }


def host_info(cpu: int) -> dict:
    cpuinfo = Path("/proc/cpuinfo")
    model = platform.processor()
    if cpuinfo.is_file():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version,
        "reserved_cpu": cpu,
        "cpu_model": model,
        "cpuinfo_sha256": sha(cpuinfo),
        "affinity": sorted(os.sched_getaffinity(0)),
        "git_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
    }


def run(source_root: Path, generator: Path, compact: Path, rho: Path,
        materialization: Path, build_receipt: Path, out: Path, cpu: int) -> dict:
    cfg = config()
    assert sys.platform == "linux" and platform.machine() == "x86_64"
    assert cpu in os.sched_getaffinity(0)
    assert not out.exists(), "never overwrite a measured window screen"
    frozen = checked_inputs(cfg)
    source = checked_source(
        cfg, frozen, source_root, generator, compact, rho, materialization, build_receipt
    )
    plan = schedule(cfg)
    assert len(plan) == cfg["blocks"] * len(cfg["arm_order_before_rotation"]) == 45
    out.mkdir(parents=True)
    (out / "materialization.json").write_bytes(materialization.read_bytes())
    (out / "BUILD_RECEIPT.json").write_bytes(build_receipt.read_bytes())
    report = {
        "schema": "ecc2k130-base-window-screen-run-v1",
        "status": "RUNNING",
        "config": cfg,
        "frozen_sha256": sha(FROZEN),
        "input_receipt_sha256": sha(INPUT_RECEIPT),
        "runner_sha256": sha(Path(__file__)),
        "source": source,
        "host": host_info(cpu),
        "limits": {
            "generator_wall_seconds": cfg["generator_wall_limit_seconds"],
            "generator_rss_bytes": cfg["generator_rss_limit_bytes"],
            "arm_wall_seconds": cfg["arm_wall_limit_seconds"],
            "arm_rss_bytes": cfg["arm_rss_limit_bytes"],
            "cell_wall_seconds": cfg["cell_wall_limit_seconds"],
        },
        "plan": [{"block": block, "arm": arm} for block, arm in plan],
        "runs": [],
    }
    report_path = out / "screen_run.json"
    write_json(report_path, report)
    cell_started = time.monotonic()
    try:
        for block, arm in plan:
            remaining = cfg["cell_wall_limit_seconds"] - (time.monotonic() - cell_started)
            if remaining <= 0:
                report["status"] = "CENSORED"
                report["failure"] = {"reason": "cell_wall_limit", "block": block, "arm": arm}
                break
            block_spec = frozen["blocks"][block]
            points = (HERE / block_spec["points_file"]).resolve()
            prefix = out / f"b{block:02d}_{arm}"
            env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
            env.update({"RAYON_NUM_THREADS": "1", "LC_ALL": "C"})
            record: dict = {
                "block": block, "arm": arm,
                "points_file": block_spec["points_file"],
                "points_sha256": block_spec["points_sha256"],
                "children": {},
            }
            report["runs"].append(record)
            if arm == "rho":
                env.update({
                    "KIC_RHO_POINT_INPUT": str(points),
                    "KIC_RHO_BATCH_CORPUS": block_spec["corpus"],
                    "KIC_RHO_DP_BITS": str(cfg["rho_dp_bits"]),
                    "KIC_RHO_CANON_BACKEND": cfg["rho_canonical_backend"],
                })
                command = [
                    "taskset", "-c", str(cpu), str(rho),
                    str(cfg["n"]), str(cfg["a"]), "signed_frobenius",
                    str(cfg["public_targets_per_block"]), str(block_spec["rho_seed"]),
                ]
                stdout = prefix.with_suffix(".rho.stdout.jsonl")
                stderr = prefix.with_suffix(".rho.stderr.txt")
                result = child(
                    command, env, stdout, stderr,
                    min(cfg["arm_wall_limit_seconds"], remaining),
                    cfg["arm_rss_limit_bytes"],
                )
                record["children"]["rho"] = {
                    "command": command,
                    "environment": {k: env[k] for k in sorted(env)
                                    if k.startswith("KIC_") or k in ("LC_ALL", "RAYON_NUM_THREADS")},
                    "files": file_meta({"stdout": stdout, "stderr": stderr}),
                    **result,
                }
                failed = result["exit_code"] != 0 or result["stopped_for"] is not None
            else:
                window = int(arm[1])
                assert arm in (f"w{window}_a", f"w{window}_b") and window in range(4)
                header = prefix.with_suffix(".header.jsonl")
                scan = prefix.with_suffix(".scan.jsonl")
                receipt = prefix.with_suffix(".generator.json")
                gen_command = [
                    "taskset", "-c", str(cpu), str(generator),
                    str(cfg["n"]), str(cfg["a"]),
                    str(cfg["useful_orbit_columns"]), str(window),
                    str(cfg["raw_x_trial_cap"]), str(header), str(scan), str(receipt),
                ]
                gen_stdout = prefix.with_suffix(".generator.stdout.txt")
                gen_stderr = prefix.with_suffix(".generator.stderr.txt")
                gen = child(
                    gen_command, env, gen_stdout, gen_stderr,
                    min(cfg["generator_wall_limit_seconds"], remaining),
                    cfg["generator_rss_limit_bytes"],
                )
                record["children"]["generator"] = {
                    "command": gen_command,
                    "environment": {k: env[k] for k in sorted(env)
                                    if k in ("LC_ALL", "RAYON_NUM_THREADS")},
                    "files": file_meta({
                        "stdout": gen_stdout, "stderr": gen_stderr,
                        "header": header, "scan": scan, "receipt": receipt,
                    }),
                    **gen,
                }
                write_json(report_path, report)
                failed = gen["exit_code"] != 0 or gen["stopped_for"] is not None
                if not failed:
                    gen_receipt = json.loads(receipt.read_text())
                    assert gen_receipt["status"] == "complete"
                    assert gen_receipt["selected_orbits"] == cfg["useful_orbit_columns"]
                    remaining = cfg["cell_wall_limit_seconds"] - (time.monotonic() - cell_started)
                    if remaining <= 0:
                        report["status"] = "CENSORED"
                        report["failure"] = {
                            "reason": "cell_wall_limit_after_generator",
                            "block": block, "arm": arm,
                        }
                        break
                    env.update({
                        "KIC_DUMP_BASE": str(prefix.with_suffix(".base.jsonl")),
                        "KIC_DUMP_RANK": str(prefix.with_suffix(".rank.jsonl")),
                        "KIC_S3_BATCH_WINDOW": str(cfg["compact_batch_window"]),
                        "KIC_S3_PREFILTER": cfg["compact_prefilter"],
                    })
                    target = prefix.with_suffix(".target.jsonl")
                    compact_command = [
                        "taskset", "-c", str(cpu), str(compact),
                        str(header), str(points), str(cfg["rank_seed"]), str(target),
                    ]
                    compact_stdout = prefix.with_suffix(".compact.stdout.jsonl")
                    compact_stderr = prefix.with_suffix(".compact.stderr.txt")
                    compact_result = child(
                        compact_command, env, compact_stdout, compact_stderr,
                        min(cfg["arm_wall_limit_seconds"], remaining),
                        cfg["arm_rss_limit_bytes"],
                    )
                    record["children"]["compact"] = {
                        "command": compact_command,
                        "environment": {k: env[k] for k in sorted(env)
                                        if k.startswith("KIC_") or k in ("LC_ALL", "RAYON_NUM_THREADS")},
                        "files": file_meta({
                            "stdout": compact_stdout, "stderr": compact_stderr,
                            "base": prefix.with_suffix(".base.jsonl"),
                            "rank": prefix.with_suffix(".rank.jsonl"),
                            "targets": target,
                        }),
                        **compact_result,
                    }
                    failed = (compact_result["exit_code"] != 0 or
                              compact_result["stopped_for"] is not None)
            write_json(report_path, report)
            print(json.dumps({
                "block": block, "arm": arm, "failed": failed,
                "children": sorted(record["children"]),
            }), flush=True)
            if failed:
                report["status"] = "CENSORED"
                report["failure"] = {"reason": "child_failure", "block": block, "arm": arm}
                break
        else:
            report["status"] = "PASS"
    except BaseException as error:
        report["status"] = "FAIL"
        report["error_type"] = type(error).__name__
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
        write_json(report_path, report)
        raise
    write_json(report_path, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--compact", required=True, type=Path)
    parser.add_argument("--rho", required=True, type=Path)
    parser.add_argument("--materialization", required=True, type=Path)
    parser.add_argument("--build-receipt", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    args = parser.parse_args()
    result = run(
        args.source_root.resolve(), args.generator.resolve(), args.compact.resolve(),
        args.rho.resolve(), args.materialization.resolve(),
        args.build_receipt.resolve(), args.out.resolve(), args.cpu,
    )
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
