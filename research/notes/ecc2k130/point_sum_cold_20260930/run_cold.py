#!/usr/bin/env python3
"""Fresh-process, same-Q point-query/S3/rho cold panel with wait4 CPU."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from run_panel import run_child, sha, write_json  # noqa: E402
from prepare import CELLS, COMPACT_SHA, RHO_SHA, SOURCE_COMMIT  # noqa: E402

FROZEN = HERE / "FROZEN.json"
RUN_CHILD_SHA = "bb57f12a4f57b6921674af026a5993c661751bdcbc47d0984b6f32af64b014b8"
LIMIT_SECONDS = 900
LIMIT_BYTES = 5 * 1024**3
ARM_ORDER = ("control_a", "point_sum", "rho", "control_b")


def checked_spec(cell: str) -> dict:
    frozen = json.loads(FROZEN.read_text())
    assert frozen["schema"] == "ecc2k130-point-sum-cold-freeze-v1"
    assert frozen["source_commit"] == SOURCE_COMMIT
    specs = frozen["specs"]
    assert set(specs) == {f"n{n}_L{length}" for n, length, *_ in CELLS}
    spec = specs[cell]
    assert cell == f"n{spec['n']}_L{spec['L']}"
    assert any((spec["n"], spec["L"], spec["K"], spec["prefilter"], spec["blocks"])
               == item for item in CELLS)
    points = HERE / spec["points_file"]
    assert sha(points) == spec["points_sha256"]
    assert len(points.read_text().splitlines()) == spec["L"]
    return spec


def checked_source(source: Path, compact: Path, rho: Path) -> dict:
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source,
                                   text=True).strip()
    assert head == SOURCE_COMMIT, "build must come from the exact merged source"
    assert sha(source / "examples/koblitz_orbit_dlp_s3_batch.rs") == COMPACT_SHA
    assert sha(source / "examples/koblitz_rho_batch_ks_v3.rs") == RHO_SHA
    lock = HERE / "CARGO_LOCK_SNAPSHOT.txt"
    assert sha(source / "Cargo.lock") == sha(lock)
    assert sha(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930/run_panel.py") == RUN_CHILD_SHA
    assert compact.is_file() and rho.is_file()
    return {"commit": head, "compact_source_sha256": COMPACT_SHA,
            "rho_source_sha256": RHO_SHA, "cargo_lock_sha256": sha(lock),
            "compact_binary_sha256": sha(compact), "rho_binary_sha256": sha(rho),
            "rustc_version_verbose": subprocess.check_output(
                ["rustc", "--version", "--verbose"], text=True).strip(),
            "cargo_version": subprocess.check_output(
                ["cargo", "--version"], text=True).strip()}


def schedule(blocks: int, mode: str) -> list[tuple[int, str]]:
    if mode == "smoke":
        return [(0, "point_sum")]
    assert mode == "measure"
    return [(block, arm) for block in range(blocks)
            for arm in ARM_ORDER[block % 4:] + ARM_ORDER[:block % 4]]


def run(cell: str, source: Path, compact: Path, rho: Path, out: Path,
        cpu: int | None, mode: str) -> dict:
    spec = checked_spec(cell)
    source_info = checked_source(source, compact, rho)
    assert not out.exists(), "never overwrite a cold-process run"
    if mode == "measure":
        assert sys.platform == "linux" and cpu is not None
        assert cpu in os.sched_getaffinity(0)
        assert platform.machine() == "x86_64"
    plan = schedule(spec["blocks"], mode)
    out.mkdir(parents=True)
    report = {
        "schema": "ecc2k130-point-sum-cold-run-v1",
        "status": "RUNNING", "mode": mode, "cell": cell,
        "spec": spec, "frozen_sha256": sha(FROZEN),
        "runner_sha256": sha(Path(__file__)),
        "run_child_source_sha256": sha(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930/run_panel.py"),
        "source": source_info,
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "python": sys.version, "reserved_cpu": cpu,
                 "processor": platform.processor()},
        "limits": {"arm_timeout_seconds": LIMIT_SECONDS,
                   "arm_rss_bytes": LIMIT_BYTES, "arm_address_space_bytes": LIMIT_BYTES},
        "plan": [{"block": block, "arm": arm} for block, arm in plan],
        "runs": [],
    }
    write_json(out / "cold_run.json", report)
    points = HERE / spec["points_file"]
    try:
        for block, arm in plan:
            prefix = out / f"b{block}_{arm}"
            env = {key: value for key, value in os.environ.items()
                   if not key.startswith("KIC_")}
            env.update({"RAYON_NUM_THREADS": "1", "LC_ALL": "C"})
            if arm == "rho":
                env.update({"KIC_RHO_POINT_INPUT": str(points.resolve()),
                            "KIC_RHO_BATCH_CORPUS": spec["corpus"],
                            "KIC_RHO_DP_BITS": "4",
                            "KIC_RHO_CANON_BACKEND": "normal_basis"})
                program = [str(rho), str(spec["n"]), "0", "signed_frobenius",
                           str(spec["L"]), str(spec["seed"])]
            else:
                env.update({"KIC_DUMP_BASE": str(prefix.with_suffix(".base.jsonl")),
                            "KIC_DUMP_RANK": str(prefix.with_suffix(".rank.jsonl")),
                            "KIC_S3_BATCH_WINDOW": "64",
                            "KIC_S3_PREFILTER": spec["prefilter"],
                            "KIC_QUERY_BACKEND": "point_sum" if arm == "point_sum" else "s3"})
                program = [str(compact),
                           f"construct:{spec['n']}:0:{spec['K']}",
                           str(points.resolve()), "7",
                           str(prefix.with_suffix(".target.jsonl"))]
            command = (["taskset", "-c", str(cpu), *program]
                       if mode == "measure" else program)
            assert not any(".fixture.jsonl" in value for value in
                           command + list(env.values()))
            stdout = prefix.with_suffix(".stdout.jsonl")
            stderr = prefix.with_suffix(".stderr.txt")
            result = run_child(command, env, stdout, stderr, LIMIT_SECONDS,
                               LIMIT_BYTES, LIMIT_BYTES if sys.platform == "linux" else None)
            files = {"stdout": stdout, "stderr": stderr}
            if arm != "rho":
                files.update({"base": prefix.with_suffix(".base.jsonl"),
                              "rank": prefix.with_suffix(".rank.jsonl"),
                              "targets": prefix.with_suffix(".target.jsonl")})
            record = {
                "block": block, "arm": arm, "command": command,
                "environment": {key: env[key] for key in sorted(env)
                                if key.startswith("KIC_") or key in ("RAYON_NUM_THREADS", "LC_ALL")},
                "files": {name: {"name": path.name, "sha256": sha(path),
                                 "bytes": path.stat().st_size}
                          for name, path in files.items() if path.exists()},
                **result,
            }
            report["runs"].append(record)
            write_json(out / "cold_run.json", report)
            print(json.dumps({"cell": cell, "block": block, "arm": arm,
                              "exit_code": result["exit_code"],
                              "stopped_for": result["stopped_for"]}), flush=True)
            if result["exit_code"] != 0 or result["stopped_for"] is not None:
                report["status"] = "FAIL"
                break
        else:
            report["status"] = "SMOKE_PASS" if mode == "smoke" else "PASS"
    except BaseException as error:
        report["status"] = "FAIL"
        report["error_type"] = type(error).__name__
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
        write_json(out / "cold_run.json", report)
        raise
    write_json(out / "cold_run.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--compact", type=Path, required=True)
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--mode", choices=("smoke", "measure"), required=True)
    args = parser.parse_args()
    report = run(args.cell, args.source_root.resolve(), args.compact.resolve(),
                 args.rho.resolve(), args.out.resolve(), args.cpu, args.mode)
    if report["status"] not in ("PASS", "SMOKE_PASS"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
