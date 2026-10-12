#!/usr/bin/env python3
"""Fresh-process orbit-disjoint v2 W64 compact/rho cold CPU panel."""
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
from prepare import (CELLS, COMPACT_SHA, LOCK_SHA, RHO_SHA, SOURCE_FREEZE,
                     SOURCE_FREEZE_SHA, cell_name)  # noqa: E402

FROZEN = HERE / "FROZEN.json"
INPUT_RECEIPT = HERE / "INPUT_RECEIPT.json"
HELPER_SHA = "bb57f12a4f57b6921674af026a5993c661751bdcbc47d0984b6f32af64b014b8"
LIMIT_SECONDS = 900
LIMIT_BYTES = 5 * 1024**3
ARM_ORDER = ("ic_a", "rho", "ic_b")


def checked_spec(cell: str) -> dict:
    frozen = json.loads(FROZEN.read_text())
    assert frozen["schema"] == "ecc2k130-disjoint-cold-v2-freeze-v1"
    input_receipt = json.loads(INPUT_RECEIPT.read_text())
    assert input_receipt["status"] == "PASS"
    assert input_receipt["point_equations_verified"] == 15390
    assert input_receipt["frozen_sha256"] == sha(FROZEN)
    assert set(frozen["specs"]) == {cell_name(n, length) for n, length, *_ in CELLS}
    spec = frozen["specs"][cell]
    assert cell == cell_name(spec["n"], spec["L"])
    assert (spec["n"], spec["L"], spec["K"], spec["prefilter"], spec["blocks"]) in CELLS
    assert len(spec["block_specs"]) == spec["blocks"]
    for block, entry in enumerate(spec["block_specs"]):
        assert entry["block"] == block
        point_path = HERE / entry["points_file"]
        assert sha(point_path) == entry["points_sha256"]
        assert len(point_path.read_text().splitlines()) == spec["L"]
    return spec


def checked_source(source_root: Path, compact: Path, rho: Path,
                   materialization: Path) -> dict:
    assert sha(SOURCE_FREEZE) == SOURCE_FREEZE_SHA
    freeze = json.loads(SOURCE_FREEZE.read_text())
    receipt = json.loads(materialization.read_text())
    assert receipt["schema"] == "compact-frozen-source-materialization-v1"
    assert receipt["pinned_files"] == len(freeze["source_sha256"]) == 20
    assert receipt["freezes"] == [str(SOURCE_FREEZE.relative_to(ROOT))]
    for path, expected in freeze["source_sha256"].items():
        assert sha(source_root / path) == expected, path
    assert sha(source_root / "examples/koblitz_orbit_dlp_s3_batch.rs") == COMPACT_SHA
    assert sha(source_root / "examples/koblitz_rho_batch_ks_v3.rs") == RHO_SHA
    assert sha(source_root / "Cargo.lock") == LOCK_SHA
    assert sha(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930/run_panel.py") == HELPER_SHA
    assert compact.is_file() and rho.is_file()
    return {"source_freeze_sha256": SOURCE_FREEZE_SHA,
            "materialization_sha256": sha(materialization),
            "materialization": receipt,
            "compact_source_sha256": COMPACT_SHA,
            "rho_source_sha256": RHO_SHA,
            "cargo_lock_sha256": LOCK_SHA,
            "compact_binary_sha256": sha(compact),
            "rho_binary_sha256": sha(rho),
            "rustc_version_verbose": subprocess.check_output(
                ["rustc", "--version", "--verbose"], text=True).strip(),
            "cargo_version": subprocess.check_output(["cargo", "--version"], text=True).strip()}


def schedule(blocks: int, mode: str) -> list[tuple[int, str]]:
    assert blocks in (5, 20)
    assert mode in ("smoke", "measure")
    count = 1 if mode == "smoke" else blocks
    return [(block, arm) for block in range(count)
            for arm in ARM_ORDER[block % 3:] + ARM_ORDER[:block % 3]]


def host_info(cpu: int | None) -> dict:
    model = platform.processor()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    return {"platform": platform.platform(), "machine": platform.machine(),
            "python": sys.version, "reserved_cpu": cpu, "cpu_model": model,
            "cpuinfo_sha256": sha(cpuinfo) if cpuinfo.is_file() else None,
            "affinity": sorted(os.sched_getaffinity(0)) if sys.platform == "linux" else None,
            "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                cwd=ROOT, text=True).strip()}


def run(cell: str, source_root: Path, compact: Path, rho: Path,
        materialization: Path, out: Path, cpu: int | None, mode: str) -> dict:
    spec = checked_spec(cell)
    source = checked_source(source_root, compact, rho, materialization)
    assert not out.exists(), "never overwrite a cold-process run"
    if mode == "measure":
        assert sys.platform == "linux" and platform.machine() == "x86_64"
        assert cpu is not None and cpu in os.sched_getaffinity(0)
    plan = schedule(spec["blocks"], mode)
    out.mkdir(parents=True)
    (out / "materialization.json").write_bytes(materialization.read_bytes())
    report = {"schema": "ecc2k130-disjoint-cold-v2-run-v1",
              "status": "RUNNING", "mode": mode, "cell": cell,
              "spec": spec, "frozen_sha256": sha(FROZEN),
              "input_receipt_sha256": sha(INPUT_RECEIPT),
              "runner_sha256": sha(Path(__file__)),
              "run_child_source_sha256": HELPER_SHA,
              "source": source, "host": host_info(cpu),
              "limits": {"arm_timeout_seconds": LIMIT_SECONDS,
                         "arm_rss_bytes": LIMIT_BYTES,
                         "arm_address_space_bytes": LIMIT_BYTES},
              "plan": [{"block": b, "arm": arm} for b, arm in plan],
              "runs": []}
    write_json(out / "cold_run.json", report)
    try:
        for block, arm in plan:
            block_spec = spec["block_specs"][block]
            points = (HERE / block_spec["points_file"]).resolve()
            prefix = out / f"b{block:02d}_{arm}"
            env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
            env.update({"RAYON_NUM_THREADS": "1", "LC_ALL": "C"})
            if arm == "rho":
                env.update({"KIC_RHO_POINT_INPUT": str(points),
                            "KIC_RHO_BATCH_CORPUS": block_spec["corpus"],
                            "KIC_RHO_DP_BITS": "4",
                            "KIC_RHO_CANON_BACKEND": "normal_basis"})
                program = [str(rho), str(spec["n"]), "0", "signed_frobenius",
                           str(spec["L"]), str(block_spec["seed"])]
            else:
                env.update({"KIC_DUMP_BASE": str(prefix.with_suffix(".base.jsonl")),
                            "KIC_DUMP_RANK": str(prefix.with_suffix(".rank.jsonl")),
                            "KIC_S3_BATCH_WINDOW": "64",
                            "KIC_S3_PREFILTER": spec["prefilter"]})
                program = [str(compact), f"construct:{spec['n']}:0:{spec['K']}",
                           str(points), "7", str(prefix.with_suffix(".target.jsonl"))]
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
            record = {"block": block, "arm": arm,
                      "points_file": block_spec["points_file"],
                      "points_sha256": block_spec["points_sha256"],
                      "command": command,
                      "environment": {key: env[key] for key in sorted(env)
                                      if key.startswith("KIC_") or
                                      key in ("RAYON_NUM_THREADS", "LC_ALL")},
                      "files": {name: {"name": path.name, "sha256": sha(path),
                                       "bytes": path.stat().st_size}
                                for name, path in files.items() if path.exists()},
                      **result}
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
    parser.add_argument("--materialization", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--mode", choices=("smoke", "measure"), required=True)
    args = parser.parse_args()
    report = run(args.cell, args.source_root.resolve(), args.compact.resolve(),
                 args.rho.resolve(), args.materialization.resolve(),
                 args.out.resolve(), args.cpu, args.mode)
    if report["status"] not in ("PASS", "SMOKE_PASS"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
