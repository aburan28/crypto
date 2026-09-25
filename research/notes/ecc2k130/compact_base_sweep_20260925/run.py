#!/usr/bin/env python3
"""Run one fresh, charged compact or same-Q rho arm with immutable receipts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import resource
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SPEC = HERE / "input_spec.json"
COMPACT = REPO / "target/release/examples/koblitz_s5_sat_instance"
RHO = REPO / "target/release/examples/koblitz_rho_fixture"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def run(args: argparse.Namespace) -> None:
    spec_bytes = SPEC.read_bytes()
    spec = json.loads(spec_bytes)
    n = args.n
    assert n in (37, 41)
    assert args.mode == "rho" or args.r in [arm["R"] for arm in spec["arms"][str(n)]]
    if args.mode == "train":
        assert args.count in (64, 512, 4096)
    else:
        assert args.count is None
    if args.mode == "rho":
        assert args.seed_index in (0, 1, 2)
    else:
        assert args.seed_index is None
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = {k: v for k, v in os.environ.items() if not k.startswith("KIC_")}
    env["RAYON_NUM_THREADS"] = "1"
    input_file = None
    input_hash = None
    if args.mode == "rho":
        exe = RHO
        target = spec["holdouts"][str(n)][args.seed_index]
        command = [str(exe), str(n), "0", "signed_frobenius", "1", "packed",
                   "13737", f"hash:{target['seed']}"]
        source = REPO / "examples/koblitz_rho_fixture.rs"
    else:
        exe = COMPACT
        arm = next(arm for arm in spec["arms"][str(n)] if arm["R"] == args.r)
        command = [str(exe), str(n), "0", str(arm["eta_numerator"]),
                   str(arm["eta_denominator"]), "natural", "1", "2000", "1", "internal"]
        source = REPO / "examples/koblitz_s5_sat_instance.rs"
        env.update({"KIC_ALGEBRA_ENCODING": "orbit_factorized",
                    "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
                    "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
                    "KIC_ORBIT_REP_ENCODING": "one_hot",
                    "KIC_ORBIT_BATCH_ONLY": "1",
                    "KIC_TASK_ID": "TASK-IC-COMPACT-BASE-SWEEP-20260925"})
        if args.mode == "train":
            source_lines = (HERE / f"training_scalars_n{n}.txt").read_bytes().splitlines()
            input_data = b"\n".join(source_lines[:args.count]) + b"\n"
            input_file = out / "target_scalars.txt"
            env["KIC_ORBIT_TARGET_SCALARS"] = str(input_file)
        else:
            input_data = (HERE / f"holdout_points_n{n}.jsonl").read_bytes()
            input_file = out / "target_points.jsonl"
            env["KIC_ORBIT_TARGET_POINTS_JSONL"] = str(input_file)
        input_file.write_bytes(input_data)
        input_hash = sha(input_data)
    source_commit = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=REPO, text=True).strip()
    manifest = {"schema_version": "1.0", "n": n, "mode": args.mode,
                "R": args.r, "count": args.count, "seed_index": args.seed_index,
                "input_spec_sha256": sha(spec_bytes),
                "source_commit": source_commit,
                "source_path": str(source.relative_to(REPO)),
                "source_sha256": sha(source.read_bytes()),
                "executable_sha256": sha(exe.read_bytes()),
                "input_file": input_file.name if input_file else None,
                "input_sha256": input_hash,
                "command": command,
                "environment": {k: env[k] for k in sorted(env) if k.startswith("KIC_")
                                or k == "RAYON_NUM_THREADS"},
                "timeout_seconds": args.timeout,
                "host": platform.node(), "platform": platform.platform()}
    write_json(out / "manifest.json", manifest)
    stdout = out / "producer.stdout.jsonl"
    stderr = out / "producer.stderr.txt"
    started = time.monotonic()
    with stdout.open("wb") as out_stream, stderr.open("wb") as err_stream:
        process = subprocess.Popen(command, cwd=REPO, env=env,
                                   stdout=out_stream, stderr=err_stream)
        deadline = started + args.timeout
        timed_out = False
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid == process.pid:
                break
            if time.monotonic() >= deadline:
                timed_out = True
                process.kill()
                pid, status, usage = os.wait4(process.pid, 0)
                assert pid == process.pid
                break
            time.sleep(0.05)
        process.returncode = os.waitstatus_to_exitcode(status)
    wall_ms = (time.monotonic() - started) * 1000
    peak_rss = int(usage.ru_maxrss * (1 if platform.system() == "Darwin" else 1024))
    receipt = {"schema_version": "1.0", "returncode": process.returncode,
               "timed_out": timed_out, "wall_ms": wall_ms,
               "user_cpu_s": usage.ru_utime, "system_cpu_s": usage.ru_stime,
               "peak_rss_bytes": peak_rss,
               "stdout_sha256": sha(stdout.read_bytes()),
               "stderr_sha256": sha(stderr.read_bytes()),
               "manifest_sha256": sha((out / "manifest.json").read_bytes())}
    write_json(out / "receipt.json", receipt)
    if process.returncode == 0:
        rows = [json.loads(line) for line in stdout.read_text().splitlines() if line.strip()]
        assert len(rows) == 1, len(rows)
        observed = rows[0]
        assert observed["n"] == n and observed["a"] == 0
        assert int(observed["subgroup_order"]) == int(spec["arms"][str(n)][0]["subgroup_order"])
        if args.mode == "rho":
            expected = spec["holdouts"][str(n)][args.seed_index]
            assert observed["published_q"] == expected["q"]
            assert observed["recovered_fixture_scalar"] == expected["scalar_validator_only"]
        else:
            assert observed["orbit_columns"] == args.r
            batch = observed["compact_orbit_batch" if args.mode == "train"
                             else "compact_orbit_point_batch"]
            assert batch["targets_requested"] == (args.count if args.mode == "train" else 3)
            assert len(batch["query_observations"]) == batch["targets_requested"]
            assert batch["sat_verification_included"] is False
    print(json.dumps({"out": str(out), **receipt}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("train", "holdout", "rho"), required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--r", type=int)
    parser.add_argument("--count", type=int)
    parser.add_argument("--seed-index", type=int)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    run(parser.parse_args())
