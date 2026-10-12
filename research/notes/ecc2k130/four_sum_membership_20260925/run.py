#!/usr/bin/env python3
"""Run one preregistered four-sum oracle or unchanged compact extractor arm."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import psutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SWEEP = HERE.parent / "compact_base_sweep_20260925"
COMPACT = REPO / "target/release/examples/koblitz_s5_sat_instance"
ORACLE = REPO / "target/release/examples/koblitz_four_sum_membership"
ORACLE_SOURCE_SHA = "6185018d17e9541245456ffbfec2df36630fe5a4d7e47e8f638ff59f072bdd8a"
COMPACT_SOURCE_SHA = "c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09"
BASE_HASHES = {
    (37, 3): "2722f3c7271ea410e47ff9d20b67496a9a28120804f8b790e0729d8c86ea6ddb",
    (41, 8): "cbdd871504ebac2f36e16ce652bc6dfe9623d3fd62723b05f7112d3d94614520",
    (41, 12): "a3856771a5ea1ad3a1e2e262eebcc1b2ccde8293e53af4360268cdcf72ab7fa8",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def run(args: argparse.Namespace) -> None:
    n, r = args.n, args.r
    assert (n, r) in BASE_HASHES
    count = 128 if (n, r) == (41, 12) else 512
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    base_file = HERE / f"base_n{n}_R{r}.json"
    full_targets_file = HERE / f"target_points_n{n}.jsonl"
    targets_file = (HERE / "target_points_n41_R12.jsonl"
                    if (n, r) == (41, 12) else full_targets_file)
    header = json.loads(base_file.read_bytes())
    assert sha(json.dumps(header["factor_base_point_coordinates"],
                          separators=(",", ":")).encode()) == BASE_HASHES[(n, r)]
    full_target_lines = full_targets_file.read_bytes().splitlines(keepends=True)
    assert len(full_target_lines) == 512
    assert len(targets_file.read_bytes().splitlines()) == count
    if (n, r) == (41, 12):
        assert targets_file.read_bytes() == b"".join(full_target_lines[:128])
        assert sha(targets_file.read_bytes()) == "1b154bdd8aa9dabdb37d2dd5a7bfee69a1fa8284ac5db678ef73eaf2c2f297a5"
    if args.mode == "oracle":
        exe = ORACLE
        source = REPO / "examples/koblitz_four_sum_membership.rs"
        expected_source_sha = ORACLE_SOURCE_SHA
        command = [str(exe), str(base_file), str(targets_file), str(count)]
        env = {k: v for k, v in os.environ.items() if not k.startswith("KIC_")}
        env["RAYON_NUM_THREADS"] = "1"
    else:
        exe = COMPACT
        source = REPO / "examples/koblitz_s5_sat_instance.rs"
        expected_source_sha = COMPACT_SOURCE_SHA
        spec = json.loads((SWEEP / "input_spec.json").read_bytes())
        arm = next(a for a in spec["arms"][str(n)] if a["R"] == r)
        command = [str(exe), str(n), "0", str(arm["eta_numerator"]),
                   str(arm["eta_denominator"]), "natural", "1", "2000", "1", "internal"]
        env = {k: v for k, v in os.environ.items() if not k.startswith("KIC_")}
        env.update({"RAYON_NUM_THREADS": "1", "KIC_ALGEBRA_ENCODING": "orbit_factorized",
                    "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
                    "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
                    "KIC_ORBIT_REP_ENCODING": "one_hot",
                    "KIC_ORBIT_BATCH_ONLY": "1",
                    "KIC_ORBIT_INCLUDE_BASE_HEADER": "1",
                    "KIC_ORBIT_TARGET_POINTS_JSONL": str(targets_file),
                    "KIC_TASK_ID": "TASK-IC-COMPACT-FOUR-SUM-ORACLE-20260925"})
    assert sha(source.read_bytes()) == expected_source_sha
    manifest = {"schema_version": "1.0", "mode": args.mode,
                "n": n, "R": r, "count": count,
                "source_pr747_commit": "fc27150df3238b6863ed5618c721e7fd8b6ce403",
                "source_input_freeze_commit": "12d2e489de77dc633f02df47de2debf865a3f301",
                "budget_amendment_commit": "4d1b48634d559e38369e85936f5c7a2190e24323",
                "runner_sha256": sha(Path(__file__).read_bytes()),
                "checkout_commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                           cwd=REPO, text=True).strip(),
                "protocol_sha256": sha((HERE / "PROTOCOL.md").read_bytes()),
                "input_manifest_sha256": sha((HERE / "input_manifest.json").read_bytes()),
                "base_file": str(base_file.relative_to(REPO)), "base_file_sha256": sha(base_file.read_bytes()),
                "point_set_sha256": BASE_HASHES[(n, r)],
                "targets_file": str(targets_file.relative_to(REPO)),
                "targets_sha256": sha(targets_file.read_bytes()),
                "full_target_stream_sha256": sha(full_targets_file.read_bytes()),
                "input_amendment_sha256": (sha((HERE / "input_amendment.json").read_bytes())
                                           if (n, r) == (41, 12) else None),
                "source_path": str(source.relative_to(REPO)),
                "source_sha256": expected_source_sha,
                "executable_sha256": sha(exe.read_bytes()),
                "command": command,
                "environment": {k: env[k] for k in sorted(env)
                                if k.startswith("KIC_") or k == "RAYON_NUM_THREADS"},
                "timeout_seconds": 900, "sampled_rss_stop_bytes": 2 * 1024**3,
                "rss_poll_interval_ms": 50, "psutil_version": psutil.__version__,
                "host": platform.node(), "platform": platform.platform()}
    write(out / "manifest.json", manifest)
    start = time.monotonic()
    timed_out = False
    rss_stop = False
    peak_sampled_rss = 0
    with (out / "producer.stdout.jsonl").open("wb") as stdout, \
         (out / "producer.stderr.txt").open("wb") as stderr:
        process = subprocess.Popen(command, cwd=REPO, env=env, stdout=stdout,
                                   stderr=stderr)
        child = psutil.Process(process.pid)
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid == process.pid:
                break
            try:
                sampled_rss = child.memory_info().rss
                peak_sampled_rss = max(peak_sampled_rss, sampled_rss)
            except psutil.NoSuchProcess:
                sampled_rss = 0  # wait4 will retrieve the exited process next loop
            if time.monotonic() - start >= 900 or sampled_rss > 2 * 1024**3:
                timed_out = time.monotonic() - start >= 900
                rss_stop = sampled_rss > 2 * 1024**3
                process.kill()
                pid, status, usage = os.wait4(process.pid, 0)
                assert pid == process.pid
                break
            time.sleep(0.05)
    returncode = os.waitstatus_to_exitcode(status)
    peak_rss = int(usage.ru_maxrss * (1 if platform.system() == "Darwin" else 1024))
    receipt = {"schema_version": "1.0", "mode": args.mode, "n": n, "R": r,
               "returncode": returncode, "timed_out": timed_out,
               "sampled_rss_stop": rss_stop, "peak_sampled_rss_bytes": peak_sampled_rss,
               "wall_ms": (time.monotonic() - start) * 1000,
               "user_cpu_s": usage.ru_utime, "system_cpu_s": usage.ru_stime,
               "peak_rss_bytes": peak_rss,
               "stdout_sha256": sha((out / "producer.stdout.jsonl").read_bytes()),
               "stderr_sha256": sha((out / "producer.stderr.txt").read_bytes()),
               "manifest_sha256": sha((out / "manifest.json").read_bytes())}
    write(out / "receipt.json", receipt)
    if returncode == 0:
        rows = [json.loads(line) for line in (out / "producer.stdout.jsonl").read_text().splitlines()]
        assert len(rows) == (count + 1 if args.mode == "oracle" else 1)
        if args.mode == "oracle":
            assert rows[0]["unique_pair_sums"] > 0
            assert all(row["unique_sum_probes"] == rows[0]["unique_pair_sums"] for row in rows[1:])
        else:
            assert rows[0]["compact_orbit_point_batch"]["targets_requested"] == count
            observed_points = rows[0]["compact_orbit_base_header"]["factor_base_point_coordinates"]
            assert sha(json.dumps(observed_points, separators=(",", ":")).encode()) == BASE_HASHES[(n, r)]
    print(json.dumps({"out": str(out), **receipt}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("oracle", "extractor"), required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--r", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    run(parser.parse_args())
