#!/usr/bin/env python3
"""Run one preregistered cold point-only signed-Frobenius rho block."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import threading
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PROTOCOL = HERE / "protocol.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args):
    protocol = json.loads(PROTOCOL.read_text())
    assert [args.seed, args.targets] in protocol["run_order"]
    assert args.timeout == protocol["external_timeout_s"]
    source = REPO / "examples/koblitz_rho_batch_ks.rs"
    assert sha(source) == protocol["algorithm"]["source_sha256"]
    points = HERE / f"target_points_{args.targets}.jsonl"
    assert sha(points) == protocol["input"][f"point_file_{args.targets}_sha256"]
    assert len(points.read_bytes().splitlines()) == args.targets
    assert not args.out.exists(), "Each cold run needs a new output directory"
    args.out.mkdir(parents=True)
    command = [str(args.exe), "53", "0", "signed_frobenius",
               str(args.targets), str(args.seed)]
    env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
    env.update({
        "KIC_RHO_TARGET_POINTS_JSONL": str(points.resolve()),
        "KIC_RHO_BATCH_CORPUS": protocol["algorithm"]["corpus"],
        "KIC_RHO_DP_BITS": str(protocol["algorithm"]["dp_bits"]),
        "KIC_RHO_PRECOMPUTE_WALKS": "0",
    })
    manifest = {
        "schema_version": "1.0",
        "command": command,
        "environment": {key: value for key, value in env.items() if key.startswith("KIC_")},
        "protocol_sha256": sha(PROTOCOL),
        "source_sha256": sha(source),
        "executable_sha256": sha(args.exe),
        "point_file_sha256": sha(points),
        "seed": args.seed,
        "targets": args.targets,
        "host": platform.platform(),
        "machine": platform.machine(),
        "timeout_s": args.timeout,
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    before_load = os.getloadavg()
    started = time.perf_counter()
    timed_out = threading.Event()
    with (args.out / "rho.stdout.jsonl").open("w") as stdout, (args.out / "rho.stderr.txt").open("w") as stderr:
        child = subprocess.Popen(command, cwd=REPO, env=env, stdout=stdout, stderr=stderr)
        def kill_on_timeout():
            timed_out.set()
            try:
                child.kill()
            except ProcessLookupError:
                pass
        timer = threading.Timer(args.timeout, kill_on_timeout)
        timer.start()
        try:
            while True:
                try:
                    pid, status, usage = os.wait4(child.pid, 0)
                    break
                except InterruptedError:
                    continue
        finally:
            timer.cancel()
    child.returncode = os.waitstatus_to_exitcode(status)
    wall_ms = (time.perf_counter() - started) * 1000
    raw = (args.out / "rho.stdout.jsonl").read_bytes()
    rows = [json.loads(line) for line in raw.splitlines() if line]
    fixture_rows = [row for row in rows if row.get("kind") == "rho_ks_batch_fixture"]
    summaries = [row for row in rows if row.get("kind") == "rho_ks_batch_summary"]
    complete = (child.returncode == 0 and not timed_out.is_set()
                and len(fixture_rows) == args.targets and len(summaries) == 1)
    if complete:
        assert summaries[0]["all_verified"] and summaries[0]["target_source"] == "explicit_public_points"
        for index, row in enumerate(fixture_rows):
            assert row["fixture_index"] == index
            assert row["target_source"] == "explicit_public_points"
            assert row["published_fixture_scalar"] is None
            assert row["published_q"] == json.loads(points.read_text().splitlines()[index])
    receipt = {
        "schema_version": "1.0",
        "seed": args.seed,
        "targets": args.targets,
        "returncode": child.returncode,
        "timed_out": timed_out.is_set(),
        "complete_producer_output": complete,
        "process_wall_ms": wall_ms,
        "child_user_cpu_ms": usage.ru_utime * 1000,
        "child_system_cpu_ms": usage.ru_stime * 1000,
        "child_peak_rss_raw": usage.ru_maxrss,
        "rss_unit": "bytes" if platform.system() == "Darwin" else "kibibytes",
        "load_average_before": before_load,
        "load_average_after": os.getloadavg(),
        "stdout_sha256": sha(args.out / "rho.stdout.jsonl"),
        "stderr_sha256": sha(args.out / "rho.stderr.txt"),
        "fixtures_reported": len(fixture_rows),
        "summary_reported": len(summaries) == 1,
    }
    (args.out / "resource_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    if not complete:
        raise SystemExit("Rho run failed, timed out, or returned incomplete output; receipt preserved")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--targets", type=int, choices=[1, 8], required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--exe", type=Path, default=REPO / "target/release/examples/koblitz_rho_batch_ks")
    parser.add_argument("--timeout", type=int, default=180)
    run(parser.parse_args())
