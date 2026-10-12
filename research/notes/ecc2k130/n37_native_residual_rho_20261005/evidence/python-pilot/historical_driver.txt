#!/usr/bin/env python3
"""Run the frozen five-block n37 native IC versus strong batch-rho panel."""

import argparse
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
ROOT = HERE.parents[3]
POINTS = ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b01.points.jsonl"
IC = ROOT / "target/release/examples/n37_native_m6_residual"
REPLAY = ROOT / "target/release/examples/n37_native_m6_residual_replay"
RHO = ROOT / "target/release/examples/koblitz_rho_batch_ks_v3"
WALL_CAP_SECONDS = 180
RSS_CAP_BYTES = 4 * 1024**3
CORPUS = "compact-disjoint-cold-v2-n37-L1024-b01-20261001"
SCHEDULE = [["ic", "rho", "rho", "ic"] if block % 2 == 0
            else ["rho", "ic", "ic", "rho"] for block in range(5)]
EXPECTED = {
    "examples/n37_native_m6_residual.rs": "086a6fbdd12e2b104ec18f809d4048453d747ad76c5a74f695befd08d4efffa8",
    "examples/n37_native_m6_residual_replay.rs": "68a4ce03029b82821fdabb9854a14559dde243bd5ac0f409de884df032f7394a",
    "examples/koblitz_rho_batch_ks_v3.rs": "98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c",
    "Cargo.toml": "f88ac8c9527bc2403cc550b8b00aaa778630911a4c9d1de47f67dbbff055d6a2",
    "Cargo.lock": "b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365",
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b01.points.jsonl": "78553fdff5ae66521d3a1052978962258e78d48c92285df6e1df962034b43717",
}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def host_record():
    cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                         capture_output=True, text=True, check=False)
    return {"platform": platform.platform(), "machine": platform.machine(),
            "python": platform.python_version(), "cpu_model": cpu.stdout.strip(),
            "cpu_model_status": cpu.returncode, "cpu_count": os.cpu_count(),
            "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
            "cargo": subprocess.check_output(["cargo", "--version"], text=True).strip(),
            "isolation": "unverified_unisolated_local_host"}


def run_arm(block, position, method, out):
    stem = f"b{block}_{position}_{method}"
    stdout_path = out / f"{stem}.stdout"
    stderr_path = out / f"{stem}.stderr"
    output_path = out / f"{stem}.result.json" if method == "ic" else stdout_path
    if method == "ic":
        command = [str(IC), str(output_path)]
        environment = {}
    else:
        command = [str(RHO), "37", "0", "signed_frobenius", "1024", "2026100110101"]
        environment = {
            "KIC_RHO_POINT_INPUT": str(POINTS),
            "KIC_RHO_BATCH_CORPUS": CORPUS,
            "KIC_RHO_CANON_BACKEND": "normal_basis",
            "KIC_RHO_DP_BITS": "8",
        }
    env = dict(os.environ)
    env.update(environment)
    before_load = os.getloadavg()
    start = time.perf_counter_ns()
    timed_out = False
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        child = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stdout,
                                 stderr=stderr, start_new_session=True)
        while True:
            pid, status, usage = os.wait4(child.pid, os.WNOHANG)
            if pid:
                child.returncode = os.waitstatus_to_exitcode(status)
                break
            if (time.perf_counter_ns() - start) / 1e9 > WALL_CAP_SECONDS:
                timed_out = True
                os.killpg(child.pid, signal.SIGKILL)
                _, status, usage = os.wait4(child.pid, 0)
                child.returncode = os.waitstatus_to_exitcode(status)
                break
            time.sleep(0.01)
    stop = time.perf_counter_ns()
    peak = usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024
    row = {
        "schema": "n37-native-batch-rho-arm-v1", "block": block,
        "position": position, "method": method,
        "command": command, "environment": environment,
        "exit_code": child.returncode, "timed_out": timed_out,
        "wall_ns": stop - start,
        "user_ns": round(usage.ru_utime * 1e9),
        "system_ns": round(usage.ru_stime * 1e9),
        "peak_rss_bytes": peak, "loadavg_before": before_load,
        "loadavg_after": os.getloadavg(),
        "stdout": stdout_path.name, "stdout_sha256": digest(stdout_path),
        "stderr": stderr_path.name, "stderr_sha256": digest(stderr_path),
        "result": output_path.name if output_path.exists() else None,
        "result_sha256": digest(output_path) if output_path.exists() else None,
    }
    row["resource_accepted"] = (child.returncode == 0 and not timed_out and
                                peak <= RSS_CAP_BYTES)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists()
    # The protocol commit is the first allowed source snapshot. Subsequent
    # commits may only add this runner or evidence; exact executed files are
    # checked below regardless of branch-head movement.
    assert subprocess.run(["git", "merge-base", "--is-ancestor", "e8c01bf5f", "HEAD"],
                          cwd=ROOT, check=False).returncode == 0
    for relative, expected in EXPECTED.items():
        assert digest(ROOT / relative) == expected, relative
    assert all(path.is_file() for path in (IC, REPLAY, RHO))
    args.out.mkdir(parents=True)
    manifest = {
        "schema": "n37-native-batch-rho-panel-v1",
        "status": "running", "schedule": SCHEDULE,
        "wall_cap_seconds": WALL_CAP_SECONDS, "rss_cap_bytes": RSS_CAP_BYTES,
        "host": host_record(),
        "source_sha256": EXPECTED,
        "protocol_sha256": digest(HERE / "PROTOCOL.md"),
        "runner_sha256": digest(Path(__file__)),
        "binary_sha256": {"ic": digest(IC), "ic_replay": digest(REPLAY),
                          "rho": digest(RHO)},
        "point_file_sha256": digest(POINTS),
        "fixture_read_by_producers": False,
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2,
                                                       sort_keys=True) + "\n")
    with (args.out / "arms.jsonl").open("x", encoding="utf-8") as stream:
        for block, order in enumerate(SCHEDULE):
            for position, method in enumerate(order):
                row = run_arm(block, position, method, args.out)
                stream.write(json.dumps(row, sort_keys=True) + "\n")
                stream.flush()
                os.fsync(stream.fileno())
                print(json.dumps({key: row[key] for key in
                                  ("block", "position", "method", "exit_code",
                                   "timed_out", "wall_ns", "resource_accepted")}),
                      flush=True)
    manifest["status"] = "completed_schedule"
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2,
                                                       sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
