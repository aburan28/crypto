#!/usr/bin/env python3
"""Run cold point-only IC and signed-Frobenius batched rho on identical Q files."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
LIMIT_SECONDS = 900
TUNE_LIMIT_SECONDS = 180
LIMIT_BYTES = 5 * 1024**3


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def limits() -> None:
    if sys.platform == "linux":
        resource.setrlimit(resource.RLIMIT_AS, (LIMIT_BYTES, LIMIT_BYTES))


def child(command: list[str], env_update: dict[str, str], prefix: Path,
          limit_seconds: int = LIMIT_SECONDS) -> dict:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    stdout = prefix.with_suffix(".stdout")
    stderr = prefix.with_suffix(".stderr")
    assert not stdout.exists() and not stderr.exists()
    env = dict(os.environ)
    for name in ("KIC_RHO_GENERATE_ONLY", "KIC_RHO_POINT_INPUT",
                 "KIC_RHO_BATCH_CORPUS", "KIC_RHO_DP_BITS",
                 "KIC_DUMP_BASE", "KIC_DUMP_RANK"):
        env.pop(name, None)
    env.update(env_update)
    before = os.getloadavg()
    with stdout.open("wb") as output, stderr.open("wb") as error:
        started = time.monotonic()
        proc = subprocess.Popen(command, stdout=output, stderr=error, env=env,
                                cwd=prefix.parent,
                                preexec_fn=limits if sys.platform == "linux" else None)
        timed_out = False
        while True:
            pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if time.monotonic() - started >= limit_seconds:
                timed_out = True
                proc.send_signal(signal.SIGKILL)
                _, status, usage = os.wait4(proc.pid, 0)
                break
            time.sleep(0.005)
        proc.returncode = os.waitstatus_to_exitcode(status)
        wall = time.monotonic() - started
    return {
        "command": command, "environment": env_update, "exit_code": proc.returncode,
        "limit_seconds": limit_seconds,
        "timeout": timed_out, "wall_seconds": wall, "user_seconds": usage.ru_utime,
        "sys_seconds": usage.ru_stime, "max_rss_kib_linux": usage.ru_maxrss
        if sys.platform == "linux" else None,
        "load_before": before, "load_after": os.getloadavg(),
        "stdout": stdout.name, "stderr": stderr.name,
        "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr),
    }


def jsonlines(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def valid_tune(prefix: Path, expected_scalar: int, k: int) -> bool:
    summary_path = prefix.with_suffix(".stdout")
    target_path = prefix.with_suffix(".target.jsonl")
    if not summary_path.exists() or not target_path.exists():
        return False
    try:
        summary, = jsonlines(summary_path)
        target, = jsonlines(target_path)
    except (ValueError, json.JSONDecodeError):
        return False
    return (summary["kind"] == "compact_orbit_dlp_summary"
            and summary["rank"] == k
            and summary["targets_solved"] == 1
            and target["recovered_scalar"] == expected_scalar
            and target["group_verified"] is True
            and target["published_fixture_scalar"] is None)


def run(n: int, length: int, rho: Path, ic: Path, out: Path) -> None:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["tuning_timeout_seconds"] == TUNE_LIMIT_SECONDS
    assert frozen["per_arm_timeout_seconds"] == LIMIT_SECONDS
    assert frozen["per_arm_address_space_limit_bytes"] == LIMIT_BYTES
    for filename, expected in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == expected, filename
    spec = frozen["specs"][f"n{n}_L{length}_eval"]
    points = HERE / spec["points_file"]
    fixture = HERE / spec["fixture_file"]
    assert sha(points) == spec["points_sha256"]
    assert sha(fixture) == spec["fixture_sha256"]
    assert not out.exists(), "run directory exists; never overwrite evidence"
    out.mkdir(parents=True)
    host = {
        "platform": platform.platform(), "uname": tuple(platform.uname()),
        "python": sys.version, "cpu_count": os.cpu_count(),
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "rho_binary_sha256": sha(rho), "ic_binary_sha256": sha(ic),
        "load_at_start": os.getloadavg(), "pid": os.getpid(),
    }
    (out / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    expected = [record["published_fixture_scalar"] for record in jsonlines(fixture)]
    k = frozen["batch_k"][str(n)]
    tune = []
    if length == 1:
        tuning = frozen["specs"][f"n{n}_L1_tune"]
        tune_points = HERE / tuning["points_file"]
        tune_fixture = HERE / tuning["fixture_file"]
        assert sha(tune_points) == tuning["points_sha256"]
        assert sha(tune_fixture) == tuning["fixture_sha256"]
        tune_scalar = jsonlines(tune_fixture)[0]["published_fixture_scalar"]
        scores = {}
        for candidate in frozen["single_k_grid"][str(n)]:
            accepted = []
            for rep in range(frozen["tuning_repetitions"]):
                prefix = out / "tune" / f"k{candidate}_r{rep}"
                target_path = prefix.with_suffix(".target.jsonl")
                command = [str(ic), f"construct:{n}:0:{candidate}",
                           str(tune_points), "7", str(target_path)]
                item = child(command, {}, prefix,
                             limit_seconds=frozen["tuning_timeout_seconds"])
                item.update({"k": candidate, "repeat": rep,
                             "valid": item["exit_code"] == 0
                             and not item["timeout"]
                             and valid_tune(prefix, tune_scalar, candidate)})
                tune.append(item)
                (out / "tune_runs.json").write_text(json.dumps(
                    tune, indent=2, sort_keys=True) + "\n")
                if item["valid"]:
                    accepted.append(item["user_seconds"] + item["sys_seconds"])
            if len(accepted) == frozen["tuning_repetitions"]:
                scores[candidate] = sorted(accepted)[len(accepted) // 2]
        assert scores, "no candidate completed all tuning repetitions"
        k = min(scores, key=lambda candidate: (scores[candidate], candidate))
        (out / "tune_scores.json").write_text(json.dumps(
            {"scores_cpu_seconds": scores, "chosen_k": k, "runs": tune},
            indent=2, sort_keys=True) + "\n")
    (out / "chosen_k.json").write_text(json.dumps(
        {"n": n, "L": length, "k": k, "selection": "disjoint_L1_tune_cpu_median"
         if length == 1 else "historical_disjoint_L1024_tune"}, sort_keys=True) + "\n")
    runs = []
    for block in range(frozen["evaluation_repetitions"]):
        for arm in (("ic", "rho") if block % 2 == 0 else ("rho", "ic")):
            prefix = out / f"b{block}_{arm}"
            if arm == "ic":
                env = {
                    "KIC_DUMP_BASE": str(prefix.with_suffix(".base.jsonl")),
                    "KIC_DUMP_RANK": str(prefix.with_suffix(".rank.jsonl")),
                }
                command = [str(ic), f"construct:{n}:0:{k}", str(points),
                           "7", str(prefix.with_suffix(".target.jsonl"))]
            else:
                env = {
                    "KIC_RHO_POINT_INPUT": str(points),
                    "KIC_RHO_BATCH_CORPUS": spec["corpus"],
                    "KIC_RHO_DP_BITS": str(frozen["rho_dp_bits"]),
                }
                command = [str(rho), str(n), "0", "signed_frobenius",
                           str(length), str(spec["seed"])]
            item = child(command, env, prefix)
            item.update({"block": block, "arm": arm})
            runs.append(item)
            (out / "runs.json").write_text(json.dumps(runs, indent=2, sort_keys=True) + "\n")
    (out / "host.json").write_text(json.dumps(
        dict(host, load_at_end=os.getloadavg()), indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(37, 41, 53), required=True)
    parser.add_argument("--L", type=int, choices=(1, 1024), required=True)
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--ic", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.n, args.L, args.rho.resolve(), args.ic.resolve(), args.out.resolve())


if __name__ == "__main__":
    main()
