#!/usr/bin/env python3
"""Execute the preregistered n37/n41 shared-log panel on a clear host."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from run import GiB, HERE, REPO, SPEC, run_one, sha, write_json

CURVE_CAP_S = 5400


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def replay(mode: str, run: Path, training: Path | None) -> dict:
    report = run / "independent_validation.json"
    cmd = [sys.executable, str(HERE / "verify.py"), "--mode", mode,
           "--run", str(run), "--out", str(report)]
    if mode == "compact":
        assert training is not None
        cmd += ["--training", str(training),
                "--training-report", str(training / "independent_validation.json")]
    else:
        assert training is None or mode == "rho"
    started = time.monotonic_ns()
    with (run / "replay.stdout.txt").open("wb") as stdout, (run / "replay.stderr.txt").open("wb") as stderr:
        completed = subprocess.run(cmd, cwd=REPO, stdout=stdout, stderr=stderr,
                                   env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    result = {"command": cmd, "returncode": completed.returncode,
              "wall_ms": (time.monotonic_ns()-started)/1e6,
              "stdout_sha256": sha((run / "replay.stdout.txt").read_bytes()),
              "stderr_sha256": sha((run / "replay.stderr.txt").read_bytes()),
              "report_sha256": sha(report.read_bytes()) if report.exists() else None}
    write_json(run / "replay_receipt.json", result)
    return result


def expected_classification(mode: str) -> str:
    return "FULL_RANK" if mode == "train" else "COMPLETE"


def run_child(panel: dict, out: Path, mode: str, n: int, block: int | None,
              length: int | None, timeout: int, training: Path | None,
              started: float) -> dict:
    name = "train" if mode == "train" else f"L{length}-b{block}-{mode}"
    record = {"name": name, "mode": mode, "block": block, "length": length,
              "child_timeout_seconds": timeout, "status": "UNRUN", "started_utc": None}
    panel["attempts"].append(record)
    remaining = CURVE_CAP_S - (time.monotonic()-started)
    if remaining < timeout:
        record["status"] = "CENSORED_CURVE_WALL_BUDGET"
        record["remaining_budget_s"] = remaining
        write_json(out / "panel_summary.json", panel)
        return record
    record["started_utc"] = utc_now()
    run_dir = out / name
    receipt = run_one(mode, n, block, length, run_dir, timeout)
    record["run_dir"] = name
    record["producer_receipt_sha256"] = sha((run_dir / "receipt.json").read_bytes())
    record["producer_wall_ms"] = receipt["wall_ms"]
    record["producer_cpu_s"] = receipt["user_cpu_s"] + receipt["system_cpu_s"]
    record["peak_rss_bytes"] = receipt["peak_rss_bytes"]
    if receipt["termination"] is not None or receipt["returncode"] != 0:
        record["status"] = "CENSORED_PRODUCER_FAILURE"
        record["failure_reason"] = receipt["termination"] or f"exit_{receipt['returncode']}"
    elif receipt["peak_rss_bytes"] >= 2*GiB:
        record["status"] = "CENSORED_RSS_CAP_OBSERVED"
    else:
        replay_receipt = replay(mode, run_dir, training)
        record["replay_receipt_sha256"] = sha((run_dir / "replay_receipt.json").read_bytes())
        record["replay_wall_ms"] = replay_receipt["wall_ms"]
        if replay_receipt["returncode"]:
            record["status"] = "CENSORED_REPLAY_FAILURE"
        else:
            report = json.loads((run_dir / "independent_validation.json").read_bytes())
            record["classification"] = report["classification"]
            record["status"] = ("COMPLETE" if report["classification"] == expected_classification(mode)
                                else "CENSORED_INCOMPLETE_LOGS")
    record["finished_utc"] = utc_now()
    write_json(out / "panel_summary.json", panel)
    print(json.dumps(record, sort_keys=True), flush=True)
    return record


def pair(panel: dict, out: Path, n: int, block: int, length: int,
         training: Path, started: float) -> bool:
    compact_first = (block % 2 == 0) == (n == 37)
    order = ("compact", "rho") if compact_first else ("rho", "compact")
    completed = []
    for mode in order:
        record = run_child(panel, out, mode, n, block, length,
                           180 if length == 8 else 600, training, started)
        completed.append(record["status"] == "COMPLETE")
        if record["status"] == "CENSORED_CURVE_WALL_BUDGET":
            break
    panel["stages"].append({"length": length, "block": block, "order": order,
                            "classification": "COMPLETE" if len(completed) == 2 and all(completed)
                            else "CENSORED"})
    write_json(out / "panel_summary.json", panel)
    return len(completed) == 2 and all(completed)


def mark_unrun(panel: dict, out: Path, reason: str, stages: list[tuple[int, int]]) -> None:
    for length, block in stages:
        panel["stages"].append({"length": length, "block": block,
                                "classification": "UNRUN_CENSORED", "reason": reason})
    write_json(out / "panel_summary.json", panel)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, choices=(37, 41), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--host-window-confirmed", action="store_true",
                        help="assert no concurrent #754/#757 timed host measurements")
    args = parser.parse_args()
    assert args.host_window_confirmed, "wait for root's clear-host coordination"
    spec_bytes = SPEC.read_bytes()
    assert (HERE / "input_spec.json").is_file()
    out = args.out.resolve()
    assert not out.is_relative_to(REPO.resolve()), "raw timed outputs must live outside the checkout"
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    panel = {"schema_version": "1.0", "n": args.n,
             "input_spec_sha256": hashlib.sha256(spec_bytes).hexdigest(),
             "checkout_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                      cwd=REPO, text=True).strip(),
             "host": platform.node(), "platform": platform.platform(),
             "started_utc": utc_now(), "curve_wall_cap_seconds": CURVE_CAP_S,
             "attempts": [], "stages": [], "extension_rule": "FROZEN_PROTOCOL"}
    write_json(out / "panel_summary.json", panel)
    training = out / "train"
    first = run_child(panel, out, "train", args.n, None, None, 1800, None, started)
    if first["status"] != "COMPLETE":
        mark_unrun(panel, out, "TRAINING_CENSORED", [(length, block)
                    for length in (8, 32) for block in range(3)])
        panel["classification"] = "TRAINING_CENSORED"
    else:
        l8_complete = True
        for block in range(3):
            if not pair(panel, out, args.n, block, 8, training, started):
                mark_unrun(panel, out, "PRIOR_L8_PAIR_CENSORED",
                           [(8, later) for later in range(block+1, 3)] +
                           [(32, later) for later in range(3)])
                l8_complete = False
                break
        if not l8_complete:
            panel["classification"] = "L8_CENSORED"
        else:
            compact = [next(row for row in panel["attempts"]
                            if row["name"] == f"L8-b{block}-compact") for block in range(3)]
            l8_walls = [row["producer_wall_ms"]/1000 for row in compact]
            all_l8_rss = [row["peak_rss_bytes"] for row in panel["attempts"]
                          if row["length"] == 8]
            can_extend = 4*max(l8_walls) < 600 and all(value < 1.5*GiB for value in all_l8_rss)
            panel["L32_entry"] = {"rule_met": can_extend,
                                  "four_times_max_L8_compact_wall_s": 4*max(l8_walls),
                                  "max_L8_rss_bytes": max(all_l8_rss)}
            if not can_extend:
                mark_unrun(panel, out, "L32_RESOURCE_ENTRY_RULE_FALSE",
                           [(32, block) for block in range(3)])
                panel["classification"] = "L8_COMPLETE_L32_NOT_ENTERED"
            else:
                panel["classification"] = "L8_COMPLETE_L32_CENSORED"
                for block in range(3):
                    if not pair(panel, out, args.n, block, 32, training, started):
                        mark_unrun(panel, out, "PRIOR_L32_PAIR_CENSORED",
                                   [(32, later) for later in range(block+1, 3)])
                        break
                else:
                    panel["classification"] = "L8_AND_L32_COMPLETE"
    panel["finished_utc"] = utc_now()
    panel["observed_panel_wall_s"] = time.monotonic()-started
    write_json(out / "panel_summary.json", panel)
    print(json.dumps({"n": args.n, "classification": panel["classification"],
                      "observed_panel_wall_s": panel["observed_panel_wall_s"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
