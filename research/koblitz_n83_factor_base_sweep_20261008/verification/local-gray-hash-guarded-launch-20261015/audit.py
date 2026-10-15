#!/usr/bin/env python3
"""Independently replay the bounded N83 quiet-host no-launch receipt."""

from __future__ import annotations

from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
AUTH = HERE.parent / "local-gray-hash-paired-20261015" / "authorization.json"

if not __debug__:
    raise RuntimeError("the receipt audit requires Python assertions to be enabled")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_bytes(commit: str, name: str) -> bytes:
    return subprocess.run(["git", "show", f"{commit}:{name}"], cwd=ROOT,
                          capture_output=True, check=True).stdout


def audit(receipts: Path = HERE, check_live_output: bool = True) -> dict:
    config_path = receipts / "config.json"
    samples_path = receipts / "samples.jsonl"
    result_path = receipts / "result.json"
    config = json.loads(config_path.read_bytes())
    result = json.loads(result_path.read_bytes())
    assert config["schema"] == "n83.local-paired-quiet-watch-config/v1"
    assert result["schema"] == "n83.local-paired-quiet-watch-result/v1"
    assert config["max_wait_seconds"] == 900
    assert config["quiet_samples_required"] == 3
    assert config["sample_interval_seconds"] == 45
    assert config["watcher_sha256"] == digest((HERE / "watch.py").read_bytes())
    assert config["authorization_sha256"] == digest(AUTH.read_bytes())
    auth = json.loads(AUTH.read_bytes())
    assert str(Path(auth["output_dir"]).resolve()) == config["authorized_output"]
    commit = config["source_commit"]
    assert re.fullmatch(r"[0-9a-f]{40}", commit)
    for name, expected in config["source_sha256"].items():
        if name == "src:git-tree":
            actual = subprocess.run(["git", "rev-parse", f"{commit}:src"],
                                    cwd=ROOT, capture_output=True, text=True,
                                    check=True).stdout.strip()
        else:
            actual = digest(git_bytes(commit, name))
        assert actual == expected, f"frozen source mismatch: {name}"

    rows = [json.loads(line) for line in samples_path.read_text().splitlines()]
    assert rows, "missing host samples"
    times = []
    quiet_streak = 0
    max_quiet_streak = 0
    loads = []
    idle_values = []
    for index, row in enumerate(rows):
        when = datetime.fromisoformat(row["observed_at_utc"])
        assert when.tzinfo is not None and (not times or when > times[-1])
        times.append(when)
        for key in ("top", "load", "cores", "docker"):
            assert row[key]["exit_code"] == 0, f"failed {key} probe in sample {index}"
        idle_match = re.search(r"CPU usage: .*?([0-9]+(?:\.[0-9]+)?)% idle",
                               row["top"]["stdout"])
        assert idle_match, f"unparseable top sample {index}"
        idle = float(idle_match.group(1))
        load_match = re.search(r"\{\s*([0-9]+(?:\.[0-9]+)?)", row["load"]["stdout"])
        assert load_match, f"unparseable load sample {index}"
        load = float(load_match.group(1))
        cores = int(row["cores"]["stdout"].strip())
        others = [name for name in row["docker"]["stdout"].splitlines()
                  if name and name != "conductor-db"]
        assert row["idle_percent"] == idle
        assert row["load_1m"] == load
        assert row["logical_cpus"] == cores
        assert row["other_containers"] == others
        quiet = (idle >= 30.0 and load <= cores and not others
                 and row["system_free_bytes"] >= 10 * (1 << 30))
        assert row["status"] == ("QUIET" if quiet else "BUSY")
        quiet_streak = quiet_streak + 1 if quiet else 0
        max_quiet_streak = max(max_quiet_streak, quiet_streak)
        loads.append(load)
        idle_values.append(idle)
    assert max_quiet_streak < config["quiet_samples_required"]
    assert 900 <= result["elapsed_watch_seconds_L0"] <= 930
    assert 0 < (times[-1] - times[0]).total_seconds() <= result["elapsed_watch_seconds_L0"]
    assert result == {
        "schema": "n83.local-paired-quiet-watch-result/v1",
        "status": "WAIT_EXPIRED_NO_SOLVER_LAUNCH",
        "launched": False,
        "launch_error": None,
        "supervisor_exit_code": None,
        "admission_sha256": None,
        "output_exists": False,
        "elapsed_watch_seconds_L0": result["elapsed_watch_seconds_L0"],
        "conservative_search_charge_seconds": 0,
    }
    assert not (receipts / "admission.json").exists()
    assert not (receipts / "supervisor.stdout").exists()
    assert not (receipts / "supervisor.stderr").exists()
    if check_live_output:
        assert not Path(config["authorized_output"]).exists()
    return {
        "schema": "n83.local-paired-quiet-watch-audit/v1",
        "status": "PASS_BOUNDED_NO_LAUNCH_ONLY",
        "source_commit": commit,
        "sample_count": len(rows),
        "quiet_sample_count": sum(row["status"] == "QUIET" for row in rows),
        "max_quiet_streak": max_quiet_streak,
        "idle_percent_min_max": [min(idle_values), max(idle_values)],
        "load_1m_min_max": [min(loads), max(loads)],
        "first_sample_utc": times[0].isoformat(),
        "last_sample_utc": times[-1].isoformat(),
        "elapsed_watch_seconds_L0": result["elapsed_watch_seconds_L0"],
        "conservative_search_charge_seconds": 0,
        "receipt_sha256": {p.name: digest(p.read_bytes()) for p in
                           (config_path, samples_path, result_path)},
    }


if __name__ == "__main__":
    receipt = audit()
    with (HERE / "audit.json").open("x", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(receipt, sort_keys=True))
