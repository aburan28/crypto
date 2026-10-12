#!/usr/bin/env python3
"""Run the preregistered nine-cell n53 fixed-base restart pilot once."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
EXPECTED_BASE_HASH = "7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973"
ORDER = (
    (530053, "control"), (530053, "cap200000"), (530053, "cap400000"),
    (530054, "cap200000"), (530054, "cap400000"), (530054, "control"),
    (530055, "cap400000"), (530055, "control"), (530055, "cap200000"),
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one_jsonl(path: Path) -> dict:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    assert len(rows) == 1, f"expected one JSONL row: {path}"
    return rows[0]


def save(path: Path, value: object) -> None:
    assert not path.exists(), f"refusing to replace {path}"
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def run_capped(
    command: list[str], environment: dict[str, str],
    stdout_path: Path, stderr_path: Path, wall_cap_seconds: int,
) -> dict:
    """Wait for one child and retain its own wait4 peak-RSS receipt."""
    started = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, env=environment, stdout=stdout, stderr=stderr)
        timed_out = False
        signal_sent_at = None
        while True:
            pid, raw_status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid == process.pid:
                exit_code = os.waitstatus_to_exitcode(raw_status)
                process.returncode = exit_code
                peak_rss = usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
                return {
                    "exit_code": exit_code,
                    "timed_out": timed_out,
                    "external_wall_ms": (time.monotonic() - started) * 1000,
                    "peak_rss_bytes": peak_rss,
                }
            now = time.monotonic()
            if not timed_out and now - started >= wall_cap_seconds:
                timed_out = True
                signal_sent_at = now
                process.terminate()
            elif timed_out and signal_sent_at is not None and now - signal_sent_at >= 2:
                process.kill()
                signal_sent_at = None
            time.sleep(0.02)


def verify_outputs(directory: Path, label: str, cap: int | None, seed: int, frozen: dict) -> dict:
    summary = one_jsonl(directory / "summary.jsonl")
    target = one_jsonl(directory / "targets.jsonl")
    base = one_jsonl(directory / "base.jsonl")
    trace = [json.loads(line) for line in (directory / "rank.jsonl").read_text().splitlines()]
    assert len(trace) >= 3
    header, solution = trace[0], trace[-1]
    attempts = trace[1:-1]
    assert base["base_hash"] == summary["base_hash"] == header["base_hash"] == EXPECTED_BASE_HASH
    assert sha(directory / "base.jsonl") == frozen["base_file_sha256"]
    assert (base["n"], base["a"], base["factor_base_points"], base["orbit_columns"]) == (
        53, 0, 23_320, 220
    )
    assert summary["rank"] == solution["rank"] == 220
    assert summary["targets"] == summary["targets_solved"] == 1
    assert summary["targets_failed"] == 0 and target["group_verified"] is True
    assert target["target"] == json.loads((HERE / "inputs/development_q.jsonl").read_text())
    assert target["published_fixture_scalar"] is None
    assert summary["rank_probe_cap"] == header["rank_probe_cap"] == solution["rank_probe_cap"] == cap
    assert len(attempts) == summary["rank_attempts"] == solution["attempts"]
    assert sum(int(row["found"]) for row in attempts) == summary["rank_relations"] == solution["relations"]
    assert sum(not row["found"] for row in attempts) == summary["rank_failures"] == solution["failures"]
    assert sum(row["probes"] for row in attempts) == summary["rank_probes_total"] == solution["probes_total"]
    capped = [row for row in attempts if row["capped"]]
    assert len(capped) == summary["rank_capped_attempts"] == solution["capped_attempts"]
    assert sum(row["probes"] for row in capped) == summary["rank_capped_probes"] == solution["capped_probes"]
    for index, row in enumerate(attempts):
        assert row["attempt_index"] == index
        if row["capped"]:
            assert cap is not None and row["probes"] == cap and not row["found"]
            assert row["rank_after"] == row["rank_before"]
        if cap is not None:
            assert row["probes"] <= cap
    assert target["online_start_event"] == "target_query_begin"
    assert target["online_stop_event"] == "recovery_check_end"
    assert abs(target["online_ms"] - target["target_phase_sum_ms"]) < 1e-6
    assert abs(sum(summary["cold_phase_ms"].values()) - summary["cold_in_process_ms"]) < 1e-6

    rank_command = [
        sys.executable,
        str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929/verify_rank.py"),
        "--trace", str(directory / "rank.jsonl"),
        "--base", str(directory / "base.jsonl"),
        "--summary", str(directory / "summary.jsonl"),
        "--out", str(directory / "rank_replay.json"),
    ]
    target_command = [
        sys.executable, str(HERE / "verify_public_target.py"),
        "--base", str(directory / "base.jsonl"),
        "--target", str(directory / "targets.jsonl"),
        "--expected-q", str(HERE / "inputs/development_q.jsonl"),
        "--out", str(directory / "target_replay.json"),
    ]
    for name, command in (("rank", rank_command), ("target", target_command)):
        with (directory / f"{name}_replay.stdout").open("wb") as stdout, (
            directory / f"{name}_replay.stderr"
        ).open("wb") as stderr:
            replay = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=600)
        assert replay.returncode == 0, f"{name} replay exit {replay.returncode}"
        assert json.loads((directory / f"{name}_replay.json").read_text())["status"] == "PASS"
    return {
        "candidate_label": label,
        "rank_seed": seed,
        "rank": summary["rank"],
        "rank_attempts": summary["rank_attempts"],
        "rank_relations": summary["rank_relations"],
        "rank_failures": summary["rank_failures"],
        "rank_capped_attempts": summary["rank_capped_attempts"],
        "rank_capped_probes": summary["rank_capped_probes"],
        "rank_probes_total": summary["rank_probes_total"],
        "cold_in_process_ms": summary["cold_in_process_ms"],
        "target_online_ms": target["online_ms"],
        "recovered_scalar": target["recovered_scalar"],
        "base_sha256": sha(directory / "base.jsonl"),
        "rank_trace_sha256": sha(directory / "rank.jsonl"),
        "summary_sha256": sha(directory / "summary.jsonl"),
        "target_sha256": sha(directory / "targets.jsonl"),
        "rank_replay_sha256": sha(directory / "rank_replay.json"),
        "target_replay_sha256": sha(directory / "target_replay.json"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--one", nargs=2, metavar=("SEED", "LABEL"))
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(ROOT / "examples/koblitz_orbit_dlp_fast_online.rs") == frozen["source_sha256"]
    assert sha(args.binary) == frozen["binary_sha256"]
    assert sha(HERE / "inputs/base_n53_k220.jsonl") == frozen["base_file_sha256"]
    assert sha(HERE / "inputs/development_q.jsonl") == frozen["development_q_sha256"]
    assert sha(HERE / "PROTOCOL.md") == frozen["protocol_sha256"]
    envelope = frozen["resource_envelope"]
    selected = ORDER if args.one is None else ((int(args.one[0]), args.one[1]),)
    for seed, label in selected:
        assert (seed, label) in ORDER
        cap = None if label == "control" else int(label.removeprefix("cap"))
        candidate_id = frozen["candidate_ids"][label]
        workload_id = frozen["workload_ids"][str(seed)]
        run_id = f"{candidate_id}{workload_id}R1"
        directory = HERE / "runs/pilot" / str(seed) / label
        command = [
            str(args.binary.resolve()), "construct:53:0:220",
            str((HERE / "inputs/development_q.jsonl").resolve()),
            str(seed), str((directory / "targets.jsonl").resolve()),
        ]
        if args.dry_run:
            print(json.dumps({"seed": seed, "label": label, "run_id": run_id, "command": command}))
            continue
        assert not directory.exists(), f"refusing to replace run {directory}"
        directory.mkdir(parents=True)
        environment = os.environ.copy()
        environment["KIC_DUMP_BASE"] = str((directory / "base.jsonl").resolve())
        environment["KIC_DUMP_RANK"] = str((directory / "rank.jsonl").resolve())
        if cap is None:
            environment.pop("KIC_RANK_PROBE_CAP", None)
        else:
            environment["KIC_RANK_PROBE_CAP"] = str(cap)
        intent = {
            "run_id": run_id,
            "candidate_id": candidate_id,
            "workload_id": workload_id,
            "rank_seed": seed,
            "rank_probe_cap": cap,
            "command": command,
            "environment_flags": {
                "KIC_RANK_PROBE_CAP": cap,
                "KIC_DUMP_BASE": True,
                "KIC_DUMP_RANK": True,
            },
            "resource_envelope": envelope,
            "source_sha256": frozen["source_sha256"],
            "binary_sha256": frozen["binary_sha256"],
            "base_input_sha256": frozen["base_file_sha256"],
            "target_input_sha256": frozen["development_q_sha256"],
        }
        save(directory / "intent.json", intent)
        process = run_capped(
            command, environment, directory / "summary.jsonl", directory / "producer.stderr",
            envelope["wall_cap_seconds"],
        )
        status = {**intent, **process}
        try:
            assert not process["timed_out"], "producer wall cap"
            assert process["exit_code"] == 0, f"producer exit {process['exit_code']}"
            assert process["peak_rss_bytes"] <= envelope["observed_rss_cap_bytes"], "RSS cap"
            status.update(verify_outputs(directory, label, cap, seed, frozen))
            status["status"] = "VERIFIED"
        except BaseException as error:
            status["status"] = "INVALID_OR_INCOMPLETE"
            status["error_type"] = type(error).__name__
            status["error"] = str(error)
            status["traceback"] = traceback.format_exc()
        save(directory / "status.json", status)
        print(json.dumps({"seed": seed, "label": label, "status": status["status"],
                          "cold_ms": status.get("cold_in_process_ms"),
                          "probes": status.get("rank_probes_total")}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
