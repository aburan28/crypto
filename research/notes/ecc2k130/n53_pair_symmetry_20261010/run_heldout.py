#!/usr/bin/env python3
"""Run the frozen six-seed pair-swap IC and same-point rho comparison once."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = ROOT / "research/notes/ecc2k130/n53_fixed_base_rank_restarts_20261010"
sys.path.insert(0, str(PARENT))
from run_pilot import one_jsonl, run_capped, save, sha  # noqa: E402

ORDER = (
    (532053, "control"), (532053, "symmetry"), (532053, "rho"),
    (532054, "symmetry"), (532054, "rho"), (532054, "control"),
    (532055, "rho"), (532055, "control"), (532055, "symmetry"),
    (532056, "control"), (532056, "rho"), (532056, "symmetry"),
    (532057, "symmetry"), (532057, "control"), (532057, "rho"),
    (532058, "rho"), (532058, "symmetry"), (532058, "control"),
)


def committed_and_pushed(relative: str) -> None:
    committed = subprocess.run(["git", "show", f"HEAD:{relative}"], cwd=ROOT,
                               capture_output=True, check=True)
    assert committed.stdout == (ROOT / relative).read_bytes(), f"not committed: {relative}"


def replay(command: list[str], directory: Path, name: str) -> dict:
    with (directory / f"{name}.stdout").open("wb") as stdout, (
        directory / f"{name}.stderr"
    ).open("wb") as stderr:
        process = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr, timeout=600)
    assert process.returncode == 0, f"{name} exit {process.returncode}"
    receipt = json.loads((directory / f"{name}.json").read_text())
    assert receipt["status"] == "PASS"
    return receipt


def verify_ic(directory: Path, label: str, seed: int, frozen: dict) -> dict:
    base = one_jsonl(directory / "base.jsonl")
    summary = one_jsonl(directory / "summary.jsonl")
    target = one_jsonl(directory / "targets.jsonl")
    q = json.loads((HERE / "inputs/heldout_q.jsonl").read_text())
    pair_symmetry = label == "symmetry"
    trace = [json.loads(line) for line in (directory / "rank.jsonl").read_text().splitlines()]
    attempts = trace[1:-1]
    assert trace[0]["rank_probe_cap"] == trace[-1]["rank_probe_cap"] is None
    assert trace[0]["pair_symmetry"] is pair_symmetry
    assert base["base_hash"] == frozen["base_hash_blake3"]
    assert sha(directory / "base.jsonl") == frozen["base_file_sha256"]
    assert (base["n"], base["a"], base["factor_base_points"], base["orbit_columns"]) == (
        53, 0, 23320, 220
    )
    assert summary["rank"] == trace[-1]["rank"] == 220
    assert summary["rank_probe_cap"] is None
    assert summary["pair_symmetry"] is pair_symmetry
    assert summary["regular_states"] == 2_565_200
    assert summary["index_pair_candidates"] == (1_282_710 if pair_symmetry else 2_565_200)
    assert summary["root_table_entries"] == 2_564_528
    assert summary["rank_attempts"] == len(attempts)
    assert summary["rank_probes_total"] == sum(row["probes"] for row in attempts)
    assert summary["rank_capped_attempts"] == 0
    assert not any(row["capped"] for row in attempts)
    assert summary["rank_root_calls_total"] == sum(row["root_calls"] for row in attempts)
    assert summary["rank_skipped_symmetric_states"] == sum(
        row["skipped_symmetric_states"] for row in attempts
    )
    assert summary["targets"] == summary["targets_solved"] == 1
    assert summary["targets_failed"] == 0
    assert target["target"] == q and target["published_fixture_scalar"] is None
    assert target["group_verified"] is True
    assert abs(sum(summary["cold_phase_ms"].values()) - summary["cold_in_process_ms"]) < 1e-6
    assert abs(target["online_ms"] - target["target_phase_sum_ms"]) < 1e-6
    rank_receipt = replay([
        sys.executable,
        str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929/verify_rank.py"),
        "--trace", str(directory / "rank.jsonl"), "--base", str(directory / "base.jsonl"),
        "--summary", str(directory / "summary.jsonl"),
        "--out", str(directory / "rank_replay.json"),
    ], directory, "rank_replay")
    target_receipt = replay([
        sys.executable, str(PARENT / "verify_public_target.py"),
        "--base", str(directory / "base.jsonl"),
        "--target", str(directory / "targets.jsonl"),
        "--expected-q", str(HERE / "inputs/heldout_q.jsonl"),
        "--out", str(directory / "target_replay.json"),
    ], directory, "target_replay")
    assert rank_receipt["rank"] == 220
    assert target_receipt["target"] == q
    return {
        "rank_seed": seed,
        "rank": summary["rank"],
        "rank_attempts": summary["rank_attempts"],
        "index_pair_candidates": summary["index_pair_candidates"],
        "regular_states": summary["regular_states"],
        "root_table_entries": summary["root_table_entries"],
        "rank_root_calls_total": summary["rank_root_calls_total"],
        "rank_skipped_symmetric_states": summary["rank_skipped_symmetric_states"],
        "cold_phase_ms": summary["cold_phase_ms"],
        "rank_probes_total": summary["rank_probes_total"],
        "cold_in_process_ms": summary["cold_in_process_ms"],
        "target_online_ms": target["online_ms"],
        "recovered_scalar": target_receipt["recovered_scalar"],
        "rank_trace_sha256": sha(directory / "rank.jsonl"),
        "target_sha256": sha(directory / "targets.jsonl"),
        "rank_replay_sha256": sha(directory / "rank_replay.json"),
        "target_replay_sha256": sha(directory / "target_replay.json"),
    }


def verify_rho(directory: Path, seed: int) -> dict:
    receipt = replay([
        sys.executable, str(PARENT / "verify_rho_public.py"),
        "--rho", str(directory / "rho.jsonl"),
        "--public-point", str(HERE / "inputs/heldout_q.jsonl"),
        "--verifier-fixture", str(HERE / "inputs/heldout_fixture_verifier_only.jsonl"),
        "--seed", str(seed),
        "--out", str(directory / "rho_replay.json"),
    ], directory, "rho_replay")
    return {
        "rho_seed": seed,
        "rho_online_ms": receipt["online_ms"],
        "rho_walk_steps": receipt["walk_steps"],
        "recovered_scalar": receipt["recovered_scalar"],
        "rho_sha256": sha(directory / "rho.jsonl"),
        "rho_replay_sha256": sha(directory / "rho_replay.json"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ic-binary", type=Path, required=True)
    parser.add_argument("--rho-binary", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN_PREP.json").read_text())
    heldout = json.loads((HERE / "HELDOUT_FROZEN.json").read_text())
    assert sha(args.ic_binary) == heldout["ic_binary_sha256"] == frozen["ic_binary_sha256"]
    assert sha(args.rho_binary) == heldout["rho_binary_sha256"]
    assert sha(ROOT / "examples/koblitz_orbit_dlp_fast_online.rs") == heldout["ic_source_sha256"]
    assert sha(ROOT / "examples/koblitz_rho_batch_ks_strong_online.rs") == heldout["rho_source_sha256"]
    assert sha(HERE / "inputs/heldout_q.jsonl") == heldout["public_point_sha256"]
    assert sha(HERE / "run_heldout.py") == heldout["runner_sha256"]
    assert sha(PARENT / "verify_rho_public.py") == heldout["rho_verifier_sha256"]
    assert sha(PARENT / "verify_public_target.py") == heldout["ic_verifier_sha256"]
    for path in ["FROZEN_PREP.json", "HELDOUT_FROZEN.json", "HELDOUT_GENERATION.json",
                 "inputs/heldout_q.jsonl", "inputs/heldout_fixture_verifier_only.jsonl",
                 "run_heldout.py", "candidate_control.json", "candidate_symmetry.json"]:
        committed_and_pushed(f"research/notes/ecc2k130/n53_pair_symmetry_20261010/{path}")
    for seed in heldout["rank_seeds"]:
        committed_and_pushed(f"research/notes/ecc2k130/n53_pair_symmetry_20261010/heldout_workload_{seed}.json")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
    upstream = subprocess.check_output(["git", "rev-parse", "@{u}"], cwd=ROOT)
    assert head == upstream, "held-out source and inputs must be pushed before timing"
    q = heldout["public_point"]
    for seed, label in ORDER:
        offset = seed - heldout["rank_seeds"][0]
        rho_seed = heldout["rho_seeds"][offset]
        workload_id = heldout["workload_ids"][str(seed)]
        run_id = (f"{heldout['ic_candidate_ids'][label]}{workload_id}R1"
                  if label != "rho" else f"RHO1N53Ce0{workload_id}R1")
        directory = HERE / "runs/heldout" / str(seed) / label
        environment = os.environ.copy()
        for key in ("KIC_PAIR_SYMMETRY", "KIC_RANK_PROBE_CAP", "KIC_RHO_TARGET_POINT", "KIC_RHO_EXPLICIT_SCALAR",
                    "KIC_RHO_BATCH_CORPUS", "KIC_RHO_FIXTURE_OFFSET"):
            environment.pop(key, None)
        if label == "rho":
            command = [str(args.rho_binary.resolve()), "53", "0", "signed_frobenius", "1",
                       str(rho_seed)]
            environment.update({"KIC_RHO_TARGET_POINT": f"{q[0]},{q[1]}",
                                "KIC_RHO_RUNG": "3", "KIC_RHO_LANES": "32",
                                "KIC_RHO_DP_BITS": "4"})
            stdout_name = "rho.jsonl"
        else:
            command = [str(args.ic_binary.resolve()), "construct:53:0:220",
                       str((HERE / "inputs/heldout_q.jsonl").resolve()), str(seed),
                       str((directory / "targets.jsonl").resolve())]
            environment["KIC_DUMP_BASE"] = str((directory / "base.jsonl").resolve())
            environment["KIC_DUMP_RANK"] = str((directory / "rank.jsonl").resolve())
            environment["KIC_PAIR_SYMMETRY"] = "1" if label == "symmetry" else "0"
            stdout_name = "summary.jsonl"
        if args.dry_run:
            print(json.dumps({"seed": seed, "label": label, "run_id": run_id,
                              "command": command}, sort_keys=True))
            continue
        assert not directory.exists(), f"refusing to replace {directory}"
        directory.mkdir(parents=True)
        intent = {
            "run_id": run_id, "label": label, "rank_seed": seed, "rho_seed": rho_seed,
            "workload_id": workload_id, "command": command,
            "rank_probe_cap": None,
            "pair_symmetry": label == "symmetry" if label != "rho" else None,
            "public_point_sha256": heldout["public_point_sha256"],
            "source_sha256": heldout["rho_source_sha256"] if label == "rho" else heldout["ic_source_sha256"],
            "binary_sha256": heldout["rho_binary_sha256"] if label == "rho" else heldout["ic_binary_sha256"],
            "runner_sha256": sha(HERE / "run_heldout.py"),
            "resource_envelope": heldout["resource_envelope"],
            "environment_flags": {
                "KIC_RANK_PROBE_CAP": environment.get("KIC_RANK_PROBE_CAP"),
                "KIC_PAIR_SYMMETRY": environment.get("KIC_PAIR_SYMMETRY"),
                "KIC_RHO_TARGET_POINT": environment.get("KIC_RHO_TARGET_POINT"),
                "KIC_RHO_RUNG": environment.get("KIC_RHO_RUNG"),
                "KIC_RHO_LANES": environment.get("KIC_RHO_LANES"),
                "KIC_RHO_DP_BITS": environment.get("KIC_RHO_DP_BITS"),
            },
        }
        save(directory / "intent.json", intent)
        process = run_capped(command, environment, directory / stdout_name,
                             directory / "producer.stderr",
                             heldout["resource_envelope"]["wall_cap_seconds"])
        status = {**intent, **process}
        try:
            assert not process["timed_out"], "producer wall cap"
            assert process["exit_code"] == 0, f"producer exit {process['exit_code']}"
            assert process["peak_rss_bytes"] <= heldout["resource_envelope"]["observed_rss_cap_bytes"]
            status.update(verify_rho(directory, rho_seed) if label == "rho"
                          else verify_ic(directory, label, seed, frozen))
            status["status"] = "VERIFIED"
        except BaseException as error:
            status["status"] = "INVALID_OR_INCOMPLETE"
            status["error_type"] = type(error).__name__
            status["error"] = str(error)
            status["traceback"] = traceback.format_exc()
        save(directory / "status.json", status)
        print(json.dumps({"seed": seed, "label": label, "status": status["status"],
                          "online_ms": status.get("rho_online_ms", status.get("target_online_ms")),
                          "cold_ms": status.get("cold_in_process_ms"),
                          "probes": status.get("rank_probes_total")}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
