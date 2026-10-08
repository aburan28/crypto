#!/usr/bin/env python3
"""Recover the n=85 pairs after the original driver was killed.

R1 is PRODUCERS_COMPLETE (untouched).  R2's rho arm completed and its
rows are retained; only its IC arm is rerun here.  R3 runs both arms.
Same receipts as run_n85_pairs.py.
"""
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_n85_pairs as base  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
TARGETS = HERE / "frozen/target_points.jsonl"


def main():
    target = json.loads(TARGETS.read_text().splitlines()[0])
    rho_env = os.environ.copy()
    for key in ("KIC_RHO_FIXED_TARGET_SCALAR", "KIC_RHO_BATCH_CORPUS", "KIC_RHO_FIXTURE_OFFSET"):
        rho_env.pop(key, None)
    rho_env["KIC_RHO_PUBLIC_TARGET_POINT"] = json.dumps(target, separators=(",", ":"))
    rho_env["KIC_RHO_WALK_SEED"] = "202610072"
    ic_env = os.environ.copy()
    for key in ("KIC_RHO_FIXED_TARGET_SCALAR", "KIC_RHO_PUBLIC_TARGET_POINT", "KIC_RHO_BATCH_CORPUS"):
        ic_env.pop(key, None)
    ic_env["KIC_RANK_THREADS"] = str(base.RANK_THREADS)

    # R2: rerun only the IC arm (rho rows retained).
    r2_dir = RUNS / "N85A0K600W12We202610072R32R2"
    run_record = json.loads((r2_dir / "run.json").read_text())
    assert run_record["status"] == "RUNNING_IC", run_record["status"]
    ic_argv = [base.IC_BINARY, base.BASE, TARGETS, "32", r2_dir / "ic.jsonl"]
    ic_status = base.run_arm("ic", ic_argv, ic_env, r2_dir, run_record)
    run_record["launch_finished_at_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    run_record["status"] = (
        "PRODUCERS_COMPLETE" if run_record.get("rho_return_code") == 0 and ic_status == 0
        else "PRODUCER_FAILURE"
    )
    base.write_json(r2_dir / "run.json", run_record)
    print(json.dumps({"run": "R2-recovered", "status": run_record["status"]}))

    # R3: full pair.
    seeds = {"rho_walk_seed": 202610073, "ic_rank_seed": 33}
    run_id = f"N85A0K600W{base.RANK_THREADS}We202610073R33R3"
    r3_dir = RUNS / run_id
    assert not r3_dir.exists(), f"{r3_dir} already exists"
    r3_dir.mkdir(parents=True)
    fixture = json.loads((HERE / "frozen/fixture.json").read_text())
    run_record = {
        "active_arm": None,
        "candidate_id": "IC1N85A0Ckb1fb102000PDP4rootRCguidedLAgaussTDdirectISO0parallel",
        "curve": {
            "n": 85,
            "a": 0,
            "b": 1,
            "field_modulus": "x^85 + x^8 + x^2 + x + 1",
            "subgroup_order": fixture["subgroup_order"],
            "cofactor": fixture["cofactor"],
        },
        "ic_binary_sha256": base.sha256(base.IC_BINARY),
        "ic_rank_seed": 33,
        "ic_rank_threads": base.RANK_THREADS,
        "ic_return_code": None,
        "known_answer_sent_to_ic": False,
        "known_answer_sent_to_rho": False,
        "launch_started_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "memory_limit_mechanism": "No finite kernel cap; per-process peak RSS recorded with wait4(2)",
        "paired_run_order": ["rho", "ic"],
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "resource_cap_bytes_per_arm": None,
        "rho_binary_sha256": base.sha256(base.RHO_BINARY),
        "rho_return_code": None,
        "rho_walk_seed": 202610073,
        "run_id": run_id,
        "sidecar_validation_after_run": True,
        "status": "FROZEN_BEFORE_RUN",
        "target": [str(v) for v in target],
        "target_count": 1,
        "validation_scalar_sidecar_path": "frozen/fixture.json",
    }
    base.write_json(r3_dir / "run.json", run_record)
    rho_env["KIC_RHO_WALK_SEED"] = "202610073"
    rho_argv = [base.RHO_BINARY, "85", "0", "signed_frobenius", "1", "packed"]
    ic_argv = [base.IC_BINARY, base.BASE, TARGETS, "33", r3_dir / "ic.jsonl"]
    rho_status = base.run_arm("rho", rho_argv, rho_env, r3_dir, run_record)
    ic_status = base.run_arm("ic", ic_argv, ic_env, r3_dir, run_record)
    run_record["launch_finished_at_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    run_record["status"] = (
        "PRODUCERS_COMPLETE" if rho_status == 0 and ic_status == 0 else "PRODUCER_FAILURE"
    )
    base.write_json(r3_dir / "run.json", run_record)
    print(json.dumps({"run": "R3", "status": run_record["status"]}))


if __name__ == "__main__":
    main()
