#!/usr/bin/env python3
"""Audit the final-binary n13 cap=1 rank-restart correctness control."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
CONTROL = HERE / "controls/n13_cap1_seed12"


def one(path: Path) -> dict:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    assert len(rows) == 1, path
    return rows[0]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    summary = one(CONTROL / "summary.jsonl")
    target = one(CONTROL / "target.jsonl")
    base = one(CONTROL / "base.jsonl")
    trace = [json.loads(line) for line in (CONTROL / "rank.jsonl").read_text().splitlines()]
    attempts = trace[1:-1]
    assert (base["n"], base["a"], base["orbit_columns"]) == (13, 0, 2)
    assert summary["rank_probe_cap"] == trace[0]["rank_probe_cap"] == 1
    assert summary["rank"] == trace[-1]["rank"] == 2
    assert summary["rank_attempts"] == len(attempts) == 5
    assert summary["rank_capped_attempts"] == sum(row["capped"] for row in attempts) == 3
    assert summary["rank_probes_total"] == sum(row["probes"] for row in attempts) == 5
    assert all(row["probes"] <= 1 for row in attempts)
    assert target["group_verified"] is True
    assert target["recovered_scalar"] == target["published_fixture_scalar"] == 7
    assert abs(target["online_ms"] - target["target_phase_sum_ms"]) < 1e-6
    assert abs(sum(summary["cold_phase_ms"].values()) - summary["cold_in_process_ms"]) < 1e-6
    assert json.loads((CONTROL / "rank_replay.json").read_text())["status"] == "PASS"
    assert json.loads((CONTROL / "target_replay.json").read_text())["status"] == "PASS"
    assert not (CONTROL / "stderr.txt").read_text()
    receipt = {
        "status": "PASS",
        "schema": "n13-rank-restart-control-v1",
        "source_sha256": frozen["source_sha256"],
        "binary_sha256": frozen["binary_sha256"],
        "rank_attempts": len(attempts),
        "rank_capped_attempts": summary["rank_capped_attempts"],
        "rank_probes_total": summary["rank_probes_total"],
        "rank": summary["rank"],
        "recovered_scalar": target["recovered_scalar"],
        "cold_in_process_ms": summary["cold_in_process_ms"],
        "target_online_ms": target["online_ms"],
        "base_sha256": sha(CONTROL / "base.jsonl"),
        "rank_trace_sha256": sha(CONTROL / "rank.jsonl"),
        "target_sha256": sha(CONTROL / "target.jsonl"),
    }
    output = CONTROL / "final_control_receipt.json"
    encoded = json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    if output.exists():
        assert output.read_text() == encoded, "control receipt changed"
    else:
        output.write_text(encoded)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
