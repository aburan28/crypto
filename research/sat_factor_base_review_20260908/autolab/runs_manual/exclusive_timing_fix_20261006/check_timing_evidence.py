"""Replay the archived timing smoke and a legacy overlapping online row."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from boundary_autolab import exclusive_online_timing_ok, exclusive_precomputation_ms


def close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1e-9, abs_tol=1e-6)


def main() -> None:
    raw = gzip.decompress((HERE / "smoke-n13.jsonl.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == (
        "09caf64810a7d3d40200d611bfaf28c1e0c606e5ded27e017da45c0e096ba417"
    )
    rows = [json.loads(line) for line in raw.splitlines()]
    summaries = [row for row in rows if row.get("kind") == "relation_rank_summary"]
    assert len(summaries) == 2
    precompute, online = summaries
    assert precompute["precomputation_fixture"] is True
    assert precompute["linear_solution_verified"] is True
    assert close(exclusive_precomputation_ms(precompute), precompute["charged_total_ms"])
    assert online["status"] == "SHARED_FACTOR_LOG_ONE_RELATION"
    assert online["linear_solution_verified"] is True
    assert online["all_relations_group_verified"] is True
    assert online["recovered_fixture_scalar"] == online["published_fixture_scalar"]
    assert exclusive_online_timing_ok(online)
    phases = online["target_online_phase_ms"]
    breakdown = online["timing_breakdown_ms"]
    assert close(sum(phases.values()), online["target_online_wall_ms"])
    assert close(phases["target_query"],
                 online["fixture_setup_ms"] + breakdown["target_generation"])
    assert close(phases["target_relation_check"],
                 online["reference_validation_ms"] + breakdown["packed_verification"])
    assert close(phases["target_pdp"],
                 online["collection_ms"] - breakdown["target_generation"]
                 - breakdown["packed_verification"] - online["linear_solve_ms"]
                 - online["solution_validation_ms"])
    batch = next(row for row in rows
                 if row.get("kind") == "retained_support_batch_summary")
    assert close(batch["online_charged_ms"], sum(
        row["fixture_setup_ms"] + row["collection_ms"]
        + row["reference_validation_ms"] for row in summaries
    ))

    legacy_bytes = (HERE / "legacy-n41-online-extract.json").read_bytes()
    assert hashlib.sha256(legacy_bytes).hexdigest() == (
        "687f94b04fe783acec8d1f0a9ca3f6c4eec0ddbb1c5ab90e382c186e33d3a62a"
    )
    legacy = json.loads(legacy_bytes)
    assert legacy["linear_solution_verified"] is True
    assert close(sum(legacy["target_online_phase_ms"].values()),
                 legacy["target_online_wall_ms"])
    assert not exclusive_online_timing_ok(legacy)
    overlap = legacy["linear_solve_ms"] + legacy["solution_validation_ms"]
    assert close(legacy["target_online_wall_ms"] - overlap, 4.613833)
    print(json.dumps({
        "status": "PASS_EXCLUSIVE_NEW_REJECT_OVERLAPPING_LEGACY",
        "new_n13_online_ms": online["target_online_wall_ms"],
        "new_n13_precomputation_ms": precompute["charged_total_ms"],
        "legacy_n41_published_online_ms": legacy["target_online_wall_ms"],
        "legacy_n41_overlap_ms": overlap,
        "legacy_n41_overlap_adjusted_ms": legacy["target_online_wall_ms"] - overlap,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
