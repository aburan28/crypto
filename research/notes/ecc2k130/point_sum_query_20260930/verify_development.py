#!/usr/bin/env python3
"""Independent replay and comparison of the preregistered public-Q pilot."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PANEL = ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"
RANK = ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"
sys.path.insert(0, str(RANK))
sys.path.insert(0, str(PANEL))
from verify_panel import check_fixture, check_target  # noqa: E402
from verify_rank import verify as verify_rank  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lines(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def verify_cell(n: int, length: int, frozen: dict) -> dict:
    spec = frozen["specs"][f"n{n}_L{length}_eval"]
    known, curve, generator = check_fixture(spec)
    raw: dict[str, dict] = {}
    for arm in ("control", "point"):
        stem = f"n{n}_{arm}_release"
        paths = {
            name: HERE / "development" / f"{stem}_{name}.{suffix}"
            for name, suffix in (
                ("base", "jsonl"),
                ("rank", "jsonl"),
                ("summary", "json"),
                ("targets", "jsonl"),
            )
        }
        rank = verify_rank(paths["rank"], paths["base"], paths["summary"])
        assert rank["status"] == "PASS" and rank["rank"] == spec.get("chosen_k", rank["rank"])
        base, = lines(paths["base"])
        summary, = lines(paths["summary"])
        targets = lines(paths["targets"])
        assert len(targets) == length
        assert summary["rank"] == summary["orbit_columns"]
        assert summary["targets_solved"] == length and summary["targets_failed"] == 0
        assert summary["query_backend"] == ("s3" if arm == "control" else "point_sum")
        checked_labels: set[tuple[str, int]] = set()
        for target, fixture in zip(targets, known):
            check_target(target, fixture, base, lines(paths["rank"])[-1]["logs"],
                         curve, generator, checked_labels)
        raw[arm] = {"base": base, "summary": summary, "targets": targets,
                    "rank": rank, "hashes": {name: sha(path) for name, path in paths.items()}}

    control = raw["control"]
    point = raw["point"]
    a = control["summary"]
    b = point["summary"]
    for key in ("base_hash", "orbit_columns", "factor_base_points", "regular_states",
                "root_table_entries", "root_table_slots", "rank_attempts", "rank_relations",
                "rank_failures", "rank", "targets_solved", "targets_failed"):
        assert a[key] == b[key], key
    assert control["hashes"]["base"] == point["hashes"]["base"]
    assert a["index_s3_counts"] == b["index_s3_counts"]
    assert a["rank_s3_counts"]["calls"] == b["rank_s3_counts"]["calls"]
    assert a["target_s3_counts"]["calls"] == b["target_s3_counts"]["calls"]
    assert [r["recovered_scalar"] for r in control["targets"]] == [
        r["recovered_scalar"] for r in point["targets"]]
    assert b["point_index_y_bytes"] > 0
    assert b["point_index_group_additions"] == 2 * b["regular_states"]
    return {
        "n": n, "targets": length, "k": a["orbit_columns"],
        "fixture_sha256": sha(PANEL / spec["fixture_file"]),
        "points_sha256": sha(PANEL / spec["points_file"]),
        "control": {"hashes": control["hashes"], "rank": control["rank"],
                    "process_ms": a["timing_ms"]["process_total"],
                    "index_ms": a["timing_ms"]["index_build"],
                    "rank_ms": a["timing_ms"]["rank_stage"],
                    "targets_ms": a["timing_ms"]["targets_total"],
                    "rank_calls": a["rank_s3_counts"]["calls"],
                    "target_calls": a["target_s3_counts"]["calls"]},
        "point": {"hashes": point["hashes"], "rank": point["rank"],
                  "process_ms": b["timing_ms"]["process_total"],
                  "index_ms": b["timing_ms"]["index_build"],
                  "rank_ms": b["timing_ms"]["rank_stage"],
                  "targets_ms": b["timing_ms"]["targets_total"],
                  "rank_calls": b["rank_s3_counts"]["calls"],
                  "target_calls": b["target_s3_counts"]["calls"],
                  "point_multiplications": b["rank_s3_counts"]["point_query_multiplications"]
                                           + b["target_s3_counts"]["point_query_multiplications"],
                  "fallback_calls": b["rank_s3_counts"]["point_query_fallback_calls"]
                                    + b["target_s3_counts"]["point_query_fallback_calls"],
                  "point_index_y_bytes": b["point_index_y_bytes"],
                  "point_index_fallback_states": b["point_index_fallback_states"]},
        "target_witness_differences": sum(
            x["point_indices"] != y["point_indices"]
            for x, y in zip(control["targets"], point["targets"])),
        "target_probe_differences": sum(
            x["probes"] != y["probes"]
            for x, y in zip(control["targets"], point["targets"])),
        "status": "PASS",
    }


def main() -> None:
    frozen = json.loads((PANEL / "FROZEN.json").read_text())
    cells = [verify_cell(n, length, frozen) for n, length in ((37, 1024), (41, 1), (53, 1))]
    receipt = {
        "schema": "compact-point-sum-development-replay-v1",
        "status": "PASS", "scope": "public-Q development; no eligible timing claim",
        "platform": platform.platform(), "machine": platform.machine(),
        "source_sha256": sha(ROOT / "examples/koblitz_orbit_dlp_s3_batch.rs"),
        "verifier_sha256": sha(Path(__file__)),
        "rank_verifier_sha256": sha(RANK / "verify_rank.py"),
        "point_verifier_sha256": sha(PANEL / "verify_panel.py"),
        "independent_arithmetic_sha256": sha(
            ROOT / "research/sat_factor_base_review_20260908/"
            "autolab_orbit_extract_20260924/independent_replay.py"),
        "cells": cells,
    }
    out = HERE / "DEVELOPMENT_RECEIPT.json"
    assert not out.exists(), "refusing to overwrite replay receipt"
    out.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "cells": [(c["n"], c["targets"]) for c in cells]}))


if __name__ == "__main__":
    main()
