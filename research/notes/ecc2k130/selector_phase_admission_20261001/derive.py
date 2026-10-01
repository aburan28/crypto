#!/usr/bin/env python3
"""Replay v2 phase receipts and price two-index selector counterfactuals."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import tarfile

HERE = Path(__file__).resolve().parent
OUTCOME = HERE.parent / "disjoint_cold_v2_outcome_20261001"
ARCHIVE = OUTCOME / "evidence_run_36803331080"
CELLS = ("n41_L1024", "n53_L1024")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def median_range(values: list[float]) -> dict[str, float]:
    assert len(values) == 5 and all(math.isfinite(x) for x in values)
    return {"median": statistics.median(values), "min": min(values), "max": max(values)}


def member(tf: tarfile.TarFile, checksums: dict[str, str], cell: str,
           relative: str) -> bytes:
    assert relative in checksums, relative
    info = tf.getmember(f"disjoint-cold-v2-{cell}/{relative}")
    assert info.isfile(), relative
    data = tf.extractfile(info).read()
    assert sha(data) == checksums[relative], relative
    return data


def derive_cell(cell: str, manifest: dict, analysis: dict) -> dict:
    entry = manifest["cases"][cell]
    assert entry["status"] == "PASS" and entry["second_host_matches_hosted"]
    assert analysis["cells"][cell]["timing_eligible"]
    path = ARCHIVE / entry["raw_path"]
    data = path.read_bytes()
    assert len(data) == entry["raw_bytes"] and sha(data) == entry["raw_sha256"]
    checksums = entry["member_sha256"]
    with tarfile.open(path, "r:gz") as tf:
        run_bytes = member(tf, checksums, cell, f"{cell}/cold_run.json")
        assert sha(run_bytes) == entry["run_json_sha256"]
        run = json.loads(run_bytes)
        assert run["status"] == "PASS" and run["cell"] == cell
        assert run["host"]["git_head"] == manifest["source_head"]
        lookup = {(int(r["block"]), r["arm"]): r for r in run["runs"]}
        assert len(lookup) == 15 and len(run["runs"]) == 15
        blocks = []
        for b in range(5):
            t = {}
            cpu = {}
            bases = {}
            for arm in ("ic_a", "ic_b", "rho"):
                raw = lookup[(b, arm)]
                assert raw["exit_code"] == 0 and raw["stopped_for"] is None
                cpu[arm] = raw["child_user_cpu_seconds"] + raw["child_system_cpu_seconds"]
                assert cpu[arm] > 0
                if arm.startswith("ic"):
                    name = f"{cell}/b{b:02d}_{arm}.stdout.jsonl"
                    lines = member(tf, checksums, cell, name).splitlines()
                    assert len(lines) == 1
                    obs = json.loads(lines[0])
                    assert obs["rank"] == analysis["cells"][cell]["K"]
                    assert obs["targets_solved"] == 1024 and obs["targets_failed"] == 0
                    assert obs["orbit_columns"] == analysis["cells"][cell]["K"]
                    assert obs["rank_attempts"] == obs["rank_relations"] == obs["rank"]
                    assert obs["rank_failures"] == obs["rank_rows_without_gain"] == 0
                    bases[arm] = obs["base_hash"]
                    t[arm] = obs["timing_ms"]
                    assert all(t[arm][k] > 0 for k in
                               ("process_total", "index_build", "rank_stage", "targets_total"))
            assert bases["ic_a"] == bases["ic_b"]
            ic_cpu = math.sqrt(cpu["ic_a"] * cpu["ic_b"])
            ratio = ic_cpu / cpu["rho"]
            total = math.sqrt(t["ic_a"]["process_total"] * t["ic_b"]["process_total"])
            index = math.sqrt(t["ic_a"]["index_build"] * t["ic_b"]["index_build"])
            variable = math.sqrt(
                (t["ic_a"]["rank_stage"] + t["ic_a"]["targets_total"])
                * (t["ic_b"]["rank_stage"] + t["ic_b"]["targets_total"]))
            index_share = index / total
            variable_share = variable / total
            assert 0 < index_share < 1 and 0 < variable_share < 1
            assert index_share + variable_share < 1
            # Scenario assumption: internal elapsed-phase shares proxy charged CPU shares.
            # The second full index is additional; scan, ranking and allocation are free.
            blocks.append({
                "block": b, "paired_ic_over_rho_cpu": round(ratio, 12),
                "index_elapsed_share": round(index_share, 12),
                "rank_plus_targets_elapsed_share": round(variable_share, 12),
                "zero_score_extra_required_variable_saving_fraction":
                    round((1 - 1 / ratio) / variable_share, 12),
                "second_full_index_required_variable_saving_fraction":
                    round((1 + index_share - 1 / ratio) / variable_share, 12),
            })
    observed = statistics.median(x["paired_ic_over_rho_cpu"] for x in blocks)
    expected = analysis["cells"][cell]["ic_over_rho_cpu"]["median"]
    assert math.isclose(observed, expected, rel_tol=1e-12), (cell, observed, expected)
    keys = [k for k in blocks[0] if k != "block"]
    return {"cell": cell, "K": analysis["cells"][cell]["K"],
            "L": 1024, "raw_sha256": entry["raw_sha256"],
            "source_head": manifest["source_head"],
            "blocks": blocks,
            "five_block_summary": {k: median_range([b[k] for b in blocks]) for k in keys}}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest_path = ARCHIVE / "MANIFEST.json"
    analysis_path = OUTCOME / "ANALYSIS.json"
    manifest = json.loads(manifest_path.read_bytes())
    analysis = json.loads(analysis_path.read_bytes())
    assert manifest["source_head"] == analysis["source_head"]
    result = {"schema": "ecc2k130-selector-phase-admission-v1",
              "classification": "counterfactual planning diagnostic, not measured selector cost",
              "assumption": "within-compact elapsed-phase fractions proxy charged CPU fractions; one extra full index costs the observed index fraction; other selector costs are zero",
              "manifest_sha256": sha(manifest_path.read_bytes()),
              "analysis_sha256": sha(analysis_path.read_bytes()),
              "source_sha256": sha(Path(__file__).read_bytes()),
              "cells": {cell: derive_cell(cell, manifest, analysis) for cell in CELLS}}
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({cell: {k: v["median"] for k, v in row["five_block_summary"].items()}
                      for cell, row in result["cells"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
