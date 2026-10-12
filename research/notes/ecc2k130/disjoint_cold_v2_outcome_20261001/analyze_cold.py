#!/usr/bin/env python3
"""Recompute v2 CPU decisions only from sealed, independently replayed bytes."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent  # noqa: E402
from run_panel import sha  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001"))
from prepare import CELLS, cell_name  # noqa: E402
from run_cold import ARM_ORDER, schedule  # noqa: E402
from verify_cold import observed_isolation_ok, paired  # noqa: E402
from archive_cold import (INPUT_FREEZE, INPUT_RECEIPT, RUN_URL, SOURCE_HEAD,
                          VERIFIER)  # noqa: E402


def sealed_members(raw: Path, prefix: str, expected: dict[str, str],
                   selected: set[str]) -> dict[str, dict]:
    assert raw.is_file()
    found = set()
    documents = {}
    with tarfile.open(raw, "r:gz") as bundle:
        for member in bundle:
            assert member.isfile() and member.name.startswith(prefix + "/")
            relative = member.name.removeprefix(prefix + "/")
            assert relative and ".." not in Path(relative).parts
            assert relative in expected and relative not in found
            stream = bundle.extractfile(member)
            assert stream is not None
            digest = hashlib.sha256()
            chunks = [] if relative in selected else None
            while chunk := stream.read(1 << 20):
                digest.update(chunk)
                if chunks is not None:
                    chunks.append(chunk)
            assert digest.hexdigest() == expected[relative], relative
            found.add(relative)
            if chunks is not None:
                documents[relative] = json.loads(b"".join(chunks))
    assert found == set(expected)
    return documents


def analyze(archive: Path) -> dict:
    manifest_path = archive / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema"] == "ecc2k130-disjoint-cold-v2-hosted-archive-v1"
    assert manifest["source_head"] == SOURCE_HEAD
    assert manifest["github_run_url"] == RUN_URL
    assert manifest["input_freeze_sha256"] == sha(INPUT_FREEZE)
    assert manifest["input_receipt_sha256"] == sha(INPUT_RECEIPT)
    assert manifest["verifier_sha256"] == sha(VERIFIER)
    host_path = archive / manifest["second_host_path"]
    assert sha(host_path) == manifest["second_host_sha256"]
    host = json.loads(host_path.read_text())
    assert host["run_url"] == RUN_URL and host["source_head"] == SOURCE_HEAD
    assert host["verifier_sha256"] == sha(VERIFIER)
    validation = manifest["validation"]
    if validation is not None:
        raw = archive / validation["raw_path"]
        assert raw.stat().st_size == validation["raw_bytes"]
        assert sha(raw) == validation["raw_sha256"]
        sealed_members(raw, "disjoint-cold-v2-validation",
                       validation["member_sha256"], set())

    result = {"schema": "ecc2k130-disjoint-cold-v2-analysis-v1",
              "source_head": SOURCE_HEAD, "github_run_url": RUN_URL,
              "manifest_sha256": sha(manifest_path),
              "input_freeze_sha256": sha(INPUT_FREEZE),
              "unit": "wait4 complete-process child user+system CPU seconds",
              "method_crossover": None, "common_group_operation_S": None,
              "n83_or_n131_transfer": None, "cells": {}}
    expected_cells = {cell_name(n, length) for n, length, *_ in CELLS}
    assert set(manifest["cases"]) == expected_cells
    for n, length, k, prefilter, blocks in CELLS:
        cell = cell_name(n, length)
        entry = manifest["cases"][cell]
        if entry["status"] == "MISSING_ARTIFACT":
            result["cells"][cell] = {"status": "MISSING_ARTIFACT",
                                     "timing_eligible": False,
                                     "decision": "CENSORED",
                                     "reason": entry["reason"]}
            continue
        raw = archive / entry["raw_path"]
        assert raw.stat().st_size == entry["raw_bytes"]
        assert sha(raw) == entry["raw_sha256"]
        sealed = sealed_members(raw, f"disjoint-cold-v2-{cell}",
                                entry["member_sha256"],
                                {f"{cell}/cold_run.json", f"{cell}/receipt.json"})
        report = sealed.get(f"{cell}/cold_run.json")
        embedded = sealed.get(f"{cell}/receipt.json")
        if entry["run_json_sha256"] is not None:
            assert report is not None
            assert report["schema"] == "ecc2k130-disjoint-cold-v2-run-v1"
            assert report["host"]["git_head"] == SOURCE_HEAD
            assert report["cell"] == cell and report["mode"] == "measure"
            assert report["frozen_sha256"] == sha(INPUT_FREEZE)
            assert (report["spec"]["K"], report["spec"]["blocks"],
                    report["spec"]["prefilter"]) == (k, blocks, prefilter)
        if entry["receipt_path"] is not None:
            hosted_path = archive / entry["receipt_path"]
            assert sha(hosted_path) == entry["receipt_sha256"]
            hosted = json.loads(hosted_path.read_text())
            assert equivalent(embedded, hosted)
            assert hosted["status"] == entry["status"]
        else:
            assert embedded is None
            hosted = None
        replay_path = entry["second_host_replay_path"]
        if replay_path is not None:
            second_path = archive / replay_path
            assert sha(second_path) == entry["second_host_replay_sha256"]
            assert host["receipts_sha256"][cell] == sha(second_path)
            second = json.loads(second_path.read_text())
        else:
            second = None
        if entry["status"] != "PASS":
            result["cells"][cell] = {"status": entry["status"],
                                     "n": n, "L": length, "K": k,
                                     "blocks": blocks,
                                     "completed_children": len(report["runs"]) if report else 0,
                                     "timing_eligible": False,
                                     "decision": "CENSORED",
                                     "hosted_failure": hosted,
                                     "second_host_status": second["status"] if second else None,
                                     "raw_sha256": entry["raw_sha256"]}
            continue
        assert report is not None and report["status"] == "PASS"
        assert hosted is not None and len(hosted["checks"]) == 3 * blocks
        assert len(report["runs"]) == 3 * blocks
        assert report["plan"] == [{"block": b, "arm": arm}
                                  for b, arm in schedule(blocks, "measure")]
        independently_verified = (second is not None and
                                  equivalent(hosted, second) and
                                  entry["second_host_matches_hosted"] is True)
        by_block = {block: {} for block in range(blocks)}
        for item in report["runs"]:
            assert item["exit_code"] == 0 and item["stopped_for"] is None
            assert item["points_sha256"] == report["spec"]["block_specs"][item["block"]]["points_sha256"]
            cpu = item["child_user_cpu_seconds"] + item["child_system_cpu_seconds"]
            assert cpu > 0
            by_block[item["block"]][item["arm"]] = (cpu, item)
        aa, ratios, ic_cpu, rho_cpu = [], [], [], []
        for block in range(blocks):
            arms = by_block[block]
            assert set(arms) == set(ARM_ORDER)
            a, rho, b = (arms[arm][0] for arm in ARM_ORDER)
            aa.append(b / a)
            pair_ic = math.sqrt(a * b)
            ic_cpu.append(pair_ic)
            rho_cpu.append(rho)
            ratios.append(pair_ic / rho)
        aa_stats, ratio_stats = paired(aa), paired(ratios)
        assert equivalent(aa_stats, hosted["paired"]["ic_b_over_a"])
        assert equivalent(ratio_stats, hosted["paired"]["ic_over_rho"])
        aa_valid = (0.9 <= aa_stats["median"] <= 1.1 and
                    aa_stats["interval_95pct"][0] <= 1 <= aa_stats["interval_95pct"][1])
        assert aa_valid == hosted["aa_valid"]
        isolation = hosted["isolation"]
        uncontended = observed_isolation_ok(
            isolation, f"disjoint-cold-v2-{cell}", report["host"]["reserved_cpu"])
        assert uncontended == hosted["uncontended"]
        assert hosted["host_timing_candidate"] == (aa_valid and uncontended)
        eligible = aa_valid and uncontended and independently_verified
        if not eligible:
            decision = "CENSORED"
        elif ratio_stats["interval_95pct"][1] < 1:
            decision = "CROSSOVER_CANDIDATE_REQUIRES_SECOND_HOST_TIMING"
        elif ratio_stats["interval_95pct"][0] > 1:
            decision = "QUANTITATIVE_NO_GO_AT_THIS_CELL"
        else:
            decision = "UNRESOLVED_INTERVAL"
        result["cells"][cell] = {
            "status": "PASS", "n": n, "L": length, "K": k,
            "prefilter": prefilter, "blocks": blocks,
            "children": 3 * blocks,
            "verified_target_logs": 3 * blocks * length if independently_verified else 0,
            "full_rank_compact_traces": 2 * blocks if independently_verified else 0,
            "ic_cpu_seconds_median": statistics.median(ic_cpu),
            "rho_cpu_seconds_median": statistics.median(rho_cpu),
            "ic_over_rho_cpu": ratio_stats,
            "ic_b_over_a": aa_stats,
            "uncontended": uncontended, "aa_valid": aa_valid,
            "independent_correctness_replay": independently_verified,
            "timing_eligible": eligible, "decision": decision,
            "required_complete_cpu_reduction_to_parity": (
                1 - 1 / ratio_stats["median"] if decision ==
                "QUANTITATIVE_NO_GO_AT_THIS_CELL" else None),
            "max_ic_rss_mib": max(item["child_max_rss_kib_linux"]
                                   for arms in by_block.values()
                                   for arm, (_, item) in arms.items()
                                   if arm != "rho") / 1024,
            "max_rho_rss_mib": max(arms["rho"][1]["child_max_rss_kib_linux"]
                                     for arms in by_block.values()) / 1024,
            "raw_sha256": entry["raw_sha256"]}
    result["all_cells_timing_eligible"] = all(
        row["timing_eligible"] for row in result["cells"].values())
    result["all_cells_fixed_policy_no_go"] = (
        result["all_cells_timing_eligible"] and all(
            row["decision"] == "QUANTITATIVE_NO_GO_AT_THIS_CELL"
            for row in result["cells"].values()))
    result["total_verified_target_logs"] = sum(
        row.get("verified_target_logs", 0) for row in result["cells"].values())
    result["total_full_rank_compact_traces"] = sum(
        row.get("full_rank_compact_traces", 0) for row in result["cells"].values())
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a cold analysis"
    result = analyze(args.archive.resolve())
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS",
                      "all_cells_timing_eligible": result["all_cells_timing_eligible"],
                      "all_cells_fixed_policy_no_go": result["all_cells_fixed_policy_no_go"],
                      "verified_target_logs": result["total_verified_target_logs"]},
                     sort_keys=True))
