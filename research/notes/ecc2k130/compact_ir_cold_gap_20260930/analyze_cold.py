#!/usr/bin/env python3
"""Recompute the preregistered cold CPU table from sealed raw artifacts."""
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
from run_cold import CONFIG, load_cell, schedule  # noqa: E402

T_CRIT_95 = {5: 2.7764451051977987, 20: 2.093024054408263}
IR_ANALYSIS = ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930/ANALYSIS.json"


def paired(values: list[float]) -> dict:
    assert len(values) in T_CRIT_95 and all(value > 0 for value in values)
    logs = [math.log(value) for value in values]
    midpoint = statistics.mean(logs)
    half = T_CRIT_95[len(values)] * statistics.stdev(logs) / math.sqrt(len(values))
    return {"median": statistics.median(values),
            "interval_95pct": [math.exp(midpoint - half), math.exp(midpoint + half)]}


def read_sealed_cell(archive: Path, cell_id: str,
                     hashes: dict[str, str]) -> tuple[dict, dict, dict[str, bytes]]:
    assert sha(archive) and archive.is_file()
    found = {}
    needed = {f"{cell_id}/cold_run.json", f"{cell_id}/receipt.json"}
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            assert member.isfile() and member.name.startswith(f"{cell_id}/")
            relative = member.name.removeprefix(f"{cell_id}/")
            assert relative and ".." not in Path(relative).parts
            assert relative in hashes and relative not in found
            stream = bundle.extractfile(member)
            assert stream is not None
            data = stream.read()
            assert hashlib.sha256(data).hexdigest() == hashes[relative]
            found[relative] = data
    assert found.keys() == hashes.keys() and needed <= found.keys()
    return (json.loads(found[f"{cell_id}/cold_run.json"]),
            json.loads(found[f"{cell_id}/receipt.json"]), found)


def analyze(archive: Path) -> dict:
    manifest = json.loads((archive / "MANIFEST.json").read_text())
    assert manifest["schema"] == "ecc2k130-compact-ir-cold-hosted-archive-v1"
    assert manifest["config_sha256"] == sha(CONFIG)
    host_path = archive / manifest["second_host_path"]
    assert sha(host_path) == manifest["second_host_sha256"]
    host = json.loads(host_path.read_text())
    assert host["measured_main_head"] == manifest["source_head"]
    ir = json.loads(IR_ANALYSIS.read_text())
    config = json.loads(CONFIG.read_text())
    assert ir["config_sha256"] == config["instruction_config_sha256"]
    assert manifest["source_freeze_sha256"] == config["source_freeze_sha256"]
    result = {"schema": "ecc2k130-compact-ir-cold-analysis-v1",
              "source_head": manifest["source_head"],
              "github_run_url": manifest["github_run_url"],
              "manifest_sha256": sha(archive / "MANIFEST.json"),
              "instruction_analysis_sha256": sha(IR_ANALYSIS),
              "all_hosted_cells_pass": True,
              "all_cells_second_host_replayed": True,
              "all_cells_timing_eligible": True,
              "method_crossover": None,
              "cells": {}}
    binaries = None
    for cell in config["cells"]:
        cell_id = cell["id"]
        entry = manifest["cases"][cell_id]
        assert entry["status"] == "PASS"
        raw = archive / entry["raw_path"]
        assert sha(raw) == entry["raw_sha256"] and raw.stat().st_size == entry["raw_bytes"]
        report, embedded_receipt, members = read_sealed_cell(
            raw, cell_id, entry["member_sha256"])
        assert report["host"]["git_head"] == manifest["source_head"]
        assert report["cell"] == cell and report["mode"] == "measure"
        assert report["status"] == "PASS" and report["config_sha256"] == sha(CONFIG)
        assert [(run["block"], run["arm"]) for run in report["runs"]] == schedule(
            config, cell, "measure")
        assert report["plan"] == [{"block": block, "arm": arm}
                                  for block, arm in schedule(config, cell, "measure")]
        if binaries is None:
            binaries = report["binaries"]
        else:
            assert report["binaries"] == binaries
        assert sha(archive / entry["receipt_path"]) == entry["receipt_sha256"]
        hosted_receipt = json.loads((archive / entry["receipt_path"]).read_text())
        assert equivalent(embedded_receipt, hosted_receipt)
        second = archive / entry["second_host_replay_path"]
        assert sha(second) == entry["second_host_replay_sha256"]
        assert host["receipts_sha256"][cell_id] == sha(second)
        replay_receipt = json.loads(second.read_text())
        assert equivalent(hosted_receipt, replay_receipt)
        assert hosted_receipt["status"] == "PASS"
        assert len(hosted_receipt["checked"]) == 3 * cell["blocks"]
        assert all(check["targets_verified"] == cell["L"]
                   for check in hosted_receipt["checked"])
        measurements = {block: {} for block in range(cell["blocks"])}
        compact_summaries = []
        for run in report["runs"]:
            assert run["status"] == "SMOKE_PASS" and run["exit_code"] == 0
            assert run["stopped_for"] is None and run["cpu_seconds"] > 0
            measurements[run["block"]][run["arm"]] = run
            if run["policy"] == "off":
                stdout = members[f"{cell_id}/{run['subdir']}/off.stdout.jsonl"]
                summary = json.loads(stdout.splitlines()[-1])
                assert summary["kind"] == "compact_orbit_dlp_summary"
                assert summary["n"] == cell["n"] and summary["rank"] == cell["K"]
                assert summary["targets_solved"] == cell["L"]
                compact_summaries.append(summary)
        assert len(compact_summaries) == 2 * cell["blocks"]
        aa, cpu, wall, off_cpu, rho_cpu = [], [], [], [], []
        for block in range(cell["blocks"]):
            data = measurements[block]
            assert set(data) == {"off_a", "rho", "off_b"}
            a, rho, b = data["off_a"], data["rho"], data["off_b"]
            aa.append(b["cpu_seconds"] / a["cpu_seconds"])
            compact = math.sqrt(a["cpu_seconds"] * b["cpu_seconds"])
            off_cpu.append(compact)
            rho_cpu.append(rho["cpu_seconds"])
            cpu.append(compact / rho["cpu_seconds"])
            wall.append(math.sqrt(a["wall_seconds"] * b["wall_seconds"])
                        / rho["wall_seconds"])
        aa_summary, cpu_summary, wall_summary = paired(aa), paired(cpu), paired(wall)
        assert equivalent(aa_summary, hosted_receipt["paired"]["off_b_over_off_a"])
        assert equivalent(cpu_summary, hosted_receipt["paired"]["off_geo_over_rho_cpu"])
        assert equivalent(wall_summary, hosted_receipt["paired"]["off_geo_over_rho_wall"])
        aa_valid = 0.9 <= aa_summary["median"] <= 1.1 and (
            aa_summary["interval_95pct"][0] <= 1 <= aa_summary["interval_95pct"][1])
        assert aa_valid == hosted_receipt["aa_valid"]
        assert hosted_receipt["uncontended"] and hosted_receipt["timing_eligible"]
        ir_row, = (row for row in ir["cells"][cell_id]["rows"] if row["arm"] == "off")
        def median_phase_fraction(phase: str) -> float:
            return statistics.median(summary["timing_ms"][phase] /
                                     summary["timing_ms"]["process_total"]
                                     for summary in compact_summaries)

        phase_diagnostics = {
            "class": "in_process_wall_timers_not_child_CPU_accounting",
            "off_arms": len(compact_summaries),
            "median_rank_fraction": median_phase_fraction("rank_stage"),
            "median_target_fraction": median_phase_fraction("targets_total"),
            "median_index_fraction": median_phase_fraction("index_build"),
            "median_rank_s3_calls": statistics.median(
                summary["rank_s3_counts"]["calls"] for summary in compact_summaries),
            "median_target_s3_calls": statistics.median(
                summary["target_s3_counts"]["calls"] for summary in compact_summaries),
            "median_rank_probes_per_relation": statistics.median(
                summary["rank_probes_mean"] for summary in compact_summaries),
        }
        row = {"status": "PASS", "n": cell["n"], "L": cell["L"],
               "K": cell["K"], "blocks": cell["blocks"],
               "children": len(report["runs"]),
               "verified_target_logs": len(report["runs"]) * cell["L"],
               "full_rank_compact_traces": 2 * cell["blocks"],
               "off_cpu_seconds_median": statistics.median(off_cpu),
               "rho_cpu_seconds_median": statistics.median(rho_cpu),
               "off_geo_over_rho_cpu": cpu_summary,
               "off_b_over_off_a": aa_summary,
               "coarse_off_geo_over_rho_wall": wall_summary,
               "wall_eligible": False,
               "wall_reason": "0.2-second wait4 polling coarsens short reference arms",
               "max_compact_rss_mib": max(run["max_rss_kib_linux"] for run in
                                           report["runs"] if run["policy"] == "off") / 1024,
               "max_rho_rss_mib": max(run["max_rss_kib_linux"] for run in
                                       report["runs"] if run["policy"] == "rho") / 1024,
               "complete_ir_over_rho": ir_row["same_Q_Ir_ratio_to_rho"],
               "phase_diagnostics": phase_diagnostics,
               "timing_eligible": aa_valid and hosted_receipt["uncontended"],
               "method_crossover": None}
        result["cells"][cell_id] = row
    assert set(result["cells"]) == set(manifest["cases"])
    result["binaries_sha256"] = binaries
    result["total_children"] = sum(row["children"] for row in result["cells"].values())
    result["total_verified_target_logs"] = sum(
        row["verified_target_logs"] for row in result["cells"].values())
    result["total_full_rank_compact_traces"] = sum(
        row["full_rank_compact_traces"] for row in result["cells"].values())
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a cold analysis"
    result = analyze(args.archive.resolve())
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "children": result["total_children"],
                      "target_logs": result["total_verified_target_logs"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
