#!/usr/bin/env python3
"""Seal a replayed same-instance selector before five fresh points are run."""

import argparse
import importlib.util
import json
import math
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = "3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03"
TARGETS_SHA256 = "4976679ac6c9e3858cd26312ec2e50b8eb30971041ea38286c9582dd9217a1c1"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


dynamic = module("n83_tier1_selector_candidate",
                 ROOT / "research/n83_root11_dynamic_20261016/materialize.py")
panel = dynamic.panel
runner = module("n83_tier1_selector_runner", HERE / "selector_run.py")


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"sealed record differs: {path}")
    path.write_text(blob)


def select(results, receipt_path):
    analysis_path = results / "analysis.json"
    analysis = panel.read(analysis_path)
    source = panel.read(results / "source.json")
    host = panel.read(results / "host.json")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                     cwd=ROOT, text=True).strip()
    if (analysis["schema"] != "n83.tier1-selector-analysis/v1"
            or analysis["status"] != "SCREEN_SELECTED"
            or analysis["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or analysis["execution_order_status"] != "PASS"
            or analysis["rho"]["status"] != "PASS_INDEPENDENT_REPLAY"
            or set(analysis["online_aa_controls"]) != {"k3072", "k2543"}
            or any(control["status"] != "PASS_10_PERCENT" for control in
                   analysis["online_aa_controls"].values())
            or host["git_sha"] != commit
            or source["source_sha256"] != SOURCE_SHA256
            or source["runner_sha256"] != panel.file_sha(
                HERE / "selector_run.py")
            or source["arm_order"] != [name for name, _, _, _ in runner.ARMS]
            or panel.file_sha(ROOT / "examples/koblitz_orbit_dlp_fast.rs")
            != SOURCE_SHA256):
        raise ValueError("selector replay, A/A or executed source differs")
    by_arm = {row["arm"]: row for row in analysis["arms"]}
    controls = analysis["online_aa_controls"]
    noise = 0.0
    for prefix in ("k3072", "k2543"):
        control = controls[prefix]
        times = [by_arm[f"{prefix}_{suffix}"]["one_target_online_ms"]
                 for suffix in "ab"]
        difference = abs(times[0] - times[1])
        if (control["online_ms"] != times
                or control["base_hash"] != by_arm[f"{prefix}_a"]["base_hash"]
                or control["base_hash"] != by_arm[f"{prefix}_b"]["base_hash"]
                or difference / (sum(times) / 2) > 0.10
                or not math.isclose(
                    control["pair_absolute_differences_ms"][0], difference,
                    rel_tol=0, abs_tol=1e-9)):
            raise ValueError(f"{prefix}: A/A control differs from replayed arms")
        noise = max(noise, difference)
    t3072 = controls["k3072"]["online_ms"]
    t2543 = controls["k2543"]["online_ms"]
    if min(t2543) - max(t3072) > noise:
        winner = "k3072root11"
        separation = min(t2543) - max(t3072)
    elif min(t3072) - max(t2543) > noise:
        winner = "k2543root11"
        separation = min(t3072) - max(t2543)
    else:
        raise ValueError("no strict online separation beyond observed A/A noise")
    if (analysis["selected_pipeline"] != winner
            or not math.isclose(analysis["winner_separation_ms"], separation,
                                rel_tol=0, abs_tol=1e-9)
            or not math.isclose(analysis["max_adjacent_aa_difference_ms"], noise,
                                rel_tol=0, abs_tol=1e-9)):
        raise ValueError("recorded winner or noise bound differs")
    k = {"k3072root11": 3072, "k2543root11": 2543}.get(winner)
    if k is None:
        raise ValueError("no admitted finalist")
    field, curve = panel.curve_record()
    records = []
    for suffix in "ab":
        name = f"k{k}_{suffix}"
        item = by_arm[name]
        folder = results / "runs" / name
        summary = panel.read(folder / "stdout.jsonl")
        if (item["status"] != "PASS_INDEPENDENT_REPLAY"
                or item["resources"]["status"] != "PASS"
                or summary["base_hash"] != item["base_hash"]
                or item["actual_factor_base_points"] != 166 * k):
            raise ValueError(f"{name}: selector arm differs")
        panel.SOURCE_COMMIT = commit
        panel.SOURCE_SHA256 = SOURCE_SHA256
        records.append(dynamic.candidate(
            k, "root11_dynamic", summary, folder / "base_dump.jsonl",
            field, curve))
    if records[0] != records[1]:
        raise ValueError("A/B exact candidate manifests differ")
    candidate = records[0]
    candidate_id = candidate["candidate_id"]
    candidate_path = HERE / f"candidate_{candidate_id}.json"
    write_once(candidate_path, candidate)
    selection = {
        "schema": "n83.tier1-online-selection/v1",
        "candidate_id": candidate_id,
        "candidate_manifest": str(candidate_path.relative_to(ROOT)),
        "candidate_manifest_sha256": panel.file_sha(candidate_path),
        "k": k, "layout": "root11_dynamic",
        "source_sha256": SOURCE_SHA256,
        "targets_sha256": TARGETS_SHA256,
        "selection_route": "same_instance_four_arm_aa_screen",
        "screen_status": analysis["status"],
        "screen_analysis_sha256": panel.file_sha(analysis_path),
        "screen_source_sha256": panel.file_sha(results / "source.json"),
        "screen_host_sha256": panel.file_sha(results / "host.json"),
        "screen_commit": commit,
        "screen_winner_separation_ms": analysis["winner_separation_ms"],
        "screen_max_aa_difference_ms":
            analysis["max_adjacent_aa_difference_ms"],
    }
    selection_path = HERE / "selection.json"
    write_once(selection_path, selection)
    receipt = {
        "schema": "n83.tier1-selector-seal/v1",
        "selection_sha256": panel.file_sha(selection_path),
        "candidate_manifest_sha256": panel.file_sha(candidate_path),
        "candidate_id": candidate_id,
        "screen_analysis_sha256": panel.file_sha(analysis_path),
        "screen_source_sha256": panel.file_sha(results / "source.json"),
        "screen_host_sha256": panel.file_sha(results / "host.json"),
        "source_commit": commit,
    }
    write_once(receipt_path, receipt)
    print(json.dumps(receipt, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    select(args.results.resolve(), args.receipt.resolve())


if __name__ == "__main__":
    main()
