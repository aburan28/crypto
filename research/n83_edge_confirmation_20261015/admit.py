#!/usr/bin/env python3
"""Freeze the unresolved finalists from the two sealed source screens."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
TOURNAMENT = REPO / "research/n83_final_tournament_results_20261015/selection.json"
EDGE = REPO / "research/n83_root11_edgek_results_20261015/selection.json"
SOURCE_ADMISSION = REPO / "research/n83_final_tournament_20261015/admission.json"
OUTPUT = HERE / "admission.json"
TARGET_SHA256 = "66c4683e252f1cd321b4b6c7c9d06744c82140222ddce160faa4e7a0e0984b7e"
ANSWER_SHA256 = "e1207d1be9a7689c42209ed50f048325812436d5f7e6f313b5e5e8a04a5730b1"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def admission():
    final, edge, source = read(TOURNAMENT), read(EDGE), read(SOURCE_ADMISSION)
    comparison = final["balanced_comparison"]
    if (final["schema"] != "n83.final-tournament-selection/v1"
            or final["status"] != "BOUNDED_UNKNOWN"
            or final["selected_pipeline"] is not None
            or final["artifact_run_id"] != 38793258972
            or final["artifact_sha256"] !=
            "a0762097a9f0ac62adda606ff9d4d8cbd94e585c1b4e5bb7dbc37f04e4ed80f6"
            or final["analysis_sha256"] !=
            "3a474962a4c8e58d7db416ca90b9746c4f3d05a047647308ab2937d46cfc147f"
            or edge["schema"] != "n83.root11-edgek-selection/v1"
            or edge["analysis_sha256"] !=
            "f50192442355a64fd1c4cfd07147f5b35f8aef4ae77d783f1c7614035992c3eb"
            or edge["no_new_challenger_under_frozen_confirmation_rule"] is not True
            or source["schema"] != "n83.final-tournament-admission/v1"
            or final["admission_sha256"] != digest(SOURCE_ADMISSION)):
        raise ValueError("sealed source comparison or edge screen differs")
    means = comparison["pipeline_mean_ms"]
    leader = comparison["descriptive_lowest_pipeline"]
    if leader != min(means, key=lambda name: (means[name], name)):
        raise ValueError("source comparison leader differs")
    accepted = {leader}
    for name in means:
        if name == leader:
            continue
        upper = comparison["pairwise"][f"{leader}_minus_{name}"]["contiguous_block_bootstrap_99pct_ms"][1]
        if upper >= 0:
            accepted.add(name)
    if not 2 <= len(accepted) <= 6:
        raise ValueError("no open multi-candidate comparison")
    candidates = [{**row, "position": index,
                   "arm_names": [f"{row['pipeline']}_a", f"{row['pipeline']}_b"]}
                  for index, row in enumerate(source["candidate_order"])
                  if row["pipeline"] in accepted]
    if len(candidates) != len(accepted):
        raise ValueError("source candidate identity missing")
    for index, row in enumerate(candidates):
        row["position"] = index
    return {
        "schema": "n83.edge-confirmation-admission/v1",
        "status": "FROZEN_ADMISSION",
        "candidate_order": candidates,
        "mirrored_arm_order": [row["pipeline"] for row in candidates] +
                              [row["pipeline"] for row in reversed(candidates)],
        "target_sha256": TARGET_SHA256,
        "answer_sidecar_sha256": ANSWER_SHA256,
        "source_tournament_selection_sha256": digest(TOURNAMENT),
        "source_edge_selection_sha256": digest(EDGE),
        "source_tournament_admission_sha256": digest(SOURCE_ADMISSION),
        "source_tournament_artifact_id": final["artifact_id"],
        "source_edge_artifact_id": edge["artifact_id"],
        "source_descriptive_leader": leader,
        "edge_challenger": None,
    }


def main():
    blob = json.dumps(admission(), sort_keys=True, indent=2).encode() + b"\n"
    if OUTPUT.exists() and OUTPUT.read_bytes() != blob:
        raise ValueError("frozen admission differs")
    if not OUTPUT.exists():
        OUTPUT.write_bytes(blob)
    print(json.dumps({"admission_sha256": digest(OUTPUT),
                      "candidates": [row["pipeline"] for row in read(OUTPUT)["candidate_order"]]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
