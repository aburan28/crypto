#!/usr/bin/env python3
"""Freeze a selected fresh-confirmation candidate for five same-point pairs."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SOURCE = REPO / "research/n83_edge_confirmation_results_20261015/selection.json"
ADMISSION = REPO / "research/n83_edge_confirmation_20261015/admission.json"
TOURNAMENT_ADMISSION = REPO / "research/n83_final_tournament_20261015/admission.json"
TOURNAMENT_SELECTION = REPO / "research/n83_final_tournament_results_20261015/selection.json"
OUTPUT = HERE / "selection.json"
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected():
    source = json.loads(SOURCE.read_text())
    admission = json.loads(ADMISSION.read_text())
    tournament = json.loads(TOURNAMENT_SELECTION.read_text())
    if (source["schema"] != "n83.edge-confirmation-selection/v1"
            or source["status"] != "SELECTED"
            or source["admission_sha256"] != digest(ADMISSION)
            or source["artifact_id"] <= 0
            or source["artifact_run_id"] <= 0
            or source["source_commit"] !=
               "a369212178ffd28ec30e0e20eaa4729ec7b81e8d"
            or admission["schema"] != "n83.edge-confirmation-admission/v1"
            or tournament["schema"] != "n83.final-tournament-selection/v1"
            or tournament["status"] != "BOUNDED_UNKNOWN"
            or tournament["artifact_id"] != admission["source_tournament_artifact_id"]
            or tournament["admission_sha256"] != digest(TOURNAMENT_ADMISSION)):
        raise ValueError("source confirmation or tournament seal differs")
    rows = [row for row in admission["candidate_order"]
            if row["pipeline"] == source["selected_pipeline"]]
    if len(rows) != 1 or source["selected_candidate_id"] != rows[0]["candidate_id"]:
        raise ValueError("source winner differs from admitted candidate")
    row = rows[0]
    return {
        "schema": "n83.final-tier1-selection/v2",
        "pipeline": row["pipeline"],
        "candidate_id": row["candidate_id"],
        "candidate_manifest_sha256": row["candidate_manifest_sha256"],
        "k": row["k"], "variant": row["variant"],
        "bloom_profile": row["bloom_profile"],
        "source_base_hash_blake3": row["base_hash"],
        "root_table_entries": row["root_table_entries"],
        "source_sha256": SOURCE_SHA256,
        "confirmation_artifact_run_id": source["artifact_run_id"],
        "confirmation_artifact_id": source["artifact_id"],
        "confirmation_artifact_sha256": source["artifact_sha256"],
        "confirmation_analysis_sha256": source["analysis_sha256"],
        "confirmation_selection_sha256": digest(SOURCE),
        "tournament_artifact_id": tournament["artifact_id"],
        "tournament_artifact_sha256": tournament["artifact_sha256"],
        "tournament_analysis_sha256": tournament["analysis_sha256"],
        "admission_sha256": digest(ADMISSION),
        "tournament_admission_sha256": digest(TOURNAMENT_ADMISSION),
        "rss_cap_bytes": 29 * 1024**3 // 2,
        "ic_child_cap_seconds": 1500,
        "rho_child_cap_seconds": 120,
    }


def main():
    blob = json.dumps(selected(), sort_keys=True, indent=2).encode() + b"\n"
    if OUTPUT.exists() and OUTPUT.read_bytes() != blob:
        raise ValueError("frozen Tier 1 selection differs")
    if not OUTPUT.exists():
        OUTPUT.write_bytes(blob)
    print(json.dumps({"selection_sha256": digest(OUTPUT),
                      "candidate_id": json.loads(blob)["candidate_id"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
