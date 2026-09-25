#!/usr/bin/env python3
"""Post-outcome diagnosis of the archived first-attempt base mismatch.

This is a read-only analysis of already sealed bytes, not a replacement for
the frozen verifier or a revised admission rule.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tempfile

import verify_archive

HERE = Path(__file__).resolve().parent
FIXTURE = (HERE.parent / "autolab_orbit_extract_20260924" /
           "independent_replay_20260924_codex/base_header.jsonl.gz")


def digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle
    manifest = json.loads((bundle / "archive_manifest.json").read_text())
    archive = bundle / "evidence.tar.gz"
    assert verify_archive.sha(archive) == manifest["archive_sha256"]
    with tempfile.TemporaryDirectory() as temp:
        unpacked = Path(temp)
        assert verify_archive.unpack_checked(archive, unpacked) == manifest["files"]
        panel = unpacked / "panel"
        summary = json.loads((panel / "panel.json").read_text())
        assert verify_archive.sha(panel / "panel.json") == manifest["panel_sha256"]
        assert summary["classification"] == "CENSORED_OR_INVALID_FRESH_BASE"
        assert set(summary["steps"]) == {"training_producer", "base_materialization"}
        assert summary["steps"]["training_producer"]["returncode"] == 0
        assert summary["steps"]["base_materialization"]["returncode"] != 0
        raw = json.loads((panel / "training/producer.stdout.jsonl").read_text())
        fresh = raw["compact_orbit_base_header"]
        certified = json.loads(gzip.decompress(FIXTURE.read_bytes()))
        fresh_points = set(map(tuple, fresh["factor_base_point_coordinates"]))
        certified_points = set(map(tuple, certified["factor_base_point_coordinates"]))
        batch = raw["compact_orbit_batch"]
        report = {
            "classification": summary["classification"],
            "archive_sha256": manifest["archive_sha256"],
            "checkout_head": summary["checkout_head"],
            "github_run_id": summary["github_run_id"],
            "github_run_attempt": summary["github_run_attempt"],
            "training_targets_requested": batch["targets_requested"],
            "training_targets_extracted": batch["targets_extracted"],
            "training_failed_targets": len(batch["failed_target_scalars"]),
            "training_loaded_prior_base": any(raw.get(key) is not None for key in (
                "factor_base_input_path", "factor_base_input_hash", "factor_base_input_blake3")),
            "fresh_cofactor": fresh["cofactor"],
            "certified_cofactor": certified["cofactor"],
            "fresh_scanned_x": fresh["scanned_x"],
            "certified_scanned_x": certified["field_x_values_scanned"],
            "fresh_points": len(fresh_points),
            "certified_points": len(certified_points),
            "coordinate_intersection": len(fresh_points & certified_points),
            "fresh_coordinates_sha256": digest(fresh["factor_base_point_coordinates"]),
            "certified_coordinates_sha256": digest(certified["factor_base_point_coordinates"]),
            "labels_equal": fresh["factor_base_point_labels"] == certified["factor_base_point_labels"],
            "representatives_equal": fresh["factor_base_representatives"] == certified["factor_base_representatives"],
        }
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
