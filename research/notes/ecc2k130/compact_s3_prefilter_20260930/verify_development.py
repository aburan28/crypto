#!/usr/bin/env python3
"""Replay same-binary off/filter runs on the prior published development Q."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import subprocess

from verify_panel import Curve, Field, ROOT, check_target, rows, sha, target_identity, verify_rank

OLD = ROOT / "research/notes/ecc2k130/compact_s3_batch_20260929"


def verify(raw: Path, baseline_binary: Path | None, candidate_binary: Path | None) -> dict:
    old_freeze = json.loads((OLD / "FROZEN.json").read_text())
    result = {"status": "PASS", "kind": "s3_prefilter_development_replay",
              "source_sha256": sha(ROOT / "examples/koblitz_orbit_dlp_s3_batch.rs"),
              "baseline_source_sha256": "e998120842faf393910a6b4a9740b3af9880eb626578fa88c575a32e5254f2c7",
              "baseline_binary_sha256": sha(baseline_binary) if baseline_binary else None,
              "candidate_binary_sha256": sha(candidate_binary) if candidate_binary else None,
              "host": platform.platform(),
              "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
              "cells": {}}
    for n, k in ((41, 255), (53, 440)):
        spec = old_freeze["specs"][f"n{n}_L1024_eval"]
        fixture_path = OLD / spec["fixture_file"]
        assert sha(fixture_path) == spec["fixture_sha256"]
        fixtures = rows(fixture_path)
        assert len(fixtures) == 1024
        curve = Curve(Field(n, spec["field_modulus_low_terms"]), 0)
        generator = tuple(spec["generator"])
        cell = {"n": n, "k": k, "fixture_sha256": sha(fixture_path), "arms": {}}
        identities = {}
        for mode in ("original", "off", "filtered"):
            stem = raw / f"n{n}_{mode}"
            base_path = Path(f"{stem}.base.jsonl")
            rank_path = Path(f"{stem}.rank.jsonl")
            target_path = Path(f"{stem}.target.jsonl")
            stdout = Path(f"{stem}.stdout")
            rank_receipt = verify_rank(rank_path, base_path, stdout)
            assert rank_receipt["status"] == "PASS" and rank_receipt["rank"] == k
            base, = rows(base_path)
            rank_rows = rows(rank_path)
            logs = rank_rows[-1]["logs"]
            targets = rows(target_path)
            summary, = rows(stdout)
            assert len(targets) == 1024
            assert summary["rank"] == k and summary["targets_solved"] == 1024
            if mode != "original":
                assert summary["root_prefilter_policy"] == (
                    "blocked_bloom_512_3hash" if mode == "filtered" else "off")
            checked_labels = set()
            for record, fixture in zip(targets, fixtures):
                check_target(record, fixture, base, logs, curve, generator, checked_labels)
                assert record["probes"] == record["s3_counts"][
                    "table_lookups" if mode == "original" else "root_keys_considered"]
            identities[mode] = (base, rank_rows, tuple(target_identity(record) for record in targets))
            cell["arms"][mode] = {
                "rank_verified": k, "targets_verified": len(targets),
                "root_prefilter_bytes": summary.get("root_prefilter_bytes", 0),
                "rank_s3_counts": summary["rank_s3_counts"],
                "target_s3_counts": summary["target_s3_counts"],
                "raw_sha256": {path.name: sha(path) for path in
                    (base_path, rank_path, target_path, stdout)},
            }
        assert identities["original"] == identities["off"] == identities["filtered"], n
        cell["all_original_off_filter_base_rank_witness_and_scalar_identities_equal"] = True
        result["cells"][str(n)] = cell
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--baseline-binary", type=Path)
    parser.add_argument("--candidate-binary", type=Path)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a development receipt"
    receipt = verify(args.raw.resolve(), args.baseline_binary, args.candidate_binary)
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": receipt["status"], "cells": list(receipt["cells"])}, sort_keys=True))


if __name__ == "__main__":
    main()
