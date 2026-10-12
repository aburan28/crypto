#!/usr/bin/env python3
"""Independent read-only diagnosis of the sealed first-attempt training stage.

This runs after the frozen base-identity gate failed. It does not change that
gate, complete the panel, or turn diagnostic time into a measured IC cost.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import importlib.util
import json
from pathlib import Path
import tempfile

import verify_archive

HERE = Path(__file__).resolve().parent
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
RANK_SOURCE = ORBIT / "cold_batch_rank.py"
GROUP_SOURCE = ORBIT / "independent_replay_20260924_codex/replay.py"
ARCHIVE_SHA256 = "2ff3d4587dd9a576de393abacabfc390c2f838997b5dcbcc9b1e2e4554c5c23d"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def replay(bundle: Path) -> dict:
    manifest = json.loads((bundle / "archive_manifest.json").read_text())
    archive = bundle / "evidence.tar.gz"
    assert manifest["archive_sha256"] == ARCHIVE_SHA256
    assert verify_archive.sha(archive) == ARCHIVE_SHA256
    with tempfile.TemporaryDirectory() as temp:
        unpacked = Path(temp)
        assert verify_archive.unpack_checked(archive, unpacked) == manifest["files"]
        panel_dir = unpacked / "panel"
        panel_file = panel_dir / "panel.json"
        assert verify_archive.sha(panel_file) == manifest["panel_sha256"]
        panel = json.loads(panel_file.read_text())
        assert panel["classification"] == "CENSORED_OR_INVALID_FRESH_BASE"
        assert set(panel["steps"]) == {"training_producer", "base_materialization"}
        assert panel["steps"]["training_producer"]["returncode"] == 0
        assert panel["steps"]["base_materialization"]["returncode"] != 0
        assert verify_archive.sha(RANK_SOURCE) == panel["source_sha256"]["training_schedule_reference"]
        assert verify_archive.sha(GROUP_SOURCE) == panel["source_sha256"]["independent_group_reference"]

        rank = load(RANK_SOURCE, "post_outcome_rank")
        verifier = rank.load_verifier()
        raw_file = panel_dir / "training/producer.stdout.jsonl"
        raw = json.loads(raw_file.read_text())
        assert all(raw.get(key) is None for key in (
            "factor_base_input_path", "factor_base_input_hash", "factor_base_input_blake3"))
        header = raw["compact_orbit_base_header"]
        curve = verifier.Curve(header)
        assert curve.on_curve(curve.generator)
        assert curve.scalar(curve.generator, curve.order) is None
        by_point, representatives, _ = rank.verify_orbit_labels(curve, header)
        assert len(by_point) == 23320 and len(representatives) == 220
        by_x = defaultdict(list)
        for point in by_point:
            by_x[point[0]].append(point)
        assert len(by_x) == 11660
        assert all(len(points) == 2 and
                   by_point[points[0]][0] == by_point[points[1]][0]
                   for points in by_x.values())

        batch = raw["compact_orbit_batch"]
        relations = batch["relations"]
        queries = batch["query_observations"]
        assert batch["targets_requested"] == batch["targets_extracted"] == 512
        assert len(relations) == len(queries) == 512
        assert not batch["failed_target_scalars"]
        scalars = rank.target_schedule(512, curve.order)
        archived_scalars = [int(line) for line in
                            (panel_dir / "training/target_scalars.txt").read_text().splitlines()]
        assert scalars == archived_scalars

        matrix = rank.Echelon(220, curve.order)
        column_occurrences = Counter()
        gains = []
        for index, (relation, query, scalar) in enumerate(
                zip(relations, queries, scalars), start=1):
            assert relation["scalar"] == query["scalar"] == scalar
            assert query["hit"]
            assert 0 <= query["partner_roots"] <= 2 * query["s3_calls"]
            assert 0 <= query["group_lift_attempts"] <= query["indexed_partner_hits"] <= query["partner_roots"]
            codes = relation["x_codes"]
            assert len(codes) == 4
            for code in codes:
                assert code in by_x
                column_occurrences[by_point[by_x[code][0]][0]] += 1
            target = curve.scalar(curve.generator, scalar)
            assert target is not None
            lift = verifier.check_witness(
                curve, by_x, target, codes, relation["pinned_intermediates"])
            row = {}
            for point in map(tuple, lift):
                column, coefficient = by_point[point]
                row[column] = (row.get(column, 0) + coefficient) % curve.order
            if matrix.insert(row, scalar):
                gains.append(index)

        assert sum(column_occurrences.values()) == 2048
        missing_columns = sorted(set(range(220)) - set(column_occurrences))
        missing_pivots = sorted(set(range(220)) - set(matrix.pivots))
        assert missing_columns == missing_pivots == [102]
        assert len(matrix.pivots) == len(gains) == 219
        assert matrix.solution() is None
        return {
            "classification": "POST_OUTCOME_DIAGNOSTIC_UNCHARGED",
            "archive_sha256": ARCHIVE_SHA256,
            "diagnostic_source_sha256": verify_archive.sha(Path(__file__)),
            "checkout_head": panel["checkout_head"],
            "github_run_id": panel["github_run_id"],
            "github_run_attempt": panel["github_run_attempt"],
            "training_raw_sha256": verify_archive.sha(raw_file),
            "independent_rank_source_sha256": verify_archive.sha(RANK_SOURCE),
            "independent_group_source_sha256": verify_archive.sha(GROUP_SOURCE),
            "orbit_labels_replayed": len(by_point),
            "representatives_replayed": len(representatives),
            "group_witnesses_replayed": len(relations),
            "relation_x_slots": sum(column_occurrences.values()),
            "columns_visited": len(column_occurrences),
            "zero_occurrence_columns": missing_columns,
            "least_positive_column_occurrences": min(column_occurrences.values()),
            "most_column_occurrences": max(column_occurrences.values()),
            "rank": len(matrix.pivots),
            "columns": 220,
            "missing_pivots": missing_pivots,
            "last_rank_gain_relation": gains[-1],
            "first_full_rank_relation": None,
            "full_rank_solution": None,
            "scope": "Read-only analysis of already sealed training bytes; no rho, point IC or cost comparison",
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = replay(args.bundle)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out is not None:
        args.out.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
