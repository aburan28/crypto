#!/usr/bin/env python3
"""Post-failure diagnostic of the censored n37 batch; never makes timing eligible."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field, verify as verify_rank  # noqa: E402
from verify_archive import verify as verify_archive  # noqa: E402

CELL = "n37_L1024"


def rows(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def matching_pinned_lifts(curve: Curve, points: list[tuple[int, int]],
                          by_x: dict[int, list[int]], record: dict) -> list[list[int]]:
    options = [by_x[code] for code in record["x_codes"]]
    assert all(0 < len(group) <= 4 for group in options)
    matches = []
    for indices in itertools.product(*options):
        chosen = [points[index] for index in indices]
        left = curve.add(chosen[0], chosen[1])
        right = curve.add(chosen[2], chosen[3])
        if left is None or right is None:
            continue
        if curve.add(left, right) != tuple(record["published_q"]):
            continue
        if set(record["pinned_intermediates"]) == {left[0], right[0]}:
            matches.append(list(indices))
    return matches


def diagnose(evidence: Path) -> dict:
    assert verify_archive(evidence)["status"] == "PASS"
    cold_dir = ROOT / "research/notes/ecc2k130/point_sum_cold_20260930"
    frozen = json.loads((cold_dir / "FROZEN.json").read_text())
    spec = frozen["specs"][CELL]
    labels = rows(cold_dir / spec["fixture_file"])
    curve = Curve(Field(37, spec["field_modulus_low_terms"]), 0)
    generator = tuple(spec["generator"])
    archive_path = evidence / "artifacts" / f"point-sum-cold-{CELL}.tar.gz"
    with tempfile.TemporaryDirectory(prefix="point-sum-n37-diagnostic-") as temp:
        root = Path(temp)
        with tarfile.open(archive_path, "r:gz") as archive:
            for member in archive:
                assert member.isfile() and not member.issym() and not member.islnk()
                name = Path(member.name)
                assert not name.is_absolute() and ".." not in name.parts
                target = root / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.extractfile(member).read())
        run_dir = root / CELL
        report = json.loads((run_dir / "cold_run.json").read_text())
        assert report["status"] == "PASS" and len(report["runs"]) == 20
        assert all(item["exit_code"] == 0 and item["stopped_for"] is None
                   for item in report["runs"])
        checked = {"rank_runs": 0, "compact_target_logs": 0, "rho_target_logs": 0}
        mismatches = []
        for item in report["runs"]:
            block, arm = item["block"], item["arm"]
            prefix = f"b{block}_{arm}"
            if arm == "rho":
                data = rows(run_dir / f"{prefix}.stdout.jsonl")
                assert len(data) == spec["L"] + 1 and data[-1]["all_verified"]
                for index, (record, label) in enumerate(zip(data[:-1], labels)):
                    assert record["fixture_index"] == index
                    assert record["published_q"] == label["published_q"]
                    assert record["recovered_fixture_scalar"] == label["published_fixture_scalar"]
                    assert curve.mul(record["recovered_fixture_scalar"], generator) == tuple(record["published_q"])
                    checked["rho_target_logs"] += 1
                continue
            base_path = run_dir / f"{prefix}.base.jsonl"
            rank_path = run_dir / f"{prefix}.rank.jsonl"
            summary_path = run_dir / f"{prefix}.stdout.jsonl"
            rank = verify_rank(rank_path, base_path, summary_path)
            assert rank["status"] == "PASS" and rank["rank"] == spec["K"]
            checked["rank_runs"] += 1
            base, = rows(base_path)
            logs = rows(rank_path)[-1]["logs"]
            points = [tuple(point) for point in base["factor_base_point_coordinates"]]
            by_x = {}
            for point_index, point in enumerate(points):
                by_x.setdefault(point[0], []).append(point_index)
            targets = rows(run_dir / f"{prefix}.target.jsonl")
            assert len(targets) == spec["L"]
            for index, (record, label) in enumerate(zip(targets, labels)):
                assert record["fixture_index"] == index
                assert record["published_q"] == label["published_q"]
                scalar = record["recovered_scalar"]
                assert scalar == label["published_fixture_scalar"]
                assert curve.mul(scalar, generator) == tuple(record["published_q"])
                indices = record["point_indices"]
                chosen = [points[point_index] for point_index in indices]
                assert [point[0] for point in chosen] == record["x_codes"]
                left = curve.add(chosen[0], chosen[1])
                right = curve.add(chosen[2], chosen[3])
                assert left is not None and right is not None
                assert curve.add(left, right) == tuple(record["published_q"])
                recovered = sum(base["factor_base_point_labels"][point_index][1]
                                * logs[base["factor_base_point_labels"][point_index][0]]
                                for point_index in indices) % base["subgroup_order"]
                assert recovered == scalar
                checked["compact_target_logs"] += 1
                if set(record["pinned_intermediates"]) != {left[0], right[0]}:
                    alternatives = matching_pinned_lifts(curve, points, by_x, record)
                    assert alternatives
                    mismatches.append({"block": block, "arm": arm, "fixture_index": index,
                                       "selected_indices": indices,
                                       "pinned_compatible_indices": alternatives,
                                       "x_codes": record["x_codes"]})
    assert checked == {"rank_runs": 15, "compact_target_logs": 15 * 1024,
                       "rho_target_logs": 5 * 1024}
    assert len(mismatches) == 15
    assert {item["fixture_index"] for item in mismatches} == {763}
    assert all(item["x_codes"][0] == item["x_codes"][2] for item in mismatches)
    return {"status": "DIAGNOSTIC_PASS_TIMING_INELIGIBLE", "cell": CELL,
            "checked": checked, "mismatch_count": len(mismatches),
            "distinct_mismatched_q": 1, "mismatches": mismatches,
            "frozen_verifier_status": "FAIL", "timing_eligible": False}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a diagnostic receipt"
    result = diagnose(args.evidence.resolve())
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "mismatches"},
                     sort_keys=True))


if __name__ == "__main__":
    main()
