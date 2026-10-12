#!/usr/bin/env python3
"""Replay the frozen n37 batch failure under a future x-only witness contract.

This is a correctness regression fixture, never a timing admission for the
2026-09-30 run.  It reads only hashed bytes from that run's committed archive.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field, verify as verify_rank  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"))
from verify_panel import check_target as old_check_target  # noqa: E402
from check_witness import check_target  # noqa: E402

OLD = ROOT / "research/notes/ecc2k130/point_sum_cold_20260930"
EVIDENCE = ROOT / "research/notes/ecc2k130/point_sum_cold_outcome_20260930/evidence/run_36764654520"
CELL = "n37_L1024"
ARTIFACT = f"point-sum-cold-{CELL}"
TARGETS = (0, 763)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def lines(data: bytes) -> list[dict]:
    return [json.loads(line) for line in data.decode().splitlines() if line.strip()]


def rejected(call) -> bool:
    try:
        call()
    except AssertionError:
        return True
    return False


def replay(evidence: Path = EVIDENCE) -> dict:
    manifest = json.loads((evidence / "MANIFEST.json").read_text())
    assert manifest["schema"] == "ecc2k130-point-sum-cold-outcome-archive-v1"
    assert manifest["run_id"] == 36764654520
    metadata = manifest["artifacts"][ARTIFACT]
    archive_path = evidence / "artifacts" / f"{ARTIFACT}.tar.gz"
    archive_bytes = archive_path.read_bytes()
    assert len(archive_bytes) == metadata["archive_bytes"]
    assert sha(archive_bytes) == metadata["archive_sha256"]

    frozen = json.loads((OLD / "FROZEN.json").read_text())
    spec = frozen["specs"][CELL]
    fixture_bytes = (OLD / spec["fixture_file"]).read_bytes()
    assert sha(fixture_bytes) == spec["fixture_sha256"]
    fixtures = lines(fixture_bytes)
    assert len(fixtures) == spec["L"]
    curve = Curve(Field(spec["n"], spec["field_modulus_low_terms"]), spec["a"])
    generator = tuple(spec["generator"])
    for index in TARGETS:
        fixture = fixtures[index]
        assert fixture["fixture_index"] == index
        assert curve.mul(fixture["published_fixture_scalar"], generator) == tuple(fixture["published_q"])

    needed = {f"{CELL}/receipt.json"}
    for block in range(spec["blocks"]):
        for arm in ("control_a", "point_sum", "control_b"):
            prefix = f"{CELL}/b{block}_{arm}"
            needed.update(f"{prefix}.{suffix}.jsonl" for suffix in
                          ("base", "rank", "target", "stdout"))
    checked_bytes = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        assert needed <= set(archive.getnames())
        for name in sorted(needed):
            member = archive.getmember(name)
            assert member.isfile() and not member.issym() and not member.islnk()
            data = archive.extractfile(member).read()
            expected = metadata["files"][name]
            assert len(data) == expected["bytes"] and sha(data) == expected["sha256"]
            checked_bytes[name] = data
    old_receipt = json.loads(checked_bytes[f"{CELL}/receipt.json"])
    assert old_receipt["status"] == "FAIL"
    assert 'set(record["pinned_intermediates"])' in old_receipt["traceback"]

    target_checks = mismatches = rank_checks = 0
    checked_labels: set[tuple[str, int]] = set()
    one_failure = None
    for block in range(spec["blocks"]):
        for arm in ("control_a", "point_sum", "control_b"):
            prefix = f"{CELL}/b{block}_{arm}"
            with tempfile.TemporaryDirectory(prefix="pinned-x-rank-") as temp:
                paths = {}
                for suffix in ("base", "rank", "stdout"):
                    paths[suffix] = Path(temp) / f"trace.{suffix}.jsonl"
                    paths[suffix].write_bytes(checked_bytes[f"{prefix}.{suffix}.jsonl"])
                rank = verify_rank(paths["rank"], paths["base"], paths["stdout"])
                assert rank["status"] == "PASS" and rank["rank"] == spec["K"]
                rank_checks += 1
            base, = lines(checked_bytes[f"{prefix}.base.jsonl"])
            logs = lines(checked_bytes[f"{prefix}.rank.jsonl"])[-1]["logs"]
            targets = lines(checked_bytes[f"{prefix}.target.jsonl"])
            assert len(targets) == spec["L"]
            for index in TARGETS:
                record = targets[index]
                result = check_target(record, fixtures[index], base, logs,
                                      curve, generator, checked_labels)
                target_checks += 1
                if index == 763:
                    assert result["pinned_matches_selected"] is False
                    mismatches += 1
                    if one_failure is None:
                        one_failure = (record, fixtures[index], base, logs, result)
                else:
                    assert result["pinned_matches_selected"] is True

    assert target_checks == 30 and mismatches == 15 and one_failure is not None
    record, fixture, base, logs, witness = one_failure
    assert rejected(lambda: old_check_target(record, fixture, base, logs,
                                             curve, generator, set()))
    bad_pin = copy.deepcopy(record)
    bad_pin["pinned_intermediates"][0] = 1 << spec["n"]
    assert rejected(lambda: check_target(bad_pin, fixture, base, logs,
                                         curve, generator, set()))
    wrong_root = copy.deepcopy(record)
    wrong_root["pinned_intermediates"][0] ^= 1
    assert rejected(lambda: check_target(wrong_root, fixture, base, logs,
                                         curve, generator, set()))
    bad_log = copy.deepcopy(record)
    bad_log["recovered_scalar"] = (record["recovered_scalar"] + 1) % spec["subgroup_order"]
    assert rejected(lambda: check_target(bad_log, fixture, base, logs,
                                         curve, generator, set()))
    bad_point = copy.deepcopy(record)
    bad_point["point_indices"][0] = 0
    assert rejected(lambda: check_target(bad_point, fixture, base, logs,
                                         curve, generator, set()))
    return {"status": "PASS_CORRECTNESS_ONLY", "cell": CELL,
            "old_frozen_status": "FAIL", "timing_eligible": False,
            "archive_sha256": metadata["archive_sha256"],
            "fixture_sha256": spec["fixture_sha256"],
            "archive_members_rehashed": len(checked_bytes),
            "full_rank_traces_replayed": rank_checks,
            "selected_targets_replayed": target_checks,
            "selected_pinned_mismatches": mismatches,
            "example_new_witness": witness,
            "tamper_rejected": ["out_of_field_pin", "in_field_wrong_root",
                                "wrong_recovered_scalar", "wrong_point_index"]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, default=EVIDENCE)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = replay(args.evidence.resolve())
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        assert not args.out.exists(), "never overwrite a replay receipt"
        args.out.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
