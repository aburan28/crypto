"""Verify T110 striped-search controls and compare deterministic work."""

import gzip
import hashlib
import json
from pathlib import Path
import tarfile


HERE = Path(__file__).resolve().parent
FREEZE = HERE / "DIAGNOSTIC_FREEZE.json"
LEDGER = HERE / "local/raw_diagnostic2.jsonl.gz"
ARCHIVE = HERE / "local/raw_reports.tar.gz"
OLD_ANALYSIS = HERE.parent / "koblitz-s3-c1800-prefix-screen-20261009/analysis.json"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    plan = json.loads(FREEZE.read_bytes())
    old = json.loads(OLD_ANALYSIS.read_bytes())
    old_t110 = next(target for target in old["targets"] if target["label"] == "T110")
    raw_ledger = gzip.decompress(LEDGER.read_bytes())
    rows = [json.loads(line) for line in raw_ledger.splitlines()]
    assert len(rows) == 2
    by_arm = {}
    with tarfile.open(ARCHIVE, "r:gz") as archive:
        assert len(archive.getnames()) == 2
        for row in rows:
            arm = row["arm"]
            assert arm not in by_arm and arm in plan["arms"]
            assert row["label"] == "T110" and row["repetition"] == 1
            assert row["candidate_id"] == plan["arms"][arm]["candidate_id"]
            assert row["source_sha256"] == plan["source_sha256"]
            assert row["binary_sha256"] == plan["binary_sha256"]
            assert row["verified_native_and_fixture"] and row["phase_ok"]
            report = row["report"]
            assert report["target_search_schedule"] == plan["arms"][arm]["target_search_schedule"]
            assert report["target_search_lead_states"] == plan["arms"][arm]["target_search_lead_states"]
            assert report["target_search_stripes"] == plan["arms"][arm]["target_search_stripes"]
            assert report["factor_base_digest"] == plan["arms"][arm]["factor_base_digest"]
            assert report["prepared_rank"] == plan["arms"][arm]["prepare_rank"]
            name = f"T110_{arm}_R1.json"
            with archive.extractfile(name) as member:
                raw = member.read()
            assert sha(raw) == row["raw_sha256"] and json.loads(raw) == report
            by_arm[arm] = row
    baseline = by_arm["baseline"]["report"]
    stripe = by_arm["expanded"]["report"]
    assert baseline["recovered_scalar"] == stripe["recovered_scalar"]
    assert baseline["target_s3_calls"] == old_t110["baseline_target_s3_calls"]
    result = {
        "kind": "s3_n53_c1800stripe_T110_diagnostic_analysis_v1",
        "freeze_sha256": sha(FREEZE.read_bytes()),
        "ledger_decompressed_sha256": sha(raw_ledger),
        "raw_reports_archive_sha256": sha(ARCHIVE.read_bytes()),
        "verified_runs": 2,
        "baseline_target_s3_calls": baseline["target_s3_calls"],
        "old_c1800_prefix_target_s3_calls": old_t110["expanded_target_s3_calls"],
        "striped_c1800_target_s3_calls": stripe["target_s3_calls"],
        "old_c1800_over_striped_call_ratio": old_t110["expanded_target_s3_calls"] / stripe["target_s3_calls"],
        "c1200_over_striped_call_ratio": baseline["target_s3_calls"] / stripe["target_s3_calls"],
        "controlled_speedup": None,
        "fresh_target_inference": None,
    }
    (HERE / "diagnostic_analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
