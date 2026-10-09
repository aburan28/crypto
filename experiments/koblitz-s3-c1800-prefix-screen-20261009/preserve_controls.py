"""Preserve the frozen T025 correctness controls with exact source receipts."""

import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TMP = Path("/private/tmp")
OLD = ROOT / "experiments/koblitz-s3-compact-state-20261009/local/n53_report.json.gz"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PLAN = ROOT / "experiments/koblitz-s3-c1800-prefix-20261009/SCREEN_FREEZE.json"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    plan = json.loads(PLAN.read_bytes())
    targets = json.loads(PARENT.read_bytes())["targets"]
    target = next(item for item in targets if item["label"] == "T025" and item["n"] == 53)
    old = json.loads(gzip.decompress(OLD.read_bytes()))
    reports = {}
    artifacts = {}
    for arm, columns in (("baseline", 1200), ("expanded", 1800)):
        name = f"s3-c1800-prefix-T025-c{columns}"
        raw = (TMP / (name + ".json")).read_bytes()
        assert raw == (TMP / (name + ".stdout.json")).read_bytes()
        report = json.loads(raw)
        assert report["group_verified"] is True
        assert report["recovered_scalar"] == target["fixture_scalar"]
        assert report["target"] == target["public_point"]
        assert report["prepared_rank"] == columns
        assert report["target_count"] == 1
        assert report["factor_base_digest"] == plan["arms"][arm]["factor_base_digest"]
        assert report["target_search_schedule"] == plan["arms"][arm]["target_search_schedule"]
        assert report["target_search_prefix_columns"] == 1200
        reports[arm] = report
        output = HERE / "local" / (f"T025_{arm}.json.gz")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(gzip.compress(raw, mtime=0))
        artifacts[arm] = {"raw_sha256": sha(raw), "compressed_sha256": sha(output.read_bytes()),
                          "target_s3_calls": report["target_s3_calls"],
                          "online_ms": report["timing_ms"]["target_online_after_reusable_setup"]}
    skip = {"timing_ms", "target_search_schedule", "target_search_prefix_columns"}
    assert {k: v for k, v in old.items() if k not in skip} == {
        k: v for k, v in reports["baseline"].items() if k not in skip}
    receipt = {
        "kind": "s3_n53_c1800prefix_T025_correctness_control_v1",
        "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "point_sha256": target["point_sha256"],
        "recovered_scalar": target["fixture_scalar"],
        "baseline_exact_nontiming_reference_match": True,
        "native_and_fixture_verified": True,
        "arms": artifacts,
        "controlled_speedup": None,
    }
    (HERE / "control_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
