"""Post hoc audit of paired run order, faults, and target phase timings."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run(root: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    ledger_raw = (root / "ledger.jsonl").read_bytes()
    analysis = json.loads((root / "analysis.json").read_bytes())
    assert analysis["ledger_sha256"] == sha(ledger_raw)
    rows = [json.loads(line) for line in ledger_raw.splitlines()]
    assert len(rows) == 16 and analysis["all_complete"]
    pairs = []
    for offset in range(0, len(rows), 2):
        pair = []
        for entry in rows[offset:offset + 2]:
            receipt_raw = (root / entry["relative_receipt"]).read_bytes()
            assert sha(receipt_raw) == entry["receipt_sha256"]
            receipt = json.loads(receipt_raw)
            packed = (root / entry["relative_receipt"]).parent / "report.json.gz"
            packed_raw = packed.read_bytes()
            meta = receipt["archive"]["report.json"]
            assert sha(packed_raw) == meta["gzip_sha256"]
            report_raw = gzip.decompress(packed_raw)
            assert sha(report_raw) == meta["raw_sha256"]
            timing = json.loads(report_raw)["timing_ms"]
            pair.append({
                "arm": entry["arm"],
                "online_ms": entry["online_ms"],
                "target_query_ms": timing["target_query"],
                "target_point_validation_ms": timing["target_point_validation"],
                "target_pdp_charged_ms": timing["target_pdp_charged"],
                "target_descent_ms": timing["target_descent"],
                "major_page_faults_entire_process":
                    receipt["resource_counters"]["major_page_faults"],
                "receipt_sha256": entry["receipt_sha256"],
            })
        assert rows[offset]["label"] == rows[offset + 1]["label"]
        assert {item["arm"] for item in pair} == {"baseline", "candidate"}
        pairs.append({"label": rows[offset]["label"], "first": pair[0],
                      "second": pair[1],
                      "second_over_first_online":
                          pair[1]["online_ms"] / pair[0]["online_ms"]})
    result = {
        "kind": "s3_n53_c3200_prefetch_posthoc_order_audit_v1",
        "audit_source_sha256": sha(Path(__file__).read_bytes()),
        "analysis_sha256": sha((root / "analysis.json").read_bytes()),
        "ledger_sha256": sha(ledger_raw),
        "pairs": pairs,
        "second_slower_count": sum(p["second_over_first_online"] > 1 for p in pairs),
        "first_zero_major_fault_count": sum(
            p["first"]["major_page_faults_entire_process"] == 0 for p in pairs),
        "second_positive_major_fault_count": sum(
            p["second"]["major_page_faults_entire_process"] > 0 for p in pairs),
        "second_over_first_online_geomean": math.exp(statistics.mean(
            math.log(p["second_over_first_online"]) for p in pairs)),
        "posthoc_diagnostic": True,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in (
        "second_slower_count", "first_zero_major_fault_count",
        "second_positive_major_fault_count", "second_over_first_online_geomean")},
        sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.run_dir, args.output)
