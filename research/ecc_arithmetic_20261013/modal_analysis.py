"""Validate two immutable Modal results and export all paired observations."""

import csv
import hashlib
import io
import json
import statistics
from pathlib import Path


HERE = Path(__file__).resolve().parent
FIELDS = ("f2m113", "f2m131", "f2m163", "f2m233")
BASELINE_SHA = "2bc2a1f55887ef7db1fdb1ce5a6cd05c86b7c75e9eb8ade36497ff5686b37fd3"
CANDIDATE_SHA = "8b32c5467ee8c87752b9cb86594fb742aa9053fd5918172686c4270c9e936364"
CSV_FIELDS = ("run", "field", "round", "order", "arm", "iterations", "nanoseconds", "checksum")


def read_run(number):
    path = HERE / f"modal_result_attempt2_run{number}.json"
    raw = path.read_bytes()
    record = json.loads(raw)
    assert record["baseline_sha256"] == BASELINE_SHA
    assert record["candidate_sha256"] == CANDIDATE_SHA
    assert record["test_exit"] == 0 and record["bench_exit"] == 0
    assert "34 passed; 0 failed" in record["test_stdout"]
    rows = list(csv.DictReader(io.StringIO(record["bench_stdout"])))
    assert len(rows) == 192
    assert tuple(rows[0]) == CSV_FIELDS[1:]
    for row in rows:
        assert row["field"] in FIELDS
        assert int(row["iterations"]) == 500_000
        assert int(row["nanoseconds"]) > 0
        round_no = int(row["round"])
        assert 0 <= round_no < 24
        assert row["order"] == ("AB" if round_no % 2 == 0 else "BA")
        assert row["arm"] in ("baseline", "candidate")
        assert len(row["checksum"]) == 16
    summary = {
        "run": number,
        "raw_path": path.name,
        "raw_sha256": hashlib.sha256(raw).hexdigest(),
        "platform": record["platform"],
        "cpu_model": record["cpu_model"],
        "cpu_flags_relevant": record["cpu_flags_relevant"],
        "affinity": record["affinity"],
        "cargo_version": record["cargo_version"],
        "rustc_version": record["rustc_version"],
        "test_exit": record["test_exit"],
        "bench_exit": record["bench_exit"],
        "fields": {},
    }
    for field in FIELDS:
        pairs = {}
        for row in rows:
            if row["field"] == field:
                arm_rows = pairs.setdefault(int(row["round"]), {})
                assert row["arm"] not in arm_rows
                arm_rows[row["arm"]] = row
        assert set(pairs) == set(range(24))
        assert all(set(p) == {"baseline", "candidate"} for p in pairs.values())
        assert all(p["baseline"]["checksum"] == p["candidate"]["checksum"] for p in pairs.values())
        baseline = [int(pairs[i]["baseline"]["nanoseconds"]) for i in range(24)]
        candidate = [int(pairs[i]["candidate"]["nanoseconds"]) for i in range(24)]
        ratios = [b / c for b, c in zip(baseline, candidate)]
        summary["fields"][field] = {
            "pairs": 24,
            "baseline_median_ns": statistics.median(baseline),
            "candidate_median_ns": statistics.median(candidate),
            "baseline_min_ns": min(baseline),
            "baseline_max_ns": max(baseline),
            "candidate_min_ns": min(candidate),
            "candidate_max_ns": max(candidate),
            "paired_baseline_over_candidate_median": statistics.median(ratios),
            "paired_ratio_min": min(ratios),
            "paired_ratio_max": max(ratios),
        }
    return [{"run": number, **row} for row in rows], summary


def main():
    observations = []
    summaries = []
    for number in (1, 2):
        rows, summary = read_run(number)
        observations.extend(rows)
        summaries.append(summary)
    with (HERE / "modal_paired.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(observations)
    (HERE / "modal_summary.json").write_text(json.dumps({"runs": summaries}, indent=2) + "\n")
    print("verified 2 runs, 384 arm rows, 192 matching pairs, 34 tests per run")


if __name__ == "__main__":
    main()
