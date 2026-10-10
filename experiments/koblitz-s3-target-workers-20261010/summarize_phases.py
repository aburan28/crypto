"""Compute paired phase ratios from the target-worker screen."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "analysis.json"
OUTPUT = HERE / "phase_summary.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def geomean(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(math.log(value) for value in values))


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    raw = SOURCE.read_bytes()
    analysis = json.loads(raw)
    assert analysis["verified_run_count"] == 16
    runs = {(row["label"], row["arm"]): row for row in analysis["runs"]}
    labels = [f"T{i}" for i in range(205, 213)]
    ratios = {}
    for metric in (*PHASES, "target_query_plus_pdp_charged", "online_ms"):
        paired = []
        for label in labels:
            baseline = runs[(label, "baseline")]
            candidate = runs[(label, "candidate")]
            if metric == "online_ms":
                earlier, later = baseline[metric], candidate[metric]
            elif metric == "target_query_plus_pdp_charged":
                earlier = sum(baseline["target_phase_ms"][phase]
                              for phase in ("target_query", "target_pdp_charged"))
                later = sum(candidate["target_phase_ms"][phase]
                            for phase in ("target_query", "target_pdp_charged"))
            else:
                earlier = baseline["target_phase_ms"][metric]
                later = candidate["target_phase_ms"][metric]
            paired.append(earlier / later)
        ratios[metric] = {"geometric_baseline_over_candidate": geomean(paired),
                          "per_target_baseline_over_candidate": paired}
    result = {
        "kind": "s3_n53_c1800_target_workers_pair_phase_summary_v1",
        "analysis_sha256": hashlib.sha256(raw).hexdigest(),
        "target_labels": labels,
        "phase_ratios": ratios,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: ratios[key]["geometric_baseline_over_candidate"]
                      for key in ("online_ms", "target_query_plus_pdp_charged",
                                  "target_descent")}, sort_keys=True))


if __name__ == "__main__":
    main()
