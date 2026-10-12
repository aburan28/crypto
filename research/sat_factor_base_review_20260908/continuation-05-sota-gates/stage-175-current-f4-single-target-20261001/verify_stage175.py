#!/usr/bin/env python3
"""Fail closed on the compact Stage 175 result and its raw artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
RESULT = json.loads((STAGE / "result.json").read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


checks: list[tuple[str, bool]] = []


def check(name: str, condition: bool) -> None:
    checks.append((name, bool(condition)))


check("schema", RESULT["schema"] == "koblitz_stage175_current_f4_single_target.v1")
check("decision_is_regression", RESULT["decision"]["status"] == "REJECTED_REGRESSION")
check("all_correctness_gates", all(item["pass"] for item in RESULT["correctness_checks"]))
check("algebraic_factor_base", RESULT["input"]["factor_base_algebraic"] is True)
check("no_subgroup_enumeration", RESULT["input"]["target_subgroup_enumerated"] is False)
check("no_known_log_labels", RESULT["input"]["discrete_log_labels_used"] is False)
check("conflicts_not_applicable", RESULT["candidate"]["conflicts"] is None)
check("single_core_not_measured", RESULT["candidate"]["single_core_seconds"] is None)

for section in ("build",):
    for receipt in RESULT[section]["artifacts"].values():
        path = STAGE / receipt["path"]
        check(f"artifact_exists:{receipt['path']}", path.is_file())
        check(f"artifact_size:{receipt['path']}", path.stat().st_size == receipt["bytes"])
        check(f"artifact_sha256:{receipt['path']}", sha256(path) == receipt["sha256"])
for group in ("candidate", "direct_mitm"):
    for run in RESULT[group]["runs"]:
        for receipt in run["artifacts"].values():
            path = STAGE / receipt["path"]
            check(f"artifact_exists:{receipt['path']}", path.is_file())
            check(f"artifact_size:{receipt['path']}", path.stat().st_size == receipt["bytes"])
            check(f"artifact_sha256:{receipt['path']}", sha256(path) == receipt["sha256"])

for group in ("candidate", "direct_mitm"):
    med = RESULT[group]["median"]
    runs = RESULT[group]["runs"]
    for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes"):
        observed = statistics.median(run["process"]["metrics"][field] for run in runs)
        check(f"median:{group}:{field}", observed == med[field])

prior = RESULT["stage174_reference_median"]
current = RESULT["candidate"]["median"]
for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes"):
    check(
        f"ratio:stage174:{field}",
        RESULT["ratios"]["candidate_over_stage174"][field] == current[field] / prior[field],
    )
check("wall_regressed", RESULT["ratios"]["candidate_over_stage174"]["wall_seconds"] > 1)
check("cpu_regressed", RESULT["ratios"]["candidate_over_stage174"]["total_core_seconds"] > 1)

verification = {
    "schema": "koblitz_stage175_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(STAGE / "result.json"),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
