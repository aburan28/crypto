#!/usr/bin/env python3
"""Verify the Stage 187 audit against all referenced stage artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
GATES = STAGE.parent
audit_path = STAGE / "audit.json"
audit = json.loads(audit_path.read_text())
checks: list[tuple[str, bool]] = []


def check(name: str, ok: bool) -> None:
    checks.append((name, bool(ok)))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


check("schema", audit["schema"] == "koblitz_stage187_current_gate_audit.v1")
check("all_gates_false", audit["gates"]["all_seven_gates_passed"] is False)
check("full_cost_false", audit["full_cost_gate_passed"] is False)
check("external_false", audit["independent_reproduction_and_novelty_review_passed"] is False)
check("sota_false", audit["sota_claim"] is False)
check("dense_default", audit["selected_implementation"]["pair_selector_default"] is True)
check("no_subgroup_enumeration", audit["selected_implementation"]["target_subgroup_enumerated"] is False)
check("no_known_logs", audit["selected_implementation"]["discrete_log_labels_used"] is False)
for number, record in audit["stage_records"].items():
    root = GATES / record["directory"]
    result = root / "result.json"
    verification = root / "verification.json"
    check(f"result_sha:{number}", sha256(result) == record["result_sha256"])
    check(f"verification_sha:{number}", sha256(verification) == record["verification_sha256"])
    certificate = json.loads(verification.read_text())
    check(f"verification_clean:{number}", not certificate.get("failed"))
comments = STAGE / "external-comments.json"
check("external_snapshot_sha", sha256(comments) == audit["external_review_snapshot"]["snapshot_sha256"])
check("zero_external_comments", audit["external_review_snapshot"]["unaffiliated_comments"] == 0)
base = json.loads((GATES / audit["stage_records"]["174"]["directory"] / "result.json").read_text())
base_charge = base["campaign_accounting"]["measured_lower_bound_through_stage174"]
inc = audit["stage175_through_186_increment"]
cum = audit["measured_campaign_lower_bound_through_stage186"]
check("component_sum", cum["resource_components"] == base_charge["resource_components"] + inc["resource_components"])
check("wall_sum", cum["summed_wall_seconds"] == base_charge["summed_wall_seconds"] + inc["summed_wall_seconds"])
check("core_sum", cum["total_core_seconds"] == base_charge["total_core_seconds"] + inc["total_core_seconds"])
check("campaign_null", cum["complete_campaign_cost"] is None)
verification = {
    "schema": "koblitz_stage187_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "audit_sha256": sha256(audit_path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
