#!/usr/bin/env python3
"""Exact point-multiset support screen from already audited N83 receipts.

This performs integer arithmetic only. It does not open a factor-base object,
run a solver, or spend any retained-N83 search allowance.
"""

import hashlib
import json
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
RUN = STUDY / "verification/mixed-orbit-search-15000-20261012"
ORDER = 2417851639230796216685689
SOURCE_COMMIT = "f0a1e55f7cb6d965b0791ef68cb636ec8b8a3a1a"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def screen() -> dict:
    imported_path = RUN / "evidence-01/run/verify.stdout"
    audit_path = RUN / "audit.json"
    result_path = RUN / "evidence-01/run/run/result.json"
    imported = json.loads(imported_path.read_text())
    audit = json.loads(audit_path.read_text())
    result = json.loads(result_path.read_text())
    if not (
        imported["schema"] == "n83.mixed-orbit-import/v1"
        and imported["status"] == "PASS"
        and imported["curve_slug"] == "icv1-f2m83-tm6151469093347-debefd74"
        and imported["small_points"] == 196212
        and imported["residual_points"] == 143756
        and imported["small_signed_orbits"] == 1182
        and imported["residual_signed_orbits"] == 866
        and audit["status"] == "PASS"
        and audit["source_commit"] == SOURCE_COMMIT
        and audit["search_status"] == "UNKNOWN_search_cap"
        and audit["conservative_charge_seconds"] == 2400
        and result["status"] == "UNKNOWN_search_cap"
        and result["solver_conflicts"] == 5000
        and result["solver_models"] == 0
    ):
        raise ValueError("pinned mixed-orbit receipts disagree")
    small, residual = imported["small_points"], imported["residual_points"]
    cases = []
    for summands in (5, 6):
        for residual_summands in (0, 1, 2):
            small_summands = summands - residual_summands
            point_multisets = (
                comb(small + small_summands - 1, small_summands)
                * comb(residual + residual_summands - 1, residual_summands)
            )
            cases.append({
                "summands": summands,
                "small_summands": small_summands,
                "residual_summands": residual_summands,
                "point_multisets": str(point_multisets),
                "uniform_nonidentity_multiplicity_ceiling_numerator": str(point_multisets),
                "uniform_nonidentity_multiplicity_ceiling_denominator": str(ORDER - 1),
                "markov_hit_probability_ceiling": "1" if point_multisets >= ORDER - 1
                    else f"{point_multisets}/{ORDER - 1}",
                "worker_admitted": residual_summands in (1, 2),
            })
    return {
        "schema": "n83.mixed-orbit-support-screen/v1",
        "status": "ANALYTIC_NON_DISCRIMINATING",
        "scope": "point-multiset counts, not observed solver yield",
        "source_commit": SOURCE_COMMIT,
        "subgroup_order": str(ORDER),
        "small_points": small,
        "residual_points": residual,
        "small_signed_orbits": imported["small_signed_orbits"],
        "residual_signed_orbits": imported["residual_signed_orbits"],
        "import_receipt_sha256": sha(imported_path),
        "run_audit_sha256": sha(audit_path),
        "search_result_sha256": sha(result_path),
        "formula": "C(B+m-j-1,m-j)*C(L+j-1,j)",
        "cases": cases,
        "retained_n83_solver_search_executed": False,
        "additional_retained_n83_charge_seconds": 0,
        "selected_best_total_runtime": None,
    }


if __name__ == "__main__":
    print(json.dumps(screen(), indent=2, sort_keys=True))
