#!/usr/bin/env python3
"""Check the pre-registered gates G1-G3 for each rung's native run.

Reads rung{0..3}_native.jsonl (one JSON object per line: 1,024 fixture records
and one summary) and the committed PR #955 R0 reference
../matched_rho_orbit_dlp_20260928_run/matched_arith_n53_L1024.jsonl.

  G1  every fixture verified and recovered == published (planted) scalar.
  G2  (rung 0 and 1) total_walk_steps == 19,103,507 and every per-fixture
      walk_steps equal to the PR #955 reference.
  G3  (rungs 2, 3) total_walk_steps within [0.85, 1.25] x 19,103,507.

Prints one JSON object per rung and exits non-zero if a gate fails.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(
    HERE, "..", "matched_rho_orbit_dlp_20260928_run", "matched_arith_n53_L1024.jsonl"
)
REF_TOTAL = 19_103_507
FIXTURES = 1024


def load(path):
    fixtures, summary = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            if obj.get("kind") == "rho_ks_batch_fixture":
                fixtures.append(obj)
            elif obj.get("kind") == "rho_ks_batch_summary":
                summary = obj
    return fixtures, summary


def main():
    ref_fix, ref_sum = load(REF)
    assert ref_sum["total_walk_steps"] == REF_TOTAL, "reference total changed"
    ref_steps = [f["walk_steps"] for f in ref_fix]
    ok = True
    for rung in range(4):
        path = os.path.join(HERE, f"rung{rung}_native.jsonl")
        if not os.path.exists(path):
            print(json.dumps({"rung": rung, "status": "missing"}))
            ok = False
            continue
        fix, summ = load(path)
        steps = [f["walk_steps"] for f in fix]
        g1 = (
            len(fix) == FIXTURES
            and summ is not None
            and summ.get("all_verified") is True
            and all(
                f["verified"] is True
                and f["recovered_fixture_scalar"] == f["published_fixture_scalar"]
                for f in fix
            )
        )
        total = summ["total_walk_steps"] if summ else None
        if rung in (0, 1):
            g2 = total == REF_TOTAL and steps == ref_steps
            g3 = None
        else:
            g2 = None
            g3 = total is not None and 0.85 * REF_TOTAL <= total <= 1.25 * REF_TOTAL
        gates = [g for g in (g1, g2, g3) if g is not None]
        passed = all(gates)
        ok &= passed
        per_target_ratio = (sum(steps) / sum(ref_steps)) if steps else None
        print(
            json.dumps(
                {
                    "rung": rung,
                    "fixtures": len(fix),
                    "G1_all_verified_and_equal_planted": g1,
                    "G2_bit_identical_to_pr955_walk": g2,
                    "G3_steps_within_band": g3,
                    "total_walk_steps": total,
                    "steps_over_reference": per_target_ratio,
                    "table_entries": summ.get("table_entries") if summ else None,
                    "cross_target_solves": summ.get("cross_target_solves") if summ else None,
                    "in_process_ms": summ.get("in_process_ms") if summ else None,
                    "rho_field_path": summ.get("rho_field_path") if summ else None,
                    "all_gates_passed": passed,
                }
            )
        )
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
