#!/usr/bin/env python3
"""Promote the n=83 a=1 vs_rho rung into the boundary ledger (fail-closed)."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLAIM = ROOT / "experiments/koblitz-single-target-n85-20261006/claim_report_vs_rho.json"
LEDGER = ROOT / "docs/ic/boundary_targets.json"


def main():
    if not CLAIM.exists():
        raise SystemExit(f"missing {CLAIM}")
    claim = json.loads(CLAIM.read_text())
    if claim.get("verdict") != "N85_A0_COMPACT_ORBIT_SINGLE_TARGET_ONLINE_IC_OVER_RHO_CROSSOVER":
        raise SystemExit("claim report verdict mismatch")
    if claim.get("median_online_speedup", 0) < 1.2:
        raise SystemExit("median online speedup below the win threshold")

    d = json.loads(LEDGER.read_text())
    records = d["regimes"]["koblitz"]["records"]
    vs = records["vs_rho"]
    old = vs["current"]
    vs.setdefault("history", []).insert(0, {
        "summary": old["summary"],
        "claim_boundary": old.get("claim_boundary"),
        "metrics": old.get("metrics", {}),
        "verdict": old.get("verdict"),
        "recorded": "2026-10-05-pre-n83-rung",
    })
    runs = claim["paired_runs"]
    vs["current"] = {
        "summary": (
            "Seventh primary single-target online rung and the a=0 arm past n=73: n=85 a=0 "
            "compact-orbit IC on u128 words (600 signed-Frobenius orbit columns, 102,000 base "
            "points, zero pair-table entries) beats automorphism-discounted rho on the identical "
            "frozen public point; three paired observations on one frozen target, median "
            f"{claim['median_online_speedup']:.1f}x online (IC 400-500 ms vs rho 114-278 s, "
            "range 227.4x-1,109.7x); deterministic 215,254-probe relation on all runs; both arms "
            "verified every run; standalone Python GF(2^85) replay PASS on all three. Precompute "
            "via the parallel guided rank (KIC_RANK_THREADS=12): ~9 min wall, mean ~4.2M "
            "probes/relation, 0 failures. Operational note: one launch was interrupted by an "
            "unrelated process-group kill after its rho arm completed; the recovery driver reran "
            "only the missing IC arm and launched R3, all receipts retained."
        ),
        "metrics": {
            "n": 85,
            "a": 0,
            "target_count": 1,
            "timing_class": "single_target_online",
            "subgroup_order": claim["subgroup_order"],
            "subgroup_log2_order": 53.7,
            "curve": claim["curve"],
            "paired_same_public_point": True,
            "public_generator": claim["public_generator"],
            "public_target_q": claim["public_target_q"],
            "fixture_scalar": claim["fixture_scalar"],
            "fixture_scalar_note": "validation-only sidecar; never supplied to either solver input",
            "base_point_set_digest": "e1b8e4a61202962d5c881775551dc9bd9dea72ce15a2e70a96aa68ea04336bd1",
            "factor_base_points": claim["precompute_detail"]["factor_base_points"],
            "orbit_columns": 600,
            "matrix_columns": 600,
            "domain": {
                "kind": "compact_orbit_u128",
                "pair_table_entries": 0,
                "edge_selectors": 0,
                "root_table_entries": 30598167,
                "regular_states": 30600000,
                "rank": 600,
                "rank_attempts": 600,
                "rank_failures": 0,
                "rank_policy": "guided parallel (KIC_RANK_THREADS=12): per-column scheduling-independent scalars; logs identical to the sequential path",
            },
            "automorphism_discount": "sqrt(2n) with A=2*n=170; rho producer uses signed_frobenius quotient",
            "paired_runs": runs,
            "median_online_speedup": claim["median_online_speedup"],
            "online_speedup_range_all_runs": claim["online_speedup_range_all_runs"],
            "ic_online_interval": (
                "after reusable target-independent base/index/log preparation (base load, S3 root "
                "index build, guided rank-600 log table) through target relation extraction, "
                "group-lift check, log recovery, and [d]G=Q check; one public Q supplied as "
                "solver input"
            ),
            "rho_online_interval": "after Q construction: per-target jump setup, walk, and in-process verification",
            "ic_precomputation_ms_excluded_from_online": claim["precompute_detail"]["rank_stage_ms_range"][1],
            "ic_verified_all_runs": True,
            "rho_verified_all_runs": True,
            "ic_and_rho_scalar_identical": True,
            "deterministic_ic_relation": {"probes": 215254, "same_all_runs": True},
            "independent_scalar_replay": {
                "method": "standalone Python GF(2^85) polynomial/affine-curve arithmetic independent of the Rust producers",
                "validated_runs": [r["run"] for r in runs],
                "status": "PASS",
            },
            "reproducibility": {
                "note": "rho walk seeds 202610071/072/073, ic rank seeds 31/32/33, KIC_RANK_THREADS=12; binaries pinned by sha256 in run.json",
            },
            "resource_caps": claim["resource_caps"],
            "asymptotic_sub_rho": False,
            "single_target_crossover": True,
        },
        "verdict": claim["verdict"],
        "claim_boundary": claim["claim_boundary"],
        "primary_single_target_status": "achieved_n41_n53_n61_n71_n73_n83a1_n85a0_20261008",
    }
    vs["next_target"] = {
        "summary": (
            "The n=59..131 scan shows the remaining rho-runnable IC-feasible rungs are n=89 a=0 "
            "(r ~ 2^58.x, cofactor 1,405,114,916; probes/relation ~ 2e8 at K=600, so K~850 with "
            "6.4e7 states ~ 2.2x the n=83 table — marginal online margin, land with care) and "
            "nothing wider: every other rung's r/(n^2 K^2) exceeds the materialized-table "
            "regime. Past that: a rho-projection-class rung at n=97 a=0 (r ~ 2^95) or the "
            "u128-ceiling a=1 rungs (n=107/109/113/127), then ECC2K-130 itself (n=131 a=0) "
            "which needs u256 words, an implicit S3 index, and 130-bit sparse LA; the 6-sum "
            "relation shape (probes/relation ~ 45r/(n^4 K^4), RESEARCH_ECC2K130_IC_FEASIBILITY.md "
            "section 4.7) is the top unexplored lever."
        ),
        "acceptance_gates": [
            "exactly one previously unseen online target",
            "IC and rho use identical frozen public point and matching resource envelope",
            "IC online clock starts after factor base, support index, and reusable logs are ready",
            "rho online clock starts at first target-dependent walk computation",
            "exclude process launch, input loading, and target fixture generation from both online intervals",
            "compact-orbit u128 domain: no explicit pair table",
            "guided rank accumulation (sequential or verified-parallel) establishes the reusable factor-log table",
            "record exclusive online phases, memory peak, and independent scalar replay for both arms",
            "report rho_online_ms / IC_online_ms; a win requires at least 1.20 and a verified answer",
            "multi-target or amortized rows remain secondary",
        ],
        "historical_secondary_note": (
            "The n=85 R2 launch was interrupted by an unrelated process-group kill after its rho "
            "arm completed; the recovery driver reran only the missing IC arm with identical "
            "frozen inputs and seeds, receipts retained per arm."
        ),
    }
    vs["evidence"] = [
        "experiments/koblitz-single-target-n85-20261006/claim_report_vs_rho.json",
        "experiments/koblitz-single-target-n83-20261006/frozen/fixture.json",
        "experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json",
        "experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json",
        "research/sat_factor_base_review_20260908/autolab/runs_manual/koblitz_parallel_rank_n73_20261005/",
    ]

    for priority in d.get("agent_priorities", []):
        if priority.get("regime") == "koblitz" and priority.get("stage") == "vs_rho":
            priority["beat"] = (
    "n=85 a=0 landed 2026-10-08 (median 693.4x online, three paired runs, Python replay PASS, "
    "deterministic 215,254-probe relation, a=0 arm past n=73, parallel-rank precompute). "
    "Remaining rho-runnable IC-feasible rungs: n=89 a=0 only (K~850, marginal). Next: a "
    "rho-projection-class rung at n=97 a=0 or the u128-ceiling a=1 rungs; ECC2K-130 itself "
    "needs u256 words + implicit S3 index + 130-bit sparse LA, and 4-sum total-work parity "
    "with rho is not reachable — the 6-sum relation shape is the top unexplored lever "
    "(RESEARCH_ECC2K130_IC_FEASIBILITY.md)."
)
    d["updated"] = "2026-10-06"
    LEDGER.write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\n")
    print("ledger koblitz vs_rho row promoted to n=83 a=1")


if __name__ == "__main__":
    main()
