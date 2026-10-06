#!/usr/bin/env python3
"""Promote the n=73 vs_rho rung into the boundary ledger (fail-closed).

Requires: experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json
assembled from three PRODUCERS_COMPLETE runs with PASSing independent
replays.  Moves the current n=71 row to history and installs n=73.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLAIM = ROOT / "experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json"
LEDGER = ROOT / "docs/ic/boundary_targets.json"


def main():
    if not CLAIM.exists():
        raise SystemExit(f"missing {CLAIM}")
    claim = json.loads(CLAIM.read_text())
    if claim.get("verdict") != "N73_U128_COMPACT_ORBIT_SINGLE_TARGET_ONLINE_IC_OVER_RHO_CROSSOVER":
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
        "recorded": "2026-10-02-pre-n73-rung",
    })
    runs = claim["paired_runs"]
    vs["current"] = {
        "summary": (
            "Fifth primary single-target online rung: n=73 compact-orbit IC on u128 words "
            "(600 signed-Frobenius orbit columns, 87,600 base points, zero pair-table entries) beats "
            "automorphism-discounted rho on the identical frozen public point; three paired observations "
            "on one frozen target, median "
            f"{claim['median_online_speedup']:.1f}x online (IC {min(r['ic_online_ms'] for r in runs):.0f}-"
            f"{max(r['ic_online_ms'] for r in runs):.0f} ms vs rho "
            f"{min(r['rho_online_ms'] for r in runs)/1000:.0f}-{max(r['rho_online_ms'] for r in runs)/1000:.0f} s); "
            "all arms verified every run; standalone Python GF(2^73) replay PASS on all three; "
            "precompute (guided rank, mean 57-70M probes/relation, 5.1 h sequential) excluded from the "
            "online clock and fully logged"
        ),
        "metrics": {
            "n": 73,
            "target_count": 1,
            "timing_class": "single_target_online",
            "subgroup_order": claim["subgroup_order"],
            "subgroup_log2_order": 56.3,
            "curve": claim["curve"],
            "paired_same_public_point": True,
            "public_generator": claim["public_generator"],
            "public_target_q": claim["public_target_q"],
            "fixture_scalar": claim["fixture_scalar"],
            "fixture_scalar_note": "validation-only sidecar; never supplied to either solver input",
            "base_point_set_digest": "a65776025e4b9ebc4a1e074f9304d834ac12cd533afc7d89a1d3ee72a2a766cc",
            "factor_base_points": claim["precompute_detail"]["factor_base_points"],
            "orbit_columns": claim["precompute_detail"]["orbit_columns"],
            "matrix_columns": claim["precompute_detail"]["orbit_columns"],
            "domain": {
                "kind": "compact_orbit_u128",
                "pair_table_entries": 0,
                "edge_selectors": 0,
                "root_table_entries": claim["precompute_detail"]["root_table_entries"],
                "regular_states": claim["precompute_detail"]["regular_states"],
                "rank": claim["precompute_detail"]["rank"],
                "rank_attempts": claim["precompute_detail"]["rank"],
                "rank_failures": 0,
                "rank_policy": "guided: decompose [a]G - R_j for the first pivotless column j (sequential frozen arms)",
            },
            "automorphism_discount": "sqrt(2n) with A=2*n=146; rho producer uses signed_frobenius quotient",
            "paired_runs": runs,
            "median_online_speedup": claim["median_online_speedup"],
            "online_speedup_range_all_runs": claim["online_speedup_range_all_runs"],
            "ic_online_interval": (
                "after reusable target-independent base/index/log preparation (base load, S3 root index "
                "build, guided rank-600 log table) through target relation extraction, group-lift check, "
                "log recovery, and [d]G=Q check; one public Q supplied as solver input"
            ),
            "rho_online_interval": "after Q construction: per-target jump setup, walk, and in-process verification",
            "ic_precomputation_ms_excluded_from_online": claim["precompute_detail"]["rank_stage_ms_range"][1],
            "ic_verified_all_runs": True,
            "rho_verified_all_runs": True,
            "ic_and_rho_scalar_identical": True,
            "independent_scalar_replay": {
                "method": "standalone Python GF(2^73) polynomial/affine-curve arithmetic independent of the Rust producers",
                "validated_runs": [r["run"] for r in runs],
                "status": "PASS",
            },
            "reproducibility": {
                "note": "rho walk seeds 202610041/042/043, ic rank seeds 7/8/9; frozen binaries pinned by sha256 in run.json",
            },
            "resource_caps": claim["resource_caps"],
            "asymptotic_sub_rho": False,
            "single_target_crossover": True,
        },
        "verdict": claim["verdict"],
        "claim_boundary": (
            "Public synthetic Koblitz n=73 (a=0, b=1 over GF(2^73), 56.3-bit subgroup); exactly one "
            "previously unseen public target point solved online by IC and rho on the identical frozen point "
            "(sequential single-process arms, same host, no finite cap; peak RSS recorded per process); "
            "constant-factor win only. Not ECC2K-130 evidence, not an asymptotic sub-rho claim, not key "
            "recovery, not deployed-curve security impact. Multi-target amortized results remain secondary."
        ),
        "primary_single_target_status": "achieved_n41_n53_n61_n71_n73_20261005",
    }
    vs["next_target"] = {
        "summary": (
            "Extend the verified single-target online ladder past n=73 on the a=1 arm of the family: "
            "n=83, a=1 (r ~ 2^52.7, cofactor 1128547018, rho online ~1 min on this producer) with the "
            "compact-orbit u128 domain and the new parallel guided rank (KIC_RANK_THREADS), whose solved "
            "logs are identical to the sequential path (verified at n=73, 2026-10-05). Farther rungs: "
            "n=97 a=0 (r ~ 2^95, cofactor 4) needs a rho projection class; n=131 a=0 is the ECC2K-130 "
            "curve itself and needs u256 words plus an implicit S3 index (see RESEARCH_ECC2K130_IC_FEASIBILITY.md)."
        ),
        "acceptance_gates": [
            "exactly one previously unseen online target",
            "IC and rho use identical frozen public point and matching resource envelope",
            "IC online clock starts after factor base, support index, and reusable logs are ready",
            "rho online clock starts at first target-dependent walk computation",
            "exclude process launch, input loading, and target fixture generation from both online intervals",
            "compact-orbit u128 domain: no explicit pair table; peak RSS under the common cap",
            "guided rank accumulation establishes the reusable factor-log table in the precompute fixture",
            "record exclusive online phases, memory peak, and independent scalar replay for both arms",
            "report rho_online_ms / IC_online_ms; a win requires at least 1.20 and a verified answer",
            "multi-target or amortized rows remain secondary and cannot satisfy this gate",
        ],
        "historical_secondary_note": (
            "The n=73 rung's sequential arms ran while an unrelated verification process shared the host; "
            "the online walls are therefore conservative (inflated), not favorable. Precompute rank-stage "
            "walls are fully logged; the parallel guided rank is a precompute-path change with the online "
            "target row invariant (reproduces the frozen row exactly)."
        ),
    }
    vs["evidence"] = [
        "experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json",
        "experiments/koblitz-single-target-n73-20261003/frozen/fixture.json",
        "experiments/koblitz-single-target-n73-20261003/runs/IC1N73Ckb1fb87600PDP4rootRCguidedLAgaussTDdirectISO0h980b4cf20ad4We79c49b134eaR1/",
        "experiments/koblitz-single-target-n71-20261002/claim_report_vs_rho.json",
        "research/sat_factor_base_review_20260908/autolab/runs_manual/koblitz_parallel_rank_n73_20261005/",
    ]

    # Global agent priorities: rank-1 rung is now n=73-complete.
    for priority in d.get("agent_priorities", []):
        if priority.get("regime") == "koblitz" and priority.get("stage") == "vs_rho":
            priority["beat"] = (
                "n=73 landed 2026-10-05 (median "
                f"{claim['median_online_speedup']:.1f}x online, three paired runs, Python replay PASS). "
                "Next: n=83, a=1 (r ~ 2^52.7, cofactor 1128547018) on the compact-orbit u128 domain with "
                "the parallel guided rank; then a rho-projection-class rung at n=97 a=0 (r ~ 2^95, cofactor 4). "
                "ECC2K-130 itself (n=131, a=0) needs u256 words, an implicit S3 index, and 130-bit sparse LA; "
                "4-sum total-work parity with rho is not reachable (see RESEARCH_ECC2K130_IC_FEASIBILITY.md) — "
                "the 5-sum relation shape is the top unexplored lever."
            )
    d["updated"] = "2026-10-05"
    LEDGER.write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\n")
    print("ledger koblitz vs_rho row promoted to n=73")


if __name__ == "__main__":
    main()
