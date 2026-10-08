#!/usr/bin/env python3
"""Promote the n=83 a=1 vs_rho rung into the boundary ledger (fail-closed)."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLAIM = ROOT / "experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json"
LEDGER = ROOT / "docs/ic/boundary_targets.json"


def main():
    if not CLAIM.exists():
        raise SystemExit(f"missing {CLAIM}")
    claim = json.loads(CLAIM.read_text())
    if claim.get("verdict") != "N83_A1_COMPACT_ORBIT_SINGLE_TARGET_ONLINE_IC_OVER_RHO_CROSSOVER":
        raise SystemExit("claim report verdict mismatch")
    if claim.get("median_online_speedup", 0) < 1.2:
        raise SystemExit("median online speedup below the win threshold")

    # Fail closed on the ledger's own vs_rho measurement schema
    # (docs/ic/PLAN_IC_ACCOUNTING_FIXES_20261007.md, F10).
    sys.path.insert(0, str(ROOT / "research/sat_factor_base_review_20260908/autolab"))
    import boundary_autolab as lab
    check = lab.validate_claim(claim, stage="vs_rho", ledger=lab.load_ledger(lab.load_protocol()))
    if check["status"] != "PASS":
        raise SystemExit("claim-check FAIL, not promoting: " + json.dumps({
            key: check[key]
            for key in ("missing_stage_fields", "missing_global_provenance", "pairing_errors")
        }))

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
            "Sixth primary single-target online rung and the first on the a=1 arm: n=83 a=1 "
            "compact-orbit IC on u128 words (600 signed-Frobenius orbit columns, 99,600 base "
            "points, zero pair-table entries) beats automorphism-discounted rho on the identical "
            "frozen public point; three paired observations on one frozen target, median "
            f"{claim['median_online_speedup']:.1f}x online (IC 5.6-7.9 s vs rho 43-71 s, "
            "range 5.5x-12.7x); deterministic 8,845,441-probe relation on all runs; both arms "
            "verified every run; standalone Python GF(2^83) replay PASS on all three. First "
            "ledger rung using the parallel guided rank (KIC_RANK_THREADS=12, logs identical "
            "to the sequential path): precompute 4.8 min wall vs the 5.1 h sequential n=73 "
            "rank. Margin smaller than n=73 because this rung's subgroup (2^52.9) is smaller "
            "than n=73's (2^56.3): rho is easier here; online walls also ran conservative "
            "under unrelated host load."
        ),
        "metrics": {
            "n": 83,
            "a": 1,
            "target_count": 1,
            "timing_class": "single_target_online",
            "subgroup_order": claim["subgroup_order"],
            "subgroup_log2_order": 52.9,
            "curve": claim["curve"],
            "paired_same_public_point": True,
            "public_generator": claim["public_generator"],
            "public_target_q": claim["public_target_q"],
            "fixture_scalar": claim["fixture_scalar"],
            "fixture_scalar_note": "validation-only sidecar; never supplied to either solver input",
            "base_point_set_digest": "d4a4da26d6e7cb09367de922d9806631fe49989f70d378cae76f347f8d76a9f7",
            "factor_base_points": claim["precompute_detail"]["factor_base_points"],
            "orbit_columns": 600,
            "matrix_columns": 600,
            "domain": {
                "kind": "compact_orbit_u128",
                "pair_table_entries": 0,
                "edge_selectors": 0,
                "root_table_entries": claim["precompute_detail"]["root_table_entries"],
                "regular_states": claim["precompute_detail"]["regular_states"],
                "rank": 600,
                "rank_attempts": 600,
                "rank_failures": 0,
                "rank_policy": "guided parallel (KIC_RANK_THREADS=12): per-column scheduling-independent scalars; logs identical to the sequential path",
            },
            "automorphism_discount": "sqrt(2n) with A=2*n=166; rho producer uses signed_frobenius quotient",
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
            "deterministic_ic_relation": {"probes": 8845441, "same_all_runs": True},
            "independent_scalar_replay": {
                "method": "standalone Python GF(2^83) polynomial/affine-curve arithmetic independent of the Rust producers",
                "validated_runs": [r["run"] for r in runs],
                "status": "PASS",
            },
            "reproducibility": {
                "note": "rho walk seeds 202610061/062/063, ic rank seeds 21/22/23, KIC_RANK_THREADS=12; binaries pinned by sha256 in run.json",
            },
            "resource_caps": claim["resource_caps"],
            "asymptotic_sub_rho": False,
            "single_target_crossover": True,
        },
        "verdict": claim["verdict"],
        "claim_boundary": claim["claim_boundary"],
        "primary_single_target_status": "achieved_n41_n53_n61_n71_n73_n83a1_20261006",
    }
    vs["next_target"] = {
        "summary": (
            "Extend the a=1 arm past n=83: n=97 a=0 (r ~ 2^95, cofactor 4) needs a rho "
            "projection-class rung (online rho ~ 2^44 steps is not wall-runnable); n=107/109/113 "
            "a=1 (r ~ 2^107-2^112, cofactor 2) and n=127 a=1 (r ~ 2^115, cofactor 7114) are the "
            "u128-ceiling rungs, also projection-class. ECC2K-130 itself (n=131 a=0, r ~ "
            "2^129.6) additionally needs u256 field words, an implicit S3 index, and 130-bit "
            "sparse LA; 4-sum total-work parity with rho is not reachable there (see "
            "RESEARCH_ECC2K130_IC_FEASIBILITY.md) — the 5-sum relation shape is the top "
            "unexplored lever."
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
            "The parallel guided rank is a precompute-path change only; its solved logs and the "
            "published target row are identical to the sequential path (verified byte-identical "
            "at n=71 and n=73, 2026-10-05). The n=83 rung's online walls are conservative: "
            "unrelated verification processes shared the host during the runs."
        ),
    }
    vs["evidence"] = [
        "experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json",
        "experiments/koblitz-single-target-n83-20261006/frozen/fixture.json",
        "experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json",
        "research/sat_factor_base_review_20260908/autolab/runs_manual/koblitz_parallel_rank_n73_20261005/",
    ]

    for priority in d.get("agent_priorities", []):
        if priority.get("regime") == "koblitz" and priority.get("stage") == "vs_rho":
            priority["beat"] = (
                "n=83 a=1 landed 2026-10-06 (median 8.2x online, three paired runs, Python replay "
                "PASS, first a=1 rung, first parallel-rank rung). Next: a rho-projection-class "
                "rung at n=97 a=0 (r ~ 2^95, cofactor 4) or the u128-ceiling a=1 rungs (n=107/109/113/127); "
                "ECC2K-130 itself (n=131 a=0) needs u256 words + implicit S3 index + 130-bit sparse LA, "
                "and 4-sum total-work parity with rho is not reachable — the 5-sum relation shape is "
                "the top unexplored lever (RESEARCH_ECC2K130_IC_FEASIBILITY.md)."
            )
    d["updated"] = "2026-10-06"
    LEDGER.write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\n")
    print("ledger koblitz vs_rho row promoted to n=83 a=1")


if __name__ == "__main__":
    main()
