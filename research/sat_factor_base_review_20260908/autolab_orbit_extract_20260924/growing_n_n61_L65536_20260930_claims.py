#!/usr/bin/env python3
"""Write vs_rho and end_to_end_dlp claim drafts for the n=61 L=65536 panel.

Only measured fields are filled.  The single-target online fields that the
vs_rho schema requires (candidate/workload identities, one-target online
intervals, IC phase costs, rho replay certificate, matched resource envelope)
are not produced by these batch producers and stay absent, so claim-check
fails closed for vs_rho by design.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = (Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "growing_n_n61_L65536_20260930").resolve()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    summary = json.loads((OUT / "growing_n_summary.json").read_text())
    host = json.loads((OUT / "host.json").read_text())
    replay = summary["replay"]
    blocks = summary["blocks"]
    scalars = OUT / "scalars_n61-ks-growing-65536-v1.txt"
    common = {
        "n_or_bits": {"n": 61, "a": 0, "subgroup_order_bits": 48},
        "fixture_hash": {
            "eval_corpus": "n61-ks-growing-65536-v1",
            "eval_scalars_sha256": sha256(scalars),
            "tune_corpus": "n61-ks-growing-tune-65536-v1",
            "tune_scalars_sha256": sha256(OUT / "scalars_n61-ks-growing-tune-65536-v1.txt"),
            "base_hash": blocks[0]["ic"]["summary"]["base_hash"],
        },
        "executable_or_source_hash": {
            "koblitz_orbit_dlp_fast_rs_sha256": sha256(ROOT / "examples/koblitz_orbit_dlp_fast.rs"),
            "koblitz_rho_batch_ks_v2_n61_rs_sha256": sha256(ROOT / "examples/koblitz_rho_batch_ks_v2_n61.rs"),
            "binaries_sha256sums": (OUT / "SHA256SUMS").read_text().strip(),
            "git_head": host["git_head"],
            "rustc": host["rustc"],
        },
        "host_id": f"{host['cpu']}, {host['logical_cpus']} logical CPUs, {host['memory_bytes']} bytes, {host['os']}",
        "resource_caps": {
            "threads": 1,
            "memory": "no cap; K candidates skipped when estimated IC RSS exceeded available memory",
            "isolation": "none: tools/isolated_bench.py is unsupported on macOS; per-run load recorded",
        },
        "seeds": {"batch_seed": 531310, "ic_seed": 7, "rho_dp_bits": 4},
        "claim_boundary_non_claims": [
            "public synthetic known-answer fixtures only",
            "no key recovery",
            "no asymptotic sub-rho claim; both arms scale as sqrt(L*r/n) at optimal K",
            "multi-target L=65536 batch: historical diagnostic, not the one-target primary comparison",
            "comparator is the historical KS v2 batched rho, superseded by the 2026-09-29 strong-rho ladder",
            "walls are not isolated (AGENTS.md section 10) and not admissible as a speedup",
        ],
    }
    end_to_end = {
        "stage": "end_to_end_dlp",
        "status": "PENDING_INDEPENDENT_VALIDATION",
        **common,
        "recovered_d_verified": {
            "targets_per_block": 65536,
            "blocks": len(blocks),
            "ic_recovered_matches_published_all": summary["all_targets_verified"],
            "rho_recovered_matches_published_all": all(b["ks"]["all_targets_verified"] for b in blocks),
            "same_target_points_all_blocks": all(b["same_target_points"] for b in blocks),
            "independent_replay": replay,
        },
        "stage_timers": {f"b{b['block']}": b["ic"]["summary"]["timing_ms"] for b in blocks},
        "claim_boundary": "Known-answer shared-log DLP on public synthetic a=0 n=61 fixtures, L=65536 per batch; every recovered log checked as [d]G = Q in-process and by the pure-Python replay.",
        "independent_replay_pointer": [
            str((OUT / f"replay_b{b['block']}/independent_replay.json").relative_to(ROOT)) for b in blocks
        ],
    }
    vs_rho = {
        "stage": "vs_rho",
        "status": "PENDING_INDEPENDENT_VALIDATION",
        **common,
        "comparator": summary["comparator"],
        "target_count": 65536,
        "timing_class": "whole_process_wall",
        "K": summary["K"],
        "wall_ratio_compact_over_rho": summary["wall_ratio"],
        "user_cpu_ratio_compact_over_rho": summary["user_ratio"],
        "instructions_retired_ratio_compact_over_rho": summary["instructions_ratio"],
        "ic_scalar_verified": summary["all_targets_verified"],
        "rho_scalar_verified": all(b["ks"]["all_targets_verified"] for b in blocks),
        "independent_validation": False,
        "verdict": "BATCH_DIAGNOSTIC_ONLY_NOT_SINGLE_TARGET",
        "claim_boundary": "L=65536 batched whole-process diagnostic against the historical KS v2 rho; ineligible for the single-target online vs_rho schema.",
        "independent_replay_pointer": end_to_end["independent_replay_pointer"],
    }
    (OUT / "claim_end_to_end_dlp.json").write_text(json.dumps(end_to_end, indent=1) + "\n")
    (OUT / "claim_vs_rho.json").write_text(json.dumps(vs_rho, indent=1) + "\n")
    print("wrote", OUT / "claim_end_to_end_dlp.json", OUT / "claim_vs_rho.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
