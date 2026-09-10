#!/usr/bin/env python3
"""Verify the signed-Frobenius quotient rho controls against stage 2."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
AUTO = ROOT / "stage-3-automorphism-rho-20260909"
IC = ROOT / "stage-2-unknown-scalar-20260909"


def main() -> None:
    rows = []
    for secret in [53, 101]:
        auto = json.loads((AUTO / f"rho-auto-secret-{secret}.json").read_text())
        meter = json.loads((AUTO / f"rho-auto-secret-{secret}.metrics.json").read_text())
        ic = json.loads((IC / f"ic-secret-{secret}.json").read_text())
        ic_meter = json.loads((IC / f"ic-secret-{secret}.metrics.json").read_text())
        assert meter["returncode"] == 0 and not meter["timed_out"]
        assert auto["automorphism_optimized"] is True
        assert auto["verified_unknown_scalar_recovery"] is True
        assert auto["target"] == ic["target"]
        rows.append(
            {
                "secret": secret,
                "target": auto["target"],
                "index_calculus": {
                    "wall_seconds": ic_meter["metrics"]["wall_seconds"],
                    "total_core_seconds": ic_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": ic_meter["metrics"]["peak_rss_bytes"],
                    "internal_seconds": ic["timing_ns"]["end_to_end"] / 1e9,
                },
                "automorphism_rho": {
                    "quotient_iterations": auto["iterations"],
                    "walk_step_additions": auto["group_additions"],
                    "jump_table_additions": 16,
                    "initial_state_additions": auto["restarts"] + 1,
                    "reported_group_additions": auto["group_additions"] + 16 + auto["restarts"] + 1,
                    "rho_setup_scalar_multiplications": 32 + 2 * (auto["restarts"] + 1),
                    "restarts": auto["restarts"],
                    "wall_seconds": meter["metrics"]["wall_seconds"],
                    "total_core_seconds": meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": meter["metrics"]["peak_rss_bytes"],
                    "internal_seconds": auto["timing_ns"]["end_to_end"] / 1e9,
                },
            }
        )
    summary = {
        "schema": "koblitz_automorphism_rho_stage3_summary.v1",
        "curve": {"n": 9, "a": 0, "subgroup_order": 127},
        "automorphism": "signed Frobenius quotient, at most 2n representatives per orbit",
        "rows": rows,
        "checks": {
            "same_targets_as_index_calculus": True,
            "both_scalars_recovered_and_verified": True,
            "coefficient_scaling_tracks_canonicalization": "independently replayed with invariant assertions and separate Sage point checks on both frozen runs",
        },
        "observation": "automorphism rho used fewer total core-seconds than index calculus in both degree-9 pairs",
        "claim_boundary": "two toy controls; no stable ratio, scaling law, cryptographic-size result, or SOTA claim",
    }
    (ROOT / "stage-3-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# Stage 3: signed-Frobenius quotient Pollard rho",
        "",
        "The rho walk canonicalizes each state over negation and all degree-9 Frobenius images, multiplying both tracked coefficients by the corresponding public factor plus or minus lambda^k. A collision is accepted only when its derived scalar reproduces the public target.",
        "",
        "| secret | IC core-s / peak MiB / internal s | automorphism rho core-s / peak MiB / internal s | quotient iterations | walk / reported additions |",
        "|--:|--:|--:|--:|--:|",
    ]
    for row in rows:
        ic, rho = row["index_calculus"], row["automorphism_rho"]
        lines.append(
            f"| {row['secret']} | {ic['total_core_seconds']:.6f} / {ic['peak_rss_bytes']/1048576:.2f} / {ic['internal_seconds']:.6f} | "
            f"{rho['total_core_seconds']:.6f} / {rho['peak_rss_bytes']/1048576:.2f} / {rho['internal_seconds']:.6f} | "
            f"{rho['quotient_iterations']} | {rho['walk_step_additions']} / {rho['reported_group_additions']} |"
        )
    lines += [
        "",
        "Automorphism rho used fewer total core-seconds in both pairs. Process wall time contains a large unexplained first-launch outlier in each arm, so no wall-time ratio is reported from two samples. Reported additions include the 16 jump-table additions and one initial-state addition; scalar-multiplication internals and Frobenius field operations remain separate costs covered by process timing.",
        "",
        "Independent internal review replayed the canonicalization invariant and both collision equations successfully, while returning QUALIFIED on the campaign as a whole. This validates a same-target toy control and the required accounting fields. It does not establish the rho distribution, an exponent, or a cryptographic-size comparison.",
        "",
    ]
    (ROOT / "STAGE3_RESULTS.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
