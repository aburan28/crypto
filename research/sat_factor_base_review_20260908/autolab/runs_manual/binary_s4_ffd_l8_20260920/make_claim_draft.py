#!/usr/bin/env python3
"""Assemble the claim draft for the l=8 subspace S4 FFD measurement run.

Reads the run artifacts in this directory and emits
claim_draft_s4_subspace_ffd.json with every field required by the
ledger measurement schema (stage=decomposition) plus honest
measurement-only framing. Run after ffd_s4_subspace finishes.
"""
import hashlib
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[4]  # .../crypto
SUMMARY = HERE / "ffd_s4_subspace_summary.json"


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    s = json.loads(SUMMARY.read_text())
    summary = s["summary"]
    draws = s["draws"]
    cfg = s["config"]

    fix = summary["ffd_fix"]
    full = summary["ffd_full"]

    # D=7 deficit invariants for the x1-fixed variant (target-independence
    # diagnostic: geometric syzygies vs target-specific signal).
    deficits = []
    for d in draws:
        for m in d["fix"]["measured"]:
            if m["degree"] == 7:
                deficits.append(
                    {
                        "draw": d["draw"],
                        "decomposable": d["decomposable"],
                        "rows": m["rows"],
                        "rank": m["rank"],
                        "deficit": m["rows"] - m["rank"],
                    }
                )
    deficit_values = sorted({x["deficit"] for x in deficits})
    decomp_deficits = sorted({x["deficit"] for x in deficits if x["decomposable"]})
    nondecomp_deficits = sorted({x["deficit"] for x in deficits if not x["decomposable"]})

    hashes = (HERE / "source_and_binary_hashes.txt").read_text().splitlines()
    hash_map = {}
    for line in hashes:
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 64:
            hash_map[pathlib.Path(parts[1]).name] = parts[0]

    host_lines = (HERE / "host.txt").read_text().splitlines()
    host_id = {
        "node": host_lines[0] if host_lines else "unknown",
        "platform": " ".join(host_lines[1:3]) if len(host_lines) > 2 else "unknown",
        "ram_bytes": int(host_lines[3]) if len(host_lines) > 3 else None,
        "cpu": host_lines[4] if len(host_lines) > 4 else "unknown",
    }

    fixture_hash = sha256(SUMMARY)
    median_ms = summary["median_decompose_ms_pairs_baseline_recheck"]

    # Largest dense Macaulay matrix actually built (bytes).
    max_bytes = 0
    for d in draws:
        for variant in ("full", "fix"):
            for m in d[variant]["measured"]:
                words = (m["cols"] + 63) // 64
                max_bytes = max(max_bytes, m["rows"] * words * 8)

    draft = {
        "beat_id": "manual.binary.decomposition.l8_s4_subspace_ffd",
        "claim_boundary": "Public synthetic binary Semaev l=8 subspace FFD measurement only.",
        "claim_boundary_non_claims": [
            "not key recovery",
            "not sub-2^(2l) oracle claim",
            "not ledger promotion",
            "no oracle speedup measured",
            "no asymptotic or extrapolated statement",
        ],
        "created_at": s["created_at"],
        "dimension_l_or_dim": cfg["l"],
        "eq_var_ratio": "24 eqs / 24 x-vars (full); 24 eqs / 16 x-vars (x1-fixed)",
        "evidence_paths": [
            str((HERE / "ffd_s4_subspace_summary.json").relative_to(REPO)),
            str((HERE / "ffd_s4_subspace.stdout.txt").relative_to(REPO)),
            str((HERE / "ffd_s4_subspace.stderr.txt").relative_to(REPO)),
            str((HERE / "source_and_binary_hashes.txt").relative_to(REPO)),
            str((HERE / "host.txt").relative_to(REPO)),
            "examples/ffd_s4_subspace.rs",
            "src/cryptanalysis/ffd_harness.rs",
        ],
        "executable_or_source_hash": {
            "source_example": hash_map.get("ffd_s4_subspace.rs"),
            "source_harness": hash_map.get("ffd_harness.rs"),
            "binary": hash_map.get("ffd_s4_subspace"),
            "rustc": [l for l in hashes if l.startswith("rustc")],
        },
        "ffd_or_degree_of_regularity": {
            "status": "measured",
            "scope": "subspace_restricted_weil_descended_S4_l8_after_exact_e_elimination",
            "draws": cfg["draws"],
            "seed": cfg["seed"],
            "convention": "operational FFD: smallest D with rank<rows AND rank<cols on the multilinear truncated Macaulay matrix; identical convention to ffd_harness quadratic path (sanity-reproduced against documented S3 FFD=3, n=5..7, rank-for-rank)",
            "x1_fixed": {
                "vars": 16,
                "eqs": 24,
                "system_degree": 4,
                "ffd_min": fix["min"],
                "ffd_max": fix["max"],
                "ffd_mean": fix["mean"],
                "draws_with_fall": fix["draws_with_fall"],
                "draws_censored": fix["draws_censored"],
                "measured_degrees": [cfg["d_min"], cfg["dmax_fix"]],
                "koszul_degree": 8,
                "fall_structural": "all observed falls at D=7 < 2*deg_min=8, i.e. below the Koszul syzygy degree",
                "deficit_at_fall_values": deficit_values,
                "deficit_decomposable_targets": decomp_deficits,
                "deficit_nondecomposable_targets": nondecomp_deficits,
                "target_independence_note": (
                    "deficit invariant across targets including decomposable ones => syzygy space looks like target-independent image geometry, not a decomposability signal"
                    if deficit_values and decomp_deficits == nondecomp_deficits == deficit_values
                    else "deficit varies across draws; see per-draw records"
                ),
            },
            "full": {
                "vars": 24,
                "eqs": 24,
                "system_degree": 6,
                "fall": None,
                "measured_degrees": [cfg["d_min"], cfg["dmax_full"]],
                "censored_note": "no fall up to D=8 (rank==rows at D=6,7,8); D=9 dense estimate ~18GB skipped under mem cap",
                "koszul_degree": 12,
            },
            "sanity_gate": s["sanity_s3"],
            "verification_gates": {
                "elimination_exactness": "64 random V^3 points per draw, zero mismatches vs symmetrised_s4_eval",
                "witness_vanishing": "every decomposable draw's witness satisfies both variants; extra draws until >=1 witness check",
                "extra_witness_draws": len(s.get("extra_witness_draws", [])),
                "decomposable_draws": summary["decomposable_draws"],
            },
            "source": "examples/ffd_s4_subspace.rs",
        },
        "fixture_hash": fixture_hash,
        "fixture_hash_semantics": "sha256 of ffd_s4_subspace_summary.json (pins targets, seeds, per-draw measurements)",
        "host_id": host_id,
        "independent_replay_pointer": str((HERE / "claim_draft_s4_subspace_ffd.json").relative_to(REPO)),
        "independent_replay_command": "./target/release/examples/ffd_s4_subspace --draws 16 --seed 0x54FFD518 --dmax-full 8 --dmax-fix 7 --mem-cap-gb 2.0 --out <fresh-path>.json",
        "largest_solvable": {
            "l": cfg["l"],
            "n": cfg["n"],
            "budget": "measurement_run_no_oracle",
            "pairs_solve_recheck_s": round(median_ms / 1000.0, 4),
            "note": "no new oracle attempted; pairs-and-solve remains the decision procedure at this rung",
        },
        "ledger_gate_assessment": {
            "sub_2_to_2l_oracle": False,
            "ffd_logged": True,
            "median_le_half_pairs": "n/a_no_new_oracle",
            "targets_le_64": True,
            "note": (
                "FFD gate datum now measured (x1-fixed: min/max/mean over "
                + str(cfg["draws"])
                + " draws; full: censored no-fall record). Oracle gate unmet: per-X1 Macaulay at the fall degree costs ~30s (=> ~2h/target over the 2^8 sweep) and full-system D=8 LA costs ~135s/target, versus the 68ms ledger-recorded pairs baseline; the D=7 syzygy space does not currently translate into a decision or speedup."
            ),
        },
        "m_summands": 3,
        "median_ms_per_target": median_ms,
        "median_ms_per_target_note": "pairs-and-solve recheck on THIS host (warm, this binary); the ledger-recorded baseline arm measured 68.0 ms — runs must not be mixed",
        "n_or_bits": cfg["n"],
        "oracle_class": "none_ffd_measurement_only",
        "regime": "binary",
        "resource_caps": {
            "mem_cap_gb": cfg["mem_cap_gb"],
            "wall_clock_s_total": s["timing_s"]["total"],
            "max_dense_matrix_bytes_observed": max_bytes,
            "d9_full_skipped_estimate_bytes": "~1.8e10",
        },
        "schema_version": 2,
        "seeds": {
            "draw_seed": cfg["seed"],
            "extra_witness_seed": cfg.get("extra_witness_seed"),
            "s3_sanity_seed": s["sanity_s3"]["seed"],
            "rng": "splitmix64 (draws); rand StdRng (S3 sanity, matching run_sweep)",
            "decomp_bench_baseline_internal": "0x9E3779B9 (baseline arm, referenced not rerun)",
        },
        "stage": "decomposition",
        "system_degree": "degree<=6 full / degree<=4 x1-fixed, after exact e-elimination (semaev quadratic in e-space; correspondence linear/quadratic/cubic in x-bits)",
        "task_id": "TASK-IC-BOUNDARY-AUTOLAB-20260910",
        "unknowns": "24 x-bits (full); 16 x-bits (x1-fixed)",
        "verdict": summary["verdict"],
        "verdict_extended": (
            "X1-fixed system falls at D=7 (below Koszul degree 8) with deficit "
            + str(deficit_values)
            + " across draws; full system shows no fall up to D=8. Measurement only: no oracle, no speedup, no promotion."
        ),
        "forward_guidance": [
            "extract the D=7 syzygy module (kernel of the x1-fixed D=7 Macaulay transpose) and identify its origin; expected: target-independent relations induced by the 16-dim sigma-image inside the 37-dim (e2,e3) space",
            "if syzygies are target-independent geometry (as the deficit invariance suggests), they cannot decide decomposability directly; a sub-2^(2l) oracle of the per-X1 Groebner shape is cost-obstructed at l=8 (~10^3-10^5x baseline)",
            "the beat's sub-2^(2l) oracle remains open; remaining mechanism classes: structural compression of the pair sweep (resultants/batched gcd over X2), Frobenius-invariant factor bases (constant-factor, changes the comparison object), or higher summation polynomial shapes",
        ],
    }

    out = HERE / "claim_draft_s4_subspace_ffd.json"
    out.write_text(json.dumps(draft, indent=2, sort_keys=True) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
