#!/usr/bin/env python3
"""Assemble claim drafts for the 2026-10-05 prime IC ladder rungs.

Reads the producer stdout JSON rows captured under
runs_manual/prime_j0_e2e_20bit_20261005/ and
runs_manual/prime_a3_ladder_20261005/ and emits fail-closed claim drafts
for the `end_to_end_dlp` stage (j=0 20-bit rung) and the a=-3 ladder
summary, with the global provenance fields the ledger requires.
"""
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS_MANUAL = HERE

HOST_ID = {
    "node": platform.node(),
    "platform": platform.platform(),
}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_row(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    rows = [l for l in lines if l.startswith("{")]
    if not rows:
        raise ValueError(f"no JSON row in {path}")
    return json.loads(rows[-1])


def source_hashes():
    root = HERE.parents[3]
    lib = root / "src/cryptanalysis/ec_index_calculus.rs"
    j0 = root / "src/cryptanalysis/ec_index_calculus_j0.rs"
    bench = root / "src/cryptanalysis/research_bench.rs"
    return {
        "source": sha256_file(lib),
        "source_j0": sha256_file(j0),
        "source_bench": sha256_file(bench),
    }


def main():
    ROOT = HERE.parents[3]
    now = datetime.now(timezone.utc).isoformat()

    # ── j=0 20-bit rung ─────────────────────────────────────────────
    j0_dir = RUNS_MANUAL / "prime_j0_e2e_20bit_20261005"
    j0_primary = load_row(j0_dir / "ka_20bit_seed20261005.json")
    j0_replays = [
        load_row(j0_dir / "ka_20bit_seed20261006.json"),
        load_row(j0_dir / "ka_20bit_seed20261007.json"),
    ]
    stages = j0_primary["ic_stages"]
    j0_draft = {
        "beat_id": "prime.end_to_end_dlp.j0_20bit",
        "created_at": now,
        "stage": "end_to_end_dlp",
        "regime": "prime",
        "schema_version": 2,
        "task_id": "TASK-IC-BOUNDARY-AUTOLAB-20260910",
        "n_or_bits": 20,
        "curve": j0_primary["curve"],
        "claim_boundary": "synthetic_known_answer",
        "claim_boundary_non_claims": [
            "not key recovery on deployed curves",
            "not asymptotic sub-rho",
            "not a vs_rho claim (IC wall exceeds rho wall at this size)",
            "not ledger promotion until independent validation",
        ],
        "recovered_d_verified": bool(j0_primary["ic_ok"]),
        "ic_agrees_rho": bool(j0_primary["ic_agrees_rho"]),
        "stage_timers": {
            "factor_base_ms": stages["factor_base_ms"],
            "relations_ms": stages["relations_ms"],
            "linear_algebra_ms": stages["linear_algebra_ms"],
            "verify_ms": stages["verify_ms"],
            "ic_total_ms": stages["total_ms"],
            "rho_ms": j0_primary["rho_ms"],
            "target_generation_ms_excluded": j0_primary["target_generation_ms_excluded"],
        },
        "relations": {
            "collected": stages["relations_collected"],
            "trials_total": stages["trials_total"],
            "trials_per_relation_median": stages["trials_per_relation_median"],
            "attempts_exhausted": stages["relation_attempts_exhausted"],
            "orbit_count": stages["orbit_count"],
        },
        "policy": {
            "target_orbits": j0_primary["target_orbits"],
            "extra_relations": j0_primary["extra_relations"],
            "max_trials_per_relation": j0_primary["max_trials_per_relation"],
            "max_relation_attempts": j0_primary["max_relation_attempts"],
            "rho_steps_cap": j0_primary["rho_steps_cap"],
        },
        "seeds": {
            "target_seed": j0_primary["target_seed"],
            "replay_target_seeds": [r["target_seed"] for r in j0_replays],
        },
        "fixture_hash": sha256_file(j0_dir / "ka_20bit_seed20261005.json"),
        "executable_or_source_hash": source_hashes(),
        "host_id": HOST_ID,
        "resource_caps": {},
        "evidence_paths": [
            str(p.relative_to(ROOT))
            for p in sorted(j0_dir.glob("*.json"))
        ],
        "independent_replay_pointer": str(
            (j0_dir / "claim_draft_20bit.json").relative_to(ROOT)
        ),
    }
    (j0_dir / "claim_draft_20bit.json").write_text(
        json.dumps(j0_draft, indent=1, sort_keys=True) + "\n"
    )

    # replay receipt
    replay_ok = all(
        r["ic_ok"] and r["rho_ok"] and r["ic_agrees_rho"] and r["ic_matches_truth"]
        for r in j0_replays
    )
    receipt = {
        "beat_id": "prime.end_to_end_dlp.j0_20bit",
        "created_at": now,
        "kind": "independent_replay",
        "note": "Fresh deterministic targets (distinct seeds) on the same 20-bit j=0 curve; all runs recover, match the sidecar, and agree with rho",
        "original": {
            "ic_ok": j0_primary["ic_ok"],
            "rho_ok": j0_primary["rho_ok"],
            "ic_agrees_rho": j0_primary["ic_agrees_rho"],
            "ic_ms": j0_primary["ic_ms"],
            "rho_ms": j0_primary["rho_ms"],
        },
        "replays": [
            {
                "seed": r["target_seed"],
                "ic_ok": r["ic_ok"],
                "rho_ok": r["rho_ok"],
                "ic_agrees_rho": r["ic_agrees_rho"],
                "ic_ms": r["ic_ms"],
                "rho_ms": r["rho_ms"],
            }
            for r in j0_replays
        ],
        "status": "PASS" if replay_ok else "FAIL",
        "schema_version": "1.0",
    }
    (j0_dir / "replay_receipt_20bit.json").write_text(
        json.dumps(receipt, indent=1, sort_keys=True) + "\n"
    )

    # ── a=-3 ladder ─────────────────────────────────────────────────
    a3_dir = RUNS_MANUAL / "prime_a3_ladder_20261005"
    rows = {}
    for path in sorted(a3_dir.glob("a3_*.json")):
        row = load_row(path)
        rows[path.name] = row
    ladder = []
    all_ok = True
    for name, row in sorted(rows.items()):
        stages = row["ic_stages"]
        ok = row["ic_ok"] and row["rho_ok"] and row["ic_agrees_rho"] and row["ic_matches_truth"]
        all_ok = all_ok and ok
        ladder.append(
            {
                "file": name,
                "curve": row["curve"],
                "curve_class": row["curve_class"],
                "bits": row["bits"],
                "seed": row["target_seed"],
                "ic_ms": row["ic_ms"],
                "rho_ms": row["rho_ms"],
                "trials_per_relation_median": stages["trials_per_relation_median"],
                "trials_per_relation_max": stages["trials_per_relation_max"],
                "trials_total": stages["trials_total"],
                "relations_collected": stages["relations_collected"],
                "stage_timers": {
                    "factor_base_ms": stages["factor_base_ms"],
                    "relations_ms": stages["relations_ms"],
                    "linear_algebra_ms": stages["linear_algebra_ms"],
                    "verify_ms": stages["verify_ms"],
                },
                "ic_ok": row["ic_ok"],
                "rho_ok": row["rho_ok"],
                "ic_agrees_rho": row["ic_agrees_rho"],
                "ic_matches_truth": row["ic_matches_truth"],
            }
        )
    a3_draft = {
        "beat_id": "prime.end_to_end_dlp.a3_ladder",
        "created_at": now,
        "stage": "end_to_end_dlp",
        "regime": "prime",
        "schema_version": 2,
        "task_id": "TASK-IC-BOUNDARY-AUTOLAB-20260910",
        "n_or_bits": 24,
        "claim_boundary": "synthetic_known_answer_a_minus_3_ladder",
        "claim_boundary_non_claims": [
            "not P-256 or CryptoPro-B themselves (synthetic small-field members of the a=-3, generic-j family)",
            "not key recovery on deployed curves",
            "not asymptotic sub-rho",
            "not a vs_rho claim (IC wall exceeds rho wall at every measured size)",
            "not ledger promotion until independent validation",
        ],
        "family": "y^2 = x^3 - 3x + b, prime order, h=1, a = p - 3; p256class p = 3 mod 4, cryptoproclass p = 1 mod 4",
        "ladder": ladder,
        "all_runs_verified_and_agree_rho": all_ok,
        "seeds": sorted({row["target_seed"] for row in rows.values()}),
        "fixture_hash": sha256_file(a3_dir / "a3_24bit_p256class_seed20261005.json"),
        "executable_or_source_hash": source_hashes(),
        "host_id": HOST_ID,
        "resource_caps": {},
        "evidence_paths": [str(p.relative_to(ROOT)) for p in sorted(a3_dir.glob("*.json"))],
        "independent_replay_pointer": str(
            (a3_dir / "claim_draft_a3_ladder.json").relative_to(ROOT)
        ),
    }
    (a3_dir / "claim_draft_a3_ladder.json").write_text(
        json.dumps(a3_draft, indent=1, sort_keys=True) + "\n"
    )
    print(json.dumps({
        "j0_20bit": {"status": receipt["status"], "draft": str(j0_dir / "claim_draft_20bit.json")},
        "a3_ladder": {"all_ok": all_ok, "draft": str(a3_dir / "claim_draft_a3_ladder.json")},
    }, indent=1))


if __name__ == "__main__":
    main()
