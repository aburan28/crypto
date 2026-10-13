#!/usr/bin/env python3
"""Exact source-only N83 mixed-orbit SAT-cap screen for v3 paired sizes."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = HERE.parents[1]
DESIGN = STUDY / "paired-frontier-v3.json"
OUTPUT = HERE / "screen.json"
N = 83
PHASE_BITS = 7
VARIABLE_CAP = 150_000
DOMAIN_CLAUSE_CAP = 12_000_000
SHAPES = ((5, 1), (5, 2), (6, 1), (6, 2))
DESIGN_SHA256 = "f7b860c9a652154541c4e6cf7d52ad13bd5a7436a41d38855d206cde36d4b0ac"
FROZEN_SOURCE_COMMIT = "f3d27f68b68cda907e0668042ac7fa13ae6b579d"
SOURCES = {
    "src/cryptanalysis/binary_semaev_chain_sat.rs": "0031972aab435d963af54a5f28f0510e45e5301f00c49466661191815becc06e",
    "src/cryptanalysis/koblitz_large_prime.rs": "accd5cf0ef7c738aedce38a21d172eeda1b309eaf9c8f4848030dac9d4421ddd",
    "examples/koblitz_n83_mixed_orbit_worker.rs": "e8781e8ba5c06528b1a8967312c5621084c0c3f8e3705e0882f849f4b1b9a152",
}
EVIDENCE = {
    "K600": {
        "audit": ("verification/k600-mixed-orbit-first-2400-20261012/audit.json",
                  "387d393b9805f35ee1ca6328e8602133a1b9bc036088b1be1a3ba5a29bb8469c"),
        "event": ("verification/k600-mixed-orbit-first-2400-20261012/evidence-01/run/run/events.jsonl",
                  "75552115ca9b34994553e6223969b331908131af212d5e5dfaa35b226f393e01"),
        "small": 600, "envelope": 2048,
    },
    "K1182": {
        "audit": ("verification/mixed-orbit-search-15000-20261012/audit.json",
                  "32420a62bb6ccf8f9df7f03dea9040e518606740968110dc8c0ac47cb50cafd5"),
        "event": ("verification/mixed-orbit-search-15000-20261012/evidence-01/run/run/events.jsonl",
                  "d9e2a43a232afe9e0d7d02dc0d32428395dd870897be772d5ca5b678c8a6f891"),
        "small": 1182, "envelope": 2048,
    },
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_blob_sha256(relative: str) -> str:
    blob = subprocess.run(
        ["git", "show", f"{FROZEN_SOURCE_COMMIT}:{relative}"],
        cwd=ROOT, capture_output=True, check=True,
    ).stdout
    return hashlib.sha256(blob).hexdigest()


def canonical(row: dict) -> bytes:
    return (json.dumps(row, indent=2, sort_keys=True) + "\n").encode()


def requirements(small: int, envelope: int, summands: int,
                 residual: int) -> dict:
    require(0 < small < envelope and 1 <= residual <= 2
            and 5 <= summands <= 6 and summands - residual >= 3,
            "unsupported mixed-orbit shape")
    small_terms = summands - residual
    suffix = envelope - small
    # These are the exact formulas in the frozen chain and orbit-domain
    # builders.  The first-pair ordering adds n-1 auxiliary variables.
    chain_variables = (N * (2 * summands - 1)
                       + 4 * N * (summands - 1)
                       + 3 * N * N * (summands - 1))
    constant_variables = -1 + N + PHASE_BITS + 2 * N * PHASE_BITS
    domain_variables = (small_terms * (2 * small + constant_variables)
                        + residual * (2 * suffix + constant_variables))
    constant_clauses = -3 + ((1 << PHASE_BITS) - N) + 4 * N * PHASE_BITS
    domain_clauses = (small_terms * ((N + 3) * small + constant_clauses)
                      + residual * ((N + 3) * suffix + constant_clauses))
    xor_rows = (5 * N * (summands - 1)
                + (PHASE_BITS + 1) * N * summands)
    variables = chain_variables + domain_variables + (N - 1)
    return {
        "small_orbits": small, "envelope_orbits": envelope,
        "summands": summands, "residual_summands": residual,
        "small_summands": small_terms,
        "chain_variables": chain_variables,
        "domain_variables": domain_variables,
        "first_pair_order_variables": N - 1,
        "total_variables": variables,
        "domain_clauses": domain_clauses,
        "xor_rows": xor_rows,
        "variable_cap_pass": variables <= VARIABLE_CAP,
        "domain_clause_cap_pass": domain_clauses <= DOMAIN_CLAUSE_CAP,
        "current_worker_cap_pass": (variables <= VARIABLE_CAP
                                    and domain_clauses <= DOMAIN_CLAUSE_CAP),
    }


def retained_check(label: str, spec: dict) -> dict:
    audit_path, audit_hash = spec["audit"]
    event_path, event_hash = spec["event"]
    require(sha256(STUDY / audit_path) == audit_hash, f"{label} audit changed")
    require(sha256(STUDY / event_path) == event_hash, f"{label} event changed")
    audit = json.loads((STUDY / audit_path).read_text())
    lines = (STUDY / event_path).read_text().splitlines()
    require(len(lines) == 1, f"{label} event count")
    event = json.loads(lines[0])
    stats = event["search"]
    calculated = requirements(spec["small"], spec["envelope"], 5, 2)
    require(audit["status"] == "PASS"
            and stats["patterns_attempted"] == 1
            and stats["solver_calls"] == 1
            and stats["peak_variables"] == calculated["total_variables"]
            and stats["domain_clauses"] == calculated["domain_clauses"]
            and stats["peak_xors"] == calculated["xor_rows"]
            and calculated["current_worker_cap_pass"],
            f"{label} exact retained counters differ")
    return {
        "audit_sha256": audit_hash, "event_sha256": event_hash,
        "small_orbits": spec["small"], "envelope_orbits": spec["envelope"],
        "observed_total_variables": stats["peak_variables"],
        "observed_domain_clauses": stats["domain_clauses"],
        "observed_xor_rows": stats["peak_xors"],
        "formula_matches_all_three_counters": True,
    }


def build() -> dict:
    require(sha256(DESIGN) == DESIGN_SHA256, "v3 design bytes changed")
    for relative, expected in SOURCES.items():
        require(frozen_blob_sha256(relative) == expected,
                f"frozen source changed: {relative}")
    design = json.loads(DESIGN.read_text())
    require(design["schema"] == "n83.factor-base-design/v3-paired-frontier"
            and design["pair_spec_count"] == 2970,
            "v3 pair population changed")
    signed = [pair for pair in design["pair_specs"]
              if pair["closure"] == "signed_frobenius"]
    require(len(signed) == 990, "signed-Frobenius pair population changed")
    repetitions = Counter((pair["small_orbits"], pair["envelope_orbits"])
                          for pair in signed)
    require(len(repetitions) == 55 and set(repetitions.values()) == {18},
            "signed pair size multiplicity changed")
    rows = [requirements(small, envelope, summands, residual)
            for small, envelope in sorted(repetitions)
            for summands, residual in SHAPES]
    summaries = []
    for summands, residual in SHAPES:
        shape_rows = [row for row in rows
                      if (row["summands"], row["residual_summands"])
                      == (summands, residual)]
        admitted = sum(row["current_worker_cap_pass"] for row in shape_rows)
        variable_failures = sum(not row["variable_cap_pass"] for row in shape_rows)
        clause_failures = sum(not row["domain_clause_cap_pass"] for row in shape_rows)
        summaries.append({
            "summands": summands, "residual_summands": residual,
            "size_pairs": len(shape_rows),
            "size_pairs_within_current_worker_caps": admitted,
            "signed_pair_specs_within_current_worker_caps": admitted * 18,
            "size_pairs_over_variable_cap": variable_failures,
            "size_pairs_over_domain_clause_cap": clause_failures,
        })
    observed = {label: retained_check(label, spec)
                for label, spec in EVIDENCE.items()}
    for small in (600, 1182):
        for summands, residual in SHAPES:
            require(requirements(small, 2048, summands, residual)
                    ["current_worker_cap_pass"],
                    "source-available mixed-orbit recipe exceeds frozen cap")
    return {
        "schema": "n83.paired-capacity-screen/v1",
        "status": "PASS_exact_source_cap_screen",
        "scope": "source_only_capacity_not_runtime",
        "v3_design_sha256": DESIGN_SHA256,
        "frozen_source_commit": FROZEN_SOURCE_COMMIT,
        "frozen_source_sha256": SOURCES,
        "screen_source_sha256": sha256(Path(__file__)),
        "limits": {"max_variables": VARIABLE_CAP,
                   "max_domain_clauses": DOMAIN_CLAUSE_CAP,
                   "order_first_pair": True},
        "signed_pair_spec_count": len(signed),
        "unique_signed_size_pairs": len(repetitions),
        "identical_specs_per_size_pair": 18,
        "shape_summaries": summaries,
        "size_pair_shape_rows": rows,
        "retained_counter_replay": observed,
        "sat_model_constructed": False,
        "sat_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = canonical(build())
    if args.write:
        OUTPUT.write_bytes(rendered)
        print(f"wrote {OUTPUT}")
    else:
        require(OUTPUT.read_bytes() == rendered, "saved screen differs")
        print("PASS: frozen paired N83 source-cap screen")


if __name__ == "__main__":
    main()
