#!/usr/bin/env python3
"""Exact coverage screen for fixed three-point prefixes and a free pair oracle."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STUDY = "koblitz_n83_factor_base_sweep_20261008"
PRIMARY_ORDER = 2417851639230796216685689
INPUTS = {
    "design": ROOT / "design.json",
    "size_frontier": ROOT / "size-frontier-v2.json",
    "k2048_manifest": ROOT / "verification/remaining-budget-k2048-20261010/manifest.json",
    "k2048_upload": ROOT / "verification/remaining-budget-k2048-20261010/upload-receipt.json",
}
OUTPUT = ROOT / "pilot-01/fixed-prefix-pair-screen.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def ceil_div(numerator: int, denominator: int) -> int:
    return (numerator + denominator - 1) // denominator


def thresholds(order: int, points: int, columns: int) -> dict:
    if order <= 1 or points <= 0 or columns <= 0:
        raise ValueError("order, point count and column count must be positive")
    pairs = comb(points + 1, 2)
    return {
        "pair_multisets_with_repetition": str(pairs),
        "single_fixed_prefix_hit_ceiling": {
            "numerator": str(min(pairs, order)), "denominator": str(order),
        },
        "minimum_prefix_checks_for_half_single_hit_ceiling":
            ceil_div(order, 2 * pairs),
        "minimum_prefix_checks_for_half_rank_ceiling":
            max(columns, ceil_div(columns * order, 2 * pairs)),
    }


def calculate() -> dict:
    design, frontier, manifest, upload = (
        json.loads(INPUTS[name].read_text(encoding="utf-8"))
        for name in ("design", "size_frontier", "k2048_manifest", "k2048_upload")
    )
    if (design.get("schema") != "n83.factor-base-design/v1" or
            frontier.get("schema") != "n83.factor-base-design/v2-size-frontier" or
            manifest.get("schema") != "n83.factor-base-panel/v2-size-frontier" or
            upload.get("schema") != "n83.factor-base-storage/v2-size-frontier" or
            any(item.get("study") != STUDY for item in (design, frontier, manifest)) or
            upload.get("status") != "PASS" or
            frontier.get("parent_design_sha256") != sha256(INPUTS["design"]) or
            frontier.get("signed_orbit_point_count_conditional_on_successful_construction") != 166):
        raise ValueError("frozen design or retained-base receipt disagrees with this screen")
    order = int(frontier["subgroup_orders"]["0"])
    if order != PRIMARY_ORDER:
        raise ValueError("primary subgroup order changed")
    old_sizes = {item["parameter"] for item in design["base_specs"]
                 if item["a"] == 0 and item["closure"] == "signed_frobenius"}
    sizes = sorted(old_sizes | set(frontier["new_sizes"]))
    frontier_sizes = {item["parameter"] for item in frontier["base_specs"]
                      if item["a"] == 0 and item["closure"] == "signed_frobenius"}
    if (len(sizes) != 21 or sizes[-1] != 16627 or 2048 not in sizes or
            frontier_sizes != set(sizes)):
        raise ValueError("the retained primary signed-orbit size grid changed")
    rows = [row for row in manifest["bases"] if row.get("a") == 0 and
            row.get("columns") == 2048 and row.get("policy") == "public_x_hash" and
            row.get("seed") == 2026100801]
    if len(rows) != 1 or len(upload.get("objects", [])) != 1:
        raise ValueError("the retained K2048 object is absent or ambiguous")
    retained = rows[0]
    stored = upload["objects"][0]
    if (retained["points"] != 166 * 2048 or
            retained["compressed_blake3"] != stored.get("compressed_blake3") or
            retained["s3_uri"] != stored.get("s3_uri") or
            stored.get("downloaded_hash_matches") is not True):
        raise ValueError("the retained K2048 point count or S3 round trip disagrees")
    cases = []
    for columns in sizes:
        points = 166 * columns
        cases.append({
            "orbit_columns": columns,
            "point_count": points,
            "point_count_basis": "retained_k2048_manifest" if columns == 2048
                                 else "conditional_full_signed_orbit_construction",
            **thresholds(order, points, columns),
        })
    return {
        "schema": "n83.fixed-prefix-pair-screen/v1", "study": STUDY,
        "status": "EXACT_CONDITIONAL_BOUND", "curve_a": 0,
        "subgroup_order": str(order),
        "input_sha256": {name: sha256(path) for name, path in INPUTS.items()},
        "source_sha256": sha256(Path(__file__)),
        "model": {
            "query_law": "uniform whole prime-order subgroup, conditional on each fixed schedule",
            "prefix_policy": "three factor-base points fixed before observing the query",
            "oracle": "complete two-point membership and witness at zero cost",
            "pair_multisets": "C(B+1,2), including repetitions",
            "rank_rows_per_query": "at most one full-smooth relation row",
            "excluded": ["target-adaptive prefixes", "cross-prefix shared algebra",
                         "partial large-prime relations", "multiple rows per query"],
        },
        "proof": {
            "single_prefix": "Pr(hit) <= min(1,C(B+1,2)/r)",
            "A_checks_single_target": "Pr(any hit) <= min(1,A*C(B+1,2)/r)",
            "A_checks_across_queries": "Pr(rank >= K) <= min(1,A*C(B+1,2)/(r*K))",
            "half_threshold": "ceil(r/(2*C(B+1,2))) for one hit; max(K,ceil(K*r/(2*C(B+1,2)))) for K rank",
        },
        "cases": cases,
        "measurement_scope": "exact counting screen only; no solver search, relation yield or cold runtime",
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="create the immutable receipt")
    args = parser.parse_args()
    result = calculate()
    if args.write:
        with OUTPUT.open("x", encoding="utf-8") as target:
            json.dump(result, target, sort_keys=True, indent=2)
            target.write("\n")
    elif json.loads(OUTPUT.read_text(encoding="utf-8")) != result:
        raise ValueError("the frozen fixed-prefix receipt disagrees with current inputs or source")
    print(json.dumps({"status": result["status"], "cases": len(result["cases"]),
                      "receipt_sha256": sha256(OUTPUT)}, sort_keys=True))


if __name__ == "__main__":
    main()
