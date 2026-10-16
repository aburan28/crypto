#!/usr/bin/env python3
"""Source-only capacity bounds for every v4 cross-stream N83 pair recipe."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY))
sys.path.insert(0, str(HERE.parent / "paired-capacity-screen-20261013"))
import paired_frontier_v4 as v4  # noqa: E402
import screen  # noqa: E402


DESIGN = STUDY / "paired-frontier-v4.json"
PILOT = STUDY / "verification/cross-stream-pilot-import-source-20261013/census.json"
AUDIT = HERE.parent / "cross-stream-first-gate-source-20261013/audit.json"
RETAINED = HERE / "retained-first-result.json"
OUTPUT = HERE / "census.json"
FROZEN_SHA256 = {
    "paired-frontier-v4.json": "8e512b8dbc6577093f8e31f4afb625894cbdedd9f9f3e87dc4f36a7df52344c2",
    "paired_frontier_v4.py": "745d4b123fe1901b73eddcd59831d7ec44a266cf7f4da89b7a989814ebf6c516",
    "paired_frontier_v3.py": "edfe3be0be5d6fbcfde877d1955d35cc1c1b4db1718503d5f05b85cb0b3d14a8",
    "cross-stream-pilot/census.json": "d1a490e2e73552d96a5ee327e34367c5bb32fbc6f6d591eb9c0e5f029118a285",
    "paired-capacity-screen/screen.py": "116617e9de7552e48b548e632e80121508515e7079878305de0b975c8a924c5d",
    "cross-stream-first-gate/audit.json": "dfc96f11e4a399a42559847dd39e94b753f37911a54d806e0e54e830e7735478",
    "retained-first-result.json": "4a7a0d79417c868d71bf95e116eab7b92da1d87a911456d38ef23098512429dd",
}
SHAPES = screen.SHAPES
CAP_CLASSES = (
    "fits_even_at_zero_overlap",
    "depends_on_exact_overlap",
    "over_cap_for_every_nonempty_residual",
)
NOT_APPLICABLE = "signed_orbit_adapter_unavailable"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def canonical(row: dict) -> bytes:
    return (json.dumps(row, sort_keys=True, indent=2) + "\n").encode()


def bounds(base_orbits: int, pool_orbits: int, summands: int,
           residual_summands: int) -> dict:
    """Bound an unknown exact set difference with a nonempty residual domain.

    The pool can share at most min(base, pool) signed orbit columns. A residual
    domain of zero would not support a positive-residual mixed-orbit shape.
    Both worker formulas increase with its orbit count, so these endpoint
    calculations cover every possible nonempty exact import.
    """
    minimum = max(1, pool_orbits - base_orbits)
    maximum = pool_orbits
    low = screen.requirements(base_orbits, base_orbits + minimum,
                              summands, residual_summands)
    high = screen.requirements(base_orbits, base_orbits + maximum,
                               summands, residual_summands)
    if high["current_worker_cap_pass"]:
        cap_class = CAP_CLASSES[0]
    elif low["current_worker_cap_pass"]:
        cap_class = CAP_CLASSES[1]
    else:
        cap_class = CAP_CLASSES[2]
    return {
        "base_orbits": base_orbits,
        "pool_orbits": pool_orbits,
        "summands": summands,
        "residual_summands": residual_summands,
        "minimum_nonempty_residual_orbits": minimum,
        "maximum_residual_orbits": maximum,
        "minimum_total_variables": low["total_variables"],
        "maximum_total_variables": high["total_variables"],
        "minimum_domain_clauses": low["domain_clauses"],
        "maximum_domain_clauses": high["domain_clauses"],
        "xor_rows": high["xor_rows"],
        "cap_class": cap_class,
    }


def retained_check(design: dict, source_keys: set[tuple]) -> dict:
    audit = json.loads(AUDIT.read_text())
    result = json.loads(RETAINED.read_text())
    require(audit["raw_sha256"]["run/search/result.json"] == sha256(RETAINED),
            "retained result is not the audited raw receipt")
    require(audit["status"] == "PASS_receipts_with_censored_search"
            and audit["search_status"] == result["status"] == "UNKNOWN_search_cap"
            and audit["shared_signed_orbits"] == result["shared_signed_orbits"] == 0
            and audit["residual_signed_orbits"] == result["residual_signed_orbits"] == 600
            and audit["base_object_blake3"] == result["base_object_blake3"]
            and audit["pool_object_blake3"] == result["pool_object_blake3"],
            "retained audited pair fields changed")
    matches = []
    for index in range(design["cross_stream_pair_count"]):
        pair = v4.pair(design, index)
        base, pool = pair["base"], pair["pool"]
        if (pair["curve_a"] == 0
                and base["closure"] == pool["closure"] == "signed_frobenius"
                and base["policy"] == result["base_policy"]
                and pool["policy"] == result["pool_policy"]
                and base["seed"] == result["base_seed"]
                and pool["seed"] == result["pool_seed"]
                and base["orbits"] == result["base_signed_orbits"]
                and pool["orbits"] == result["pool_signed_orbits"]):
            matches.append((index, pair))
    require(len(matches) == 1, "retained pair identity is not unique in v4")
    index, pair = matches[0]
    require(v4.source_class(pair, source_keys)
            == "pilot_exact_import_gate_available_unexecuted",
            "retained pair source class changed")
    row = bounds(600, 600, 5, 2)
    require(row["cap_class"] == CAP_CLASSES[0]
            and result["search"]["peak_variables"] == row["maximum_total_variables"]
            and result["search"]["domain_clauses"] == row["maximum_domain_clauses"]
            and result["search"]["peak_xors"] == row["xor_rows"]
            and result["search"]["models"] == 0
            and result["search"]["conflicts"] == 5000,
            "retained zero-overlap SAT counters do not replay")
    return {
        "pair_index": index,
        "audit_sha256": FROZEN_SHA256["cross-stream-first-gate/audit.json"],
        "raw_result_sha256": FROZEN_SHA256["retained-first-result.json"],
        "shared_signed_orbits": 0,
        "observed_residual_orbits": 600,
        "observed_total_variables": result["search"]["peak_variables"],
        "observed_domain_clauses": result["search"]["domain_clauses"],
        "observed_xor_rows": result["search"]["peak_xors"],
        "cap_class": row["cap_class"],
        "search_status": result["status"],
    }


def build() -> dict:
    paths = {
        "paired-frontier-v4.json": DESIGN,
        "paired_frontier_v4.py": STUDY / "paired_frontier_v4.py",
        "paired_frontier_v3.py": STUDY / "paired_frontier_v3.py",
        "cross-stream-pilot/census.json": PILOT,
        "paired-capacity-screen/screen.py": Path(screen.__file__),
        "cross-stream-first-gate/audit.json": AUDIT,
        "retained-first-result.json": RETAINED,
    }
    for label, path in paths.items():
        require(sha256(path) == FROZEN_SHA256[label],
                f"frozen input changed: {label}")
    screen.build()  # Recheck its pinned Rust blobs and retained nested counters.
    design = json.loads(DESIGN.read_text())
    require(design["schema"] == "n83.factor-base-design/v4-cross-stream-frontier"
            and design["cross_stream_pair_count"] == 169884,
            "v4 pair population changed")
    source_keys = v4.pilot_keys(json.loads(PILOT.read_text()))
    # A size/shape calculation can be reused across distinct streams, but
    # every ordered pair is classified and counted individually below.
    size_rows = {
        (base, pool, summands, residual): bounds(base, pool, summands, residual)
        for base in v4.v3.EXPECTED_SIZES
        for pool in v4.v3.EXPECTED_SIZES
        for summands, residual in SHAPES
    }
    by_class_shape: dict[tuple, Counter] = defaultdict(Counter)
    source_counts: Counter = Counter()
    totals: dict[tuple, Counter] = defaultdict(Counter)
    for index in range(design["cross_stream_pair_count"]):
        pair = v4.pair(design, index)
        source_class = v4.source_class(pair, source_keys)
        source_counts[source_class] += 1
        for summands, residual in SHAPES:
            if source_class in (
                    "cross_closure_domain_unimplemented",
                    "non_signed_closure_domain_unimplemented"):
                cap_class = NOT_APPLICABLE
            else:
                row = size_rows[(pair["base"]["orbits"], pair["pool"]["orbits"],
                                 summands, residual)]
                cap_class = row["cap_class"]
            by_class_shape[(source_class, summands, residual)][cap_class] += 1
            totals[(summands, residual)][cap_class] += 1
    require(dict(source_counts) == design["source_class_pair_counts"],
            "source class population differs from frozen v4")
    require(sum(source_counts.values()) == 169884,
            "cross-stream pair census incomplete")
    require(all(sum(totals[shape].values()) == 169884 for shape in SHAPES),
            "shape census incomplete")
    summaries = [
        {
            "source_class": category,
            "summands": summands,
            "residual_summands": residual,
            "pair_count": source_counts[category],
            "cap_classes": {name: by_class_shape[(category, summands, residual)][name]
                            for name in (*CAP_CLASSES, NOT_APPLICABLE)},
        }
        for category in sorted(source_counts)
        for summands, residual in SHAPES
    ]
    return {
        "schema": "n83.cross-stream-source-capacity-census/v1",
        "status": "PASS_source_only_capacity_bounds",
        "scope": "all_v4_ordered_cross_stream_pair_recipes_four_current_worker_shapes",
        "frozen_input_sha256": FROZEN_SHA256,
        "current_worker_caps": {
            "variables": screen.VARIABLE_CAP,
            "domain_clauses": screen.DOMAIN_CLAUSE_CAP,
            "order_first_pair": True,
        },
        "ordered_pair_count": sum(source_counts.values()),
        "shape_count_per_pair": len(SHAPES),
        "source_class_pair_counts": dict(sorted(source_counts.items())),
        "shape_cap_totals": [
            {"summands": summands, "residual_summands": residual,
             "cap_classes": {name: totals[(summands, residual)][name]
                             for name in (*CAP_CLASSES, NOT_APPLICABLE)}}
            for summands, residual in SHAPES
        ],
        "source_class_shape_counts": summaries,
        "size_shape_bounds": [size_rows[key] for key in sorted(size_rows)],
        "retained_first_pair_replay": retained_check(design, source_keys),
        "exact_overlap_measured_for_all_pairs": False,
        "sat_search_executed_by_census": False,
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
        require(OUTPUT.read_bytes() == rendered, "saved capacity census changed")
        print("PASS: frozen v4 cross-stream source-capacity census")


if __name__ == "__main__":
    main()
