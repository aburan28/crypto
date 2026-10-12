#!/usr/bin/env python3
"""Verify nested retained orbit prefixes and count finite large-prime supports.

This counts possible point multisets, not discovered partials, graph cycles,
independent matrix rows, or index-calculus runtime. The bounded outer command
charges even a timed-out or failed header check to the existing pilot budget.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
from decimal import Decimal

from support_moments import ORDERS, ORBIT_POINTS_PER_COLUMN, compact_ratio


ROOT = Path(__file__).resolve().parent
PANEL = ROOT / "pilot-01"
SCREEN = PANEL / "large-prime-pool-screen.json"
BUDGET = PANEL / "large-prime-pool-budget.json"
PAIRS = ((64, 256), (64, 600), (256, 600))
SUMMANDS = range(2, 7)


def multiset_count(points: int, size: int) -> int:
    if points < 0 or size < 0:
        raise ValueError("negative pool or summand count")
    if size == 0:
        return 1
    if points == 0:
        return 0
    return math.comb(points + size - 1, size)


def partial_support(base_points: int, residual_points: int, summands: int, residuals: int) -> int:
    """Unordered sums with exactly residuals from the disjoint residual pool."""
    if not 0 <= residuals <= summands:
        raise ValueError("residual count exceeds total summands")
    return multiset_count(base_points, summands - residuals) * multiset_count(
        residual_points, residuals
    )


def assert_prefix(small: dict, large: dict) -> None:
    for field in ("curve_a", "policy", "seed", "subgroup_order", "closure"):
        assert small[field] == large[field], f"nested header mismatch: {field}"
    assert small["orbit_columns"] < large["orbit_columns"]
    assert small["point_count"] == ORBIT_POINTS_PER_COLUMN * small["orbit_columns"]
    assert large["point_count"] == ORBIT_POINTS_PER_COLUMN * large["orbit_columns"]
    for field in ("representatives", "accepted_candidate_indices"):
        assert len(small[field]) == small["orbit_columns"]
        assert len(large[field]) == large["orbit_columns"]
        assert large[field][: small["orbit_columns"]] == small[field], (
            f"retained orbits are not a prefix: {field}"
        )


def header_for(row: dict) -> tuple[dict, str]:
    path = (PANEL / row["object"]).resolve(strict=True)
    assert PANEL.resolve() in path.parents, "object escaped panel directory"
    assert path.stat().st_size == row["bytes"], "compressed byte count changed"
    with gzip.open(path, "rb") as stream:
        header_bytes = stream.readline()
    header = json.loads(header_bytes)
    for key, expected in (
        ("curve_a", row["a"]),
        ("policy", row["policy"]),
        ("seed", row["seed"]),
        ("orbit_columns", row["columns"]),
        ("point_count", row["points"]),
    ):
        assert header[key] == expected, f"header/manifest mismatch: {key}"
    assert int(header["subgroup_order"]) == ORDERS[row["a"]]
    return header, hashlib.sha256(header_bytes).hexdigest()


def summarize() -> dict:
    manifest_bytes = (PANEL / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    replay = json.loads((PANEL / "replay.json").read_bytes())
    upload = json.loads((PANEL / "upload-receipt.json").read_bytes())
    prior_budget = json.loads((PANEL / "primary-cold-budget-audit.json").read_bytes())
    assert manifest["status"] == "completed_factor_base_panel"
    assert manifest["completed_base_count"] == len(manifest["bases"]) == 54
    assert replay["status"] == upload["status"] == "PASS"
    assert len(replay["checks"]) == len(upload["objects"]) == 54
    rows = {row["object"]: row for row in manifest["bases"]}
    assert len(rows) == 54
    assert {check["object"] for check in replay["checks"]} == set(rows)
    assert all(check["status"] == "PASS" for check in replay["checks"])
    assert {item["s3_uri"] for item in upload["objects"]} == {
        row["s3_uri"] for row in rows.values()
    }
    assert all(item["downloaded_hash_matches"] is True for item in upload["objects"])
    assert prior_budget["schema"] == "n83.primary-cold-active-budget-audit/v1"

    groups: dict[tuple[int, str, int], dict[int, tuple[dict, dict, str]]] = {}
    for row in rows.values():
        assert row["a"] in ORDERS and row["columns"] in (64, 256, 600)
        assert row["points"] == ORBIT_POINTS_PER_COLUMN * row["columns"]
        header, header_hash = header_for(row)
        group = (row["a"], row["policy"], row["seed"])
        sizes = groups.setdefault(group, {})
        assert row["columns"] not in sizes, "duplicate retained size"
        sizes[row["columns"]] = (row, header, header_hash)
    assert len(groups) == 18

    bindings = []
    for (arm, policy, seed), sizes in sorted(groups.items()):
        assert set(sizes) == {64, 256, 600}
        for smaller, larger in ((64, 256), (256, 600)):
            assert_prefix(sizes[smaller][1], sizes[larger][1])
        assert len({tuple(p) for p in sizes[600][1]["representatives"]}) == 600
        bindings.append(
            {
                "curve_a": arm,
                "policy": policy,
                "seed": seed,
                "objects": {
                    str(columns): {
                        "object": sizes[columns][0]["object"],
                        "s3_uri": sizes[columns][0]["s3_uri"],
                        "compressed_blake3_from_manifest": sizes[columns][0][
                            "compressed_blake3"
                        ],
                        "header_sha256": sizes[columns][2],
                        "point_set_blake3": sizes[columns][0]["point_set_blake3"],
                    }
                    for columns in (64, 256, 600)
                },
                "prefix_check": "PASS",
            }
        )

    cases = []
    for arm in sorted(ORDERS):
        order = ORDERS[arm]
        for base_columns, envelope_columns in PAIRS:
            base_points = ORBIT_POINTS_PER_COLUMN * base_columns
            residual_points = ORBIT_POINTS_PER_COLUMN * (
                envelope_columns - base_columns
            )
            for summands in SUMMANDS:
                cumulative = 0
                for residuals in range(3):
                    exact = partial_support(
                        base_points, residual_points, summands, residuals
                    )
                    cumulative += exact
                    cases.append(
                        {
                            "curve_a": arm,
                            "subgroup_order": str(order),
                            "base_columns": base_columns,
                            "envelope_columns": envelope_columns,
                            "base_points": base_points,
                            "residual_points": residual_points,
                            "summands": summands,
                            "exact_residuals": residuals,
                            "unordered_multisets_exact_residuals": str(exact),
                            "unordered_multisets_at_most_residuals": str(cumulative),
                            "uniform_nonidentity_target_hit_ceiling": {
                                "numerator": str(min(cumulative, order - 1)),
                                "denominator": str(order - 1),
                                "decimal": compact_ratio(
                                    min(cumulative, order - 1), order - 1
                                ),
                            },
                        }
                    )
    assert len(cases) == 2 * len(PAIRS) * len(SUMMANDS) * 3
    return {
        "schema": "n83.large-prime-pool-support-screen/v1",
        "study": manifest["study"],
        "status": "PASS_prefix_and_exact_support",
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "replay_panel_manifest_blake3": replay["panel_manifest_blake3"],
        "source_constructor_blake3": manifest["source_blake3"],
        "screen_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "target_model": "uniform nonidentity point of the stated prime-order subgroup",
        "pool_policy": "same arm, policy and seed; residual pool is envelope signed-Frobenius orbits after the base prefix",
        "counting_model": "unordered point multisets with repetition; zero, one or two points from a disjoint finite residual pool",
        "scope_limit": "first-moment support ceiling only; not partial discovery, graph cycle yield, independent rank, solver runtime or a factor-base winner",
        "prior_active_seconds_approx": prior_budget["total_active_seconds_approx"],
        "group_bindings": bindings,
        "cases": cases,
        "selected_best_total_runtime": None,
    }


def write_new(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as output:
        json.dump(value, output, indent=2, sort_keys=True)
        output.write("\n")


def bounded_main(wall_seconds: float) -> None:
    if SCREEN.exists() or BUDGET.exists():
        raise FileExistsError("bounded screen requires new immutable output paths")
    if not math.isfinite(wall_seconds):
        raise ValueError("screen wall cap must be finite")
    prior = json.loads((PANEL / "primary-cold-budget-audit.json").read_bytes())
    remaining = Decimal(prior["remaining_active_seconds_approx"])
    overhead_allowance = Decimal("1")
    wall_cap = Decimal(str(wall_seconds))
    if wall_cap <= 0 or wall_cap + overhead_allowance > remaining:
        raise ValueError("screen wall cap exceeds remaining authorized pilot budget")
    started = time.monotonic()
    status = "PRODUCER_FAILURE"
    exit_code = None
    error = None
    try:
        child = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--compute"],
            capture_output=True,
            text=True,
            timeout=wall_seconds,
            check=False,
        )
        exit_code = child.returncode
        if (
            child.returncode == 0
            and SCREEN.is_file()
            and json.loads(SCREEN.read_bytes()).get("status")
            == "PASS_prefix_and_exact_support"
        ):
            status = "PASS_prefix_and_exact_support"
        else:
            error = (child.stderr or child.stdout)[-2000:]
    except subprocess.TimeoutExpired:
        status = "UNKNOWN_wall_cap"
        error = "header screen exceeded the child process wall cap"
    except (OSError, ValueError) as failure:
        status = "PRODUCER_FAILURE_receipt"
        error = f"{type(failure).__name__}: {failure}"
    elapsed = Decimal(str(time.monotonic() - started))
    charged = elapsed + overhead_allowance
    total = Decimal(prior["total_active_seconds_approx"]) + charged
    if total > Decimal(prior["authorized_active_seconds"]):
        status = "PRODUCER_FAILURE_budget_overshoot"
    write_new(
        BUDGET,
        {
            "schema": "n83.large-prime-pool-active-budget/v1",
            "status": status,
            "prior_active_seconds_approx": prior["total_active_seconds_approx"],
            "child_wall_cap_seconds": str(wall_cap),
            "measured_outer_interval_seconds": str(elapsed),
            "conservative_process_overhead_allowance_seconds": str(
                overhead_allowance
            ),
            "charged_active_seconds": str(charged),
            "total_active_seconds_approx": str(total),
            "authorized_active_seconds": prior["authorized_active_seconds"],
            "child_exit_code": exit_code,
            "error": error,
            "screen_sha256": (
                hashlib.sha256(SCREEN.read_bytes()).hexdigest()
                if status == "PASS_prefix_and_exact_support"
                else None
            ),
            "selected_best_total_runtime": None,
        },
    )
    print(BUDGET)
    if status != "PASS_prefix_and_exact_support":
        raise RuntimeError(status)


def main() -> None:
    if sys.argv[1:] == ["--compute"]:
        write_new(SCREEN, summarize())
        print(SCREEN)
    elif len(sys.argv) == 3 and sys.argv[1] == "--bounded-seconds":
        bounded_main(float(sys.argv[2]))
    else:
        raise SystemExit(
            "usage: large_prime_pool_screen.py --bounded-seconds WALL | --compute"
        )


if __name__ == "__main__":
    main()
