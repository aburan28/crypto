#!/usr/bin/env python3
"""Derive per-rung cost breakdowns from autolab run logs into a committed artifact.

`boundary_autolab.py` writes raw producer stdout to `autolab/runs/<id>/logs/`,
which `autolab/.gitignore` excludes: it is tens of megabytes and the per-relation
records carry factor-base point coordinates, target point keys and walk
coefficients. The cost components inside those logs are still worth citing, so
this script reduces them to aggregates and writes `rung_breakdown.json` next to
the promoted bundles.

Run it while the source logs are still on disk. Once they are cleaned the output
file is the only surviving record, which is the point of committing it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent.parent / "runs"

DIRECT_SUMMARY_KIND = "retained_support_batch_summary"
DIRECT_TARGET_KIND = "relation_rank_summary"
RHO_KIND = "rho_public_fixture"

DIRECT_TARGET_FIELDS = (
    "collection_ms",
    "solution_validation_ms",
    "linear_solve_ms",
    "setup_ms",
    "fixture_setup_ms",
    "charged_total_ms",
)
RHO_FIELDS = ("setup_ms", "walk_ms", "validation_ms", "total_ms", "walk_steps")


def fold(path: Path, kind: str, fields: tuple[str, ...]) -> tuple[int, dict, int]:
    totals = {field: 0.0 for field in fields}
    count = 0
    verified = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get("kind") != kind:
                continue
            count += 1
            if row.get("verified"):
                verified += 1
            for field in fields:
                if field in row:
                    totals[field] += float(row[field])
    return count, totals, verified


def first_of_kind(path: Path, kind: str) -> dict | None:
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get("kind") == kind:
                return row
    return None


def describe(run_dir: Path) -> dict:
    draft = json.loads((run_dir / "artifacts/claim_draft.json").read_text())
    logs = run_dir / "logs"
    direct_log = logs / "direct.stdout.jsonl"
    rho_log = logs / "rho.stdout.jsonl"

    targets, direct_totals, _ = fold(direct_log, DIRECT_TARGET_KIND, DIRECT_TARGET_FIELDS)
    batch = first_of_kind(direct_log, DIRECT_SUMMARY_KIND)
    rho_rows, rho_totals, rho_verified = fold(rho_log, RHO_KIND, RHO_FIELDS)

    entry: dict = {
        "run_id": run_dir.name,
        "beat_id": draft["beat_id"],
        "n": draft["n"],
        "timing_class": draft["timing_class"],
        "targets": targets,
        "ic_cost_ms": draft["ic_cost"],
        "rho_cost_ms": draft["rho_cost"],
        "rho_over_ic": draft["rho_cost"] / draft["ic_cost"],
        "direct_ms_per_target": {
            field.removesuffix("_ms"): value / targets
            for field, value in direct_totals.items()
        },
        "rho_ms_per_target": {
            field.removesuffix("_ms"): value / rho_rows
            for field, value in rho_totals.items()
        },
        "rho_targets_verified": f"{rho_verified}/{rho_rows}",
    }
    if batch is not None:
        entry["direct_batch_totals_ms"] = {
            field: batch[field]
            for field in (
                "curve_setup_ms",
                "support_setup_ms",
                "online_charged_ms",
                "full_algorithm_charged_total_ms",
            )
            if field in batch
        }
        entry["direct_batch_totals_ms"]["target_mode"] = batch.get("target_mode")
    else:
        # Single-fixture beats do not emit a batch summary row. The per-target
        # row's own setup_ms is then the whole support-table build, charged to
        # the one target rather than spread over a batch.
        entry["direct_batch_summary_absent"] = (
            "single-fixture run; per-target setup_ms is the full support-table "
            "build with nothing to amortize it against"
        )
    return entry


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE / "rung_breakdown.json")
    args = parser.parse_args()

    promoted = sorted((HERE / "autolab_runs").iterdir())
    entries = []
    missing = []
    for bundle in promoted:
        if not bundle.is_dir():
            continue
        source = RUNS / bundle.name
        if not (source / "logs/direct.stdout.jsonl").exists():
            missing.append(bundle.name)
            continue
        entries.append(describe(source))

    if missing:
        print(
            f"warning: no logs on disk for {', '.join(missing)}; "
            "these rungs are omitted",
            file=sys.stderr,
        )
    if not entries:
        raise SystemExit("no run logs found; nothing to derive")

    payload = {
        "derived_from": "research/sat_factor_base_review_20260908/autolab/runs/<run_id>/logs/*.jsonl",
        "note": (
            "Aggregates only. The source logs are gitignored because they are tens "
            "of megabytes and carry point-level records; this file is the committed "
            "reduction of their cost fields."
        ),
        "direct_charged_definition": (
            "Prefer direct_batch_totals_ms.full_algorithm_charged_total_ms, which "
            "charges curve and support-table setup once for the batch. The "
            "per-target charged field re-adds the shared setup for every target and "
            "must not be summed across a batch."
        ),
        "rho_over_ic_note": "above 1.0 would be an index-calculus win",
        "rungs": entries,
    }
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.out} ({len(entries)} rungs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
