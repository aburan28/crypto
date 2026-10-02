#!/usr/bin/env python3
"""Compose the Stage 186 dense-pair elimination-policy screen."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def load_run(root: Path, variant: str, order: int) -> dict:
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    process = json.loads(metrics_path.read_text())
    report = json.loads(stdout_path.read_text())
    extra = report["cost"]["extra"]
    full = variant == "full"
    correct = (
        process["returncode"] == 0
        and not process["timed_out"]
        and report["status"] == "unsat"
        and report["exhaustive"] is True
        and report["source_instance_verified"] is True
        and report["regenerated_source_exact"] is True
        and report["solver_equations_blake3"] == EXPECTED_EQUATIONS
        and report["fixed_x1_masks_visited"] == 512
        and report["fixed_x1_systems_completed"] == 242
        and extra["pair_dense_select_calls"] == 1_011_275
        and extra["pair_quadratic_select_calls"] == 0
        and extra["pair_candidate_visits"] == 1_137_001_812
        and extra["pair_lcm_groups"] == 594_604_504
        and extra["pair_cover_lookups"] == 698_657_372
        and extra["pair_dense_scratch_bytes_max"] == 4_472_832
        and report["cost"]["ops"] == (318_635_818_320 if full else 319_313_687_585)
        and extra["word_xors_performed"] == (99_192_937_526 if full else 147_794_583_858)
        and extra["full_m4ri_matrices"] == (723 if full else 0)
        and extra["basis_len"] == 1
        and extra["pairs_left"] == 0
        and report["conflicts"] is None
    )
    return {
        "order": order,
        "variant": variant,
        "correct": correct,
        "process": process,
        "report": report,
        "artifacts": {
            "metrics": receipt(metrics_path),
            "stdout": receipt(stdout_path),
            "stderr": receipt(stderr_path),
        },
    }


root = STAGE / "development" / "screen"
blocktables = load_run(root / "01-blocktables", "blocktables", 1)
full = load_run(root / "02-full", "full", 2)
if not blocktables["correct"] or not full["correct"]:
    raise SystemExit("screen correctness failure")
bm = blocktables["process"]["metrics"]
fm = full["process"]["metrics"]
ratio = {
    "wall": fm["wall_seconds"] / bm["wall_seconds"],
    "core": fm["total_core_seconds"] / bm["total_core_seconds"],
    "rss": fm["peak_rss_bytes"] / bm["peak_rss_bytes"],
    "performed_xors": 99_192_937_526 / 147_794_583_858,
}
continued = ratio["core"] < 1 and ratio["performed_xors"] < 1
new_processes = [blocktables["process"], full["process"]]
new_charge = {
    "components": 2,
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in new_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in new_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in new_processes),
}
stage185 = json.loads(
    (STAGE.parent / "stage-185-dense-pair-default-replay-20261001" / "result.json").read_text()
)
result = {
    "schema": "koblitz_stage186_dense_pair_elimination_policy.v1",
    "claim_boundary": "One opened n=59 target and a selected-stack elimination screen only; not a full index-calculus run or SOTA.",
    "source": {
        "selection_commit": "8014149a2b55cd2cca202237a302643e39a50f6e",
        "binary_sha256": "3aefd5cd0daf73602bf8a38e7fdcb363e630ffd57bfebe47a37123b1d03429b1",
        "supplied_lock_sha256": "4f17b356fa7bac392b6d801d1c74fb9e36b6517f9465c8ebc19bb9a2792a84c5",
    },
    "inherited_stage185_result_sha256": sha256(
        STAGE.parent / "stage-185-dense-pair-default-replay-20261001" / "result.json"
    ),
    "screen": {
        "runs": [blocktables, full],
        "full_over_blocktables": ratio,
        "continued": continued,
    },
    "decision": {
        "status": "REJECTED_SCREEN_CPU_REGRESSION",
        "confirmation_required": continued,
        "selected_phase_b_elimination": "current_five_column_blocktables",
        "selected_pair_selector": "dense_exact",
        "reason": "Full M4RI reduced performed XORs but increased both wall and total core-seconds after dense pair selection, so confirmation was not run.",
    },
    "new_campaign_charge": new_charge,
    "cumulative_with_stage185": {
        "components": stage185["campaign_charge"]["components"] + 2,
        "wall_seconds_sum": stage185["campaign_charge"]["wall_seconds_sum"] + new_charge["wall_seconds_sum"],
        "total_core_seconds_sum": stage185["campaign_charge"]["total_core_seconds_sum"] + new_charge["total_core_seconds_sum"],
        "peak_rss_bytes_max": max(stage185["campaign_charge"]["peak_rss_bytes_max"], new_charge["peak_rss_bytes_max"]),
    },
    "single_core_seconds": None,
    "conflicts": None,
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 186: elimination policy after dense-pair selection",
    "",
    "Both arms used the selected dense pair updates, authenticated the same target, and returned exhaustive UNSAT with their exact elimination counters.",
    "",
    "| candidate / control | wall | core | RSS | performed XORs |",
    "|---|---:|---:|---:|---:|",
    f"| full M4RI / current BlockTables | {ratio['wall']:.6f} | {ratio['core']:.6f} | {ratio['rss']:.6f} | {ratio['performed_xors']:.6f} |",
    "",
    "Full M4RI reduces actual XORs and RSS but increases both wall and CPU on the selected dense-pair stack. The screen continuation condition fails, so no confirmation panel is run.",
    "",
    f"The two screen processes charge {new_charge['wall_seconds_sum']:.6f} wall seconds, {new_charge['total_core_seconds_sum']:.6f} core-seconds, and {new_charge['peak_rss_bytes_max']} bytes maximum RSS.",
    "",
    "Current five-column BlockTables plus dense exact pair selection remains the Phase B target-specific and repository-default configuration. This is not a full attack or SOTA result.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"ratio": ratio, "decision": result["decision"], "new_campaign_charge": new_charge}, indent=2))
