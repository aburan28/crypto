#!/usr/bin/env python3
"""Measure Koblitz index calculus against rho at n=37 across all three target modes.

The autolab beat `koblitz.vs_rho.n37_wall` only exercises `target_mode=independent`.
This sweep runs the other two modes under identical parameters so the ledger can
quote the cheapest mode rather than the one the beat happens to pick.

Charged cost is taken from the producer's own batch summary row
(`full_algorithm_charged_total_ms`), which charges curve setup and support-table
setup once for the whole batch. Do not derive it by summing the per-target
`charged_total_ms` field: that field re-adds the shared `setup_ms` for every
target, which inflates the n=37 batch by roughly 14 ms per target.

Producer stdout is aggregated in-flight and never written to disk. The per-relation
records carry factor-base point coordinates, target point keys and walk
coefficients; only aggregates belong in a committed evidence directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

N = 37
A = 0
ETA = (1, 2)
PAIR_MODE = "signed_expanded"
QUERY_MODE = "pair_pair_16"
TARGET_MODES = ("independent", "coefficient_walk", "partition_walk")
RHO_QUOTIENT_MODE = "signed_frobenius"
RHO_BACKEND = "packed"

# Reused from autolab run 20260912T163230Z-3f7e563ca0 so the `independent` arm
# reproduces that recorded bundle exactly and acts as a cross-check.
DIRECT_SEED = 8866588946129576708
RHO_SEED = 14412094330580221154

DIRECT_SUMMARY_KIND = "retained_support_batch_summary"
RHO_ROW_KIND = "rho_public_fixture"

# Per-target fields aggregated from the direct producer's `relation_rank_summary`
# rows. These are cost components only; nothing here identifies a point.
DIRECT_ROW_FIELDS = (
    "collection_ms",
    "linear_solve_ms",
    "solution_validation_ms",
    "fixture_setup_ms",
    "fixture_generation_ms",
)
RHO_ROW_FIELDS = ("total_ms", "setup_ms", "walk_ms", "validation_ms", "walk_steps")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


class RowAccumulator:
    """Sums cost fields over rows of one kind without retaining the rows.

    The direct producer emits one record per accepted relation (~46k at n=37,
    ~53 MB of JSON), so rows are folded in as they stream rather than buffered.
    """

    def __init__(self, kind: str, fields: tuple[str, ...]) -> None:
        self.kind = kind
        self.fields = fields
        self.rows = 0
        self.totals = {field: 0.0 for field in fields}
        self.verified = 0

    def offer(self, row: dict) -> None:
        if row.get("kind") != self.kind:
            return
        self.rows += 1
        if row.get("verified"):
            self.verified += 1
        for field in self.fields:
            if field in row:
                self.totals[field] += float(row[field])

    def as_dict(self) -> dict:
        return {"rows": self.rows, "totals": self.totals}


def stream_rows(
    command: list[str], accumulators: list[RowAccumulator], keep_kinds: set[str]
) -> tuple[list[dict], float]:
    """Run `command`, folding JSONL stdout into `accumulators`.

    Only rows whose `kind` is in `keep_kinds` are returned; everything else is
    discarded once the accumulators have seen it.
    """
    started = time.perf_counter()
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=REPO_ROOT,
    )
    kept: list[dict] = []
    assert process.stdout is not None
    for line in process.stdout:
        line = line.strip()
        if not line:
            continue
        row = json.loads(line)
        for accumulator in accumulators:
            accumulator.offer(row)
        if row.get("kind") in keep_kinds:
            kept.append(row)
    stderr = process.stderr.read() if process.stderr else ""
    exit_code = process.wait()
    wall_ms = (time.perf_counter() - started) * 1000.0
    if exit_code != 0:
        raise SystemExit(f"producer failed ({exit_code}): {' '.join(command)}\n{stderr}")
    return kept, wall_ms


def run_direct(binary: Path, target_mode: str, fixtures: int) -> dict:
    command = [
        str(binary),
        str(N),
        str(A),
        str(ETA[0]),
        str(ETA[1]),
        str(DIRECT_SEED),
        PAIR_MODE,
        target_mode,
        QUERY_MODE,
        str(fixtures),
    ]
    per_target = RowAccumulator("relation_rank_summary", DIRECT_ROW_FIELDS)
    kept, wall_ms = stream_rows(command, [per_target], {DIRECT_SUMMARY_KIND})
    if len(kept) != 1:
        raise SystemExit(f"expected one batch summary, got {len(kept)}")
    summary = kept[0]
    charged_total_ms = float(summary["full_algorithm_charged_total_ms"])
    return {
        "arm": "direct",
        "target_mode": target_mode,
        "command": command,
        "fixtures": fixtures,
        "whole_process_wall_ms": wall_ms,
        "charged_total_ms": charged_total_ms,
        "charged_ms_per_target": charged_total_ms / fixtures,
        "batch_summary": summary,
        "per_target_totals": per_target.as_dict(),
        "all_fixtures_rank_plus_32": bool(summary["all_fixtures_rank_plus_32"]),
        "all_relations_group_verified": bool(summary["all_relations_group_verified"]),
    }


def run_rho(binary: Path, fixtures: int) -> dict:
    command = [
        str(binary),
        str(N),
        str(A),
        RHO_QUOTIENT_MODE,
        str(fixtures),
        RHO_BACKEND,
        str(RHO_SEED),
    ]
    per_target = RowAccumulator(RHO_ROW_KIND, RHO_ROW_FIELDS)
    _, wall_ms = stream_rows(command, [per_target], set())
    charged_total_ms = per_target.totals["total_ms"]
    return {
        "arm": "rho",
        "command": command,
        "fixtures": fixtures,
        "whole_process_wall_ms": wall_ms,
        "charged_total_ms": charged_total_ms,
        "charged_ms_per_target": charged_total_ms / fixtures,
        "per_target_totals": per_target.as_dict(),
        "targets_verified": per_target.verified,
        "all_targets_verified": per_target.verified == fixtures,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=int, default=1024)
    parser.add_argument(
        "--target-modes", nargs="+", default=list(TARGET_MODES), choices=TARGET_MODES
    )
    parser.add_argument("--out", type=Path, default=HERE)
    parser.add_argument(
        "--from-summaries",
        action="store_true",
        help="rebuild target_mode_sweep.json from summaries/ without rerunning producers",
    )
    args = parser.parse_args()

    direct_binary = REPO_ROOT / "target/release/examples/koblitz_rank_fixture"
    rho_binary = REPO_ROOT / "target/release/examples/koblitz_rho_fixture"
    summaries_dir = args.out / "summaries"

    if args.from_summaries:
        results = [
            json.loads((summaries_dir / f"direct.{mode}.json").read_text())
            for mode in args.target_modes
        ]
        rho = json.loads((summaries_dir / "rho.json").read_text())
        executable_sha256 = {
            "direct": sha256_file(direct_binary) if direct_binary.exists() else None,
            "rho": sha256_file(rho_binary) if rho_binary.exists() else None,
        }
    else:
        for binary in (direct_binary, rho_binary):
            if not binary.exists():
                raise SystemExit(
                    f"missing {binary}; build with "
                    "`cargo build --release --example koblitz_rank_fixture "
                    "--example koblitz_rho_fixture`"
                )
        summaries_dir.mkdir(parents=True, exist_ok=True)

        results = []
        for target_mode in args.target_modes:
            print(f"[sweep] direct target_mode={target_mode} ...", flush=True)
            result = run_direct(direct_binary, target_mode, args.fixtures)
            results.append(result)
            (summaries_dir / f"direct.{target_mode}.json").write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n"
            )
            print(
                f"[sweep]   charged {result['charged_ms_per_target']:.4f} ms/target, "
                f"wall {result['whole_process_wall_ms'] / 1000:.1f} s",
                flush=True,
            )

        print("[sweep] rho ...", flush=True)
        rho = run_rho(rho_binary, args.fixtures)
        (summaries_dir / "rho.json").write_text(
            json.dumps(rho, indent=2, sort_keys=True) + "\n"
        )
        print(
            f"[sweep]   charged {rho['charged_ms_per_target']:.4f} ms/target, "
            f"{rho['targets_verified']}/{args.fixtures} verified",
            flush=True,
        )
        executable_sha256 = {
            "direct": sha256_file(direct_binary),
            "rho": sha256_file(rho_binary),
        }

    fixtures = rho["fixtures"]
    rho_per_target = rho["charged_ms_per_target"]
    comparison = {
        "task_id": "TASK-IC-BOUNDARY-AUTOLAB-20260910",
        "n": N,
        "a": A,
        "eta": {"numerator": ETA[0], "denominator": ETA[1]},
        "fixtures": fixtures,
        "timing_class": "algorithmic_charged",
        "charged_definition": (
            "direct: full_algorithm_charged_total_ms from the producer batch summary "
            "(curve setup + support-table setup charged once, plus per-target fixture "
            "setup and relation collection). rho: sum of per-target total_ms "
            "(setup + walk + validation)."
        ),
        "claim_boundary": (
            "Public synthetic Koblitz fixtures only. Not key recovery, not asymptotic "
            "sub-rho, not imported points. A negative result for this implementation "
            "at this rung, not a proof that no crossover exists."
        ),
        "host_id": {
            "machine": platform.machine(),
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "executable_sha256": executable_sha256,
        "seeds": {"direct": DIRECT_SEED, "rho": RHO_SEED},
        "rho": {
            "charged_ms_per_target": rho_per_target,
            "whole_process_wall_ms": rho["whole_process_wall_ms"],
            "targets_verified": rho["targets_verified"],
            "all_targets_verified": rho["all_targets_verified"],
        },
        "direct_by_target_mode": {
            result["target_mode"]: {
                "charged_ms_per_target": result["charged_ms_per_target"],
                "whole_process_wall_ms": result["whole_process_wall_ms"],
                "rho_over_ic": rho_per_target / result["charged_ms_per_target"],
                "ic_times_behind_rho": result["charged_ms_per_target"] / rho_per_target,
                "all_fixtures_rank_plus_32": result["all_fixtures_rank_plus_32"],
            }
            for result in results
        },
        "ratio_note": (
            "rho_over_ic above 1.0 would be an index-calculus win. Every mode measured "
            "here is below 1.0."
        ),
        "cost_decomposition": {
            "note": (
                "Diagnostic only. The headline ratios above charge each arm its own "
                "full cost as the producers define it; nothing here is subtracted from "
                "a quoted ratio."
            ),
            "rho_ms_per_target": {
                field.removesuffix("_ms"): rho["per_target_totals"]["totals"][field]
                / fixtures
                for field in ("setup_ms", "walk_ms", "validation_ms", "total_ms")
            },
            "rho_setup_includes_target_construction": (
                "rho charges setup_ms for building the target point Q = d0*G and the "
                "jump table. Constructing Q is instance generation, which the direct "
                "arm reports separately as fixture_generation_ms and does not charge. "
                "This asymmetry runs against rho."
            ),
            "direct_ms_per_target": {
                result["target_mode"]: {
                    "batch_setup_amortized": (
                        float(result["batch_summary"]["curve_setup_ms"])
                        + float(result["batch_summary"]["support_setup_ms"])
                    )
                    / fixtures,
                    "collection": result["per_target_totals"]["totals"]["collection_ms"]
                    / fixtures,
                    "of_which_solution_validation": result["per_target_totals"][
                        "totals"
                    ]["solution_validation_ms"]
                    / fixtures,
                    "of_which_linear_solve": result["per_target_totals"]["totals"][
                        "linear_solve_ms"
                    ]
                    / fixtures,
                    "fixture_setup": result["per_target_totals"]["totals"][
                        "fixture_setup_ms"
                    ]
                    / fixtures,
                }
                for result in results
            },
            "direct_solution_validation_is_a_correctness_assertion": (
                "solution_validation_ms re-derives the discrete log of every "
                "factor-base representative with a scalar multiplication and replays "
                "every collected relation against the solved vector, under assert_eq!. "
                "It confirms an answer the linear solve already produced rather than "
                "producing it. It is the largest single component of the direct arm's "
                "charged cost and is the obvious place to look next; it is not "
                "subtracted from any ratio quoted here."
            ),
            "setup_amortization_is_not_the_lever": (
                "Curve plus support-table setup is a few tens of milliseconds for the "
                "whole batch, so it is well under 0.1 ms per target at 1024 targets. "
                "Amortizing the base across more targets cannot close the remaining "
                "gap; the cost is per-target relation collection."
            ),
        },
    }
    best = max(
        comparison["direct_by_target_mode"].items(),
        key=lambda item: item[1]["rho_over_ic"],
    )
    comparison["best_target_mode"] = best[0]
    comparison["verdict"] = (
        "NO_CHARGED_CROSSOVER_N37" if best[1]["rho_over_ic"] < 1.0 else "CHARGED_CROSSOVER_N37"
    )

    out_path = args.out / "target_mode_sweep.json"
    out_path.write_text(json.dumps(comparison, indent=2, sort_keys=True) + "\n")
    print(f"[sweep] wrote {out_path}")
    print(f"[sweep] best mode {best[0]}: rho/ic = {best[1]['rho_over_ic']:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
