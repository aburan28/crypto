#!/usr/bin/env python3
"""Operation-count accounting for the Koblitz single-target vs_rho rungs.

The retained single-target rungs report exploratory wall-clock ratios. The IC probe
loop and the rho fixture run at different per-operation speeds, so a wall
ratio mixes the algorithmic comparison with implementation speed.  This
tool reads the raw run rows that every rung already retains and writes an
``operation_accounting`` block into the rung's ``claim_report_vs_rho.json``:

- IC target probes per paired run (deterministic per target) and the
  rank-stage mean probes per relation on its separate guided-query law;
- rho walk steps per paired run and the expected walk
  ``sqrt(pi*r/(4n))`` on signed-Frobenius classes;
- native counter quotients (rho steps per IC probe) next to the wall
  ratios, each arm's native counter rate, the rank-to-target probe ratio,
  thread counts, and the host-contention flags each rung's README states.

Probes and walk steps remain distinct native units.  The block keeps the
counter quotients as diagnostics; a calibrated operation speedup and total
work S remain unknown until a common boundary is measured.

See docs/ic/PLAN_IC_ACCOUNTING_FIXES_20261007.md, findings F2-F4 and steps
0.3, 1.1, 1.4, 1.5.

Usage:
    python3 op_accounting.py                 # print the table for all rungs
    python3 op_accounting.py --write         # also update the claim reports
    python3 op_accounting.py --rung n83 --write
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[3]
SCHEMA_VERSION = "2.0"
GENERATED_BY = "research/sat_factor_base_review_20260908/autolab/op_accounting.py"

UNIT_ASSUMPTION = (
    "IC probes (S3 quadratic solve plus root-table lookup) and rho walk steps "
    "(group addition plus class canonicalisation) are distinct native counters; "
    "no common operation calibration or complete total-work boundary is recorded"
)

NON_CLAIMS = [
    "not compared against Pollard rho with precomputation (Bernstein-Lange "
    "distinguished-point tables) at equal precompute and memory",
    "single frozen target: IC online probes are deterministic per target, so "
    "repeated paired runs sample timing noise, not the target distribution",
    "wall-clock speedup includes implementation speed; operation counts are in "
    "operation_accounting",
]
CONTENTION_NON_CLAIM = (
    "runs flagged host_contended in operation_accounting shared the host with "
    "unrelated processes and are exploratory"
)

N73_RUN = "IC1N73Ckb1fb87600PDP4rootRCguidedLAgaussTDdirectISO0h980b4cf20ad4We79c49b134ea"

# Paired ledger runs per rung.  ``contended`` repeats what each rung's own
# README states; nothing here is inferred from timings.
RUNGS: dict[str, dict[str, Any]] = {
    "n61": {
        "dir": "experiments/koblitz-single-target-n61-20260925",
        "n": 61,
        "subgroup_order": 162888033982417,
        "runs": [
            {"run": "R1", "path": "runs/R1", "contended": False},
            {"run": "R2", "path": "runs/R2", "contended": False},
            {"run": "R3", "path": "runs/R3", "contended": False},
            {"run": "R4", "path": "runs/R4", "contended": False},
        ],
        "rank_summaries": ["runs/R1/ic-summary.json", "runs/R2/ic-summary.json", "runs/R3/ic-summary.json"],
        "contention_source": None,
    },
    "n71": {
        "dir": "experiments/koblitz-single-target-n71-20261002",
        "n": 71,
        "subgroup_order": 5513228015079457,
        "runs": [
            {"run": "ledger-R1", "path": "runs/ledger-R1", "contended": False},
            {"run": "ledger-R2", "path": "runs/ledger-R2", "contended": False},
            {"run": "ledger-R3", "path": "runs/ledger-R3", "contended": False},
        ],
        # The ledger runs reuse the frozen base and rank seed 7 but did not
        # retain the summary row; the pilot run on the same base did.
        "rank_summaries": ["runs/pilot-R1-concurrent-session/ic.stdout.txt"],
        "contention_source": None,
    },
    "n73": {
        "dir": "experiments/koblitz-single-target-n73-20261003",
        "n": 73,
        "subgroup_order": 86020738150056119,
        "runs": [
            {"run": "R1", "path": f"runs/{N73_RUN}R1", "contended": False},
            {"run": "R2", "path": f"runs/{N73_RUN}R2", "contended": True},
            {"run": "R3", "path": f"runs/{N73_RUN}R3", "contended": True},
        ],
        "rank_summaries": [
            f"runs/{N73_RUN}R1/ic_summary.jsonl",
            f"runs/{N73_RUN}R2/ic_summary.jsonl",
            f"runs/{N73_RUN}R3/ic_summary.jsonl",
        ],
        "contention_source": "README.md: R2/R3 ran while an unrelated verification process shared the host",
    },
    "n83": {
        "dir": "experiments/koblitz-single-target-n83-20261006",
        "n": 83,
        "subgroup_order": 8569786107849059,
        "runs": [
            {"run": "R1", "path": "runs/N83A1K600W12We202610061R21R1", "contended": True},
            {"run": "R2", "path": "runs/N83A1K600W12We202610062R22R2", "contended": True},
            {"run": "R3", "path": "runs/N83A1K600W12We202610063R23R3", "contended": True},
        ],
        "rank_summaries": [
            "runs/N83A1K600W12We202610061R21R1/ic_summary.jsonl",
            "runs/N83A1K600W12We202610062R22R2/ic_summary.jsonl",
            "runs/N83A1K600W12We202610063R23R3/ic_summary.jsonl",
        ],
        "contention_source": "README.md: unrelated verification processes shared the host during the runs",
    },
}


class AccountingError(RuntimeError):
    """Missing or inconsistent run data; the block is not written."""


def json_rows(path: Path) -> list[dict[str, Any]]:
    """Read a JSON object, a JSON-lines file, or JSON lines mixed with text."""
    text = path.read_text()
    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        value = None
    if isinstance(value, dict):
        return [value]
    if isinstance(value, list):
        return [row for row in value if isinstance(row, dict)]
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("{"):
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return rows


def one_row(path: Path, kind: str) -> dict[str, Any]:
    rows = [row for row in json_rows(path) if row.get("kind") == kind]
    if len(rows) != 1:
        raise AccountingError(f"{path}: expected one {kind} row, found {len(rows)}")
    return rows[0]


def optional_json(path: Path) -> dict[str, Any] | None:
    return json.loads(path.read_text()) if path.is_file() else None


def expected_rho_steps(subgroup_order: int, n: int) -> float:
    """Expected rho walk on the signed-Frobenius class space of size r/(2n)."""
    return math.sqrt(math.pi * subgroup_order / (4 * n))


def account_run(rung_dir: Path, spec: dict[str, Any]) -> dict[str, Any]:
    run_dir = rung_dir / spec["path"]
    ic = one_row(run_dir / "ic.jsonl", "compact_orbit_dlp_target")
    rho = one_row(run_dir / "rho.jsonl", "rho_public_fixture")
    probes = int(ic["probes"])
    ic_ms = float(ic["target_ms"])
    steps = int(rho["walk_steps"])
    walk_ms = float(rho["walk_ms"])
    rho_ms = float(rho["total_ms"])
    if min(probes, steps) <= 0 or min(ic_ms, walk_ms, rho_ms) <= 0:
        raise AccountingError(f"{run_dir}: non-positive probe, step or time field")
    run_json = optional_json(run_dir / "run.json") or {}
    ic_exec = optional_json(run_dir / "ic_execution.json") or {}
    rho_exec = optional_json(run_dir / "rho_execution.json") or {}
    ic_rate = probes / (ic_ms / 1000.0)
    rho_rate = steps / (walk_ms / 1000.0)
    return {
        "run": spec["run"],
        "host_contended": spec["contended"],
        "exploratory": spec["contended"],
        "ic_probes": probes,
        "ic_online_ms": ic_ms,
        "rho_walk_steps": steps,
        "rho_ideal_steps_recorded": rho.get("ideal_steps"),
        "rho_walk_ms": walk_ms,
        "rho_online_ms": rho_ms,
        "wall_speedup": rho_ms / ic_ms,
        "rho_steps_per_ic_probe_measured": steps / probes,
        "ic_probes_per_s": ic_rate,
        "rho_steps_per_s": rho_rate,
        "ic_to_rho_rate_ratio": ic_rate / rho_rate,
        "ic_rank_threads": run_json.get("ic_rank_threads"),
        "ic_target_threads": run_json.get("ic_target_threads"),
        "rho_threads": run_json.get("rho_threads"),
        "ic_peak_rss_bytes": ic_exec.get("peak_rss_bytes"),
        "rho_peak_rss_bytes": rho_exec.get("peak_rss_bytes"),
    }


def account_rung(name: str) -> dict[str, Any]:
    cfg = RUNGS[name]
    rung_dir = REPO / cfg["dir"]
    n = int(cfg["n"])
    r = int(cfg["subgroup_order"])
    expected = expected_rho_steps(r, n)
    runs = [account_run(rung_dir, spec) for spec in cfg["runs"]]

    for row in runs:
        recorded = row["rho_ideal_steps_recorded"]
        if recorded is not None and abs(recorded - expected) > expected * 1e-6:
            raise AccountingError(
                f"{name} {row['run']}: rho fixture ideal_steps {recorded} != sqrt(pi r/(4n)) {expected}"
            )
    target_probes = {row["ic_probes"] for row in runs}
    if len(target_probes) != 1:
        raise AccountingError(f"{name}: frozen-target probes differ across runs: {sorted(target_probes)}")
    frozen_probes = target_probes.pop()

    claim = json.loads((rung_dir / "claim_report_vs_rho.json").read_text())
    base_digest = (claim.get("fixture_hash") or {}).get("base_point_set_digest")
    rank_means, rank_attempts, rank_threads = [], set(), set()
    for rel in cfg["rank_summaries"]:
        summary = one_row(rung_dir / rel, "compact_orbit_dlp_summary")
        if base_digest and summary.get("base_hash") != base_digest:
            raise AccountingError(f"{name}: {rel} base_hash differs from the claim's base digest")
        rank_means.append(float(summary["rank_probes_mean"]))
        rank_attempts.add(int(summary["rank_attempts"]))
        rank_threads.add(summary.get("threads"))
    if len(rank_attempts) != 1:
        raise AccountingError(f"{name}: rank_attempts differ across summaries")
    attempts = rank_attempts.pop()
    rank_mean = statistics.fmean(rank_means)
    precompute_probes = rank_mean * attempts

    walls = [row["wall_speedup"] for row in runs]
    clean_walls = [row["wall_speedup"] for row in runs if not row["host_contended"]]
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_by": GENERATED_BY,
        "unit_assumption": UNIT_ASSUMPTION,
        "comparison_status": "native_counters_only",
        "operation_units": {"ic": "target probes", "rho": "walk steps"},
        "ic_online_operations": frozen_probes,
        "ic_online_operations_basis": "deterministic probes to the first relation for the frozen target (identical on every run)",
        "rho_online_operations": expected,
        "rho_online_operations_basis": "expected walk sqrt(pi*r/(4n)) on signed-Frobenius classes; measured walks are under runs",
        "rho_per_ic_native_counter": expected / frozen_probes,
        "ops_speedup_online": None,
        "total_work_S": None,
        "ic_rank_probes_mean_per_relation": rank_means,
        "ic_rank_probes_mean_pooled": rank_mean,
        "ic_rank_probes_source": cfg["rank_summaries"],
        "ic_rank_probes_note": (
            "guided rank queries and target extraction use different input laws and scan starts; "
            "this mean is a rank-stage diagnostic, not an estimate of unseen-target cost; "
            "per-relation counts are not retained, so no confidence interval is available"
        ),
        "rho_expected_steps_per_rank_mean_probe": expected / rank_mean,
        "rank_mean_to_frozen_probe_ratio": rank_mean / frozen_probes,
        "precompute_rank_probes": precompute_probes,
        "precompute_rank_probes_note": "rank_probes_mean x rank_attempts; the S3 root-index build is not counted in probes",
        "rho_expected_steps_per_rank_plus_target_probe": expected / (precompute_probes + rank_mean),
        "wall_speedup_median_all_runs": statistics.median(walls),
        "wall_speedup_median_uncontended_runs": statistics.median(clean_walls) if clean_walls else None,
        "rho_steps_per_ic_probe_measured_median": statistics.median(row["rho_steps_per_ic_probe_measured"] for row in runs),
        "ic_to_rho_rate_ratio_median": statistics.median(row["ic_to_rho_rate_ratio"] for row in runs),
        "ic_rank_threads": sorted(t for t in rank_threads if t is not None) or None,
        "thread_note": "ic_target_threads and rho_threads were not recorded by these run drivers (null)",
        "host_contention": {
            "contended_runs": [row["run"] for row in runs if row["host_contended"]],
            "source": cfg["contention_source"],
        },
        "runs": runs,
    }


def apply_to_claim(name: str, block: dict[str, Any]) -> Path:
    path = REPO / RUNGS[name]["dir"] / "claim_report_vs_rho.json"
    text = path.read_text()
    claim = json.loads(text)
    # Keep the file's own layout so the diff shows only the accounting change.
    second = text.splitlines()[1] if text.count("\n") > 1 else "  "
    indent = len(second) - len(second.lstrip()) or 2
    sort_keys = list(claim) == sorted(claim)
    claim["operation_accounting"] = block
    non_claims = list(claim.get("claim_boundary_non_claims") or [])
    wanted = list(NON_CLAIMS)
    if block["host_contention"]["contended_runs"]:
        wanted.append(CONTENTION_NON_CLAIM)
    for item in wanted:
        if item not in non_claims:
            non_claims.append(item)
    claim["claim_boundary_non_claims"] = non_claims
    path.write_text(json.dumps(claim, indent=indent, sort_keys=sort_keys) + "\n")
    return path


def table(blocks: dict[str, dict[str, Any]]) -> str:
    header = (
        "| rung | IC probes, frozen | IC probes, rank mean | rank/frozen probes | rho expected steps | "
        "rho steps/IC probe | rho steps/rank probe | rho steps/(rank+target probe) | wall x median | wall x uncontended | native rate ratio |"
    )
    lines = [header, "|" + "---|" * 11]
    for name, b in blocks.items():
        clean = b["wall_speedup_median_uncontended_runs"]
        lines.append(
            f"| {name} | {b['ic_online_operations']:,} | {b['ic_rank_probes_mean_pooled']:,.0f} | "
            f"{b['rank_mean_to_frozen_probe_ratio']:.2f} | {b['rho_online_operations']:,.0f} | "
            f"{b['rho_per_ic_native_counter']:.3g} | {b['rho_expected_steps_per_rank_mean_probe']:.3g} | "
            f"{b['rho_expected_steps_per_rank_plus_target_probe']:.3g} | {b['wall_speedup_median_all_runs']:.4g} | "
            f"{'-' if clean is None else format(clean, '.4g')} | {b['ic_to_rho_rate_ratio_median']:.3g} |"
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rung", choices=sorted(RUNGS), action="append")
    parser.add_argument("--write", action="store_true", help="update claim_report_vs_rho.json")
    args = parser.parse_args()
    names = args.rung or list(RUNGS)
    blocks = {name: account_rung(name) for name in names}
    print(table(blocks))
    if args.write:
        for name, block in blocks.items():
            print("wrote", apply_to_claim(name, block).relative_to(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
