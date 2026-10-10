#!/usr/bin/env python3
"""Compose the immutable Stage 175 result from raw process-meter artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics
import subprocess


STAGE = Path(__file__).resolve().parent
REPO = next(parent for parent in STAGE.parents if (parent / "Cargo.toml").is_file())
CAMPAIGN = STAGE / "development" / "campaign"
BUILD = STAGE / "development" / "build"
EXPECTED_COMMIT = "c38626a5e2674dbfdff796474ea83af0c863917a"
EXPECTED_INSTANCE = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7"
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
STAGE174_MEDIAN = {
    "wall_seconds": 26.358217916989815,
    "total_core_seconds": 147.841771,
    "peak_rss_bytes": 2_756_362_240,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict:
    return {
        "path": str(path.relative_to(STAGE)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def load_run(backend: str, repeat: int) -> dict:
    root = CAMPAIGN / f"{backend}-r{repeat}"
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    metrics = json.loads(metrics_path.read_text())
    report = json.loads(stdout_path.read_text())
    return {
        "backend": backend,
        "repeat": repeat,
        "process": metrics,
        "report": report,
        "artifacts": {
            "metrics": artifact(metrics_path),
            "stdout": artifact(stdout_path),
            "stderr": artifact(stderr_path),
        },
    }


def medians(runs: list[dict]) -> dict:
    return {
        field: statistics.median(run["process"]["metrics"][field] for run in runs)
        for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }


def ratio(numerator: dict, denominator: dict) -> dict:
    return {
        field: numerator[field] / denominator[field]
        for field in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }


def main() -> int:
    provenance_path = CAMPAIGN / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    if provenance["source_commit"] != EXPECTED_COMMIT:
        raise SystemExit("campaign source commit mismatch")
    current = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    if current != EXPECTED_COMMIT:
        raise SystemExit("compose from the exact measured commit")

    f4 = [load_run("native-f4", repeat) for repeat in range(1, 4)]
    direct = [load_run("direct-mitm", repeat) for repeat in range(1, 4)]
    checks: list[dict] = []
    structural = []
    for run in f4:
        process = run["process"]
        report = run["report"]
        contract = report["factor_base_contract"]
        ok = (
            process["returncode"] == 0
            and not process["timed_out"]
            and report["status"] == "unsat"
            and report["source_instance_id"] == EXPECTED_INSTANCE
            and report["source_instance_verified"] is True
            and report["regenerated_source_exact"] is True
            and contract["target_subgroup_enumerated"] is False
            and contract["discrete_log_labels_used"] is False
            and report["fixed_x1_masks_visited"] == 512
            and report["fixed_x1_systems_constructed"] == 242
            and report["fixed_x1_systems_completed"] == 242
            and report["exhaustive"] is True
            and report["algebraic_roots"] == 0
            and report["witness_points"] is None
            and report["conflicts"] is None
            and report["solver_equations_blake3"] == EXPECTED_EQUATIONS
        )
        checks.append({"name": f"native-f4-r{run['repeat']}", "pass": ok})
        extra = report["cost"]["extra"]
        structural.append(
            {
                "ops": report["cost"]["ops"],
                "word_xors_performed": extra["word_xors_performed"],
                "steps": extra["steps"],
                "pairs_reduced": extra["pairs_reduced"],
                "field_pairs_reduced": extra["field_pairs_reduced"],
                "reducer_rows": extra["reducer_rows"],
                "matrix_rows_sum": extra["matrix_rows_sum"],
                "new_elements": extra["new_elements"],
            }
        )
    checks.append({"name": "native-f4-structural-repeatability", "pass": len({json.dumps(x, sort_keys=True) for x in structural}) == 1})
    for run in direct:
        process = run["process"]
        report = run["report"]
        checks.append(
            {
                "name": f"direct-mitm-r{run['repeat']}",
                "pass": process["returncode"] == 0
                and not process["timed_out"]
                and report["status"] == "unsat"
                and report["source_instance_id"] == EXPECTED_INSTANCE
                and report["source_instance_verified"] is True
                and report["exhaustive"] is True,
            }
        )
    if not all(check["pass"] for check in checks):
        raise SystemExit("one or more correctness gates failed")

    f4_median = medians(f4)
    direct_median = medians(direct)
    vs_stage174 = ratio(f4_median, STAGE174_MEDIAN)
    vs_direct = ratio(f4_median, direct_median)
    improvement_pass = (
        vs_stage174["wall_seconds"] < 1 and vs_stage174["total_core_seconds"] < 1
    )
    build_metrics_path = BUILD / "metrics.json"
    build_metrics = json.loads(build_metrics_path.read_text())
    all_processes = [build_metrics, *(run["process"] for run in f4), *(run["process"] for run in direct)]
    campaign_charge = {
        "wall_seconds_sum": sum(p["metrics"]["wall_seconds"] for p in all_processes),
        "total_core_seconds_sum": sum(p["metrics"]["total_core_seconds"] for p in all_processes),
        "peak_rss_bytes_max": max(p["metrics"]["peak_rss_bytes"] for p in all_processes),
        "components": len(all_processes),
    }

    result = {
        "schema": "koblitz_stage175_current_f4_single_target.v1",
        "claim_boundary": "One already-opened public n=59 true-negative decomposition target; implementation engineering only, not relation yield, a completed index-calculus DLP, a rho crossover, independent reproduction, novelty review, or SOTA.",
        "source": {
            "commit": EXPECTED_COMMIT,
            "parents": subprocess.check_output(
                ["git", "show", "-s", "--format=%P", EXPECTED_COMMIT], cwd=REPO, text=True
            ).strip().split(),
            "binary_sha256": provenance["binary_sha256"],
            "engine": "current repository Boolean F4 with adaptive leading-block BlockTables",
        },
        "input": {
            "n": 59,
            "ell": 9,
            "m": 3,
            "blind_instance_id": "b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4",
            "source_instance_id": EXPECTED_INSTANCE,
            "manifest": artifact(STAGE / "input" / "manifest.json"),
            "equation_fingerprint": EXPECTED_EQUATIONS,
            "factor_base_algebraic": True,
            "target_subgroup_enumerated": False,
            "discrete_log_labels_used": False,
        },
        "build": {
            "process": build_metrics,
            "artifacts": {
                "metrics": artifact(build_metrics_path),
                "stdout": artifact(BUILD / "stdout.txt"),
                "stderr": artifact(BUILD / "stderr.txt"),
            },
        },
        "candidate": {
            "configuration": provenance["environment_controls"],
            "runs": f4,
            "median": f4_median,
            "structural_counts": structural[0],
            "conflicts": None,
            "single_core_seconds": None,
        },
        "direct_mitm": {"runs": direct, "median": direct_median},
        "stage174_reference_median": STAGE174_MEDIAN,
        "ratios": {
            "candidate_over_stage174": vs_stage174,
            "candidate_over_same_binary_direct_mitm": vs_direct,
        },
        "correctness_checks": checks,
        "decision": {
            "engineering_improvement_pass": improvement_pass,
            "status": "REJECTED_REGRESSION" if not improvement_pass else "ACCEPTED_ENGINEERING_IMPROVEMENT",
            "reason": "Both wall and total core-seconds had to beat Stage 174; the current-main engine regressed both." if not improvement_pass else "Both preregistered timing conditions passed.",
        },
        "campaign_charge": campaign_charge,
        "provenance": artifact(provenance_path),
        "validation": {
            "backend_tests": "3 passed",
            "current_f4_tests": "8 passed",
            "fixed_x1_specialisation_test": "1 passed",
        },
    }
    (STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    def f(value: float) -> str:
        return f"{value:.6f}"

    markdown = f"""# Stage 175: current repository F4 on one frozen target

The current repository `BlockTables` F4 engine was run three times on the exact
opened `n=59, ell=9, m=3` true-negative used by Stage 174. All repeats
authenticated the source, visited all 512 fixed-X1 masks, completed all 242
rational systems, reproduced equation fingerprint `{EXPECTED_EQUATIONS}`, and
returned exhaustive UNSAT. The algebraic factor-base contract still records no
target-subgroup enumeration and no discrete-log labels.

| backend | wall median (s) | core median (s) | peak RSS median (bytes) | status |
|---|---:|---:|---:|---|
| Stage 174 selected F4 | {f(STAGE174_MEDIAN['wall_seconds'])} | {f(STAGE174_MEDIAN['total_core_seconds'])} | {STAGE174_MEDIAN['peak_rss_bytes']} | historical exhaustive UNSAT |
| current repository F4 | {f(f4_median['wall_seconds'])} | {f(f4_median['total_core_seconds'])} | {int(f4_median['peak_rss_bytes'])} | exhaustive UNSAT |
| same-binary direct MITM | {f(direct_median['wall_seconds'])} | {f(direct_median['total_core_seconds'])} | {int(direct_median['peak_rss_bytes'])} | exhaustive UNSAT |

The current F4 is `{vs_stage174['wall_seconds']:.3f}x` the Stage 174 wall,
`{vs_stage174['total_core_seconds']:.3f}x` its CPU, and
`{vs_stage174['peak_rss_bytes']:.3f}x` its RSS. It is therefore rejected by the
preregistered rule. Against same-binary direct MITM it is
`{vs_direct['wall_seconds']:.2f}x` wall, `{vs_direct['total_core_seconds']:.2f}x`
CPU, and `{vs_direct['peak_rss_bytes']:.2f}x` RSS.

The repeated current-F4 structural counts are exact: 242 calls, 1,204 steps,
`319313687585` row-equivalent word XORs, `147794583858` actually performed
table-assisted word XORs, and no conflicts (not applicable to F4). The clean
build cost {f(build_metrics['metrics']['wall_seconds'])} wall seconds,
{f(build_metrics['metrics']['total_core_seconds'])} core-seconds, and
{build_metrics['metrics']['peak_rss_bytes']} bytes peak RSS. The complete
seven-process Stage 175 charge is {f(campaign_charge['wall_seconds_sum'])} wall
seconds, {f(campaign_charge['total_core_seconds_sum'])} core-seconds, and
{campaign_charge['peak_rss_bytes_max']} bytes maximum RSS.

This is a negative single-target engineering result. It does not change any
SOTA gate and does not measure natural relation yield or a full unknown-scalar
index-calculus run.
"""
    (STAGE / "RESULTS.md").write_text(markdown)
    print(json.dumps({"decision": result["decision"], "ratios": result["ratios"], "campaign_charge": campaign_charge}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
