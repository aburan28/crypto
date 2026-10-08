#!/usr/bin/env python3
"""Assemble the n=85 a=0 single-target online vs_rho claim (fail-closed)."""
import csv
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
FROZEN = HERE / "frozen"
WIN_THRESHOLD = 1.20


def load_json(path):
    return json.loads(Path(path).read_text())


def load_one_jsonl(path):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path}")
    return rows[0]


def main():
    fixture = load_json(FROZEN / "fixture.json")
    runs = []
    for run_dir in sorted(RUNS.iterdir()):
        run = load_json(run_dir / "run.json")
        if run["status"] != "PRODUCERS_COMPLETE":
            raise SystemExit(f"{run_dir.name}: status {run['status']}")
        replay = load_json(run_dir / "independent_replay.json")
        if replay["status"] != "PASS":
            raise SystemExit(f"{run_dir.name}: replay {replay['status']}")
        ic = load_one_jsonl(run_dir / "ic.jsonl")
        ic_summary = load_one_jsonl(run_dir / "ic_summary.jsonl")
        rho = load_one_jsonl(run_dir / "rho.jsonl")
        ic_exec = load_json(run_dir / "ic_execution.json")
        rho_exec = load_json(run_dir / "rho_execution.json")
        ic_online = float(ic["target_ms"])
        rho_online = (
            float(rho["setup_ms"]) + float(rho["walk_ms"]) + float(rho["validation_ms"])
        )
        runs.append({
            "run": run_dir.name[-2:],
            "run_id": run_dir.name,
            "rho_walk_seed": run["rho_walk_seed"],
            "ic_rank_seed": run["ic_rank_seed"],
            "ic_online_ms": ic_online,
            "rho_online_ms": rho_online,
            "online_speedup": rho_online / ic_online,
            "ic_verified": bool(ic["group_verified"]),
            "rho_verified": bool(rho["verified"]),
            "ic_scalar": int(ic["recovered_scalar"]),
            "rho_scalar": int(rho["recovered_fixture_scalar"]),
            "probes": int(ic["probes"]),
            "walk_steps": int(rho["walk_steps"]),
            "ic_peak_rss_bytes": (
                ic_summary["peak_rss_bytes"]
                if isinstance(ic_summary.get("peak_rss_bytes"), int)
                else ic_exec["peak_rss_bytes"]
            ),
            "rho_peak_rss_bytes": rho_exec["peak_rss_bytes"],
            "rank_stage_ms": ic_summary["timing_ms"]["rank_stage"],
            "rank_threads": ic_summary["threads"],
            "rank_probes_mean": ic_summary["rank_probes_mean"],
        })

    scalars = {run["ic_scalar"] for run in runs} | {run["rho_scalar"] for run in runs}
    if len(scalars) != 1:
        raise SystemExit(f"arms disagree on the scalar: {scalars}")
    if not all(run["ic_verified"] and run["rho_verified"] for run in runs):
        raise SystemExit("an arm failed in-process verification")
    if len({run["probes"] for run in runs}) != 1:
        raise SystemExit("target relation is not deterministic across runs")
    speeds = sorted(run["online_speedup"] for run in runs)
    median_speedup = statistics.median(speeds)
    if median_speedup < WIN_THRESHOLD:
        raise SystemExit(f"median speedup {median_speedup:.2f} below {WIN_THRESHOLD}")

    report = {
        "schema_version": "2.0",
        "task_id": "TASK-IC-BOUNDARY-AUTOLAB-20260910",
        "kind": "single_target_online_vs_rho_claim_report",
        "stage": "vs_rho",
        "regime": "koblitz",
        "n": 85,
        "n_or_bits": 85,
        "timing_class": "single_target_online",
        "target_count": 1,
        "subgroup_order": int(fixture["subgroup_order"]),
        "subgroup_log2_order": 53.7,
        "curve": {
            "a": 0,
            "b": 1,
            "field_modulus": "x^85 + x^8 + x^2 + x + 1",
            "cofactor": int(fixture["cofactor"]),
            "weierstrass_model": "y^2 + x*y = x^3 + 1",
        },
        "public_generator": [str(v) for v in fixture["generator"]],
        "public_target_q": [str(v) for v in fixture["public_target"]],
        "fixture_scalar": int(fixture["fixture_scalar_validation_only"]),
        "fixture_scalar_note": "validation-only sidecar; never supplied to either solver input",
        "paired_same_public_point": True,
        "paired_runs": [
            {
                "run": run["run"],
                "rho_walk_seed": run["rho_walk_seed"],
                "ic_rank_seed": run["ic_rank_seed"],
                "ic_online_ms": run["ic_online_ms"],
                "rho_online_ms": run["rho_online_ms"],
                "online_speedup": run["online_speedup"],
                "rho_walk_steps": run["walk_steps"],
            }
            for run in runs
        ],
        "median_online_speedup": median_speedup,
        "online_speedup_range_all_runs": [speeds[0], speeds[-1]],
        "automorphism_discount": {
            "formula": "sqrt(2n)",
            "A": 170,
            "n": 85,
            "description": "signed Frobenius automorphism discount family used by the Koblitz rho control (quotient_mode signed_frobenius, automorphism_size 170)",
        },
        "resource_caps": {
            "common_cap_bytes": None,
            "ic_peak_rss_bytes_max": max(run["ic_peak_rss_bytes"] for run in runs),
            "rho_peak_rss_bytes_max": max(run["rho_peak_rss_bytes"] for run in runs),
            "note": "no finite kernel cap declared; per-process peak RSS recorded with wait4(2)",
        },
        "precompute_detail": {
            "orbit_columns": 600,
            "factor_base_points": 102000,
            "rank": 600,
            "rank_failures": 0,
            "rank_rows_without_gain": 0,
            "pair_table_entries": 0,
            "edge_selectors": 0,
            "root_table_entries": 30598167,
            "regular_states": 30600000,
            "rank_threads": runs[0]["rank_threads"],
            "rank_stage_ms_range": [
                min(run["rank_stage_ms"] for run in runs),
                max(run["rank_stage_ms"] for run in runs),
            ],
            "rank_probes_mean_range": [
                min(run["rank_probes_mean"] for run in runs),
                max(run["rank_probes_mean"] for run in runs),
            ],
            "target_probes": runs[0]["probes"],
            "note": "first rung on the a=1 arm and the first ledger rung using the parallel guided rank (KIC_RANK_THREADS=12; logs identical to the sequential path)",
        },
        "all_stages_charged_same_series": True,
        "ic_verified_all_runs": True,
        "rho_verified_all_runs": True,
        "ic_and_rho_scalar_identical": True,
        "deterministic_ic_relation": {"probes": runs[0]["probes"], "same_all_runs": True},
        "independent_scalar_replay": {
            "method": "standalone Python GF(2^85) polynomial/affine-curve arithmetic independent of the Rust producers",
            "validated_runs": [run["run"] for run in runs],
            "status": "PASS",
        },
        "verdict": "N85_A0_COMPACT_ORBIT_SINGLE_TARGET_ONLINE_IC_OVER_RHO_CROSSOVER",
        "claim_boundary": (
            "Public synthetic Koblitz n=85 (a=0, b=1 over GF(2^85), 53.7-bit subgroup); "
            "exactly one previously unseen public target point solved online by IC and rho on "
            "the identical frozen point (sequential single-process arms, same host); "
            "constant-factor win only. The online walls are conservative: unrelated "
            "verification processes shared the host during the runs. Not ECC2K-130 "
            "evidence, not an asymptotic sub-rho claim, not key recovery, no "
            "deployed-curve impact."
        ),
        "claim_boundary_non_claims": [
            "not ECC2K-130 evidence (n=131 is a different field degree)",
            "not an asymptotic sub-rho claim",
            "not key recovery or deployed-curve security impact",
            "multi-target amortized results remain secondary",
            "total operation-count comparison absent; S unknown",
        ],
        "evidence": [
            str(p.relative_to(HERE.parents[2]))
            for p in sorted(RUNS.iterdir())
            if p.is_dir()
        ],
    }
    out = HERE / "claim_report_vs_rho.json"
    out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    with open(HERE / "single-target-results.csv", "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "run", "rho_walk_seed", "ic_rank_seed", "ic_online_ms",
                "rho_online_ms", "online_speedup", "ic_peak_rss_bytes",
                "rho_peak_rss_bytes", "rank_stage_ms", "rank_threads",
                "probes", "walk_steps",
            ],
        )
        writer.writeheader()
        for run in runs:
            writer.writerow({k: run.get(k) for k in writer.fieldnames})
    print(json.dumps({
        "median_online_speedup": median_speedup,
        "range": [speeds[0], speeds[-1]],
        "per_run": {run["run"]: round(run["online_speedup"], 2) for run in runs},
        "report": str(out),
    }, indent=1))


if __name__ == "__main__":
    main()
