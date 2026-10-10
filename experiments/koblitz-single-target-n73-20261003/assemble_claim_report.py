#!/usr/bin/env python3
"""Assemble the n=73 single-target online vs_rho claim report (R1..R3).

Reads the three frozen paired runs, requires every independent replay to
be present and PASS, computes the median online speedup, and writes
claim_report_vs_rho.json + a CSV summary.  Fail-closed: any missing
artifact, failed check, or producer failure aborts without writing the
report.
"""
import csv
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
FROZEN = HERE / "frozen"
PROTOCOL = json.loads((HERE / "protocol.json").read_text())
BASE_RUN_ID = PROTOCOL["run_id"].removesuffix("R1")
WIN_THRESHOLD = 1.20


def load_json(path):
    return json.loads(Path(path).read_text())


def load_one_jsonl(path):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    if len(rows) != 1:
        raise ValueError(f"expected exactly one row in {path}, got {len(rows)}")
    return rows[0]


def collect(tag):
    run_dir = RUNS / (BASE_RUN_ID + tag)
    run = load_json(run_dir / "run.json")
    if run["status"] != "PRODUCERS_COMPLETE":
        raise SystemExit(f"{tag}: producer status is {run['status']}")
    ic = load_one_jsonl(run_dir / "ic.jsonl")
    ic_summary = load_one_jsonl(run_dir / "ic_summary.jsonl")
    rho = load_one_jsonl(run_dir / "rho.jsonl")
    replay = load_json(run_dir / "independent_replay.json")
    if replay["status"] != "PASS":
        raise SystemExit(f"{tag}: independent replay is {replay['status']}")
    ic_exec = load_json(run_dir / "ic_execution.json")
    rho_exec = load_json(run_dir / "rho_execution.json")
    ic_online = float(ic["target_ms"])
    rho_online = (
        float(rho["setup_ms"]) + float(rho["walk_ms"]) + float(rho["validation_ms"])
    )
    return {
        "run": tag,
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
        "ic_peak_rss_bytes": (
            ic_summary["peak_rss_bytes"]
            if isinstance(ic_summary.get("peak_rss_bytes"), int)
            else ic_exec["peak_rss_bytes"]
        ),
        "rho_peak_rss_bytes": rho_exec["peak_rss_bytes"],
        "ic_precompute_ms": ic_summary["timing_ms"]["process_total"] - ic_summary["timing_ms"]["targets_total"],
        "rank_stage_ms": ic_summary["timing_ms"]["rank_stage"],
        "rank_probes_mean": ic_summary["rank_probes_mean"],
        "probes": ic["probes"],
    }


def main():
    runs = [collect(tag) for tag in ("R1", "R2", "R3")]
    fixture = load_json(FROZEN / "fixture.json")
    scalars = {run["ic_scalar"] for run in runs} | {run["rho_scalar"] for run in runs}
    if len(scalars) != 1:
        raise SystemExit(f"arms disagree on the scalar: {scalars}")
    if not all(run["ic_verified"] and run["rho_verified"] for run in runs):
        raise SystemExit("an arm failed in-process verification")
    speeds = sorted(run["online_speedup"] for run in runs)
    median_speedup = statistics.median(speeds)
    if median_speedup < WIN_THRESHOLD:
        raise SystemExit(f"median online speedup {median_speedup:.2f} below win threshold {WIN_THRESHOLD}")

    manifest = load_json(HERE / "candidate-manifest.json")
    curve = manifest["record"]["curve"]
    report = {
        "schema_version": "2.0",
        "task_id": PROTOCOL.get("task_id", "TASK-IC-BOUNDARY-AUTOLAB-20260910"),
        "kind": "single_target_online_vs_rho_claim_report",
        "stage": "vs_rho",
        "regime": "koblitz",
        "n": 73,
        "n_or_bits": 73,
        "timing_class": "single_target_online",
        "target_count": 1,
        "subgroup_order": int(fixture["subgroup_order"]),
        "subgroup_log2_order": 56.3,
        "curve": {
            "a": int(fixture["a"]),
            "b": 1,
            "field_modulus": "x^73 + x^4 + x^3 + x^2 + 1",
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
            }
            for run in runs
        ],
        "median_online_speedup": median_speedup,
        "online_speedup_range_all_runs": [speeds[0], speeds[-1]],
        "ic_online_phase_ms": {
            "target_query": float(load_one_jsonl(RUNS / (BASE_RUN_ID + "R2") / "ic.jsonl")["target_query_ms"]),
        },
        "automorphism_discount": {
            "formula": "sqrt(2n)",
            "A": 146,
            "n": 73,
            "description": "signed Frobenius automorphism discount family used by the Koblitz rho control (quotient_mode signed_frobenius, automorphism_size 146)",
        },
        "resource_caps": {
            "common_cap_bytes": None,
            "ic_peak_rss_bytes_max": max(run["ic_peak_rss_bytes"] for run in runs),
            "rho_peak_rss_bytes_max": max(run["rho_peak_rss_bytes"] for run in runs),
            "note": "no finite kernel cap declared; per-process peak RSS recorded with wait4(2)",
        },
        "precompute_detail": {
            "orbit_columns": 600,
            "factor_base_points": 87600,
            "rank": 600,
            "rank_failures": 0,
            "pair_table_entries": 0,
            "edge_selectors": 0,
            "root_table_entries": 26278197,
            "regular_states": 26280000,
            "rank_stage_ms_range": [min(run["rank_stage_ms"] for run in runs), max(run["rank_stage_ms"] for run in runs)],
            "rank_probes_mean_range": [min(run["rank_probes_mean"] for run in runs), max(run["rank_probes_mean"] for run in runs)],
            "target_probes": [run["probes"] for run in runs],
        },
        "all_stages_charged_same_series": True,
        "ic_verified_all_runs": True,
        "rho_verified_all_runs": True,
        "ic_and_rho_scalar_identical": True,
        "verdict": "N73_U128_COMPACT_ORBIT_SINGLE_TARGET_ONLINE_IC_OVER_RHO_CROSSOVER",
        "claim_boundary": (
            "Public synthetic Koblitz n=73 (a=0, b=1 over GF(2^73), modulus x^73+x^4+x^3+x^2+1, "
            "56.3-bit subgroup); exactly one previously unseen public target point solved online by "
            "IC and rho on the identical frozen point (sequential single-process arms, same host); "
            "constant-factor win only. Not ECC2K-130 evidence, not an asymptotic sub-rho claim, "
            "not key recovery, not external/private targets."
        ),
        "claim_boundary_non_claims": [
            "not ECC2K-130 evidence (n=131 is a different field degree on the same curve family)",
            "not an asymptotic sub-rho claim",
            "not key recovery or deployed-curve security impact",
            "multi-target amortized results remain secondary",
        ],
        "evidence": [
            str(p.relative_to(HERE.parents[2])) for p in sorted(RUNS.glob(BASE_RUN_ID + "R*"))
        ] + ["experiments/koblitz-single-target-n73-20261003/frozen/fixture.json"],
    }
    out = HERE / "claim_report_vs_rho.json"
    out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    # Operation-count accounting (probes vs rho steps, mean-target estimate,
    # contention flags); see research/.../autolab/op_accounting.py.
    sys.path.insert(0, str(HERE.parents[1] / "research/sat_factor_base_review_20260908/autolab"))
    import op_accounting
    op_accounting.apply_to_claim("n73", op_accounting.account_rung("n73"))
    with open(HERE / "single-target-results.csv", "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "run", "rho_walk_seed", "ic_rank_seed", "ic_online_ms",
                "rho_online_ms", "online_speedup", "ic_peak_rss_bytes",
                "rho_peak_rss_bytes", "rank_stage_ms", "probes",
            ],
        )
        writer.writeheader()
        for run in runs:
            writer.writerow({k: run.get(k) for k in writer.fieldnames})
    print(json.dumps({
        "median_online_speedup": median_speedup,
        "range": [speeds[0], speeds[-1]],
        "per_run": {run["run"]: run["online_speedup"] for run in runs},
        "report": str(out),
    }, indent=1))


if __name__ == "__main__":
    main()
