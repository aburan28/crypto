"""Describe the frozen N41 C800 and N53 C1200 full-rank pilot gains.

These targets were frozen in separate experiments. The independent-curve
contrast describes these complete pipelines, not a causal field-size effect.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


ROOT = Path(__file__).resolve().parents[3]
INPUTS = {
    "N41": {
        "analysis": "experiments/koblitz-s3-n41-c800-fullrank-20261009/local/analysis.json",
        "ledger": "experiments/koblitz-s3-n41-c800-fullrank-20261009/local/raw_holdout96.jsonl.gz",
        "replay": "experiments/koblitz-s3-n41-c800-fullrank-20261009/local/sage_replay.json",
        "freeze": "experiments/koblitz-s3-n41-c800-fullrank-20261009/FREEZE.json",
    },
    "N53": {
        "analysis": "experiments/koblitz-s3-target-conjugate-20261009/local/analysis_c1200_fullrank.json",
        "ledger": "experiments/koblitz-s3-target-conjugate-20261009/local/raw_c1200_fullrank_paired32.jsonl.gz",
        "replay": "experiments/koblitz-s3-target-conjugate-20261009/local/sage_replay_c1200_fullrank.json",
        "freeze": "experiments/koblitz-s3-target-conjugate-20261009/C1200_FULLRANK_FREEZE.json",
    },
}
ARMS = ("baseline", "expanded")
BOOTSTRAP_SAMPLES = 20_000
Z_SUM_SQUARED = (1.96 + 0.84) ** 2


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def percentile_95(values: list[float]) -> list[float]:
    ordered = sorted(values)
    assert len(ordered) == BOOTSTRAP_SAMPLES
    return [ordered[499], ordered[19_499]]


def load_curve(name: str) -> dict:
    paths = {key: ROOT / path for key, path in INPUTS[name].items()}
    analysis_bytes = paths["analysis"].read_bytes()
    analysis = json.loads(analysis_bytes)
    ledger_bytes = paths["ledger"].read_bytes()
    raw = gzip.decompress(ledger_bytes)
    replay_bytes = paths["replay"].read_bytes()
    replay = json.loads(replay_bytes)
    freeze_bytes = paths["freeze"].read_bytes()
    assert analysis["all_native_verified"] and not analysis["failure_rows"]
    assert replay["all_scalars_independently_recovered"]
    assert analysis["controlled_speedup"] is None
    assert analysis["raw_ledger_sha256"] == replay["ledger_sha256"] == digest(raw)
    assert analysis["freeze_sha256"] == replay["freeze_sha256"] == digest(freeze_bytes)
    assert analysis["run_count"] == replay["run_count"]
    blocks = []
    for label, block in sorted(analysis["target_blocks"].items()):
        logs = {}
        for arm in ARMS:
            times = block[arm]["online_ms"]
            assert len(times) == 2 and all(value > 0 for value in times)
            logs[arm] = [math.log(value) for value in times]
        blocks.append({"label": label, "logs": logs})
    assert len(blocks) * 4 == analysis["run_count"]
    effects = [statistics.mean(block["logs"]["baseline"]) -
               statistics.mean(block["logs"]["expanded"]) for block in blocks]
    mean_effect = statistics.mean(effects)
    assert math.isclose(math.exp(mean_effect),
                        analysis["exploratory_baseline_over_expanded_online_geomean"],
                        rel_tol=1e-12)
    within_sum = sum(statistics.mean(statistics.variance(block["logs"][arm])
                                     for block in blocks) for arm in ARMS)
    observed_var = statistics.variance(effects)
    aa = {}
    for arm in ARMS:
        values = sorted(abs(block["logs"][arm][0] - block["logs"][arm][1])
                        for block in blocks)
        aa[arm] = values[math.ceil(0.95 * len(values)) - 1]
        assert math.isclose(aa[arm],
                            analysis["arms"][arm]["aa_p95_abs_log_online_ratio"],
                            rel_tol=1e-12)
    return {
        "input_sha256": {key: digest(paths[key].read_bytes()) for key in paths},
        "target_count": len(blocks),
        "run_count": analysis["run_count"],
        "candidate_ids": {arm: analysis["arms"][arm]["candidate_id"] for arm in ARMS},
        "geometric_online_gain": math.exp(mean_effect),
        "mean_log_gain": mean_effect,
        "observed_target_log_variance": observed_var,
        "mean_within_target_log_variance_sum": within_sum,
        "estimated_between_target_log_variance": max(0, observed_var - within_sum / 2),
        "aa_p95_abs_log_online_ratio": aa,
        "_blocks": blocks,
    }


def sampled_mean_log_gain(rng: random.Random, blocks: list[dict]) -> float:
    chosen = rng.choices(blocks, k=len(blocks))
    effects = []
    for block in chosen:
        logs = block["logs"]
        effects.append(statistics.mean(rng.choices(logs["baseline"], k=2)) -
                       statistics.mean(rng.choices(logs["expanded"], k=2)))
    return statistics.mean(effects)


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    curves = {name: load_curve(name) for name in ("N41", "N53")}
    rng = random.Random(20261009)
    samples = {name: [] for name in curves}
    contrasts = []
    for _ in range(BOOTSTRAP_SAMPLES):
        draws = {name: sampled_mean_log_gain(rng, curve["_blocks"])
                 for name, curve in curves.items()}
        for name, value in draws.items():
            samples[name].append(math.exp(value))
        contrasts.append(math.exp(draws["N53"] - draws["N41"]))
    for name, curve in curves.items():
        curve.pop("_blocks")
        curve["hierarchical_bootstrap_95"] = percentile_95(samples[name])
        curve["aa_noise_gate_pass"] = all(value <= math.log(1.10)
                                          for value in curve["aa_p95_abs_log_online_ratio"].values())
    delta = curves["N53"]["mean_log_gain"] - curves["N41"]["mean_log_gain"]
    plan = {}
    for repeats in (2, 4, 8):
        variance_sum = sum(curve["estimated_between_target_log_variance"] +
                           curve["mean_within_target_log_variance_sum"] / repeats
                           for curve in curves.values())
        required = math.ceil(Z_SUM_SQUARED * variance_sum / delta**2) if delta else None
        plan[str(repeats)] = {
            "equal_independent_targets_per_curve": required,
            "one_target_solves_both_curves": 4 * repeats * required if required else None,
        }
    result = {
        "kind": "s3_fullrank_crosscurve_descriptive_pilot_v1",
        "curves": curves,
        "N53_over_N41_geometric_gain_ratio": math.exp(delta),
        "N53_over_N41_gain_ratio_hierarchical_bootstrap_95": percentile_95(contrasts),
        "independent_target_equal_allocation_plan": plan,
        "power_method": "ceil((1.96+0.84)^2 * (v41(k)+v53(k)) / delta^2), with v(k)=estimated_between_target_variance+mean_within_target_variance_sum/k",
        "power_assumptions": "two-sided alpha 0.05, power 0.80, observed positive log contrast persists, independent targets, k repeats per arm, stable host and unchanged target laws",
        "bootstrap_method": "20,000 deterministic independent target-block draws per curve, followed by two within-arm timing draws with replacement",
        "controlled_speedup": None,
        "curve_causal_effect": None,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.output)
    print(json.dumps({key: summary[key] for key in
                      ("N53_over_N41_geometric_gain_ratio",
                       "N53_over_N41_gain_ratio_hierarchical_bootstrap_95",
                       "independent_target_equal_allocation_plan")}, sort_keys=True))
