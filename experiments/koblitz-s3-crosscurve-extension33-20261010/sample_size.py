"""Estimate distinct matched-target counts for curve-gain differences."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "local/analysis.json"
BLOCKS = HERE / "local/prior_analysis.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def estimate(sd: float, log_effect: float, z_total: float) -> int:
    return math.ceil((z_total * sd / abs(log_effect)) ** 2)


def make() -> dict:
    analysis_raw, blocks_raw = ANALYSIS.read_bytes(), BLOCKS.read_bytes()
    analysis, prior = json.loads(analysis_raw), json.loads(blocks_raw)
    assert analysis["complete"] and analysis["matched_index_count"] == 33
    assert analysis["prior_analysis_sha256"] == sha(blocks_raw)
    assert prior["matched_index_count"] == 33 and prior["run_count"] == 264
    blocks = {(item["n"], item["label"]): item for item in prior["blocks"]}
    assert len(blocks) == 66 and all(item["all_verified"] for item in blocks.values())
    normal = statistics.NormalDist()
    z_total = normal.inv_cdf(.975) + normal.inv_cdf(.8)
    metrics = {}
    for name, source, expected in (
        ("online_root_call_ratio", "operation_count",
         analysis["operation_count"]["n41_over_n53_gain"]["geomean_ratio"]),
        ("exploratory_online_wall_ratio", "online_wall_exploratory",
         analysis["online_wall_exploratory"]["n41_over_n53_gain"]["geomean_ratio"]),
    ):
        key = "online_root_call_ratio" if source == "operation_count" else "exploratory_online_ratio"
        logs = [math.log(blocks[41, f"T{index:03d}"][key] /
                         blocks[53, f"T{index:03d}"][key])
                for index in range(1, 34)]
        mean, sd = statistics.mean(logs), statistics.stdev(logs)
        assert math.isclose(math.exp(mean), expected, rel_tol=1e-12)
        metrics[name] = {
            "observed_n41_over_n53_gain": math.exp(mean),
            "mean_paired_log_gain_difference": mean,
            "sample_sd_paired_log_gain_difference": sd,
            "provisional_total_distinct_pairs_for_observed_effect":
                estimate(sd, mean, z_total),
            "provisional_total_distinct_pairs_by_effect": {
                "1.25x": estimate(sd, math.log(1.25), z_total),
                "1.5x": estimate(sd, math.log(1.5), z_total),
                "2x": estimate(sd, math.log(2), z_total),
            },
        }
    return {
        "kind": "s3_crosscurve_extension33_exploratory_sample_size_v1",
        "analysis_script_sha256": sha(Path(__file__).read_bytes()),
        "analysis_sha256": sha(analysis_raw),
        "prior_analysis_sha256": sha(blocks_raw),
        "current_distinct_matched_indices": 33,
        "frozen_distinct_matched_indices_available": 200,
        "method": "normal_approximation_to_mean_of_paired_log_gain_differences",
        "sidedness": "two_sided", "alpha": "0.05", "target_power": "0.80",
        "z975_plus_z80": z_total,
        "formula": "ceil(((z_0.975+z_0.80)*sample_sd/log(effect_ratio))**2)",
        "status": "exploratory_planning_estimate",
        "interpretation": "More distinct matched target indices estimate a curve-wide difference; same-target repetitions assess wall-time noise. Re-estimate wall variance on an isolated host before powering a controlled comparison.",
        "metrics": metrics,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = make()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["metrics"], sort_keys=True))
