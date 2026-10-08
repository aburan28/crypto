#!/usr/bin/env python3
"""Describe target-to-target spread without promoting a CPU speedup claim."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
HOLDOUT = HERE / "holdout"
PAIRS = HOLDOUT / "paired_summary_exploratory.json"
REPLAY = HOLDOUT / "independent_sage_replay.json"
OUT = HOLDOUT / "uncertainty_exploratory.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quantile(sorted_values: list[float], fraction: float) -> float:
    position = fraction * (len(sorted_values) - 1)
    low = int(position)
    high = min(low + 1, len(sorted_values) - 1)
    return sorted_values[low] * (high - position) + sorted_values[high] * (position - low)


def main() -> None:
    if OUT.exists():
        raise FileExistsError("refusing to replace uncertainty artifact")
    panel = json.loads(PAIRS.read_text())
    replay = json.loads(REPLAY.read_text())
    assert replay["target_count"] == len(panel["pairs"]) == 12
    assert replay["all_target_scalars_independently_recovered"] is True
    ratios = [pair["exploratory_online_ratio"] for pair in panel["pairs"]]
    deltas = [pair["baseline_online_ms"] - pair["candidate_online_ms"]
              for pair in panel["pairs"]]
    rng = random.Random(20261007)
    bootstrap = sorted(statistics.median(rng.choices(ratios, k=len(ratios)))
                       for _ in range(20000))
    record = {
        "kind": "conditional_fresh_target_exploratory_spread",
        "target_count": len(ratios),
        "paired_ratio_median": statistics.median(ratios),
        "paired_ratio_range": [min(ratios), max(ratios)],
        "paired_online_delta_ms_median": statistics.median(deltas),
        "paired_online_delta_ms_range": [min(deltas), max(deltas)],
        "bootstrap_median_ratio_percentile_95_interval": [
            quantile(bootstrap, 0.025), quantile(bootstrap, 0.975)
        ],
        "bootstrap_seed": 20261007,
        "bootstrap_resamples": 20000,
        "interpretation": (
            "Descriptive resampling of twelve fixed paired targets on an "
            "unisolated host. This interval is not an auditable CPU speedup "
            "claim or a population guarantee."
        ),
        "paired_summary_sha256": sha(PAIRS),
        "independent_sage_replay_sha256": sha(REPLAY),
    }
    OUT.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, sort_keys=True))


if __name__ == "__main__":
    main()
