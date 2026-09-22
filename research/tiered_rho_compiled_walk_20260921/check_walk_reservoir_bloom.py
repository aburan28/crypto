"""Bloom equivalence audit for the single-register reservoir-guided walk."""

from hashlib import sha256
from pathlib import Path
import json

from check_charged_polynomial_construction import Transcript, add, direction_profile
from check_multiline_stripes import prefix_metrics
from check_coverage_guided_chain import SlopeBloom
from check_reservoir_guided_chain import prf, slope
from check_walk_reservoir_guided import (
    BUDGET,
    CONFIGS as ALL_CONFIGS,
    addend_indices,
    bottom_k_indices,
    run as run_exact,
)


CONFIGS = (ALL_CONFIGS[0], ALL_CONFIGS[-1])
SEEDS = tuple(20261400 + i for i in range(32))


def run_bloom(choices, reservoir_size, seed):
    transcript = Transcript()
    materialized = set(transcript.events)
    bloom = SlopeBloom()
    exact_sketch = set()
    false_positives = 0
    false_negatives = 0
    changed_choices = 0
    membership_queries = 0
    selected_addends = []
    for step in range(BUDGET):
        current_index = len(transcript.events) - 1
        witnesses = [
            transcript.events[index]
            for index in bottom_k_indices(transcript.events, reservoir_size, seed)
        ]
        bloom_ranked = []
        exact_ranked = []
        for addend_index in addend_indices(
            len(transcript.events), choices, seed, step
        ):
            candidate = add(
                transcript.events[current_index], transcript.events[addend_index]
            )
            slopes = set()
            for witness in witnesses:
                value = slope(candidate, witness)
                if value is not None:
                    slopes.add(value)
            bloom_gain = 0
            for value in slopes:
                present = bloom.contains(value)
                exact = value in exact_sketch
                membership_queries += 1
                false_positives += int(present and not exact)
                false_negatives += int((not present) and exact)
                bloom_gain += int(not present)
            exact_gain = len(slopes - exact_sketch)
            exploration = prf(
                "walk-explore", seed, step, candidate[0], candidate[1]
            )
            tie = (
                candidate not in materialized,
                exploration,
                -addend_index,
                addend_index,
                candidate,
                slopes,
            )
            bloom_ranked.append((bloom_gain, *tie))
            exact_ranked.append((exact_gain, *tie))
        bloom_best = max(bloom_ranked)
        exact_best = max(exact_ranked)
        changed_choices += int(bloom_best[-3] != exact_best[-3])
        addend_index, candidate, selected_slopes = (
            bloom_best[-3], bloom_best[-2], bloom_best[-1]
        )
        result = transcript.charge_add(
            transcript.events[current_index], transcript.events[addend_index],
            "walk_reservoir_bloom", f"step_{step}_addend_{addend_index}",
        )
        assert result == candidate
        materialized.add(result)
        selected_addends.append(addend_index)
        exact_sketch.update(selected_slopes)
        for value in selected_slopes:
            bloom.insert(value)
    assert false_negatives == 0
    profile = direction_profile(transcript.events)
    prefix = prefix_metrics(transcript.events)
    return {
        "choices": choices,
        "reservoir_size": reservoir_size,
        "seed": seed,
        "profile": profile,
        "prefix_metrics": prefix,
        "bloom": {
            "bits": len(bloom.data) * 8,
            "hashes": 7,
            "membership_queries": membership_queries,
            "hash_probes": bloom.hash_probes,
            "false_positives": false_positives,
            "false_negatives": false_negatives,
            "changed_choices_vs_exact_online_score": changed_choices,
            "bits_set": bloom.bit_sets,
        },
        "selected_addends": selected_addends,
        "events_sha256": sha256(
            json.dumps(transcript.events, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def main():
    rows = []
    summaries = []
    for config in CONFIGS:
        current = []
        for seed in SEEDS:
            bloom = run_bloom(config["choices"], config["reservoir_size"], seed)
            exact = run_exact(
                True, config["choices"], config["reservoir_size"], seed
            )
            row = {
                **config,
                "seed": seed,
                "events_match": bloom["events_sha256"] == exact["events_sha256"],
                "bloom_endpoint": bloom["profile"]["unsigned_useful_finite_slopes"],
                "exact_endpoint": exact["profile"]["unsigned_useful_finite_slopes"],
                "bloom_prefix_numerator": bloom["prefix_metrics"]["censored_mean_numerator"],
                "exact_prefix_numerator": exact["prefix_metrics"]["censored_mean_numerator"],
                "bloom": bloom["bloom"],
            }
            rows.append(row)
            current.append(row)
        summaries.append({
            **config,
            "trials": len(current),
            "matching_event_transcripts": sum(r["events_match"] for r in current),
            "false_positives": sum(r["bloom"]["false_positives"] for r in current),
            "false_negatives": sum(r["bloom"]["false_negatives"] for r in current),
            "changed_choices": sum(r["bloom"]["changed_choices_vs_exact_online_score"] for r in current),
            "mean_endpoint_difference_bloom_minus_exact": sum(
                r["bloom_endpoint"] - r["exact_endpoint"] for r in current
            ) / len(current),
            "mean_censored_addition_difference_bloom_minus_exact": sum(
                r["bloom_prefix_numerator"] - r["exact_prefix_numerator"] for r in current
            ) / (len(current) * 65537),
        })
    report = {
        "evidence_type": "finite_walk_reservoir_bloom_equivalence_audit",
        "field_prime": 65537,
        "budget": BUDGET,
        "configs": CONFIGS,
        "public_seeds": SEEDS,
        "summaries": summaries,
        "rows": rows,
        "limitations": [
            "The one-megabit finite Bloom sketch is lightly loaded.",
            "Projected large-run false positives require independent-uniform-hash assumptions.",
            "Bloom guidance does not certify collisions; exact archive verification remains mandatory.",
            "No asymptotic, novelty, hardware, or ECDLP speedup claim follows.",
        ],
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("walk_reservoir_bloom_checks.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"summaries": summaries, "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
