"""Derive rank-243 online S3-call counts from frozen relation traces.

The original 100 N53 targets were first run with rank 237. Three were later
rerun with rank 243 as stress controls; the other 97 rank-243 entries here
are exact trace-derived counterfactuals. A 27-target rank-243 ledger checks
the formula against executed runs.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tarfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARCHIVES = {
    "initial": ("raw_panel.tar.gz", "88d726f4508d6083311115324c311cae5a1f3c4af7ab46a3298da690ee5f3e79"),
    "extension": ("raw_extension.tar.gz", "60fec46ea00db4a2ad7c8a7e76685b305855c75a1568b4865f2d9e33481c4e37"),
    "rank243": ("raw_rank243.tar.gz", "62ce737b1d1e3c37bfcba63a0403d8fc7ef7c60240ae674a6e962b229ac3770b"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_ledger(name: str) -> list[dict]:
    filename, expected_hash = ARCHIVES[name]
    path = HERE / "local" / filename
    assert sha256(path) == expected_hash, f"archive changed: {path}"
    with tarfile.open(path, "r:gz") as archive:
        raw = archive.extractfile("./raw.jsonl")
        assert raw is not None
        return [json.loads(line) for line in raw]


def by_target(rows: list[dict], expected_count: int) -> dict[str, dict[str, dict]]:
    result: dict[str, dict[str, dict]] = {}
    for row in rows:
        if row["n"] != 53:
            continue
        arms = result.setdefault(row["label"], {})
        assert row["arm"] not in arms, (row["label"], row["arm"])
        assert row["verified_native_and_fixture"] and row["phase_ok"]
        arms[row["arm"]] = row
    assert len(result) == expected_count, len(result)
    return result


def fraction_record(value: Fraction) -> dict:
    return {"exact": f"{value.numerator}/{value.denominator}",
            "decimal": float(value)}


def main(output: Path) -> None:
    plan = json.loads((HERE / "rank243" / "FREEZE.json").read_text())
    assert plan["rank_policies"]["deep"] == 243
    assert plan["workers"] == 14 and plan["relation_seed"] == 20260928
    assert sha256(ROOT / "examples/s3_partial_prepared_targetspan.rs") == plan["partial_source_sha256"]
    assert sha256(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs") == plan["reference_source_sha256"]

    measured = by_target(read_ledger("rank243"), 27)
    master = measured["S067"]["reference"]["report"]["rank_relation_witnesses"]
    assert len(master) == 560
    cumulative = [0]
    for witness in master:
        cumulative.append(cumulative[-1] + witness["target_s3_calls"])
    assert cumulative[238] == 7_699_066
    assert cumulative[462] == 15_323_984
    assert cumulative[560] == 18_819_864

    regular_states = measured["S067"]["reference"]["report"]["regular_states"]
    assert regular_states == 3_155_408
    direct_cap = 2 * regular_states  # extract() visits each state at most once.

    def check_reference(row: dict) -> dict:
        report = row["report"]
        attempts = report["rank_attempts"]
        assert 238 <= attempts <= 560
        assert report["rank_attempts_completed"] == attempts
        assert report["rank_relation_witnesses"] == master[:attempts]
        assert report["rank_target_s3_calls"] == cumulative[attempts]
        assert report["target_relation_found"] and report["group_verified"]
        assert report["target_s3_calls"] == 2 * report["target_state_probes"]
        assert 0 < report["target_s3_calls"] <= direct_cap
        assert report["regular_states"] == regular_states
        assert row["binary_sha256"] == plan["reference_binary_sha256"]
        return report

    def predicted_online_calls(reference: dict) -> tuple[int, int]:
        attempts = reference["rank_attempts"]
        tail = cumulative[attempts] - cumulative[462] if attempts > 462 else 0
        return reference["target_s3_calls"] + tail, tail

    validation = []
    for label, arms in sorted(measured.items()):
        assert set(arms) == {"reference", "shallow", "deep"}
        reference = check_reference(arms["reference"])
        deep = arms["deep"]["report"]
        assert arms["deep"]["binary_sha256"] == plan["partial_binary_sha256"]
        assert deep["prepare_rank_requested"] == deep["prepared_rank"] == 243
        assert deep["prepared_rank_reached"] and deep["rank_attempts"] == 462
        assert deep["rank_relation_witnesses"] == master[:462]
        assert deep["target"] == reference["target"]
        assert deep["target_relation_indices"] == reference["target_relation_indices"]
        assert deep["target_s3_calls"] == reference["target_s3_calls"]
        predicted, tail = predicted_online_calls(reference)
        assert deep["fallback_attempts"] == max(0, reference["rank_attempts"] - 462)
        assert deep["fallback_relation_witnesses"] == master[462:reference["rank_attempts"]]
        assert deep["fallback_s3_calls"] == tail
        assert predicted == deep["target_s3_calls"] + deep["fallback_s3_calls"]
        validation.append({"label": label, "family": arms["deep"]["family"],
                           "reference_attempts": reference["rank_attempts"],
                           "measured_rank243_online_s3_calls": predicted,
                           "ratio": fraction_record(Fraction(
                               reference["rank_target_s3_calls"] + reference["target_s3_calls"],
                               predicted))})

    original = by_target(read_ledger("initial") + read_ledger("extension"), 100)
    stress_overlap = {f"T{label[1:]}": label for label, arms in measured.items()
                      if arms["deep"]["family"] == "full_rank_stress_control"}
    assert stress_overlap == {"T067": "S067", "T085": "S085", "T096": "S096"}
    counterfactual = []
    for label, arms in sorted(original.items()):
        assert set(arms) == {"reference", "partial"}
        reference = check_reference(arms["reference"])
        partial = arms["partial"]["report"]
        assert arms["partial"]["binary_sha256"] == plan["partial_binary_sha256"]
        assert partial["prepared_rank"] == 237 and partial["rank_attempts"] == 238
        assert partial["rank_relation_witnesses"] == master[:238]
        assert partial["target"] == reference["target"]
        assert partial["target_relation_indices"] == reference["target_relation_indices"]
        assert partial["target_s3_calls"] == reference["target_s3_calls"]
        assert partial["group_verified"]
        predicted, tail = predicted_online_calls(reference)
        executed_label = stress_overlap.get(label)
        if executed_label is not None:
            executed = measured[executed_label]["deep"]
            assert executed["target_sha256"] == arms["reference"]["target_sha256"]
            assert executed["report"]["target"] == reference["target"]
            assert (executed["report"]["target_s3_calls"] +
                    executed["report"]["fallback_s3_calls"]) == predicted
        ratio = Fraction(reference["rank_target_s3_calls"] + reference["target_s3_calls"],
                         predicted)
        counterfactual.append({
            "label": label,
            "point_sha256": arms["reference"]["target_sha256"],
            "reference_attempts": reference["rank_attempts"],
            "direct_target_s3_calls": reference["target_s3_calls"],
            "reference_online_s3_calls": reference["rank_target_s3_calls"] + reference["target_s3_calls"],
            "rank243_predicted_fallback_attempts": max(0, reference["rank_attempts"] - 462),
            "rank243_predicted_fallback_s3_calls": tail,
            "rank243_predicted_online_s3_calls": predicted,
            "rank243_executed_stress_label": executed_label,
            "predicted_ratio": fraction_record(ratio),
        })

    # These bounds use the loop's structural direct-target cap, rather than
    # the observed maximum over the particular target list.
    min_attempts = min(row["reference_attempts"] for row in counterfactual)
    max_attempts = max(row["reference_attempts"] for row in counterfactual)
    assert min_attempts == 238 and max_attempts == 560
    early_bound = Fraction(cumulative[min_attempts] + direct_cap, direct_cap)
    late_tail_max = cumulative[max_attempts] - cumulative[462]
    late_bound = Fraction(cumulative[462] + late_tail_max + direct_cap,
                          late_tail_max + direct_cap)
    panel_bound = min(early_bound, late_bound)
    assert panel_bound >= 2
    worst = min(counterfactual, key=lambda row: Fraction(row["predicted_ratio"]["exact"]))
    measured_worst = min(validation, key=lambda row: Fraction(row["ratio"]["exact"]))

    result = {
        "kind": "s3_rank243_n53_frozen_panel_root_call_bound_v1",
        "metric": "online S3 root evaluations; two charged per visited index state",
        "claim_type": "exact trace projection for 100 prior targets, including 3 executed rank243 controls; 27 executed rank243 targets validate the formula",
        "controlled_speedup": None,
        "curve_id": "EC1N53Ckb1h738a572b91f2",
        "reference_candidate_id": "IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h5ed3e1d8954b",
        "rank243_candidate_id": "IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd278c9cf111",
        "rank243_plan_sha256": sha256(HERE / "rank243" / "FREEZE.json"),
        "archive_sha256": {name: expected for name, (_, expected) in ARCHIVES.items()},
        "relation_seed": plan["relation_seed"],
        "workers": plan["workers"],
        "regular_states": regular_states,
        "direct_target_s3_call_cap": direct_cap,
        "relation_prefix_s3_calls": {"at_238": cumulative[238], "at_462": cumulative[462],
                                     "at_560": cumulative[560]},
        "relation_tail_463_through_560_s3_calls": late_tail_max,
        "frozen_panel": {
            "prior_target_count": len(counterfactual),
            "trace_only_rank243_target_count": len(counterfactual) - len(stress_overlap),
            "executed_rank243_stress_overlap_count": len(stress_overlap),
            "executed_rank243_stress_labels": stress_overlap,
            "reference_attempt_distribution": dict(sorted(Counter(
                row["reference_attempts"] for row in counterfactual).items())),
            "observed_max_direct_target_s3_calls": max(
                row["direct_target_s3_calls"] for row in counterfactual),
            "early_stop_lower_bound": fraction_record(early_bound),
            "post_rank243_stop_lower_bound": fraction_record(late_bound),
            "uniform_lower_bound": fraction_record(panel_bound),
            "predicted_at_least_2_count": sum(
                Fraction(row["predicted_ratio"]["exact"]) >= 2 for row in counterfactual),
            "minimum_predicted_ratio": worst["predicted_ratio"],
            "minimum_predicted_ratio_label": worst["label"],
            "targets": counterfactual,
        },
        "executed_rank243_validation": {
            "target_count": len(validation),
            "exact_online_call_match_count": len(validation),
            "minimum_measured_ratio": measured_worst["ratio"],
            "minimum_measured_ratio_label": measured_worst["label"],
            "targets": validation,
        },
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "panel_lower_bound": float(panel_bound),
                      "prior_panel_minimum_predicted_ratio": worst["predicted_ratio"]["decimal"],
                      "executed_validation_count": len(validation)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "local" / "rank243_root_call_bound.json")
    main(parser.parse_args().output)
