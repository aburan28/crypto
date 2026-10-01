#!/usr/bin/env python3
"""Independent recomputation of the frozen PDP difficulty score.

This file does not import the producer. It repeats the custody check, the
training choice, both holdout cohorts, and the permutation hashes, then
compares those scientific fields to the archived choice and outcome.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import statistics
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "disjoint_cold_v2_outcome_20261001/evidence_run_36803331080"
FEATURES = ("x_weight", "signed_y_weight", "frobenius_x_min_weight")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def ranks(values: list[int]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    start = 0
    while start < len(values):
        stop = start
        while stop + 1 < len(values) and values[order[stop + 1]] == values[order[start]]:
            stop += 1
        mid = (start + stop + 2) / 2.0
        for cursor in range(start, stop + 1):
            out[order[cursor]] = mid
        start = stop + 1
    return out


def correlation(left: list[float], right: list[float]) -> float:
    count = len(left)
    mean_left = sum(left) / count
    mean_right = sum(right) / count
    left_ss = right_ss = product = 0.0
    for a, b in zip(left, right):
        da = a - mean_left
        db = b - mean_right
        left_ss += da * da
        right_ss += db * db
        product += da * db
    if left_ss == 0.0 or right_ss == 0.0:
        return 0.0
    return product / (left_ss ** 0.5 * right_ss ** 0.5)


def rank_correlation(counts: list[int], probes: list[int]) -> float:
    return correlation(ranks(counts), ranks(probes))


def step_square(element: int, degree: int, low_terms: list[int]) -> int:
    expanded = 0
    cursor = element
    index = 0
    while cursor:
        if cursor & 1:
            expanded |= 1 << (2 * index)
        cursor >>= 1
        index += 1
    for exponent in range(2 * (degree - 1), degree - 1, -1):
        if (expanded >> exponent) & 1:
            expanded ^= 1 << exponent
            shift = exponent - degree
            for term in low_terms:
                expanded ^= 1 << (term + shift)
    if expanded >> degree:
        raise AssertionError("unreduced square")
    return expanded


def min_frobenius_weight(element: int, degree: int, low_terms: list[int]) -> int:
    image = element
    lightest = image.bit_count()
    for _ in range(degree - 1):
        image = step_square(image, degree, low_terms)
        lightest = min(lightest, image.bit_count())
    if step_square(image, degree, low_terms) != element:
        raise AssertionError("Frobenius did not close")
    return lightest


def parse_lines(payload: bytes) -> list[dict]:
    return [json.loads(line) for line in payload.splitlines() if line]


def verified_members(config: dict) -> dict[str, dict[str, bytes]]:
    manifest_bytes = (ARCHIVE / "MANIFEST.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    assert digest(manifest_bytes) == config["manifest_sha256"]
    assert manifest["source_head"] == config["source_head"]
    stored: dict[str, dict[str, bytes]] = {}
    for cell, spec in config["cells"].items():
        entry = manifest["cases"][cell]
        assert entry["status"] == "PASS" and entry["second_host_matches_hosted"]
        assert spec["raw_sha256"] == entry["raw_sha256"]
        blob = (ARCHIVE / entry["raw_path"]).read_bytes()
        assert len(blob) == entry["raw_bytes"] and digest(blob) == spec["raw_sha256"]
        prefix = f"disjoint-cold-v2-{cell}/"
        found: set[str] = set()
        selected: dict[str, bytes] = {}
        with tarfile.open(fileobj=__import__("io").BytesIO(blob), mode="r:gz") as bundle:
            for member in bundle:
                assert member.isfile() and member.name.startswith(prefix)
                relative = member.name[len(prefix):]
                assert relative not in found and relative in entry["member_sha256"]
                payload = bundle.extractfile(member).read()
                assert digest(payload) == entry["member_sha256"][relative]
                found.add(relative)
                if relative.endswith((".target.jsonl", ".stdout.jsonl", ".base.jsonl")):
                    selected[relative] = payload
        assert found == set(entry["member_sha256"])
        stored[cell] = selected
    return stored


def paired_rows(stored: dict[str, bytes], cell: str, block: int, spec: dict) -> list[dict]:
    paired = []
    decoded = {}
    for arm in ("a", "b"):
        stem = f"{cell}/b{block:02d}_ic_{arm}"
        targets = parse_lines(stored[stem + ".target.jsonl"])
        summary = parse_lines(stored[stem + ".stdout.jsonl"])
        base = parse_lines(stored[stem + ".base.jsonl"])
        assert len(targets) == 1024 and len(summary) == 1 and len(base) == 1
        assert summary[0]["targets_solved"] == 1024 and summary[0]["targets_failed"] == 0
        assert summary[0]["base_hash"] == base[0]["base_hash"]
        assert base[0]["n"] == spec["n"] and base[0]["orbit_columns"] == spec["K"]
        assert base[0]["field_modulus_low_terms"] == spec["field_modulus_low_terms"]
        assert int(base[0]["subgroup_order"]) == spec["subgroup_order"]
        by_index = {}
        for record in targets:
            assert record["kind"] == "compact_orbit_dlp_target"
            assert record["exit_code"] == 0 and record["group_verified"] is True
            assert record["published_q"] == record["target"]
            assert type(record["probes"]) is int and type(record["recovered_scalar"]) is int
            assert record["fixture_index"] not in by_index
            by_index[record["fixture_index"]] = record
        assert set(by_index) == set(range(1024))
        decoded[arm] = (by_index, summary[0]["base_hash"])
    assert decoded["a"][1] == decoded["b"][1]
    for index in range(1024):
        left = decoded["a"][0][index]
        right = decoded["b"][0][index]
        assert left["published_q"] == right["published_q"]
        assert left["probes"] == right["probes"]
        assert left["recovered_scalar"] == right["recovered_scalar"]
        x_coord, y_coord = left["published_q"]
        paired.append((x_coord, y_coord, left["probes"]))
    return paired


def count_for(point: tuple[int, int, int], name: str, degree: int, low_terms: list[int]) -> int:
    x_coord, y_coord, _probes = point
    if name == "x_weight":
        return x_coord.bit_count()
    if name == "signed_y_weight":
        return min(y_coord.bit_count(), (x_coord ^ y_coord).bit_count())
    return min_frobenius_weight(x_coord, degree, low_terms)


def evaluate(points: list[tuple[int, int, int]], name: str, sign: int, degree: int, low_terms: list[int], seed: int, draws: int) -> dict:
    counts = [count_for(point, name, degree, low_terms) for point in points]
    probes = [point[2] for point in points]
    observed = rank_correlation(counts, probes)
    signed = sign * observed
    shuffled = probes.copy()
    generator = random.Random(seed)
    wins = 0
    rendered: list[str] = []
    for _ in range(draws):
        generator.shuffle(shuffled)
        draw = sign * rank_correlation(counts, shuffled)
        wins += draw >= signed
        rendered.append(float.hex(draw))
    order = sorted(range(len(points)), key=lambda i: ((-sign * counts[i]), points[i][0], points[i][1]))
    width = len(points) // 4
    chosen_probes = [probes[i] for i in order[:width]]
    all_median = statistics.median(probes)
    part_median = statistics.median(chosen_probes)
    return {
        "n": len(points),
        "selected_n": width,
        "observed_spearman": observed,
        "signed_spearman": signed,
        "permutation_p": (1 + wins) / (draws + 1),
        "permutation_scores_sha256": digest(("\n".join(rendered) + "\n").encode()),
        "median_probes_all": all_median,
        "median_probes_selected": part_median,
        "selected_quartile_median_probes_ratio": part_median / all_median,
        "feature_evaluations": len(points),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--choice", type=Path, required=True)
    parser.add_argument("--outcome", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    assert abs(rank_correlation([1, 2, 3], [1, 2, 3]) - 1.0) < 1e-12
    assert abs(rank_correlation([1, 2, 3], [3, 2, 1]) + 1.0) < 1e-12
    assert ranks([1, 1, 2]) == [1.5, 1.5, 3.0]
    config = json.loads((HERE / "CONFIG.json").read_bytes())
    stored = verified_members(config)
    n41 = config["cells"]["n41_L1024"]
    training: list[tuple[int, int, int]] = []
    for block in n41["train_blocks"]:
        training.extend(paired_rows(stored["n41_L1024"], "n41_L1024", block, n41))
    assert len(training) == 3072
    found = {}
    for name in FEATURES:
        counts = [count_for(point, name, n41["n"], n41["field_modulus_low_terms"]) for point in training]
        found[name] = rank_correlation(counts, [point[2] for point in training])
    if all(value == 0.0 for value in found.values()):
        status = "NO_PREDICTOR"
        chosen = FEATURES[0]
    else:
        status = "CHOSEN"
        chosen = max(FEATURES, key=lambda name: (abs(found[name]), -FEATURES.index(name)))
    sign = 1 if found[chosen] >= 0.0 else -1
    recorded_choice = json.loads(args.choice.read_text())
    assert recorded_choice["holdout_parsed"] is False
    assert recorded_choice["predictor_status"] == status
    assert recorded_choice["chosen_feature"] == chosen
    assert recorded_choice["training_sign"] == sign
    assert recorded_choice["training_correlations"] == found
    archived = json.loads(args.outcome.read_text())
    assert archived["decision"] in {"NEGATIVE", "MIXED", "TRANSFER_WORTHY_PREDICTOR"}
    assert archived["chosen_feature"] == chosen
    assert archived["training_sign"] == sign
    assert archived["training_correlations"] == found
    assert archived["predictor_status"] == status
    evaluations = 3 * len(training)
    gate = config["positive_gate"]
    passes = []
    for cell_name, label in (("n41_L1024", "n41_holdout"), ("n53_L1024", "n53_holdout")):
        spec = config["cells"][cell_name]
        points: list[tuple[int, int, int]] = []
        for block in spec["evaluation_blocks"]:
            points.extend(paired_rows(stored[cell_name], cell_name, block, spec))
        metrics = evaluate(
            points, chosen, sign, spec["n"], spec["field_modulus_low_terms"],
            spec["permutation_seed"], config["permutations"],
        )
        metrics["passes_gate"] = (
            metrics["signed_spearman"] >= gate["minimum_signed_spearman"]
            and metrics["permutation_p"] <= gate["maximum_permutation_p"]
            and metrics["selected_quartile_median_probes_ratio"] <= gate["maximum_selected_quartile_median_probes_ratio"]
        )
        passes.append(metrics["passes_gate"])
        evaluations += metrics["feature_evaluations"]
        saved = archived["cohorts"][label]
        for key, value in metrics.items():
            assert saved[key] == value, (label, key, saved[key], value)
    if all(passes):
        decision = "TRANSFER_WORTHY_PREDICTOR"
    elif any(passes):
        decision = "MIXED"
    else:
        decision = "NEGATIVE"
    assert archived["decision"] == decision
    assert archived["feature_evaluations"] == evaluations
    assert archived["status"] == "PASS"
    assert archived["classification"] == "diagnostic_only_no_speed_or_n131_claim"
    elapsed = time.perf_counter() - started
    if elapsed > config["checker_limit_seconds"]:
        raise SystemExit("CENSORED checker wall limit")
    print(json.dumps({"status": "PASS", "decision": decision, "elapsed_seconds": elapsed}, sort_keys=True))


if __name__ == "__main__":
    main()
