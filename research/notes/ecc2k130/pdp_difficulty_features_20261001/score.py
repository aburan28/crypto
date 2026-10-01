#!/usr/bin/env python3
"""Score the frozen public-point PDP difficulty protocol.

Training choice is written to disk before any holdout target file is parsed.
Holdout bytes may already have been hashed; hashing is the custody check and
does not compute a feature or a correlation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import statistics
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "disjoint_cold_v2_outcome_20261001/evidence_run_36803331080"
FEATURES = ("x_weight", "signed_y_weight", "frobenius_x_min_weight")
PERMUTATION_ALGORITHM = (
    "Python random.Random(seed).shuffle, 2000 successive Fisher-Yates "
    "shuffles of the archive-order probe list; average-tie Spearman"
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def average_ranks(values: list[int]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    i = 0
    n = len(values)
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(a: list[float], b: list[float]) -> float:
    n = len(a)
    ma = sum(a) / n
    mb = sum(b) / n
    var_a = var_b = cov = 0.0
    for x, y in zip(a, b):
        dx = x - ma
        dy = y - mb
        var_a += dx * dx
        var_b += dy * dy
        cov += dx * dy
    if var_a == 0.0 or var_b == 0.0:
        return 0.0
    return cov / (var_a ** 0.5 * var_b ** 0.5)


def spearman(feature_counts: list[int], probes: list[int]) -> float:
    return pearson(average_ranks(feature_counts), average_ranks(probes))


def square(value: int, degree: int, low_terms: list[int]) -> int:
    squared = 0
    bit = value
    position = 0
    while bit:
        if bit & 1:
            squared |= 1 << (2 * position)
        bit >>= 1
        position += 1
    for exponent in range(2 * (degree - 1), degree - 1, -1):
        if (squared >> exponent) & 1:
            squared ^= 1 << exponent
            shift = exponent - degree
            for term in low_terms:
                squared ^= 1 << (term + shift)
    if squared >> degree:
        raise AssertionError("reduction left a bit at or above the degree")
    return squared


def frobenius_min_weight(x: int, degree: int, low_terms: list[int]) -> int:
    image = x
    best = image.bit_count()
    for _ in range(degree - 1):
        image = square(image, degree, low_terms)
        weight = image.bit_count()
        if weight < best:
            best = weight
    if square(image, degree, low_terms) != x:
        raise AssertionError("Frobenius orbit did not return after n squares")
    return best


def self_test() -> None:
    assert abs(spearman([1, 2, 3], [1, 2, 3]) - 1.0) < 1e-12
    assert abs(spearman([1, 2, 3], [3, 2, 1]) + 1.0) < 1e-12
    assert average_ranks([1, 1, 2]) == [1.5, 1.5, 3.0]
    assert square(1, 41, [0, 3]) == 1
    assert frobenius_min_weight(1, 41, [0, 3]) == 1
    value = (1 << 40) | (1 << 3) | 1
    image = value
    for _ in range(41):
        image = square(image, 41, [0, 3])
    assert image == value


def load_jsonl(data: bytes) -> list[dict]:
    return [json.loads(line) for line in data.splitlines() if line]


def custody(config: dict) -> dict[str, dict[str, bytes]]:
    manifest_bytes = (ARCHIVE / "MANIFEST.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if sha(manifest_bytes) != config["manifest_sha256"]:
        raise SystemExit("FAIL manifest hash")
    if manifest["source_head"] != config["source_head"]:
        raise SystemExit("FAIL source head")
    kept: dict[str, dict[str, bytes]] = {}
    for cell, spec in config["cells"].items():
        entry = manifest["cases"][cell]
        if entry["status"] != "PASS" or not entry["second_host_matches_hosted"]:
            raise SystemExit(f"FAIL {cell} receipt")
        if spec["raw_sha256"] != entry["raw_sha256"]:
            raise SystemExit(f"FAIL {cell} raw pin")
        path = ARCHIVE / entry["raw_path"]
        data = path.read_bytes()
        if len(data) != entry["raw_bytes"] or sha(data) != spec["raw_sha256"]:
            raise SystemExit(f"FAIL {cell} raw bytes")
        prefix = f"disjoint-cold-v2-{cell}/"
        seen: set[str] = set()
        needed: dict[str, bytes] = {}
        with tarfile.open(path, "r:gz") as archive:
            for member in archive:
                if not member.isfile() or not member.name.startswith(prefix):
                    raise SystemExit(f"FAIL member {member.name}")
                relative = member.name[len(prefix):]
                if relative in seen or relative not in entry["member_sha256"]:
                    raise SystemExit(f"FAIL member set {relative}")
                payload = archive.extractfile(member).read()
                if sha(payload) != entry["member_sha256"][relative]:
                    raise SystemExit(f"FAIL member hash {relative}")
                seen.add(relative)
                if relative.endswith((".target.jsonl", ".stdout.jsonl", ".base.jsonl")) and "/b" in relative:
                    needed[relative] = payload
        if seen != set(entry["member_sha256"]):
            raise SystemExit(f"FAIL {cell} member coverage")
        kept[cell] = needed
    return kept


def rows_for(kept: dict[str, bytes], cell: str, block: int, spec: dict) -> list[dict]:
    stem = f"{cell}/b{block:02d}_ic_"
    arms = {}
    for arm in ("a", "b"):
        targets = load_jsonl(kept[stem + arm + ".target.jsonl"])
        summary = load_jsonl(kept[stem + arm + ".stdout.jsonl"])
        base = load_jsonl(kept[stem + arm + ".base.jsonl"])
        if len(summary) != 1 or len(base) != 1 or len(targets) != 1024:
            raise SystemExit(f"FAIL {cell} block {block} {arm} shape")
        if summary[0]["targets_solved"] != 1024 or summary[0]["targets_failed"] != 0:
            raise SystemExit(f"FAIL {cell} block {block} {arm} solves")
        if summary[0]["base_hash"] != base[0]["base_hash"]:
            raise SystemExit(f"FAIL {cell} block {block} {arm} base hash")
        if base[0]["n"] != spec["n"] or base[0]["orbit_columns"] != spec["K"]:
            raise SystemExit(f"FAIL {cell} block {block} {arm} curve header")
        if base[0]["field_modulus_low_terms"] != spec["field_modulus_low_terms"]:
            raise SystemExit(f"FAIL {cell} block {block} {arm} modulus")
        if int(base[0]["subgroup_order"]) != spec["subgroup_order"]:
            raise SystemExit(f"FAIL {cell} block {block} {arm} subgroup")
        by_index = {}
        for record in targets:
            if record["kind"] != "compact_orbit_dlp_target" or record["exit_code"] != 0:
                raise SystemExit(f"FAIL {cell} block {block} {arm} target status")
            if record["group_verified"] is not True:
                raise SystemExit(f"FAIL {cell} block {block} {arm} verification")
            if record["published_q"] != record["target"]:
                raise SystemExit(f"FAIL {cell} block {block} {arm} published point")
            if not isinstance(record["probes"], int) or isinstance(record["probes"], bool):
                raise SystemExit(f"FAIL {cell} block {block} {arm} probes")
            if not isinstance(record["recovered_scalar"], int) or isinstance(record["recovered_scalar"], bool):
                raise SystemExit(f"FAIL {cell} block {block} {arm} scalar")
            index = record["fixture_index"]
            if index in by_index:
                raise SystemExit(f"FAIL {cell} block {block} {arm} duplicate")
            by_index[index] = record
        if set(by_index) != set(range(1024)):
            raise SystemExit(f"FAIL {cell} block {block} {arm} coverage")
        arms[arm] = (by_index, summary[0]["base_hash"])
    if arms["a"][1] != arms["b"][1]:
        raise SystemExit(f"FAIL {cell} block {block} arm base hash")
    rows = []
    for index in range(1024):
        left = arms["a"][0][index]
        right = arms["b"][0][index]
        if left["published_q"] != right["published_q"] or left["probes"] != right["probes"]:
            raise SystemExit(f"FAIL {cell} block {block} A/B point or probes")
        if left["recovered_scalar"] != right["recovered_scalar"]:
            raise SystemExit(f"FAIL {cell} block {block} A/B scalar")
        x, y = left["published_q"]
        rows.append({"x": x, "y": y, "probes": left["probes"]})
    return rows


def feature_count(row: dict, name: str, degree: int, low_terms: list[int]) -> int:
    if name == "x_weight":
        return row["x"].bit_count()
    if name == "signed_y_weight":
        return min(row["y"].bit_count(), (row["x"] ^ row["y"]).bit_count())
    if name == "frobenius_x_min_weight":
        return frobenius_min_weight(row["x"], degree, low_terms)
    raise AssertionError(name)


def cohort_metrics(rows: list[dict], name: str, sign: int, degree: int, low_terms: list[int], seed: int, permutations: int) -> dict:
    counts = [feature_count(row, name, degree, low_terms) for row in rows]
    probes = [row["probes"] for row in rows]
    observed = spearman(counts, probes)
    signed = sign * observed
    labels = probes.copy()
    rng = random.Random(seed)
    extremes = 0
    score_hex: list[str] = []
    for _ in range(permutations):
        rng.shuffle(labels)
        value = sign * spearman(counts, labels)
        if value >= signed:
            extremes += 1
        score_hex.append(float.hex(value))
    ordered = sorted(range(len(rows)), key=lambda i: ((sign * counts[i]), rows[i]["x"], rows[i]["y"]))
    take = len(rows) // 4
    selected = [probes[i] for i in ordered[:take]]
    median_all = statistics.median(probes)
    median_selected = statistics.median(selected)
    return {
        "n": len(rows),
        "selected_n": take,
        "observed_spearman": observed,
        "signed_spearman": signed,
        "permutation_p": (1 + extremes) / (permutations + 1),
        "permutation_scores_sha256": sha(("\n".join(score_hex) + "\n").encode()),
        "median_probes_all": median_all,
        "median_probes_selected": median_selected,
        "selected_quartile_median_probes_ratio": median_selected / median_all,
        "feature_evaluations": len(rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--choice-out", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.choice_out.exists() or args.out.exists():
        raise SystemExit("refuse to overwrite a scored outcome")
    started = time.perf_counter()
    cpu = time.process_time()
    self_test()
    config_bytes = (HERE / "CONFIG.json").read_bytes()
    config = json.loads(config_bytes)
    kept = custody(config)
    n41 = config["cells"]["n41_L1024"]
    train: list[dict] = []
    for block in n41["train_blocks"]:
        train.extend(rows_for(kept["n41_L1024"], "n41_L1024", block, n41))
    if len(train) != 3072:
        raise SystemExit("FAIL training row count")
    correlations = {}
    evaluations = 0
    for name in FEATURES:
        counts = [feature_count(row, name, n41["n"], n41["field_modulus_low_terms"]) for row in train]
        correlations[name] = spearman(counts, [row["probes"] for row in train])
        evaluations += len(train)
    if all(value == 0.0 for value in correlations.values()):
        predictor_status = "NO_PREDICTOR"
        chosen = FEATURES[0]
    else:
        predictor_status = "CHOSEN"
        chosen = max(FEATURES, key=lambda name: (abs(correlations[name]), -FEATURES.index(name)))
    sign = 1 if correlations[chosen] >= 0.0 else -1
    choice = {
        "schema": "ecc2k130-pdp-difficulty-train-choice-v1",
        "predictor_status": predictor_status,
        "chosen_feature": chosen,
        "training_sign": sign,
        "training_correlations": correlations,
        "training_rows": len(train),
        "holdout_parsed": False,
    }
    args.choice_out.write_text(json.dumps(choice, indent=2, sort_keys=True) + "\n")
    args.choice_out.chmod(0o644)
    cohorts = {}
    for cell_name, label in (("n41_L1024", "n41_holdout"), ("n53_L1024", "n53_holdout")):
        spec = config["cells"][cell_name]
        rows: list[dict] = []
        for block in spec["evaluation_blocks"]:
            rows.extend(rows_for(kept[cell_name], cell_name, block, spec))
        metrics = cohort_metrics(
            rows, chosen, sign, spec["n"], spec["field_modulus_low_terms"],
            spec["permutation_seed"], config["permutations"],
        )
        evaluations += metrics["feature_evaluations"]
        gate = config["positive_gate"]
        metrics["passes_gate"] = (
            metrics["signed_spearman"] >= gate["minimum_signed_spearman"]
            and metrics["permutation_p"] <= gate["maximum_permutation_p"]
            and metrics["selected_quartile_median_probes_ratio"] <= gate["maximum_selected_quartile_median_probes_ratio"]
        )
        cohorts[label] = metrics
    passed = [row["passes_gate"] for row in cohorts.values()]
    if all(passed):
        decision = "TRANSFER_WORTHY_PREDICTOR"
    elif any(passed):
        decision = "MIXED"
    else:
        decision = "NEGATIVE"
    elapsed = time.perf_counter() - started
    cpu_elapsed = time.process_time() - cpu
    if elapsed > config["producer_limit_seconds"]:
        raise SystemExit("CENSORED producer wall limit")
    outcome = {
        "schema": "ecc2k130-pdp-difficulty-outcome-v1",
        "status": "PASS",
        "classification": config["classification"],
        "decision": decision,
        "predictor_status": predictor_status,
        "chosen_feature": chosen,
        "training_sign": sign,
        "training_correlations": correlations,
        "cohorts": cohorts,
        "feature_evaluations": evaluations,
        "permutation_algorithm": PERMUTATION_ALGORITHM,
        "median_definition": "statistics.median: average of the two central values when the count is even",
        "tie_break": "predicted ease, then ascending public (x, y)",
        "zero_variance_spearman": 0.0,
        "config_sha256": sha(config_bytes),
        "source_head": config["source_head"],
        "host": {
            "system": platform.system(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "elapsed_seconds": elapsed,
            "cpu_seconds": cpu_elapsed,
        },
    }
    args.out.write_text(json.dumps(outcome, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"decision": decision, "chosen_feature": chosen, "predictor_status": predictor_status}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
