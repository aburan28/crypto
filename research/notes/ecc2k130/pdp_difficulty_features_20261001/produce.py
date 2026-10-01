#!/usr/bin/env python3
"""One-shot, held-out public-point predictor experiment; never use scalar labels."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import resource
import signal
import statistics
import subprocess
import sys
import tarfile
import time
import traceback

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "disjoint_cold_v2_outcome_20261001/evidence_run_36803331080"
FEATURES = ("x_weight", "signed_y_weight", "frobenius_x_min_weight")
FEATURE_METER = {"evaluations": 0, "field_squarings": 0, "cpu_seconds": 0.0}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: object) -> None:
    data = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    with path.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())


def write_jsonl(path: Path, rows: list[dict]) -> None:
    data = "".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in rows).encode()
    with path.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())


def square(value: int, n: int, low_terms: list[int]) -> int:
    spread = 0
    for bit in range(n):
        if (value >> bit) & 1:
            spread ^= 1 << (2 * bit)
    modulus = 1 << n
    for bit in low_terms:
        modulus |= 1 << bit
    for bit in range(2 * n - 2, n - 1, -1):
        if (spread >> bit) & 1:
            spread ^= modulus << (bit - n)
    assert spread < 1 << n
    return spread


def features(x: int, y: int, n: int, low_terms: list[int]) -> dict[str, int]:
    started = time.process_time()
    assert 0 <= x < 1 << n and 0 <= y < 1 << n
    value = x
    smallest = n + 1
    for _ in range(n):
        smallest = min(smallest, value.bit_count())
        value = square(value, n, low_terms)
    assert value == x, "Frobenius cycle mismatch"
    value = {"x_weight": x.bit_count(),
             "signed_y_weight": min(y.bit_count(), (x ^ y).bit_count()),
             "frobenius_x_min_weight": smallest}
    FEATURE_METER["evaluations"] += 1
    FEATURE_METER["field_squarings"] += n
    FEATURE_METER["cpu_seconds"] += time.process_time() - started
    return value


def tied_ranks(values: list[int]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    left = 0
    while left < len(order):
        right = left + 1
        while right < len(order) and values[order[right]] == values[order[left]]:
            right += 1
        average = (left + 1 + right) / 2
        for position in range(left, right):
            ranks[order[position]] = average
        left = right
    return ranks


def correlation(x: list[float], y: list[float]) -> float:
    assert len(x) == len(y) and len(x) > 1
    mx = statistics.fmean(x)
    my = statistics.fmean(y)
    sx = sum((v - mx) ** 2 for v in x)
    sy = sum((v - my) ** 2 for v in y)
    if sx == 0 or sy == 0:
        return 0.0
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sx * sy)


def spearman(features_: list[int], probes: list[int]) -> float:
    return correlation(tied_ranks(features_), tied_ranks(probes))


def checked_member(archive: tarfile.TarFile, manifest: dict, cell: str,
                   relative: str) -> bytes:
    expected = manifest["cases"][cell]["member_sha256"][relative]
    member = archive.getmember(f"disjoint-cold-v2-{cell}/{relative}")
    assert member.isfile()
    data = archive.extractfile(member).read()
    assert sha(data) == expected, relative
    return data


def load_rows(cell: str, blocks: list[int], config: dict, manifest: dict) -> list[dict]:
    spec = config["cells"][cell]
    entry = manifest["cases"][cell]
    path = ARCHIVE / entry["raw_path"]
    raw = path.read_bytes()
    assert sha(raw) == spec["raw_sha256"] == entry["raw_sha256"]
    assert len(raw) == entry["raw_bytes"]
    observations = []
    with tarfile.open(path, "r:gz") as archive:
        run = json.loads(checked_member(archive, manifest, cell, f"{cell}/cold_run.json"))
        assert run["status"] == "PASS" and run["host"]["git_head"] == config["source_head"]
        for block in blocks:
            base = []
            traces = []
            for arm in ("ic_a", "ic_b"):
                prefix = f"{cell}/b{block:02d}_{arm}"
                summary_lines = checked_member(archive, manifest, cell, prefix + ".stdout.jsonl").splitlines()
                assert len(summary_lines) == 1
                summary = json.loads(summary_lines[0])
                assert summary["n"] == spec["n"] and summary["orbit_columns"] == spec["K"]
                assert summary["rank"] == spec["K"]
                assert summary["targets_solved"] == 1024 and summary["targets_failed"] == 0
                base.append(summary["base_hash"])
                lines = checked_member(archive, manifest, cell, prefix + ".target.jsonl").splitlines()
                assert len(lines) == 1024
                traces.append([json.loads(line) for line in lines])
            assert base[0] == base[1]
            for index, (a, b) in enumerate(zip(*traces)):
                assert a["fixture_index"] == b["fixture_index"] == index
                assert a["exit_code"] == b["exit_code"] == 0
                assert a["group_verified"] and b["group_verified"]
                assert a["published_fixture_scalar"] is None and b["published_fixture_scalar"] is None
                assert a["published_q"] == b["published_q"] == a["target"] == b["target"]
                assert a["probes"] == b["probes"]
                assert a["recovered_scalar"] == b["recovered_scalar"]
                assert a["n"] == b["n"] == spec["n"]
                x, y = a["published_q"]
                observations.append({"cell": cell, "block": block, "fixture_index": index,
                                     "x": x, "y": y, "probes": a["probes"],
                                     "features": features(x, y, spec["n"], spec["field_modulus_low_terms"])})
    assert len(observations) == 1024 * len(blocks)
    return observations


def permutation_p(rows: list[dict], feature: str, sign: int,
                  seed: int, count: int) -> tuple[float, list[float]]:
    xranks = tied_ranks([row["features"][feature] for row in rows])
    yranks = tied_ranks([row["probes"] for row in rows])
    observed = sign * correlation(xranks, yranks)
    rng = random.Random(seed)
    scores = []
    exceed = 0
    for _ in range(count):
        shuffled = yranks.copy()
        rng.shuffle(shuffled)
        score = sign * correlation(xranks, shuffled)
        exceed += score >= observed
        scores.append(round(score, 12))
    p = (1 + exceed) / (count + 1)
    return p, scores


def evaluate(rows: list[dict], choice: dict, spec: dict, count: int) -> tuple[dict, list[float]]:
    feature = choice["feature"]
    sign = choice["sign"]
    signed_rho = sign * spearman([r["features"][feature] for r in rows],
                                 [r["probes"] for r in rows])
    ordered = sorted(rows, key=lambda r: (sign * r["features"][feature], r["x"], r["y"]))
    selected = ordered[:len(rows) // 4]
    ratio = statistics.median(r["probes"] for r in selected) / statistics.median(r["probes"] for r in rows)
    p, scores = permutation_p(rows, feature, sign, spec["permutation_seed"], count)
    return {"rows": len(rows), "selected": len(selected), "signed_spearman": round(signed_rho, 12),
            "selected_quartile_median_probes_ratio": round(ratio, 12),
            "selected_median_probes": statistics.median(r["probes"] for r in selected),
            "all_median_probes": statistics.median(r["probes"] for r in rows),
            "permutation_p": p, "permutation_seed": spec["permutation_seed"],
            "permutation_count": count,
            "permutation_scores_sha256": sha((json.dumps(scores, separators=(",", ":")) + "\n").encode())}, scores


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    def timeout(_signum: int, _frame: object) -> None:
        raise TimeoutError("producer wall cap")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(120)
    try:
        config_bytes = (HERE / "CONFIG.json").read_bytes()
        config = json.loads(config_bytes)
        assert tuple(config["feature_tie_order"]) == FEATURES
        assert config["producer_limit_seconds"] == 120
        preflight = [sys.executable, str(HERE / "check_config.py"), "--out", str(out / "input_check.json")]
        with out.joinpath("input_check.stdout.txt").open("xb") as log:
            subprocess.run(preflight, check=True, stdout=log)
        assert (out / "input_check.json").read_bytes() == (HERE / "INPUT_CHECK.json").read_bytes()
        source_preflight = [sys.executable, str(HERE / "source_preflight.py"),
                            "--out", str(out / "source_check.json")]
        with out.joinpath("source_check.stdout.txt").open("xb") as log:
            subprocess.run(source_preflight, check=True, stdout=log)
        assert (out / "source_check.json").read_bytes() == (HERE / "SOURCE_CHECK.json").read_bytes()
        manifest_bytes = (ARCHIVE / "MANIFEST.json").read_bytes()
        assert sha(manifest_bytes) == config["manifest_sha256"]
        manifest = json.loads(manifest_bytes)
        train = load_rows("n41_L1024", config["cells"]["n41_L1024"]["train_blocks"], config, manifest)
        write_jsonl(out / "training.jsonl", train)
        correlations = {f: spearman([row["features"][f] for row in train],
                                    [row["probes"] for row in train]) for f in FEATURES}
        selected_feature = max(FEATURES, key=lambda f: abs(correlations[f]))
        selected_rho = correlations[selected_feature]
        choice = {"feature": selected_feature, "sign": 1 if selected_rho >= 0 else -1,
                  "training_spearman": {f: round(correlations[f], 12) for f in FEATURES},
                  "training_rows": len(train),
                  "classification": "NO_PREDICTOR" if all(v == 0 for v in correlations.values()) else "SELECTED"}
        write_json(out / "train_choice.json", choice)
        evaluations = {}
        gate = config["positive_gate"]
        passed = []
        for cell in ("n41_L1024", "n53_L1024"):
            spec = config["cells"][cell]
            rows = load_rows(cell, spec["evaluation_blocks"], config, manifest)
            write_jsonl(out / f"{cell}_evaluation.jsonl", rows)
            summary, scores = evaluate(rows, choice, spec, config["permutations"])
            write_json(out / f"{cell}_permutations.json", scores)
            evaluations[cell] = summary
            passed.append(summary["signed_spearman"] >= gate["minimum_signed_spearman"]
                          and summary["permutation_p"] <= gate["maximum_permutation_p"]
                          and summary["selected_quartile_median_probes_ratio"]
                          <= gate["maximum_selected_quartile_median_probes_ratio"])
        after = resource.getrusage(resource.RUSAGE_SELF)
        after_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        rss = max(after.ru_maxrss, after_children.ru_maxrss) * (1 if platform.system() == "Darwin" else 1024)
        evidence_names = ("input_check.json", "input_check.stdout.txt",
                          "source_check.json", "source_check.stdout.txt", "training.jsonl",
                          "train_choice.json", "n41_L1024_evaluation.jsonl",
                          "n53_L1024_evaluation.jsonl", "n41_L1024_permutations.json",
                          "n53_L1024_permutations.json")
        evidence_sha256 = {name: sha((out / name).read_bytes()) for name in evidence_names}
        receipt = {"schema": "ecc2k130-pdp-difficulty-producer-v1",
                   "evidence_sha256": evidence_sha256,
                   "status": "PASS", "decision": "TRANSFER_WORTHY_DIAGNOSTIC" if all(passed) else "NO_TRANSFER_WORTHY_PREDICTOR",
                   "source_sha256": sha(Path(__file__).read_bytes()), "config_sha256": sha(config_bytes),
                   "feature_meter": FEATURE_METER.copy(),
                   "manifest_sha256": sha(manifest_bytes), "choice": choice,
                   "evaluation": evaluations, "host": platform.node(),
                   "python": sys.version, "cpu_seconds": ((after.ru_utime + after.ru_stime) - (before.ru_utime + before.ru_stime)
                   + (after_children.ru_utime + after_children.ru_stime) - (before_children.ru_utime + before_children.ru_stime)),
                   "wall_seconds": time.monotonic() - started, "peak_rss_bytes": rss,
                   "limit_seconds": config["producer_limit_seconds"],
                   "memory_limit_bytes": config["memory_limit_bytes"]}
        assert receipt["wall_seconds"] <= config["producer_limit_seconds"]
        assert rss <= config["memory_limit_bytes"]
        write_json(out / "RESULT.json", receipt)
        print(json.dumps({"status": receipt["status"], "decision": receipt["decision"]}, sort_keys=True))
    except BaseException as exc:
        write_json(out / "FAILURE.json", {"type": type(exc).__name__, "error": str(exc),
                                         "traceback": traceback.format_exc(),
                                         "elapsed_seconds": time.monotonic() - started})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()
