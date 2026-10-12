#!/usr/bin/env python3
"""Independent field, rank and archive replay of the frozen PDP feature test."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import random
import signal
import statistics
import tarfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ARCHIVE = HERE.parent / "disjoint_cold_v2_outcome_20261001/evidence_run_36803331080"
FEATURES = ("x_weight", "signed_y_weight", "frobenius_x_min_weight")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def point_features(x: int, y: int, degree: int, terms: list[int]) -> dict[str, int]:
    modulus = (1 << degree) | sum(1 << t for t in terms)
    mask = (1 << degree) - 1
    def multiply(a: int, b: int) -> int:
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & (1 << degree):
                a ^= modulus
        return result & mask
    assert 0 <= x <= mask and 0 <= y <= mask
    orbit = [x]
    for _ in range(degree):
        orbit.append(multiply(orbit[-1], orbit[-1]))
    assert orbit[-1] == x
    return {"x_weight": x.bit_count(),
            "signed_y_weight": min(y.bit_count(), (y ^ x).bit_count()),
            "frobenius_x_min_weight": min(z.bit_count() for z in orbit[:-1])}


def average_ranks(values: list[int]) -> list[float]:
    counts = Counter(values)
    start = 1
    by_value = {}
    for value in sorted(counts):
        count = counts[value]
        by_value[value] = (start + start + count - 1) / 2
        start += count
    return [by_value[v] for v in values]


def rho(x: list[float], y: list[float]) -> float:
    assert len(x) == len(y) and len(x) > 1
    ax, ay = statistics.fmean(x), statistics.fmean(y)
    cov = math.fsum((a - ax) * (b - ay) for a, b in zip(x, y))
    vx = math.fsum((a - ax) ** 2 for a in x)
    vy = math.fsum((b - ay) ** 2 for b in y)
    return cov / math.sqrt(vx * vy) if vx and vy else 0.0


def spearman(x: list[int], y: list[int]) -> float:
    return rho(average_ranks(x), average_ranks(y))


def read_member(archive: tarfile.TarFile, cell: str, name: str,
                checksums: dict[str, str]) -> bytes:
    path = f"{cell}/{name}"
    data = archive.extractfile(archive.getmember(f"disjoint-cold-v2-{cell}/{path}")).read()
    assert sha(data) == checksums[path], path
    return data


def check_cohort(cell: str, blocks: list[int], evidence: Path,
                 config: dict, manifest: dict, file_name: str) -> list[dict]:
    path = evidence / file_name
    rows = [json.loads(line) for line in path.read_bytes().splitlines()]
    assert len(rows) == 1024 * len(blocks)
    by_key = {(r["block"], r["fixture_index"]): r for r in rows}
    assert len(by_key) == len(rows)
    spec = config["cells"][cell]
    entry = manifest["cases"][cell]
    packed = (ARCHIVE / entry["raw_path"]).read_bytes()
    assert sha(packed) == entry["raw_sha256"] == spec["raw_sha256"]
    assert len(packed) == entry["raw_bytes"]
    with tarfile.open(ARCHIVE / entry["raw_path"], "r:gz") as archive:
        for block in blocks:
            summaries = []
            target_rows = []
            for arm in ("ic_a", "ic_b"):
                head = f"b{block:02d}_{arm}"
                summaries.append(json.loads(read_member(archive, cell, head + ".stdout.jsonl",
                                                       entry["member_sha256"])))
                target_rows.append([json.loads(line) for line in
                                    read_member(archive, cell, head + ".target.jsonl",
                                                entry["member_sha256"]).splitlines()])
            assert summaries[0]["base_hash"] == summaries[1]["base_hash"]
            assert all(s["rank"] == spec["K"] and s["targets_solved"] == 1024
                       and s["targets_failed"] == 0 for s in summaries)
            assert len(target_rows[0]) == len(target_rows[1]) == 1024
            for index in range(1024):
                a, b = target_rows[0][index], target_rows[1][index]
                assert a["fixture_index"] == b["fixture_index"] == index
                assert a["published_q"] == b["published_q"] == a["target"] == b["target"]
                assert a["probes"] == b["probes"]
                assert a["recovered_scalar"] == b["recovered_scalar"]
                assert a["group_verified"] is True and b["group_verified"] is True
                assert a["exit_code"] == b["exit_code"] == 0
                assert a["published_fixture_scalar"] is None and b["published_fixture_scalar"] is None
                x, y = a["published_q"]
                observed = by_key[(block, index)]
                assert (observed["cell"], observed["x"], observed["y"], observed["probes"]) == (cell, x, y, a["probes"])
                assert observed["features"] == point_features(x, y, spec["n"], spec["field_modulus_low_terms"])
    replay_path = ARCHIVE / entry["second_host_replay_path"]
    assert sha(replay_path.read_bytes()) == entry["second_host_replay_sha256"]
    replay = json.loads(replay_path.read_bytes())
    assert replay["status"] == "PASS" and replay["cell"] == cell
    for check in replay["checks"]:
        if check["arm"] in ("ic_a", "ic_b"):
            assert check["rank"]["rank"] == spec["K"]
            assert check["target_logs_verified"] == 1024
    return rows


def check_evaluation(rows: list[dict], choice: dict, spec: dict,
                     config: dict, evidence: Path, cell: str, saved: dict) -> None:
    feature, sign = choice["feature"], choice["sign"]
    x = [r["features"][feature] for r in rows]
    y = [r["probes"] for r in rows]
    observed_rho = sign * spearman(x, y)
    assert abs(saved["signed_spearman"] - observed_rho) <= 1e-10
    selected = sorted(rows, key=lambda r: (sign * r["features"][feature], r["x"], r["y"]))[:len(rows) // 4]
    ratio = statistics.median(r["probes"] for r in selected) / statistics.median(y)
    assert saved["rows"] == len(rows) and saved["selected"] == len(selected)
    assert abs(saved["selected_quartile_median_probes_ratio"] - ratio) <= 1e-10
    assert saved["selected_median_probes"] == statistics.median(r["probes"] for r in selected)
    assert saved["all_median_probes"] == statistics.median(y)
    rng = random.Random(spec["permutation_seed"])
    rx, ry = average_ranks(x), average_ranks(y)
    scores = []
    exceed = 0
    for _ in range(config["permutations"]):
        perm = ry.copy()
        rng.shuffle(perm)
        score = sign * rho(rx, perm)
        exceed += score >= observed_rho
        scores.append(round(score, 12))
    p = (1 + exceed) / (config["permutations"] + 1)
    assert saved["permutation_seed"] == spec["permutation_seed"]
    assert saved["permutation_count"] == config["permutations"]
    assert saved["permutation_p"] == p
    recorded_scores = json.loads((evidence / f"{cell}_permutations.json").read_bytes())
    assert len(recorded_scores) == len(scores)
    assert all(abs(a - b) <= 1e-10 for a, b in zip(recorded_scores, scores))
    assert saved["permutation_scores_sha256"] == sha((json.dumps(recorded_scores, separators=(",", ":")) + "\n").encode())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    evidence = args.evidence.resolve()
    started = time.monotonic()
    def timeout(_signum: int, _frame: object) -> None:
        raise TimeoutError("independent checker wall cap")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(120)
    config_bytes = (HERE / "CONFIG.json").read_bytes()
    config = json.loads(config_bytes)
    assert tuple(config["feature_tie_order"]) == FEATURES
    manifest_bytes = (ARCHIVE / "MANIFEST.json").read_bytes()
    assert sha(manifest_bytes) == config["manifest_sha256"]
    manifest = json.loads(manifest_bytes)
    assert manifest["source_head"] == config["source_head"]
    assert (evidence / "input_check.json").read_bytes() == (HERE / "INPUT_CHECK.json").read_bytes()
    assert (evidence / "source_check.json").read_bytes() == (HERE / "SOURCE_CHECK.json").read_bytes()
    lock = json.loads((HERE / "SOURCE_LOCK.json").read_bytes())
    for relative, expected in lock["file_sha256"].items():
        assert sha((ROOT / relative).read_bytes()) == expected, relative
    result = json.loads((evidence / "RESULT.json").read_bytes())
    assert result["status"] == "PASS"
    assert set(result["evidence_sha256"]) == {
        "input_check.json", "input_check.stdout.txt",
        "source_check.json", "source_check.stdout.txt", "training.jsonl",
        "train_choice.json", "n41_L1024_evaluation.jsonl",
        "n53_L1024_evaluation.jsonl", "n41_L1024_permutations.json",
        "n53_L1024_permutations.json"}
    for name, expected in result["evidence_sha256"].items():
        assert sha((evidence / name).read_bytes()) == expected, name
    assert result["source_sha256"] == sha((HERE / "produce.py").read_bytes())
    assert result["config_sha256"] == sha(config_bytes)
    assert result["manifest_sha256"] == sha(manifest_bytes)
    assert result["wall_seconds"] <= config["producer_limit_seconds"]
    assert result["peak_rss_bytes"] <= config["memory_limit_bytes"]
    assert result["feature_meter"]["evaluations"] == 10240
    assert result["feature_meter"]["field_squarings"] == 5120 * 41 + 5120 * 53
    assert result["feature_meter"]["cpu_seconds"] >= 0
    train = check_cohort("n41_L1024", [0, 1, 2], evidence, config, manifest, "training.jsonl")
    correlations = {name: spearman([r["features"][name] for r in train],
                                    [r["probes"] for r in train]) for name in FEATURES}
    feature = max(FEATURES, key=lambda f: abs(correlations[f]))
    sign = 1 if correlations[feature] >= 0 else -1
    choice = json.loads((evidence / "train_choice.json").read_bytes())
    assert choice["feature"] == feature and choice["sign"] == sign
    assert choice["training_rows"] == len(train)
    assert choice["classification"] == ("NO_PREDICTOR" if all(v == 0 for v in correlations.values()) else "SELECTED")
    assert choice["training_spearman"] == {f: round(correlations[f], 12) for f in FEATURES}
    assert result["choice"] == choice
    gates = []
    for cell in ("n41_L1024", "n53_L1024"):
        spec = config["cells"][cell]
        rows = check_cohort(cell, spec["evaluation_blocks"], evidence, config, manifest,
                            f"{cell}_evaluation.jsonl")
        saved = result["evaluation"][cell]
        check_evaluation(rows, choice, spec, config, evidence, cell, saved)
        target = config["positive_gate"]
        gates.append(saved["signed_spearman"] >= target["minimum_signed_spearman"]
                     and saved["permutation_p"] <= target["maximum_permutation_p"]
                     and saved["selected_quartile_median_probes_ratio"]
                     <= target["maximum_selected_quartile_median_probes_ratio"])
    all_rows = list(train)
    for cell in ("n41_L1024", "n53_L1024"):
        all_rows += [json.loads(line) for line in (evidence / f"{cell}_evaluation.jsonl").read_bytes().splitlines()]
    assert len({(r["cell"], r["x"], r["y"]) for r in all_rows}) == len(all_rows)
    decision = "TRANSFER_WORTHY_DIAGNOSTIC" if all(gates) else "NO_TRANSFER_WORTHY_PREDICTOR"
    assert result["decision"] == decision
    elapsed = time.monotonic() - started
    usage = resource.getrusage(resource.RUSAGE_SELF)
    rss = usage.ru_maxrss * (1 if platform.system() == "Darwin" else 1024)
    assert elapsed <= config["checker_limit_seconds"]
    assert rss <= config["memory_limit_bytes"]
    receipt = {"schema": "ecc2k130-pdp-difficulty-independent-replay-v1", "status": "PASS",
               "wall_seconds": elapsed, "peak_rss_bytes": rss,
               "decision": decision, "evidence_result_sha256": sha((evidence / "RESULT.json").read_bytes()),
               "config_sha256": sha(config_bytes), "manifest_sha256": sha(manifest_bytes),
               "verified_training_rows": len(train),
               "verified_evaluation_rows": {"n41_L1024": 2048, "n53_L1024": 5120},
               "checker_sha256": sha(Path(__file__).read_bytes())}
    with args.out.open("x") as stream:
        stream.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    signal.alarm(0)
    print(json.dumps({"status": receipt["status"], "decision": receipt["decision"]}, sort_keys=True))


if __name__ == "__main__":
    main()
