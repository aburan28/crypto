#!/usr/bin/env python3
"""Apply frozen three-block same-Q wall gates to an archived panel."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"


def read_run(n: int, name: str) -> tuple[dict, dict]:
    run = EVIDENCE / "runs" / f"n{n}-{name}"
    return (json.loads((run / "receipt.json").read_bytes()),
            json.loads((run / "independent_validation.json").read_bytes()))


def evaluate(n: int) -> dict:
    panel = json.loads((EVIDENCE / f"panel_n{n}.json").read_bytes())
    attempts = {row["name"]: row for row in panel["attempts"]}
    train = attempts.get("train")
    result = {"n": n, "panel_classification": panel["classification"],
              "panel_sha256": sha256((EVIDENCE / f"panel_n{n}.json").read_bytes()).hexdigest(),
              "training": None, "levels": {}}
    if train is None or train["status"] != "COMPLETE":
        for length in (8, 32):
            result["levels"][str(length)] = {"classification": "CENSORED_TRAINING"}
        return result
    train_receipt, train_report = read_run(n, "train")
    training_child_ms = train_receipt["wall_ms"]
    training_replay_ms = train["replay_wall_ms"]
    result["training"] = {"child_wall_ms": training_child_ms,
                          "independent_rank_replay_wall_ms": training_replay_ms,
                          "rank": train_report["rank"], "R": train_report["R"],
                          "relations": train_report["relations"],
                          "root_occupancy": train_report["root_occupancy"],
                          "total_s3_calls": train_report["total_s3_calls"],
                          "peak_rss_bytes": train_receipt["peak_rss_bytes"]}
    for length in (8, 32):
        block_rows = []
        for block in range(3):
            pair = {}
            for mode in ("compact", "rho"):
                name = f"L{length}-b{block}-{mode}"
                attempt = attempts.get(name)
                pair[mode] = {"status": attempt["status"] if attempt else "UNRUN"}
                if attempt and attempt["status"] == "COMPLETE":
                    receipt, report = read_run(n, name)
                    pair[mode].update({"child_wall_ms": receipt["wall_ms"],
                                       "peak_rss_bytes": receipt["peak_rss_bytes"],
                                       "report": report})
                    if mode == "compact":
                        pair[mode]["independent_recovery_wall_ms"] = attempt["replay_wall_ms"]
            row = {"block": block, "compact": pair["compact"], "rho": pair["rho"]}
            if all(pair[mode]["status"] == "COMPLETE" for mode in pair):
                ic_child = pair["compact"]["child_wall_ms"]
                ic_replay = pair["compact"]["independent_recovery_wall_ms"]
                rho_wall = pair["rho"]["child_wall_ms"]
                assert rho_wall > 0
                lower = training_child_ms + ic_child
                upper = training_child_ms + training_replay_ms + ic_child + ic_replay
                row["charged_wall_ms"] = {"ic_lower": lower, "ic_upper": upper,
                                           "rho_cold_batch": rho_wall,
                                           "ic_lower_over_rho": lower/rho_wall,
                                           "ic_upper_over_rho": upper/rho_wall}
            block_rows.append(row)
        complete = all("charged_wall_ms" in row for row in block_rows)
        if not complete:
            classification = "CENSORED_FIXED_STREAM"
        elif all(row["charged_wall_ms"]["ic_upper"] <
                 row["charged_wall_ms"]["rho_cold_batch"] for row in block_rows):
            classification = "FAVORABLE_FULLY_CHARGED_SAME_HOST_WALL_OBSERVATION"
        elif all(row["charged_wall_ms"]["ic_lower"] >
                   row["charged_wall_ms"]["rho_cold_batch"] for row in block_rows):
            classification = "NEGATIVE_FULLY_CHARGED_SAME_HOST_WALL_OBSERVATION"
        else:
            classification = "UNRESOLVED_FIXED_STREAM"
        level = {"classification": classification, "all_three_pairs_complete": complete,
                 "blocks": block_rows}
        if complete:
            lower_ratios = [row["charged_wall_ms"]["ic_lower_over_rho"] for row in block_rows]
            upper_ratios = [row["charged_wall_ms"]["ic_upper_over_rho"] for row in block_rows]
            level["paired_ratio_descriptive"] = {
                "lower_min_median_max": [min(lower_ratios), median(lower_ratios), max(lower_ratios)],
                "upper_min_median_max": [min(upper_ratios), median(upper_ratios), max(upper_ratios)],
                "inferential_95_percent_interval": None,
                "note": "Three fixed blocks cannot exclude parity with a distribution-free 95% sign/permutation interval"}
            level["three_block_portfolio_wall_ms"] = {
                "ic_lower_training_once": training_child_ms + sum(
                    row["compact"]["child_wall_ms"] for row in block_rows),
                "ic_upper_training_and_rank_once": training_child_ms + training_replay_ms + sum(
                    row["compact"]["child_wall_ms"] +
                    row["compact"]["independent_recovery_wall_ms"] for row in block_rows),
                "rho_three_cold_batches": sum(row["rho"]["child_wall_ms"] for row in block_rows)}
        result["levels"][str(length)] = level
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    assert EVIDENCE.is_dir()
    report = {"schema_version": "1.0", "classification": "TOY_RUNG_SHARED_LOG_CONTROL",
              "n131_transfer": None,
              "common_group_addition_equivalent_S": None,
              "arms": {str(n): evaluate(n) for n in (37, 41)}}
    data = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(data)
    print(json.dumps({n: report["arms"][str(n)]["levels"]["8"]["classification"]
                      for n in (37, 41)}, sort_keys=True))


if __name__ == "__main__":
    main()
