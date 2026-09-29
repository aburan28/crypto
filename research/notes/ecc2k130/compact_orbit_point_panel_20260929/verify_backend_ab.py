#!/usr/bin/env python3
"""Check same-Q, same-work rho backend A/B output and full-process costs."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics
import traceback

from verify_panel import HERE, ROOT, check_fixture, rows, sha


def interval(ratios: list[float]) -> list[float]:
    logs = [math.log(value) for value in ratios]
    mean = statistics.mean(logs)
    half = 2.776445105 * statistics.stdev(logs) / math.sqrt(len(logs))
    return [math.exp(mean - half), math.exp(mean + half)]


def check(n: int, length: int, run_dir: Path) -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert {"n": n, "L": length} in frozen["backend_ab_cases"]
    for filename, expected in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == expected, filename
    spec = frozen["specs"][f"n{n}_L{length}_eval"]
    fixtures, curve, generator = check_fixture(spec)
    host = json.loads((run_dir / "host.json").read_text())
    assert "x86_64" in host["platform"]
    assert "pclmulqdq" in host["cpu_flags"]
    runs = json.loads((run_dir / "runs.json").read_text())
    assert len(runs) == 2 * frozen["backend_ab_repetitions"]
    measurements = {}
    for block in range(frozen["backend_ab_repetitions"]):
        expected = ("baseline", "candidate") if block % 2 == 0 else ("candidate", "baseline")
        pair = runs[2 * block:2 * block + 2]
        assert tuple(item["arm"] for item in pair) == expected
        outputs = {}
        measurements[block] = {}
        for item in pair:
            arm = item["arm"]
            assert item["block"] == block
            assert item["exit_code"] == 0 and not item["timeout"]
            stdout = run_dir / item["stdout"]
            stderr = run_dir / item["stderr"]
            assert sha(stdout) == item["stdout_sha256"]
            assert sha(stderr) == item["stderr_sha256"]
            assert item["environment"]["KIC_RHO_BATCH_CORPUS"] == spec["corpus"]
            assert item["environment"]["KIC_RHO_DP_BITS"] == str(frozen["rho_dp_bits"])
            assert Path(item["environment"]["KIC_RHO_POINT_INPUT"]).name == Path(spec["points_file"]).name
            assert not any(".fixture.jsonl" in value
                           for value in item["command"] + list(item["environment"].values()))
            data = rows(stdout)
            assert len(data) == length + 1
            summary = data[-1]
            assert summary["kind"] == "rho_ks_batch_summary"
            assert summary["target_source"] == "public_point_jsonl"
            assert summary["quotient_mode"] == "signed_frobenius"
            assert summary["fixtures"] == length and summary["all_verified"] is True
            assert summary["corpus"] == spec["corpus"]
            if arm == "candidate":
                assert summary["field_product_backend"] == "pclmulqdq"
                assert summary["inversion_backend"] == "itoh_tsujii"
            else:
                assert "field_product_backend" not in summary
            for index, (record, fixture) in enumerate(zip(data[:-1], fixtures)):
                assert record["kind"] == "rho_ks_batch_fixture"
                assert record["fixture_index"] == index
                assert record["published_fixture_scalar"] is None
                assert record["target_source"] == "public_point_jsonl"
                assert record["published_q"] == fixture["published_q"]
                assert record["recovered_fixture_scalar"] == fixture["published_fixture_scalar"]
                # check_fixture already independently proved this scalar times G is Q.
            outputs[arm] = data
            measurements[block][arm] = {
                "wall_seconds": item["wall_seconds"],
                "cpu_seconds": item["user_seconds"] + item["sys_seconds"],
                "rss_kib": item["max_rss_kib_linux"],
            }
        baseline, candidate = outputs["baseline"], outputs["candidate"]
        for old, new in zip(baseline[:-1], candidate[:-1]):
            for key in ("published_q", "recovered_fixture_scalar", "walk_steps",
                        "walks", "fruitless_two_cycles", "capped_walks",
                        "wasted_merges", "table_entries_before", "table_entries_after",
                        "solved_via_target", "cross_target_solve"):
                assert old[key] == new[key], (block, key)
        for key in ("total_walk_steps", "table_entries", "cross_target_solves", "charges"):
            assert baseline[-1][key] == candidate[-1][key], (block, key)
    wall = [measurements[b]["candidate"]["wall_seconds"] /
            measurements[b]["baseline"]["wall_seconds"] for b in measurements]
    cpu = [measurements[b]["candidate"]["cpu_seconds"] /
           measurements[b]["baseline"]["cpu_seconds"] for b in measurements]
    return {
        "status": "PASS", "n": n, "L": length, "blocks": len(measurements),
        "same_Q_scalar_walk_and_charges": True,
        "targets_verified_per_arm": length * len(measurements),
        "candidate_over_baseline_wall_median": statistics.median(wall),
        "candidate_over_baseline_wall_range": [min(wall), max(wall)],
        "candidate_over_baseline_wall_95pct_log_t_interval": interval(wall),
        "candidate_over_baseline_cpu_median": statistics.median(cpu),
        "candidate_over_baseline_cpu_range": [min(cpu), max(cpu)],
        "candidate_over_baseline_cpu_95pct_log_t_interval": interval(cpu),
        "measurements": measurements,
        "baseline_binary_sha256": host["baseline_binary_sha256"],
        "candidate_binary_sha256": host["candidate_binary_sha256"],
        "host_cpu_model": host["cpu_model"],
        "host_cpu_flags_sha256": __import__("hashlib").sha256(
            " ".join(host["cpu_flags"]).encode()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--L", type=int, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    try:
        result = check(args.n, args.L, args.run_dir.resolve())
    except BaseException as error:
        result = {"status": "FAIL", "error_type": type(error).__name__,
                  "error": str(error), "traceback": traceback.format_exc()}
        args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "measurements"}, sort_keys=True))


if __name__ == "__main__":
    main()
