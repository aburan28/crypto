#!/usr/bin/env python3
"""Second-machine arithmetic replay of an extracted hosted panel archive."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

from verify_panel import (ARMS, check_counts, check_fixture, check_target, freeze,
                          rows, sha, target_identity, verify_rank)


def identity_digest(values: tuple) -> str:
    return hashlib.sha256(json.dumps(values, separators=(",", ":"),
                                     sort_keys=True).encode()).hexdigest()


def replay(n: int, run_dir: Path) -> dict:
    frozen = freeze()
    spec = frozen["specs"][f"n{n}_L1024_eval"]
    known, curve, generator = check_fixture(spec)
    assert run_dir.name == f"n{n}_L1024"
    hosted = json.loads((run_dir / "receipt.json").read_text())
    assert hosted["status"] == "PASS" and hosted["n"] == n
    runs = json.loads((run_dir / "runs.json").read_text())
    assert len(runs) == 5 * len(ARMS)
    checks = []
    verified_labels: set[tuple[str, int]] = set()
    for block in range(5):
        identities = {}
        block_runs = [item for item in runs if item["block"] == block]
        assert set(item["arm"] for item in block_runs) == set(ARMS)
        for item in block_runs:
            arm = item["arm"]
            assert item["exit_code"] == 0 and not item["timeout"]
            stdout = run_dir / item["stdout"]
            stderr = run_dir / item["stderr"]
            assert sha(stdout) == item["stdout_sha256"]
            assert sha(stderr) == item["stderr_sha256"]
            assert not stderr.read_text().strip()
            if arm == "rho_normal":
                rho_rows = rows(stdout)
                assert len(rho_rows) == 1025
                summary = rho_rows[-1]
                assert summary["kind"] == "rho_ks_batch_summary"
                assert summary["target_source"] == "public_point_jsonl"
                assert summary["canonicalization_backend"] == "normal_basis"
                assert summary["all_verified"] and summary["fixtures"] == 1024
                for fixture, output in zip(known, rho_rows[:-1]):
                    assert output["kind"] == "rho_ks_batch_fixture"
                    assert output["published_q"] == fixture["published_q"]
                    assert output["recovered_fixture_scalar"] == fixture["published_fixture_scalar"]
                    assert curve.mul(output["recovered_fixture_scalar"], generator) == tuple(
                        output["published_q"])
                checks.append({"block": block, "arm": arm, "logs_verified": 1024,
                               "stdout_sha256": sha(stdout)})
                continue
            filtered = arm == "filter"
            base_path = run_dir / f"b{block}_{arm}.base.jsonl"
            rank_path = run_dir / f"b{block}_{arm}.rank.jsonl"
            target_path = run_dir / f"b{block}_{arm}.target.jsonl"
            rank = verify_rank(rank_path, base_path, stdout)
            assert rank["status"] == "PASS" and rank["rank"] == frozen["k_by_n"][str(n)]
            base, = rows(base_path)
            rank_rows = rows(rank_path)
            logs = rank_rows[-1]["logs"]
            targets = rows(target_path)
            summary, = rows(stdout)
            assert len(targets) == summary["targets_solved"] == 1024
            assert summary["targets_failed"] == 0
            assert summary["root_prefilter_policy"] == (
                "blocked_bloom_512_3hash" if filtered else "off")
            for phase in ("index", "rank", "target"):
                check_counts(summary[f"{phase}_s3_counts"], filtered)
            for fixture, output in zip(known, targets):
                check_target(output, fixture, base, logs, curve, generator, verified_labels)
                check_counts(output["s3_counts"], filtered)
                assert output["probes"] == output["s3_counts"]["root_keys_considered"]
            identities[arm] = (base, rank_rows,
                               tuple(target_identity(output) for output in targets))
            checks.append({"block": block, "arm": arm, "rank_verified": rank["rank"],
                           "logs_verified": 1024, "base_sha256": sha(base_path),
                           "rank_sha256": sha(rank_path),
                           "target_sha256": sha(target_path),
                           "identity_sha256": identity_digest(identities[arm])})
        assert identities["off_a"] == identities["filter"] == identities["off_b"]
    return {
        "status": "PASS", "n": n, "L": 1024,
        "hosted_receipt_sha256": sha(run_dir / "receipt.json"),
        "replay_host": platform.platform(), "python": sys.version,
        "ic_rank_runs_verified": 15, "ic_logs_verified": 15 * 1024,
        "rho_logs_verified": 5 * 1024,
        "same_base_rank_first_witness_and_scalar_all_blocks": True,
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(41, 53), required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a replay receipt"
    receipt = replay(args.n, args.run_dir.resolve())
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in (
        "status", "n", "ic_rank_runs_verified", "ic_logs_verified", "rho_logs_verified")},
        sort_keys=True))


if __name__ == "__main__":
    main()
