#!/usr/bin/env python3
"""Freeze and run one public-point IC/rho pair under the strict isolated service.

Run ``manifest`` on the prepared physical Linux host after compiling both
online examples from one committed source tree. The manifest is consumed by
``cryptanalysis/scripts/isolated_bench.py``. Each ``arm`` invocation performs
target-independent setup inside its producer, reads that producer's internal
online interval, independently replays the answer, and emits the service's
``online_ms=... verified=1`` line. The service retains the full stdout and
host/noise receipts for both arms.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUTOLAB = ROOT / "research/sat_factor_base_review_20260908/autolab"
sys.path.insert(0, str(AUTOLAB))
import single_target_panel as panel  # noqa: E402

BEAT_ID = "koblitz.compact_orbit.n61_single_target"
PROTOCOL = AUTOLAB / "protocol.json"
IC_PHASES = ("target_query_ms", "target_pdp_ms", "target_relation_check_ms",
             "target_descent_ms", "target_recovery_check_ms")


def beat_record() -> dict:
    return json.loads(PROTOCOL.read_text())["beats"][BEAT_ID]


def point_hash(n: int, a: int, point: list[int]) -> str:
    record = {"n": n, "a": a, "target": point}
    return hashlib.sha256(json.dumps(record, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def rows(text: str) -> list[dict]:
    return [json.loads(line) for line in text.splitlines() if line.startswith("{")]


def require_one(records: list[dict], kind: str) -> dict:
    found = [record for record in records if record.get("kind") == kind]
    if len(found) != 1:
        raise ValueError(f"expected one {kind} record, found {len(found)}")
    return found[0]


def equal_sum(values: list[float], total: float) -> bool:
    return math.isclose(math.fsum(values), total, rel_tol=1e-6, abs_tol=1e-3)


def run_arm(args: argparse.Namespace) -> None:
    beat = beat_record()
    n, a = int(beat["n"]), int(beat["a"])
    target = [args.x, args.y]
    fixture = beat["curve_fixture"]
    curve = panel.cached_curve(fixture)
    target_point = curve.decode(target)
    if target_point is None or curve.mul(target_point, curve.r) is not None:
        raise ValueError("target is not a point in the declared subgroup")
    env = os.environ.copy()
    env["RAYON_NUM_THREADS"] = "1"
    with tempfile.TemporaryDirectory(prefix="ic-one-target-") as directory:
        temp = Path(directory)
        if args.arm == "ic":
            input_path = temp / "target.jsonl"
            output_path = temp / "ic.jsonl"
            base_path = temp / "base.jsonl"
            input_path.write_text(json.dumps(target) + "\n")
            env["KIC_DUMP_BASE"] = str(base_path)
            command = [str(args.binary), f"construct:{n}:{a}:{args.k}", str(input_path),
                       str(beat["ic_rank_seed"]), str(output_path)]
        else:
            env.update({key: str(value) for key, value in beat["panel_producers"]["rho"]["env"].items()})
            env["KIC_RHO_TARGET_POINT"] = f"{args.x},{args.y}"
            command = [str(args.binary), str(n), str(a), "signed_frobenius", "1",
                       str(int(beat["rho_seed_base"]) + args.target_seed)]
        result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True)
        if result.returncode:
            print(result.stdout, file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            raise RuntimeError(f"{args.arm} producer exited {result.returncode}")

        if args.arm == "ic":
            record = require_one(rows(output_path.read_text()), "compact_orbit_dlp_target")
            summary = require_one(rows(result.stdout), "compact_orbit_dlp_summary")
            base = rows(base_path.read_text())[0]
            indices = record["point_indices"]
            points = [base["factor_base_point_coordinates"][index] for index in indices]
            relation_sum = None
            for point in points:
                relation_sum = curve.add(relation_sum, curve.decode(point))
            valid = (record["group_verified"] is True and record["target"] == target
                     and summary["rank"] == summary["orbit_columns"]
                     and relation_sum == target_point
                     and [point[0] for point in points] == record["x_codes"]
                     and equal_sum([float(record[key]) for key in IC_PHASES], float(record["online_ms"])))
            scalar = int(record["recovered_scalar"])
            certificate = {"base_hash": summary["base_hash"], "point_indices": indices,
                           "relation_points": points, "x_codes": record["x_codes"]}
        else:
            record = require_one(rows(result.stdout), "rho_ks_batch_fixture")
            summary = require_one(rows(result.stdout), "rho_ks_batch_summary")
            valid = (record["verified"] is True and record["published_q"] == target
                     and record["target_source"] == "public_point"
                     and record["table_entries_before"] == 0
                     and equal_sum([float(record["walk_and_collision_ms"]),
                                    float(record["recovery_check_ms"])], float(record["online_ms"])))
            scalar = int(record["recovered_fixture_scalar"])
            certificate = {"walk_steps": record["walk_steps"], "rho_policy": {
                "rung": summary["rung"], "lanes": summary["lanes"],
                "distinguished_point_bits": summary["dp_bits"],
                "automorphism_size": summary["automorphism_size"]}}
        valid = valid and 0 <= scalar < curve.r and curve.mul(curve.g, scalar) == target_point
        if not valid:
            raise ValueError(f"{args.arm} producer record failed independent replay")
        online_ms = float(record["online_ms"])
        if not math.isfinite(online_ms) or online_ms <= 0:
            raise ValueError("invalid online interval")
        print("producer_json=" + json.dumps(record, sort_keys=True, separators=(",", ":")))
        print("certificate_json=" + json.dumps(certificate, sort_keys=True, separators=(",", ":")))
        print("summary_json=" + json.dumps(summary, sort_keys=True, separators=(",", ":")))
        print(f"online_ms={online_ms:.9f} verified=1 target_hash={point_hash(n, a, target)} "
              f"scalar={scalar} arm={args.arm} n={n} a={a}")


def make_manifest(args: argparse.Namespace) -> None:
    beat = beat_record()
    n, a = int(beat["n"]), int(beat["a"])
    curve = panel.cached_curve(beat["curve_fixture"])
    target = list(panel.public_point(curve, args.domain, "isolated", args.target_seed))
    script = Path(__file__).resolve()
    common = ["--k", str(args.k), "--x", str(target[0]), "--y", str(target[1]),
              "--target-seed", str(args.target_seed)]
    sources = set(beat["ic_sources"] + beat["rho_sources"])
    sources.update(("Cargo.lock", "Cargo.toml",
                    "research/sat_factor_base_review_20260908/autolab/protocol.json",
                    "research/sat_factor_base_review_20260908/autolab/single_target_panel.py",
                    "research/ic_candidate_tournament_20260915/oracle.py",
                    "research/ic_candidate_tournament_20260915/identity.py"))
    artifacts = [str(ROOT / name) for name in sorted(sources)]
    artifacts.extend((str(script), str(args.ic_binary), str(args.rho_binary)))
    manifest = {
        "schema": 1, "name": "koblitz-n61-one-public-target-ic-vs-strong-rho",
        "workdir": str(ROOT), "isolation": {"cgroup": str(args.cgroup), "cpus": args.cpus,
                                          "execution_cpu": args.execution_cpu,
                                          "mem_nodes": args.mem_nodes},
        "artifacts": artifacts, "timeout_s": args.timeout_s,
        "repetitions": args.repetitions,
        "measurement_boundary": "producer online_ms: first target query to scalar replay for IC; first target-dependent rho walk operation to scalar replay for rho; reusable IC setup and point generation excluded",
        "pair_fields": ["n", "a", "target_hash", "scalar"],
        "target_law": {"domain": args.domain, "role": "isolated", "seed": args.target_seed,
                       "public_point": target, "target_hash": point_hash(n, a, target)},
        "candidate_configuration": {"K": args.k, "rank_seed": beat["ic_rank_seed"],
                                    "rho_seed": int(beat["rho_seed_base"]) + args.target_seed},
        "cases": [{"id": "n61-public-target-0",
                   "candidate": [str(script), "arm", "--arm", "ic", "--binary", str(args.ic_binary), *common],
                   "reference": [str(script), "arm", "--arm", "rho", "--binary", str(args.rho_binary), *common]}],
    }
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(args.output), "target": target,
                      "target_hash": manifest["target_law"]["target_hash"]}, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    arm = sub.add_parser("arm")
    arm.add_argument("--arm", choices=("ic", "rho"), required=True)
    arm.add_argument("--binary", type=Path, required=True)
    arm.add_argument("--k", type=int, required=True)
    arm.add_argument("--x", type=int, required=True)
    arm.add_argument("--y", type=int, required=True)
    arm.add_argument("--target-seed", type=int, required=True)
    frozen = sub.add_parser("manifest")
    frozen.add_argument("--ic-binary", type=Path, required=True)
    frozen.add_argument("--rho-binary", type=Path, required=True)
    frozen.add_argument("--cgroup", type=Path, required=True)
    frozen.add_argument("--cpus", required=True)
    frozen.add_argument("--execution-cpu", type=int, required=True)
    frozen.add_argument("--mem-nodes", required=True)
    frozen.add_argument("--output", type=Path, required=True)
    frozen.add_argument("--k", type=int, default=400)
    frozen.add_argument("--domain", default="kic-n61-isolated-one-target-20261009-v1")
    frozen.add_argument("--target-seed", type=int, default=0)
    frozen.add_argument("--repetitions", type=int, default=1)
    frozen.add_argument("--timeout-s", type=int, default=86400)
    args = parser.parse_args()
    if args.action == "arm":
        run_arm(args)
    else:
        make_manifest(args)


if __name__ == "__main__":
    main()
