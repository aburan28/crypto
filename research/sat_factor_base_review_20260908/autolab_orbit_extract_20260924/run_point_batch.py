#!/usr/bin/env python3
"""Freeze and run point-only n=53 compact-orbit target probes."""

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RANK = HERE / "cold_batch_rank.py"
DOMAIN = b"ECC2K53-COMPACT-POINT-BATCH-20260924-v1/"


def load_rank():
    spec = importlib.util.spec_from_file_location("cold_batch_rank", RANK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run(args):
    rank = load_rank()
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    with gzip.open(rank.BASE_GZ, "rb") as stream:
        base = stream.read()
    header = json.loads(base)
    assert header["base_hash"] == rank.BASE_HASH and header["n"] == 53
    curve = rank.load_verifier().Curve(header)
    order = curve.order
    scalars = [1 + int.from_bytes(hashlib.sha256(DOMAIN + str(i).encode()).digest(), "big")
               % (order - 1) for i in range(args.targets)]
    assert len(set(scalars)) == args.targets
    training_schedule = set(rank.target_schedule(
        json.loads((args.training / "manifest.json").read_text())["targets"], order
    ))
    assert not training_schedule.intersection(scalars)
    points = [curve.scalar(curve.generator, scalar) for scalar in scalars]
    assert all(point is not None for point in points)
    point_bytes = "".join(json.dumps(point, separators=(",", ":")) + "\n"
                          for point in points).encode()
    (out / "base_header.jsonl").write_bytes(base)
    (out / "target_points.jsonl").write_bytes(point_bytes)
    exe = REPO / "target/release/examples/koblitz_s5_sat_instance"
    source = REPO / "examples/koblitz_s5_sat_instance.rs"
    command = [str(exe), "53", "0", "1", "10", "natural", "1", "2000", "1", "internal"]
    manifest = {
        "schema_version": "1.0", "task": "point_input_log_recovery_probe",
        "scope": "Public synthetic Q points; validator scalar labels withheld from producer",
        "domain": DOMAIN.decode(), "point_file_sha256": digest(point_bytes),
        "base_hash": header["base_hash"],
        "base_gzip_sha256": digest(rank.BASE_GZ.read_bytes()),
        "source_sha256": digest(source.read_bytes()), "exe_sha256": digest(exe.read_bytes()),
        "training_validation_sha256": digest((args.training / "validation.json").read_bytes()),
        "training_manifest_sha256": digest((args.training / "manifest.json").read_bytes()),
        "validation_scalars": scalars, "command": command,
        "stopping_rule": "One cold process on the first N precommitted points; count all failures",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    env = {key: val for key, val in os.environ.items() if not key.startswith("KIC_")}
    env.update({
        "KIC_ALGEBRA_ENCODING": "orbit_factorized",
        "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
        "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
        "KIC_ORBIT_REP_ENCODING": "one_hot",
        "KIC_FACTOR_BASE_JSONL": str(out / "base_header.jsonl"),
        "KIC_TASK_ID": "TASK-IC-COMPACT-POINT-RECOVERY-20260924",
        "KIC_ORBIT_TARGET_POINTS_JSONL": str(out / "target_points.jsonl"),
        "KIC_ORBIT_BATCH_ONLY": "1",
    })
    start = time.perf_counter()
    process = subprocess.run(command, cwd=REPO, env=env, capture_output=True,
                             text=True, timeout=args.timeout)
    wall_ms = (time.perf_counter() - start) * 1000
    peak_rss_raw = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    (out / "producer.stdout.jsonl").write_text(process.stdout)
    (out / "producer.stderr.txt").write_text(process.stderr)
    assert process.returncode == 0, (process.returncode, process.stderr[-1000:])
    batch = json.loads(process.stdout)["compact_orbit_point_batch"]
    assert batch["targets_requested"] == args.targets
    assert batch["targets_extracted"] + len(batch["failed_target_points"]) == args.targets
    receipt = {
        "schema_version": "1.0", "positive_process_wall_ms": wall_ms,
        "positive_process_peak_rss_raw": peak_rss_raw,
        "producer_stdout_sha256": digest(process.stdout.encode()),
        "producer_stderr_sha256": digest(process.stderr.encode()),
        "targets_requested": args.targets, "targets_extracted": batch["targets_extracted"],
        "negative_controls": [],
    }
    for name, point, expected in (
        ("off_curve", [0, 0], "point target must be on the curve"),
        ("wrong_subgroup", [0, 1], "point target must belong to the prime-order subgroup"),
        ("missing_compact_mode", list(points[0]),
         "compact batch target files require KIC_ORBIT_LAZY_RELATIVE_SUPPORT=1"),
    ):
        path = out / f"{name}.jsonl"
        path.write_text(json.dumps(point) + "\n")
        negative_env = env.copy()
        negative_env["KIC_ORBIT_TARGET_POINTS_JSONL"] = str(path)
        if name == "missing_compact_mode":
            negative_env.pop("KIC_ORBIT_LAZY_RELATIVE_SUPPORT")
        result = subprocess.run(command, cwd=REPO, env=negative_env, capture_output=True,
                                text=True, timeout=args.timeout)
        (out / f"{name}.stderr.txt").write_text(result.stderr)
        assert result.returncode != 0 and expected in result.stderr
        receipt["negative_controls"].append({
            "case": name, "point": point, "returncode": result.returncode,
            "expected_rejection": expected, "stderr_sha256": digest(result.stderr.encode()),
        })
    (out / "resource_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in (
        "positive_process_wall_ms", "positive_process_peak_rss_raw",
        "targets_requested", "targets_extracted")}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--targets", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=90)
    arguments = parser.parse_args()
    assert 1 <= arguments.targets <= 32
    run(arguments)
