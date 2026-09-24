#!/usr/bin/env python3
"""Independently check natural seeds 9..16 after fixing compact-orbit traversal."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from replay import Curve, check_witness

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "examples/koblitz_s5_sat_instance.rs"
EXE = ROOT / "target/release/examples/koblitz_s5_sat_instance"
BASE = HERE / "base_header.jsonl.gz"
BASE_HASH = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"
SEEDS = tuple(range(9, 17))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--variant", choices=("deterministic", "unsorted"), required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    args.out.mkdir(parents=True)
    source_sha, exe_sha = sha(SOURCE), sha(EXE)
    with gzip.open(BASE, "rb") as stream:
        header_bytes = stream.readline()
        assert not stream.read()
    header = json.loads(header_bytes)
    assert header["base_hash"] == BASE_HASH
    curve = Curve(header)
    assert curve.scalar(curve.generator, curve.order) is None
    points = [tuple(point) for point in header["factor_base_point_coordinates"]]
    by_x = {}
    for point in points:
        assert curve.on_curve(point)
        by_x.setdefault(point[0], []).append(point)
    assert len(points) == 23320 and len(by_x) == 11660
    fixture = args.out / "base_header.jsonl"
    fixture.write_bytes(header_bytes)
    env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
    env.update({
        "KIC_ALGEBRA_ENCODING": "orbit_factorized",
        "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
        "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
        "KIC_ORBIT_REP_ENCODING": "one_hot",
        "KIC_FACTOR_BASE_JSONL": str(fixture),
        "KIC_TASK_ID": "TASK-IC-COMPACT-ORBIT-HOLDOUT-20260924",
    })
    rows = []
    for seed in SEEDS:
        command = [str(EXE), "53", "0", "1", "10", "natural", str(seed), "2000", "1", "internal"]
        started = time.perf_counter()
        p = subprocess.run(command, capture_output=True, text=True, env=env, timeout=30)
        wall = time.perf_counter() - started
        assert p.returncode == 0, (seed, p.stderr[-1000:])
        assert len(p.stdout.splitlines()) == 1
        record = json.loads(p.stdout)
        assert record["factor_base_input_hash"] == BASE_HASH
        assert record["target_class"] == "natural" and record["seed"] == seed
        assert record["pair_table_entries"] == record["pair_selector_variables"] == 0
        assert record["invalid_group_lifts"] == 0
        label = record["published_scalar_validator_label"]
        assert isinstance(label, int) and 1 <= label < curve.order
        target = curve.scalar(curve.generator, label)
        assert target is not None and curve.on_curve(target)
        extraction = record["compact_orbit_extraction"]
        assert extraction["enabled"] and extraction["pair_table_entries"] == extraction["edge_selectors"] == 0
        lift = None
        if extraction["group_valid"]:
            assert record["decomposition_verdict"] == "SAT"
            lift = check_witness(curve, by_x, target, extraction["x_codes"], extraction["pinned_intermediates"])
            for point in lift:
                assert curve.scalar(tuple(point), curve.order) is None
        rows.append({
            "seed": seed, "target": list(target), "published_scalar": label,
            "group_valid": extraction["group_valid"], "independent_lift": lift,
            "extract_ms": extraction["extract_ms"], "trials": extraction["trials"],
            "index_entries": extraction["index_entries"],
            "process_wall_seconds": wall, "raw_producer": record,
            "stderr": p.stderr,
        })
        print(json.dumps({"seed": seed, "group_valid": extraction["group_valid"], "wall_seconds": wall}), flush=True)
    fixture.unlink()
    assert sha(SOURCE) == source_sha and sha(EXE) == exe_sha
    report = {
        "schema": "compact-orbit-natural-holdout-v1",
        "status": "PASS",
        "variant": args.variant,
        "seeds": list(SEEDS),
        "host": platform.platform(),
        "python": sys.version,
        "verifier_sha256": sha(__file__),
        "independent_arithmetic_sha256": sha(HERE / "replay.py"),
        "source_sha256": source_sha,
        "executable_sha256": exe_sha,
        "base_gzip_sha256": sha(BASE),
        "base_header_sha256": hashlib.sha256(header_bytes).hexdigest(),
        "base_hash": BASE_HASH,
        "base_points_checked_on_curve": len(points),
        "verified_hits": sum(row["group_valid"] for row in rows),
        "rows": rows,
        "claim_boundary": "Held-out positive relation check and stage diagnostics only; no exhaustive negative certificate, full-rank matrix, scalar recovery, or rho comparison.",
    }
    (args.out / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "hits": report["verified_hits"], "seeds": list(SEEDS)}))


if __name__ == "__main__":
    main()
