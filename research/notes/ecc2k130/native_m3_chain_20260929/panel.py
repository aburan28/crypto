#!/usr/bin/env python3
"""Frozen toy SAT and exact-leaf representation producer; one arm per child."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import resource
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

from chain import (EXACT_MODULUS, NodeCapExceeded, build_native_chain,
                   exact_dimacs_size, leaf_slot_bases)
from reference import (TOY_A, TOY_B, TOY_BASES, TOY_MODULUS, TOY_N,
                       branch, independent_leaf_bases, parse_cnf_and_model,
                       primary_model, decode_model, sha, toy_oracle)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
NOTES = HERE.parent
sys.path.insert(0, str(NOTES / "symbolic_dag_dimacs_gate_20260925"))
from export import write_cnf  # noqa: E402

EXACT = ROOT / "research/ecc2k130_factor_base_replication_20260925/exact_smoke.json"
EXACT_REPLAY = ROOT / "research/ecc2k130_factor_base_replication_20260925/exact_replay.json"
DOMAIN = "ecc2k130-native-m3-chain-20260929-v1"


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def child_usage() -> tuple[float, int]:
    row = resource.getrusage(resource.RUSAGE_CHILDREN)
    return row.ru_utime + row.ru_stime, row.ru_maxrss if sys.platform == "darwin" else row.ru_maxrss * 1024


def target_units(chain, point):
    o, x, y = point
    nodes = chain.roles["SUM"]
    values = [o] + [(x >> i) & 1 for i in range(chain.n)] + [(y >> i) & 1 for i in range(chain.n)]
    flat = [nodes[0], *nodes[1], *nodes[2]]
    assert len(flat) == len(values)
    return [node + 1 if value else -(node + 1)
            for node, value in zip(flat, values, strict=True)]


def produce_toy(out: Path) -> dict:
    oracle = toy_oracle()
    curve = oracle["curve"]
    chain = build_native_chain(TOY_N, TOY_MODULUS, TOY_A, TOY_B,
                               TOY_BASES, 20000)
    positive, generic_positive, negative = (oracle["positive"],
                                             oracle["generic_positive"],
                                             oracle["negative"])
    branch_counts = Counter()
    triple_rows = out / "triples.jsonl"
    with triple_rows.open("w") as stream:
        for triple in itertools.product(*oracle["factors"]):
            s2 = curve.add(triple[0], triple[1])
            total = curve.add(s2, triple[2])
            wrong = next(point for point in oracle["points"] if point != total)
            slopes = (curve.slope(triple[0], triple[1]), curve.slope(s2, triple[2]))
            masks = [int(point[1] == basis[0])
                     for point, basis in zip(triple, TOY_BASES, strict=True)]
            correct_model = primary_model(chain, triple, s2, total, slopes, masks)
            wrong_model = primary_model(chain, triple, s2, wrong, slopes, masks)
            accepted = chain.dag.evaluate(correct_model, chain.output)
            rejected = not chain.dag.evaluate(wrong_model, chain.output)
            kinds = (branch(triple[0], triple[1]), branch(s2, triple[2]))
            branch_counts[f"edge0_{kinds[0]}"] += 1
            branch_counts[f"edge1_{kinds[1]}"] += 1
            row = {"factors": [list(point) for point in triple],
                   "masks": masks, "s2": list(s2), "target": list(total),
                   "wrong_target": list(wrong), "slopes": list(slopes),
                   "branches": list(kinds), "accepted": accepted,
                   "wrong_rejected": rejected}
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
            assert accepted and rejected, row
    assert sum(count for key, count in branch_counts.items() if key.startswith("edge0_")) == 27
    assert branch_counts["edge0_generic"] == branch_counts["edge1_generic"] == 24
    triple, s2 = oracle["support"][positive]
    slopes = (curve.slope(triple[0], triple[1]), curve.slope(s2, triple[2]))
    masks = [int(point[1] == basis[0]) for point, basis in zip(triple, TOY_BASES, strict=True)]
    witness = primary_model(chain, triple, s2, positive, slopes, masks)
    assert chain.dag.evaluate(witness, chain.output)
    assert all(chain.dag.evaluate(witness, node) for node in chain.edge_outputs)
    solver = shutil.which("cryptominisat5")
    if not solver:
        raise RuntimeError("CryptoMiniSat binary unavailable")
    solver_path = Path(solver).resolve()
    version = subprocess.run([str(solver_path), "--version"], capture_output=True,
                             text=True, timeout=10)
    assert version.returncode == 0 and "CryptoMiniSat version 5." in version.stdout
    (out / "solver.version.txt").write_text(version.stdout)
    runs = []
    for label, target in (("positive", positive),
                          ("generic_positive", generic_positive),
                          ("negative", negative)):
        cnf = out / f"{label}.cnf"
        exported = write_cnf(chain, cnf, units=target_units(chain, target),
                             byte_cap=16 * 1024 * 1024)
        assert sha(cnf) == exported["sha256"]
        command = [str(solver_path), "--verb", "0", "--threads", "1", str(cnf)]
        started = time.monotonic()
        cpu_before, _ = child_usage()
        try:
            job = subprocess.run(command, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired as exc:
            (out / f"{label}.stdout.txt").write_text(
                exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or ""))
            (out / f"{label}.stderr.txt").write_text(
                exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or ""))
            raise TimeoutError(f"{label} CryptoMiniSat exceeded 120 seconds") from exc
        wall = time.monotonic() - started
        child_cpu, child_rss = child_usage()
        stdout, stderr = out / f"{label}.stdout.txt", out / f"{label}.stderr.txt"
        stdout.write_text(job.stdout)
        stderr.write_text(job.stderr)
        status, bits, variables, clauses = parse_cnf_and_model(cnf, job.stdout, job.returncode)
        assert (variables, clauses) == (exported["variables"], exported["clauses"])
        assert status == ("UNSAT" if label == "negative" else "SAT")
        lift = decode_model(chain, bits, curve, target) if bits is not None else None
        assert (lift is not None) == (label != "negative")
        if label == "generic_positive":
            lifted_factors = [tuple(point) for point in lift["factors"]]
            assert branch(lifted_factors[0], lifted_factors[1]) == "generic"
            assert branch(tuple(lift["s2"]), lifted_factors[2]) == "generic"
        runs.append({"label": label, "target": list(target), "status": status,
                     "classification": ("SAT_LIFTED" if status == "SAT" else
                                        "SOLVER_UNSAT_ORACLE_CONFIRMED_TOY"),
                     "command": command, "exit_code": job.returncode,
                     "wall_seconds": wall, "child_cpu_seconds": child_cpu - cpu_before,
                     "child_peak_rss_bytes": child_rss,
                     "cnf": exported, "stdout_sha256": sha(stdout),
                     "stderr_sha256": sha(stderr), "lift": lift})
    return {"arm": "toy", "decision": "PASS", "field_degree": TOY_N,
            "modulus": hex(TOY_MODULUS), "curve_a": TOY_A, "curve_b": TOY_B,
            "bases": [list(slot) for slot in TOY_BASES],
            "factor_sizes": list(map(len, oracle["factors"])),
            "factor_triples": oracle["triples"],
            "distinct_supported_targets": len(oracle["support"]),
            "rational_curve_points": len(oracle["points"]),
            "positive": list(positive), "generic_positive": list(generic_positive),
            "negative": list(negative),
            "triple_rows": oracle["triples"],
            "triple_rows_sha256": sha(triple_rows),
            "branch_counts": dict(sorted(branch_counts.items())),
            "dag": chain.counts(), "prefix_sha256": chain.dag.prefix_sha256(),
            "solver_binary_sha256": sha(solver_path),
            "solver_version_sha256": sha(out / "solver.version.txt"),
            "runs": runs}


def produce_leaf(out: Path, arm: str) -> dict:
    archived = json.loads(EXACT.read_text())
    replay = json.loads(EXACT_REPLAY.read_text())
    assert archived["status"] == replay["status"] == "PASS"
    assert replay["input_sha256"] == sha(EXACT)
    line = (1, 0) if arm == "leaf10" else (1, 4)
    record = next(item for item in archived["lines"] if tuple(item["line"]) == line)
    b = int(record["codomain_b"], 16)
    bases = leaf_slot_bases()
    assert bases == independent_leaf_bases()
    progress = out / "progress.jsonl"
    def checkpoint(value):
        with progress.open("a") as stream:
            stream.write(json.dumps(value, sort_keys=True) + "\n")
            stream.flush()
    try:
        chain = build_native_chain(131, EXACT_MODULUS, 0, b, bases,
                                   2_000_000, checkpoint=checkpoint)
    except NodeCapExceeded as exc:
        return {"arm": arm, "decision": "CAPACITY_CENSORED",
                "censor": "DAG_NODE_CAP", "stage": exc.stage,
                "partial_counts": exc.counts,
                "partial_prefix_sha256": exc.prefix_sha256,
                "line": list(line), "curve_b": hex(b),
                "dimensions": list(map(len, bases)), "basis_rank": 131,
                "exact_smoke_sha256": sha(EXACT),
                "exact_replay_sha256": sha(EXACT_REPLAY)}
    generic = exact_dimacs_size(chain, fixed_target_O=False)
    target_o = exact_dimacs_size(chain, fixed_target_O=True)
    decision = ("PASS" if target_o["bytes"] <= 256 * 1024 * 1024
                else "CAPACITY_CENSORED")
    return {"arm": arm, "decision": decision,
            "censor": None if decision == "PASS" else "CNF_BYTE_CAP",
            "line": list(line), "curve_b": hex(b), "dimensions": list(map(len, bases)),
            "basis_rank": 131, "edge_count": len(chain.edges),
            "dag": chain.counts(), "local_edge_dag": chain.local_relation_counts,
            "prefix_sha256": chain.dag.prefix_sha256(),
            "cnf_generic": generic, "cnf_target_O": target_o,
            "exact_smoke_sha256": sha(EXACT),
            "exact_replay_sha256": sha(EXACT_REPLAY)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=("toy", "leaf10", "leaf14"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started, cpu_started = time.monotonic(), time.process_time()
    result = {"domain": DOMAIN, "arm": args.arm, "decision": "STOP"}
    try:
        result.update(produce_toy(args.out) if args.arm == "toy" else
                      produce_leaf(args.out, args.arm))
    except BaseException as exc:
        result.update({"decision": "STOP", "error_type": type(exc).__name__,
                       "error": str(exc), "traceback": traceback.format_exc(limit=10)})
    result.update({"wall_seconds": time.monotonic() - started,
                   "cpu_seconds": time.process_time() - cpu_started,
                   "peak_rss_bytes": rss_bytes()})
    (args.out / "result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"arm": args.arm, "decision": result["decision"],
                      "wall_seconds": result["wall_seconds"],
                      "peak_rss_bytes": result["peak_rss_bytes"]}, sort_keys=True))
    return 0 if result["decision"] in ("PASS", "CAPACITY_CENSORED") else 1


if __name__ == "__main__":
    raise SystemExit(main())
