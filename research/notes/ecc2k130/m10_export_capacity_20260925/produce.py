#!/usr/bin/env python3
"""Held n131 m10 complete-chain DAG/CNF size producer; never invokes a solver."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

from basis_verify import replay as replay_basis
from capacity import NodeCapExceeded, build_chain, exact_dimacs_size, slot_bases_n131

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clock():
    rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return (time.monotonic(), time.process_time(),
            rss if sys.platform == "darwin" else rss * 1024)


def timeout(_signum, _frame):
    raise TimeoutError("frozen per-arm wall cap reached")


def check_hashes(frozen):
    for rel, expected in frozen["input_sha256"].items():
        assert sha(ROOT / rel) == expected, rel
    for name, expected in frozen["source_sha256"].items():
        assert sha(HERE / name) == expected, name


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=("balanced", "unequal"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    spec = json.loads((HERE / "INPUT.json").read_text())
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    caps = spec["caps"]
    started = clock()
    phase = "release_gate"
    result = {"schema": "ecc2k130-m10-capacity-arm-result-v1",
              "status": "STOP", "arm": args.arm, "phase": phase,
              "freeze_sha256": sha(HERE / "FROZEN.json"),
              "release_main_head": frozen["release_main_head"],
              "scope": "complete-chain representation size only; no solver or PDP"}
    try:
        if frozen["release_main_head"] is None:
            raise RuntimeError("NOT_ADMITTED: #802, #804 and #784 must merge, then re-freeze")
        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(caps["per_arm_external_wall_seconds"])
        try:
            resource.setrlimit(resource.RLIMIT_AS,
                               (caps["per_arm_address_space_and_rss_bytes"],
                                caps["per_arm_address_space_and_rss_bytes"]))
        except (AttributeError, OSError, ValueError) as exc:
            raise RuntimeError("frozen memory cap could not be enforced") from exc
        phase = "frozen_hashes_and_basis"
        check_hashes(frozen)
        basis = replay_basis()
        assert basis["status"] == "PASS"
        result["basis_replay"] = basis["arms"][args.arm]
        source_name = f"{args.arm}_m10_result.json"
        source = json.loads((HERE / "inputs" / source_name).read_text())
        assert source["beta"] == spec["normal_beta"]
        assert source["field_poly"] == int(spec["field_modulus_hex"], 16)
        assert source["q"] == 1 << spec["field_degree"]
        expected_dims = next(row["dimensions"] for row in spec["arms"]
                             if row["name"] == args.arm)
        if args.arm == "balanced":
            assert source["arm"] == {"d": 13, "m": 10}
            assert source["physical_f0_points"] == 7977
            assert source["nonzero_signed_columns"] == 3988
        else:
            assert source["slot_dimensions"] == expected_dims
            assert source["normalized_global_nonzero_signed_columns"] == 8062
            assert source["ordered_physical_tuples"] == 2108900315408742840629516059380574909125
        bases = slot_bases_n131(args.arm)
        assert list(map(len, bases)) == expected_dims
        result["source_result_sha256"] = sha(HERE / "inputs" / source_name)
        result["dimensions"] = expected_dims
        phase = "complete_chain_dag"
        progress_file = args.out / "progress.jsonl"
        def checkpoint(row):
            with progress_file.open("a") as stream:
                stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
                stream.flush()
                os.fsync(stream.fileno())
        try:
            chain = build_chain(spec["field_degree"], int(spec["field_modulus_hex"], 16),
                                bases, caps["dag_nodes"], checkpoint=checkpoint)
        except NodeCapExceeded as exc:
            result["status"] = "CENSORED_DAG_NODE_CAP"
            result["complete_chain"] = False
            result["cap_stage"] = exc.stage
            result["partial_dag_counts"] = exc.counts
            result["partial_prefix_sha256"] = exc.prefix_sha256
        else:
            result["complete_chain"] = True
            result["dag_counts"] = chain.counts()
            result["dag_prefix_sha256"] = chain.dag.prefix_sha256()
            result["local_relation_counts"] = chain.local_relation_counts
            result["edge_count"] = len(chain.edges)
            result["edge_labels"] = [list(edge) for edge in chain.edges]
            result["primary_inputs_expected"] = (
                sum(expected_dims) + 10 * 131 + 9 * (3 * 131 + 1))
            assert chain.counts()["variables"] == result["primary_inputs_expected"]
            phase = "exact_dimacs_size"
            result["cnf_generic"] = exact_dimacs_size(chain, fixed_target_O=False)
            result["cnf_target_O"] = exact_dimacs_size(chain, fixed_target_O=True)
            assert (result["cnf_target_O"]["clauses"] ==
                    result["cnf_generic"]["clauses"] + 263)
            if result["cnf_target_O"]["bytes"] > caps["dimacs_bytes"]:
                result["status"] = "CENSORED_CNF_BYTE_CAP"
            else:
                result["status"] = "WITHIN_FROZEN_REPRESENTATION_CAPS"
            result["node_cap_ratio"] = chain.counts()["total_nodes"] / caps["dag_nodes"]
            result["cnf_byte_cap_ratio"] = result["cnf_target_O"]["bytes"] / caps["dimacs_bytes"]
    except BaseException as exc:
        result["status"] = "STOP"
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
        result["traceback"] = traceback.format_exc(limit=12)
    finally:
        end = clock()
        result["phase"] = phase
        result["resources"] = {"wall_seconds": end[0] - started[0],
                               "cpu_seconds": end[1] - started[1],
                               "peak_rss_bytes": end[2]}
        if end[0] - started[0] > caps["per_arm_external_wall_seconds"]:
            result["status"] = "STOP"
            result.setdefault("error", "wall cap exceeded")
        if end[2] > caps["per_arm_address_space_and_rss_bytes"]:
            result["status"] = "STOP"
            result.setdefault("error", "memory cap exceeded")
        (args.out / "result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
        signal.alarm(0)
    return 0 if result["status"] in ("WITHIN_FROZEN_REPRESENTATION_CAPS",
                                      "CENSORED_CNF_BYTE_CAP",
                                      "CENSORED_DAG_NODE_CAP") else 1


if __name__ == "__main__":
    sys.exit(main())
