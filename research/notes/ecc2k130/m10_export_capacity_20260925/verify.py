#!/usr/bin/env python3
"""Independent m10 capacity receipt replay; no SAT, PDP or DLP work."""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import signal
import sys
from pathlib import Path

from basis_verify import replay as replay_basis
from capacity import NodeCapExceeded, build_chain, slot_bases_n131

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_dimacs_size(chain, target_O: bool) -> dict:
    """Byte arithmetic independent of capacity.exact_dimacs_size and export.py."""
    d = chain.dag
    count = d.counts()
    variables = count["total_nodes"]
    unit_ids = []
    if target_O:
        o, x, y = chain.roles["SUM"]
        unit_ids = [o + 1] + [-(node + 1) for node in (*x, *y)]
    clauses = 3 + 4 * count["xor"] + 3 * count["and"] + len(unit_ids)
    size = len(f"p cnf {variables} {clauses}\n") + len("-1 0\n") + len("2 0\n")
    for node_id, (op, a, b) in enumerate(d.nodes[2:], 2):
        if op == "var":
            continue
        da, db, dz = len(str(a + 1)), len(str(b + 1)), len(str(node_id + 1))
        if op == "xor":
            size += 4 * (da + db + dz) + 26
        elif op == "and":
            size += 2 * da + 2 * db + 3 * dz + 17
        else:
            raise AssertionError("unknown node")
    size += len(f"{chain.output + 1} 0\n")
    size += sum(len(str(unit)) + 3 for unit in unit_ids)
    return {"variables": variables, "clauses": clauses, "bytes": size,
            "unit_count": len(unit_ids), "target": "O" if target_O else "generic",
            "dag": count, "dag_prefix_sha256": d.prefix_sha256()}


def read_progress(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def check_manifest(out: Path, receipt: dict):
    manifest_path = out / "MANIFEST.json"
    assert sha(manifest_path) == receipt["manifest_sha256"]
    rows = json.loads(manifest_path.read_text())
    expected_paths = {row["path"] for row in rows}
    actual_paths = {p.relative_to(out).as_posix()
                    for p in out.rglob("*") if p.is_file()
                    and p.name not in ("MANIFEST.json", "receipt.json")}
    assert expected_paths == actual_paths
    for row in rows:
        path = out / row["path"]
        assert path.is_file() and path.stat().st_size == row["bytes"]
        assert sha(path) == row["sha256"]


def timeout(_signum, _frame):
    raise TimeoutError("independent replay wall cap reached")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite independent replay receipt")
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    spec = json.loads((HERE / "INPUT.json").read_text())
    assert frozen["release_main_head"] is not None
    resource.setrlimit(resource.RLIMIT_AS,
                       (spec["caps"]["per_arm_address_space_and_rss_bytes"],
                        spec["caps"]["per_arm_address_space_and_rss_bytes"]))
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(2 * spec["caps"]["per_arm_external_wall_seconds"])
    out = args.evidence.resolve()
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["status"] == "PRODUCER_ARCHIVED"
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["release"]["main_head"] == frozen["release_main_head"]
    assert set(receipt["release"]["merged_parent_commits"]) == {"802", "804", "784"}
    for merge_oid in receipt["release"]["merged_parent_commits"].values():
        subprocess.run(["git", "merge-base", "--is-ancestor", merge_oid,
                        frozen["release_main_head"]], cwd=ROOT, check=True)
    check_manifest(out, receipt)
    for rel, expected in frozen["input_sha256"].items():
        assert sha(ROOT / rel) == expected, rel
    for name, expected in frozen["source_sha256"].items():
        assert sha(HERE / name) == expected, name
    basis = replay_basis()
    assert basis["status"] == "PASS"
    arms = []
    for attempt in receipt["attempts"]:
        name = attempt["arm"]
        assert name in ("balanced", "unequal")
        assert attempt["exit_code"] == 0 and attempt["stop_reason"] is None
        source = json.loads((out / name / "result.json").read_text())
        assert sha(out / name / "result.json") == attempt["result_sha256"]
        assert source["arm"] == name and source["freeze_sha256"] == receipt["freeze_sha256"]
        assert source["basis_replay"] == basis["arms"][name]
        assert source["dimensions"] == next(row["dimensions"] for row in spec["arms"]
                                             if row["name"] == name)
        checkpoints = []
        expected_bases = slot_bases_n131(name)
        try:
            chain = build_chain(spec["field_degree"], int(spec["field_modulus_hex"], 16),
                                expected_bases, spec["caps"]["dag_nodes"],
                                checkpoint=checkpoints.append)
        except NodeCapExceeded as exc:
            assert source["status"] == "CENSORED_DAG_NODE_CAP"
            assert source["complete_chain"] is False
            assert source["cap_stage"] == exc.stage
            assert source["partial_dag_counts"] == exc.counts
            assert source["partial_prefix_sha256"] == exc.prefix_sha256
            replay_status = "REPRODUCED_DAG_CAP"
        else:
            assert source["complete_chain"] is True
            assert source["dag_counts"] == chain.counts()
            assert source["dag_prefix_sha256"] == chain.dag.prefix_sha256()
            assert source["local_relation_counts"] == chain.local_relation_counts
            assert source["edge_count"] == 9
            assert source["edge_labels"] == [list(edge) for edge in chain.edges]
            assert source["primary_inputs_expected"] == chain.counts()["variables"]
            generic = independent_dimacs_size(chain, False)
            target_o = independent_dimacs_size(chain, True)
            assert generic == source["cnf_generic"]
            assert target_o == source["cnf_target_O"]
            assert target_o["clauses"] == generic["clauses"] + 263
            expected_status = ("CENSORED_CNF_BYTE_CAP" if target_o["bytes"] >
                               spec["caps"]["dimacs_bytes"]
                               else "WITHIN_FROZEN_REPRESENTATION_CAPS")
            assert source["status"] == expected_status
            replay_status = expected_status
        assert read_progress(out / name / "progress.jsonl") == checkpoints
        arms.append({"arm": name, "replay_status": replay_status,
                     "source_result_sha256": attempt["result_sha256"]})
    assert [row["arm"] for row in arms] == ["balanced", "unequal"]
    result = {"schema": "ecc2k130-m10-capacity-independent-replay-v1",
              "status": "PASS", "scope": "representation size only; no solver or PDP",
              "producer_receipt_sha256": sha(out / "receipt.json"),
              "freeze_sha256": receipt["freeze_sha256"],
              "arms": arms,
              "decision": ("BOTH_WITHIN_FROZEN_CAPS" if all(
                  row["replay_status"] == "WITHIN_FROZEN_REPRESENTATION_CAPS"
                  for row in arms) else "CAPACITY_CENSORED")}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    signal.alarm(0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
