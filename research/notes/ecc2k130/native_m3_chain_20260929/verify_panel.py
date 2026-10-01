#!/usr/bin/env python3
"""Replay raw SAT models/CNF and leaf capacity with independent arithmetic."""
from __future__ import annotations

import argparse
from collections import Counter
import importlib.util
import itertools
import json
import resource
import sys
import tempfile
import time
from pathlib import Path

from chain import (EXACT_MODULUS, NodeCapExceeded, build_native_chain,
                   leaf_slot_bases)
from reference import (TOY_A, TOY_B, TOY_BASES, TOY_MODULUS, TOY_N,
                       PolynomialCurve, decode_model, independent_leaf_bases,
                       parse_cnf_and_model, primary_model, sha)
from panel import DOMAIN, EXACT, EXACT_REPLAY, target_units

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
sys.path.insert(0, str(HERE.parent / "native_fullpoint_edge_20260929"))
sys.path.insert(0, str(HERE.parent / "symbolic_dag_dimacs_gate_20260925"))
from relations import Koblitz  # noqa: E402
from export import write_cnf  # noqa: E402

_edge_producer_path = HERE.parent / "native_fullpoint_edge_20260929/produce.py"
_edge_spec = importlib.util.spec_from_file_location("frozen_native_edge_producer", _edge_producer_path)
_edge_module = importlib.util.module_from_spec(_edge_spec)
assert _edge_spec.loader is not None
_edge_spec.loader.exec_module(_edge_module)
BitField = _edge_module.BitField

O = (1, 0, 0)


def alternate_branch(p, q):
    if p == O:
        return "copy_q"
    if q == O:
        return "copy_p"
    if p[1] == q[1] and p[2] ^ q[2] == p[1]:
        return "inverse"
    if p[1] == q[1]:
        return "double"
    return "generic"


def alternate_slope(curve, p, q):
    kind = alternate_branch(p, q)
    field = curve.F
    if kind in ("copy_q", "copy_p", "inverse"):
        return 0
    if kind == "double":
        return p[1] ^ field.mul(p[2], field.inv(p[1]))
    return field.mul(p[2] ^ q[2], field.inv(p[1] ^ q[1]))


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def alternate_toy_support():
    field = BitField(TOY_N, TOY_MODULUS)
    curve = Koblitz(field, a=TOY_A, b=TOY_B)
    points = [O] + [(0, x, y) for x in range(1 << TOY_N)
                    for y in range(1 << TOY_N) if curve.on_curve((x, y))]
    factors = [tuple(p for p in points if p != O and p[1] in (0, basis[0]))
               for basis in TOY_BASES]
    support = set()
    witnesses = {}
    for p, q, r in itertools.product(*factors):
        intermediate = curve.add(p[1:], q[1:])
        result = curve.add(intermediate, r[1:])
        total = O if result is None else (0, *result)
        s2 = O if intermediate is None else (0, *intermediate)
        support.add(total)
        witnesses.setdefault(total, []).append(((p, q, r), s2))
    positive = min(p for p in support if p != O)
    generic = min(point for point, paths in witnesses.items()
                  if point != O and all(alternate_branch(path[0][0], path[0][1]) == "generic"
                                        and alternate_branch(path[1], path[0][2]) == "generic"
                                        for path in paths))
    negative = min(p for p in points if p not in support)
    return {"curve": curve, "points": points, "factors": factors,
            "support": support, "positive": positive,
            "generic_positive": generic, "negative": negative}


def verify_toy(producer: Path, summary: dict):
    alt = alternate_toy_support()
    assert summary["factor_sizes"] == list(map(len, alt["factors"]))
    assert summary["factor_triples"] == __import__("math").prod(summary["factor_sizes"])
    assert summary["distinct_supported_targets"] == len(alt["support"])
    assert summary["rational_curve_points"] == len(alt["points"])
    assert summary["positive"] == list(alt["positive"])
    assert summary["generic_positive"] == list(alt["generic_positive"])
    assert summary["negative"] == list(alt["negative"])
    chain = build_native_chain(TOY_N, TOY_MODULUS, TOY_A, TOY_B,
                               TOY_BASES, 20000)
    assert chain.counts() == summary["dag"]
    assert chain.dag.prefix_sha256() == summary["prefix_sha256"]
    assert sha(producer / "solver.version.txt") == summary["solver_version_sha256"]
    assert "CryptoMiniSat version 5." in (producer / "solver.version.txt").read_text()
    rows_path = producer / "triples.jsonl"
    assert sha(rows_path) == summary["triple_rows_sha256"]
    branch_counts = Counter()
    with rows_path.open() as stream:
        row_count = 0
        for triple in itertools.product(*alt["factors"]):
            p, q, r = triple
            intermediate = alt["curve"].add(p[1:], q[1:])
            s2 = O if intermediate is None else (0, *intermediate)
            got = alt["curve"].add(intermediate, r[1:])
            total = O if got is None else (0, *got)
            wrong = next(point for point in sorted(alt["points"]) if point != total)
            slopes = (alternate_slope(alt["curve"], p, q),
                      alternate_slope(alt["curve"], s2, r))
            masks = [int(point[1] == basis[0])
                     for point, basis in zip(triple, TOY_BASES, strict=True)]
            kinds = (alternate_branch(p, q), alternate_branch(s2, r))
            branch_counts[f"edge0_{kinds[0]}"] += 1
            branch_counts[f"edge1_{kinds[1]}"] += 1
            correct_model = primary_model(chain, triple, s2, total, slopes, masks)
            wrong_model = primary_model(chain, triple, s2, wrong, slopes, masks)
            accepted = chain.dag.evaluate(correct_model, chain.output)
            rejected = not chain.dag.evaluate(wrong_model, chain.output)
            expected = {"factors": [list(point) for point in triple],
                        "masks": masks, "s2": list(s2), "target": list(total),
                        "wrong_target": list(wrong), "slopes": list(slopes),
                        "branches": list(kinds), "accepted": accepted,
                        "wrong_rejected": rejected}
            line = stream.readline()
            assert line and json.loads(line) == expected
            assert accepted and rejected
            row_count += 1
        assert stream.readline() == ""
    assert row_count == summary["triple_rows"] == 27
    assert dict(sorted(branch_counts.items())) == summary["branch_counts"]
    assert branch_counts["edge0_generic"] == branch_counts["edge1_generic"] == 24
    checks = []
    for label, target in (("positive", alt["positive"]),
                          ("generic_positive", alt["generic_positive"]),
                          ("negative", alt["negative"])):
        run = next(row for row in summary["runs"] if row["label"] == label)
        assert run["target"] == list(target)
        cnf = producer / f"{label}.cnf"
        stdout = producer / f"{label}.stdout.txt"
        stderr = producer / f"{label}.stderr.txt"
        assert sha(cnf) == run["cnf"]["sha256"]
        assert sha(stdout) == run["stdout_sha256"]
        assert sha(stderr) == run["stderr_sha256"]
        with tempfile.TemporaryDirectory(prefix="native-m3-cnf-") as temp:
            duplicate = Path(temp) / "duplicate.cnf"
            fresh = write_cnf(chain, duplicate, units=target_units(chain, target),
                              byte_cap=16 * 1024 * 1024)
            assert fresh == run["cnf"] and sha(duplicate) == sha(cnf)
        # This parser truth-tables every gate's local CNF template.
        verifier_path = HERE.parent / "symbolic_dag_dimacs_gate_20260925/verify.py"
        spec = importlib.util.spec_from_file_location("frozen_dimacs_verifier", verifier_path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        variables, clauses, rows = module.parse_cnf(cnf)
        assert (variables, clauses) == (run["cnf"]["variables"], run["cnf"]["clauses"])
        module.check_relation_cnf(chain, rows, target_units(chain, target))
        status, bits, _, _ = parse_cnf_and_model(cnf, stdout.read_text(), run["exit_code"])
        assert status == run["status"] == ("UNSAT" if label == "negative" else "SAT")
        assert (target in alt["support"]) == (status == "SAT")
        if status == "SAT":
            lifted = decode_model(chain, bits, PolynomialCurve(TOY_N, TOY_MODULUS,
                                                                TOY_A, TOY_B), target)
            assert lifted == run["lift"]
            factors = [tuple(point) for point in lifted["factors"]]
            assert all(point in alt["factors"][i] for i, point in enumerate(factors))
            got = alt["curve"].add(alt["curve"].add(factors[0][1:], factors[1][1:]),
                                   factors[2][1:])
            assert (O if got is None else (0, *got)) == target
            assert bits[chain.output] == 1
            assert all(bits[node] == 1 for node in chain.edge_outputs)
            if label == "generic_positive":
                assert alternate_branch(factors[0], factors[1]) == "generic"
                assert alternate_branch(tuple(lifted["s2"]), factors[2]) == "generic"
        else:
            assert run["classification"] == "SOLVER_UNSAT_ORACLE_CONFIRMED_TOY"
            assert bits is None and run["lift"] is None
        checks.append({"label": label, "status": status,
                       "cnf_sha256": sha(cnf), "variables": variables,
                       "clauses": clauses})
    return {"arm": "toy", "decision": "PASS", "checks": checks,
            "triple_rows": row_count, "triple_rows_sha256": sha(rows_path),
            "branch_counts": dict(sorted(branch_counts.items())),
            "factor_triples": summary["factor_triples"],
            "supported_targets": len(alt["support"])}


def independent_dimacs_size(chain, target_o: bool) -> dict:
    d = chain.dag
    count = d.counts()
    variables = count["total_nodes"]
    units = []
    if target_o:
        o, x, y = chain.roles["SUM"]
        units = [o + 1] + [-(node + 1) for node in (*x, *y)]
    clauses = 3 + 4 * count["xor"] + 3 * count["and"] + len(units)
    size = len(f"p cnf {variables} {clauses}\n") + len("-1 0\n") + len("2 0\n")
    for i, (op, a, b) in enumerate(d.nodes[2:], 2):
        if op == "var":
            continue
        da, db, dz = len(str(a + 1)), len(str(b + 1)), len(str(i + 1))
        if op == "xor":
            size += 4 * (da + db + dz) + 26
        elif op == "and":
            size += 2 * da + 2 * db + 3 * dz + 17
        else:
            raise AssertionError("invalid DAG node")
    size += len(f"{chain.output + 1} 0\n")
    size += sum(len(str(unit)) + 3 for unit in units)
    return {"variables": variables, "clauses": clauses, "bytes": size,
            "unit_count": len(units), "target": "O" if target_o else "generic",
            "dag": count, "dag_prefix_sha256": d.prefix_sha256()}


def verify_leaf(producer: Path, summary: dict, arm: str):
    archived = json.loads(EXACT.read_text())
    replay = json.loads(EXACT_REPLAY.read_text())
    assert archived["status"] == replay["status"] == "PASS"
    assert replay["input_sha256"] == sha(EXACT)
    line = [1, 0] if arm == "leaf10" else [1, 4]
    record = next(item for item in archived["lines"] if item["line"] == line)
    b = int(record["codomain_b"], 16)
    assert summary["line"] == line and summary["curve_b"] == hex(b)
    assert summary["exact_smoke_sha256"] == sha(EXACT)
    assert summary["exact_replay_sha256"] == sha(EXACT_REPLAY)
    bases = independent_leaf_bases()
    assert bases == leaf_slot_bases()
    try:
        chain = build_native_chain(131, EXACT_MODULUS, 0, b, bases, 2_000_000)
    except NodeCapExceeded as exc:
        assert summary["decision"] == "CAPACITY_CENSORED"
        assert summary["censor"] == "DAG_NODE_CAP"
        assert summary["stage"] == exc.stage
        assert summary["partial_counts"] == exc.counts
        assert summary["partial_prefix_sha256"] == exc.prefix_sha256
        return {"arm": arm, "decision": "CAPACITY_CENSORED",
                "censor": "DAG_NODE_CAP", "partial_counts": exc.counts}
    assert chain.counts() == summary["dag"]
    assert summary["prefix_sha256"] == chain.dag.prefix_sha256()
    assert summary["edge_count"] == 2 and summary["dimensions"] == [44, 44, 43]
    assert summary["cnf_generic"] == independent_dimacs_size(chain, False)
    assert summary["cnf_target_O"] == independent_dimacs_size(chain, True)
    fits = summary["cnf_target_O"]["bytes"] <= 256 * 1024 * 1024
    assert summary["decision"] == ("PASS" if fits else "CAPACITY_CENSORED")
    assert summary["censor"] == (None if fits else "CNF_BYTE_CAP")
    progress = [json.loads(line) for line in (producer / "progress.jsonl").read_text().splitlines()]
    assert [row["stage"] for row in progress] == ["primary_wires", "edge_0", "edge_1", "chain_output"]
    assert progress[-1]["prefix_sha256"] == summary["prefix_sha256"]
    assert progress[-1]["counts"] == summary["dag"]
    return {"arm": arm, "decision": summary["decision"], "line": line,
            "dag": summary["dag"], "cnf_target_O": summary["cnf_target_O"]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=("toy", "leaf10", "leaf14"), required=True)
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started, cpu_started = time.monotonic(), time.process_time()
    summary = json.loads((args.producer / "result.json").read_text())
    assert summary["domain"] == DOMAIN and summary["arm"] == args.arm
    assert summary["decision"] in ("PASS", "CAPACITY_CENSORED")
    result = (verify_toy(args.producer, summary) if args.arm == "toy" else
              verify_leaf(args.producer, summary, args.arm))
    result.update({"domain": DOMAIN, "producer_sha256": sha(args.producer / "result.json"),
                   "wall_seconds": time.monotonic() - started,
                   "cpu_seconds": time.process_time() - cpu_started,
                   "peak_rss_bytes": rss_bytes()})
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"arm": args.arm, "decision": result["decision"],
                      "wall_seconds": result["wall_seconds"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
