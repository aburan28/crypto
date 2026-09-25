#!/usr/bin/env python3
"""Capped n19 explicit-state sparse S3 CNF producer. No run before parent merge."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import resource
import signal
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
DENSE = NOTES / "rotated_s3_o_branch_20260925/export.py"
SPARSE = NOTES / "rotated_s3_sparse_cnf_20260925/export.py"
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
TARGETS = NOTES / "rotated_s3_candidate_20260925/evidence/n19-m6/producer/result.json"
N, M, POLY = 19, 6, 0x80027
STATE_CAP = 100_000
PAIR_CAP = 500_000
PATH_CAP = 250_000
BYTE_CAP = 20_000_000
WALL_CAP = 180
RSS_CAP = 512 * 1024 * 1024


class Censored(RuntimeError):
    pass


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def check_freeze(frozen: dict) -> None:
    for key, path in (("protocol_sha256", HERE / "PROTOCOL.md"),
                      ("export_sha256", Path(__file__)),
                      ("dense_export_sha256", DENSE),
                      ("sparse_export_sha256", SPARSE),
                      ("corpus_sha256", CORPUS),
                      ("target_panel_sha256", TARGETS)):
        assert sha(path) == frozen[key], key


def input_rows():
    with tarfile.open(CORPUS, "r:gz") as tar:
        factors = json.load(tar.extractfile("raw/n19-m6/factors.json"))
    archive = json.loads(TARGETS.read_text())
    assert archive["arm"] == "n19-m6" and len(archive["targets"]) == 32
    targets = [{"id": f"Q{row['Q_index']}T{row['T_index']}",
                "point": row["target"], "class": row["target_class"],
                "archived_exact_tuple_count": row["true_point_tuple_count"]}
               for row in archive["targets"]]
    assert [(row["Q_index"], row["T_index"]) for row in archive["targets"]] == [
        (i, j) for i in range(8) for j in range(4)]
    targets.append({"id": "O", "point": None, "class": "identity",
                    "archived_exact_tuple_count": None})
    return factors, targets


def capped_growth(dense, field, domains, targets, out: Path):
    reachable = None
    states = []
    growth = []
    total_pairs = 0
    for length in range(2, M + 1):
        left = domains[0] if length == 2 else states[-1]
        right = domains[1] if length == 2 else domains[length - 1]
        pairs = len(left) * len(right)
        total_pairs += pairs
        if total_pairs > PAIR_CAP:
            save(out / "growth.json", {"completed": growth, "censored_stage": length,
                                       "prospective_stage_pairs": pairs,
                                       "prospective_total_pairs": total_pairs})
            raise Censored(f"transition input-pair cap at S{length}")
        if length == 2:
            next_reachable = {v for a in domains[0] for b in domains[1]
                              for v in dense.local(field, a, b)}
        else:
            next_reachable = {v for u in reachable for a in domains[length - 1]
                              for v in dense.local(field, u, a)}
        admitted = set(next_reachable) | {None}
        if length == M:
            admitted.update(None if row["point"] is None else row["point"][0]
                            for row in targets)
        values = sorted(admitted, key=dense.state_sort)
        row = {"state_group": f"S{length}", "reachable_states": len(next_reachable),
               "admitted_states": len(values), "reachable_O": None in next_reachable,
               "admitted_O": None in admitted, "stage_input_pairs": pairs,
               "cumulative_input_pairs": total_pairs, "peak_rss_bytes": rss()}
        growth.append(row)
        save(out / "growth.json", {"completed": growth})
        if len(values) > STATE_CAP:
            raise Censored(f"admitted state cap at S{length}")
        if rss() > RSS_CAP:
            raise Censored(f"producer RSS cap at S{length}")
        states.append(values)
        reachable = next_reachable
    return states, growth, total_pairs


def groups_for(domains, states):
    groups = []
    next_var = 1
    for name, values in ([(f"F{i}", values) for i, values in enumerate(domains)] +
                         [(f"S{i}", values) for i, values in enumerate(states, start=2)]):
        variables = list(range(next_var, next_var + len(values)))
        groups.append({"name": name, "values": values, "vars": variables})
        next_var += len(values)
    return groups, next_var - 1


def enumerate_paths(dense, field, domains, out: Path):
    rows = []
    for tuple_index, xs in enumerate(itertools.product(*domains), start=1):
        paths = [(state,) for state in dense.local(field, xs[0], xs[1])]
        for factor in xs[2:]:
            paths = [path + (state,) for path in paths
                     for state in dense.local(field, path[-1], factor)]
        for path in paths:
            rows.append({"factor_x": list(xs), "states": list(path)})
        if len(rows) > PATH_CAP:
            save(out / "path_growth.json", {"completed_factor_x_tuples":
                 tuple_index, "primary_paths_so_far": len(rows), "limit": PATH_CAP})
            raise Censored("primary path cap")
    rows.sort(key=lambda row: (row["factor_x"],
                               [dense.state_sort(v) for v in row["states"]]))
    with (out / "paths.jsonl.gz").open("wb") as stream:
        with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as zipped:
            for row in rows:
                zipped.write((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode())
    return rows


def write_cnf(dense, sparse, field, groups, stage_rows, variables, out: Path):
    by_name = {row["name"]: row for row in groups}
    next_var = variables + 1
    aux = []
    for group in groups:
        variables_row = list(range(next_var, next_var + len(group["vars"]) - 1))
        aux.append({"name": group["name"], "vars": variables_row})
        next_var += len(variables_row)
    onehot = sum(1 if len(group["vars"]) == 1 else 3 * len(group["vars"]) - 3
                 for group in groups)
    pairs = sum(len(by_name[row["left"]]["vars"]) * len(by_name[row["right"]]["vars"])
                for row in stage_rows)
    count = onehot + pairs
    total_variables = next_var - 1
    case_counts = []
    written = 0
    path = out / "base.cnf"
    with path.open("w") as stream:
        def line(value: str):
            nonlocal written
            written += stream.write(value)
            if written > BYTE_CAP:
                raise Censored("DIMACS byte cap")
        line("c frozen n19 sparse O-aware recursive S3 base\n")
        line(f"p cnf {total_variables} {count}\n")
        for group, auxiliaries in zip(groups, aux):
            for clause in sparse.sequential_group(group["vars"], auxiliaries["vars"]):
                line(" ".join(map(str, clause)) + " 0\n")
        for stage in stage_rows:
            left, right, end = (by_name[stage[key]] for key in ("left", "right", "out"))
            cases = Counter()
            end_values = set(end["values"])
            for u, vu in zip(left["values"], left["vars"]):
                for a, va in zip(right["values"], right["vars"]):
                    permitted = set(dense.local(field, u, a))
                    outputs = [vv for v, vv in zip(end["values"], end["vars"])
                               if v in permitted]
                    assert len(outputs) == len(permitted & end_values)
                    cases[len(outputs)] += 1
                    line(" ".join(map(str, [-vu, -va, *outputs])) + " 0\n")
            case_counts.append({"stage": stage["out"],
                                "admitted_output_size_counts": dict(sorted(cases.items()))})
    assert path.stat().st_size == written and sum(sum(x["admitted_output_size_counts"].values())
                                                     for x in case_counts) == pairs
    return aux, total_variables, onehot, pairs, count, case_counts


def run(out: Path):
    started = time.perf_counter()
    cpu = time.process_time()
    out.mkdir(parents=True, exist_ok=False)
    def expired(_signum, _frame):
        raise Censored("exporter wall cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, WALL_CAP)
    try:
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        check_freeze(frozen)
        dense = module(DENSE, "n19_dense_root_source")
        sparse = module(SPARSE, "n19_sparse_block_source")
        field = dense.GF(N, POLY)
        factors, targets = input_rows()
        domains = dense.validate_factors(field, [[tuple(p) for p in slot] for slot in factors])
        assert len(domains) == M and [len(slot) for slot in factors] == [7] * M
        assert [len(slot) for slot in domains] == [4] * M
        assert all(row["point"] is None or tuple(row["point"]) in field.fibre(row["point"][0])
                   for row in targets)
        states, growth, total_pairs = capped_growth(dense, field, domains, targets, out)
        groups, primary = groups_for(domains, states)
        stage_rows = dense.stages(M)
        aux, variables, onehot, pairs, clauses, case_counts = write_cnf(
            dense, sparse, field, groups, stage_rows, primary, out)
        assert pairs == total_pairs
        target_literals = {v: literal for v, literal in zip(groups[-1]["values"],
                                                             groups[-1]["vars"])}
        targets = [{**row, "assumption_literal": target_literals[
            None if row["point"] is None else row["point"][0]]} for row in targets]
        schema = {"panel": "n19-m6", "field_degree": N, "field_poly": POLY, "m": M,
                  "factor_points": factors, "targets": targets, "groups": groups,
                  "stages": stage_rows, "auxiliary_groups": aux,
                  "encoding": "sinz-sequential-amo-and-allowed-output-v1",
                  "primary_variables": primary, "variables": variables,
                  "clauses": clauses, "onehot_clauses": onehot,
                  "implication_clauses": pairs}
        save(out / "schema.json", schema)
        paths = enumerate_paths(dense, field, domains, out)
        result = {"decision": "PRODUCED", "panel": "n19-m6", "growth": growth,
                  "transition_cases": case_counts, "factor_x_tuples": math.prod(map(len, domains)),
                  "signed_point_tuples": math.prod(map(len, factors)),
                  "primary_paths": len(paths), "target_labels": len(targets),
                  "primary_variables": primary, "auxiliary_variables": variables - primary,
                  "variables": variables, "clauses": clauses, "onehot_clauses": onehot,
                  "implication_clauses": pairs, "bytes": (out / "base.cnf").stat().st_size,
                  "schema_sha256": sha(out / "schema.json"),
                  "base_sha256": sha(out / "base.cnf"),
                  "paths_sha256": sha(out / "paths.jsonl.gz"),
                  "field_operations": dict(field.ops),
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu,
                  "peak_rss_bytes": rss()}
        if result["peak_rss_bytes"] > RSS_CAP:
            raise Censored("producer RSS cap")
        save(out / "result.json", result)
    except Exception as error:
        save(out / "failure.json", {"decision": "CENSORED" if isinstance(error, Censored)
                                     else "FAILED", "error": repr(error),
                                    "wall_seconds": time.perf_counter() - started,
                                    "cpu_seconds": time.process_time() - cpu,
                                    "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)


if __name__ == "__main__":
    main()
