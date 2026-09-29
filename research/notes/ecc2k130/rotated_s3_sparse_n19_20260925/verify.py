#!/usr/bin/env python3
"""Independent parsed-CNF and full signed-point replay of n19 sparse growth."""
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
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
TARGETS = NOTES / "rotated_s3_candidate_20260925/evidence/n19-m6/producer/result.json"
INDEPENDENT = NOTES / "rotated_pdp_corpus_20260925/verify.py"
SPARSE_VERIFY = NOTES / "rotated_s3_sparse_cnf_20260925/verify.py"
N, M, POLY = 19, 6, 0x80027
WALL_CAP = 600
RSS_CAP = 512 * 1024 * 1024
PATH_CAP = 250_000


class SemanticError(AssertionError):
    pass


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


def state_sort(value):
    return -1 if value is None else value


def check_freeze(frozen):
    for key, path in (("verify_sha256", Path(__file__)),
                      ("independent_sha256", INDEPENDENT),
                      ("sparse_verify_sha256", SPARSE_VERIFY),
                      ("corpus_sha256", CORPUS),
                      ("target_panel_sha256", TARGETS)):
        assert sha(path) == frozen[key], key


def input_rows(independent):
    with tarfile.open(CORPUS, "r:gz") as tar:
        factors = [[tuple(p) for p in slot]
                   for slot in json.load(tar.extractfile("raw/n19-m6/factors.json"))]
        corpus_targets = json.load(tar.extractfile("raw/n19-m6/targets.json"))
    archive = json.loads(TARGETS.read_text())["targets"]
    assert len(factors) == 6 and all(len(slot) == 7 for slot in factors)
    assert len(corpus_targets) == 8 and len(archive) == 32
    f = independent.parent_verify.GF(N, POLY)
    curve = independent.parent_verify.E(f)
    torsion = [None, (0, 1), (1, 0), (1, 1)]
    assert all(curve.on(t) for t in torsion)
    targets = []
    for i, row in enumerate(corpus_targets):
        Q = tuple(row["Q"])
        assert curve.on(Q)
        for j, T in enumerate(torsion):
            point = curve.add(Q, T)
            expected = archive[4 * i + j]
            if (expected["Q_index"], expected["T_index"], expected["target"],
                expected["target_class"]) != (i, j, list(point) if point else None,
                                                row["class"]):
                raise SemanticError("frozen Q+T target mismatch")
            targets.append({"id": f"Q{i}T{j}", "point": point,
                            "class": row["class"],
                            "archived_exact_tuple_count": expected["true_point_tuple_count"]})
    targets.append({"id": "O", "point": None, "class": "identity",
                    "archived_exact_tuple_count": None})
    return f, curve, factors, targets


def complete_fibres(independent, curve, factors):
    domains = []
    for slot in factors:
        actual = set(slot)
        xs = sorted({p[0] for p in slot})
        for x in xs:
            expected = set(independent.lifts(curve, x))
            if not expected or {p for p in actual if p[0] == x} != expected:
                raise SemanticError(f"sign-incomplete or nonlift factor x={x}")
        domains.append(xs)
    return domains


def make_geometry(independent, curve):
    @lru_cache(maxsize=1 << N)
    def fibre(x):
        return (None,) if x is None else tuple(independent.lifts(curve, x))

    @lru_cache(maxsize=500_000)
    def local(u, a):
        return frozenset(None if point is None else point[0]
                         for left in fibre(u) for right in fibre(a)
                         for point in (curve.add(left, right),))
    return fibre, local


def check_target_row(actual, expected, final_group):
    state = None if expected["point"] is None else expected["point"][0]
    literal = final_group["vars"][final_group["values"].index(state)]
    reference = {**expected, "point": None if expected["point"] is None
                 else list(expected["point"]), "assumption_literal": literal}
    if actual != reference:
        raise SemanticError(f"target point/assumption {expected['id']}")


def validate_schema(schema, factors, targets, domains, fibre):
    if (schema["panel"], schema["field_degree"], schema["field_poly"], schema["m"]) != (
            "n19-m6", N, POLY, M):
        raise SemanticError("panel/field metadata")
    if schema["factor_points"] != [[list(p) for p in slot] for slot in factors]:
        raise SemanticError("factor point mismatch")
    if [g["name"] for g in schema["groups"]] != ([f"F{i}" for i in range(M)] +
                                                  [f"S{i}" for i in range(2, M + 1)]):
        raise SemanticError("group order")
    if schema["stages"] != ([{"left": "F0", "right": "F1", "out": "S2"}] +
                            [{"left": f"S{i}", "right": f"F{i}", "out": f"S{i+1}"}
                             for i in range(2, M)]):
        raise SemanticError("stage order")
    next_var = 1
    for i, group in enumerate(schema["groups"]):
        values = group["values"]
        if (not values or len(values) != len(set(values)) or
            group["vars"] != list(range(next_var, next_var + len(values)))):
            raise SemanticError("primary group values or variables")
        if i < M:
            if values != domains[i]:
                raise SemanticError("factor x domain mismatch")
        else:
            if values != sorted(values, key=state_sort) or values[0] is not None:
                raise SemanticError("prefix domain order/O")
            if any(not fibre(v) for v in values):
                raise SemanticError("nonrational prefix x")
        next_var += len(values)
    if schema["primary_variables"] != next_var - 1:
        raise SemanticError("primary variable count")
    if len(schema["auxiliary_groups"]) != len(schema["groups"]):
        raise SemanticError("auxiliary group count")
    aux_next = next_var
    for group, aux in zip(schema["groups"], schema["auxiliary_groups"]):
        if aux["name"] != group["name"] or aux["vars"] != list(
                range(aux_next, aux_next + len(group["vars"]) - 1)):
            raise SemanticError("auxiliary group numbering")
        aux_next += len(aux["vars"])
    if schema["variables"] != aux_next - 1:
        raise SemanticError("variable count")
    if len(schema["targets"]) != 33:
        raise SemanticError("target count")
    terminal = schema["groups"][-1]
    for actual, expected in zip(schema["targets"], targets):
        check_target_row(actual, expected, terminal)
    return {row["name"]: row for row in schema["groups"]}


def check_pair(row, vu, va, end, u, a, local, label):
    expected_values = local(u, a)
    expected = tuple([-vu, -va, *[vv for v, vv in zip(end["values"], end["vars"])
                                  if v in expected_values]])
    if row != expected:
        raise SemanticError(f"point-law implication mismatch at {label}: {(u, a)}")
    return len(expected) - 2


def parse_cnf(sv, path, schema, by_name, local):
    reader = sv.CNFReader(path)
    if (reader.variables, reader.declared) != (schema["variables"], schema["clauses"]):
        raise SemanticError("DIMACS header")
    onehot = 0
    for group, aux in zip(schema["groups"], schema["auxiliary_groups"]):
        expected = sv.expected_block(group["vars"], aux["vars"])
        actual = [reader.take() for _ in expected]
        if actual != expected:
            raise SemanticError(f"one-hot clause mismatch {group['name']}")
        sv.check_onehot_block(group["vars"], aux["vars"], actual)
        onehot += len(actual)
    if onehot != schema["onehot_clauses"]:
        raise SemanticError("one-hot clause count")
    tables = []
    cases = []
    located = {}
    geometry_cases = Counter()
    for stage in schema["stages"]:
        left, factor, end = (by_name[stage[key]] for key in ("left", "right", "out"))
        table = {}
        counts = Counter()
        for u, vu in zip(left["values"], left["vars"]):
            for a, va in zip(factor["values"], factor["vars"]):
                row = reader.take()
                size = check_pair(row, vu, va, end, u, a, local, stage["out"])
                geometry_cases["O_plus_factor" if u is None else
                               "equal_x_zero" if u == a == 0 else
                               "equal_x_nonzero" if u == a else "distinct_x"] += 1
                table[u, a] = tuple(v for v, vv in zip(end["values"], end["vars"])
                                    if vv in row[2:])
                counts[size] += 1
                located[stage["out"], u, a] = (row, vu, va, end)
        tables.append(table)
        cases.append({"stage": stage["out"],
                      "admitted_output_size_counts": {str(k): v for k, v in sorted(counts.items())}})
    reader.finish()
    if onehot + sum(sum(c["admitted_output_size_counts"].values())
                    for c in cases) != schema["clauses"]:
        raise SemanticError("total clause count")
    return tables, cases, located, dict(geometry_cases)


def expected_prefix_domains(domains, targets, local):
    reachable = None
    states = []
    growth = []
    cumulative_pairs = 0
    for length in range(2, M + 1):
        left = domains[0] if length == 2 else states[-1]
        right = domains[1] if length == 2 else domains[length - 1]
        pairs = len(left) * len(right)
        cumulative_pairs += pairs
        if length == 2:
            next_reachable = {v for a in domains[0] for b in domains[1]
                              for v in local(a, b)}
        else:
            next_reachable = {v for u in reachable for a in domains[length - 1]
                              for v in local(u, a)}
        admitted = set(next_reachable) | {None}
        if length == M:
            admitted.update(None if row["point"] is None else row["point"][0]
                            for row in targets)
        states.append(sorted(admitted, key=state_sort))
        growth.append({"state_group": f"S{length}",
                       "reachable_states": len(next_reachable),
                       "admitted_states": len(admitted),
                       "reachable_O": None in next_reachable,
                       "admitted_O": None in admitted,
                       "stage_input_pairs": pairs,
                       "cumulative_input_pairs": cumulative_pairs})
        reachable = next_reachable
    return states, growth


def parsed_paths(domains, tables):
    rows = []
    for xs in itertools.product(*domains):
        paths = [(v,) for v in tables[0][xs[0], xs[1]]]
        for stage, x in enumerate(xs[2:], start=1):
            paths = [path + (v,) for path in paths
                     for v in tables[stage][path[-1], x]]
        rows.extend((tuple(xs), path) for path in paths)
        if len(rows) > PATH_CAP:
            raise Censored("independent primary path cap")
    rows.sort(key=lambda item: (item[0], tuple(map(state_sort, item[1]))))
    return rows


def signed_oracle(curve, factors):
    geometry = defaultdict(set)
    target_counts = Counter()
    prefixes = Counter()
    for points in itertools.product(*factors):
        summed = curve.add(points[0], points[1])
        states = [None if summed is None else summed[0]]
        for point in points[2:]:
            summed = curve.add(summed, point)
            states.append(None if summed is None else summed[0])
        key = (tuple(p[0] for p in points), tuple(states))
        geometry[key].add(summed)
        target_counts[summed] += 1
        prefixes["O_prefix" if None in states[:-1] else "all_affine"] += 1
        prefixes["O_terminal" if summed is None else "finite_terminal"] += 1
    return geometry, target_counts, prefixes


def controls(schema, factors, independent, curve, located, local):
    results = {}
    incomplete = [list(slot) for slot in factors]
    incomplete[0].remove(next(p for p in incomplete[0] if p[0] == 3))
    try:
        complete_fibres(independent, curve, incomplete)
    except SemanticError as error:
        results["sign_incomplete"] = str(error)
    else:
        raise SemanticError("incomplete factor fibre accepted")
    cases = (("zero_zero_to_O", "S2", 0, 0, False),
             ("O_plus_zero", "S3", None, 0, True),
             ("zero_double_to_O", "S4", 0, 0, False))
    for label, stage, u, a, wrong_O in cases:
        row, vu, va, end = located[stage, u, a]
        wanted = None if wrong_O else next(v for v in end["values"] if v is not None)
        wrong = (-vu, -va, end["vars"][end["values"].index(wanted)])
        if wrong == row:
            raise SemanticError(f"mutation did not change {label}")
        try:
            check_pair(wrong, vu, va, end, u, a, local, stage)
        except SemanticError as error:
            results[label] = {"original": row, "mutated": wrong, "rejection": str(error)}
        else:
            raise SemanticError(f"wrong {label} transition accepted")
    terminal = schema["groups"][-1]
    identity = schema["targets"][-1]
    wrong_literal = next(v for value, v in zip(terminal["values"], terminal["vars"])
                         if value is not None)
    expected_literal = terminal["vars"][terminal["values"].index(None)]
    if identity["assumption_literal"] != expected_literal or wrong_literal == expected_literal:
        raise SemanticError("O target setup")
    changed = dict(identity, assumption_literal=wrong_literal)
    expected = {"id": "O", "point": None, "class": "identity",
                "archived_exact_tuple_count": None}
    try:
        check_target_row(changed, expected, terminal)
    except SemanticError as error:
        results["wrong_O_target"] = {"original": expected_literal,
                                      "mutated": wrong_literal,
                                      "rejection": str(error)}
    else:
        raise SemanticError("wrong O target accepted")
    assert set(results) == {"sign_incomplete", "zero_zero_to_O", "O_plus_zero",
                            "zero_double_to_O", "wrong_O_target"}
    return results


def run(produced_dir: Path, out: Path):
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_signum, _frame):
        raise Censored("independent verifier wall cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, WALL_CAP)
    try:
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        check_freeze(frozen)
        produced = json.loads((produced_dir / "result.json").read_text())
        schema = json.loads((produced_dir / "schema.json").read_text())
        for file, field in (("base.cnf", "base_sha256"),
                            ("schema.json", "schema_sha256"),
                            ("paths.jsonl.gz", "paths_sha256")):
            if sha(produced_dir / file) != produced[field]:
                raise SemanticError(f"producer artifact hash {file}")
        independent = module(INDEPENDENT, "n19_independent_corpus_verify")
        sv = module(SPARSE_VERIFY, "n19_pinned_sparse_parser")
        f, curve, factors, targets = input_rows(independent)
        domains = complete_fibres(independent, curve, factors)
        fibre, local = make_geometry(independent, curve)
        by_name = validate_schema(schema, factors, targets, domains, fibre)
        tables, cases, located, geometry_cases = parse_cnf(
            sv, produced_dir / "base.cnf", schema, by_name, local)
        independent_states, independent_growth = expected_prefix_domains(
            domains, targets, local)
        for length, values in enumerate(independent_states, start=2):
            if values != by_name[f"S{length}"]["values"]:
                raise SemanticError(f"prefix state group S{length}")
        if len(produced["growth"]) != len(independent_growth) or any(
                {key: value for key, value in measured.items() if key != "peak_rss_bytes"}
                != expected for measured, expected in zip(produced["growth"], independent_growth)):
            raise SemanticError("prefix growth ledger")
        if (produced["variables"], produced["clauses"], produced["bytes"]) != (
                schema["variables"], schema["clauses"],
                (produced_dir / "base.cnf").stat().st_size):
            raise SemanticError("producer CNF size ledger")
        paths = parsed_paths(domains, tables)
        with gzip.open(produced_dir / "paths.jsonl.gz", "rt") as stream:
            archived = [(tuple(row["factor_x"]), tuple(row["states"]))
                        for row in (json.loads(line) for line in stream)]
        if paths != archived:
            raise SemanticError("producer/parsed primary paths")
        geometry, target_counts, point_cases = signed_oracle(curve, factors)
        if set(paths) != set(geometry) or len(paths) != len(geometry):
            raise SemanticError("parsed path/full signed-tuple oracle")
        for (_, states), point_set in geometry.items():
            if point_set != set(fibre(states[-1])):
                raise SemanticError("terminal sign completeness")
        target_rows = []
        for row in schema["targets"]:
            point = None if row["point"] is None else tuple(row["point"])
            state = None if point is None else point[0]
            selected = [key for key in paths if key[1][-1] == state]
            if any(point not in geometry[key] for key in selected):
                raise SemanticError(f"target terminal sign {row['id']}")
            count = target_counts[point]
            if (row["archived_exact_tuple_count"] is not None and
                    count != row["archived_exact_tuple_count"]):
                raise SemanticError(f"#774 exact target count {row['id']}")
            target_rows.append({"id": row["id"], "exact_signed_tuple_count": count,
                                "cnf_primary_path_count": len(selected),
                                "assumption_literal": row["assumption_literal"]})
        mutation_controls = controls(schema, factors, independent, curve, located, local)
        if cases != produced["transition_cases"]:
            raise SemanticError("transition case totals")
        if produced["primary_paths"] != len(paths) or produced["target_labels"] != 33:
            raise SemanticError("producer path/target totals")
        if produced["signed_point_tuples"] != math.prod(map(len, factors)):
            raise SemanticError("signed point tuple total")
        if produced["factor_x_tuples"] != math.prod(map(len, domains)):
            raise SemanticError("factor x tuple total")
        result = {"decision": "PASS", "primary_paths": len(paths),
                  "signed_point_tuples": sum(target_counts.values()),
                  "target_rows": target_rows, "transition_cases": cases,
                  "point_cases": dict(point_cases), "local_geometry_cases": geometry_cases,
                  "negative_controls": mutation_controls,
                  "producer_sha256": sha(produced_dir / "result.json"),
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu,
                  "peak_rss_bytes": rss()}
        if result["peak_rss_bytes"] > RSS_CAP:
            raise Censored("independent verifier RSS cap")
        save(out, result)
    except Exception as error:
        save(out, {"decision": "CENSORED" if isinstance(error, (Censored, MemoryError)) else "FAILED",
                   "error": repr(error), "wall_seconds": time.perf_counter() - started,
                   "cpu_seconds": time.process_time() - cpu,
                   "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--produced", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.produced, args.out)


if __name__ == "__main__":
    main()
