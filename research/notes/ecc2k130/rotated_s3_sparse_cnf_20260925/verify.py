#!/usr/bin/env python3
"""Independent parsed-CNF, dense-reference and signed-point replay for sparse S3."""
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
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DENSE = HERE.parent / "rotated_s3_o_branch_20260925"
PANELS = (("n2-m4", 2, 4, 0x7), ("n3-m5", 3, 5, 0xb),
          ("n4-m4", 4, 4, 0x13), ("n13-m5", 13, 5, 0x201b))
CAP_SECONDS = 180
CAP_RSS = 512 * 1024 * 1024


class ClauseError(AssertionError):
    pass


class TargetError(AssertionError):
    pass


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def independent_module():
    spec = importlib.util.spec_from_file_location("sparse_independent_dense_verifier", DENSE / "verify.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CNFReader:
    """Streaming DIMACS parser, also used for the million-clause dense reference."""
    def __init__(self, path: Path):
        self.stream = path.open("r")
        self.variables = self.declared = self.seen = 0
        while True:
            line = self.stream.readline()
            if not line:
                raise ClauseError("missing DIMACS header")
            if not line.strip() or line.startswith("c"):
                continue
            parts = line.split()
            if len(parts) != 4 or parts[:2] != ["p", "cnf"]:
                raise ClauseError("invalid DIMACS header")
            self.variables, self.declared = map(int, parts[2:])
            break
        self.pending = self._read()

    def _read(self):
        for line in self.stream:
            if not line.strip() or line.startswith("c"):
                continue
            if line.startswith("p "):
                raise ClauseError("duplicate DIMACS header")
            try:
                values = tuple(map(int, line.split()))
            except ValueError as error:
                raise ClauseError("noninteger DIMACS literal") from error
            if not values or values[-1] != 0 or 0 in values[:-1]:
                raise ClauseError("invalid DIMACS clause terminator")
            row = values[:-1]
            if any(not 1 <= abs(value) <= self.variables for value in row):
                raise ClauseError("literal outside declared range")
            return row
        return None

    def take(self):
        if self.pending is None:
            raise ClauseError("missing DIMACS clause")
        row = self.pending
        self.seen += 1
        self.pending = self._read()
        return row

    def finish(self):
        if self.pending is not None or self.seen != self.declared:
            raise ClauseError("DIMACS clause count or order mismatch")
        self.stream.close()


def require(actual, expected, label):
    if actual != expected:
        raise ClauseError(f"{label}: {actual!r} != {expected!r}")


def expected_block(primary, auxiliary):
    """Literal order is fixed separately from the producer's loop."""
    k = len(primary)
    result = [tuple(primary)]
    if k == 1:
        return result
    result.append((-primary[0], auxiliary[0]))
    result.extend(row for j in range(1, k - 1)
                  for row in ((-primary[j], auxiliary[j]),
                              (-auxiliary[j - 1], auxiliary[j]),
                              (-primary[j], -auxiliary[j - 1])))
    result.append((-primary[-1], -auxiliary[-1]))
    return result


def satisfies(clause, true_variables: set[int]):
    return any((literal > 0) == (abs(literal) in true_variables) for literal in clause)


def check_onehot_block(primary, auxiliary, block):
    k = len(primary)
    assert block == expected_block(primary, auxiliary)
    if k <= 5:
        for chosen in range(1 << k):
            primary_truth = {primary[j] for j in range(k) if chosen >> j & 1}
            sat_extensions = 0
            for aux_mask in range(1 << len(auxiliary)):
                truth = primary_truth | {auxiliary[j] for j in range(len(auxiliary))
                                         if aux_mask >> j & 1}
                sat_extensions += all(satisfies(row, truth) for row in block)
            assert sat_extensions == (1 if chosen.bit_count() == 1 else 0)
    else:
        for chosen in range(k):
            truth = {primary[chosen]} | {auxiliary[j] for j in range(k - 1)
                                         if chosen <= j}
            assert all(satisfies(row, truth) for row in block)


def check_schema(schema, baseline, frozen, name, iv, factors, targets, domains, f, by_h):
    for key, value in baseline.items():
        if key not in ("variables", "clauses", "targets"):
            assert schema[key] == value, (name, key)
    assert len(schema["targets"]) == len(baseline["targets"])
    for target, reference in zip(schema["targets"], baseline["targets"]):
        assert {key: value for key, value in target.items() if key != "assumption_literal"} == {key: value for key, value in reference.items() if key != "assumption_literal"}
    assert schema["panel"] == name
    assert schema["encoding"] == "sinz-sequential-amo-and-allowed-output-v1"
    assert schema["primary_variables"] == baseline["variables"]
    assert schema["dense_schema_sha256"] == frozen["panels"][name]["schema.json"]["sha256"]
    assert schema["dense_base_sha256"] == frozen["panels"][name]["base.cnf"]["sha256"]
    assert schema["dense_paths_sha256"] == frozen["panels"][name]["paths.jsonl.gz"]["sha256"]
    aux = schema["auxiliary_groups"]
    assert [row["name"] for row in aux] == [group["name"] for group in baseline["groups"]]
    next_var = baseline["variables"] + 1
    for group, row in zip(baseline["groups"], aux):
        count = len(group["vars"]) - 1
        assert row["vars"] == list(range(next_var, next_var + count))
        next_var += count
    assert schema["variables"] == next_var - 1
    by_name = iv.check_schema(baseline, name, baseline["field_degree"],
                              baseline["m"], baseline["field_poly"],
                              factors, targets, domains, f, by_h)
    for target in schema["targets"]:
        state = None if target["point"] is None else target["point"][0]
        final = by_name[baseline["stages"][-1]["out"]]
        expected = final["vars"][final["values"].index(state)]
        if target["assumption_literal"] != expected:
            raise TargetError(f"wrong target assumption {target['id']}")
        assert target["assumption_literal"] <= baseline["variables"]
    return by_name


def parse_dense(path: Path, baseline):
    reader = CNFReader(path)
    require((reader.variables, reader.declared),
            (baseline["variables"], baseline["clauses"]), "dense header")
    for group in baseline["groups"]:
        xs = group["vars"]
        require(reader.take(), tuple(xs), "dense at least one")
        for a, b in itertools.combinations(xs, 2):
            require(reader.take(), (-a, -b), "dense at most one")
    by_name = {group["name"]: group for group in baseline["groups"]}
    tables = []
    for stage in baseline["stages"]:
        left, right, end = (by_name[stage[key]] for key in ("left", "right", "out"))
        table = {}
        output_set = set(end["vars"])
        for u, vu in zip(left["values"], left["vars"]):
            for a, va in zip(right["values"], right["vars"]):
                forbidden = set()
                previous = 0
                while reader.pending is not None and len(reader.pending) == 3 and \
                        reader.pending[:2] == (-vu, -va):
                    row = reader.take()
                    vv = -row[2]
                    if vv not in output_set or vv <= previous:
                        raise ClauseError("dense forbidden output variable/order")
                    forbidden.add(vv)
                    previous = vv
                table[u, a] = tuple(v for v, vv in zip(end["values"], end["vars"])
                                    if vv not in forbidden)
        tables.append(table)
    reader.finish()
    return tables


def parse_sparse(path: Path, schema, by_name, iv, curve, by_h):
    reader = CNFReader(path)
    require((reader.variables, reader.declared),
            (schema["variables"], schema["clauses"]), "sparse header")
    onehot = 0
    blocks = {}
    for group, aux in zip(schema["groups"], schema["auxiliary_groups"]):
        expected = expected_block(group["vars"], aux["vars"])
        actual = [reader.take() for _ in expected]
        if actual != expected:
            first = next(i for i, pair in enumerate(zip(actual, expected)) if pair[0] != pair[1])
            raise ClauseError(f"sequential one-hot clause {group['name']}:{first}")
        check_onehot_block(group["vars"], aux["vars"], actual)
        blocks[group["name"]] = actual
        onehot += len(actual)
    require(onehot, schema["onehot_clauses"], "onehot count")
    tables = []
    implication = 0
    for stage in schema["stages"]:
        left, right, end = (by_name[stage[key]] for key in ("left", "right", "out"))
        table = {}
        end_by_literal = dict(zip(end["vars"], end["values"]))
        for u, vu in zip(left["values"], left["vars"]):
            for a, va in zip(right["values"], right["vars"]):
                row = reader.take()
                if row[:2] != (-vu, -va) or len(row) < 3:
                    raise ClauseError("sparse input-pair implication")
                outputs = row[2:]
                if len(set(outputs)) != len(outputs) or any(v not in end_by_literal for v in outputs):
                    raise ClauseError("sparse output literals")
                if tuple(sorted(outputs)) != outputs:
                    raise ClauseError("sparse output order")
                selected = tuple(end_by_literal[v] for v in outputs)
                table[u, a] = selected
                expected_set = iv.local_geometry(curve, by_h, u, a)
                expected = tuple(v for v in end["values"] if v in expected_set)
                if selected != expected:
                    raise ClauseError(f"point-law output mismatch {stage} {(u, a)}")
                implication += 1
        tables.append(table)
    require(implication, schema["implication_clauses"], "implication count")
    reader.finish()
    return tables, blocks


def path_truth(schema, path):
    xs, states = path
    m = len(xs)
    choices = list(xs) + list(states)
    assert len(choices) == len(schema["groups"]) and len(states) == m - 1
    truth = set()
    for group, aux, value in zip(schema["groups"], schema["auxiliary_groups"], choices):
        index = group["values"].index(value)
        truth.add(group["vars"][index])
        truth.update(aux["vars"][j] for j in range(len(aux["vars"])) if index <= j)
    return truth


def read_paths(path: Path):
    with gzip.open(path, "rt") as file:
        return [json.loads(line) for line in file]


def verify_panel(panel, sparse_dir: Path, produced, frozen, iv, parent):
    name, n, m, poly = panel
    dense_dir = DENSE / "evidence/producer" / name
    for file in ("base.cnf", "schema.json", "paths.jsonl.gz"):
        assert sha(dense_dir / file) == frozen["panels"][name][file]["sha256"]
    base = json.loads((dense_dir / "schema.json").read_text())
    schema = json.loads((sparse_dir / "schema.json").read_text())
    assert sha(sparse_dir / "schema.json") == produced["schema_sha256"]
    assert sha(sparse_dir / "base.cnf") == produced["base_sha256"]
    assert sha(sparse_dir / "paths.jsonl.gz") == produced["paths_sha256"]
    f = parent.Field(n, [bit for bit in range(n) if poly >> bit & 1])
    iv.verify_field(f)
    curve = parent.Curve(f)
    by_h = iv.lift_map(f)
    factors, targets = iv.input_rows(curve, by_h, n, m)
    domains = iv.complete_fibres(factors, f, by_h)
    by_name = check_schema(schema, base, frozen, name, iv,
                           factors, targets, domains, f, by_h)
    dense_tables = parse_dense(dense_dir / "base.cnf", base)
    sparse_tables, blocks = parse_sparse(sparse_dir / "base.cnf", schema,
                                         by_name, iv, curve, by_h)
    if sparse_tables != dense_tables:
        raise ClauseError(f"sparse/dense local truth mismatch: {name}")
    sparse_paths = iv.candidate_paths(domains, sparse_tables)
    dense_paths = iv.candidate_paths(domains, dense_tables)
    assert sparse_paths == dense_paths
    dense_archive = read_paths(dense_dir / "paths.jsonl.gz")
    normalized = [{"factor_x": list(xs), "states": list(states)}
                  for xs, states in sparse_paths]
    assert normalized == dense_archive == read_paths(sparse_dir / "paths.jsonl.gz")
    geometry, oracle_counts = iv.point_oracle(curve, factors)
    assert set(sparse_paths) == set(geometry) and len(sparse_paths) == len(geometry)
    for (_, states), points in geometry.items():
        assert points == set(iv.value_point_set(states[-1], f, by_h))
    mask_by_x = iv.verify_mask_map(schema, parent, factors, n)
    target_rows = iv.target_statistics(targets, sparse_paths, geometry, schema, mask_by_x)
    if n == 13:
        assert [0, 0, 0, 2, 1] in target_rows[12]["exceptional_only_masks"]
        witness = ((0, 0, 0, 6433, 217),)
        assert any(xs == witness[0] and states[0] is None and states[-1] == 7256
                   and (7256, 3272) in geometry[(xs, states)]
                   for xs, states in sparse_paths)
        assert (produced["variables"], produced["clauses"]) == (2517, 4731)
        assert produced["bytes"] <= 250000 and produced["clauses"] * 100 <= produced["dense_clauses"]
    sparse_clauses = []
    reader = CNFReader(sparse_dir / "base.cnf")
    while reader.pending is not None:
        sparse_clauses.append(reader.take())
    reader.finish()
    for path in sparse_paths:
        truth = path_truth(schema, path)
        assert all(satisfies(row, truth) for row in sparse_clauses)
    assert produced["model_paths"] == len(sparse_paths)
    assert produced["target_labels"] == len(targets)
    assert produced["factor_x_tuples"] == math.prod(map(len, domains))
    assert produced["variables"] == schema["variables"]
    assert produced["clauses"] == schema["clauses"]
    assert produced["bytes"] == (sparse_dir / "base.cnf").stat().st_size
    assert produced["dense_bytes"] == (dense_dir / "base.cnf").stat().st_size
    return {"panel": name, "candidate_paths": len(sparse_paths),
            "signed_point_tuples": oracle_counts["signed_point_tuples"],
            "target_labels": len(targets), "variables": schema["variables"],
            "clauses": schema["clauses"], "bytes": produced["bytes"],
            "dense_variables": base["variables"], "dense_clauses": base["clauses"],
            "dense_bytes": produced["dense_bytes"], "target_rows": target_rows,
            "point_cases": oracle_counts["point_addition_cases"]}, schema, sparse_clauses, blocks


def assert_exact_clauses(actual, expected):
    if len(actual) != len(expected):
        raise ClauseError(f"clause count {len(actual)} != {len(expected)}")
    for index, (row, reference) in enumerate(zip(actual, expected)):
        if row != reference:
            raise ClauseError(f"clause mismatch at index {index}")


def negative_controls(iv, parent, schema, clauses, blocks):
    by_name = {g["name"]: g for g in schema["groups"]}
    f0, f1, s2, final = (by_name[name] for name in ("F0", "F1", "S2", "S4"))
    v0 = f0["vars"][f0["values"].index(0)]
    v1 = f1["vars"][f1["values"].index(0)]
    o = s2["vars"][s2["values"].index(None)]
    zero = s2["vars"][s2["values"].index(0)]
    allowed = (-v0, -v1, o)
    assert clauses.count(allowed) == 1
    archived = read_paths(DENSE / "evidence/producer/n2-m4/paths.jsonl.gz")
    o_prefix = next(row for row in archived if row["factor_x"][:2] == [0, 0]
                    and row["states"][0] is None)
    o_prefix_truth = path_truth(schema, (tuple(o_prefix["factor_x"]),
                                         tuple(o_prefix["states"])))
    assert all(satisfies(row, o_prefix_truth) for row in clauses)
    results = {}
    mutation = (-v0, -v1, -o)
    assert mutation not in clauses and not satisfies(mutation, o_prefix_truth)
    try:
        assert_exact_clauses(clauses + [mutation], clauses)
    except ClauseError as error:
        results["forbid_zero_zero_O"] = str(error)
    else:
        raise AssertionError("extra O-forbidding clause accepted")
    altered = clauses[:]
    altered[altered.index(allowed)] = (-v0, -v1, zero)
    assert not satisfies((-v0, -v1, zero), o_prefix_truth)
    try:
        assert_exact_clauses(altered, clauses)
    except ClauseError as error:
        results["replace_O_output"] = str(error)
    else:
        raise AssertionError("changed O output accepted")
    group = f0
    first = (-group["vars"][0], schema["auxiliary_groups"][0]["vars"][0])
    block = blocks["F0"]
    assert first in block
    shortened = block[:]
    shortened.remove(first)
    witness_found = False
    for pair in itertools.combinations(range(len(group["vars"])), 2):
        primary_truth = {group["vars"][index] for index in pair}
        for mask in range(1 << len(schema["auxiliary_groups"][0]["vars"])):
            truth = primary_truth | {v for j, v in enumerate(schema["auxiliary_groups"][0]["vars"])
                                     if mask >> j & 1}
            witness_found |= all(satisfies(row, truth) for row in shortened)
    assert witness_found
    try:
        assert_exact_clauses(shortened, block)
    except ClauseError as error:
        results["drop_sequential_clause"] = str(error)
    else:
        raise AssertionError("dropped AMO clause accepted")
    changed = json.loads(json.dumps(schema))
    o_target = next(row for row in changed["targets"] if row["point"] is None)
    o_target["assumption_literal"] = final["vars"][final["values"].index(0)]
    o_terminal = next(row for row in archived if row["states"][-1] is None)
    o_terminal_truth = path_truth(schema, (tuple(o_terminal["factor_x"]),
                                           tuple(o_terminal["states"])))
    assert final["vars"][final["values"].index(None)] in o_terminal_truth
    assert o_target["assumption_literal"] not in o_terminal_truth
    f = parent.Field(2, [0, 1])
    by_h = iv.lift_map(f)
    factors, targets = iv.input_rows(parent.Curve(f), by_h, 2, 4)
    domains = iv.complete_fibres(factors, f, by_h)
    base = json.loads((DENSE / "evidence/producer/n2-m4/schema.json").read_text())
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    try:
        check_schema(changed, base, frozen, "n2-m4", iv,
                     factors, targets, domains, f, by_h)
    except (TargetError, AssertionError) as error:
        results["wrong_O_target"] = f"{type(error).__name__}: {error}"
    else:
        raise AssertionError("wrong O target accepted")
    return results


def run(sparse_dir: Path, out: Path):
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"{CAP_SECONDS}s sparse independent verifier cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        assert sha(Path(__file__)) == frozen["verify_sha256"]
        assert sha(DENSE / "verify.py") == frozen["dense_verify_sha256"]
        assert sha(DENSE / "evidence/receipt.json") == frozen["dense_receipt_sha256"]
        assert sha(DENSE / "evidence/producer/result.json") == frozen["dense_result_sha256"]
        for relative, digest in frozen["dense_dependency_sha256"].items():
            assert sha(HERE.parent / relative) == digest, relative
        produced = json.loads((sparse_dir / "result.json").read_text())
        assert produced["domain"] == frozen["domain"]
        iv = independent_module()
        parent = iv.parent_module()
        panels = []
        control_data = None
        for panel, row in zip(PANELS, produced["panels"]):
            assert row["panel"] == panel[0]
            checked, schema, clauses, blocks = verify_panel(panel, sparse_dir / panel[0],
                                                              row, frozen, iv, parent)
            panels.append(checked)
            if panel[0] == "n2-m4":
                control_data = (schema, clauses, blocks)
        assert len(panels) == len(PANELS) and control_data is not None
        controls = negative_controls(iv, parent, *control_data)
        assert set(controls) == {"forbid_zero_zero_O", "replace_O_output",
                                 "drop_sequential_clause", "wrong_O_target"}
        result = {"decision": "PASS", "domain": frozen["domain"],
                  "producer_sha256": sha(sparse_dir / "result.json"),
                  "panels": panels, "negative_controls": controls,
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu,
                  "peak_rss_bytes": rss()}
        assert result["wall_seconds"] <= CAP_SECONDS and result["peak_rss_bytes"] <= CAP_RSS
        save(out, result)
    except Exception as error:
        save(out, {"decision": "FAIL", "error": repr(error),
                   "wall_seconds": time.perf_counter() - started,
                   "cpu_seconds": time.process_time() - cpu,
                   "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sparse", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.sparse, args.out)


if __name__ == "__main__":
    main()
