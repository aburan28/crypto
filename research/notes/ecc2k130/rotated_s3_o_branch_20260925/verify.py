#!/usr/bin/env python3
"""Independent full-point and parsed-DIMACS replay of O-aware S3 export."""
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
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
PRIOR = NOTES / "rotated_s3_candidate_20260925/evidence/n13-m5/producer/result.json"
PARENT = NOTES / "rotated_subspace_support_20260925/gate.py"
GATE_INPUTS = NOTES / "rotated_m56_export_gate_20260925/INPUTS.json"
GATE_FROZEN = NOTES / "rotated_m56_export_gate_20260925/FROZEN.json"
GATE_N13 = NOTES / "rotated_m56_export_gate_20260925/evidence/n13-m5.json"
PANELS = (("n2-m4", 2, 4, 0x7), ("n3-m5", 3, 5, 0xb),
          ("n4-m4", 4, 4, 0x13), ("n13-m5", 13, 5, 0x201b))
CAP_SECONDS = 180
CAP_RSS = 512 * 1024 * 1024


class ClauseError(AssertionError):
    pass


class TargetError(AssertionError):
    pass


class DomainError(AssertionError):
    pass


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def parent_module():
    spec = importlib.util.spec_from_file_location("s3_o_independent_curve", PARENT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def poly_rem(a: int, b: int) -> int:
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_rem(a, b)
    return a


def verify_field(f):
    n, poly = f.n, f.poly
    assert poly.bit_length() == n + 1 and poly & 1
    x = poly_rem(2, poly)
    values = [x]
    for _ in range(n):
        values.append(f.square(values[-1]))
    assert values[n] == x
    primes = [p for p in range(2, n + 1) if n % p == 0 and
              all(p % d for d in range(2, int(p ** .5) + 1))]
    assert all(poly_gcd(poly, values[n // p] ^ x) == 1 for p in primes)


def expected_order(n: int) -> int:
    trace = [2, -1]
    for _ in range(2, n + 1):
        trace.append(-trace[-1] - 2 * trace[-2])
    return (1 << n) + 1 - trace[n]


def lift_map(f):
    mapping = {}
    for z in range(1 << f.n):
        mapping.setdefault(f.square(z) ^ z, []).append(z)
    assert all(len(row) == 2 for row in mapping.values())
    return mapping


def fibre(f, x: int, by_h: dict):
    assert 0 <= x < 1 << f.n
    if x == 0:
        return ((0, 1),)
    h = x ^ f.square(f.inverse(x))
    return tuple(sorted((x, f.mul(x, z)) for z in by_h.get(h, ())))


def input_rows(curve, by_h, n: int, m: int):
    f = curve.f
    if n != 13:
        points = sorted((x, y) for x in range(1 << n) for y in range(1 << n)
                        if curve.on_curve((x, y)))
        assert points == sorted(p for x in range(1 << n) for p in fibre(f, x, by_h))
        assert len(points) + 1 == expected_order(n)
        factors = [points[:] for _ in range(m)]
        targets = [{"id": "O", "point": None, "class": "identity"}]
        targets += [{"id": f"p-{x}-{y}", "point": [x, y], "class": "all_points"}
                    for x, y in points]
        return factors, targets
    with tarfile.open(CORPUS, "r:gz") as tar:
        content = tar.extractfile("raw/n13-m5/factors.json").read()
    factors = [[tuple(p) for p in slot] for slot in json.loads(content)]
    prior = json.loads(PRIOR.read_text())
    targets = [{"id": f"Q{row['Q_index']}T{row['T_index']}",
                "point": row["target"], "class": row["target_class"]}
               for row in prior["targets"]]
    assert len(targets) == 32 and targets[12]["point"] == [7256, 3272]
    targets.append({"id": "O", "point": None, "class": "identity"})
    return factors, targets


def complete_fibres(factors, f, by_h):
    if not factors or any(not slot for slot in factors):
        raise DomainError("empty factor slot")
    domains = []
    for slot in factors:
        actual = {tuple(p) for p in slot}
        xs = sorted({p[0] for p in actual})
        for x in xs:
            if not (0 <= x < 1 << f.n):
                raise DomainError("factor x outside field")
            expected = set(fibre(f, x, by_h))
            if not expected:
                raise DomainError(f"nonlift factor x={x}")
            if {p for p in actual if p[0] == x} != expected:
                raise DomainError(f"incomplete fibre x={x}")
        domains.append(xs)
    return domains


def value_point_set(value, f, by_h):
    return (None,) if value is None else fibre(f, value, by_h)


def local_geometry(curve, by_h, left, right):
    first = value_point_set(left, curve.f, by_h)
    second = fibre(curve.f, right, by_h)
    assert first and second
    return {None if p is None else p[0] for p in
            (curve.add(a, b) for a in first for b in second)}


def parse_cnf(path: Path):
    header = None
    clauses = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("c"):
            continue
        if line.startswith("p "):
            assert header is None
            parts = line.split()
            assert len(parts) == 4 and parts[1] == "cnf"
            header = (int(parts[2]), int(parts[3]))
            continue
        values = [int(v) for v in line.split()]
        assert values and values[-1] == 0 and 0 not in values[:-1]
        clauses.append(tuple(values[:-1]))
    assert header is not None and len(clauses) == header[1]
    assert all(0 < abs(v) <= header[0] for row in clauses for v in row)
    return header, clauses


def expected_stages(m: int):
    return ([{"left": "F0", "right": "F1", "out": "S2"}] +
            [{"left": f"S{i}", "right": f"F{i}", "out": f"S{i+1}"}
             for i in range(2, m)])


def check_schema(schema, name: str, n: int, m: int, poly: int, factors,
                 targets, domains, f, by_h):
    assert (schema["panel"], schema["field_degree"], schema["m"], schema["field_poly"]) == (name, n, m, poly)
    assert schema["factor_points"] == [[list(p) for p in slot] for slot in factors]
    assert schema["stages"] == expected_stages(m)
    groups = schema["groups"]
    assert [g["name"] for g in groups] == ([f"F{i}" for i in range(m)] +
                                          [f"S{i}" for i in range(2, m + 1)])
    next_var = 1
    for i, group in enumerate(groups):
        values = group["values"]
        assert values and len(values) == len(set(values))
        if i < m:
            assert values == domains[i]
        else:
            assert values == sorted(values, key=lambda x: -1 if x is None else x)
            assert values[0] is None
            assert all(value is None or fibre(f, value, by_h) for value in values)
        assert group["vars"] == list(range(next_var, next_var + len(values)))
        next_var += len(values)
    assert schema["variables"] == next_var - 1
    assert len(schema["targets"]) == len(targets)
    assert all({key: row[key] for key in ("id", "point", "class")} == target
               for row, target in zip(schema["targets"], targets))
    assert all(target["point"] is None or tuple(target["point"]) in
               fibre(f, target["point"][0], by_h) for target in targets)
    final = {value: var for value, var in zip(groups[-1]["values"], groups[-1]["vars"])}
    for target in schema["targets"]:
        state = None if target["point"] is None else target["point"][0]
        if target["assumption_literal"] != final[state]:
            raise TargetError(f"wrong target assumption {target['id']}")
    return {g["name"]: g for g in groups}


def expected_clauses(schema, by_name, curve, by_h):
    output = []
    local_cases = Counter()
    for group in schema["groups"]:
        variables = group["vars"]
        output.append(tuple(variables))
        output.extend((-a, -b) for a, b in itertools.combinations(variables, 2))
    for stage in schema["stages"]:
        left, right, end = (by_name[stage[key]] for key in ("left", "right", "out"))
        for u, vu in zip(left["values"], left["vars"]):
            for a, va in zip(right["values"], right["vars"]):
                permitted = local_geometry(curve, by_h, u, a)
                local_cases["O_plus_factor" if u is None else
                            "equal_x" if u == a else "distinct_x"] += 1
                for v, vv in zip(end["values"], end["vars"]):
                    if v not in permitted:
                        output.append((-vu, -va, -vv))
    return output, local_cases


def assert_clauses(actual, expected):
    if actual != expected:
        first = next((i for i, (a, e) in enumerate(zip(actual, expected)) if a != e),
                     min(len(actual), len(expected)))
        raise ClauseError(f"clause mismatch at index {first}")


def parsed_relation(schema, by_name, clauses):
    forbidden = set(row for row in clauses if len(row) == 3 and all(v < 0 for v in row))
    tables = []
    for stage in schema["stages"]:
        left, right, end = (by_name[stage[key]] for key in ("left", "right", "out"))
        table = {}
        for u, vu in zip(left["values"], left["vars"]):
            for a, va in zip(right["values"], right["vars"]):
                table[u, a] = tuple(v for v, vv in zip(end["values"], end["vars"])
                                    if (-vu, -va, -vv) not in forbidden)
        tables.append(table)
    return tables


def candidate_paths(domains, tables):
    rows = []
    for factors in itertools.product(*domains):
        paths = [(state,) for state in tables[0][factors[0], factors[1]]]
        for stage, factor in enumerate(factors[2:], start=1):
            paths = [path + (state,) for path in paths
                     for state in tables[stage][path[-1], factor]]
        rows.extend((factors, path) for path in paths)
    return rows


def point_oracle(curve, factors):
    geometry = {}
    total = o_prefix = o_terminal = 0
    branch_counts = Counter()
    for choice in itertools.product(*factors):
        total += 1
        running = choice[0]
        states = []
        for factor in choice[1:]:
            if running is None:
                branch_counts["O_plus_factor"] += 1
            elif running[0] == factor[0] == 0:
                branch_counts["zero_zero_to_O"] += 1
            elif running[0] == factor[0]:
                branch_counts["equal_x_doubling_or_inverse"] += 1
            else:
                branch_counts["distinct_x"] += 1
            running = curve.add(running, factor)
            states.append(None if running is None else running[0])
        o_prefix += None in states[:-1]
        o_terminal += running is None
        key = (tuple(p[0] for p in choice), tuple(states))
        geometry.setdefault(key, set()).add(running)
    return geometry, {"signed_point_tuples": total,
                      "signed_tuples_with_O_prefix": o_prefix,
                      "signed_tuples_with_O_terminal": o_terminal,
                      "point_addition_cases": dict(sorted(branch_counts.items()))}


def read_rows(path: Path):
    with gzip.open(path, "rt") as stream:
        return [json.loads(line) for line in stream]


def target_statistics(targets, candidate, geometry, schema, mask_by_x):
    output = []
    for target, source in zip(schema["targets"], targets):
        point = None if source["point"] is None else tuple(source["point"])
        state = None if point is None else point[0]
        selected = [key for key in candidate if key[1][-1] == state]
        exact = [key for key, points in geometry.items() if point in points]
        assert set(selected) == set(exact), source["id"]
        masks = {key[0] for key in selected}
        affine = {key[0] for key in selected if None not in key[1][:-1]}
        exceptional = sorted(masks - affine)
        expected = {"model_path_count": len(selected), "model_x_tuple_count": len(masks),
                    "o_prefix_path_count": sum(None in key[1][:-1] for key in selected),
                    "exceptional_only_x_tuples": [list(xs) for xs in exceptional]}
        if mask_by_x is not None:
            expected["exceptional_only_masks"] = [[mask_by_x[i][x]
                                                     for i, x in enumerate(xs)]
                                                    for xs in exceptional]
        assert all(target[key] == value for key, value in expected.items()), source["id"]
        output.append({"id": source["id"], "target": source["point"],
                       **expected})
    return output


def verify_mask_map(schema, parent, factors, n):
    if n != 13:
        assert schema["mask_by_x"] is None
        return None
    field = parent.Field(13, [0, 1, 3, 4])
    bases = parent.subspace_basis(parent.normal_conjugates(field, 3), 5, 2)
    maps = []
    for basis in bases:
        mapping = {}
        for mask in range(4):
            x = 0
            for bit, element in enumerate(basis):
                if mask >> bit & 1:
                    x ^= element
            assert x not in mapping
            mapping[x] = mask
        maps.append(mapping)
    assert schema["mask_by_x"] == [{str(x): mask for x, mask in row.items()} for row in maps]
    assert all({p[0] for p in slot} <= set(row) for slot, row in zip(factors, maps))
    return maps


def verify_panel(name: str, n: int, m: int, poly: int, output: Path,
                 producer_row: dict, parent):
    f = parent.Field(n, [bit for bit in range(n) if poly >> bit & 1])
    assert f.poly == poly
    verify_field(f)
    curve = parent.Curve(f)
    by_h = lift_map(f)
    factors, targets = input_rows(curve, by_h, n, m)
    domains = complete_fibres(factors, f, by_h)
    schema = json.loads((output / "schema.json").read_text())
    assert sha(output / "schema.json") == producer_row["schema_sha256"]
    assert sha(output / "base.cnf") == producer_row["base_sha256"]
    assert sha(output / "paths.jsonl.gz") == producer_row["paths_sha256"]
    by_name = check_schema(schema, name, n, m, poly, factors, targets, domains, f, by_h)
    mask_by_x = verify_mask_map(schema, parent, factors, n)
    header, actual_clauses = parse_cnf(output / "base.cnf")
    assert header == (schema["variables"], schema["clauses"])
    expected, local_cases = expected_clauses(schema, by_name, curve, by_h)
    assert_clauses(actual_clauses, expected)
    tables = parsed_relation(schema, by_name, actual_clauses)
    candidate = candidate_paths(domains, tables)
    raw_rows = read_rows(output / "paths.jsonl.gz")
    assert [{"factor_x": list(xs), "states": list(states)} for xs, states in candidate] == raw_rows
    geometry, oracle_counts = point_oracle(curve, factors)
    assert set(candidate) == set(geometry)
    assert len(candidate) == len(set(candidate))
    for (_, states), points in geometry.items():
        terminal = states[-1]
        assert points == set(value_point_set(terminal, f, by_h))
    target_rows = target_statistics(targets, candidate, geometry, schema, mask_by_x)
    if n == 13:
        assert [0, 0, 0, 2, 1] in target_rows[12]["exceptional_only_masks"]
        witness_x = (0, 0, 0, 6433, 217)
        assert any(xs == witness_x and states[0] is None and states[-1] == 7256
                   and (7256, 3272) in geometry[(xs, states)]
                   for xs, states in candidate)
    else:
        assert target_rows[0]["model_path_count"] > 0
    assert producer_row["factor_x_tuples"] == math.prod(map(len, domains))
    assert producer_row["signed_point_tuples"] == oracle_counts["signed_point_tuples"]
    assert producer_row["model_paths"] == len(candidate)
    assert producer_row["variables"] == header[0]
    assert producer_row["clauses"] == header[1]
    assert producer_row["targets"] == len(targets)
    assert producer_row["target_rows"] == schema["targets"]
    assert producer_row["transition_cases"] == dict(sorted(local_cases.items()))
    return {"panel": name, "candidate_paths": len(candidate),
            "target_labels": len(targets), "variables": header[0],
            "clauses": header[1], "oracle": oracle_counts,
            "target_rows": target_rows,
            "field_operations": dict(sorted(f.operations.items())),
            "curve_operations": dict(sorted(curve.operations.items()))}


def negative_controls(parent, producer_dir: Path):
    f3 = parent.Field(3, [0, 1])
    f2 = parent.Field(2, [0, 1])
    try:
        complete_fibres([[(2, 0)]], f3, lift_map(f3))
    except DomainError as error:
        nonlift = f"{type(error).__name__}: {error}"
    else:
        raise AssertionError("nonlift control accepted")
    try:
        complete_fibres([[fibre(f2, 1, lift_map(f2))[0]]], f2, lift_map(f2))
    except DomainError as error:
        sign = f"{type(error).__name__}: {error}"
    else:
        raise AssertionError("sign-incomplete control accepted")
    schema = json.loads((producer_dir / "n2-m4/schema.json").read_text())
    header, clauses = parse_cnf(producer_dir / "n2-m4/base.cnf")
    by_name = {g["name"]: g for g in schema["groups"]}
    f0 = by_name["F0"]["vars"][by_name["F0"]["values"].index(0)]
    f1 = by_name["F1"]["vars"][by_name["F1"]["values"].index(0)]
    o = by_name["S2"]["vars"][by_name["S2"]["values"].index(None)]
    mutation = (-f0, -f1, -o)
    assert mutation not in clauses
    try:
        assert_clauses(clauses + [mutation], clauses)
    except ClauseError as error:
        mutated_o = f"{type(error).__name__}: {error}"
    else:
        raise AssertionError("mutated O clause accepted")
    o_target = next(row for row in schema["targets"] if row["point"] is None)
    final = by_name["S4"]
    finite_zero = final["vars"][final["values"].index(0)]
    assert finite_zero != o_target["assumption_literal"]
    altered = json.loads(json.dumps(schema))
    next(row for row in altered["targets"] if row["point"] is None)["assumption_literal"] = finite_zero
    factors, targets = input_rows(parent.Curve(f2), lift_map(f2), 2, 4)
    try:
        check_schema(altered, "n2-m4", 2, 4, 0x7, factors, targets,
                     complete_fibres(factors, f2, lift_map(f2)), f2, lift_map(f2))
    except TargetError as error:
        wrong_target = f"{type(error).__name__}: {error}"
    else:
        raise AssertionError("wrong O target literal accepted")
    return {"nonlift": nonlift, "sign_incomplete": sign,
            "mutated_O_clause": mutated_o, "wrong_O_target_literal": wrong_target}


def replay(producer_dir: Path):
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == frozen["verify_sha256"]
    assert sha(PARENT) == frozen["parent_sha256"]
    assert sha(CORPUS) == frozen["corpus_sha256"]
    assert sha(PRIOR) == frozen["prior_sha256"]
    assert sha(GATE_INPUTS) == frozen["gate_inputs_sha256"]
    assert sha(GATE_FROZEN) == frozen["gate_frozen_sha256"]
    assert sha(GATE_N13) == frozen["gate_n13_sha256"]
    producer = json.loads((producer_dir / "result.json").read_text())
    assert producer["domain"] == frozen["domain"]
    parent = parent_module()
    panels = []
    for (name, n, m, poly), produced in zip(PANELS, producer["panels"]):
        assert produced["panel"] == name
        panels.append(verify_panel(name, n, m, poly, producer_dir / name,
                                   produced, parent))
    assert len(panels) == len(PANELS)
    controls = negative_controls(parent, producer_dir)
    assert producer["negative_input_controls"] == {
        "nonlift": "nonlift factor x=2",
        "sign_incomplete": "sign-incomplete or invalid factor x=1"}
    return {"decision": "PASS", "domain": frozen["domain"],
            "producer_sha256": sha(producer_dir / "result.json"),
            "panels": panels, "negative_controls": controls}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"{CAP_SECONDS}s O-branch independent replay cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        result = replay(args.producer)
        result.update({"wall_seconds": time.perf_counter() - started,
                       "cpu_seconds": time.process_time() - cpu,
                       "peak_rss_bytes": rss()})
        assert result["wall_seconds"] <= CAP_SECONDS and result["peak_rss_bytes"] <= CAP_RSS
        save(args.out, result)
    except Exception as error:
        save(args.out, {"decision": "FAIL", "error": repr(error),
                        "wall_seconds": time.perf_counter() - started,
                        "cpu_seconds": time.process_time() - cpu,
                        "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    main()
